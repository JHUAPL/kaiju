module init

    use gamtypes
    use math
    use gridutils
    use gamutils
    use mhdgroup
    use xml_input
    use gioH5
    use prob
    use output
    use bcs
    use quadrature
    use background
    use ringutils
    use ringrecon
    use recon
    use multifluid
    use files    
    use step
    
    implicit none

    !Initialization defaults (2D Square)
    integer :: Nc = 64
    real(rp) :: xMin = -1.0, xMax = 1.0
    integer :: seed0 = 31337 !Default random seed
    character(len=strLen) :: reconMethod
    !Max value of CFL 
    !-Can be lower from PDMB
    !-Can be forced lower via sim/CFL
    real(rp) :: MaxCFL = 0.3

    contains
    
    !Hatch Gamera
    !Initialize main data structures
    subroutine Hatch(Model,Grid,State,oState,Solver,xmlInp,userInitFunc,endTime)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State,oState
        type(Solver_T), intent(inout) :: Solver
        !OMEGA can overrule what GAMERA has
        type(XML_Input_T), intent(inout) :: xmlInp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        real(rp), optional, intent(in) :: endTime 

        !Alwasys zero for single process job
        Grid%ijkShift(1:3) = 0

        ! call appropriate subroutines to read corner info and mesh size data
        call ReadCorners(Model,Grid,xmlInp,endTime)

        ! call appropriate subroutines to calculate all appropriate grid data from the corner data
        call CalcGridInfo(Model,Grid,State,oState,Solver,xmlInp,userInitFunc)

        !Initialization complete!
        
    end subroutine Hatch

    ! Read corner data for the mesh and set grid size variables
    subroutine ReadCorners(Model,Grid,xmlInp,endTime,childGameraOpt)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(XML_Input_T), intent(inout) :: xmlInp
        real(rp), optional, intent(in) :: endTime
        logical, optional, intent(in) ::childGameraOpt

        logical :: doH5g
        character(len=strLen) :: inH5,resId
        integer :: dotLoc,nRes
        logical :: childGamera

        if(present(childGameraOpt)) then
            childGamera = childGameraOpt
        else
            childGamera = .false.
        endif

        if(.not. childGamera) then
            !Setup OMP info unless gamera is owned by someone else
            call SetOMP(xmlInp,Model%isLoud)
        endif
        !Either way set the current number of threads
        Model%nTh = NumOMP()
   
!--------------------------------------------
        !Initalize model data structure
        call initModel(Model,xmlInp)
        if(present(endTime)) Model%tFin = endTime

!--------------------------------------------
        !Initialize grid data

        !Prepare for Grid/IC generation, either from file or internal routines
        call xmlInp%Set_Val(doH5g ,"sim/doH5g" ,.false.)

        !Get input H5 if necessary
        !Restart file overwrites doH5g
        !Always read the full mesh file
        if (doH5g) call xmlInp%Set_Val(inH5,"sim/H5Grid","gMesh.h5")
        if (Model%isRestart .and. .not. childGamera) then
            !Get restart file information
            call xmlInp%Set_Val(resID,"restart/resID","msphere")
            call xmlInp%Set_Val(nRes,"restart/nRes" ,-1)
            call getRestart(Model,Grid,resId,nRes,inH5)

        endif

        !Do grid generation
        if (doH5g .or. Model%isRestart) then
            !Use H5 input for Grid
            if (Model%isLoud) write(*,*) 'Reading grid from file: ', trim(inH5)
            call readH5Grid(Model,Grid,inH5)
        else
            !Create grid (corners) from XML info
            call genGridXML(Model,Grid,xmlInp)
        endif
    end subroutine ReadCorners

    !Get name of restart file
    subroutine getRestart(Model,Grid,resId,nRes,inH5)
        type(Model_T)    , intent(inout)   :: Model
        type(Grid_T)     , intent(in)      :: Grid
        character(len=*) , intent(in)      :: resId
        integer          , intent(in)      :: nRes
        character(len=strLen), intent(out) :: inH5

        character(len=strLen) :: bStr,nStr

        !Get filename base
        if (Grid%isTiled) then
            !If case is tiled, adjust the h5 filename for each rank
            bStr = genRunId(resID,Grid%NumRi,Grid%NumRj,Grid%NumRk,Grid%Ri+1,Grid%Rj+1,Grid%Rk+1)
        else
            bStr = trim(resID)
        endif

        !Get number string
        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') nRes
        endif
        inH5 = trim(bStr) // ".gam.Res." // trim(nStr) // ".h5"
        write(*,*) 'Assigned restart file: ', trim(inH5)
        call CheckFileOrDie(inH5,"Restart file not found ...")
    end subroutine getRestart

    subroutine CalcGridInfo(Model,Grid,State,oState,Solver,xmlInp,userInitFunc)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State,oState
        type(Solver_T), intent(inout) :: Solver
        type(XML_Input_T), intent(inout) :: xmlInp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc

        character(len=strLen) :: inH5, FileCode,resId
        logical :: fExist, doReset
        real(rp) :: tReset
        integer :: dotLoc,nRes

        !Set default domains (needs to be done after grid generation/reading)
        call SetDomain(Model,Grid)

        !Figure out ring avg (if using) before C2G
        call SetRings(Model,Grid,xmlInp)

        !Turn corners into grid
        call Corners2Grid(Model,Grid)

        !Set default BCs (after C2G)
        call DefaultBCs(Model,Grid)

!--------------------------------------------
        !Initialize state data
        !First prep state, then do restart if necessary, then finish state

        call PrepState(Model,Grid,oState,State,xmlInp,userInitFunc)

        !Do background stuff
        !Incorporate background field, B0, if necessary
        if (Model%doBackground .and. Grid%doB0Init) then
            call AddB0(Model,Grid,Model%B0)
        endif
        !Incorporate gravity if necessary
        if (Model%doGrav .and. Grid%doG0Init) then
            call AddGrav(Model,Grid,Model%Phi)
        endif

        if (Model%isRestart) then
            !If restart replace State variable w/ restart file
            !Make sure inH5 is set to restart
            call xmlInp%Set_Val(resID,"restart/resID","msphere")
            call xmlInp%Set_Val(nRes,"restart/nRes" ,-1)
            call getRestart(Model,Grid,resId,nRes,inH5)

            !Test for resetting time
            call xmlInp%Set_Val(doReset ,"restart/doReset" ,.false.)
            call xmlInp%Set_Val( tReset ,"restart/tReset"   ,0.0_rp)

            !Read restart
            call readH5Restart(Model,Grid,State,oState,inH5,doReset,tReset)
        else
            ! set initial dt0 to 0, it will be set once the case settles
            Model%dt0 = 0
        endif

        if (Model%doMHD) then
            call bFlux2Fld(Model,Grid,State%magFlux,State%Bxyz)
            if (Model%isRestart) then
                call bFlux2Fld(Model,Grid,oState%magFlux,oState%Bxyz)
            else
                oState%magFlux = State%magFlux
                oState%Bxyz    = State%Bxyz
            endif
        endif


    !Finalize setup
        !Enforce initial BC's
        call Tic("BCs")
        call EnforceBCs(Model,Grid,State)
        if (Model%isRestart) then
            call EnforceBCs(Model,Grid,oState)
        else
            oState = State
        endif
        call Toc("BCs")

        !Setup timestep and initial previous state for predictor
        Model%dt = CalcDT(Model,Grid,State)
        if (.not. Model%isRestart) then
            !If restart then oState%time already set by actual value
            oState%time = State%time-Model%dt !Initial old state
        endif

        !Initialize solver data
        call initSolver(Solver, Model, Grid)

        !Setup output file
        GamH5File   = genName (Model%RunID, Grid%NumRi, Grid%NumRj, Grid%NumRk, Grid%Ri+1, Grid%Rj+1, Grid%Rk+1)
        Model%RunID = genRunId(Model%RunID, Grid%NumRi, Grid%NumRj, Grid%NumRk, Grid%Ri+1, Grid%Rj+1, Grid%Rk+1)

        if (.not. Model%isRestart) then
            !Kill output file if it exists
            call CheckAndKill(GamH5File)    
            !Write grid to output file
            call writeH5GridInit(Model,Grid)
        endif

    end subroutine CalcGridInfo

    !Prepare state and call IC
    subroutine PrepState(Model,Grid,oState,State,xmlInp,userInitFunc)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: oState,State
        type(XML_Input_T), intent(in) :: xmlInp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc

        procedure(StateIC_T), pointer :: initState => NULL()        
        logical :: doH5ic
        character(len=strLen) :: icStr
        integer :: n

        call xmlInp%Set_Val(doH5ic,"sim/doH5ic",.false.)
        !Set IC subroutine
        if (doH5ic) then
            !Use H5 input for ICs
            write(*,*) 'H5 Init currently not implemented ...'
            !initState => initH5
            stop
        else
            call xmlInp%Set_Val(icStr,"sim/icType","OT2D")
            call setIC_T(initState,icStr,userInitFunc)
        endif           

        !Setup ICs either from routine or file
        !Alloc first
        call allocState(Model,Grid,State)
        call allocState(Model,Grid,oState)

        !Call IC function
        !Call even for restart to reset BCs, background, etc ...
        call initState(Model,Grid,State,xmlInp)
        oState = State

        !Ensure the BC objects are OK
        call ValidateBCs(Model,Grid)

        !Initialize BC objects
        do n=1,Grid%NumBC
            if (allocated(Grid%externalBCs(n)%p)) then
                call Grid%externalBCs(n)%p%doInit(Model,Grid,State,xmlInp)
            endif
        enddo

    end subroutine PrepState

    !Initialize Model data structure
    subroutine initModel(Model,xmlInp)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: xmlInp

        real(rp) :: C0,MJD0
        integer :: nSeed, icSeed
        integer, dimension(:), allocatable :: vSeed
        logical :: doFatIO,doHighPrecOut

        !Start by shutting up extra ranks
        if (.not. Model%isLoud) call xmlInp%BeQuiet()

        !Get/set model defaults (can be overwritten)
        Model%nSpc = 0
        Model%nDim = 3
        Model%nG = 4
        Model%t = 0.0
        Model%ts = 0

    !Main logicals
        !These are set by default until they're implemented
        Model%doMultiF = .false.
        Model%doBackground = .false.
        Model%doHall = .false.
        Model%doGrav = .false.

        call xmlInp%Set_Val(Model%doMHD        ,'physics/doMHD'        ,.false.)
        call xmlInp%Set_Val(Model%do25D        ,'physics/do25D'        ,.false.)
        call xmlInp%Set_Val(Model%doResistive  ,'physics/doResistive'  ,.false.)

    !Misc. algorithmic/physics options
        !Need CFL & PDMB values
        !pdmb value used in module recon
        !Start by reading PDMB, set CFL based on that but limit to max 0.3
        !Set global pdmb value used in recon
        !Set Vd0, the coefficient of the diffusive electric field
        call xmlInp%Set_Val(pdmb,'sim/pdmb',0.75_rp) !Standard msphere value
        C0 = min(0.5/(pdmb+0.5) ,MaxCFL) !Set CFL based on PDM
        call SetFloors(Model,xmlInp)

        !Set CFL from XML
        call xmlInp%Set_Val(Model%CFL ,'sim/CFL'  ,C0)
        if (Model%isLoud) then
           if (Model%CFL > C0) then
               write(*,*) '-------------------------------------'
               write(*,*) 'WARNING, CFL is above critical value!'
               write(*,*) 'CFL/Critical/PDMB = ', Model%CFL,C0,pdmb
               write(*,*) '-------------------------------------'
           else
                if (Model%isLoud) write(*,*) 'CFL # = ', Model%CFL
           endif
        endif
        call xmlInp%Set_Val(Model%Vd0,'sim/Vd0',0.5)
        
        call xmlInp%Set_Val(Model%gamma,'physics/gamma',5.0_rp/3.0)
        call xmlInp%Set_Val(Model%doHogs,'physics/doHogs',.true.)
        
        if (Model%doHogs) then
            call xmlInp%Set_Val(Model%cHogH,'physics/cHogH',0.25_rp)
            call xmlInp%Set_Val(Model%cHogM,'physics/cHogM',0.25_rp)
        else
            Model%cHogH = 0.0
            Model%cHogM = 0.0
        endif
    !Time options
        !Check both omega/sim/tFin & gamera/time/tFin
        call xmlInp%Set_Val(Model%tFin,'time/tFin',1.0_rp)
        call xmlInp%Set_Val(Model%dt,'time/fixedTimestep', -1.0_rp)
        if(Model%dt > 0) then
            Model%fixedTimestep = .true.
        else
            Model%fixedTimestep = .false.
        endif

        Model%MJD0 = 0.0 !Set this by default
    
    !Timestep stuff
        !Ratio of dt0 to die
        call xmlInp%Set_Val(Model%limDT0,'timestep/limDT0',Model%limDT0)
        !Whether to try and fix low timesteps
        call xmlInp%Set_Val(Model%doCPR,'timestep/doCPR',Model%doCPR)
        if (Model%doCPR) then
            !Ratio of dt0 to start CPR
            call xmlInp%Set_Val(Model%limCPR,'timestep/limCPR',Model%limCPR)
        endif
        
    !Output/Restart (IOCLOCK)
        call Model%IO%init(xmlInp,Model%t,Model%ts)
        call xmlInp%Set_Val(Model%doDivB ,'output/DivB'    ,.true. )
        call xmlInp%Set_Val(doFatIO      ,'output/doFatIO' ,.false.)
        if (doFatIO) then
            call SetFatIO()
        endif
        call xmlInp%Set_Val(doHighPrecOut,'output/doHighPrecOut' ,.false.)
        if(doHighPrecOut) then
            call SetHPOut()
        endif
        
        !Whether to read restart
        call xmlInp%Set_Val(Model%isRestart,'restart/doRes',.false.)

    !Boris info
        call xmlInp%Set_Val(Model%doBoris,'physics/doBoris',.false.)

        if (Model%doBoris) then
            call xmlInp%Set_Val(Model%Ca,'physics/Ca',HUGE)
        else
            Model%Ca = HUGE
        endif

    !Multi-fluid
        call xmlInp%Set_Val(Model%doMultiF,'multifluid/doMF',.false.)
        if (Model%doMultiF) then
            call InitMultiF(Model,xmlInp)
        endif
    !Source terms
        call xmlInp%Set_Val(Model%doSource,'source/doSource',.false.)

    !Get RunID
        call xmlInp%Set_Val(Model%RunID,'sim/runid',"Sim")

    !Set function pointers for algorithm details
        !Cen8LRs,Up7LRs 
        GetLRs => NULL()
        call xmlInp%Set_Val(reconMethod,"sim/rmeth","8CENT")
        
        select case (trim(toUpper(reconMethod)))
        case("7UP")
            GetLRs => Up7LRs
            if (Model%isLoud) write(*,*) 'Using 7UP Reconstruction'
        case("8CENT","8C")
            GetLRs => Cen8LRs
            if (Model%isLoud) write(*,*) 'Using 8CENT Reconstruction'
        case default
            write(*,*) "Uknown sim/rmeth, exiting ..."
            stop
        end select

        call xmlInp%Set_Val(Model%doLFMLim,"sim/lfmlim",Model%doLFMLim)

        RingLR => NULL()
        call xmlInp%Set_Val(reconMethod,"ring/rngrec","PPM")
        
        select case (trim(toUpper(reconMethod)))
        case("PLM")
            RingLR => PLM_IntLR
        case("PPM")
            RingLR => PPM_IntLR
        case("WENO")
            RingLR => WENO_IntLR
        case("PCM")
            RingLR => PCM_IntLR
        case default
            write(*,*) "Uknown sim/rmeth, exiting ..."
            stop
        end select

        !Get RNG seed and setup RNG
        call xmlInp%Set_Val(icSeed,"random/seed",seed0)

        !Setup random number generator
        call random_seed(size=nSeed)
        allocate(vSeed(nSeed))
        vSeed(:) = seed0
        call random_seed(put=vSeed)

    end subroutine initModel

    subroutine initSolver(Solver, Model, Grid)
        type(Solver_T), intent(inout) :: Solver
        type(Model_T), intent(in)     :: Model
        type(Grid_T), intent(in)      :: Grid

        ! initialize stress variables
        allocate(Solver%gFlx(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NVAR,1:NDIM,BLK:Model%nSpc))
        Solver%gFlx = 0.0
        if ( Model%doMHD  ) then
            allocate(Solver%mFlx(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM,1:NDIM))
            Solver%mFlx = 0.0
        endif

        ! initialize electric field variable
        allocate(Solver%Vf(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM))
        Solver%Vf = 0.0

        ! initialize mhd variables
        call allocState(Model,Grid,Solver%StateHf)
        allocate(Solver%dGasH(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NVAR,0:Model%nSpc) )
        if (Model%doMHD) call allocGridVec(Model,Grid,Solver%dGasM)

    end subroutine initSolver
    
    !Set default Grid domain indices (after grid generation)
    subroutine SetDomain(Model,Grid)
        type(Model_T), intent(in)   :: Model
        type(Grid_T), intent(inout) :: Grid

        !Set default domain for DT calculation
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

        if (Grid%hasLowerBC(IDIR)) Grid%isDT = Grid%is-1
        if (Grid%hasLowerBC(JDIR)) Grid%jsDT = Grid%js-1
        if (Grid%hasLowerBC(KDIR)) Grid%ksDT = Grid%ks-1

        if (Grid%hasUpperBC(IDIR)) Grid%ieDT = Grid%ie+1
        if (Grid%hasUpperBC(JDIR)) Grid%jeDT = Grid%je+1
        if (Grid%hasUpperBC(KDIR)) Grid%keDT = Grid%ke+1
    end subroutine SetDomain

    !Set default Grid domain indices (after grid generation)
    subroutine DefaultBCs(Model,Grid)
        type(Model_T), intent(in)   :: Model
        type(Grid_T), intent(inout) :: Grid

        integer :: i
        !Set default BCs to triply periodic, problem IC can over-ride
        allocate(periodicInnerIBC_T :: Grid%externalBCs(INI )%p)
        allocate(periodicOuterIBC_T :: Grid%externalBCs(OUTI)%p)
        allocate(periodicInnerJBC_T :: Grid%externalBCs(INJ )%p)
        allocate(periodicOuterJBC_T :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(OUTK)%p)
        Grid%NumBC = 6
    end subroutine DefaultBCs

    !Figure out ring avg (if using) before C2G
    subroutine SetRings(Model,Grid,xmlInp)
        type(Model_T),intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(XML_Input_T), intent(inout) :: xmlInp

        call xmlInp%Set_Val(Model%doRing,"ring/doRing",.false.)
        
        if (Model%doRing) then
            call xmlInp%Set_Val(Model%Ring%GridID,"ring/gid","NONE")
            call InitRings(Model,Grid,xmlInp)
        endif
    end subroutine SetRings

    !Generate 3d corner arrays (xyz) from XML input
    !For now only do uniform
    subroutine genGridXML(Model,Grid,xmlInp,NipIN,NjpIN,NkpIN,xyzBdsIN)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(XML_Input_T), intent(in) :: xmlInp

        !If we use MPI we will send in these params else
        ! we read them from file
        integer, optional, intent(in) :: NipIN,NjpIN,NkpIN
        real(rp),optional, intent(in) :: xyzBdsIN(6)

    
        logical :: doWarp, doCyl, doSph, doStretch
        real(rp) :: xyzBds(6)
        real(rp) :: dx,dy,dz
        integer :: i,j,k
        !Variables for warp generation
        real(rp) :: Lx,Ly,Lz, xp,yp,zp, Ax,Ay,Az,R,w0, dsp
        real(rp) :: xW,yW,zW
        real(rp) :: x1,x2,x3,x,y,z
        integer :: Nw
        

        !Start with main indices
        if(present(NipIN)) then
            Grid%Nip  = NipIN
        else 
            call xmlInp%Set_Val(Grid%Nip,"idir/N",Nc)
        end if
        if(present(NjpIN)) then
            Grid%Njp  = NjpIN
        else
            call xmlInp%Set_Val(Grid%Njp,"jdir/N",Nc) 
        end if
        if(present(NkpIN)) then
            Grid%Nkp  = NkpIN
        else 
            call xmlInp%Set_Val(Grid%Nkp,"kdir/N",1)
        end if


        if (Grid%Nkp == 1) then
            Model%do25D = .true.
        else
            Model%do25D = .false.
        endif
        
        !Get mins/maxes
        if(present(xyzBdsIN)) then
            xyzBds(:) = xyzBdsIN(:)
        else
            call xmlInp%Set_Val(xyzBds(1),"idir/min",xMin)
            call xmlInp%Set_Val(xyzBds(2),"idir/max",xMax)
            call xmlInp%Set_Val(xyzBds(3),"jdir/min",xMin)
            call xmlInp%Set_Val(xyzBds(4),"jdir/max",xMax)
            call xmlInp%Set_Val(xyzBds(5),"kdir/min",xMin)
            call xmlInp%Set_Val(xyzBds(6),"kdir/max",xMax)
        end if

        if (Model%isLoud) then
            write(*,*) 'Grid generation ...'
            write(*,*) '   Cells = ', Grid%Nip,Grid%Njp,Grid%Nkp
            write(*,*) '   xMin/xMax = ', xyzBds(1),xyzBds(2)
            write(*,*) '   yMin/yMax = ', xyzBds(3),xyzBds(4)
            write(*,*) '   zMin/zMax = ', xyzBds(5),xyzBds(6)
            write(*,*) ''
        endif
        
        !Get grid geometry options
        !Use right-handed system
        !Cyl: x1,x2,x3 -> R,phi/2pi,z
        !Sph: x1,x2,x3 -> r,theta/pi,phi/2pi

        call xmlInp%Set_Val(doCyl,"grid/doCyl" ,.false.)
        call xmlInp%Set_Val(doSph,"grid/doSph" ,.false.)

        !Derived quantities
        Grid%Ni = Grid%Nip + 2*Model%nG
        Grid%Nj = Grid%Njp + 2*Model%nG
        Grid%Nk = Grid%Nkp + 2*Model%nG

        Grid%is = 1; Grid%ie = Grid%Nip
        Grid%js = 1; Grid%je = Grid%Njp
        Grid%ks = 1; Grid%ke = Grid%Nkp

        Grid%isg = Grid%is-Model%nG
        Grid%ieg = Grid%ie+Model%nG

        Grid%jsg = Grid%js-Model%nG
        Grid%jeg = Grid%je+Model%nG

        Grid%ksg = Grid%ks-Model%nG
        Grid%keg = Grid%ke+Model%nG

        Grid%ijkShift(IDIR) = Grid%Nip*Grid%Ri
        Grid%ijkShift(JDIR) = Grid%Njp*Grid%Rj
        Grid%ijkShift(KDIR) = Grid%Nkp*Grid%Rk

        !Set dx's (uniform for now)
        dx = (xyzBds(2)-xyzBds(1))/Grid%Nip
        dy = (xyzBds(4)-xyzBds(3))/Grid%Njp
        dz = (xyzBds(6)-xyzBds(5))/Grid%Nkp

        if( Grid%isTiled) then
            !Adjust grid info to only generate data for this tile
            Grid%Nip = Grid%Nip/Grid%NumRi
            Grid%Njp = Grid%Njp/Grid%NumRj
            Grid%Nkp = Grid%Nkp/Grid%NumRk

            Grid%ijkShift(IDIR) = Grid%Nip*Grid%Ri
            Grid%ijkShift(JDIR) = Grid%Njp*Grid%Rj
            Grid%ijkShift(KDIR) = Grid%Nkp*Grid%Rk

            Grid%Ni = Grid%Nip + 2*Model%nG
            Grid%Nj = Grid%Njp + 2*Model%nG
            Grid%Nk = Grid%Nkp + 2*Model%nG

            Grid%is = 1; Grid%ie = Grid%Nip
            Grid%js = 1; Grid%je = Grid%Njp
            Grid%ks = 1; Grid%ke = Grid%Nkp

            Grid%isg = Grid%is-Model%nG
            Grid%ieg = Grid%ie+Model%nG

            Grid%jsg = Grid%js-Model%nG
            Grid%jeg = Grid%je+Model%nG

            Grid%ksg = Grid%ks-Model%nG
            Grid%keg = Grid%ke+Model%nG
        endif

        !Allocate corner grid holders
        call allocGrid(Model,Grid)

        ! cell corners - a part of the Grid structure
        do k=Grid%ksg, Grid%keg+1
            do j=Grid%jsg, Grid%jeg+1
                do i=Grid%isg, Grid%ieg+1
                    x1 = xyzBds(1)+(Grid%ijkShift(IDIR)+i-1)*dx
                    x2 = xyzBds(3)+(Grid%ijkShift(JDIR)+j-1)*dy
                    x3 = xyzBds(5)+(Grid%ijkShift(KDIR)+k-1)*dz
                    if (doCyl) then
                        !Treat x1,x2,x3 as R,phi/2pi,z
                        x = x1*cos(2*pi*x2)
                        y = x1*sin(2*pi*x2)
                        z = x3
                    else if (doSph) then
                        !Sph: x1,x2,x3 -> r,theta/pi,phi/2pi
                        x = x1*sin(pi*x2)*cos(2*pi*x3)
                        y = x1*sin(pi*x2)*sin(2*pi*x3)
                        z = x1*cos(pi*x2)
                    else !Cartesian
                        x = x1
                        y = x2
                        z = x3
                    endif

                    Grid%xyz(i,j,k,XDIR:ZDIR) = [x,y,z]
                enddo
            enddo
        enddo

        !Do warp if necessary
        call xmlInp%Set_Val(doWarp,"grid/doWarp" ,.false.)
        if (doWarp) then
            call xmlInp%Set_Val(w0,"grid/w0",0.1_rp)
            call xmlInp%Set_Val(R ,"grid/R" ,1.0_rp)
            call xmlInp%Set_Val(Nw,"grid/Nw",1)

            Lx = (xyzBds(2)-xyzBds(1))
            Ly = (xyzBds(4)-xyzBds(3))
            Lz = (xyzBds(6)-xyzBds(5))
            Ax = Nw*pi/Lx
            Ay = Nw*pi/Ly
            Az = Nw*pi/Lz

            !Add warp factor to each corner
            do k=Grid%ksg, Grid%keg+1
                do j=Grid%jsg, Grid%jeg+1
                    do i=Grid%isg, Grid%ieg+1
                        
                        xp = Grid%xyz(i,j,k,XDIR) - xyzBds(1)
                        yp = Grid%xyz(i,j,k,YDIR) - xyzBds(3)
                        zp = Grid%xyz(i,j,k,ZDIR) - xyzBds(5)

                        dsp = w0*sin(Ax*xp)*sin(Ay*yp)

                        Grid%xyz(i,j,k,XDIR) = Grid%xyz(i,j,k,XDIR) + dsp
                        Grid%xyz(i,j,k,YDIR) = Grid%xyz(i,j,k,YDIR) - dsp
                        if (.not. Model%do25D) then
                            dsp = w0*sin(Az*zp)
                            Grid%xyz(i,j,k,ZDIR) = Grid%xyz(i,j,k,ZDIR) + dsp
                        endif
                    enddo
                enddo
            enddo

        endif

        call xmlInp%Set_Val(doStretch,"grid/doStretch" ,.false.)
        if (doStretch) then
            call xmlInp%Set_Val(xW,"grid/xW",1.0_rp)
            call xmlInp%Set_Val(yW,"grid/yW",1.0_rp)
            call xmlInp%Set_Val(zW,"grid/zW",1.0_rp)
            
            Lx = (xyzBds(2)-xyzBds(1))
            Ly = (xyzBds(4)-xyzBds(3))
            Lz = (xyzBds(6)-xyzBds(5))
            !Recalculate certain corners, for now just doing yp

            do k=Grid%ksg, Grid%keg+1
                do j=Grid%jsg, Grid%jeg+1
                    do i=Grid%isg, Grid%ieg+1
                        
                        !Remap to be eta = [0,1]
                        xp = ( Grid%xyz(i,j,k,XDIR) - xyzBds(1) )/Lx
                        yp = ( Grid%xyz(i,j,k,YDIR) - xyzBds(3) )/Ly
                        zp = ( Grid%xyz(i,j,k,ZDIR) - xyzBds(5) )/Lz

                        Grid%xyz(i,j,k,XDIR) = 2*Lx*xp + xW*(0.5*Lx - Lx*xp)*(1-xp)*xp - Lx
                        Grid%xyz(i,j,k,YDIR) = 2*Ly*yp + yW*(0.5*Ly - Ly*yp)*(1-yp)*yp - Ly
                        Grid%xyz(i,j,k,ZDIR) = 2*Lz*zp + zW*(0.5*Lz - Lz*zp)*(1-zp)*zp - Lz
                    enddo
                enddo
            enddo
        endif

    end subroutine genGridXML

    !Calculates grid data structure from cell corner arrays (x,y,z)
    subroutine Corners2Grid(Model,Grid)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid

        integer :: i,j,k,d
        integer :: T1,T2,dNorm
        integer :: is,ie,js,je,ks,ke

        real(rp), dimension(NDIM) :: xi,xj,xk,xip1,xjp1,xkp1, di,dj,dk
        real(rp), dimension(NDIM) :: dT1,dT2,eT1,eT2,eN
        real(rp), dimension(NDIM) :: dEdge, eAvg
        real(rp), dimension(NDIM) :: xFace
        real(rp), dimension(8,NDIM) :: xyzC
        real(rp) :: ijkDV(NVAR)

        real(rp) :: fA
        integer :: DelI,DelJ,DelK
        procedure(GasIC_T), pointer :: CdV
        procedure(VectorField_T), pointer :: fCtr,fArea !Functions for surface area/area center
        real(rp), dimension(NDIM) :: f0,f1,f2,f3, fInt,fInt2,fIntx
        logical :: doQuadFT = .false. !Whether to do quadrature for face system
        
        !Check for min length veclen
        !NOTE: Do we need <= here?
        if ( (Grid%ie-Grid%is+1) < VECLEN ) then
            write(*,*) 'Grid I-size smaller than VECLEN!'
            stop
        endif

        !Check for divisibility by double prec
        if ( mod(Grid%ie-Grid%is+1,ALIGN/8) /= 0 ) then
            write(*,*) 'Grid divisibility error, Grid I-size not divisible by ', ALIGN/8
            stop
        endif
        
        !Set pointers needed for Gaussian integrals
        CdV => CellDV
        fCtr => rVec
        fArea => IdVec

        !------------------------------------------------
        !Calculate face-centered coordinates
        allocate(Grid%xyzcc(Grid%isg:Grid%ieg  ,Grid%jsg:Grid%jeg  ,Grid%ksg:Grid%keg  ,NDIM))
        allocate(Grid%xfc  (Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,NDIM,NDIM))
        call allocGridVec(Model,Grid,Grid%face,doP1=.true.)
        Grid%face = 1.0
        Grid%xfc = 0.0

        allocate(Grid%Tf(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,NDIM*NDIM,NDIM))
        Grid%Tf = 0.0

        do d=1,NDIM
            DelI = ijkD(d,IDIR)
            DelJ = ijkD(d,JDIR)
            DelK = ijkD(d,KDIR)

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(f0,f1,f2,f3,fA,fInt,fInt2,fIntX,eN,eT1,eT2) &
            !$OMP private(dT1,dT2)
            do k=Grid%ksg, Grid%keg+DelK
                do j=Grid%jsg, Grid%jeg+DelJ
                    do i=Grid%isg, Grid%ieg+DelI
                        !Need specific ordering of corner points
                        call faceCoords(Model,Grid,i,j,k,d,f0,f1,f2,f3)

                        !Calculate face area/face center
                        call GaussianFaceIntegral(Model,f0,f1,f2,f3,fArea,fInt,fInt2,fIntX)
                        fA = fInt(1)
                        call GaussianFaceIntegral(Model,f0,f1,f2,f3,fCtr,fInt,fInt2,fIntX)
                        Grid%xfc(i,j,k,:,d) = fInt/fA
                        Grid%face(i,j,k,d) = fA

                        !Get face coordinate system
                        if (doQuadFT) then
                            call GaussianFaceSystem(f0,f1,f2,f3,eN,eT1,eT2)
                        else
                            select case(d)
                            case(IDIR)
                                !N,T1,T2 = i,j,k
                                dT1 = 0.5*( Grid%xyz(i  ,j+1,k+1,:) - Grid%xyz(i  ,j  ,k+1,:) ) + &
                                      0.5*( Grid%xyz(i  ,j+1,k  ,:) - Grid%xyz(i  ,j  ,k  ,:) )
                                dT2 = 0.5*( Grid%xyz(i  ,j+1,k+1,:) - Grid%xyz(i  ,j+1,k  ,:) ) + &
                                      0.5*( Grid%xyz(i  ,j  ,k+1,:) - Grid%xyz(i  ,j  ,k  ,:) )
                            case(JDIR)
                                !N,T1,T2 = j,k,i
                                dT1 = 0.5*( Grid%xyz(i+1,j  ,k+1,:) - Grid%xyz(i+1,j  ,k  ,:) ) + &
                                      0.5*( Grid%xyz(i  ,j  ,k+1,:) - Grid%xyz(i  ,j  ,k  ,:) )
                                dT2 = 0.5*( Grid%xyz(i+1,j  ,k+1,:) - Grid%xyz(i  ,j  ,k+1,:) ) + &
                                      0.5*( Grid%xyz(i+1,j  ,k  ,:) - Grid%xyz(i  ,j  ,k  ,:) )
                            case(KDIR)
                                !N,T1,T2 = k,i,j
                                dT1 = 0.5*( Grid%xyz(i+1,j+1,k  ,:) - Grid%xyz(i  ,j+1,k  ,:) ) + &
                                      0.5*( Grid%xyz(i+1,j  ,k  ,:) - Grid%xyz(i  ,j  ,k  ,:) )
                                dT2 = 0.5*( Grid%xyz(i+1,j+1,k  ,:) - Grid%xyz(i+1,j  ,k  ,:) ) + &
                                      0.5*( Grid%xyz(i  ,j+1,k  ,:) - Grid%xyz(i  ,j  ,k  ,:) )
                            end select
                            !Face normal
                            eN = cross(dT1,dT2)/norm2(cross(dT1,dT2))
                            eT2 = dT2/norm2(dT2)
                            eT1 = cross(eT2,eN)/norm2(cross(eT2,eN))

                        endif
                        !Use whichever system you calculated
                        Grid%Tf(i,j,k,NORMX:NORMZ,d) = eN
                        Grid%Tf(i,j,k,TAN1X:TAN1Z,d) = eT1
                        Grid%Tf(i,j,k,TAN2X:TAN2Z,d) = eT2

                    enddo
                enddo
            enddo                        
        enddo


        !------------------------------------------------
        !Calculate cell volume, grid edges
        call allocGridVar(Model,Grid,Grid%volume)
        call allocGridVar(Model,Grid,Grid%di)
        call allocGridVar(Model,Grid,Grid%dj)
        call allocGridVar(Model,Grid,Grid%dk)

   
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyzC,di,dj,dk,ijkDV)
        do k=Grid%ksg, Grid%keg
            do j=Grid%jsg, Grid%jeg
                do i=Grid%isg, Grid%ieg
                    !Get di,dj,dk vectors across cell
                    di = Grid%xfc(i+1,j,k,:,IDIR) - Grid%xfc(i,j,k,:,IDIR)
                    dj = Grid%xfc(i,j+1,k,:,JDIR) - Grid%xfc(i,j,k,:,JDIR)
                    dk = Grid%xfc(i,j,k+1,:,KDIR) - Grid%xfc(i,j,k,:,KDIR)

                    !Calculate average cell length in each direction for timestep
                    Grid%di(i,j,k) = norm2(di)
                    Grid%dj(i,j,k) = norm2(dj)
                    Grid%dk(i,j,k) = norm2(dk)

                !Calculate volume/volume-center using quadrature
                    !Get cell corners
                    call cellCoords(Model,Grid,i,j,k,xyzC)

                    ijkDV = GaussianVolumeIntegral(xyzC,CdV)                                        
                    Grid%volume(i,j,k)   = ijkDV(1)
                    Grid%xyzcc (i,j,k,:) = ijkDV(VELX:VELZ)/ijkDV(1)

                    !Grid%volume(i,j,k) = dot_product(di,cross(dj,dk))
                enddo
            enddo
        enddo

        !Replace di,dj,dk calculations
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(fA)
        do k=Grid%ksg, Grid%keg
            do j=Grid%jsg, Grid%jeg
                do i=Grid%isg, Grid%ieg
                    fA = max(Grid%face(i+1,j,k,IDIR),Grid%face(i,j,k,IDIR))
                    Grid%di(i,j,k) = Grid%volume(i,j,k)/fA
                    fA = max(Grid%face(i,j+1,k,JDIR),Grid%face(i,j,k,JDIR))
                    Grid%dj(i,j,k) = Grid%volume(i,j,k)/fA
                    fA = max(Grid%face(i,j,k+1,KDIR),Grid%face(i,j,k,KDIR))
                    Grid%dk(i,j,k) = Grid%volume(i,j,k)/fA

                enddo
            enddo
        enddo

        !Fix transforms if necessary
        if (Model%doRing) call RingGridFix(Model,Grid)

        !------------------------------------------------
        !Calculate coordinate systems at edges for magnetic field updates (velocity)
        call allocGridVec(Model,Grid,Grid%edge,doP1=.true.)
        Grid%edge = 1.0
        
        allocate(Grid%Te(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,NDIM*NDIM,NDIM))
        Grid%Te = 0.0

        do d=1,NDIM
            !Use maximal bounds
            is = Grid%isg
            ie = Grid%ieg+1
            js = Grid%jsg
            je = Grid%jeg+1
            ks = Grid%ksg
            ke = Grid%keg+1

            !Correct bounds for direction (go maximal)
            select case(d)
            case(IDIR)
                ie = Grid%ieg
                ks = Grid%ksg+1
                ke = Grid%keg
            case(JDIR)
                je = Grid%jeg
                is = Grid%isg+1
                ie = Grid%ieg
            case(KDIR)
                ke = Grid%keg
                js = Grid%jsg+1
                je = Grid%jeg
            end select

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,dEdge,eAvg,eN,eT1,eT2)
            do k = ks, ke
                do j = js, je
                    do i = is, ie

                        !Get edge vector
                        !Get average (k,i,j) vector for (i,j,k) edge center (needed for T1 direction)
                        select case(d)

                        case (IDIR)
                            dEdge = Grid%xyz(i+1,j,k,:) - Grid%xyz(i,j,k,:)
                            eAvg = 0.5*( Grid%xyz(i  ,j  ,k+1,:) - Grid%xyz(i  ,j  ,k-1,:) &
                                       + Grid%xyz(i+1,j  ,k+1,:) - Grid%xyz(i+1,j  ,k-1,:) )

                        case (JDIR)
                            dEdge = Grid%xyz(i,j+1,k,:) - Grid%xyz(i,j,k,:)
                            eAvg = 0.5*( Grid%xyz(i+1,j  ,k  ,:) - Grid%xyz(i-1,j  ,k  ,:) &
                                       + Grid%xyz(i+1,j+1,k  ,:) - Grid%xyz(i-1,j+1,k  ,:) )

                        case (KDIR)
                            dEdge = Grid%xyz(i,j,k+1,:) - Grid%xyz(i,j,k,:)
                            eAvg = 0.5*( Grid%xyz(i  ,j+1,k  ,:) - Grid%xyz(i  ,j-1,k  ,:) &
                                       + Grid%xyz(i  ,j+1,k+1,:) - Grid%xyz(i  ,j-1,k+1,:) )

                        end select

                        !Edge vector
                        eN = dEdge/norm2(dEdge)

                        !Finish triad calculation
                        !T1 = eAvg x eN
                        !T2 = eNormal x eT1, right-handed
                        eT1 = cross(eAvg,eN)/norm2(cross(eAvg,eN))
                        eT2 = cross(eN,eT1)/norm2(cross(eN,eT1))

                        !Save edge lengths and transform
                        Grid%edge(i,j,k,d) = norm2(dEdge)
                        Grid%Te(i,j,k,NORMX:NORMZ,d) = eN
                        Grid%Te(i,j,k,TAN1X:TAN1Z,d) = eT1
                        Grid%Te(i,j,k,TAN2X:TAN2Z,d) = eT2

                    enddo
                enddo
            enddo !k loop
        enddo !Dim loop

        
        !------------------------------------------------
        !Calculate coordinate systems at edges for magnetic field updates (magnetic field)
        allocate(Grid%Teb(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,4,NDIM))
        Grid%Teb = 0.0
        !Loop over normal direction, and get edge system for plane w/ that normal
        do dNorm=1,NDIM
            !Get local triad
            select case(dNorm)
                !TODO Replace this with Levi-Cevita?
                case(IDIR)
                    T1 = JDIR; T2 = KDIR
                case(JDIR)
                    T1 = KDIR; T2 = IDIR
                case(KDIR)
                    T1 = IDIR; T2 = JDIR
            end select  
            !xnqi,ynqi <- (dNorm,dT1,dT2)
            !xnqj,ynqj <- (dNorm,dT2,dT1)

            call ebGeom(Model,Grid,Grid%Teb(:,:,:,XNQI:YNQI,dNorm),dNorm,T1,T2)
            call ebGeom(Model,Grid,Grid%Teb(:,:,:,XNQJ:YNQJ,dNorm),dNorm,T2,T1)

        enddo    

        if (Model%doSource) then
            allocate(Grid%Gas0(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NVAR,0:Model%nSpc))
            Grid%Gas0 = 0.0
        endif
        
    end subroutine Corners2Grid

    !Calculate one of the two normal vectors of the edge-centered mag field system
    subroutine ebGeom(Model,Gr,nQ,dNorm,dT1,dT2)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(inout) :: nQ(Gr%isg:Gr%ieg+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,2)
        integer, intent(in) :: dNorm,dT1,dT2
        
        !Reconstruction stencils for Nx,Ny,Nz and face area
        real(rp), dimension(vecLen,recLen,NDIM) :: NxyzB
        real(rp), dimension(vecLen,recLen) :: fArB
        integer :: i,j,k,ie
        integer :: iG,iB,iMax
        real(rp) :: dxT,dyT,dScl,fAi
        real(rp) :: Nec(NDIM),Nij(2)

        !Initialize all values (array is larger than necessary)
        nQ = 1.0
        ie = Gr%ie+1

        !K: Testing restricted interpolation (6Cent)
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,iB,j,k,iMax,iG) &
        !$OMP private(fAi,Nec,Nij,dxT,dyT,dScl,NxyzB,fArB)
        do k=Gr%ks,Gr%ke+1
            do j=Gr%js,Gr%je+1
                do iB=Gr%is,ie,vecLen
                    iMax = min(vecLen,ie-iB+1)

                    !Get stencils in the dT2 direction
                    call recLoadBlockI(Model,Gr,fArB           ,Gr%face(:,:,:,dT1)    ,iB,j,k,iMax,dT2)
                    call recLoadBlockI(Model,Gr,NxyzB(:,:,XDIR),Gr%Tf(:,:,:,NORMX,dT1),iB,j,k,iMax,dT2)
                    call recLoadBlockI(Model,Gr,NxyzB(:,:,YDIR),Gr%Tf(:,:,:,NORMY,dT1),iB,j,k,iMax,dT2)
                    call recLoadBlockI(Model,Gr,NxyzB(:,:,ZDIR),Gr%Tf(:,:,:,NORMZ,dT1),iB,j,k,iMax,dT2)

                    !Do weighted interpolation from face to corner of each component
                    do i=1,iMax
                        iG = iB+i-1 !Global index
                        fAi = dot_product(interpWgt6,fArB(i,:))
                        Nec(XDIR) = dot_product(NxyzB(i,:,XDIR)*fArB(i,:),interpWgt6)/fAi
                        Nec(YDIR) = dot_product(NxyzB(i,:,YDIR)*fArB(i,:),interpWgt6)/fAi
                        Nec(ZDIR) = dot_product(NxyzB(i,:,ZDIR)*fArB(i,:),interpWgt6)/fAi

                        !Map xyz->1,2 with edge mapping (dNorm)
                        Nij(1) = dot_product(Nec,Gr%Te(iG,j,k,TAN1X:TAN1Z,dNorm))
                        Nij(2) = dot_product(Nec,Gr%Te(iG,j,k,TAN2X:TAN2Z,dNorm))

                        dxT = Nij(1)
                        dyT = Nij(2)
                        dScl = sqrt(dxT**2.0+dyT**2.0)
                        
                        if (dScl <= TINY) then
                            dScl = 1.0
                            write(*,*) 'dScl < TINY, ijk = ', i,j,k
                        endif
                        nQ(iG,j,k,1) = dxT/dScl
                        nQ(iG,j,k,2) = dyT/dScl                        

                    enddo

                enddo !iB
            enddo
        enddo !K loop

    end subroutine ebGeom

end module init
