!Various routines to read/write HDF5 files

module gioH5
    use gamtypes
    use gamutils
    use gridutils
    use ioH5
    use multifluid
    use planethelper
    use dates
    use files
    
    implicit none

    integer, parameter, private :: MAXIOVAR = 50
    type(IOVAR_T), dimension(MAXIOVAR), private :: IOVars
    logical, private :: doRoot = .true. !Whether root variables need to be written
    logical, private :: doFat = .false. !Whether to output lots of extra data
    logical, private :: doHighPrec = .false. !Whether to output 64bit standard output

    !Necessary for IO routines
    character(len=strLen) ,public:: GamH5File

    contains

    subroutine SetFatIO()

        doFat = .true.
    end subroutine SetFatIO
    
    subroutine SetHPOut()
        doHighPrec = .true.
    end subroutine SetHPOut

    subroutine readH5Grid(Model,Grid,inH5)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Grid
        character(len=*), intent(in) :: inH5

        integer, dimension(NDIM) :: dims        

        !Reset IO chain
        call ClearIO(IOVars)

        call CheckFileOrDie(inH5,"Unable to open input mesh, exiting ...")

        dims = GridSizeH5(inH5)

        !Start with main indices, convert # of ghost corners to active
        Grid%Nip = dims(1) - 2*Model%nG - 1
        Grid%Njp = dims(2) - 2*Model%nG - 1
        Grid%Nkp = dims(3) - 2*Model%nG - 1
    
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

        if( Model%isRestart .or. .not. Grid%isTiled ) then
            !Read the entire grid

            !Setup input chain
            call AddInVar(IOVars,"X")
            call AddInVar(IOVars,"Y")
            call AddInVar(IOVars,"Z")
        else
            !Read only the part of the grid for this tile

            !Adjust the grid values for this tile
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

            !Setup input chain
            call AddInVarHyper(IOVars,"X",Grid%ijkShift,(/ Grid%Ni+1,Grid%Nj+1,Grid%Nk+1 /))
            call AddInVarHyper(IOVars,"Y",Grid%ijkShift,(/ Grid%Ni+1,Grid%Nj+1,Grid%Nk+1 /))
            call AddInVarHyper(IOVars,"Z",Grid%ijkShift,(/ Grid%Ni+1,Grid%Nj+1,Grid%Nk+1 /))
        endif

        call ReadVars(IOVars,.false.,inH5) !Don't use io precision

        call allocGrid(Model,Grid)
        
        Grid%xyz(:,:,:,XDIR) = reshape(IOVars(XDIR)%data,[Grid%Ni+1,Grid%Nj+1,Grid%Nk+1])
        Grid%xyz(:,:,:,YDIR) = reshape(IOVars(YDIR)%data,[Grid%Ni+1,Grid%Nj+1,Grid%Nk+1])
        Grid%xyz(:,:,:,ZDIR) = reshape(IOVars(ZDIR)%data,[Grid%Ni+1,Grid%Nj+1,Grid%Nk+1])
    
    end subroutine readH5Grid

    !Write initial grid info to root of H5 output file
    subroutine writeH5GridInit(Model,Gr)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr

        real(rp), dimension(:,:,:)  , allocatable :: gQ !Grid quality
        real(rp), dimension(:,:,:,:), allocatable :: gVec
        logical :: isExist

        integer :: iMin,iMax,jMin,jMax,kMin,kMax
        integer :: i,j,k
        character(len=strLen) :: vID

        !Don't call this function again
        doRoot = .false.

        !Test if root variables (grid/force info is already there)
        vID = "X" !Value to test for
        isExist = ioExist(GamH5File,vID)

        if (isExist) then
            return
        endif
        
        !Calculate grid quality
        call allocGridVar(Model,Gr,gQ)
        call GridQuality(Model,Gr,gQ) 

        !Reset IO chain
        call ClearIO(IOVars)

        if(writeGhosts) then
            ! for debugging
            iMin = Gr%isg
            iMax = Gr%ieg
            jMin = Gr%jsg
            jMax = Gr%jeg
            kMin = Gr%ksg
            kMax = Gr%keg
        else
            iMin = Gr%is
            iMax = Gr%ie
            jMin = Gr%js
            jMax = Gr%je
            kMin = Gr%ks
            kMax = Gr%ke
        endif

        !Fill IO chain, start with coordinates
        if (Gr%Nkp > 1) then
            !3D problem
            call AddOutVar(IOVars,"X",Gr%xyz(iMin:iMax+1,jMin:jMax+1,kMin:kMax+1,XDIR))
            call AddOutVar(IOVars,"Y",Gr%xyz(iMin:iMax+1,jMin:jMax+1,kMin:kMax+1,YDIR))
            call AddOutVar(IOVars,"Z",Gr%xyz(iMin:iMax+1,jMin:jMax+1,kMin:kMax+1,ZDIR))
        else
            !2D problem
            !Squash corner arrays to 2D
            call AddOutVar(IOVars,"X",reshape(Gr%xyz(iMin:iMax+1,jMin:jMax+1,kMin:kMin,XDIR),[iMax-iMin+2,jMax-jMin+2]))
            call AddOutVar(IOVars,"Y",reshape(Gr%xyz(iMin:iMax+1,jMin:jMax+1,kMin:kMin,YDIR),[iMax-iMin+2,jMax-jMin+2]))
        endif

        call AddOutVar(IOVars,"dV",Gr%volume(iMin:iMax,jMin:jMax,kMin:kMax))
        if (doFat) then
            call AddOutVar(IOVars,"gQ",       gQ(iMin:iMax,jMin:jMax,kMin:kMax))
        endif
        if (Model%doMHD .and. Model%doBackground) then
            !Write out background field and force density
            associate(gamOut=>Model%gamOut)
            call GameraOut("Bx0",gamOut%bID,gamOut%bScl,Gr%B0(iMin:iMax,jMin:jMax,kMin:kMax,XDIR))
            call GameraOut("By0",gamOut%bID,gamOut%bScl,Gr%B0(iMin:iMax,jMin:jMax,kMin:kMax,YDIR))
            call GameraOut("Bz0",gamOut%bID,gamOut%bScl,Gr%B0(iMin:iMax,jMin:jMax,kMin:kMax,ZDIR))
            end associate

            if (doFat) then
                call AddOutVar(IOVars,"dPxB0",Gr%dpB0(iMin:iMax,jMin:jMax,kMin:kMax,XDIR))
                call AddOutVar(IOVars,"dPyB0",Gr%dpB0(iMin:iMax,jMin:jMax,kMin:kMax,YDIR))
                call AddOutVar(IOVars,"dPzB0",Gr%dpB0(iMin:iMax,jMin:jMax,kMin:kMax,ZDIR))
            endif
        endif
        if (Model%doGrav .and. doFat) then
            !Write out grav accelerations
            call AddOutVar(IOVars,"gx",Gr%gxyz(iMin:iMax,jMin:jMax,kMin:kMax,XDIR))
            call AddOutVar(IOVars,"gy",Gr%gxyz(iMin:iMax,jMin:jMax,kMin:kMax,YDIR))
            call AddOutVar(IOVars,"gz",Gr%gxyz(iMin:iMax,jMin:jMax,kMin:kMax,ZDIR))

        endif

        if (Model%doMHD .and. Model%isMagsphere) then
            !Add mag moment
            call AddOutVar(IOVars,"MagM0",Model%MagM0*Model%gamOut%bScl)
            
            !Write out dipole field values
            allocate(gVec (iMin:iMax,jMin:jMax,kMin:kMax,1:NDIM))
            !Subtract dipole before calculating current
            !$OMP PARALLEL DO default(shared) collapse(2)
            do k=kMin,kMax
                do j=jMin,jMax
                    do i=iMin,iMax
                        gVec(i,j,k,:) = MagsphereDipole(Gr%xyzcc(i,j,k,:),Model%MagM0)
                    enddo
                enddo
            enddo
            associate(gamOut=>Model%gamOut)
            call GameraOut("BxD",gamOut%bID,gamOut%bScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,XDIR))
            call GameraOut("ByD",gamOut%bID,gamOut%bScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,YDIR))
            call GameraOut("BzD",gamOut%bID,gamOut%bScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,ZDIR))
            end associate
        endif

        !Add information about time scaling/units
        call AddOutVar(IOVars,"tScl"   ,Model%gamOut%tScl)
        call AddOutVar(IOVars,"timeID" ,Model%gamOut%tID)
        call AddOutVar(IOVars,"UnitsID",Model%gamOut%uID)

        !Call user routine if defined
        if (associated(Model%HackIO_0)) then
            call Model%HackIO_0(Model,Gr,IOVars)
        endif

        !Write out the chain
        if(doHighPrec) then
            call WriteVars(IOVars,.false.,GamH5File)
        else
            call WriteVars(IOVars,.true.,GamH5File)
        endif
        
    end subroutine writeH5GridInit

    !Write specific output slice
    subroutine writeSlc(Model,Gr,State,gStr)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        type(State_T), intent(in) :: State
        character(len=*), intent(in) :: gStr

        integer :: i,j,k,s
        character(len=strLen) :: dID,VxID,VyID,VzID,PID
        integer iMin,iMax,jMin,jMax,kMin,kMax

        !Fill base data
        real(rp), dimension(:,:,:),   allocatable :: gVar,DivBcc,gVar1
        real(rp), dimension(:,:,:,:), allocatable :: gVec
        real(rp), dimension(:,:,:,:), allocatable :: VecA,VecB !Full-sized arrays
        real(rp), dimension(:,:,:)  , allocatable :: dVacMask
        real(rp), dimension(:)      , allocatable :: dVacFloors

        real(rp) :: totDivB,MJD

        !Check if root variables need to be written
        if (doRoot) then
            call writeH5GridInit(Model,Gr)
            doRoot = .false.
        endif
        !Reset IO chain
        call ClearIO(IOVars)

        if(writeGhosts) then
            ! for debugging
            call AddOutVar(IOVars,'hasGhosts',1)
            iMin = Gr%isg
            iMax = Gr%ieg
            jMin = Gr%jsg
            jMax = Gr%jeg
            kMin = Gr%ksg
            kMax = Gr%keg
        else
            call AddOutVar(IOVars,'hasGhosts',0)
            iMin = Gr%is
            iMax = Gr%ie
            jMin = Gr%js
            jMax = Gr%je
            kMin = Gr%ks
            kMax = Gr%ke
        endif

        !Allocate holders
        allocate(gVar (iMin:iMax,jMin:jMax,kMin:kMax))
        allocate(gVar1(iMin:iMax,jMin:jMax,kMin:kMax))
        allocate(gVec (iMin:iMax,jMin:jMax,kMin:kMax,1:NDIM))
        allocate(dVacMask(iMin:iMax,jMin:jMax,kMin:kMax))
        allocate(dVacFloors(0:Model%nSpc))

        associate(gamOut=>Model%gamOut)

        dVacFloors(BLK) = 0.0
        if (Model%doMultiF) then
            dVacFloors(1:Model%nSpc) = Spcs(:)%dVac+TINY !For zeroing out mass output in degenerate fluids
        endif

        do s=0,Model%nSpc
            if (s == BLK) then
                dID = "D"
                VxID = "Vx"
                VyID = "Vy"
                VzID = "Vz"
                PID = "P"
            else
                write(dID ,'(A,I0)') "D" , s
                write(VxID,'(A,I0)') "Vx", s
                write(VyID,'(A,I0)') "Vy", s
                write(VzID,'(A,I0)') "Vz", s
                write(PID ,'(A,I0)') "P" , s
            endif

            !Set species density to zero if it's below threshold for that fluid
            do i=iMin,iMax
                do j=jMin,jMax
                    do k=kMin,kMax
                        if (State%Gas(i,j,k,DEN,s) < dVacFloors(s)) then
                            dVacMask(i,j,k) = 0.0
                        else
                            dVacMask(i,j,k) = State%Gas(i,j,k,DEN,s)
                        endif
                    enddo
               enddo
            enddo

            !Density
            call GameraOut(dID,gamOut%dID,gamOut%dScl,dVacMask)

            !---------------------
            !Calculate Velocities/Pressure
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k)
            do k=kMin,kMax
                do j=jMin,jMax
                    do i=iMin,iMax
                        if (State%Gas(i,j,k,DEN,s)>TINY) then
                            gVec(i,j,k,:) = State%Gas(i,j,k,MOMX:MOMZ,s)/State%Gas(i,j,k,DEN,s)
                            gVar(i,j,k) = (Model%gamma-1)*( State%Gas(i,j,k,ENERGY,s) - &
                                          0.5*(State%Gas(i,j,k,MOMX,s)**2 + &
                                               State%Gas(i,j,k,MOMY,s)**2 + &
                                               State%Gas(i,j,k,MOMZ,s)**2)/ &
                                               State%Gas(i,j,k,DEN ,s))
                        else
                            gVec(i,j,k,:) = 0.0
                            gVar(i,j,k)   = 0.0
                        endif
                        if (s == BLK) then
                            call CellPress2Cs(Model,State%Gas(i,j,k,:,s),gVar1(i,j,k))
                        endif

                    enddo
                enddo
            enddo

            !Add V/P to chain
            call GameraOut(VxID,gamOut%vID,gamOut%vScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,XDIR))
            call GameraOut(VyID,gamOut%vID,gamOut%vScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,YDIR))
            call GameraOut(VzID,gamOut%vID,gamOut%vScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,ZDIR))
            call GameraOut(PID ,gamOut%pID,gamOut%pScl,gVar(iMin:iMax,jMin:jMax,kMin:kMax))
            if (s == BLK) then
                call GameraOut("Cs",gamOut%vID,gamOut%vScl,gVar1(iMin:iMax,jMin:jMax,kMin:kMax))
            endif
        enddo !Species loop

        !---------------------
        !Write MHD variables
        if (Model%doMHD) then
            call allocGridVec(Model,Gr,VecA)
            call allocGridVec(Model,Gr,VecB)

            !For current, use VecA to hold total Bxyz, VecB for Jxyz
            if (Model%doBackground) then
                VecA = State%Bxyz + Gr%B0 !Full size
                gVec(:,:,:,:) = Gr%B0(iMin:iMax,jMin:jMax,kMin:kMax,XDIR:ZDIR) + State%Bxyz(iMin:iMax,jMin:jMax,kMin:kMax,XDIR:ZDIR)
            else
                gVec(:,:,:,:) = State%Bxyz(iMin:iMax,jMin:jMax,kMin:kMax,XDIR:ZDIR)
                VecA = State%Bxyz
            endif

            ! If writing magnetic flux for debugging
            if (writeMagFlux) then
                call AddOutVar(IOVars,'MFi',State%magFlux(iMin:iMax+1,jMin:jMax,kMin:kMax,IDIR))
                call AddOutVar(IOVars,'MFj',State%magFlux(iMin:iMax,jMin:jMax+1,kMin:kMax,JDIR))
                call AddOutVar(IOVars,'MFk',State%magFlux(iMin:iMax,jMin:jMax,kMin:kMax+1,KDIR))
            endif

            !Add magnetic fields
            call GameraOut("Bx",gamOut%bID,gamOut%bScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,XDIR))
            call GameraOut("By",gamOut%bID,gamOut%bScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,YDIR))
            call GameraOut("Bz",gamOut%bID,gamOut%bScl,gVec(iMin:iMax,jMin:jMax,kMin:kMax,ZDIR))

            !Add mag pressure
            gVar = 0.5*(gVec(:,:,:,XDIR)**2.0 + gVec(:,:,:,YDIR)**2.0 + gVec(:,:,:,ZDIR)**2.0)
            call GameraOut("Pb",gamOut%pID,gamOut%pScl,gVar(iMin:iMax,jMin:jMax,kMin:kMax))

            !Write current
            if (Model%isMagsphere) then
                !Subtract dipole before calculating current
                !$OMP PARALLEL DO default(shared) collapse(2)
                do k=Gr%ksg,Gr%keg
                    do j=Gr%jsg,Gr%jeg
                        do i=Gr%isg,Gr%ieg
                            VecA(i,j,k,:) = State%Bxyz(i,j,k,:) + Gr%B0(i,j,k,:) - MagsphereDipole(Gr%xyzcc(i,j,k,:),Model%MagM0)
                        enddo
                    enddo
                enddo

                call bFld2Jxyz(Model,Gr,VecA,VecB)
                gVec(:,:,:,:) = VecB(iMin:iMax,jMin:jMax,kMin:kMax,XDIR:ZDIR)
            else
                !Do full current
                VecA = State%Bxyz
                call bFld2Jxyz(Model,Gr,VecA,VecB)
                gVec(:,:,:,:) = VecB(iMin:iMax,jMin:jMax,kMin:kMax,XDIR:ZDIR)
            endif

            !Output current density
            call GameraOut("Jx",gamOut%jID,gamOut%jScl,gVec(:,:,:,XDIR))
            call GameraOut("Jy",gamOut%jID,gamOut%jScl,gVec(:,:,:,YDIR))
            call GameraOut("Jz",gamOut%jID,gamOut%jScl,gVec(:,:,:,ZDIR))

            !Calculate/Write xyz electric fields
            if (doFat) then
            
                !Divide by edge-length to go from potential to field
                !$OMP PARALLEL DO default(shared) collapse(2)
                do k=Gr%ksg,Gr%keg
                    do j=Gr%jsg,Gr%jeg
                        do i=Gr%isg,Gr%ieg
                            if(any(Gr%Edge(i,j,k,:) == 0)) then
                                 VecA(i,j,k,:) = 0.0
                            else
                                 VecA(i,j,k,:) = State%Efld(i,j,k,:)/Gr%Edge(i,j,k,:)
                            end if
                        enddo
                    enddo
                enddo
                call Eijk2xyz(Model,Gr,VecA,VecB)
                gVec(:,:,:,:) = VecB(iMin:iMax,jMin:jMax,kMin:kMax,XDIR:ZDIR)
                call FixRAVec(gVec(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,1:NDIM))
            
                call GameraOut("Ex",gamOut%eID,gamOut%eScl,gVec(:,:,:,XDIR))
                call GameraOut("Ey",gamOut%eID,gamOut%eScl,gVec(:,:,:,YDIR))
                call GameraOut("Ez",gamOut%eID,gamOut%eScl,gVec(:,:,:,ZDIR))

            endif

            if (Model%doSource) then
                call GameraOut("SrcD" ,gamOut%dID,gamOut%dScl,Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,DEN     ,BLK))
                call GameraOut("SrcP" ,gamOut%pID,gamOut%pScl,Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,PRESSURE,BLK))
                if (Model%isMagsphere) then
                    call GameraOut("SrcX1","DEG",1.0_rp     ,Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,VELX    ,BLK))
                    call GameraOut("SrcX2","DEG",1.0_rp     ,Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,VELY    ,BLK))
                    call GameraOut("SrcDT","s"  ,gamOut%tScl,Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,VELZ    ,BLK))
                else
                    call GameraOut("SrcVx","CODE"    ,1.0_rp     ,Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,VELX    ,BLK))
                    call GameraOut("SrcVy","CODE"    ,1.0_rp     ,Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,VELY    ,BLK))
                    call GameraOut("SrcVz","CODE"    ,1.0_rp     ,Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,VELZ    ,BLK))
                endif                    
            endif

            if(Model%doResistive) then
                !$OMP PARALLEL DO default(shared) collapse(2) &
                !$OMP private(i,j,k)
                do k=kMin,kMax
                    do j=jMin,jMax
                        do i=iMin,iMax
                            !Save cell-centered eta
                            gVar(i,j,k) = EdgeScalar2CC(Model,Gr,State%Deta,i,j,k)

                            !Save cell-centered diffusive velocity
                            gVar1(i,j,k) = gVar(i,j,k)*2.0/minval([Gr%di(i,j,k),Gr%dj(i,j,k),Gr%dk(i,j,k)])
                        enddo
                    enddo
                enddo

                !Should change this to have more meaningful scaling
                call GameraOut("Eta","CODE",1.0_rp,gVar(iMin:iMax,jMin:jMax,kMin:kMax))
                !Output diffusive velocity scaled to proper output velocity units
                call GameraOut("Vdiff",gamOut%vID,gamOut%vScl,gVar1(iMin:iMax,jMin:jMax,kMin:kMax))
            end if

            !Write divergence if necessary
            if (Model%doDivB .and. doFat) then
                call allocGridVar(Model,Gr,DivBcc)
                call DivB(Model,Gr,State,totDivB,DivBcc,doTotO=.true.)
                gVar = DivBcc(iMin:iMax,jMin:jMax,kMin:kMax)
                call AddOutVar(IOVars,"DivB",gVar)

                if (Model%doBackground) then
                    !Also calculate divergence of perturbation field
                    call DivB(Model,Gr,State,totDivB,DivBcc,doTotO=.false.)
                    gVar = DivBcc(iMin:iMax,jMin:jMax,kMin:kMax)
                    call AddOutVar(IOVars,"DivdB",gVar)
                endif !B0
                
                deallocate(DivBcc)
            endif
            deallocate(VecA,VecB)
        endif
        deallocate(gVec,gVar)

        !---------------------
        !Do attributes

        call AddOutVar(IOVars,"time",gamOut%tScl*Model%t)
        call AddOutVar(IOVars,"timestep",Model%ts)
        call AddOutVar(IOVars,"dt",gamOut%tScl*Model%dt)
        if ( Model%MJD0 >= (-TINY) ) then
            MJD = T2MJD(Model%t*Model%Units%gT0,Model%MJD0)
            call AddOutVar(IOVars,"MJD",MJD)
        endif
        call AddOutVar(IOVars,"nSpc",Model%nSpc)
        call AddOutVar(IOVars,"kzcsMHD",Model%kzcsMHD,uStr="kZCs",dStr="MHD-only kZCs")
        call AddOutVar(IOVars,"kzcsTOT",Model%kzcsTOT,uStr="kZCs",dStr="Total kZCs"   )

        !---------------------

        !Call user routine if defined
        if (associated(Model%HackIO)) then
            call Model%HackIO(Model,Gr,State,IOVars)
        endif

        !------------------
        !Finalize
        if(doHighPrec) then
            call WriteVars(IOVars,.false.,GamH5File,trim(gStr))
        else
            call WriteVars(IOVars,.true.,GamH5File,trim(gStr))
        endif

        end associate

        contains

            !Fix up cell-centered vector (like current or electric field) to deal with axis
            !Note, this is only for output purposes since we don't have proper ghost information

            !Estimate Qxyz at pole by averaging about second ring
            !Calculate around first ring to interpolate pole/ring-2 values
            subroutine FixRAVec(Qxyz)
                real(rp), intent(inout) :: Qxyz(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,1:NDIM)

                real(rp) :: Q0(NDIM)
                integer :: i,k
                real(rp) :: w1,w2
                !Check if we're doing ring avg and have either pole
                if ( Model%doRing .and. (Model%Ring%doS .or. Model%Ring%doE) ) then

                    select case (Model%Ring%GridID)
                        !------------------
                        case ("lfm")
                            !Move along axis
                            do i=Gr%is,Gr%ie
                                if (Model%Ring%doS) then
                                    call FixRAVec_S(Model,Gr,Qxyz(i,Gr%js:Gr%js+1,Gr%ks:Gr%ke,1:NDIM))
                                endif
                                if (Model%Ring%doE) then
                                    call FixRAVec_E(Model,Gr,Qxyz(i,Gr%je-1:Gr%je,Gr%ks:Gr%ke,1:NDIM))
                                endif
                            enddo !I-loop
                    end select

                endif
            end subroutine FixRAVec
    end subroutine writeSlc

    subroutine GameraOut(vID,uID,vScl,V)
        character(len=*), intent(in) :: vID,uID
        real(rp), intent(in) :: vScl
        real(rp), intent(in) :: V(:,:,:)

        integer :: n0
        call AddOutVar(IOVars,vID,V)
        n0 = FindIO(IOVars,vID)
        IOVars(n0)%scale = vScl
        IOVars(n0)%unitStr = uID
    end subroutine GameraOut

    subroutine GridQuality(Model,Gr,gQ)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        real(rp), intent(inout) :: gQ(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg)

        integer :: i,j,k
        real(rp), dimension(NDIM) :: dXcc, ddV,ddx,gdV

        gQ = 0.0
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,dXcc,ddV,ddx,gdV)
        do k=Gr%ks, Gr%ke
            do j=Gr%js, Gr%je
                do i=Gr%is, Gr%ie
                    dXcc = [Gr%di(i,j,k),Gr%dj(i,j,k),Gr%dk(i,j,k)]
                    
                    ddV = [Gr%volume(i+1,j,k)-Gr%volume(i-1,j,k), &
                           Gr%volume(i,j+1,k)-Gr%volume(i,j-1,k), &
                           Gr%volume(i,j,k+1)-Gr%volume(i,j,k-1) ]
                    ddx = [norm2( Gr%xyzcc(i+1,j,k,:)-Gr%xyzcc(i-1,j,k,:) ), &
                           norm2( Gr%xyzcc(i,j+1,k,:)-Gr%xyzcc(i,j-1,k,:) ), &
                           norm2( Gr%xyzcc(i,j,k+1,:)-Gr%xyzcc(i,j,k-1,:) ) ]
                    gdV(IDIR) = ddV(IDIR)/ddx(IDIR) !Gradient dV
                    gdV(JDIR) = ddV(JDIR)/ddx(JDIR)
                    gdV(KDIR) = ddV(KDIR)/ddx(KDIR)

                    gQ(i,j,k) = norm2(gdV)*norm2(dXcc)/Gr%volume(i,j,k)
                enddo
            enddo
        enddo
    end subroutine GridQuality

    !Write restart dump to "ResF" output file
    !NOTE: Rewritten to output full data structure including halos
    subroutine writeH5Res(Model,Gr,oState,State,ResF)
        type(Model_T), intent(inout) :: Model
        type(Grid_T),  intent(in) :: Gr
        type(State_T), intent(in) :: State,oState
        character(len=*), intent(in) :: ResF

        !Reset IO chain
        call ClearIO(IOVars)

        !Main attributes
        call AddOutVar(IOVars,"nOut",Model%IO%nOut)
        call AddOutVar(IOVars,"nRes",Model%IO%nRes)
        call AddOutVar(IOVars,"ts"  ,Model%ts)
        call AddOutVar(IOVars,"t"   ,Model%t)
        call AddOutVar(IOVars,"ot"  ,oState%time)
        call AddOutVar(IOVars,"dt"  ,Model%dt)

        if (Model%dt0 < TINY*10.0) then
            call AddOutVar(IOVars,"dt0"   ,0.0)
        else
            call AddOutVar(IOVars,"dt0"   ,Model%dt0)
        endif

        !Coordinates of corners
        call AddOutVar(IOVars,"X",Gr%xyz(:,:,:,XDIR))
        call AddOutVar(IOVars,"Y",Gr%xyz(:,:,:,YDIR))
        call AddOutVar(IOVars,"Z",Gr%xyz(:,:,:,ZDIR))

        !State variable
        call AddOutVar(IOVars, "Gas",  State%Gas(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,:))        
        call AddOutVar(IOVars,"oGas" ,oState%Gas(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,:))

        if (Model%doMHD) then
            call AddOutVar(IOVars, "magFlux", State%magFlux(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,:))
            call AddOutVar(IOVars,"omagFlux",oState%magFlux(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,:))
            call AddOutVar(IOVars, "Bxyz"   , State%Bxyz(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:  ))
            call AddOutVar(IOVars,"oBxyz"   ,oState%Bxyz(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:  ))
        endif

        if (Model%doSource) then
            !Add source terms to output
            call AddOutVar( IOVars,"Gas0",Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,:) )
        endif

        !Write out, force real precision
        call WriteVars(IOVars,.false.,ResF)

    end subroutine writeH5Res
    
    subroutine readH5Restart(Model,Gr,State,oState,inH5,doResetO,tResetO)
        type(Model_T), intent(inout) :: Model
        type(Grid_T),  intent(inout) :: Gr
        type(State_T), intent(inout) :: State,oState
        character(len=*), intent(in) :: inH5
        logical , intent(in), optional :: doResetO
        real(rp), intent(in), optional :: tResetO

        logical :: doReset,fExist,hasSrc
        real(rp) :: tReset,dt
        integer :: rSpc

        !Test for resetting
        if (present(doResetO)) then
            doReset = doResetO
            if (present(tResetO)) then
                tReset = tResetO
            else
                tReset = 0.0
            endif
        else
            doReset = .false.
            tReset = 0.0
        endif

        write(*,*) 'Reading restart from ', trim(inH5)
        inquire(file=inH5,exist=fExist)
        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open input restart file, exiting'
            stop
        endif

        !Check for oState info
        if ( (.not. ioExist(inH5,"oGas")) .or. (.not. ioExist(inH5,"omagFlux")) ) then
            write(*,*) "Restart file too old, does not include oState!"
            stop
        endif

        !Reset IO chain
        call ClearIO(IOVars)

        call AddInVar(IOVars,"Gas" )
        call AddInVar(IOVars,"oGas")
        call AddInVar(IOVars,"magFlux")
        call AddInVar(IOVars,"omagFlux")
        call AddInVar(IOVars,"nOut",vTypeO=IOINT)
        call AddInVar(IOVars,"nRes",vTypeO=IOINT)
        call AddInVar(IOVars,"ts"  ,vTypeO=IOINT)
        call AddInVar(IOVars,"t"   ,vTypeO=IOREAL)
        call AddInVar(IOVars,"ot"   ,vTypeO=IOREAL)

        !Get data
        call ReadVars(IOVars,.false.,inH5)

        !Find number of species in restart
        rSpc = IOVars(1)%dims(5)-1

        if (Model%nSpc == rSpc) then
            !Restart and State variable agree
            call IOArray5DFill(IOVars,"Gas" ,State %Gas(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,:))
            call IOArray5DFill(IOVars,"oGas",oState%Gas(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,:))

        else if (Model%nSpc > rSpc) then
            !Not enough species in restart, fill as many as possible
            call IOArray5DFill(IOVars,"Gas" ,State %Gas(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,0:rSpc))
            call IOArray5DFill(IOVars,"oGas",oState%Gas(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,0:rSpc))

            !Now initialize to empty remaining species
            State %Gas(:,:,:,1:NVAR,rSpc+1:Model%nSpc) = 0.0
            oState%Gas(:,:,:,1:NVAR,rSpc+1:Model%nSpc) = 0.0

            !Now reaccumulate
            call State2Bulk(Model,Gr,State )
            call State2Bulk(Model,Gr,oState)
        else
            !Too many species in restart, this isn't good
            write(*,*) 'Restart error, more species in restart than room in State!'
            stop
        endif

        !Now handle magnetic fields
        if (Model%doMHD) then
            call IOArray4DFill(IOVars,"magFlux" ,State %magFlux(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,:))
            call IOArray4DFill(IOVars,"omagFlux",oState%magFlux(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,:))
        endif


        !Get main attributes
        dt = GetIOReal(IOVars,"t") - GetIOReal(IOVars,"ot") !Spacing between restart states
        if (doReset) then
            Model%IO%nOut = 0
            Model%IO%nRes = GetIOInt(IOVars,"nRes") + 1
            Model%ts      = 0
            Model%t       = tReset      
        else
            Model%IO%nOut = GetIOInt(IOVars,"nOut")
            Model%IO%nRes = GetIOInt(IOVars,"nRes") + 1
            Model%ts      = GetIOInt(IOVars,"ts")
            Model%t       = GetIOReal(IOVars,"t")
        endif
        
        State %time = Model%t
        oState%time = State%time - dt !Handles tReset

        !Set back to old dt0 if possible
        if (ioExist(inH5,"dt0")) then
            call ClearIO(IOVars)
            call AddInVar(IOVars,"dt0")
            call ReadVars(IOVars,.false.,inH5)

            Model%dt0 = GetIOReal(IOVars,"dt0")

            if (Model%dt0<TINY*10) then
                Model%dt0 = 0.0
            else
                if (Model%isLoud) write(*,*) 'Found dt0, setting to ', Model%dt0*Model%Units%gT0
            endif
        else
            if (Model%isLoud) then
                write(*,*) 'No dt0 found in restart, setting to 0'
                Model%dt0 = 0.0
            endif
        endif

    !Do source term stuff if necessary
        hasSrc = ioExist(inH5,"Gas0")
        if (Model%doSource .and. hasSrc) then
            if (Model%isLoud) then
                write(*,*) 'Found MHD source term data in restart, reading ...'
            endif
            !We want source and it's got some, let's do this thing
            call ClearIO(IOVars)
            call AddInVar(IOVars,"Gas0")
            call ReadVars(IOVars,.false.,inH5)

            rSpc = IOVars(1)%dims(5)-1
            if (Model%nSpc == rSpc) then
                !Restart and Gas0 species agree, do stuff
                call IOArray5DFill(IOVars,"Gas0",Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,:))
            else
                if (Model%isLoud) write(*,*) 'Gas0 is wrong size, ignoring ...'
            endif
        else
            if (Model%isLoud) write(*,*) 'No Gas0 found in restart, starting fresh ...'
        endif !Gas0

    !Do touchup to data structures
        Model%IO%tOut = floor(Model%t/Model%IO%dtOut)*Model%IO%dtOut + Model%IO%dtOut
        Model%IO%tRes = floor(Model%t/Model%IO%dtRes)*Model%IO%dtRes + Model%IO%dtRes

    end subroutine readH5Restart

    !Output black box from crash
    subroutine WriteBlackBox(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        type(State_T), intent(in) :: State

        !Reset output file
        GamH5File = "CRASH" // trim(Model%RunID) // ".h5"
        doRoot = .true. !Make sure root vars are rewritten
        call CheckAndKill(GamH5File)

        call writeSlc(Model,Gr,State,"Step#0")

    end subroutine WriteBlackBox

end module gioH5
