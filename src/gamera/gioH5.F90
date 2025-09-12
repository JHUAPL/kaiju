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
    use volttypes

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

        integer :: i,j,k,s,nClkSteps
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
                !$OMP PARALLEL DO default(shared) collapse(2) &
                !$OMP private(i,j,k)
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
                !$OMP PARALLEL DO default(shared) collapse(2) &
                !$OMP private(i,j,k)
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

            if (Model%doSource .and. Model%isMagsphere) then
                
                !Volt variables
                call GameraOut("SrcX1" ,"deg",rad2deg,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,PROJLAT))
                call GameraOut("SrcX2" ,"deg",rad2deg,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,PROJLON))

                call GameraOut("SrcIONEx" ,gamOut%eID,gamOut%eScl,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,IONEX))
                call GameraOut("SrcIONEy" ,gamOut%eID,gamOut%eScl,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,IONEY))
                call GameraOut("SrcIONEz" ,gamOut%eID,gamOut%eScl,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,IONEZ))

                !IMAG variables
                call GameraOut("SrcD_RING" ,gamOut%dID,gamOut%dScl,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,IM_D_RING))
                call GameraOut("SrcP_RING" ,gamOut%pID,gamOut%pScl,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,IM_P_RING))
                call GameraOut("SrcD_COLD" ,gamOut%dID,gamOut%dScl,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,IM_D_COLD))
                call GameraOut("SrcP_COLD" ,gamOut%pID,gamOut%pScl,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,IM_P_COLD))
                call GameraOut("SrcDT"     ,"s"       ,gamOut%tScl,Gr%Gas0(iMin:iMax,jMin:jMax,kMin:kMax,IM_TSCL  ))
            
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
        !Performance metrics

        nClkSteps = readNCalls('Advance')
        call AddOutVar(IOVars,"_perf_stepTime",readClock(1)/nClkSteps)
        call AddOutVar(IOVars,"_perf_mathTime", readClock('Gamera')/nClkSteps)
        call AddOutVar(IOVars,"_perf_bcTime", readClock('BCs')/nClkSteps)
        call AddOutVar(IOVars,"_perf_haloTime", readClock('Halos')/nClkSteps)
        call AddOutVar(IOVars,"_perf_voltTime", readClock('VoltSync')/nClkSteps)
        call AddOutVar(IOVars,"_perf_ioTime", readClock('IO')/nClkSteps)
        call AddOutVar(IOVars,"_perf_advanceTime", readClock('Advance')/nClkSteps)

        !----------------------

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
    subroutine writeH5Res(Model,Gr,State,oState,ooState,ResF)
        type(Model_T), intent(inout) :: Model
        type(Grid_T),  intent(in) :: Gr
        type(State_T), intent(in) :: State,oState,ooState
        character(len=*), intent(in) :: ResF

        !Reset IO chain
        call ClearIO(IOVars)

        !Main attributes
        call AddOutVar(IOVars,"nOut",Model%IO%nOut)
        call AddOutVar(IOVars,"nRes",Model%IO%nRes)
        call AddOutVar(IOVars,"ts"  ,Model%ts)
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

        !States
        call AddState2IO(Model,Gr, State,"")
        call AddState2IO(Model,Gr,oState,"o")
        if (Model%doAB3) then
            call AddState2IO(Model,Gr,ooState,"oo")
        endif

        if (Model%doSource) then
            !Add source terms to output
            call AddOutVar( IOVars,"Gas0",Gr%Gas0(:,:,:,:) )
        endif

        if(associated(Model%HackSaveRestart)) then
            ! Add optional variables from IC
            call Model%HackSaveRestart(Model,Gr,State,IOVars)
        endif

        !Write out, force real precision
        call WriteVars(IOVars,.false.,ResF)

        contains
            subroutine AddState2IO(Model,Gr,xState,xID)
                type(Model_T)   , intent(in) :: Model
                type(Grid_T)    , intent(in) :: Gr
                type(State_T)   , intent(in) :: xState
                character(len=*), intent(in) :: xID

                call AddOutVar(IOVars, trim(xID) // "t"  , xState%time)
                call AddOutVar(IOVars, trim(xID) // "Gas", xState%Gas(:,:,:,:,:))
                if (Model%doMHD) then
                    call AddOutVar(IOVars, trim(xID) // "magFlux",xState%magFlux(:,:,:,:))
                    call AddOutVar(IOVars, trim(xID) // "Bxyz"   ,xState%Bxyz   (:,:,:,:)) 
                endif

            end subroutine AddState2IO
    end subroutine writeH5Res
    
    !Read H5 restart
    subroutine readH5Restart(Model,Gr,State,oState,ooState,inH5)
        type(Model_T), intent(inout) :: Model
        type(Grid_T),  intent(inout) :: Gr
        type(State_T), intent(inout) :: State,oState,ooState
        character(len=*), intent(in) :: inH5

        logical  :: hasSrc,hasO,hasOO,hasdt0,skipGhosts
        integer  :: rSpc,n0

        !NOTE: Removing fortran support for changing numfluids on restart, it's easier to just do that on python side

        if (Model%isLoud) write(*,*) 'Reading restart from ', trim(inH5)

        !Check file
        call CheckFileOrDie(inH5,"Unable to open restart file")

        !Okay, let's see what we got here
        hasO   = ioExist(inH5, "oGas")
        hasOO  = ioExist(inH5,"ooGas")
        hasSrc = ioExist(inH5,"Gas0")
        hasdt0 = ioExist(inH5,"dt0")

        !Do some error/warning checking
        if (.not. hasO) then
            write(*,*) "Restart file too old, does not include oState!"
            stop
        endif
        if (Model%doSource .and. (.not. hasSrc)) then
            if (Model%isLoud) write(*,*) 'No Gas0 found in restart, starting fresh ...'
        endif
        if (Model%doAB3 .and. (.not. hasOO)) then
            if (Model%isLoud) write(*,*) "No ooState found in restart, faking it ..."
        endif

        !Now setup chain
        call ClearIO(IOVars)
        call AddInVar(IOVars,"nOut",vTypeO=IOINT)
        call AddInVar(IOVars,"nRes",vTypeO=IOINT)
        call AddInVar(IOVars,"ts"  ,vTypeO=IOINT)

        call AddInVar(IOVars,"t"   ,vTypeO=IOREAL)
        call AddInVar(IOVars,"dt"  ,vTypeO=IOREAL)
        call AddInVar(IOVars,"Gas" )
        call AddInVar(IOVars,"magFlux")

        call AddInVar(IOVars,"ot"   ,vTypeO=IOREAL)
        call AddInVar(IOVars,"oGas")
        call AddInVar(IOVars,"omagFlux")

        if (Model%doAB3 .and. hasOO) then
            call AddInVar(IOVars,"oot"   ,vTypeO=IOREAL)
            call AddInVar(IOVars,"ooGas")
            call AddInVar(IOVars,"oomagFlux")
        endif
        if (Model%doSource .and. hasSrc) then
            call AddInVar(IOVars,"Gas0")
        endif
        if (hasdt0) call AddInVar(IOVars,"dt0")

        !Get data
        call ReadVars(IOVars,.false.,inH5)

        !Now process data
        n0 = FindIO(IOVars,"Gas")
        rSpc = IOVars(n0)%dims(5)-1 !Number of species in restart (subtract BLK)
        if (Model%nSpc /= rSpc) then
            write(*,*) "Restart species is inconsistent w/ MHD, bailing ..."
            stop
        endif

        !n0 still Gas
        if(IOVars(n0)%dims(1) .eq. size(State%Gas,1)) then
            skipGhosts = .false.
        elseif(IOVars(n0)%dims(1) .eq. (size(State%Gas,1) - 2*Model%nG)) then
            write(*,*) "Reading older version gamera restart data without ghosts"
            skipGhosts = .true.
        else
            write(*,*) "Restart data was written from a grid with a different size. Cannot read. Bailing."
            stop
        endif
        
        !Fill state/ostate
        call PullState(Model,Gr, State,"")
        call PullState(Model,Gr,oState,"o")

        !Get main attributes
        Model%IO%nOut = GetIOInt(IOVars,"nOut")
        Model%IO%nRes = GetIOInt(IOVars,"nRes") + 1
        Model%ts      = GetIOInt(IOVars,"ts")
        Model%t       = GetIOReal(IOVars,"t")
        Model%dt      = GetIOReal(IOVars,"dt")

        State %time = Model%t
        oState%time = GetIOReal(IOVars,"ot")

        !Now handle AB3 state
        if (Model%doAB3) then
            if (hasOO) then
                call PullState(Model,Gr,ooState,"oo")
            else
                !Fake it
                ooState%time    = oState%time - Model%dt
                ooState%Gas     = oState%Gas
                if (Model%doMHD) ooState%magFlux = oState%magFlux
            endif
        endif

        !Handle gas0
        if (Model%doSource .and. hasSrc) then
            if(skipGhosts) then
                call IOArray4DFill(IOVars,"Gas0",Gr%Gas0(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:))
            else
                call IOArray4DFill(IOVars,"Gas0",Gr%Gas0(:,:,:,:))
            endif
        endif

        if (hasdt0) then
            Model%dt0 = GetIOReal(IOVars,"dt0")

            if (Model%dt0<TINY*10) Model%dt0 = 0.0
        else
            Model%dt0 = 0.0
        endif

        !Do touchup to data structures
        Model%IO%tOut = floor(Model%t/Model%IO%dtOut)*Model%IO%dtOut + Model%IO%dtOut
        Model%IO%tRes = floor(Model%t/Model%IO%dtRes)*Model%IO%dtRes + Model%IO%dtRes

        if(associated(Model%HackLoadRestart)) then
            ! Load optional variables from IC
            call Model%HackLoadRestart(Model,Gr,State,inH5)
        endif


        contains
            subroutine PullState(Model,Gr,xState,xID)
                type(Model_T)   , intent(in) :: Model
                type(Grid_T)    , intent(in) :: Gr
                type(State_T)   , intent(inout) :: xState
                character(len=*), intent(in) :: xID

                xState%time = GetIOReal(IOVars,trim(xID) // "t")

                if(skipGhosts) then
                    call IOArray5DFill(IOVars,trim(xID) // "Gas",xState%Gas(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:,:))
                else
                    call IOArray5DFill(IOVars,trim(xID) // "Gas",xState%Gas(:,:,:,:,:))
                endif

                if (Model%doMHD) then
                    if (skipGhosts) then
                        call IOArray4DFill(IOVars,trim(xID) // "magFlux",xState%magFlux(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,:))
                    else
                        call IOArray4DFill(IOVars,trim(xID) // "magFlux",xState%magFlux(:,:,:,:))
                    endif
                endif

            end subroutine PullState

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
