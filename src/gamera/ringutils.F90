!Utilities to handle ring-averaging

!Notation
!nSi/nSe: Slices.  Cyl (Z), LFM (R)
!nR: Rings. Cyl (R), LFM (J)

module ringutils
    use gamtypes
    use gamutils
    use gamdebug
    use math
    use xml_input
    implicit none

    !Enumerators for singularity sides (i.e. positive/negative axis)
    !SPOLE is pole at xs, EPOLE is pole at xe
    enum, bind(C)
        enumerator :: SPOLE=1,EPOLE
    endenum

    integer, parameter :: NFTMAX = 2 !Max number of FFT modes (beyond 0th) possible

    !TODO: Handle these values better for getting them too ringav.F90
    !Information for Fourier reductions
    integer :: NFTMAG = 1 !Number of Fourier modes (beyond 0th) to remove from magflux
    integer :: NFTVEL = 0 !Number of Fourier modes (beyond 0th) to remove from velocity
    logical :: doShift = .false. !Whether to add random circular shift to ring chunking
    logical :: doVClean = .false. !Whether to clean velocity Fourier components before reconstruction
    
    !Cleaning field loops    
    logical :: doCleanLoop = .true. !Whether to remove magnetic field loops
    logical :: doFastLoop  = .true. !Whether to remove magnetic field loops instantly

    contains

    !Init ring averager, change edge length used in timestep calculation
    subroutine InitRings(Model,Grid,xmlInp)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(XML_Input_T), intent(inout) :: xmlInp

        real(rp), dimension(NDIM) :: iHatO,jHatO,kHatO,iHatG,jHatG,kHatG

        integer, dimension(NDIM) :: Nx
        integer :: NCh0, iR, Nr, dJ, iLim
        integer, dimension(:), allocatable :: NCr !Default number of chunks/ring
        character(len=strLen) :: ChunkID
        logical :: doAggressive

        Nx = [Grid%Nip,Grid%Njp,Grid%Nkp]
        Model%Ring%NumR = 1 !Fake value for now

        Model%Ring%doS = .false.
        Model%Ring%doE = .false.

        !Whether or not to do xs/xe sides of singularity
        if (Grid%isTiled) then
            !If we're tiled we need to figure out if we own either singularity
            if (Grid%hasLowerBC(JDIR)) Model%Ring%doS = .true.
            if (Grid%hasUpperBC(JDIR)) Model%Ring%doE = .true.
        else
            !Assume we have both singularities
            Model%Ring%doS = .true.
            Model%Ring%doE = .true.
        endif

    !Set some extra knobs

        !Choice of ring coordinates
        if (Model%doMultiF) then
            !Set mass coordinates as default if doing MF
            call xmlInp%Set_Val(Model%Ring%doMassRA,"ringknobs/doMassRA",.true.)
        else
            call xmlInp%Set_Val(Model%Ring%doMassRA,"ringknobs/doMassRA",.false.)
        endif
                
        call xmlInp%Set_Val(doCleanLoop,"ringknobs/doClean"    ,doCleanLoop)
        call xmlInp%Set_Val(doFastLoop ,"ringknobs/doFastClean",doFastLoop )
        call xmlInp%Set_Val(doShift    ,"ringknobs/doShift",doShift)

        call xmlInp%Set_Val(NFTMAG,"ringknobs/NMag",NFTMAG)
        call xmlInp%Set_Val(doVClean,"ringknobs/doVClean",.false.)
        if (doVClean) then
            call xmlInp%Set_Val(NFTVEL,"ringknobs/NVel",NFTVEL)
        else
            NFTVEL = -1
        endif

        if ( (NFTMAG>NFTMAX) .or. (NFTVEL>NFTMAX) ) then
            if (Model%isLoud) then
                write(*,*) 'Trying to subtract more modes than NFTMAX = ', NFTMAX
                write(*,*) 'Bailing ...'
            endif
            stop
        endif

        !Choose safe vs. dangerous defaults
        call xmlInp%Set_Val(doAggressive,"ring/doAggressive",.false.)
        if (doAggressive .and. Model%isLoud) then
            write(*,*) ANSIRED
            write(*,*) "<Using aggressive ring configuration, that could be dangerous ...>"
            write(*,'(a)',advance="no") ANSIRESET, ''
        endif

    !Set singularity information and default ring configurations
        select case (Model%Ring%GridID)
            !------------------
            case ("cyl")
                !Set pole geometry
                Model%Ring%PLDIR = JDIR
                Model%Ring%Np = Nx(Model%Ring%PLDIR)
                Model%Ring%nSi = Grid%ks
                Model%Ring%nSe = Grid%ke

                !Set chunk/ring info
                call xmlInp%Set_Val(NCh0,"ring/NCh0",8)

                Nr = log(1.0*Model%Ring%Np/NCh0)/log(2.0)
                allocate(NCr(Nr))
                do iR=1,Nr
                    NCr(iR) = NCh0*2**(iR-1)
                enddo

            !------------------
            case ("lfm")
                Model%Ring%PLDIR = KDIR
                Model%Ring%Np = Nx(Model%Ring%PLDIR)                
                Model%Ring%nSi = Grid%is
                Model%Ring%nSe = Grid%ie
                
                !Set chunk/ring info
                select case(Model%Ring%Np)
                case(8)
                    Nr = 1
                    allocate(NCr(Nr))
                    NCr = [4]
                case(16)
                    Nr = 1
                    allocate(NCr(Nr))
                    NCr = [8]
                case(32)
                    Nr = 2
                    allocate(NCr(Nr))
                    NCr = [8,16]
                case(64)
                    !DBL res
                    if (doAggressive) then
                        Nr = 8
                        allocate(NCr(Nr))
                        NCr = [8,16,32,32,32,32,32,32]
                    else
                        Nr = 4
                        allocate(NCr(Nr))
                        NCr = [8,16,32,32]
                    endif
                case(128)
                    !QUAD res
                    if (doAggressive) then
                        !Use possibly dangerous config
                        Nr = 12
                        allocate(NCr(Nr))
                        NCr = [8,16,32,32,64,64,64,64,64,64,64,64]
                    else
                        Nr = 8
                        allocate(NCr(Nr))
                        NCr = [8,16,32,32,64,64,64,64]
                    endif
                case(256)
                    !OCT res
                    if (doAggressive) then
                        !Use old LFM config
                        Nr = 16
                        allocate(NCr(Nr))
                        NCr = [8,16,32,32,64,64,64,64,128,128,128,128,128,128,128,128]
                    else
                        !Use safer gamera defaults
                        Nr = 12
                        allocate(NCr(Nr))
                        NCr = [8,16,32,32,64,64,64,64,128,128,128,128]
                    endif

                case(512)
                    !Default HEX case
                    if (doAggressive) then
                        !Use dangerous hex config
                        Nr = 22
                        allocate(NCr(Nr))
                        ! 8x1,16x1,32x2,64x4,128x6,256x8
                        NCr = [8,16,32,32,64,64,64,64,128,128,128,128,128,128,256,256,256,256,256,256,256,256]
                    else
                        !Use safer defaults
                        Nr = 16
                        allocate(NCr(Nr))
                        ! 8x1,16x1,32x2,64x4,128x4,256x4
                        NCr = [8,16,32,32,64,64,64,64,128,128,128,128,256,256,256,256]
                    endif
                case default
                    write(*,*) 'Unsupported LFM-Nk, no default ...'
                end select
                
            !------------------
            case default
                write(*,*) 'Unknown ring type, bailing ...'
                stop                
        end select

    !Let user over-ride ring config defaults via XML input
        !Set number of rings
        call xmlInp%Set_Val(Model%Ring%NumR,"ring/Nr",Nr)
        allocate(Model%Ring%Nch(Model%Ring%NumR))
        
        !Loop over number of rings and set chunking
        do iR=1,Model%Ring%NumR
            !Set default value for this ring chunk
            iLim = min(iR,Nr) !Don't go over number of rings in set defaults
            dJ = NCr(iLim)
            write(ChunkID,'(A,I0)') "ring/Nc", iR
            call xmlInp%Set_Val(Model%Ring%Nch(iR),trim(ChunkID),dJ)
        enddo
        
        !Do some sanity checks
        select case (Model%Ring%GridID)
            !------------------
            case ("lfm")
                !Test that we haven't over-refined past ring avg chunks
                if (Model%Ring%NumR>=Grid%Njp) then
                    write(*,*) 'Number of rings is larger than J cells per rank!'
                    stop
                endif
                if (.not. (Grid%hasUpperBC(KDIR) .and. Grid%hasLowerBC(KDIR)) ) then
                    write(*,*) 'Ring averaging does not work with tiling in K'
                    write(*,*) 'Somebody (else) should fix that ...'
                    stop
                endif
        end select
    end subroutine InitRings

    !Make corrections to grid data structure to account for pole
    !Correct ghost grid using values from conjugate cell
    !Things to fix:
    !Face-transforms, di,dj,dk, dV
    subroutine RingGridFix(Model,Gr)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Gr

        real(rp), dimension(NDIM) :: iHatO,jHatO,kHatO,iHatG,jHatG,kHatG
        integer :: i,j,k,ig,jg,kg,ip,jp,kp
        integer :: iR,dJ,Np2
        integer :: n,dNorm,T1,T2
        integer :: jxP,jxM

        associate(NumR=>Model%Ring%NumR,Np=>Model%Ring%Np)

        Np2 = Model%Ring%Np/2

        !Increase the size of the edge length
        select case (Model%Ring%GridID)
            case ("cyl")
                do iR=1,NumR
                    dJ = Np/Model%Ring%NCh(iR)
                    Gr%dj(iR,:,:) = dJ*Gr%dj(iR,:,:)
                enddo
            case ("lfm")
                do iR=1,NumR
                    dJ = Np/Model%Ring%NCh(iR)

                    if (Model%Ring%doS) then
                        jxP = Gr%js+iR-1
                        Gr%dk(:,jxP,:) = dJ*Gr%dk(:,jxP,:)   
                    endif
                    
                    if (Model%Ring%doE) then
                        jxM = Gr%je-iR+1
                        Gr%dk(:,jxM,:) = dJ*Gr%dk(:,jxM,:)
    
                    endif
                enddo
        end select

        !ig,jg,kg = ghost cell
        !ip,jp,kp = opposite cell

        !Fix face normal vectors in ghost regions to account for switch from right to left hand system going through pole
        !Need for magnetic flux through pole to have consistent hat vectors for reconstruction on fluxes                
        !For poleward direction, use xp+1 to avoid double-counting face on pole
        select case (Model%Ring%GridID)
        !---------------
        case ("cyl")
            do k=Gr%ksg,Gr%keg
                do j=Gr%jsg,Gr%jeg
                    do n=1,Model%Ng
                        ig = Gr%is-n
                        jg = j; kg = k

                        ip = Gr%is+n-1
                        jp = j+Np2
                        if (jp>Np) jp = jp-Np !Wrap
                        kp = k

                        Gr%volume(ig,jg,kg) = Gr%volume(ip,jp,kp)
                        Gr%di(ig,jg,kg) = Gr%di(ip,jp,kp)
                        Gr%dj(ig,jg,kg) = Gr%dj(ip,jp,kp)
                        Gr%dk(ig,jg,kg) = Gr%dk(ip,jp,kp)
                        
                        Gr%Tf(ig,jg,kg,NORMX:NORMZ,IDIR) = -Gr%Tf(ip+1,jp,kp,NORMX:NORMZ,IDIR)
                        Gr%Tf(ig,jg,kg,NORMX:NORMZ,JDIR) = -Gr%Tf(ip  ,jp,kp,NORMX:NORMZ,JDIR)
                        Gr%Tf(ig,jg,kg,NORMX:NORMZ,KDIR) =  Gr%Tf(ip  ,jp,kp,NORMX:NORMZ,KDIR)

                        ! Gr%face(ig,jg,kg,IDIR) = Gr%face(ip+1,jp,kp,IDIR)
                        ! Gr%face(ig,jg,kg,JDIR) = Gr%face(ip  ,jp,kp,JDIR)
                        ! Gr%face(ig,jg,kg,KDIR) = Gr%face(ip  ,jp,kp,KDIR)
                    enddo
                enddo
            enddo
        !---------------
        case ("lfm")
            do k=Gr%ksg,Gr%keg+1
            !do k=Gr%ks,Gr%ke+1
                do i=Gr%isg,Gr%ieg+1
                    do n=1,Model%Ng                   
                        !j-boundaries (IN)
                        if (Model%Ring%doS) then
                            ig = i
                            kg = k
                            jg = Gr%js-n
                            
                            !Do cell-centered stuff
                            if ( (kg<=Gr%keg) .and. (ig<=Gr%ieg) ) then
                                call lfmIJKcc(Model,Gr,ig,jg,kg,ip,jp,kp)
                                Gr%volume(ig,jg,kg) = Gr%volume(ip,jp,kp)
                                Gr%di(ig,jg,kg) = Gr%di(ip,jp,kp)
                                Gr%dj(ig,jg,kg) = Gr%dj(ip,jp,kp)
                                Gr%dk(ig,jg,kg) = Gr%dk(ip,jp,kp)
                            endif

                            !Do faces
                            !I face (+ signature)
                            call lfmIJKfc(Model,Gr,IDIR,ig,jg,kg,ip,jp,kp)
                            Gr%Tf  (ig,jg,kg,NORMX:NORMZ,IDIR) = Gr%Tf  (ip,jp,kp,NORMX:NORMZ,IDIR)
                            Gr%face(ig,jg,kg,            IDIR) = Gr%face(ip,jp,kp,            IDIR)
                            !J face (- signature)
                            call lfmIJKfc(Model,Gr,JDIR,ig,jg,kg,ip,jp,kp)
                            Gr%Tf  (ig,jg,kg,NORMX:NORMZ,JDIR) = -Gr%Tf  (ip,jp,kp,NORMX:NORMZ,JDIR)
                            Gr%face(ig,jg,kg,            JDIR) =  Gr%face(ip,jp,kp,            JDIR)

                            !K face (- signature)
                            call lfmIJKfc(Model,Gr,KDIR,ig,jg,kg,ip,jp,kp)
                            Gr%Tf  (ig,jg,kg,NORMX:NORMZ,KDIR) = -Gr%Tf  (ip,jp,kp,NORMX:NORMZ,KDIR)
                            Gr%face(ig,jg,kg,            KDIR) =  Gr%face(ip,jp,kp,            KDIR)

                        endif !doS

                        !j-boundaries (OUT)
                        if (Model%Ring%doE) then
                            ig = i
                            kg = k
                            jg = Gr%je+n
                            
                            !Do cell-centered stuff
                            if ( (kg<=Gr%keg) .and. (ig<=Gr%ieg) ) then
                                call lfmIJKcc(Model,Gr,ig,jg,kg,ip,jp,kp)
                                Gr%volume(ig,jg,kg) = Gr%volume(ip,jp,kp)
                                Gr%di(ig,jg,kg) = Gr%di(ip,jp,kp)
                                Gr%dj(ig,jg,kg) = Gr%dj(ip,jp,kp)
                                Gr%dk(ig,jg,kg) = Gr%dk(ip,jp,kp)
                            endif
                            !Do faces
                            !I face (+ signature)
                            call lfmIJKfc(Model,Gr,IDIR,ig,jg,kg,ip,jp,kp)
                            Gr%Tf  (ig,jg,kg,NORMX:NORMZ,IDIR) = Gr%Tf  (ip,jp,kp,NORMX:NORMZ,IDIR)
                            Gr%face(ig,jg,kg,            IDIR) = Gr%face(ip,jp,kp,            IDIR)

                            !J face (- signature), first active is je+1+1
                            call lfmIJKfc(Model,Gr,JDIR,ig,jg+1,kg,ip,jp,kp)
                            Gr%Tf  (ig,jg+1,kg,NORMX:NORMZ,JDIR) = -Gr%Tf  (ip,jp,kp,NORMX:NORMZ,JDIR)
                            Gr%face(ig,jg+1,kg,            JDIR) =  Gr%face(ip,jp,kp,            JDIR)

                            !K face (- signature)
                            call lfmIJKfc(Model,Gr,KDIR,ig,jg,kg,ip,jp,kp)
                            Gr%Tf  (ig,jg,kg,NORMX:NORMZ,KDIR) = -Gr%Tf  (ip,jp,kp,NORMX:NORMZ,KDIR)
                            Gr%face(ig,jg,kg,            KDIR) =  Gr%face(ip,jp,kp,            KDIR)
                        endif !doE

                    enddo !n (ghost) loop
                enddo
            enddo

        end select

        end associate
    end subroutine RingGridFix

    subroutine RingFlux2FieldFix(Model,Grid,Bxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        real(rp), dimension(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM), intent(inout) :: Bxyz
        integer :: n,i,k,ig,jg,kg,ip,jp,kp

        !Kick out the riff raff
        if (.not. Model%doRing) return
        if ( (.not. Model%Ring%doE) .and. (.not. Model%Ring%doS) ) return

        !Okay, let's get shit done
        select case (trim(toUpper(Model%Ring%GridID)))
            !------------------
            case ("CYL")
                write(*,*) "If you care about this configuration, you should write this routine ..."
                stop

            !------------------
            case ("LFM")
                !Loop through cells on the other side of the singularity and replace their Bxyz with the real cell they map to
                
                !$OMP PARALLEL DO default(shared) &
                !$OMP private(n,i,k,ig,jg,kg,ip,jp,kp)
                do k=Grid%ksg,Grid%keg
                    do i=Grid%isg,Grid%ieg
                        do n=1,Model%Ng
                            if (Model%Ring%doS) then
                                ig = i
                                kg = k
                                jg = Grid%js-n
                                call lfmIJKcc(Model,Grid,ig,jg,kg,ip,jp,kp)
                                Bxyz(ig,jg,kg,:) = Bxyz(ip,jp,kp,:)
                            endif

                            if (Model%Ring%doE) then
                                ig = i
                                kg = k
                                jg = Grid%je+n
                                call lfmIJKcc(Model,Grid,ig,jg,kg,ip,jp,kp)
                                Bxyz(ig,jg,kg,:) = Bxyz(ip,jp,kp,:)
                            endif
                        enddo
                    enddo
                enddo !k loop        
                        
            !------------------
            case default
                write(*,*) 'Unknown ring type, bailing ...'
                stop                
        end select        

    end subroutine RingFlux2FieldFix

    !Replace Bxyz's inside singular region with their conjugate cell value
    subroutine HackLFMPredictor(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,i,k,ig,jg,kg,ip,jp,kp

        if ( (.not. Model%Ring%doE) .and. (.not. Model%Ring%doS) ) return

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,i,k,ig,jg,kg,ip,jp,kp)
        do k=Grid%ksg,Grid%keg
            do i=Grid%isg,Grid%ieg
                do n=1,Model%Ng
                    if (Model%Ring%doS) then
                        ig = i
                        kg = k
                        jg = Grid%js-n
                        call lfmIJKcc(Model,Grid,ig,jg,kg,ip,jp,kp)
                        State%Bxyz(ig,jg,kg,:) = State%Bxyz(ip,jp,kp,:)
                    endif

                    if (Model%Ring%doE) then
                        ig = i
                        kg = k
                        jg = Grid%je+n
                        call lfmIJKcc(Model,Grid,ig,jg,kg,ip,jp,kp)
                        State%Bxyz(ig,jg,kg,:) = State%Bxyz(ip,jp,kp,:)
                    endif
                enddo
            enddo
        enddo
    end subroutine HackLFMPredictor

    !Cell-centered conjugate mapping
    subroutine lfmIJKcc(Model,Grid,i,j,k,ip,jp,kp)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        integer, intent(in) :: i,j,k
        integer, intent(out) :: ip,jp,kp

        integer :: Np
        Np  = Grid%Nkp

        !Map i to itself
        ip = i

        !Next do k, map via periodicity
        !NOTE: This is assuming you have all
        if (k < Grid%ks) then
            kp = k + Np
        elseif (k > Grid%ke) then
            kp = k - Np
        else
            kp = k
        endif

        !Finally do j
        jp = j ! default value
        if ( Model%Ring%doS .and. (j<Grid%js) ) then
            !js-1 => js
            jp = Grid%js + (Grid%js-j) - 1
            kp = WrapK(k,Np)
        endif

        if ( Model%Ring%doE .and. (j>Grid%je) ) then
            !je+1 => je
            jp = Grid%je - (j-Grid%je) + 1
            kp = WrapK(k,Np)
        endif

    end subroutine lfmIJKcc

    !d Face-centered conjugate mapping
    subroutine lfmIJKfc(Model,Grid,d,i,j,k,ip,jp,kp)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        integer, intent(in) :: i,j,k,d
        integer, intent(out) :: ip,jp,kp

        integer :: Np
        Np  = Grid%Nkp

        !Map i to itself
        ip = i

        if (d == IDIR) then
            !In i direction just use cell-centered
            call lfmIJKcc(Model,Grid,i,j,k,ip,jp,kp)
            return
        endif

        !Now do k via periodicity
        !NOTE: This is assuming you have all
        kp = k
        !k<ks same for both J/K
        !k>ke (J), k>ke+1 (K)
        if (k < Grid%ks) then
            kp = k + Np
        endif
        if ( (d == JDIR) .and. (k > Grid%ke) ) then
            kp = k - Np
        endif
        if ( (d == KDIR) .and. (k > Grid%ke+1) ) then
            kp = k - Np
        endif

        !Finally do j
        jp = j ! default value
        if (d == JDIR) then
            !Wraps differently
            !For j you offset so you only see the singularity once
            !js-1 => js+1 & wrap K
            !js-2 => js+2 & wrap K
            !For cell center, js-1 => js & wrap K
            if ( Model%Ring%doS .and. (j<Grid%js) ) then
                jp = Grid%js + (Grid%js-j)
                kp = WrapK(k,Np)
            endif
            if ( Model%Ring%doE .and. (j>Grid%je+1) ) then
                !je+1 => je+1
                !je+2 => je
                !je+3 => je-1
                jp = Grid%je + 1 - (j-Grid%je) + 1
                kp = WrapK(k,Np)
            endif
        else !KDIR
            !Wrapping like cell-center
            if ( Model%Ring%doS .and. (j<Grid%js) ) then
                jp = Grid%js + (Grid%js-j) - 1
                kp = WrapK(k,Np)
            endif

            if ( Model%Ring%doE .and. (j>Grid%je) ) then
                jp = Grid%je - (j-Grid%je) + 1
                kp = WrapK(k,Np)
            endif
        endif

    end subroutine lfmIJKfc

    !d Edge-centered conjugate mapping
    subroutine lfmIJKec(Model,Grid,d,i,j,k,ip,jp,kp)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        integer, intent(in) :: i,j,k,d
        integer, intent(out) :: ip,jp,kp

        integer :: Np
        Np  = Grid%Nkp

        !Map i to itself
        ip = i
        
        !Now do k via periodicity
        !NOTE: This is assuming you have all
        kp = k
        if ( (d == IDIR) .or. (d == JDIR) ) then
            !i/j edges in k plane
            if (k < Grid%ks) then
                kp = k + Np
            endif
            if (k > Grid%ke+1) then
                kp = k - Np
            endif
        else 
            !k edge
            if (k < Grid%ks) then
                kp = k + Np
            endif
            if (k > Grid%ke) then
                kp = k - Np
            endif
        endif

        !Finally do j
        jp = j
        if ( (d == IDIR) .or. (d==KDIR) ) then
            if ( Model%Ring%doS .and. (j<Grid%js) ) then
                jp = Grid%js + (Grid%js-j) - 1
                kp = WrapK(k,Np)
            endif
            if ( Model%Ring%doE .and. (j>Grid%je+1) ) then
                !je+2=>je+1
                jp = Grid%je + 1 - (j-Grid%je-1) + 1
                kp = WrapK(k,Np)
            endif
        else
            !J dir
            if ( Model%Ring%doS .and. (j<Grid%js) ) then
                !js-1 => js
                jp = Grid%js + (Grid%js-j) - 1
                kp = WrapK(k,Np)
            endif
            if ( Model%Ring%doE .and. (j>Grid%je) ) then
                !je+1 => je
                jp = Grid%je - (j-Grid%je) + 1
                kp = WrapK(k,Np)
            endif
        endif

    end subroutine lfmIJKec

    function WrapK(k,Np) result(kp)
        integer, intent(in) :: k,Np
        integer :: kp
        integer :: Np2

        Np2 = Np/2
        kp = k + Np2
        if (kp>Np) kp = kp-Np

    end function WrapK

    !Ensure no flux through degenerate faces
    subroutine RingFlux(Model,Gr,gFlx,mFlx)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(inout) :: gFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(inout), optional :: mFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM,1:NDIM)

        select case (Model%Ring%GridID)
        case("lfm")
            if (Model%Ring%doS) gFlx(:,Gr%js  ,:,:,JDIR,:) = 0.0
            if (Model%Ring%doE) gFlx(:,Gr%je+1,:,:,JDIR,:) = 0.0
            if (Model%doMHD .and. present(mFlx)) then
                if (Model%Ring%doS) mFlx(:,Gr%js  ,:,:,JDIR) = 0.0
                if (Model%Ring%doE) mFlx(:,Gr%je+1,:,:,JDIR) = 0.0
                call FixGasFluxLFM(Model,Gr,gFlx,mFlx)
            endif
        end select
    end subroutine RingFlux
    
!-----
!Routines for pulling/pushing and converting variables

    !Push ring back into cell-centered variables
    subroutine PushRingCC(Model,Gr,Qcc,Qr,nR,nS,NumV,XPOLE)
        type (Model_T), intent(in) :: Model
        type (Grid_T) , intent(in) :: Gr
        real(rp), intent(inout) :: Qcc(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NumV)
        integer, intent(in) :: nR,nS,NumV,XPOLE
        real(rp), dimension(Model%Ring%Np,NumV), intent(in) :: Qr

        select case (Model%Ring%GridID)
        case ("cyl")
            !Cylindrical
            if (XPOLE == SPOLE) then
                Qcc(Gr%is+nR-1,Gr%js:Gr%je,nS,:) = Qr
            endif

            if (XPOLE == EPOLE) then
                !No meaningfull data
                return
            endif
        case ("lfm")        
            !LFM-style axis
            if (XPOLE == SPOLE) then
                Qcc(nS,Gr%js+nR-1,Gr%ks:Gr%ke,:) = Qr
            endif
            if (XPOLE == EPOLE) then
                Qcc(nS,Gr%je-nR+1,Gr%ks:Gr%ke,:) = Qr
            endif

        end select
    end subroutine PushRingCC

    !Pull ring from cell-centered variables
    subroutine PullRingCC(Model,Gr,Qcc,Qr,nR,nS,NumV,XPOLE)
        type (Model_T), intent(in) :: Model
        type (Grid_T) , intent(in) :: Gr
        real(rp), intent(in) :: Qcc(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NumV)
        integer, intent(in) :: nR,nS,NumV,XPOLE
        real(rp), dimension(Model%Ring%Np,NumV), intent(inout) :: Qr

        select case (Model%Ring%GridID)
        case ("cyl")
            !Cylindrical
            if (XPOLE == SPOLE) then
                Qr = Qcc(Gr%is+nR-1,Gr%js:Gr%je,nS,:)
            endif

            if (XPOLE == EPOLE) then
                !No meaningfull data
                Qr = 0.0
            endif
        case ("lfm")        
            !LFM-style axis
            if (XPOLE == SPOLE) then
                Qr = Qcc(nS,Gr%js+nR-1,Gr%ks:Gr%ke,:)
            endif
            if (XPOLE == EPOLE) then
                Qr = Qcc(nS,Gr%je-nR+1,Gr%ks:Gr%ke,:)
            endif

        end select
    end subroutine PullRingCC

    !Pull ring from interface-centered fluxes
    subroutine PullRingI(Model,Gr,Qfc,Qr,nR,nS,XPOLE)
        type (Model_T), intent(in) :: Model
        type (Grid_T) , intent(in) :: Gr
        real(rp), intent(in) :: Qfc(Gr%isg:Gr%ieg+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,NDIM)
        integer, intent(in) :: nR,nS,XPOLE
        real(rp), dimension(Model%Ring%Np,NDIM), intent(inout) :: Qr

        select case (Model%Ring%GridID)
        case ("cyl")
            !Cylindrical
            if (XPOLE == SPOLE) then
                !Start at I-flux is+1, rest @ is
                Qr(:,IDIR) = Qfc(Gr%is+nR  ,Gr%js:Gr%je,nS,IDIR)
                Qr(:,JDIR) = Qfc(Gr%is+nR-1,Gr%js:Gr%je,nS,JDIR)
                Qr(:,KDIR) = Qfc(Gr%is+nR-1,Gr%js:Gr%je,nS,KDIR)
            endif

            if (XPOLE == EPOLE) then
                !No meaningfull data
                Qr = 0.0
            endif
        case ("lfm")
            !LFM-style singularity
            if (XPOLE == SPOLE) then
                !For +X pole
                !start at J-Flux js+1
                !start at I-Flux js
                Qr(:,IDIR) = Qfc(nS,Gr%js+nR-1,Gr%ks:Gr%ke,IDIR)
                Qr(:,JDIR) = Qfc(nS,Gr%js+nR  ,Gr%ks:Gr%ke,JDIR)
                Qr(:,KDIR) = Qfc(nS,Gr%js+nR-1,Gr%ks:Gr%ke,KDIR)
            endif
            if (XPOLE == EPOLE) then
                !For -X pole
                !start at J-Flux je
                !start at I-Flux je
                Qr(:,IDIR:KDIR) = Qfc(nS,Gr%je-nR+1,Gr%ks:Gr%ke,IDIR:KDIR)
            endif

        end select        
    end subroutine PullRingI

!-----
!Convert to ringav variables and back
    !Convert ringav variables back to hydro variables
    !Which ring variables
    !doMassRA = T => rho,rho*v,rho*Cs^2
    !doMassRA = F => rho,mom  ,inte

    !Con -> RAVars
    subroutine Gas2Ring(Model,rW)
        type (Model_T), intent(in) :: Model
        real(rp), intent(inout) :: rW(Model%Ring%Np,NVAR)

        integer :: n
        real(rp) :: D,E,P,KinE,IntE
        real(rp), dimension(NDIM) :: Mom

        do n=1,Model%Ring%Np
            D = max( rW(n,DEN), dFloor )
            Mom = rW(n,MOMX:MOMZ)
            E = rW(n,ENERGY)

            KinE = 0.5*dot_product(Mom,Mom)/D
            P = max( (Model%gamma-1)*(E-KinE) , pFloor )
            IntE = P/(Model%gamma-1)
            !Put ring variables back in (in place)
            rW(n,DEN) = D
            rW(n,MOMX:MOMZ) = Mom
            if (.not. Model%Ring%doMassRA) then
                rW(n,ENERGY) = IntE
            else
                rW(n,ENERGY) = (Model%gamma)*P
            endif

        enddo

    end subroutine Gas2Ring

    !RAVars => Con
    subroutine Ring2Gas(Model,rW)
        type (Model_T), intent(in) :: Model
        real(rp), intent(inout) :: rW(Model%Ring%Np,NVAR)

        integer :: n
        real(rp) :: D,P,KinE,IntE,Cs2
        real(rp), dimension(NDIM) :: Mom,V

        do n=1,Model%Ring%Np
            D  = max( rW(n,DEN), dFloor )
            Mom = rW(n,MOMX:MOMZ)
            if (.not. Model%Ring%doMassRA) then
                IntE = rW(n,PRESSURE)
                P = (Model%gamma-1)*IntE
            else
                Cs2 = rW(n,PRESSURE)/D
                P = Cs2*D/Model%gamma
            endif
            P = max(P,pFloor)
            KinE = 0.5*dot_product(Mom,Mom)/D
            !Recalculate IntE w/ floored pressure
            IntE = P/(Model%gamma-1)
            !Put conserved variables back
            rW(n,DEN) = D
            rW(n,MOMX:MOMZ) = Mom
            rW(n,ENERGY) = KinE+IntE

        enddo

    end subroutine Ring2Gas

end module ringutils
