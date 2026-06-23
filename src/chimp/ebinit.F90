!Various EB initialization routines
module ebinit
    use chmpdefs
    use chmpunits
    use chmpdbz
    use ebtypes
    use xml_input
    use ioH5
    use chmpfields
    use plasmaputils

    use volttypes, only : NVARVOLTSRC
    
    implicit none
    integer, parameter :: Ngm = 4 !Order for metric differencing

    contains

    ! A version of the above
    !Reads from MHD Grid received via coupling
    subroutine ebInit_fromMHDGrid(Model,ebState,inpXML,mhdGridCorners,doStaticO)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(inout)   :: ebState
        type(XML_Input_T), intent(inout) :: inpXML
        real(rp), dimension(:,:,:,:), intent(in) :: mhdGridCorners
        logical, optional, intent(in) :: doStaticO

        logical :: doStatic

        if (present(doStaticO)) then
            doStatic = doStaticO
        else
            doStatic = .true.
        endif

        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        ! get grid from MHD corners
        call getGridFromMHD(Model,ebGr,inpXML,mhdGridCorners)

        !Allocate eb data (includes space for ghosts)
        !Always do eb1
        call allocEB(Model,ebGr,eb1)
        eb1%E  = 0.0
        eb1%dB = 0.0

        if (doStatic) then
            ebState%doStatic = .true.
        else
            call allocEB(Model,ebGr,eb2)

            eb2%E  = 0.0
            eb2%dB = 0.0
            ebState%doStatic = .false.
        endif

        !Initialize and allocate plasmapause parameters
        if (Model%doPP) call initPP(ebState,inpXML,doStatic)

        end associate

    end subroutine ebInit_fromMHDGrid

    subroutine setGrid(Model,ebGr)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(inout) :: ebGr

        !Calculate augmented sizes
        ebGr%Ni = ebGr%Nip+2*Model%Ng
        ebGr%Nj = ebGr%Njp+2*Model%Ng
        ebGr%Nk = ebGr%Nkp+2*Model%Ng
        !Calculate grid bounds
        ebGr%is = 1; ebGr%ie = ebGr%Nip
        ebGr%js = 1; ebGr%je = ebGr%Njp
        ebGr%ks = 1; ebGr%ke = ebGr%Nkp
        ebGr%isg = ebGr%is-Model%nG
        ebGr%ieg = ebGr%ie+Model%nG
        ebGr%jsg = ebGr%js-Model%nG
        ebGr%jeg = ebGr%je+Model%nG
        ebGr%ksg = ebGr%ks-Model%nG
        ebGr%keg = ebGr%ke+Model%nG

        write(*,'(a,I0,a,I0,a,I0,a)') '<Grid size = (',ebGr%Nip,',',ebGr%Njp,',',ebGr%Nkp,')>'

        !Allocate data

        !NOTE: Leaving corners same size as cells b/c of redundant ghosts here
        allocate(ebGr%xyz  (ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        allocate(ebGr%xyzcc(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        allocate(ebGr%B0cc (ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        allocate(ebGr%Tix  (ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM,NDIM))
        allocate(ebGr%Txi  (ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM,NDIM))
        allocate(ebGr%dV   (ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg))

        !Zero out initial values
        ebGr%xyz = 0.0 ; ebGr%xyzcc = 0.0
        ebGr%B0cc = 0.0
        ebGr%Tix = 0.0 ; ebGr%Txi = 0.0
        ebGr%dV = 0.0

    end subroutine setGrid

    subroutine fixGrid(Model,ebGr,inpXML)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(inout) :: ebGr
        type(XML_Input_T), intent(inout) :: inpXML

        character(len=strLen) :: grStr

        !Figure out grid type & fix corners if necessary
        call inpXML%Set_Val(grStr,"fields/grType","EGG")

        select case(trim(toUpper(grStr)))
        case("EGG")
            ebGr%GrID = EGGGRID
            write(*,*) '<Grid type = EGG>'
            call FixCorners(Model,ebGr)

        case("LFM")
            ebGr%GrID = LFMGRID
            write(*,*) '<Grid type = LFM>'
            call FixCorners(Model,ebGr)
        case("SPH")
            ebGr%GrID = SPHGRID
            write(*,*) '<Grid type = SPH>'
        case default
            ebGr%GrID = EGGGRID
            write(*,*) '<Unknown grid type, assuming EGG>'
        end select

        write(*,*) 'Doing grid calculations ...'
        call Corners2ebG(Model,ebGr)
        write(*,*) 'Finished grid calculations'
    end subroutine fixGrid
        
    !Read grid from file(s)
    subroutine rdGrid(Model,ebGr,ebTab,inpXML)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(inout) :: ebGr
        type(ioTab_T), intent(in) :: ebTab
        type(XML_Input_T), intent(inout) :: inpXML

        integer :: i,j,k,is,ie,js,je,ks,ke
        character(len=strLen) :: ebFile
        integer :: dims(NDIM),dimscc(NDIM)


        ebGr%Nip = ebTab%Ri*ebTab%dNi 
        ebGr%Njp = ebTab%Rj*ebTab%dNj 
        ebGr%Nkp = ebTab%Rk*ebTab%dNk 
        !requires that the above be set (Nip,Njp,Nkp)
        call setGrid(Model,ebGr)

        dims   = [ebTab%dNi+1,ebTab%dNj+1,ebTab%dNk+1]
        dimscc = [ebTab%dNi  ,ebTab%dNj  ,ebTab%dNk  ]

        !------------
        !Loop over grid pieces and get sub-grids
        !NOTE: These reshape commands can cause seg faults if they blow up your stack size
        !Either increase the stack size or rewrite them yourself
        write(*,*) 'Reshaping grid data (may blow up stack) ...'

        do k=1,ebTab%Rk
            do j=1,ebTab%Rj
                do i=1,ebTab%Ri
                    ebFile = genName(ebTab%bStr,ebTab%Ri,ebTab%Rj,ebTab%Rk,i,j,k,Model%doOldNaming)
                    !write(*,'(3a)') '<Reading grid from ', trim(ebFile), '>'

                    !Get piece from file
                    call ClearIO(ebIOs)
                    call AddInVar(ebIOs,"X")
                    call AddInVar(ebIOs,"Y")
                    call AddInVar(ebIOs,"Z")
                    call AddInVar(ebIOs,"dV")

                    call ReadVars(ebIOs,.false.,ebFile) !Use IO precision

                    !Push piece to grid
                    is = (i-1)*ebTab%dNi + 1
                    js = (j-1)*ebTab%dNj + 1
                    ks = (k-1)*ebTab%dNk + 1
                    ie = is + ebTab%dNi - 1
                    je = js + ebTab%dNj - 1
                    ke = ks + ebTab%dNk - 1

                    ebGr%xyz(is:ie+1,js:je+1,ks:ke+1,XDIR) = reshape(ebIOs(XDIR)%data,dims)
                    ebGr%xyz(is:ie+1,js:je+1,ks:ke+1,YDIR) = reshape(ebIOs(YDIR)%data,dims)
                    ebGr%xyz(is:ie+1,js:je+1,ks:ke+1,ZDIR) = reshape(ebIOs(ZDIR)%data,dims)
                    ebGr%dV (is:ie  ,js:je  ,ks:ke       ) = reshape(ebIOs(NDIM+1)%data,dimscc)

                enddo
            enddo
        enddo
        write(*,*) 'Finished reshaping grids'

        !Pull corners and calculate derived geometric quantities
        call fixGrid(Model,ebGr,inpXML)

    end subroutine rdGrid


    subroutine getGridFromMHD(Model,ebGr,inpXML,mhdGridCorners)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(inout) :: ebGr
        type(XML_Input_T), intent(inout) :: inpXML
        real(rp), dimension(:,:,:,:), intent(in) :: mhdGridCorners

        integer :: dims(NDIM)
        integer :: is,ie,js,je,ks,ke
        integer :: isg,ieg,jsg,jeg,ksg,keg

        ! get grid info from MHD
        if (NDIM/=3) then 
           print *,'NDIM!=3 not supported for grids received from MHD'
        end if
        dims = shape(mhdGridCorners(:,:,:,1))

        !Grid metadata
        !Calculate active cells from nodes
        ebGr%Nip = dims(IDIR)-1
        ebGr%Njp = dims(JDIR)-1
        ebGr%Nkp = dims(KDIR)-1

        ! requires that the above be set (Nip,Njp,Nkp)
        call setGrid(Model,ebGr)

        !Pull grid bounds to local variables for brevity
        is = ebGr%is ; ie = ebGr%ie ; isg = ebGr%isg ; ieg = ebGr%ieg
        js = ebGr%js ; je = ebGr%je ; jsg = ebGr%jsg ; jeg = ebGr%jeg
        ks = ebGr%ks ; ke = ebGr%ke ; ksg = ebGr%ksg ; keg = ebGr%keg

        ebGr%xyz(is:ie+1,js:je+1,ks:ke+1,:) = mhdGridCorners

        !Pull corners and calculate derived geometric quantities
        call fixGrid(Model,ebGr,inpXML)
    end subroutine getGridFromMHD

    subroutine Corners2ebG(Model,ebGr)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(inout) :: ebGr

        integer :: n,i,j,k
        integer :: is,ie,js,je,ks,ke
        integer :: isg,ieg,jsg,jeg,ksg,keg
        real(rp) :: xcc(NDIM)
        !Values for LFM/EGG metric calculations
        real(rp) :: ThX,Rx,dth,dzp,dyp


        !Pull grid bounds to local variables for brevity
        is = ebGr%is ; ie = ebGr%ie ; isg = ebGr%isg ; ieg = ebGr%ieg
        js = ebGr%js ; je = ebGr%je ; jsg = ebGr%jsg ; jeg = ebGr%jeg
        ks = ebGr%ks ; ke = ebGr%ke ; ksg = ebGr%ksg ; keg = ebGr%keg

        !Calculate derived quantities over active grid
        !Cell-centers (need to get all before differencing), B0 @ cell-centers
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(xcc)
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    !Calculate cell center (simple averaging)
                    ebGr%xyzcc(i,j,k,:) = 0.125*(  ebGr%xyz(i,j,k,:)   + ebGr%xyz(i+1,j,k,:) &
                                                 + ebGr%xyz(i,j+1,k,:) + ebGr%xyz(i,j,k+1,:) &
                                                 + ebGr%xyz(i+1,j+1,k,:) + ebGr%xyz(i+1,j,k+1,:) &
                                                 + ebGr%xyz(i,j+1,k+1,:) + ebGr%xyz(i+1,j+1,k+1,:) )
                    xcc = ebGr%xyzcc(i,j,k,:)
                    ebGr%B0cc(i,j,k,:) = Model%B0(xcc)
                enddo
            enddo
        enddo

        !Now go back and difference for metric terms        
        !Calculate derivatives of spatial coordinates (@ cc's) wrt to index coordinates
        !Do centered differences except for boundaries

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,n,xcc,ThX,dth,Rx,dyp,dzp)
        do k=ks,ke
            do j=js,je
                do i=is,ie
                    do n=1,NDIM
                        !I derivatives
                        ebGr%Txi(i,j,k,n,IDIR) = Diff1D(Model,ebGr%xyzcc(:,j,k,n),is,ie,i)
                        !J derivatives
                        ebGr%Txi(i,j,k,n,JDIR) = Diff1D(Model,ebGr%xyzcc(i,:,k,n),js,je,j)
                        !K derivatives
                        ebGr%Txi(i,j,k,n,KDIR) = Diff1D(Model,ebGr%xyzcc(i,j,:,n),ks,ke,k)
                    enddo

                    !Do some fixes for geometry
                    !NOTE: This probably doesn't matter
                    if (ebGr%GrID == EGGGRID .or. ebGr%GrID == LFMGRID) then
                        xcc = ebGr%xyzcc(i,j,k,:)
                        ThX = atan2(xcc(ZDIR),xcc(YDIR))
                        dth = 2*PI/ebGr%Nkp
                        Rx = sqrt(xcc(YDIR)**2.0 + xcc(ZDIR)**2.0)
                        dyp = -Rx*sin(ThX)*dth
                        dzp =  Rx*cos(ThX)*dth
                        !write(*,*) 'oTrk = ', ebGr%Txi(i,j,k,XDIR:ZDIR,KDIR)
                        ebGr%Txi(i,j,k,XDIR:ZDIR,KDIR) = [0.0_rp,dyp,dzp]
                        !write(*,*) 'cTrk = ', ebGr%Txi(i,j,k,XDIR:ZDIR,KDIR)
                    endif
                    !Now calculate inverse matrix
                    call matinv3(ebGr%Txi(i,j,k,:,:),ebGr%Tix(i,j,k,:,:))

                enddo
            enddo
        enddo

    end subroutine Corners2ebG
    !Takes 1D differencing of variable Q(is-Model%Ng,ie+Model%Ng) at i0
    !Use 4-point stencil to calculate first derivative of coordinates
    function Diff1D(Model,Q,is,ie,i0) result(Qp)
        integer, intent(in) :: is,ie,i0
        type(chmpModel_T), intent(in) :: Model
        real(rp), intent(in) :: Q(is-Model%Ng:ie+Model%Ng)
        real(rp) :: Qp
        real(rp) :: Qblk(Ngm),c(Ngm)
        if (i0 == is) then
            !Forward
            Qblk = [Q(is),Q(is+1),Q(is+2),Q(is+3)]
            c = [-11.0,18.0,-9.0,2.0]/6.0
        else if (i0 == is+1) then
            !1 back
            Qblk = [Q(is),Q(is+1),Q(is+2),Q(is+3)]
            c = [-2.0,-3.0,6.0,-1.0]/6.0
        else if (i0 == ie) then
            Qblk = [Q(ie-3),Q(ie-2),Q(ie-1),Q(ie)]
            c = [-2.0,9.0,-18.0,11.0]/6.0
        else if (i0 == ie-1) then
            Qblk = [Q(ie-3),Q(ie-2),Q(ie-1),Q(ie)]
            c = [1.0,-6.0,3.0,2.0]/6.0
        else
            !Centered
            Qblk = [Q(i0-2),Q(i0-1),Q(i0+1),Q(i0+2)]
            c = [1.0,-8.0,8.0,-1.0]/12.0
        endif
        Qp = dot_product(Qblk,c)
    end function Diff1D

    !Fix periodicity/matching from single precision
    subroutine FixCorners(Model,ebGr)
        type(chmpModel_T), intent(in)    :: Model
        type(   ebGrid_T), intent(inout) :: ebGr

        integer :: Nkp,Nk2
        ebGr%xyz(:,:,1,ZDIR) = 0.0
        Nkp = ebGr%Nkp
        Nk2 =  Nkp/2
        ebGr%xyz(:,:,1+Nk2,ZDIR) = 0.0
        ebGr%xyz(:,:,Nkp+1,ZDIR) = 0.0

    end subroutine FixCorners

    !Allocates ebField_T given ebGr
    subroutine allocEB(Model,ebGr,ebF)
        type(chmpModel_T), intent(inout) :: Model
        type(ebGrid_T), intent(in)  :: ebGr
        type(ebField_T), intent(inout) :: ebF

        ebF%time = 0.0
        if (.not. allocated(ebF%dB)) then
            allocate(ebF%dB(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        endif

        if (.not. allocated(ebF%E)) then
            allocate(ebF%E(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        endif

        ebF%dB = 0.0
        ebF%E  = 0.0
        if (Model%doMHD .and. (.not. allocated(ebF%W))) then
            allocate(ebF%W(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NVARMHD,0:Model%nSpc))
            ebF%W = 0.0
        endif
        if (Model%doJ .and. (.not. allocated(ebF%Jxyz))) then
            allocate(ebF%Jxyz(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
            ebF%Jxyz = 0.0
        endif
        if (Model%doReadSrc .and. (.not. allocated(ebF%src))) then
            allocate(ebF%src(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NVARVOLTSRC))
            ebF%src = 0.0
        endif


    end subroutine allocEB

end module ebinit
