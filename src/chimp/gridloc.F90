!Routines to deal with grid localization
module gridloc
    use chmpdefs
    use chmpunits
    use ebtypes
    use math
    use xml_input
    use strings

    implicit none

    !Globals for in/out domain
    !Simple spherical bounds for now
    real(rp) :: DomR(2) !Rin/Rout
    real(rp) :: DomE(3) !xSun,xTail,yMax

    real(rp) :: buffR ! extent of Gap region, assume pure dipole field within this and rmin

    !Aux data for grid localization
    type LocAux_T
        real(rp) :: dPhi,dTh !Angular spacing
        real(rp), allocatable :: rrI(:,:) !R-Interface centered radius
        real(rp), allocatable :: xxC(:,:), yyC(:,:) !Cell-centered 2D cells (LFM grids)
        real(rp), allocatable, dimension(:) :: rMin,rMax,pMin,pMax
        logical :: isInit = .false. !Has been initialized
    end type LocAux_T

    !General function type for localiztion
    !Locate point (xyz) in grid cell (ijk)
    !Optionally return isIn with inDomain info
    !ijkO (Optional): Guess as to location of xyz in grid
    !NOTE: ijkO may not actually do anything
    abstract interface
        subroutine Loc_T(xyz,ijk,Model,ebGr,isInO,ijkO)
            Import :: rp,NDIM,chmpModel_T,ebGrid_T
            real(rp), intent(in) :: xyz(NDIM)
            integer, intent(out) :: ijk(NDIM)
            type(chmpModel_T), intent(in) :: Model
            type(ebGrid_T), intent(in) :: ebGr
            logical, intent(out), optional :: isInO
            integer, intent(in), optional :: ijkO(NDIM)
        end subroutine Loc_T
    end interface

    !General inDomain function
    abstract interface
        function inDom_T(xyz,Model,ebGr) result(inDom)
            import :: rp,NDIM,chmpModel_T,ebGrid_T
            real(rp), intent(in) :: xyz(NDIM)
            type(chmpModel_T), intent(in) :: Model
            type(ebGrid_T), intent(in) :: ebGr
            logical :: inDom        
        end function inDom_T
    end interface

    procedure(Loc_T)  , pointer :: locate=>NULL()
    procedure(inDom_T), pointer :: inDomain=>inDomain_Sph
    procedure(inDom_T), pointer :: inGap=>inGap_Default

    type(LocAux_T) :: locAux

    integer, parameter :: NSnake = 25
    !Snake search
    !GFEDC
    !H432B
    !I501A
    !J6789
    !KLMNO
    !                                                            1  2  3  4  5  6  7  8  9  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O
    integer, dimension(NSnake), parameter, private :: dsI = [ 0,+1,+1, 0,-1,-1,-1, 0,+1,+2,+2,+2,+2,+1, 0,-1,-2,-2,-2,-2,-2,-1, 0,+1,+2]
    integer, dimension(NSnake), parameter, private :: dsJ = [ 0, 0,+1,+1,+1,-1,-1,-1,-1,-1, 0,+1,+2,+2,+2,+2,+2,+1, 0,-1,-2,-2,-2,-2,-2]
    
    contains

    !If necessary, deallocate arrays in locAux
    subroutine CleanupLoc()
        if(allocated(locAux%rMin)) deallocate(locAux%rMin)
        if(allocated(locAux%rMax)) deallocate(locAux%rMax)
        if(allocated(locAux%pMin)) deallocate(locAux%pMin)
        if(allocated(locAux%pMax)) deallocate(locAux%pMax)
        if(allocated(locAux%rrI))  deallocate(locAux%rrI)
        if(allocated(locAux%xxC))  deallocate(locAux%xxC)
        if(allocated(locAux%yyC))  deallocate(locAux%yyC)
    end subroutine CleanupLoc

    !Initialize various things for grid localization
    subroutine InitLoc(Model,ebGr,inpXML)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in)   :: ebGr
        type(XML_Input_T), intent(in) :: inpXML

        integer :: i,j
        real(rp) :: xJ(NDIM),xJp(NDIM)
        character(len=strLen) :: domStr,gapStr,gapDefStr

        call cleanupLoc()

        !Setting default of using a dipole in the gap region between grid inner boundary and ionosphere
        ! not used unless it MAGE/Earth simulations
        gapDefStr = "NONE"

        !Initialize aux grid variables for inDomain function
        select case(ebGr%GrID)
        case(LFMGRID,EGGGRID)
            ! if MAGE/Earth run turn on gap region be default
            if (Model%isMAGE .or. trim(toUpper(Model%uID))=="EARTH") then
                gapDefStr = "LFM"

                ! Set Gap radius from sunward line
                buffR = ebGr%xyz(ebGr%is+1,ebGr%js,ebGr%ks,XDIR)

                ! Setting Rin to be the ionosphere
                DomR(1) = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880
            else
                !Take Rin from sunward line
                DomR(1) = ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,XDIR)
            endif
            !Take Rout from sunward line
            !DomR(2) = ebGr%xyz(ebGr%ie+1,ebGr%js,ebGr%ks,XDIR)
            DomR(2) = ebGr%xyz(ebGr%ie,ebGr%js,ebGr%ks,XDIR)

            !Find min/max r and phi along each line of constant i or j
            allocate(locAux%rMin(ebGr%Nip+1))
            allocate(locAux%rMax(ebGr%Nip+1))
            allocate(locAux%pMin(ebGr%Njp+1))
            allocate(locAux%pMax(ebGr%Njp+1))
            do i=1,ebGr%Nip+1
                locAux%rMin(i) = minval(norm2(ebGr%xyz(i,ebGr%js:ebGr%je,ebGr%ks,:),DIM=2))
                locAux%rMax(i) = maxval(norm2(ebGr%xyz(i,ebGr%js:ebGr%je,ebGr%ks,:),DIM=2))
            enddo
            do j=1,ebGr%Njp+1
                locAux%pMin(j) = minval(atan2(ebGr%xyz(ebGr%is:ebGr%ie,j,ebGr%ks,YDIR),ebGr%xyz(ebGr%is:ebGr%ie,j,ebGr%ks,XDIR)))
                locAux%pMax(j) = maxval(atan2(ebGr%xyz(ebGr%is:ebGr%ie,j,ebGr%ks,YDIR),ebGr%xyz(ebGr%is:ebGr%ie,j,ebGr%ks,XDIR)))
            enddo

            if (ebGr%GrID == EGGGRID) then
                write(*,*) 'Initializing EGG locator'
                locate=>Loc_Egg
                locAux%dPhi = PI/ebGr%Njp !Constant j spacing
                locAux%dTh  = 2*PI/ebGr%Nkp !Angular spacing about x-axis
                !Calculate i-interface centered radii
                allocate(locAux%rrI(ebGr%Nip+1,ebGr%Njp))
                do j=1,ebGr%Njp
                    do i=1,ebGr%Nip
                        !Use average of two corner points for this face
                        !Know that ebGr%ks is first (z=0) slice
                        xJ  = ebGr%xyz(i,j  ,ebGr%ks,:)
                        xJp = ebGr%xyz(i,j+1,ebGr%ks,:)
                        locAux%rrI(i,j) = 0.5*(norm2(xJ)+norm2(xJp))
                    enddo
                enddo
                locAux%isInit = .true.
            else
                write(*,*) 'Initializing LFM locator'
                locate=>Loc_LFM
                locAux%dTh  = 2*PI/ebGr%Nkp !Angular spacing about x-axis
                !Get cell centers of 2D grid
                allocate(locAux%xxC(ebGr%Nip,ebGr%Njp))
                allocate(locAux%yyC(ebGr%Nip,ebGr%Njp))
                do j=1,ebGr%Njp
                    do i=1,ebGr%Nip
                        locAux%xxC(i,j) = 0.25*( sum(ebGr%xyz(i:i+1,j:j+1,ebGr%ks,XDIR)) )
                        locAux%yyC(i,j) = 0.25*( sum(ebGr%xyz(i:i+1,j:j+1,ebGr%ks,YDIR)) )
                    enddo
                enddo   
                locAux%isInit = .true.             
            endif

        case(CARTGRID)
            write(*,*) 'Cartesian not implemented'
            stop
        case(SPHGRID)
            !Take Rin/Rout
            DomR(1) = norm2(ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,XDIR:ZDIR))
            DomR(2) = norm2(ebGr%xyz(ebGr%ie-1,ebGr%js,ebGr%ks,XDIR:ZDIR))
            write(*,*) 'Initializing SPH locator'
            locate=>Loc_SPH
            locAux%isInit = .true.
        end select

        !Set inDomain function here
        call inpXML%Set_Val(DomR(1),'domain/rmin',DomR(1))
        call inpXML%Set_Val(DomR(2),'domain/rmax',DomR(2))
        DomE = [DomR(2),-100.0_rp,40.0_rp]        
        call inpXML%Set_Val(domStr,'domain/dtype',"SPH")

        if (Model%isMAGE) domStr = "MAGE"

        select case (trim(toUpper(domStr)))
            case("SPH")
                write(*,*) 'Using spherical inDomain'
                inDomain=>inDomain_Sph
            case("LFM","LFMCYL")
                write(*,*) 'Using LFM inDomain'
                inDomain=>inDomain_LFM
                !Set bounds for yz-cylinder grid
                call inpXML%Set_Val(DomE(1),'domain/xSun' ,30.0_rp)
                call inpXML%Set_Val(DomE(2),'domain/xTail',-100.0_rp)
                call inpXML%Set_Val(DomE(3),'domain/yzMax',40.0_rp)
            case("MAGE")
                write(*,*) 'Using MAGE inDomain'
                inDomain=>inDomain_MAGE

                call inpXML%Set_Val(DomE(1),'domain/xSun' ,  20.0_rp)
                call inpXML%Set_Val(DomE(2),'domain/xTail', -50.0_rp)
                call inpXML%Set_Val(DomE(3),'domain/yzMax',  40.0_rp)

            case("EGG")
                write(*,*) 'Using EGG inDomain'
                inDomain=>inDomain_Egg
            case("ELL")
                write(*,*) 'Using ellipse inDomain'
                !inDomain=>inDomain_Ell
                call inpXML%Set_Val(DomE(1),'domain/xSun' ,DomE(1))
                call inpXML%Set_Val(DomE(2),'domain/xTail',DomE(2))
                call inpXML%Set_Val(DomE(3),'domain/yzMax',DomE(3))
                write(*,*) 'Elliptical inDom not yet implemented ...'
                stop

        end select
        
        ! Gap region options
        gapDefStr = trim(toUpper(gapDefStr))
        call inpXML%Set_Val(gapStr,'domain/gtype',gapDefStr)
        select case (trim(toUpper(gapStr)))
            case("NONE")
                write(*,*) 'Not using inGap'
                inGap=>inGap_Default
            case("SPH","LFM")
                write(*,*) 'Using spherical inGap'
                inGap=>inGap_Sph
                call inpXML%Set_Val(buffR,'domain/rgap',buffR)

                ! Check to make sure gap region at leaset extends to MHD grid
                if (buffR < ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,XDIR)) then
                    write(*,*) 'Gap region must extend to inner boundary of MHD grid'
                    stop
                endif

        end select
    end subroutine InitLoc

!---------------------------------------
!Locator functions

    !3D localization routine for egg grid
    subroutine Loc_Egg(xyz,ijk,Model,ebGr,isInO,ijkO)
        real(rp), intent(in) :: xyz(NDIM)
        integer, intent(out) :: ijk(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        logical, intent(out), optional :: isInO
        integer, intent(in), optional :: ijkO(NDIM)

        logical :: isIn
        real(rp) :: lfmC(NDIM)
        integer :: i0,j0,k0

        ijk = 0
        !Always check is in first
        isIn = inDomain(xyz,Model,ebGr)
        if (present(isInO)) isInO = isIn

        if (.not. isIn) then
            return
        endif

        !Calculate lfmCoords
        lfmC = lfmCoords(xyz)

        !Get k0/j0
        k0 = min(floor(lfmC(KDIR)/locAux%dTh) +1,ebGr%Nkp) !Evenly spaced k
        j0 = min(floor(lfmC(JDIR)/locAux%dPhi)+1,ebGr%Njp) !Evenly spaced j

        ijk(KDIR) = k0
        ijk(JDIR) = j0

        !Test optional guess if present for i0
        if (present(ijkO)) then
            if (norm2(1.0*ijkO)>0) then !Trap for unset ijk
                i0 = ijkO(IDIR)
                if ( lfmC(IDIR) >= locAux%rrI(i0,j0) .and. lfmC(IDIR) <= locAux%rrI(i0+1,j0) ) then
                    !Got it!
                    ijk(IDIR) = i0
                    return
                endif
            endif
        endif !present ijkO

        !If still here then guess failed
        ijk(IDIR) = maxloc( locAux%rrI(:,j0) ,dim=1,mask=lfmC(IDIR) >= locAux%rrI(:,j0))

        !write(*,*) 'Mapped xyz -> ijk ', xyz,ijk

    end subroutine Loc_Egg

    !3D localization routine for egg grid
    subroutine Loc_LFM(xyz,ijk,Model,ebGr,isInO,ijkO)
        real(rp), intent(in) :: xyz(NDIM)
        integer, intent(out) :: ijk(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        logical, intent(out), optional :: isInO
        integer, intent(in), optional :: ijkO(NDIM)

        logical :: isIn,inCell
        real(rp) :: xs,ys,lfmC(NDIM)
        real(rp) :: xp(2)
        integer :: i,j,i0,j0,k0,n
        integer :: i1,i2,j1,j2

        ijk = 0
        !Always check is in first
        isIn = inDomain(xyz,Model,ebGr)
        if (present(isInO)) isInO = isIn


        if (.not. isIn) then
            return
        endif

        !Calculate lfmCoords
        lfmC = lfmCoords(xyz)

        xs = xyz(XDIR)
        ys = sqrt(xyz(YDIR)**2.0 + xyz(ZDIR)**2.0)

        !Get k0/j0
        k0 = min(floor(lfmC(KDIR)/locAux%dTh) +1,ebGr%Nkp) !Evenly spaced k
        ijk(KDIR) = k0

        xp = [xs,ys]
        !Use provided guess if present,
        if (present(ijkO)) then
            !Do snake search around guess
            do n=1,NSnake
                i = ijkO(IDIR) + dsI(n)
                j = ijkO(JDIR) + dsJ(n)
                !Check if point is valid
                isIn = (i>=ebGr%is) .and. (i<=ebGr%ie) .and. (j>=ebGr%js) .and. (j<=ebGr%je)

                if (.not. isIn) cycle
                !Otherwise check if this is the right cell
                inCell = CheckIJ(xp,[i,j],Model,ebGr)
                if (inCell) then !Found it, let's get the hell out of here
                    ijk(IDIR:JDIR) = [i,j]
                    return
                endif
            enddo

        endif !Using guess
        

        !If we're still here, do this the hard way
        !Cut out obviously incorrect 2D indices
        call lfmChop(Model,ebGr,[xs,ys],i1,i2,j1,j2)

        !If still here, just pick the closest one
        ijk(IDIR:JDIR) = minloc( (locAux%xxC(i1:i2,j1:j2)-xs)**2.0 + (locAux%yyC(i1:i2,j1:j2)-ys)**2.0 )
        ijk(IDIR:JDIR) = ijk(IDIR:JDIR) + [i1-1,j1-1] !Correct for offset
        if (ijk(KDIR)<0) then
            write(*,*) 'xyz/ijk = ',xyz,ijk
        endif

    end subroutine Loc_LFM

    !Check whether 2D point xy is in LFM cell ijG
    function CheckIJ(xy,ijG,Model,ebGr) result(isIn)
        real(rp), intent(in) :: xy (NDIM-1)
        integer, intent(in)  :: ijG(NDIM-1)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr

        logical :: isIn

        real(rp) :: xCs(4,2)

        xCs(1,:) = ebGr%xyz(ijG(IDIR)  ,ijG(JDIR)  ,ebGr%ks,XDIR:YDIR)
        xCs(2,:) = ebGr%xyz(ijG(IDIR)+1,ijG(JDIR)  ,ebGr%ks,XDIR:YDIR)
        xCs(3,:) = ebGr%xyz(ijG(IDIR)+1,ijG(JDIR)+1,ebGr%ks,XDIR:YDIR)
        xCs(4,:) = ebGr%xyz(ijG(IDIR),  ijG(JDIR)+1,ebGr%ks,XDIR:YDIR)

        !Test guess
        isIn = inCell2D(xy,xCs)
     
    end function CheckIJ

    !3D localization routine for spherical grid
    subroutine Loc_SPH(xyz,ijk,Model,ebGr,isInO,ijkO)
        real(rp), intent(in) :: xyz(NDIM)
        integer, intent(out) :: ijk(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        logical, intent(out), optional :: isInO
        integer, intent(in), optional :: ijkO(NDIM)


        logical :: isIn
        real(rp) :: helioC(NDIM)
        !For evenly spaced grid in three dimensions
        real(rp) :: dphi,dtheta,dr, thetaMin, rMin
        integer :: i0,j0,k0

        !Add localization routine here
        ijk = 0
        isInO = .false.
  
        !Always check is in first
        isIn = inDomain(xyz,Model,ebGr)
        if (present(isInO)) isInO = isIn

        if (.not. isIn) then
            return
        endif

        !Calculate helioCoords
        helioC = helioCoords(xyz)

        !Even spacing in k
        dphi = 2*PI/ebGr%Nkp
        !dtheta = theta2 - theta1 for even spacing in j
        dtheta = acos(ebGr%xyz(ebGr%is,ebGr%js+1,ebGr%ks,ZDIR)/norm2(ebGr%xyz(ebGr%is,ebGr%js+1,ebGr%ks,:))) &
               - acos(ebGr%xyz(ebGr%is,ebGr%js,  ebGr%ks,ZDIR)/norm2(ebGr%xyz(ebGr%is,ebGr%js,  ebGr%ks,:)))
        !dr = r2-r1 for even spacing in r
        dr = norm2(ebGr%xyz(ebGr%is+1,ebGr%js,ebGr%ks,:)) - norm2(ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,:))
        
        thetaMin = acos(ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,ZDIR)/norm2(ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,:)))
        rMin = norm2(ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,:)) 
 
        ! pick the closest one
        i0 = min(floor((helioC(IDIR)-rMin)/dr) + 1,ebGr%Nip) !Evenly spaced i
        j0 = min(floor((helioC(JDIR)-thetaMin)/dtheta) + 1,ebGr%Njp) !Evenly spaced j
        k0 = min(floor(helioC(KDIR)/dphi) + 1,ebGr%Nkp) !Evenly spaced k

        ijk(IDIR) = i0
        ijk(JDIR) = j0
        ijk(KDIR) = k0
        
        !write(*,*) 'Mapped xyz -> HelioC ', xyz, helioC
        !write(*,*) 'Mapped xyz -> ijk ', xyz, ijk

    end subroutine Loc_SPH
!---------------------------------------
!inDomain functions

    !Simple spherical inDomain
    function inDomain_Sph(xyz,Model,ebGr) result(inDom)
        real(rp), intent(in) :: xyz(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr

        real(rp) :: r
        logical :: inDom
        r = norm2(xyz)
        if (r <= DomR(1) .or. r >= DomR(2)) then
            inDom = .false.
            
        else
            inDom = .true.
        endif

    end function inDomain_Sph

    function inDomain_MAGE(xyz,Model,ebGr) result(inDom)
        real(rp), intent(in) :: xyz(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr

        real(rp) :: r,xs,ys,xSun
        integer :: ij(2), ix
        integer :: i1,i2,j1,j2

        logical :: inDom

        xSun = ebGr%xyz(ebGr%ie,ebGr%js,ebGr%ks,XDIR) !Sunward line

        !Start w/ some silly short circuits
        r = norm2(xyz)
        if (r <= DomR(1)) then
            !Inner boundary
            inDom = .false.
            return
        endif

        if (r <= xSun) then
            !Within sphere bounded by dayside, you're good
            inDom = .true.
            return
        endif

        !Do some more short circuits but work in upper half plane coords
        xs = xyz(XDIR)
        ys = sqrt(xyz(YDIR)**2.0 + xyz(ZDIR)**2.0)

        !Use a nice rectangular bound
        if ( (xs <=0) .and. (xs >= -225) .and. (ys <= 80) ) then
            inDom = .true.
            return
        endif

        !All right, we've had a good time but if you're still here we're gonna have to do this the hard way
        
        !Cut out obviously incorrect 2D indices
        call lfmChop(Model,ebGr,[xs,ys],i1,i2,j1,j2)

        !If still here, just pick the closest one

        ij = minloc( (locAux%xxC(i1:i2,j1:j2)-xs)**2.0 + (locAux%yyC(i1:i2,j1:j2)-ys)**2.0 )
        ij = ij + [i1-1,j1-1] !Correct for offset
        ix = ij(IDIR)
        if ( abs(ix - ebGr%ie) <= 4 ) then
            !Closest point is pretty damn close to the outer boundary, bag that
            inDom = .false.
            return
        else
            inDom = .true.
            return
        endif

    end function inDomain_MAGE

    function inDomain_LFM(xyz,Model,ebGr) result(inDom)
        real(rp), intent(in) :: xyz(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr

        real(rp) :: r,x,yzR
        logical :: inDom


        inDom = .true.
        r = norm2(xyz)

        !Inner boundary
        if (r <= DomR(1)) then
            inDom = .false.
            return
        endif

        !X-SM bounds
        x = xyz(XDIR)
        if ( (x>=DomE(1)) .or. (x<=DomE(2)) ) then
            inDom = .false.
            return
        endif

        !Y-Z cylinder bounds
        yzR = sqrt( xyz(YDIR)**2.0 + xyz(ZDIR)**2.0 )
        if (yzR >= DomE(3)) then
            inDom = .false.
            return
        endif

    end function inDomain_LFM

    function inDomain_Egg(xyz,Model,ebGr) result(inDom)
        real(rp), intent(in) :: xyz(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        logical :: inDom

        real(rp) :: r,rOutJ,rOutJP
        real(rp) :: lfmC(NDIM)
        integer :: j0

        inDom = .true.
        r = norm2(xyz)
        !Start w/ short circuit on radius
        if ( (r>DomR(1)) .and. (r < DomR(2)) ) return
        if (r <= DomR(1)) then
            inDom = .false.
            return
        endif

        !If we're still here, then we have to do more work
        lfmC = lfmCoords(xyz)
        !For egg, phi doesn't depend on r
        !Find j cell
        j0 = maxloc( locAux%pMin,dim=1,mask=lfmC(JDIR)>=locAux%pMin )
        !For this j cell, check that 2d radius is inside
        !Get min of both interface radii
        rOutJ  = norm2(ebGr%xyz(ebGr%ie,j0+0,ebGr%ks,:))
        rOutJP = norm2(ebGr%xyz(ebGr%ie,j0+1,ebGr%ks,:))

        if (lfmC(IDIR) >= min(rOutJ,rOutJP)) then
            inDom = .false.
        else
            inDom = .true.
        endif

    end function inDomain_Egg

    !Default inGap is always F
    function inGap_Default(xyz,Model,ebGr) result(inDom)
        real(rp), intent(in) :: xyz(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr

        logical  :: inDom

        inDom = .false.

    end function inGap_Default

    function inGap_Sph(xyz,Model,ebGr) result(inDom)
        real(rp), intent(in) :: xyz(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        real(rp) :: r
        logical  :: inDom

        r = norm2(xyz)
        if (r >= DomR(1) .and. r <= buffR) then
            inDom = .true.
        else
            inDom = .false.
        endif

    end function inGap_Sph

    !Lazy hard-wired function to decide if foot-point is closed
    function isClosed(xyz,Model) result(inDom)
        real(rp), intent(in) :: xyz(NDIM)
        type(chmpModel_T), intent(in) :: Model

        real(rp) :: r
        logical :: inDom

        r = norm2(xyz)
        if (r <= rClosed) then
            inDom = .true.
        else
            inDom = .false.
        endif
    end function isClosed

    !Given xy (projected to upper half plane) return index bounds for 2D search space
    subroutine lfmChop(Model,ebGr,xy,i1,i2,j1,j2)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        real(rp), intent(in) :: xy(2)
        integer, intent(out) :: i1,i2,j1,j2

        real(rp) :: lfmC(NDIM)
        real(rp) :: r,phi
        lfmC = lfmCoords([xy(1),xy(2),0.0_rp])

        r = lfmC(IDIR)
        phi = lfmC(JDIR)

        i1 = ebGr%is
        i2 = ebGr%ie
        j1 = ebGr%js
        j2 = ebGr%je

        i1 = maxloc( locAux%rMax,dim=1,mask=locAux%rMax<r )
        i2 = maxloc( locAux%rMin,dim=1,mask=r>locAux%rMin )

        j1 = maxloc( locAux%pMax,dim=1,mask=locAux%pMax<phi)
        j2 = maxloc( locAux%pMin,dim=1,mask=phi>locAux%pMin)

        !Finish up, add extra cell for safety
        i1 = max(i1-1,ebGr%is)
        i2 = min(i2+1,ebGr%ie)
        j1 = max(j1-1,ebGr%js)
        j2 = min(j2+1,ebGr%je)

        !write(*,*) 'Chop down to i1,i2,j1,j2 = ', i1,i2,j1,j2
    end subroutine lfmChop


    !Convert xyz into LFM (xs,ys,ThX) coords
    !rs,ps are polar coords in upper half plane, ThX is angle about x axis
    function lfmCoords(xyz) result(lfmC)
        real(rp), intent(in) :: xyz(NDIM)
        real(rp) :: lfmC(NDIM)

        real(rp) :: xs,ys,ThX
        xs = xyz(XDIR)
        ys = sqrt(xyz(YDIR)**2.0 + xyz(ZDIR)**2.0)

        lfmC(IDIR) = sqrt(xs**2.0+ys**2.0)
        lfmC(JDIR) = atan2(ys,xs) !Has to be >=0 b/c y>=0

        !Calculate angle about x axis
        ThX = atan2(xyz(ZDIR),xyz(YDIR))
        !Force to 0,2pi
        if (ThX<0) ThX=ThX+2*PI
        lfmC(KDIR) = ThX

    end function lfmCoords

    !Convert xyz into Helio spherical (r, theta, phi) coords
    !theta is angle from z axis, phi is angle from x axis around z axis 
    function helioCoords(xyz) result(helioC)
        real(rp), intent(in) :: xyz(NDIM)
        real(rp) :: helioC(NDIM)

        real(rp) :: zs,rxy,phi

        zs = xyz(ZDIR)
        rxy = sqrt(xyz(XDIR)**2.0 + xyz(YDIR)**2.0)

        helioC(IDIR) = sqrt(zs**2.0 + rxy**2.0)
        !Calculate theta [0, pi]
        helioC(JDIR) = acos(zs/sqrt(zs**2.0 + rxy**2.0)) 

        !Calculate phi
        phi = atan2(xyz(YDIR),xyz(XDIR))
        !Force to 0,2pi
        if (phi<0) phi=phi+2*PI
        helioC(KDIR) = phi

    end function helioCoords
end module gridloc
