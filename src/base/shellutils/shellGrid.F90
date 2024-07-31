!> Various data structures and routines to define grids on spherical shells
module shellGrid
    use kdefs
    use math

    implicit none

    !> Identifiers for location of variable data relative to shell grid
    !> Cell center, corners, theta faces, phi faces
    !> SHGR_FACE_THETA is face with normal pointing in theta direction
    !> SHGR_FACE_PHI is face with normal pointing in phi direction
    enum, bind(C)
        enumerator :: SHGR_CC,SHGR_CORNER,SHGR_FACE_THETA,SHGR_FACE_PHI
    end enum
    
    !> Data type for holding 2D spherical shell grid
    type ShellGrid_T
        
        character(len=strLen) :: name
            !! Name assigned to this ShellGrid instance, determined by the model initializing it
        real(rp) :: radius
            !! Radius that this ShellGrid lives at, in units of Rp (planetary radii)
        integer :: Nt,Np
            !! Number of colat/lon cells (theta, phi)
        real(rp), dimension(:), allocatable :: th, ph, lat
            !! (Nt+1 or Np+1) [radians] grid corners
            !! th (theta) is colatitude and runs from north pole toward south
            !! Phi is longitude, with zero/2pi at 12 MLT
            !! Assuming lat in -pi/2,pi/2 and lon in [0,2pi]
        real(rp), dimension(:), allocatable :: thc, phc, latc  
            !! (Nt or Np) [radians] grid centers
        logical :: doSP = .false., doNP = .false. 
            !! Whether active grid contains south/north pole, no ghosts in this case
        logical :: ghostSP = .false., ghostNP = .false. 
            !! Whether the last ghost corner is the south/north pole

        real(rp) :: minTheta, maxTheta, minPhi, maxPhi
            !! Theta and phi bounds of grid excluding ghost cells
        real(rp) :: minGTheta, maxGTheta, minGPhi, maxGPhi
            !! Theta and phi bounds of grid including ghost cells

        integer :: is=1,ie,js=1,je
            !! Local indices, active region
        integer :: isg,ieg,jsg,jeg
            !! Local indices, including ghosts

        integer :: Ngn=0, Ngs=0, Nge=1, Ngw=1
            !! Ghosts for north, south, east, and west boundaries. East=>larger phi. West => smaller phi.
            !! default to 0 for N/S, 1 for E/W

        logical :: isPeriodic   
            !! Whether the low/high phi boundary is periodic
        logical :: isPhiUniform 
            !! Define this to speed up search for interpolation

        !> Subgrid information
        !> ShellGrids that are subgrids of other shellGrids store info about their parent grid
        logical :: isChild = .false.
        character(len=strLen) :: parentName = "N/A"
            !! Name of the parent grid  that this one derives from
        integer :: bndis,bndie,bndjs,bndje
            !! Indices of parent grid that bound this grid

    end type ShellGrid_T

    type ShellGridVar_T

        integer :: loc
            !! Location of data on the shellGrid (e.g. center, corner, theta of phi face)
            !! Corresponds to enum above (SHGR_CC, SHGR_CORNER, SHGR_FACE_THETA, SHGR_FACE_PHI)
        integer :: Ni, Nj
            !! Number of values in i and j direction
        integer :: isv,iev,jsv,jev
            !! Start and end indices for this variable
            !! ex: if loc=SHGR_CORNER, isv = sh%isg, iev=sh%ieg+1
            !! This is helpful for e.g. InterpShellVar_TSC_pnt determining size of dtheta and dPhi arrays
        real(rp), dimension(:,:), allocatable :: data
            !! The actual variable values
        logical, dimension(:,:), allocatable :: mask
            !! Mask indicating whether the data at a given index is valid
            !! e.g. good for interpolation, etc.


        ! Commenting out for now, I think we will ultimately not use this
        !logical, dimension(4) :: bcsApplied 
            !! Flag indicating whether BCs were applied (ghosts filled) for [n,s,e,w] boundaries

    end type ShellGridVar_T


    contains

    !> Create a shell grid data structure
    !> Takes Theta and Phi 1D arrays (uniform or not)
    !> Decides if isPeriodic and isPhiUniform based on the Phi array passed in
    subroutine GenShellGrid(shGr,Theta,Phi,name,nGhosts,radO)
        type(ShellGrid_T), intent(inout) :: shGr
        real(rp), dimension(:), intent(in) :: Theta, Phi
        character(len=*) :: name
            !! Name identifier used for this grid instance
        integer, optional, dimension(4), intent(in) :: nGhosts 
            !! How many ghosts on each side (n,s,e,w)
        real(rp), optional, intent(in) :: radO
            !! Radius in planetary radii from planet center this grid lives at
            !! WARNING: Will default to 1 if not provided! That should rarely be the case
        
        integer :: i,j
        real(rp) :: delta
        real(rp), dimension(:), allocatable :: dphi

        ! Parse optional parameters
        if (present(nGhosts)) then
            shGr%Ngn = nGhosts(NORTH) ! otherwise, always 0 as set in ShellGrid type
            shGr%Ngs = nGhosts(SOUTH) ! otherwise, always 0 as set in ShellGrid type
            shGr%Nge = nGhosts(EAST) ! otherwise, always 1 as set in ShellGrid type
            shGr%Ngw = nGhosts(WEST) ! otherwise, always 1 as set in ShellGrid type
        end if

        ! do some checks first
        if (.not.(isAscending(Theta))) then 
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Theta array must be ascending. Quitting..."
            stop
        end if

        if (.not.(isAscending(Phi))) then 
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Phi array must be ascending. Quitting..."
            stop
        end if

        if (any(Theta<0.).or.(any(Theta>PI))) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Theta array should be in the range [0,PI]. Quitting..."
            stop
        end if

        if (any(Phi<0.).or.(any(Phi>2*PI))) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Phi array should be in the range [0,2*PI]. Quitting..."
            stop
        end if

        ! decide if the grid is periodic
        if ( ( minval(Phi) <= TINY ).and.( maxval(Phi) >= 2*PI - TINY) ) then
            shGr%isPeriodic = .True.
        else
            shGr%isPeriodic = .False.
        end if
        
        if ( (shGr%isPeriodic).and.(shGr%Nge/=shGr%Ngw) ) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Periodic grid must have the same number of ghosts on low/high phi boundaries. Quitting..."
            stop
        end if

        if ( (.not.(shGr%isPeriodic)).and.(shGr%Nge/=shGr%Ngw) ) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Grids with Nge/=Ngw have not been implemented. Quitting..."
            stop
        end if
        
        ! Decide if this has north/south pole
        shGr%doNP = .false.
        shGr%doSP = .false.
        if ( maxval(Theta) >= (PI - TINY) ) then
            shGr%doSP = .true.
        endif

        if ( minval(Theta) <= (TINY) ) then
            shGr%doNP = .true.
        endif

        if ( (shGr%doNP).and.(shGr%Ngn/=0) ) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Grid containing the north pole can't have Ngn/=0. Quitting..."
            stop
        end if

        if ( (shGr%doSP).and.(shGr%Ngs/=0) ) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Grid containing the south pole can't have Ngs/=0. Quitting..."
            stop
        end if

        ! define various array indices
        shGr%Np = size(Phi)-1
        shGr%Nt = size(Theta)-1

        shGr%is = 1; shGr%ie = shGr%Nt
        shGr%js = 1; shGr%je = shGr%Np

        shGr%isg = shGr%is - shGr%Ngn
        shGr%ieg = shGr%ie + shGr%Ngs
        shGr%jsg = shGr%js - shGr%Ngw
        shGr%jeg = shGr%je + shGr%Nge

        ! decide if the grid is uniform in Phi
        ! helps to speed up interpolation search
        if (allocated(dphi)) deallocate(dphi)
        allocate(dphi(shGr%Np))   ! helper
        dphi = Phi(2:shGr%Np+1) - Phi(1:shGr%Np)
        if ( all( dphi - dphi(1) <= TINY ) ) then
            shGr%isPhiUniform = .True.
        else 
            shGr%isPhiUniform = .False.
        end if
        deallocate(dphi)

        associate(is=>shGr%is, ie=>shGr%ie, js=>shGr%js, je=>shGr%je, &
            isg=>shGr%isg, ieg=>shGr%ieg, jsg=>shGr%jsg, jeg=>shGr%jeg)

        ! Nuke arrays if already allocated
        if (allocated(shGr%th))   deallocate(shGr%th)
        if (allocated(shGr%ph))   deallocate(shGr%ph)
        if (allocated(shGr%lat))  deallocate(shGr%lat)
        if (allocated(shGr%thc))  deallocate(shGr%thc)
        if (allocated(shGr%phc))  deallocate(shGr%phc)
        if (allocated(shGr%latc)) deallocate(shGr%latc)
        ! Create new arrays
        allocate(shGr%th  (isg:ieg+1))
        allocate(shGr%thc (isg:ieg  ))
        allocate(shGr%ph  (jsg:jeg+1))
        allocate(shGr%phc (jsg:jeg  ))
        allocate(shGr%lat (isg:ieg+1))
        allocate(shGr%latc(isg:ieg  ))

        ! Set grid coordinates
        shGr%th(is:ie+1) = Theta  ! note the arrays are conformable because of the index definitions above
        shGr%ph(js:je+1) = Phi

        ! Define ghost coordinates

        ! Do north, unless it's the pole (no ghosts for poles)
        if ( (.not.shGr%doNP).and.(shGr%Ngn/=0) ) then
            ! linearly extrapolate
            do i=1,shGr%Ngn
                shGr%th(is-i) = 2*shGr%th(is) - shGr%th(is+i)
            end do

            ! check if we got into negative thetas
            ! then just fill the space with same size cells down to theta=0
            if ( shGr%th(is-shGr%Ngn) < 0 ) then
                delta = shGr%th(is)/shGr%Ngn

                do i=1,shGr%Ngn
                    shGr%th(is-i) = shGr%th(is) - delta*i
                end do

                shGr%ghostNP = .true.
            end if
        end if

        ! Do south, unless it's the pole (no ghosts for poles)
        if ( (.not.shGr%doSP).and.(shGr%Ngs/=0) ) then
            ! linearly extrapolate
            do i=1,shGr%Ngs
                shGr%th(ie+1+i) = 2*shGr%th(ie+1) - shGr%th(ie+1-i)
            end do

            ! check if we got into thetas > pi
            ! then just fill the space with same size cells up to theta=pi
            if ( shGr%th(ie+1+shGr%Ngs) > PI ) then
                delta = (PI - shGr%th(ie+1))/shGr%Ngs
                do i=1,shGr%Ngs
                    shGr%th(ie+1+i) = shGr%th(ie+1) + delta*i
                end do

                shGr%ghostSP = .true.
            end if
        end if

        ! if non-periodic grid is needed, can implement ghosts on E/W ends similarly to N/S above
        ! but we assume the grid is always periodic
        if ( (.not.shGr%isPeriodic) )  then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Non-periodic grids are not implemented. Quitting..."
            stop
        endif

        ! Note, we already made sure that for a periodic grid Nge=Ngw, phi(js)=0 and phi(je+1)=pi
        do j=1,shGr%Nge
            shGr%ph(js-j)   = shGr%ph(je+1-j) - 2*PI
            shGr%ph(je+1+j) = shGr%ph(js+j) + 2*PI
        end do

        ! Set centers
        shGr%thc = ( shGr%th(isg+1:ieg+1) + shGr%th(isg:ieg) )/2.
        shGr%phc = ( shGr%ph(jsg+1:jeg+1) + shGr%ph(jsg:jeg) )/2.

        ! Set latitude
        shGr%lat  = PI/2.0_rp - shGr%th
        shGr%latc = PI/2.0_rp - shGr%thc

        ! Note, the bounds below only include the active cells
        shGr%minTheta = minval(shGr%th(is:ie+1))
        shGr%maxTheta = maxval(shGr%th(is:ie+1))
        shGr%minPhi   = minval(shGr%ph(js:je+1))
        shGr%maxPhi   = maxval(shGr%ph(js:je+1))

        ! this includes ghosts
        shGr%minGTheta = minval(shGr%th)
        shGr%maxGTheta = maxval(shGr%th)
        shGr%minGPhi   = minval(shGr%ph)
        shGr%maxGPhi   = maxval(shGr%ph)

        end associate

        shGr%name = name

        if (present(radO)) then
            ! Try to trap for people thinking units are meters or km
            if (radO > 10) then
                write(*,*) "WARNING for ShellGrid with name: ",name
                write(*,*) "Radius being set to ",radO," planetary radii. That seems kinda big."
            endif
            shGr%radius = radO
        else
            write(*,*) "WARNING for ShellGrid with name: ",name
            write(*,*) "No radius provided, assuming Earth's ionosphere according to kdefs:"
            shGr%radius = RIonE*1.0e6/REarth
            write(*,*) shGr%radius," Rp"
        endif

    end subroutine GenShellGrid


    subroutine initShellVar(shGr, loc, shellVar, maskO)
        !! Inits a ShellGridVar that associated with provided and initialized ShellGrid
        type(ShellGrid_T), intent(in) :: shGr
            !! ShellGrid that this variable is related to
        integer, intent(in) :: loc
            !! Location of data (cell center, corner, theta or phi face)
        type(ShellGridVar_T), intent(out) :: shellVar
        logical, dimension(:,:), optional, intent(in) :: maskO
            !! Optional mask to initialize with

        integer :: iExtra, jExtra
        
        ! If you didn't want your data blown up you shouldn't have called init
        if (allocated(shellVar%data)) deallocate(shellVar%data)
        if (allocated(shellVar%data)) deallocate(shellVar%mask)
            
        shellVar%loc = loc

        ! Determine which dimensions have extra index relative to # cells based on variable's location on grid
        select case(loc)
            case(SHGR_CC)
                iExtra = 0
                jExtra = 0
            case(SHGR_CORNER)
                iExtra = 1
                jExtra = 1
            case(SHGR_FACE_THETA)
                iExtra = 1
                jExtra = 0
            case(SHGR_FACE_PHI)
                iExtra = 0
                jExtra = 1
            case default
                write(*,*) "initShellGridVar got an invalid data location:",loc
                stop
        end select

        allocate(shellVar%data(shGr%isg:shGr%ieg+iExtra, shGr%jsg:shGr%jeg+jExtra))
        allocate(shellVar%mask(shGr%isg:shGr%ieg+iExtra, shGr%jsg:shGr%jeg+jExtra))
        shellVar%Ni = shGr%Nt + shGr%Ngn + shGr%Ngs + iExtra
        shellVar%Nj = shGr%Np + shGr%Nge + shGr%Ngw + jExtra

        shellVar%isv = shGr%isg
        shellVar%iev = shGr%ieg + iExtra
        shellVar%jsv = shGr%jsg
        shellVar%jev = shGr%jeg + jExtra


        shellVar%data = 0.  ! initialize to 0

        ! Init mask, either by maskO or defaults
        if (present(maskO)) then
            if ( all(shape(shellVar%mask) == shape(maskO)) ) then
                shellVar%mask = maskO
            else
                write(*,*)"ERROR in initShellVar: maskO shape doesn't match."
                stop
            endif
        else
            shellVar%mask = .false.  ! Up to user to determine which points are valid
        endif
        
    end subroutine initShellVar


    subroutine GenChildShellGrid(pSG, cSG, name, nGhosts, sub_is, sub_ie, sub_js, sub_je)
        !! Given a parent ShellGrid, makes a child ShellGrid as a subset of parent
        type(ShellGrid_T), intent(in) :: pSG
            !! Parent ShellGrid
        type(ShellGrid_T), intent(out) :: cSG
            !! Child ShellGrid
        character(len=*) :: name
            !! Name identifier used for this grid instance
        integer, optional, dimension(4), intent(in) :: nGhosts
            !! Number of ghosts the child grid will have
        integer, optional, intent(in) :: sub_is, sub_ie, sub_js, sub_je
            !! Start and end i/j indices of parent grid that bound active domain of new child grid
            !! These are optional. If left out, we will essentially make a copy of the parent grid

        integer :: is, ie, js, je
            !! Actual bounds used

        ! If a bound is provided then we use that, if not we default to parent grid's bounds
        is = whichInd(sub_is, pSG%is,   present(sub_is))
        ie = whichInd(sub_ie, pSG%ie+1, present(sub_ie))
        js = whichInd(sub_js, pSG%js  , present(sub_js))
        je = whichInd(sub_je, pSG%je+1, present(sub_je))

        ! Check for valid bounds
        if (is < pSG%is) then
            write(*,*) "ERROR GenChildShellGrid: Invalid is bound."
            write(*,*) "Requested:",is
            write(*,*) "Mimumum:",pSG%is
            stop
        endif
        if (ie > pSG%ie+1) then
            write(*,*) "ERROR GenChildShellGrid: Invalid ie bound."
            write(*,*) "Requested:",ie
            write(*,*) "Maximum:",pSG%ie+1
            stop
        endif
        if (js < pSG%js) then
            write(*,*) "ERROR GenChildShellGrid: Invalid js bound."
            write(*,*) "Requested:",js
            write(*,*) "Mimumum:",pSG%js
            stop
        endif
        if (je > pSG%je+1) then
            write(*,*) "ERROR GenChildShellGrid: Invalid je bound."
            write(*,*) "Requested:",je
            write(*,*) "Maximum:",pSG%je+1
            stop
        endif
        if (is > ie) then
            write(*,*) "ERROR GenChildShellGrid: is > ie"
            stop
        endif
        if (js > je) then
            write(*,*) "ERROR GenChildShellGrid: js > je"
            stop
        endif

        ! Make sure name is not a copy of its parent's
        if (trim(name) == trim(pSG%name)) then
            write(*,*) "ERROR GenChildShellGrid: Child can't have same name as its parent"
            write(*,*) " Try '",trim(pSG%name)," ",trim(pSG%name),"sson'"
            stop
        endif

        ! Otherwise we are ready to make a child grid
        if ( present(nGhosts) ) then
            ! The ghosts can't overrun parent grid's ghost bounds
            if (    is - nGhosts(NORTH) < pSG%isg  ) then
                write(*,*) "ERROR GenChildShellGrid: Child ghosts overrun parent ghosts in is direction"
                stop
            elseif (ie + nGhosts(SOUTH) > pSG%ieg+1) then
                write(*,*) "ERROR GenChildShellGrid: Child ghosts overrun parent ghosts in ie direction"
                stop
            elseif (js - nGhosts(WEST ) < pSG%jsg  ) then
                write(*,*) "ERROR GenChildShellGrid: Child ghosts overrun parent ghosts in js direction"
                stop
            elseif (je + nGhosts(EAST ) > pSG%jeg+1) then
                write(*,*) "ERROR GenChildShellGrid: Child ghosts overrun parent ghosts in je direction"
                stop
            endif

            ! Otherwise, this ghost definition is okay
            call GenShellGrid(cSG, pSG%th(is:ie), pSG%ph(js:je), name, nGhosts=nGhosts, radO=pSG%radius)

        else
            ! No ghosts defined, go with default
            call GenShellGrid(cSG, pSG%th(is:ie), pSG%ph(js:je), name, radO=pSG%radius)
        endif
        

        cSG%isChild = .true.
        cSG%parentName = pSG%name
        cSG%bndis = is
        cSG%bndie = ie
        cSG%bndjs = js
        cSG%bndje = je

        contains 

        function whichInd(indyes,indno,ispresent) result(indchoice)
            integer, intent (in) :: indyes, indno
            logical, intent (in) :: ispresent
            integer              :: indchoice
        
            if (ispresent) then
                indchoice = indyes
            else    
                indchoice = indno
            end if  

            return
        end function
    end subroutine GenChildShellGrid


    subroutine deallocShellGrid(sh)
        !! Deallocates any allocated memory
        type(ShellGrid_T), intent(inout) :: sh

        deallocate(sh%th)
        deallocate(sh%ph)
        deallocate(sh%thc)
        deallocate(sh%phc)
        deallocate(sh%lat)
        deallocate(sh%latc)

    end subroutine deallocShellGrid

end module shellGrid
