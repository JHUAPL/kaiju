!> Various data structures and routines to define grids on spherical shells
module shellGrid
    use kdefs
    use math
    use planethelper

    implicit none

    !> Identifiers for location of variable data relative to shell grid
    !> Cell center, corners, theta faces, phi faces
    !> SHGR_FACE_THETA is face with normal pointing in theta direction
    !> SHGR_FACE_PHI is face with normal pointing in phi direction
    enum, bind(C)
        enumerator :: SHGR_CC,SHGR_CORNER,SHGR_FACE_THETA,SHGR_FACE_PHI
    end enum

    character(len=strLen), parameter :: def_noParentName ="N/A"
    

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
        !! True = good value, False = bad value


    end type ShellGridVar_T


    !> Data type for holding 2D spherical shell grid
    type ShellGrid_T
        
        character(len=strLen) :: name
            !! Name assigned to this ShellGrid instance, determined by the model initializing it
        real(rp) :: radius
            !! [Rp] Radius that this ShellGrid lives at, in units of Rp (planetary radii)
        integer :: Nt,Np
            !! Number of active colat/lon cells (theta, phi)
        real(rp), dimension(:), allocatable :: th, ph, lat
            !! (Nt+1 or Np+1) [radians] grid corners
            !! th (theta) is colatitude and runs from north pole toward south
            !! Phi is longitude, with zero/2pi at 12 MLT
            !! Assuming lat in -pi/2,pi/2 and lon in [0,2pi]
        real(rp), dimension(:), allocatable :: thc, phc, latc  
            !! (Nt or Np) [radians] grid centers
        real(rp), dimension(:), allocatable :: thRp
            !! (Nt+1) [rad] theta corners projected to 1 Rp assuming dipole
        real(rp), dimension(:), allocatable :: thcRp
            !! (Nt) [rad] theta centers projected to 1 Rp assuming dipole
        real(rp), dimension(:,:), allocatable :: areaCC
            !! [Rp^2] Area of each cell
        real(rp), dimension(:,:), allocatable :: lenFace_th
            !! (Nt+1, Np) [Rp] arc length of each Theta face
        real(rp), dimension(:,:), allocatable :: lenFace_ph
            !! (Nt, Np+1) [Rp] arc length of each Phi face
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
        character(len=strLen) :: parentName = def_noParentName
            !! Name of the parent grid that this one derives from
        integer :: bndis,bndie,bndjs,bndje
            !! Indices of parent grid that bound this grid

    end type ShellGrid_T


    contains

    !> Create a shell grid data structure
    !> Takes Theta and Phi 1D arrays (uniform or not)
    !> Decides if isPeriodic and isPhiUniform based on the Phi array passed in
    subroutine GenShellGrid(shGr,Theta,Phi,name,nGhosts,radO,TPIncludesGhostsO)
        type(ShellGrid_T), intent(inout) :: shGr
        real(rp), dimension(:), intent(in) :: Theta, Phi
        character(len=*) :: name
            !! Name identifier used for this grid instance
        integer, optional, dimension(4), intent(in) :: nGhosts 
            !! How many ghosts on each side (n,s,e,w)
        real(rp), optional, intent(in) :: radO
            !! Radius in planetary radii from planet center this grid lives at
            !! WARNING: Will default to Earth's ionospheric altitude if not provided! That should rarely be the case
        logical, optional, intent(in) :: TPIncludesGhostsO
        
        integer :: i,j
        integer :: phiLo, phiHi
        logical :: addGhosts
        real(rp) :: delta, L
        real(rp), dimension(:), allocatable :: dphi
        real(rp), dimension(3) :: xyz

        ! Parse optional parameters
        if (present(nGhosts)) then
            shGr%Ngn = nGhosts(NORTH) ! is ghost (default 0)
            shGr%Ngs = nGhosts(SOUTH) ! ie ghost (default 0)
            shGr%Nge = nGhosts(EAST)  ! js ghost (default 1)
            shGr%Ngw = nGhosts(WEST)  ! je ghost (default 1)
        end if

        if(present(TPIncludesGhostsO)) then
            addGhosts = .not. TPIncludesGhostsO
        else
            addGhosts = .true.
        endif

        ! do some checks first
        if (.not.(isAscending(Theta))) then 
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Theta array must be ascending. Quitting..."
            write(*,*) Theta
            stop
        end if

        if (.not.(isAscending(Phi))) then 
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Phi array must be ascending. Quitting..."
            write(*,*) Phi
            stop
        end if


        ! Decide if the grid is periodic
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

        ! Define various array indices
        ! Np and Nt are for active domain only. If ghosts are already included in Theta/Phi, subtract nghosts from them
        if (addGhosts) then
            shGr%Np = size(Phi)-1
            shGr%Nt = size(Theta)-1
        else
            shGr%Np = size(Phi  ) - shGr%Nge - shGr%Ngw - 1
            shGr%Nt = size(Theta) - shGr%Ngn - shGr%Ngs - 1
        endif

        shGr%is = 1; shGr%ie = shGr%Nt
        shGr%js = 1; shGr%je = shGr%Np

        shGr%isg = shGr%is - shGr%Ngn
        shGr%ieg = shGr%ie + shGr%Ngs
        shGr%jsg = shGr%js - shGr%Nge
        shGr%jeg = shGr%je + shGr%Ngw


        ! Do some more checks now that we have some sizes defined

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

        ! No one allowed beyond Theta=[0,PI], even ghosts
        if (any(Theta<0.).or.(any(Theta>PI))) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Theta array should be in the range [0,PI]. Quitting..."
            write(*,*) Theta
            stop
        end if

        ! Only ghosts can go outside of PHI=[0,2*PI]
        if (addGhosts) then
            if (any(Phi<0.).or.(any(Phi>2*PI))) then
                write(*,*) "Inside shell grid generator (GenShellGrid)."
                write(*,*) "Phi array should be in the range [0,2*PI]. Quitting..."
                write(*,*) Phi
                stop
            end if
        else
            phiLo = shGr%Nge + 1
            phiHi = shGr%Nge + shGr%Np + 1

            if (phiLo < lbound(Phi,1) .or. phiHi > ubound(Phi,1)) then
                write(*,*) "Inside shell grid generator (GenShellGrid)."
                write(*,*) "Phi array ghost indexing is inconsistent with active-domain bounds. Quitting..."
                stop
            end if

            if (any(Phi(phiLo:phiHi)<0.).or.(any(Phi(phiLo:phiHi)>2*PI))) then
                write(*,*) "Inside shell grid generator (GenShellGrid)."
                write(*,*) "Phi array should be in the range [0,2*PI]. Quitting..."
                write(*,*) Phi
                stop
            end if
        end if


        ! decide if the grid is uniform in Phi
        ! helps to speed up interpolation search
        if (allocated(dphi)) deallocate(dphi)
        allocate(dphi(shGr%Np))   ! helper
        dphi = Phi(2:shGr%Np+1) - Phi(1:shGr%Np)
        if ( all( dphi - dphi(1) <= TINY ) ) then
            shGr%isPhiUniform = .true.
        else 
            shGr%isPhiUniform = .false.
        end if
        deallocate(dphi)

        associate(is=>shGr%is, ie=>shGr%ie, js=>shGr%js, je=>shGr%je, &
            isg=>shGr%isg, ieg=>shGr%ieg, jsg=>shGr%jsg, jeg=>shGr%jeg)

        ! Nuke arrays if already allocated
        if (allocated(shGr%th))    deallocate(shGr%th)
        if (allocated(shGr%ph))    deallocate(shGr%ph)
        if (allocated(shGr%lat))   deallocate(shGr%lat)
        if (allocated(shGr%thc))   deallocate(shGr%thc)
        if (allocated(shGr%phc))   deallocate(shGr%phc)
        if (allocated(shGr%latc))  deallocate(shGr%latc)
        if (allocated(shGr%thRp))  deallocate(shGr%thRp)
        if (allocated(shGr%thcRp)) deallocate(shGr%thcRp)
        if (allocated(shGr%areaCC)) deallocate(shGr%areaCC)
        if (allocated(shGr%lenFace_th)) deallocate(shGr%lenFace_th)
        if (allocated(shGr%lenFace_ph)) deallocate(shGr%lenFace_ph)
        ! Create new arrays
        allocate(shGr%th   (isg:ieg+1))
        allocate(shGr%thc  (isg:ieg  ))
        allocate(shGr%ph   (jsg:jeg+1))
        allocate(shGr%phc  (jsg:jeg  ))
        allocate(shGr%lat  (isg:ieg+1))
        allocate(shGr%latc (isg:ieg  ))
        allocate(shGr%thRp (isg:ieg+1))
        allocate(shGr%thcRp(isg:ieg  ))
        allocate(shGr%areaCC(isg:ieg,jsg:jeg))
        allocate(shGr%lenFace_th(isg:ieg+1,jsg:jeg  ))
        allocate(shGr%lenFace_ph(isg:ieg  ,jsg:jeg+1))



        if (.not. addGhosts) then
            ! Ghosts already included, we have enerything we need
            shGr%th = Theta
            shGr%ph = Phi
        else
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
                write(*,*)name
                write(*,*) Phi
                stop
            endif

            ! Note, we already made sure that for a periodic grid Nge=Ngw, phi(js)=0 and phi(je+1)=pi
            do j=1,shGr%Nge
                shGr%ph(js-j)   = shGr%ph(je+1-j) - 2*PI
                shGr%ph(je+1+j) = shGr%ph(js+j) + 2*PI
            end do
        endif

        ! Set centers
        shGr%thc = ( shGr%th(isg+1:ieg+1) + shGr%th(isg:ieg) )/2.
        shGr%phc = ( shGr%ph(jsg+1:jeg+1) + shGr%ph(jsg:jeg) )/2.

        ! Set latitude
        shGr%lat  = PI/2.0_rp - shGr%th
        shGr%latc = PI/2.0_rp - shGr%thc

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

        ! Project to planet surface
        xyz = 0.0
        do i=shGr%isg, shGr%ieg+1
            ! Don't need longitude information so just assume we are at Y=0
            xyz(1) = shGr%radius*sin(shGr%th(i))
            xyz(3) = shGr%radius*cos(shGr%th(i))
            L = DipoleL(xyz)
            shGr%thRp(i) = abs(asin(sqrt(1.0_rp/L)))
        enddo
        do i=shGr%isg, shGr%ieg
            xyz(1) = shGr%radius*sin(shGr%thc(i))
            xyz(3) = shGr%radius*cos(shGr%thc(i))
            L = DipoleL(xyz)
            shGr%thcRp(i) = abs(asin(sqrt(1.0_rp/L)))
        enddo

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


        ! Do some extra geom calculations for easy reference when needed
        do i=shGr%isg,shGr%ieg
            do j=shGr%jsg,shGr%jeg
                ! r^2 * (cos(th1)-cos(th2)) * dPh
                shGr%areaCC(i,j) = shGr%radius**2 &
                                 * ( cos(shGr%th(i)) - cos(shGr%th(i+1)) ) &
                                 * (shGr%ph(j+1) - shGr%ph(j))
            enddo
        enddo

        ! Faces in theta dir are i +/- 1/2 faces, each of which spans the phi direction
        do j=shGr%jsg,shGr%jeg
            do i=shGr%isg,shGr%ieg+1
                shGr%lenFace_th(i, j) = shGr%radius &
                                      * sin(shGr%th(i)) &
                                      * (shGr%ph(j+1) - shGr%ph(j))
            enddo
        enddo
        do j=shGr%jsg,shGr%jeg+1
            do i=shGr%isg,shGr%ieg
                shGr%lenFace_ph(i,j) = shGr%radius &
                                     * (shGr%th(i+1) - shGr%th(i))
            enddo
        enddo

        end associate

        shGr%name = name


    end subroutine GenShellGrid


    subroutine initShellVar(shGr, loc, shellVar, maskO)
        !! Inits a ShellGridVar thats associated with provided and initialized ShellGrid
        !! TODO: Maybe make mask alloc optional, in case we have an array of ShellGridVars with the same mask
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
        if (allocated(shellVar%mask)) deallocate(shellVar%mask)
            
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

        shellVar%data = 0.0_rp  ! initialize to 0

        ! Init mask, either by maskO or defaults
        if (present(maskO)) then
            if ( all(shape(shellVar%mask) == shape(maskO)) ) then
                shellVar%mask(:,:) = maskO(:,:)
            else
                write(*,*)"ERROR in initShellVar: maskO shape doesn't match."
                stop
            endif
        else
            shellVar%mask = .true.
        endif
    end subroutine initShellVar


    subroutine GenChildShellGrid(pSG, cSG, name, nGhosts, sub_is, sub_ie, sub_js, sub_je, ghostsOnParentGridO)
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
            !! Start and end cell-centered i/j indices of parent grid that bound active domain of new child grid
            !! These are optional. If left out, we will essentially make a copy of the parent grid
        logical, optional, intent(in) :: ghostsOnParentGridO
            !! Whether or not the ghost cells of the child grid should exactly match points on parent grid beyond defined subdomain

        integer :: is, ie, js, je
            !! Actual bounds used
        logical :: ghostsOnParentGrid

        if (present(ghostsOnParentGridO)) then
            ghostsOnParentGrid = ghostsOnParentGridO
        else
            ghostsOnParentGrid = .true.
        endif

        ! If a bound is provided then we use that, if not we default to parent grid's bounds
        is = whichInd(sub_is, pSG%is, present(sub_is))
        ie = whichInd(sub_ie, pSG%ie, present(sub_ie))
        js = whichInd(sub_js, pSG%js, present(sub_js))
        je = whichInd(sub_je, pSG%je, present(sub_je))

        ! Check for valid bounds
        if (is < pSG%is) then
            write(*,*) "ERROR GenChildShellGrid: Invalid is bound."
            write(*,*) "Requested:",is
            write(*,*) "Mimumum:",pSG%is
            stop
        endif
        if (ie > pSG%ie) then
            write(*,*) "ERROR GenChildShellGrid: Invalid ie bound."
            write(*,*) "Requested:",ie
            write(*,*) "Maximum:",pSG%ie
            stop
        endif
        if (js < pSG%js) then
            write(*,*) "ERROR GenChildShellGrid: Invalid js bound."
            write(*,*) "Requested:",js
            write(*,*) "Mimumum:",pSG%js
            stop
        endif
        if (je > pSG%je) then
            write(*,*) "ERROR GenChildShellGrid: Invalid je bound."
            write(*,*) "Requested:",je
            write(*,*) "Maximum:",pSG%je
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
            if (ghostsOnParentGrid) then
                ! The ghosts can't overrun parent grid's ghost bounds
                if (    is - nGhosts(NORTH) < pSG%isg  ) then
                    write(*,*) "ERROR GenChildShellGrid: Child ghosts overrun parent ghosts in is direction"
                    stop
                elseif (ie + nGhosts(SOUTH) > pSG%ieg+1) then
                    write(*,*) "ERROR GenChildShellGrid: Child ghosts overrun parent ghosts in ie direction"
                    stop
                elseif (js - nGhosts(EAST ) < pSG%jsg  ) then
                    write(*,*) "ERROR GenChildShellGrid: Child ghosts overrun parent ghosts in js direction"
                    stop
                elseif (je + nGhosts(WEST ) > pSG%jeg+1) then
                    write(*,*) "ERROR GenChildShellGrid: Child ghosts overrun parent ghosts in je direction"
                    stop
                endif
                ! Ghost cells have to be on parent grid. Include those cells in the theta/phi we provide and let GenShellGrid know
                call GenShellGrid(cSG, pSG%th(is-nGhosts(NORTH):ie+nGhosts(SOUTH)+1), pSG%ph(js-nGhosts(EAST):je+nGhosts(WEST)+1), name, nGhosts=nGhosts, radO=pSG%radius, TPIncludesGhostsO=.true.)    
            else
                call GenShellGrid(cSG, pSG%th(is:ie+1)                              , pSG%ph(js:je+1)                            , name, nGhosts=nGhosts, radO=pSG%radius, TPIncludesGhostsO=.false.)
            endif

        else
            ! No ghosts defined, go with default
            call GenShellGrid(cSG, pSG%th(is:ie+1), pSG%ph(js:je+1), name, radO=pSG%radius, TPIncludesGhostsO=.false.)
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

        if (allocated(sh%th)  ) deallocate(sh%th)
        if (allocated(sh%ph)  ) deallocate(sh%ph)
        if (allocated(sh%thc) ) deallocate(sh%thc)
        if (allocated(sh%phc) ) deallocate(sh%phc)
        if (allocated(sh%lat) ) deallocate(sh%lat)
        if (allocated(sh%latc)) deallocate(sh%latc)
        if (allocated(sh%latc)) deallocate(sh%thRp)
        if (allocated(sh%latc)) deallocate(sh%thcRp)

    end subroutine deallocShellGrid


    subroutine deallocShellGridVar(sgVar)
        type(ShellGridVar_T), intent(inout) :: sgVar

        sgVar%loc = -1
        sgVar%Ni  = 0
        sgVar%Nj  = 0
        sgVar%isv = 0
        sgVar%iev = 0
        sgVar%jsv = 0
        sgVar%jev = 0

        if (allocated(sgVar%data)) deallocate(sgVar%data)
        if (allocated(sgVar%mask)) deallocate(sgVar%mask)

    end subroutine deallocShellGridVar


end module shellGrid
