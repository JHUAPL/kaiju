! Various data structures and routines to do interpolation from (and to?) a spherical shell grid
! Routines are fairly broken up to try and allow for a reasonable balance of flexibility and performance
! TODO: Add routine to take list of scattered points and interpolate to ShellGrid_T
module shellInterp
    use kdefs
    use math
    use shellgrid
    use shellutils

    implicit none

    integer, parameter, private :: NumTSC = 9

    contains


    ! InterpShellVar - TODO
    ! InterpShellVar_ChildToParent - TODO
    ! InterpShellVar_ParentToChild - TODO

    subroutine InterpShellVar(sgSource, sgVar, sgDest, varOut)
        !! This is meant to be the highest-abstraction option for interpolation
        !! Not expected to be used very frequently, but maybe it makes sense for some simple variables
        !! Basically, the source and destination are well-defined within ShellGrid's intended use, 
        !!   so we can make the best decisions here on how to interpolate this variable
        type(ShellGrid_T)   , intent(in) :: sgSource
            !! Source shellGrid that sgVar lives on
        type(ShellGridVar_T), intent(in) :: sgVar
            !! the variable we are interpolating
        type(ShellGrid_T)   , intent(in) :: sgDest
            !! The destination grid
        type(ShellGridVar_T), intent(inout) :: varOut
            !! sgVar interpolated ontp sgDest


    end subroutine InterpShellVar


    subroutine InterpShellVar_TSC_SG(sgSource, sgVar, sgDest, varOut, dThetaO, dPhiO, srcMaskO)
        !! Interpolate a ShellGridVar to another ShellGrid using the Triangle-Schaped Cloud method
        type(ShellGrid_T)   , intent(in) :: sgSource
            !! Source shellGrid that sgVar lives on
        type(ShellGridVar_T), intent(in) :: sgVar
            !! the variable we are interpolating
        type(ShellGrid_T)   , intent(in) :: sgDest
            !! The destination grid
        type(ShellGridVar_T), intent(inout) :: varOut
            !! sgVar interpolated ontp sgDest
        real(rp), dimension(:), optional, intent(in) :: dThetaO
            !! Cell width in theta direction
        real(rp), dimension(:), optional, intent(in) :: dPhiO
            !! Width in theta direction
        logical, dimension(sgVar%isv:sgVar%iev, sgVar%jsv:sgVar%jev), optional, intent(in) :: srcMaskO
            !! If provided, we use this to mask the source variable instead of its mask


        integer :: extraPnt
        integer :: i,j
        real(rp), dimension(:), allocatable :: dTheta
        real(rp), dimension(:), allocatable :: dPhi
        logical :: goodInterp

        if ( present(dThetaO) ) then
            ! First, make sure they are the same shape
            if (size(dTheta) .ne. sgVar%Ni) then
                write(*,*) "WARNING: InterpShellVar_TSC_SG got a mismatched dThetaO shape. Dying."
                stop
            endif
            ! Otherwise, we can safely copy one to the other
            allocate(dTheta(sgVar%isv:sgVar%iev))
            dTheta = dThetaO
        else
            ! If dTheta not present, we calculate it ourselves

            if (sgVar%loc == SHGR_CC .or. sgVar%loc == SHGR_FACE_PHI) then
                ! Note that the location id is the variable's 1D location w.r.t. the theta axis
                call calcdx_TSC(sgSource%th, sgSource%isg, sgSource%ieg, SHGR_CC, dTheta)
            elseif (sgVar%loc == SHGR_CORNER .or. sgVar%loc == SHGR_FACE_THETA) then
                call calcdx_TSC(sgSource%th, sgSource%isg, sgSource%ieg, SHGR_CORNER, dTheta)
            endif
            !! Note: calcdx_TSC sets the dimension of dTheta and dPhi
        endif

        if ( present(dPhiO) ) then
            ! First, make sure they are the same shape
            if (size(dPhiO) .ne. sgVar%Nj) then
                write(*,*) "WARNING: InterpShellVar_TSC_SG got a mismatched dPhiO shape. Dying."
                stop
            endif
            ! Otherwise, we can safely copy one to the other
            allocate(dPhi(sgVar%jsv:sgVar%jev))
            dPhi = dPhiO
        else
            ! If dPhi not present, we calculate it ourselves
            if (sgVar%loc == SHGR_CC .or. sgVar%loc == SHGR_FACE_THETA) then
                call calcdx_TSC(sgSource%ph, sgSource%jsg, sgSource%jeg, SHGR_CC, dPhi)
            elseif (sgVar%loc == SHGR_CORNER .or. sgVar%loc == SHGR_FACE_PHI) then
                call calcdx_TSC(sgSource%ph, sgSource%jsg, sgSource%jeg, SHGR_CORNER, dPhi)
            endif
        endif

        ! Now that we have our dTheta and dPhi, we can start interpolating

        ! Which destination grid locations we loop over depends on the destination variable location
        select case(varOut%loc)
            case(SHGR_CC)
                !$OMP PARALLEL DO default(shared) &
                !$OMP schedule(dynamic) &
                !$OMP private(i,j,goodInterp)
                do j=sgDest%jsg,sgDest%jeg
                    do i=sgDest%isg,sgDest%ieg
                        if (.not. varOut%mask(i,j)) cycle
                        ! NOTE/TODO: This is where we would do transformations of destination grid's theta dn phi to source grid
                        ! in the case where they have different coordinate systems
                        call InterpShellVar_TSC_pnt( \
                                sgSource, sgVar,\
                                sgDest%thc(i), sgDest%phc(j),\
                                varOut%data(i,j), \
                                dThetaO=dTheta, dPhiO=dPhi,\
                                srcMaskO=srcMaskO, goodInterpO=goodInterp)
                        ! TODO: Handle case where goodInterp is false here
                        ! Probably will be model dependent. Maybe we return a 2D goodInterp array if an optional array is provided to us
                    enddo
                enddo
            case(SHGR_CORNER)
                !$OMP PARALLEL DO default(shared) &
                !$OMP schedule(dynamic) &
                !$OMP private(i,j,goodinterp)
                do j=sgDest%jsg,sgDest%jeg+1
                    do i=sgDest%isg,sgDest%ieg+1
                        if (.not. varOut%mask(i,j)) cycle
                        ! NOTE/TODO: This is where we would do transformations of destination grid's theta dn phi to source grid
                        ! in the case where they have different coordinate systems
                        call InterpShellVar_TSC_pnt( \
                                sgSource, sgVar,\
                                sgDest%th(i), sgDest%ph(j),\
                                varOut%data(i,j), \
                                dThetaO=dTheta, dPhiO=dPhi,\
                                srcMaskO=srcMaskO, goodInterpO=goodInterp)
                    enddo
                enddo
            case(SHGR_FACE_THETA)
                !$OMP PARALLEL DO default(shared) &
                !$OMP schedule(dynamic) &
                !$OMP private(i,j,goodInterp)
                do j=sgDest%jsg,sgDest%jeg
                    do i=sgDest%isg,sgDest%ieg+1
                        if (.not. varOut%mask(i,j)) cycle
                        ! NOTE/TODO: This is where we would do transformations of destination grid's theta dn phi to source grid
                        ! in the case where they have different coordinate systems
                        call InterpShellVar_TSC_pnt( \
                                sgSource, sgVar,\
                                sgDest%th(i), sgDest%phc(j),\
                                varOut%data(i,j), \
                                dThetaO=dTheta, dPhiO=dPhi,\
                                srcMaskO=srcMaskO, goodInterpO=goodInterp)
                    enddo
                enddo
            case(SHGR_FACE_PHI)
                !$OMP PARALLEL DO default(shared) &
                !$OMP schedule(dynamic) &
                !$OMP private(i,j,goodInterp)
                do j=sgDest%jsg,sgDest%jeg+1
                    do i=sgDest%isg,sgDest%ieg
                        if (.not. varOut%mask(i,j)) cycle
                        ! NOTE/TODO: This is where we would do transformations of destination grid's theta dn phi to source grid
                        ! in the case where they have different coordinate systems
                        call InterpShellVar_TSC_pnt( \
                                sgSource, sgVar,\
                                sgDest%thc(i), sgDest%ph(j),\
                                varOut%data(i,j), \
                                dThetaO=dTheta, dPhiO=dPhi,\
                                srcMaskO=srcMaskO, goodInterpO=goodInterp)
                    enddo
                enddo
        end select

        ! Fill j ghosts with active cell data
        call wrapJ_SGV(sgDest, varOut)
        

    end subroutine InterpShellVar_TSC_SG


    subroutine InterpShellVar_TSC_pnt(sgSource, sgVar, th, pin, Qinterp, dThetaO, dPhiO, srcMaskO, goodInterpO)
        !! Given the source information, interpolate sgVar to point (t,pin) and return as Qout
        type(ShellGrid_T   ), intent(in) :: sgSource
            !! Source ShellGrid
        type(ShellGridVar_T), intent(in) :: sgVar
            !! Variable relative to provided ShellGrid
        real(rp), intent(in) :: th
            !! Theta coordinate with respect to source grid
        real(rp), intent(in) :: pin
            !! Phi coordinate with respect to source grid
        real(rp), intent(out) :: Qinterp
            !! Interpolated value we are returning
        real(rp), dimension(sgVar%isv:sgVar%iev), optional, intent(in) :: dThetaO
            !! Cell width in theta, centered at variable positions on source grid
        real(rp), dimension(sgVar%jsv:sgVar%jev), optional, intent(in) :: dPhiO
            !! Cell width in phi, centered at variable positions on source grid
        logical, dimension(sgVar%isv:sgVar%iev, sgVar%jsv:sgVar%jev), optional, intent(in) :: srcMaskO
            !! If provided, we use this to mask the source variable instead of its mask
        logical, optional, intent(inout) :: goodInterpO
            !! True if we are returning a meaningful interpolated value

        
        real(rp) :: ph
            !! Cleaned-up phi location we actually use
        integer :: i0_cell, j0_cell
            !! i and j locations of the cell containing point t,p, irrespective of scVar%loc
        integer :: i0, j0
            !! i and j locations of point closest to t,p
            !! Whether they are respect to corner, center, or a face depends on sgVar%loc
        real(rp) :: dTh, dPh
            !! dTheta and dPhi at location i0,j0 , assuming we are in domain
        real(rp) :: t0, p0
            !! Theta and Phi values at point i0,j0
        real(rp) :: eta, zeta
        real(rp), dimension(-1:+1) :: wE,wZ
        real(rp), dimension(NumTSC) :: Ws,Qs
        logical , dimension(NumTSC) :: isGs
        integer :: ipnt,jpnt,n,di,dj

        ! Default return values
        Qinterp = 0.0
        if (present(goodInterpO)) goodInterpO = .false.

        ! Check bounds
        ph = modulo(pin,2*PI)

        ! Also make sure the requested theta is in bounds
        if ( (th<0.).or.(th>PI) ) then
            write(*,*) "ERROR in InterpShellVar_TSC_pnt: Theta should be in the range [0,PI]. Quitting..."
            stop
        end if

        call getSGCellILoc(sgSource, th, i0_cell, t0)
        if (sgVar%loc .eq. SHGR_CORNER .or. sgVar%loc .eq. SHGR_FACE_THETA) then
            call iLocCC2Corner(sgSource, th, i0_cell, iLocCornerO=i0, tLocO=t0)
        else
            i0 = i0_cell
        endif
        call getSGCellJLoc(sgSource, ph, j0_cell, p0)
        if (sgVar%loc .eq. SHGR_CORNER .or. sgVar%loc .eq. SHGR_FACE_PHI  ) then
            call jLocCC2Corner(sgSource, ph, j0_cell, jLocCornerO=j0, pLocO=p0)
        else
            j0 = j0_cell
        endif

        if (i0 > sgVar%iev .or. i0 < sgVar%isv) then
            return
        endif

        if (j0 > sgVar%jev .or. j0 < sgVar%jsv) then
            write(*,*) "ERROR in InterpShellVar_TSC_pnt: Phi out of bounds. idx=",j0
            write(*,*) "This wasn't supposed to be possible, good job."
            return
        endif

        ! Make sure data is good at this point
        if (.not. sgVar%mask(i0,j0)) then
            return
        endif

        ! If still here we're gonna do something, so we can tell our caller we are returning a valid value
        if (present(goodInterpO)) goodInterpO = .true.

        ! Determine dTheta and dPhi for this ij point
        if (present(dThetaO)) then
            ! First, make sure dTheta and dPhi are defined at sgVar locations     
            if ( size(dThetaO) .ne. sgVar%Ni ) then
                write(*,*)"ERROR in InterpShellVar_TSC_pnt: size(dTheta) != sgVar%Ni"
                write(*,*) size(dThetaO), sgVar%Ni
                stop
            endif
            dTh = dThetaO(i0)
        else
            ! dThetaO array not provided, so we calculate it ourselves
            if (sgVar%loc == SHGR_CC .or. sgVar%loc == SHGR_FACE_PHI) then
                dTh = Diff1D_4halfh(sgSource%th, sgSource%isg, sgSource%ieg  , i0)
            else if (sgVar%loc == SHGR_CORNER .or. sgVar%loc == SHGR_FACE_THETA) then
                dTh = Diff1D_4h    (sgSource%th, sgSource%isg, sgSource%ieg+1, i0)
            endif
        endif

        if (present(dPhiO)) then
            if ( size(dPhiO) .ne. sgVar%Nj ) then
                write(*,*)"ERROR in InterpShellVar_TSC_pnt: size(dPhi) != sgVar%Nj"
                write(*,*) size(dPhiO), sgVar%Nj
                stop
            endif
            dPh = dPhiO(j0)
        else
            ! dPhiO array not provided, so we calculate it ourselves
            if (sgVar%loc == SHGR_CC .or. sgVar%loc == SHGR_FACE_THETA) then
                dPh = Diff1D_4halfh(sgSource%ph, sgSource%jsg, sgSource%jeg  , j0)
            else if (sgVar%loc == SHGR_CORNER .or. sgVar%loc == SHGR_FACE_PHI) then
                dPh = Diff1D_4h    (sgSource%ph, sgSource%jsg, sgSource%jeg+1, j0)
            endif
        endif

        ! First, check if active grid has poles
        ! note, if the destination point is closer to the top cell boundary
        ! than the bottom (for corner centered variables), then i0 by now is 2.
        ! in other words, only the half of the cell closest to the pole will be interpolated as a pole.
        if (sgSource%doNP .and. (i0_cell==sgSource%is)) then
            ! Handle north pole and return
            call interpPole(sgSource,sgVar,th,ph,Qinterp)
            return
        endif

        ! same comment as above but for south pole
        if (sgSource%doSP .and. (i0_cell==sgSource%ie)) then
            ! Handle south pole and return
            call interpPole(sgSource,sgVar,th,ph,Qinterp)
            return
        endif

        ! Now, if ghost grid has poles
        if (sgSource%ghostSP .and. (i0==sgsource%ieg)) then
            write(*,*) "Not implemented!"
            stop
        endif

        if (sgSource%ghostNP .and. (i0==sgSource%isg)) then
            write(*,*) "Not implemented!"
            stop
        endif

        ! Note: If still here we know i0 isn't on the boundary

        eta  = (th - t0)/dTh
        zeta = (ph - p0)/dPh

        call ClampMapVar(eta)
        call ClampMapVar(zeta)

        ! Calculate weights
        call TSCweight1D(eta ,wE)
        call TSCweight1D(zeta,wZ)

        ! Collect weights/values
        n = 1
        do dj=-1,+1
            do di=-1,+1
                ipnt = i0+di
                jpnt = j0+dj
                
                ! Wrap around boundary
                if (jpnt<1)           jpnt = sgSource%Np
                if (jpnt>sgSource%Np) jpnt = 1

                ! Do zero-grad for theta
                if (ipnt<sgVar%isv)   ipnt = sgVar%isv
                if (ipnt>sgVar%iev)   ipnt = sgVar%iev

                Qs(n) = sgVar%data(ipnt,jpnt)
                Ws(n) = wE(di)*wZ(dj)
                if(present(srcMaskO)) then
                    isGs(n) = srcMaskO(ipnt,jpnt)
                else
                    isGs(n) = sgVar%mask(ipnt,jpnt)
                endif

                if (.not. isGs(n)) Ws(n) = 0.0

                n = n + 1
            enddo
        enddo

        ! Renormalize
        Ws = Ws/sum(Ws)

        ! Get final value
        Qinterp = dot_product(Qs,Ws)

        ! Have some internal functions
        contains
        
        ! Clamps mapping in [-0.5,0.5]
        subroutine ClampMapVar(ez)
            REAL(rp), intent(inout) :: ez
            if (ez<-0.5) ez = -0.5
            if (ez>+0.5) ez = +0.5
        end subroutine ClampMapVar

        ! 1D triangular shaped cloud weights
        ! 1D weights for triangular shaped cloud interpolation
        ! Assuming on -1,1 reference element, dx=1
        ! Check for degenerate cases ( |eta| > 0.5 )
        subroutine TSCweight1D(eta,wE)
            real(rp), intent(in)  :: eta
            real(rp), intent(out) :: wE(-1:1)

            wE(-1) = 0.5*(0.5-eta)**2.0
            wE( 1) = 0.5*(0.5+eta)**2.0
            wE( 0) = 0.75 - eta**2.0

        end subroutine TSCweight1D

    end subroutine InterpShellVar_TSC_pnt


    subroutine InterpShellVar_ParentToChild(srcGrid, srcVar, destGrid, destVar)
        type(ShellGrid_T   ), intent(in) :: srcGrid
        type(ShellGridVar_T), intent(in) :: srcVar
        type(ShellGrid_T   ), intent(in) :: destGrid
        type(ShellGridVar_T), intent(inout) :: destVar

        integer :: subis, subie, subjs, subje, iExtra, jExtra
        type(ShellGridVar_T) :: src2destVar
            !! srcVar chopped down to dest grid's domain 

        ! Are you actually my parent or do I need to call the police?
        if (trim(destGrid%parentName) == def_noParentName .or. trim(srcGrid%name) .ne. trim(destGrid%parentName)) then
            write(*,*) " ERROR in InterpShellVar_ParentToChild: Src grid not parent of Dest grid"
            write(*,*) " Src: ", trim(srcGrid%name), ", Dest: ",trim(destGrid%parentName)
            stop
        endif

        ! Determine SGV bounds of source var that map to destination grid's bounds
        call ijExtra_SGV(srcVar%loc, iExtra, jExtra)

        subis = destGrid%is + destGrid%bndis-1 - destGrid%Ngn
        subie = destGrid%is + destGrid%bndie-1 + destGrid%Ngs + iExtra
        subjs = destGrid%js + destGrid%bndjs-1 - destGrid%Ngw
        subje = destGrid%js + destGrid%bndje-1 + destGrid%Nge + jExtra

        ! Now fill in destination variable
        if (srcVar%loc == destVar%loc) then
            destVar%data(:,:) = srcVar%data(subis:subie, subjs:subje)
            destVar%mask(:,:) = srcVar%mask(subis:subie, subjs:subje)
        else  ! Need to interpolate
            ! Chop srcVar into destination grid size
            call initShellVar(destGrid, srcVar%loc, src2destVar)
            src2destVar%data(:,:) = srcVar%data(subis:subie, subjs:subje)
            src2destVar%mask(:,:) = srcVar%mask(subis:subie, subjs:subje)
            call InterpShellVar_TSC_SG(destGrid, src2destVar, destGrid, destVar)
        endif

    end subroutine InterpShellVar_ParentToChild


    subroutine interpPole(shGr,Qin,tin,pin,Qinterp)
        type(ShellGrid_T), intent(in)     :: shGr     ! source grid
        type(ShellGridVar_T), intent(in)  :: Qin      ! source grid variable
        real(rp), intent(out)             :: Qinterp  ! interpolated variable
        real(rp), intent(in)              :: tin,pin  ! theta,phi of the destination point
        real(rp) :: f0,f1,f2,I1,I2,dphi,phi
        integer  :: ishift,j,pole
        integer  :: iinterp ! the index of the pole cell (first/last for NORTH/SOUTH)
        real(rp) :: tfactor,Qtemp

        ! Find out which pole we're at
        if ( (tin.ge.shGr%th(shGr%is)).and.(tin.le.shGr%th(shGr%is+1)) ) then
            ! note, shGr%th(shGr%is) is within TINY of 0 since we are at the pole
            pole = NORTH
            iinterp = shGr%is
        else if ( (tin.le.shGr%th(shGr%ie+1)).and.(tin.ge.shGr%th(shGr%ie)) ) then
            ! note, shGr%th(shGr%ie+1) is within TINY of PI since we are at the pole
            pole = SOUTH
            iinterp = shGr%ie
        else
            write(*,*) "Attempting to call pole interpolation for a point that's not on the pole. Quitting..."
            write(*,*)"tin=",tin
            write(*,*)"lower cell bounds: (",shGr%th(shGr%is),",",shGr%th(shGr%is+1),")"
            write(*,*)"upper cell bounds: (",shGr%th(shGr%ie),",",shGr%th(shGr%ie+1),")"
            
        end if

        ! represent the function near pole to first order in theta as
        ! f(t,p) = f0 + f1*cos(p)*t + f2*sin(p)*t 
        ! (Lewis&Bellan, J. Math. Phys. 31, 2592 (1990); 
        ! https://doi.org/10.1063/1.529009
        ! 
        ! for corner centered and theta face centered, 
        ! iinterp will be right on the pole for NORTH, so we will have to shift up
        ! otherwise, there's no shift
        ishift = 0

        ! check centering wrt theta
        select case(Qin%loc)
        case(SHGR_CC, SHGR_FACE_PHI)
            if (pole.eq.NORTH) then
                tfactor = tin/shGr%thc(iinterp)
            else
                tfactor = (PI-tin)/(PI-shGr%thc(iinterp))
            end if  
        
        case(SHGR_CORNER, SHGR_FACE_THETA)
            if (pole.eq.NORTH) then
                ishift = 1
                ! note, using th, not thc, since we're on a theta face
                tfactor = tin/shGr%th(iinterp+ishift)
            else
                tfactor = (PI-tin)/(PI-shGr%th(iinterp))
            end if  

        case default
            write(*,*) "interpPole got an invalid data location:",Qin%loc
            stop
        end select 

        ! determine f0, f1, f2 from the Fourier transform
        Qinterp = 0.0
        f0      = 0.0
        f1      = 0.0
        f2      = 0.0

        do j=1,shGr%Np
            ! check centering wrt phi
            select case(Qin%loc)
            case(SHGR_CC, SHGR_FACE_THETA)
                Qtemp = Qin%data(iinterp+ishift,j)
            case(SHGR_CORNER, SHGR_FACE_PHI)  
                Qtemp = 0.5*(Qin%data(iinterp+ishift,j)+Qin%data(iinterp+ishift,j+1))
            end select 

            dphi = shGr%ph(j+1)-shGr%ph(j)
            phi  = shGr%phc(j)
            f0   = f0 + Qtemp*dphi
            f1   = f1 + Qtemp*dphi*cos(phi)
            f2   = f2 + Qtemp*dphi*sin(phi)                
        end do

        f0 = f0/(2.*PI)
        f1 = f1/PI
        f2 = f2/PI
        
        Qinterp = f0 + ( f1*cos(pin) + f2*sin(pin) )*tfactor

    end subroutine interpPole

    ! this version of the function attempts to do the expansion to an arbitrary order
    ! however, it's incorrect in that it neglects the higher order terms in theta
    ! in the polynomials multiplying the harmonics in eqn. 9 of Lewis & Bellan, 1990
    ! in other words, it only retains f_m^(0) terms in that eqn. which is technically incorrect
    ! although the error is trivially small and in fact the result looks slightly better
    subroutine interpPoleHighOrder(shGr,Qin,tin,pin,Qinterp)
        ! note, calling this function with order=1 is equivalent to calling InterpPole above
        
        type(ShellGrid_T), intent(in)     :: shGr     ! source grid
        type(ShellGridVar_T), intent(in)  :: Qin      ! source grid variable
        real(rp), intent(out)             :: Qinterp  ! interpolated variable
        real(rp), intent(in)              :: tin,pin  ! theta,phi of the destination point
        real(rp) :: f0,I1,I2,dphi,phi
        real(rp),dimension(:,:),allocatable :: fcoef
        integer  :: ishift,j,pole
        integer  :: iinterp ! the index of the pole cell (first/last for NORTH/SOUTH)
        real(rp) :: tfactor,Qtemp
        integer :: oind, order=1 ! 12 -- I ran it up to 12th order just for funsies

        write(*,*)"WARNING: unless you are Slava or have talked to him about pole interpolation you shouldn't be seeing this message"

        ! Find out which pole we're at
        if ( (tin.ge.shGr%th(shGr%is)).and.(tin.le.shGr%th(shGr%is+1)) ) then
            ! note, shGr%th(shGr%is) is within TINY of 0 since we are at the pole
            pole = NORTH
            iinterp = shGr%is
        else if ( (tin.le.shGr%th(shGr%ie+1)).and.(tin.ge.shGr%th(shGr%ie)) ) then
            ! note, shGr%th(shGr%ie+1) is within TINY of PI since we are at the pole
            pole = SOUTH
            iinterp = shGr%ie
        else
            write(*,*) "Attempting to call pole interpolation for a point that's not on the pole. Quitting..."
        end if

        ! represent the function near pole to first order in theta as
        ! f(t,p) = f0 + f1*cos(p)*t + f2*sin(p)*t + higher order terms
        ! (Lewis&Bellan, J. Math. Phys. 31, 2592 (1990); 
        ! https://doi.org/10.1063/1.529009
        ! 
        ! for corner centered and theta face centered, 
        ! iinterp will be right on the pole for NORTH, so we will have to shift up
        ! otherwise, there's no shift
        ishift = 0

        ! check centering wrt theta
        select case(Qin%loc)
        case(SHGR_CC, SHGR_FACE_PHI)
            if (pole.eq.NORTH) then
                tfactor = tin/shGr%thc(iinterp)
            else
                tfactor = (PI-tin)/(PI-shGr%thc(iinterp))
            end if  
        
        case(SHGR_CORNER, SHGR_FACE_THETA)
            if (pole.eq.NORTH) then
                ishift = 1
                ! note, using th, not thc, since we're on a theta face
                tfactor = tin/shGr%th(iinterp+ishift)
            else
                tfactor = (PI-tin)/(PI-shGr%th(iinterp))
            end if  

        case default
            write(*,*) "interpPole got an invalid data location:",Qin%loc
            stop
        end select 

        ! determine f0, f1, f2 from the Fourier transform
        Qinterp = 0.0
        f0      = 0.0

        if (allocated(fcoef)) deallocate(fcoef)
        allocate(fcoef(order,2))
        fcoef(:,:)  = 0.0

        do j=1,shGr%Np
            ! check centering wrt phi
            select case(Qin%loc)
            case(SHGR_CC, SHGR_FACE_THETA)
                Qtemp = Qin%data(iinterp+ishift,j)
            case(SHGR_CORNER, SHGR_FACE_PHI)  
                Qtemp = 0.5*(Qin%data(iinterp+ishift,j)+Qin%data(iinterp+ishift,j+1))
            end select 

            dphi = shGr%ph(j+1)-shGr%ph(j)
            phi  = shGr%phc(j)
            f0   = f0 + Qtemp*dphi

            do oind=1,order
                fcoef(oind,1) = fcoef(oind,1) + Qtemp*dphi*cos(oind*phi)
                fcoef(oind,2) = fcoef(oind,2) + Qtemp*dphi*sin(oind*phi)                
            enddo
        end do 

        f0 = f0/(2.*PI)

        do oind=1,order
            fcoef(oind,:) = fcoef(oind,:)/PI
            Qinterp = Qinterp + ( fcoef(oind,1)*cos(oind*pin) + fcoef(oind,2)*sin(oind*pin) )*tfactor**oind
        enddo

        Qinterp = Qinterp + f0
    end subroutine interpPoleHighOrder
    

!------
! Low-level interp helpers
!------

    subroutine calcdx_TSC(x, isg, ieg, loc, dx)
        !! Calculates dx for positional array x using a 4-point stencil
        !! This is desirable for TSC so that we can have better accuracy 
        !! in the case of a non-uniformly space grid
        
        ! Note, even though we are not doing much in this function, we are 
        ! breaking it out here in case we want to do higher-order for TSC or something later
        real(rp), dimension(isg:ieg+1), intent(in) :: x
            !! 1D spatial grid, should always be cell corners
        integer, intent(in) :: isg, ieg
            !! Start/end indices of cell centers
        integer, intent(in) :: loc
            !! Location identifier, MUST BE SHGR_CC OR SHGR_CORNER
            !! This is the location of the points we are calculating dx at relative to x
        real(rp), dimension(:), allocatable, intent(out) :: dx
            !! 'cell width' we return

        integer :: i

        if (allocated(dx)) deallocate(dx)

        if (loc == SHGR_CORNER) then
            
            allocate(dx(isg:ieg+1))
            do i=isg,ieg+1
                dx(i) = Diff1D_4h(x, isg, ieg+1, i)
            enddo

        else if (loc == SHGR_CC) then

            allocate(dx(isg:ieg))
            do i=isg,ieg
                dx(i) = Diff1D_4halfh(x, isg, ieg, i)
            enddo

        else
            write(*,*) "ERROR: Invalid location id in calcdx_TSC. Must be SHGR_CC or SHGR_CORNER"
            stop
        endif

    end subroutine calcdx_TSC


    function Diff1D_4halfh(Q,is,ie,i0) result(Qp)
        !! Use 4-point stencil to calculate first derivative of coordinates
        !! This is for the case where we are using corners to calculate the difference at cell center
        !! adapted from chimp/ebinit.F90
        real(rp), intent(in) :: Q(is:ie+1)
            !! Cell corners
        integer, intent(in) :: is,ie
            !! Start and end indices of bounding grid
        integer, intent(in) :: i0
            !! Index of the point we are evaluating, offset half a cell from its lower bound Q(i0)
        real(rp) :: Qp
        real(rp) :: Qblk(4),c(4)

        ! Note that we have fewer cases than Diff1D_4h
        ! That's because even when i0 is the last cell-centered coordinate,
        ! we still have 1 usable corner value before we reach the void
        if (i0 == is) then
            ! Q coordinates at -0.5,0.5,1.5,2.5 relative to our point
            Qblk = [Q(is), Q(is+1), Q(is+2), Q(is+3)]
            c = [-23.0, 21.0, 3.0, -1.0]/25.0
        else if (i0 == ie) then
            ! Q coordinates at -2.5,-1.5,-0.5,0.5 relative to our point
            Qblk = [Q(ie-2), Q(ie-1), Q(ie), Q(ie+1)]
            c = [1.0, -3.0, -21.0, 23.0]/25.0
        else
            ! Q coordinates at -1.5,-0.5,0.5,1.5 relative to our point
            Qblk = [Q(i0-1), Q(i0), Q(i0+1), Q(i0+2)]
            c = [1.0, -27.0, 27.0, -1.0]/24.0
        endif
        Qp = dot_product(Qblk,c)

    end function Diff1D_4halfh
    

    function Diff1D_4h(Q,is,ie,i0) result(Qp)
        !! Use 4-point stencil to calculate first derivative of coordinates
        !! This is for the case where position we are calculating the difference for is at Q(i0)
        !! (In contrast to e.g. using cell corner values to calculate the difference located at the cell center)
        !! adapted from chimp/ebinit.F90
        real(rp), intent(in) :: Q(is:ie)
            !! Cell corners
        integer, intent(in) :: is,ie,i0
        real(rp) :: Qp
        real(rp) :: Qblk(4),c(4)
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
        
    end function Diff1D_4h

end module shellInterp
