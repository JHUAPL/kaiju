module raijuDomain

    use raijudefs
    use raijuTypes
    use raijugrids
    use math

    implicit none

    contains


    subroutine setActiveDomain(Model, Grid, State)
        !! Sets State%active array
        !! Determines Inactive -> buffer -> active cells
        type(raijuModel_T), intent(in   ) :: Model
        type(raijuGrid_T ), intent(in   ) :: Grid
        type(raijuState_T), intent(inout) :: State

        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: isInactive
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: isBuffer


        !call getInactiveCells(Model, Grid, State, isInactive)
        call getInactiveCells_domWeights(Model, Grid, State, isInactive, State%domWeights)
        call getBufferCells  (Model, Grid, State, isInactive, isBuffer, State%OCBDist)

        ! Set zones
        where (isInactive)
            State%active = RAIJUINACTIVE
        else where (isBuffer)
            State%active = RAIJUBUFFER
        elsewhere
            State%active = RAIJUACTIVE
        end where

        ! HOWEVER...
        ! Force all cells below X Re to be active, with enough buffer cells, and suffer the consequences if that includes bad flux tubes
        call forceActiveWithinRadius(Model, Grid, State, State%active)

    end subroutine setActiveDomain


    subroutine getInactiveCells_domWeights(Model, Grid, State, isInactive, domWeights)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: isInactive
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: domWeights

        integer :: i,j
        !real(rp) :: xyMin
        real(rp), dimension(2,2) :: rMin22
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: vaFrac_cc, bmin_cc
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: tmpWeights
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: tmpPass
        integer, dimension(Grid%shGrid%jsg:Grid%shGrid%jeg) :: bndLoc, bndR, bndL
        integer :: n_bndLim, n_contig
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: isInactive_tmp

        isInactive = .false.
        domWeights = 1.0_rp
        bndLoc     = Grid%shGrid%isg
        vaFrac_cc  = 0.0_rp
        bmin_cc    = 0.0_rp

        associate(sh=>Grid%shGrid)

        ! Start with hard requirement of closed field line
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,rMin22)
        do j=sh%jsg, sh%jeg
            do i=sh%isg, sh%ieg
                if (any(State%topo(i:i+1,j:j+1) .eq. RAIJUOPEN)) then
                    isInactive(i,j) = .true.
                    domWeights(i,j) = 0.0_rp
                    cycle
                endif
                
                bmin_cc(i,j)   = 0.25*(sum(State%Bmin  (i:i+1,j:j+1,ZDIR)))
                vaFrac_cc(i,j) = 0.25*(sum(State%vaFrac(i:i+1,j:j+1)))
                
                ! Also, apply buffer outer boundary constraint
                if ( any(State%radMax(i:i+1,j:j+1) > Model%maxTail_buffer) .or. any(State%radMax(i:i+1,j:j+1) > Model%maxSun_buffer) ) then
                    isInactive(i,j) = .true.
                    domWeights(i,j) = 0.0_rp
                endif
            enddo
        enddo

        ! Apply soft/hard limits
        call evalDomLimits(bmin_cc, Model%lim_bmin_soft    , Model%lim_bmin_hard   , &
                           .not. isInactive, weights2D=tmpWeights, pass2D=tmpPass)
        isInactive = isInactive .or. (.not. tmpPass)
        domWeights = domWeights*tmpWeights
        
        call evalDomLimits(vaFrac_cc, Model%lim_vaFrac_soft, Model%lim_vaFrac_hard, &
                           .not. isInactive, weights2D=tmpWeights, pass2D=tmpPass)
        isInactive = isInactive .or. (.not. tmpPass)
        domWeights = domWeights*tmpWeights


        ! Now that we applied all criteria, update bndLoc
        do j=sh%jsg, sh%jeg
            do i=sh%isg, sh%ieg
                if (isInactive(i,j) .and. (i > bndLoc(j))) then
                    bndLoc(j) = i
                endif
            enddo
        enddo

        ! Handle periodic boundary in j direction
        !isInactive(:, sh%jsg:sh%js-1) = isInactive(:, sh%je-sh%Ngw+1:sh%je)
        !isInactive(:, sh%je+1:sh%jeg) = isInactive(:, sh%js:sh%js+sh%Nge-1)
        call wrapJcc(sh, isInactive)
        !isInactive(:, sh%jsg:sh%js-1) = isInactive(:, sh%je-sh%Ngw+1:sh%je)
        !isInactive(:, sh%je+1:sh%jeg) = isInactive(:, sh%js:sh%js+sh%Nge-1)
        bndLoc(sh%jsg:sh%js-1) = bndLoc(sh%je-sh%Nge+1:sh%je)
        bndLoc(sh%je+1:sh%jeg) = bndLoc(sh%js:sh%js+sh%Ngw-1)

        ! Try to make boundary smooth in ionosphere
        bndL = bndLoc
        bndR = bndLoc
        n_bndLim = Model%n_bndLim

        ! Right sweep
        do j=sh%js,sh%je
            bndR(j) = max(bndR(j), bndR(j-1)-n_bndLim)
        enddo
        bndR(sh%je+1:sh%jeg) = bndR(sh%je)  ! Extend final right sweep result into j ghosts
        ! Left sweep
        do j=sh%je,sh%js, -1
            bndL(j) = max(bndL(j), bndL(j+1)-n_bndLim)
        enddo
        bndL(sh%jsg:sh%js-1) = bndL(sh%js)  ! Extend final left sweep into ghosts
        ! Now combine back into bndLoc
        do j=sh%jsg,sh%jeg
            bndLoc(j) = max(bndL(j), bndR(j))
            isInactive(sh%isg:bndLoc(j), j) = .true.
        enddo

        ! Finally, kick out any stumpy regions
        n_contig=4
        isInactive_tmp = isInactive
        do i=sh%isg, sh%ieg-n_contig
            do j=sh%js, sh%je-n_contig
                if (isInactive(i,j) .and. isInactive(i+n_contig,j)) then
                    isInactive_tmp(i:i+n_contig, j) = .true.
                endif
                if (isInactive(i,j) .and. isInactive(i,j+n_contig)) then
                    isInactive_tmp(i,j:j+n_contig) = .true.
                endif
            enddo
        enddo
        isInactive = isInactive_tmp
        call wrapJcc(sh, isInactive)

        end associate

        contains


        subroutine evalDomLimits(var, softLim, hardLim, isGood, weights2D, pass2D)
            real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: var
            real(rp), intent(in) :: softLim
            real(rp), intent(in) :: hardLim
            logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isGood
            real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(out) :: weights2D
            logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(out) :: pass2D

            real(rp) :: varSm_pnt
                !! Smoothed var at point i,j
            real(rp) wgt_pnt
                !! Weight after smoothing for point i,j

            do j=Grid%shGrid%jsg, Grid%shGrid%jeg
                do i=Grid%shGrid%isg, Grid%shGrid%ieg
                    if (.not. isGood(i,j)) then
                        weights2D(i,j) = 0.0_rp
                    else
                        weights2D(i,j) = (var(i,j) - hardLim)/(softLim - hardLim)
                        call ClampValue(weights2D(i,j), 0.0_rp, 1.0_rp)
                            !! 0 = hard fail, 1 = soft pass, in between = hard pass soft fail
                            !! Gives the behavior we want for both var < soft < hard and var > soft > hard
                    endif
                enddo
            enddo

            ! Now that we have everyone's weights, see if smoothing can save some points
            pass2D = .false.
            do j=Grid%shGrid%jsg, Grid%shGrid%jeg
                do i=Grid%shGrid%isg, Grid%shGrid%ieg
                    if (abs(weights2D(i,j) - 1.0_rp) < TINY) then
                        pass2D(i,j) = .true.
                    else if (weights2D(i,j) > TINY .and. notGhost(i,j)) then
                        ! Give it a second chance based on how good surrounding points are
                        ! (Ghosts don't get second chance, they are full pass or full fail)
                        varSm_pnt = SmoothOperator55(var(i-2:i+2,j-2:j+2), weights2D(i-2:i+2,j-2:j+2) > TINY)
                        wgt_pnt = (varSm_pnt - hardLim)/(softLim - hardLim)
                        if (abs(wgt_pnt - 1.0_rp) < TINY) then
                            pass2D(i,j) = .true.
                        endif
                    endif
                enddo
            enddo
        end subroutine evalDomLimits

        function notGhost(i,j)
            integer, intent(in) :: i,j
            logical :: notGhost
            if ( (i >= Grid%shGrid%is) .and. (i <= Grid%shGrid%ie) .and. (j >= Grid%shGrid%js) .and. (j <= Grid%shGrid%je)) then
                notGhost = .true.
            else
                notGhost = .false.
            endif
        end function notGhost


    end subroutine getInactiveCells_domWeights


!    subroutine getInactiveCells(Model, Grid, State, isInactive, isCoreInactiveO)
!        !! Applies series of criteria to determine which cells should be marked as inactive
!        type(raijuModel_T), intent(in) :: Model
!        type(raijuGrid_T ), intent(in) :: Grid
!        type(raijuState_T), intent(in) :: State
!        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: isInactive
!            !! The ultimate domain we want to call inactive, don't trust at all
!        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), &
!                optional, intent(inout) :: isCoreInactiveO
!            !! Debug option: Contians the inactive region according to only core physical constraints, no extra trimming
!
!        integer :: i,j
!        logical :: iShellHasCheck
!        real(rp) :: bndRateLim
!        integer :: n_bndLim, n_contig
!        integer, dimension(Grid%shGrid%jsg:Grid%shGrid%jeg) :: bndLoc, bndR, bndL
!        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: cornerNormAngle
!        integer, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: ocbDist
!        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: isCoreInactive, checkMask, isInactive_tmp
!
!        isInactive = .false.
!        checkMask = .false.
!        bndLoc = Grid%shGrid%isg
!        call calcCoreConstraints(Model, Grid, State, isInactive, bndLoc)
!
!        ! Prep trim criteria
!        call calcCornerNormAngles(Model, Grid, State, cornerNormAngle)
!
!        ! Now do trimming
!        associate(sh=>Grid%shGrid)
!        
!        ! Get dist from OCb
!        call CalcOCBDist(sh, isCoreInactive, 4, ocbDist)
!        ! Highlight areas to start checking for extra constraints
!        where(ocbDist > 0 .and. ocbDist < 5)
!            checkMask = .true.
!        end where
!        
!        !write(*,*)"NOTE: turning off normAngle check"
!        do i=sh%isg, sh%ie ! NOTE: Not touching upper theta ghosts, we shouldn't reach there anyways. If we do then we have bigger problems
!            iShellHasCheck = .false.
!            if ( all(State%active(i,:) == RAIJUINACTIVE) ) then
!                cycle
!            endif
!            do j=sh%js, sh%je  ! NOTE: Not touching j ghosts
!                if (checkMask(i,j)) then
!                    iShellHasCheck = .true.
!                    ! Criteria check
!                    if ( any(State%Pstd(i:i+1,j:j+1,0) > Model%PstdThresh)) then
!                    !if ( any(State%Pstd(i:i+1,j:j+1,0) > Model%PstdThresh) .or. any(cornerNormAngle(i:i+1,j:j+1) < Model%normAngThresh)) then
!                        isInactive(i,j) = .true.
!                        if (i > bndLoc(j)) then
!                            bndLoc(j) = i
!                        endif
!                        ! Flag points a little depper into the domain as places to check for criteria pass
!                        checkMask(i+1,j-2:j+2) = .true.
!                        checkMask(i+2,j-1:j+1) = .true.
!                        checkMask(i+3,j      ) = .true.
!                    endif
!                endif
!            enddo
!            if (.not. iShellHasCheck) then
!                exit
!            endif
!        enddo
!        ! Wrap in J
!        !call wrapJcc(sh, isInactive)
!        isInactive(:, sh%jsg:sh%js-1) = isInactive(:, sh%je-sh%Ngw+1:sh%je)
!        isInactive(:, sh%je+1:sh%jeg) = isInactive(:, sh%js:sh%js+sh%Nge-1)
!        bndLoc(sh%jsg:sh%js-1) = bndLoc(sh%je-sh%Ngw+1:sh%je)
!        bndLoc(sh%je+1:sh%jeg) = bndLoc(sh%js:sh%js+sh%Nge-1)
!
!        bndL = bndLoc
!        bndR = bndLoc
!        n_bndLim = Model%n_bndLim
!
!        ! Right sweep
!        do j=sh%js,sh%je
!            bndR(j) = max(bndR(j), bndR(j-1)-n_bndLim)
!        enddo
!        bndR(sh%je+1:sh%jeg) = bndR(sh%je)  ! Extend final right sweep result into j ghosts
!        ! Left sweep
!        do j=sh%je,sh%js, -1
!            bndL(j) = max(bndL(j), bndL(j+1)-n_bndLim)
!        enddo
!        bndL(sh%isg:sh%is-1) = bndL(sh%is)  ! Extend final left sweep into ghosts
!        ! Now combine back into bndLoc
!        do j=sh%jsg,sh%jeg
!            bndLoc(j) = max(bndL(j), bndR(j))
!            isInactive(sh%isg:bndLoc(j), j) = .true.
!        enddo
!
!        ! Finally, kick out any stumpy regions
!        n_contig=4
!        isInactive_tmp = isInactive
!        do i=sh%isg, sh%ieg-n_contig
!            do j=sh%js, sh%je-n_contig
!                if (isInactive(i,j) .and. isInactive(i+n_contig,j)) then
!                    isInactive_tmp(i:i+n_contig, j) = .true.
!                endif
!                if (isInactive(i,j) .and. isInactive(i,j+n_contig)) then
!                    isInactive_tmp(i,j:j+n_contig) = .true.
!                endif
!            enddo
!        enddo
!        isInactive = isInactive_tmp
!
!        end associate
!
!        contains
!
!        subroutine calcCoreConstraints(Model, Grid, State, isCoreInactive, bndLoc)
!            type(raijuModel_T), intent(in) :: Model
!            type(raijuGrid_T ), intent(in) :: Grid
!            type(raijuState_T), intent(in) :: State
!            logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: isCoreInactive
!            integer, dimension(Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: bndLoc
!
!            integer :: i,j
!            real(rp) :: xyMin
!            real(rp), dimension(2,2) :: bminSquare
!
!            isCoreInactive = .false.
!            do j=Grid%shGrid%jsg, Grid%shGrid%jeg
!                do i=Grid%shGrid%isg, Grid%shGrid%ieg
!                    
!
!                    ! Any cell with an open corner is bad
!                    if (any(State%topo(i:i+1,j:j+1) .eq. RAIJUOPEN)) then
!                        isCoreInactive(i,j) = .true.
!                    endif
!
!                    if (any(State%vaFrac(i:i+1,j:j+1) .le. Model%vaFracThresh)) then
!                        isCoreInactive(i,j) = .true.
!                    endif
!
!                    bminSquare(1,1) = norm2(State%Bmin(i  ,j  ,:))
!                    bminSquare(2,1) = norm2(State%Bmin(i+1,j  ,:))
!                    bminSquare(1,2) = norm2(State%Bmin(i  ,j+1,:))
!                    bminSquare(2,2) = norm2(State%Bmin(i+1,j+1,:))
!                    if (any( bminSquare .le. Model%bminThresh) ) then
!                    !if (any( norm2(State%Bmin(i:i+1,j:j+1,:),dim=3) .le. Model%bminThresh) ) then
!                        isCoreInactive(i,j) = .true.
!                    endif
!
!                    xyMin = sqrt(State%xyzMin(i,j,XDIR)**2 + State%xyzMin(i,j,YDIR)**2)
!                    ! Simple circle limit
!                    if ( (xyMin > Model%maxTail_buffer) .or. xyMin > Model%maxSun_buffer) then
!                        isCoreInactive(i,j) = .true.
!                    endif
!
!                    ! Now that we applied all criteria, update bndLoc
!                    if (isCoreInactive(i,j) .and. (i > bndLoc(j))) then
!                        bndLoc(j) = i
!                    endif
!                enddo
!            enddo
!        end subroutine calcCoreConstraints
!    end subroutine getInactiveCells


    subroutine getBufferCells(Model, Grid, State, isInactive, isBuffer, ocbDistO)
        !! Determines which cells should be marked as inactive
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isInactive
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: isBuffer
        integer, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), optional, intent(inout) :: ocbDistO
            !! Optional array to get ocbDist out of this function if desired

        integer :: i,j
        real(rp) :: xyMin
        integer, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: ocbDist
            !! Distance each cell is to open-closed boundary, up to nB + 1
        
        ! As a first pass, we will blanket enforce that there must always be nB buffer between any inactive and active cells
        ! We can do more complex stuff later

        call CalcOCBDist(Grid%shGrid, isInactive, Grid%nBufferCells, ocbDist)

        where (ocbDist > 0 .and. ocbDist .le. Grid%nBufferCells)
            isBuffer = .true.
        elsewhere
            isBuffer = .false.
        end where

        ! Also limit to max allowed sun/tail extent
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (ocbDist(i,j) > 0) then  !! OCBdist will be zero in inactive zone, we want to avoid editing those cells
                    ! Simple circle limit
                    if ( any(State%radMax_cc(i:i+1,j:j+1) > Model%maxTail_active) .or. any(State%radMax_cc(i:i+1,j:j+1) > Model%maxSun_active)) then
                        isBuffer(i,j) = .true.
                    endif
                endif
            enddo
        enddo

        if (present(ocbDistO)) then
            ocbDistO = ocbDist
        endif

    end subroutine getBufferCells


    subroutine CalcOCBDist(sh, isInactive, nBnd, ocbDist)
        !! Calculates distance each cell is from open/closed boundary, up to nBnd cells away
        type(ShellGrid_T), intent(in) :: sh
            !! RAIJU shell grid
        logical, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(in) :: isInactive
            !! Whether cell centers are inactive
        integer, intent(in) :: nBnd
            !! Number of desired layers between open/closed boundary and active domain
        integer, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(out) :: ocbDist


        integer :: iLayer, i, j, iL, iU, jL, jU


        where ( isInactive )
            ocbDist = 0
        elsewhere
            ocbDist = nBnd + 1
        end where


        ! Grow out from open/closed boundary, set proper distance for closed points
        do iLayer=1,nBnd
            !! Note: We can safely parallelize within eah iLayer loop if needed. But this isn't a bottleneck right now
            do i=sh%isg, sh%ieg
                do j=sh%jsg, sh%jeg
                    iL = max(i-1, sh%isg)
                    iU = min(i+1, sh%ieg)
                    jL = max(j-1, sh%jsg)
                    jU = min(j+1, sh%jeg)

                    if (isInactive(i,j)) then
                        cycle
                    !else if ( (ocbDist(i,j) .eq. nBnd+1) .and. any(ocbDist(iL:iU,jL:jU) .eq. iLayer-1) ) then  ! This version includes the diagonal cells, which are technically not adjacent to our cell
                    else if ( (ocbDist(i,j) .eq. nBnd+1) .and. ( i .eq. sh%isg .or. i .eq. sh%ieg .or. any(ocbDist(iL:iU,j) .eq. iLayer-1) .or. any(ocbDist(i,jL:jU) .eq. iLayer-1) ) ) then  ! Doesn't include diagonals
                        !! If current point's distance hasn't been decided and its bordering cell with iL-1, this point is distance iL from ocb
                        ocbDist(i,j) = iLayer
                    endif
                enddo
            enddo
        enddo

        call wrapJcc(sh, ocbDist)
        
    end subroutine CalcOCBDist


    subroutine calcMapJacNorm(Grid, bMin, jacNorm)
        type(raijuGrid_T), intent(in) :: Grid
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 3), intent(in) :: bMin
            !! xyz coordinates of bMin surface
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: jacNorm
            !! Cell-centered jacobian norm we calculate

        integer :: i,j, dim
        real(rp) :: ddTheta, ddPhi, aij
        jacNorm = 0.0

        ! Iterate through cell centers
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg

                ! Doing L 2,1 norm (https://en.wikipedia.org/wiki/Matrix_norm)
                ! Derivative by mapping cell corner values to edges, gradient across cell
                ! Deriv w.r.t. theta first:
                !ddTheta = 0.0
                !do dim = 1,3
                !    aij = 0.5*(bMin(i+1,j+1,dim) + bMin(i+1,j,dim) \
                !            - (bMin(i,j+1,dim) + bMin(i,j,dim))) /  Grid%lenFace(i,j,RAI_PH)
                !    ddTheta = ddtheta + aij**2
                !enddo
!
                !ddPhi = 0.0
                !do dim=1,3
                !    aij = 0.5*(bMin(i+1,j+1,dim) + bMin(i,j+1,dim) \
                !            - (bMin(i+1,j,dim) - bMin(i,j,dim))) / (0.5*(Grid%lenFace(i+1,j,RAI_TH) + Grid%lenFace(i,j,RAI_TH)))
                !    ddPhi = ddPhi + aij**2
                !enddo
                !jacNorm(i,j) = sqrt(ddTheta) + sqrt(ddPhi)
                
                
                do dim=1,3
                    aij = (0.5*(bMin(i+1,j+1,dim) + bMin(i+1,j,dim) \
                                - (bMin(i,j+1,dim) + bMin(i,j,dim))) /  Grid%lenFace(i,j,RAI_PH))**2
                    aij = aij + (0.5*(bMin(i+1,j+1,dim) + bMin(i,j+1,dim) \
                                - (bMin(i+1,j,dim) - bMin(i,j,dim))) / (0.5*(Grid%lenFace(i+1,j,RAI_TH) + Grid%lenFace(i,j,RAI_TH))))**2
                    jacNorm(i,j) = jacNorm(i,j) + sqrt(aij)
                    !jacNorm(i,j) = max(jacNorm(i,j), aij)
                enddo

            enddo
        enddo
    end subroutine calcMapJacNorm


    subroutine calcCornerNormAngles(Model, Grid, State, normAngle)
        !! For each cell corner, calculate the maximum angle between the normals of the 4 triangles its a part of
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(inout) :: normAngle

        integer :: i,j,d,u,v
        real(rp), dimension(3) :: v1, v2, v3, v4
        real(rp), dimension(4,3) :: crosses
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1,3) :: xyz_sm

        normAngle = 1.0

        ! First, smooth xyzmins
        xyz_sm = State%xyzMin(:,:,:)
        do i=Grid%shGrid%isg+1,Grid%shGrid%ieg  ! Ignore last row and column, no smoothing for you
            do j=Grid%shGrid%jsg+1,Grid%shGrid%jeg
                do d=XDIR,ZDIR
                    xyz_sm(i,j,d) = SmoothOperator33(State%xyzMin(i-1:i+1,j-1:j+1,d))
                enddo
            enddo
        enddo
        do i=Grid%shGrid%isg+1,Grid%shGrid%ieg
            do j=Grid%shGrid%jsg+1,Grid%shGrid%jeg
                ! Build vectors from corner i,j
                do d=XDIR,ZDIR
                    v1(d) = xyz_sm(i,j,d) - xyz_sm(i  ,j-1,d)
                    v2(d) = xyz_sm(i,j,d) - xyz_sm(i-1,j  ,d)
                    v3(d) = xyz_sm(i,j,d) - xyz_sm(i  ,j+1,d)
                    v4(d) = xyz_sm(i,j,d) - xyz_sm(i+1,j  ,d)
                enddo

                ! Do cross products
                crosses(1,:) = cross(v2, v1)
                crosses(2,:) = cross(v3, v2)
                crosses(3,:) = cross(v4, v3)
                crosses(4,:) = cross(v1, v4)
                crosses(1,:) = crosses(1,:)/max(norm2(crosses(1,:)), TINY)
                crosses(2,:) = crosses(2,:)/max(norm2(crosses(2,:)), TINY)
                crosses(3,:) = crosses(3,:)/max(norm2(crosses(3,:)), TINY)
                crosses(4,:) = crosses(4,:)/max(norm2(crosses(4,:)), TINY)

                ! Calculate min
                do u=1,3
                    do v=u+1,4
                        ! We are doing dot product, so 1 = zero angle and -1 means 180 deg angle
                        normAngle(i,j) = min(normAngle(i,j), dot_product(crosses(u,:), crosses(v,:)))
                    enddo
                enddo

            enddo
        enddo

    end subroutine calcCornerNormAngles


    subroutine forceActiveWithinRadius(Model, Grid, State, active)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: active

        integer :: i,j, iWalk
        integer :: iBufStart, iActStart
        real(rp) :: radL, radR

!        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
!            i = Grid%shGrid%is
!            
!            radL = State%radMax(i,j)
!            radR = State%radMax(i,j+1)
!            do while (i < Grid%shGrid%ie .and. (radL > Model%activeDomRad .or. radR > Model%activeDomRad))
!                i = i + 1
!                radL = State%radMax(i,j)
!                radR = State%radMax(i,j+1)
!            enddo
!
!            ! i is now the first cell that has corners under activeDomRad
!            ! If this point is active, we're good and have nothing to do for this j
!            if (State%active(i,j) .eq. RAIJUACTIVE) then
!                cycle
!            endif
!
!            ! If still here, something went wrong and active domain is very limited.
!            ! Try to save it by ignoring all but the most important constraint that the field line is closed
!            ! Move i nBuffer cells outward since we gotta add buffer cells and do the rest of our checks from there
!            i = max(Grid%shGrid%isg, i - Grid%nBufferCells)
!            do while ( i < Grid%shGrid%ie .and. any(State%topo(i:i+1,j:j+1) .eq. RAIJUOPEN) )
!                i = i + 1
!            enddo
!
!            ! Now i is a safe cell topology-wise, we can set active domain as we want from here
!            iBufStart = min(i, Grid%shGrid%ieg)
!            iActStart = min(iBufStart + Grid%nBufferCells, Grid%shGrid%ieg)
!
!            active(Grid%shGrid%isg:iBufStart-1, j) = RAIJUINACTIVE
!            if (iActStart > iBufStart) then
!                active(iBufStart:iActStart-1, j) = RAIJUBUFFER
!            endif
!            if (iActStart <= Grid%shGrid%ieg) then
!                active(iActStart:Grid%shGrid%ieg, j) = RAIJUACTIVE
!            endif
!        enddo

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,radL,radR,iWalk,iBufStart,iActStart)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%ie,Grid%shGrid%is,-1

                radL = State%radMax(i,j)
                radR = State%radMax(i,j+1)
                if (radL > Model%activeDomRad .or. radR > Model%activeDomRad) then
                    ! Out of our jurisdiction, pass on this j
                    exit
                endif

                ! Otherwise, inspect
                if (State%active(i,j) .eq. RAIJUACTIVE) then
                    ! This i,j is clear
                    cycle
                endif

                ! If still here, we are within activeDomRad and we ran into either an INACTIVE or BUFFER cell
                ! Keep walking out until either we hit bad topology or active would be placed within activeDomRad
                do while (iWalk > Grid%shGrid%is .or. radL > Model%activeDomRad .or. radR > Model%activeDomRad &
                          .or. any(State%topo(iWalk:iWalk+1,j:j+1) .eq. RAIJUINACTIVE))
                    iWalk = i - 1
                    radL = State%radMax(iWalk,j)
                    radR = State%radMax(iWalk,j+1)
                enddo

                if (iWalk == Grid%shGrid%is) then
                    ! We didn't find any good values above i,j don't rock the boat even harder
                    exit
                endif

                ! Now iWalk is one cell too far, make it the last usable cell
                iWalk = iWalk + 1
                ! Now iWalk is a safe cell topology-wise, we can set active domain as we want from here
                iBufStart = min(iWalk, Grid%shGrid%ieg)
                iActStart = min(iBufStart + Grid%nBufferCells, Grid%shGrid%ieg)
                active(Grid%shGrid%isg:iBufStart-1, j) = RAIJUINACTIVE
                if (iActStart > iBufStart) then
                    active(iBufStart:iActStart-1, j) = RAIJUBUFFER
                endif
                if (iActStart <= Grid%shGrid%ieg) then
                    active(iActStart:Grid%shGrid%ieg, j) = RAIJUACTIVE
                endif
                ! We are done for this j now
                exit
            enddo
        enddo

    end subroutine forceActiveWithinRadius

end module raijuDomain
