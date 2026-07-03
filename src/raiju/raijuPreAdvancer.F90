module raijuPreAdvancer
    use clocks
    use planethelper
    use kai2geo

    use raijudefs
    use raijutypes
    use raijuetautils
    use raijuBCs
    use raijulosses
    use raijuRecon

    use shellGrid
    use shellInterp

    use arrayutil

    implicit none

    contains


!------
! Main high-level functions
!------

    subroutine raijuPreAdvance(Model, Grid, State)
        !! Takes a state and calculates what is needed in order to advance
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: k

        ! Clear things that will be accumulated over the advance
        call fillArray(State%eta_avg, 0.0_rp)
        ! (losses handled in updateRaiLosses)

        ! Now that topo is set, we can calculate active domain
        call setActiveDomain(Model, Grid, State)

        ! Moments to etas, initial active shell calculation
        call Tic("BCs")
        call applyRaijuBCs(Model, Grid, State, doWholeDomainO=State%isFirstCpl) ! If fullEtaMap=True, mom2eta map is applied to the whole domain
        if (State%isFirstCpl) then
            call setRaijuInitPsphere(Model, Grid, State, Model%psphInitKp)
        endif
        call Toc("BCs")

        ! Handle plasmasphere refilling for the full step about to happen
        call plasmasphereRefill(Model,Grid,State)

        ! Handle edge cases that may effect the validity of information carried over from last coupling period
        call prepEtaLast(Grid%shGrid, State, State%isFirstCpl)

        call Tic("Calc cell-centered velocities")
        call calcPotGrads_cc(Model, Grid, State)
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call calcVelocity_cc(Model, Grid, State, k, State%cVel(:,:,k,:))  ! Get velocity at cell interfaces
            ! Calc sub-time step. Each channel will do this on its own, but this way we can output the step sizes everyone is using
            call reconVelocityLRs(Model, Grid, State, k, State%iVelL(:,:,k,:), State%iVelR(:,:,k,:))
            State%dtk(k) = activeDt_LR(Model, Grid, State, k)
        enddo
        call Toc("Calc cell-centered velocities")

        ! Some loss processes depends on up-to-date densities, so we should run EvalMoments first
        call Tic("Moments Eval PreAdvance")
        call EvalMoments(Grid, State)
        call Toc("Moments Eval PreAdvance")

        call Tic("Losses")
        if (Model%doLosses) then
            call Tic("Update loss states")
            call updateRaiLosses(Model, Grid, State)
            call Toc("Update loss states")
            
            call Tic("Calc loss rates")
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(k)
            do k=1,Grid%Nk
                call calcChannelLossRates(Model, Grid, State, k)
            enddo
            call Toc("Calc loss rates")
        endif
        call Toc("Losses")

    end subroutine raijuPreAdvance


!------
! Transitions between coupling chunks
!------

    subroutine prepEtaLast(sh, State, isFirstCpl)
        ! To be called before main advance loop
        ! Edits the initial eta_last state to account for edge cases
        type(ShellGrid_T), intent(in) :: sh
        type(raijuState_T), intent(inout) :: State
        logical, intent(in) :: isFirstCpl

        integer :: i,j

        ! If its the first coupling instance, we have no good previous information
        ! Set previous cpl info to current state
        if (isFirstCpl) then
            State%active_last = State%active
            State%eta_last    = State%eta
            return
        endif

        ! For now, don't try to be fancy, just reset everyone
        State%active_last = State%active
        State%eta_last    = State%eta
        return

        ! For a given cell, if it is currently active, but wasn't active last step, 
        !   we need to set its eta_last to something usable in the first eta_half calculation
        ! We do this in the buffer region as well, because this was freshly populated with MHD moments and any existing eta_last is stale
        ! After the first dt, eta_last will be set with valid information and no further adjustments are necessay
        !do j=sh%jsg,sh%jeg
        !    do i=sh%isg,sh%ieg
        !        if ( (State%active(i,j) .eq. RAIJUACTIVE .and. State%active_last(i,j) .ne. RAIJUACTIVE) &
        !        .or. (State%active(i,j) .eq. RAIJUBUFFER) ) then
        !            State%eta_last(i,j,:) = State%eta(i,j,:)
        !        endif
        !    enddo
        !enddo

        ! NOTE: Any case not addressed here should mean that eta_last and active_last are valid in those cases

    end subroutine prepEtaLast


!------
! Cell Potentials and their gradients
!------
    !! Potential variables come in pre-allocated, since they may be subsets of larger arrays
    subroutine potExB(sh, State, pExB, doSmoothO, isGCornerO)
        type(ShellGrid_T), intent(in) :: sh
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(inout) :: pExB
        logical, intent(in), optional :: doSmoothO
        logical, dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(in), optional :: isGCornerO

        logical :: doSmooth

        if(present(doSmoothO)) then
            doSmooth = doSmoothO
        else
            doSmooth = .false.
        endif
        
        pExB = State%espot * 1.e3  ! [kV -> V]

        if (doSmooth .and. present(isGCornerO)) then 
            call smoothCornerVar(sh, isGcornerO, pExB)
        endif

    end subroutine potExB

    subroutine potCorot(planet, sh, pCorot, doGeoCorotO)
        ! Simple corotation potential [V]
        type(planet_T), intent(in) :: planet
        type(ShellGrid_T), intent(in) :: sh
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(inout) :: pCorot
        logical, intent(in), optional :: doGeoCorotO
        
        integer :: i,j

        if (present(doGeoCorotO)) then
            if (doGeoCorotO) then
                pCorot = 0.0
                do j=sh%jsg,sh%jeg+1
                    do i=sh%isg,sh%ieg+1
                        call geocorotation_from_SMTP(sh%th(i), sh%ph(j), pCorot(i,j))
                    enddo
                enddo
                ! geopack corotation returns in kV, we need V
                pCorot = pCorot * 1.e3  ! [kV -> V]
                return
            endif
        endif

        ! IfdoGeoCorotO was not provided, or it was false, we default to corotation on aligned dipole and rotational axis
        do j=sh%jsg,sh%jeg+1
            pCorot(:,j) = -planet%psiCorot*(planet%rp_m/planet%ri_m)*sin(sh%th)**2 * 1.e3  ! [kV -> V]
        enddo

    end subroutine potCorot


    subroutine potGC(shGrid, alamc, bVol, pGC)
        ! Simple gradient curvature potential [V]
        type(ShellGrid_T), intent(in) :: shGrid
        real(rp), intent(in) :: alamc
        real(rp), dimension(shGrid%isg:shGrid%ieg+1,&
                            shGrid%jsg:shGrid%jeg+1), intent(in) :: bVol

        real(rp), dimension(shGrid%isg:shGrid%ieg+1,&
                            shGrid%jsg:shGrid%jeg+1), intent(inout) :: pGC

        pGC = alamc*bVol**(-2./3.)

    end subroutine potGC


    subroutine calcEffectivePotential(Model, Grid, State, pEff)
        !! Note: This is not used to calculate velocities in the solver itself
        !! Calculates effective potential [V] for all lambda channels as the sum of pExB, pCorot, and pGC
        type(raijuModel_T) , intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1,&
                            Grid%Nk), intent(inout) :: pEff
        
        integer :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: pExB, pCorot, pGC
        
        pEff = 0.0

        ! Grab 2D potentials
        call potExB(Grid%shGrid, State, pExB)
        call potCorot(Model%planet, Grid%shGrid, pCorot, Model%doGeoCorot)
        
        ! Build 3D effective potential
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(k, pGC)
        do k=1,Grid%Nk
            call potGC(Grid%shGrid, Grid%alamc(k), State%bVol, pGC)
            pEff(:,:,k) = pExB + pCorot + pGC
        enddo

    end subroutine calcEffectivePotential


    subroutine calcPotGrads_cc(Model, Grid, State)
        !! Calculate cell-centered potential gradients using Green-Gauss method
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State


        integer :: i,j
        real(rp) :: bvol_cc
        ! Cell corners
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: isGCorner
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: pExB, pCorot
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradVM

        State%cVel = 0.0

        where (State%topo .eq. RAIJUCLOSED)
            isGCorner = .true.
        elsewhere
            isGCorner = .false.
        end where

        associate(sh=>Grid%shGrid)
        ! Gauss-Green calculation of cell-averaged gradients
        call potExB(Grid%shGrid, State, pExB, doSmoothO=Model%doSmoothGrads, isGCornerO=isGCorner)  ! [V]
        call potCorot(Model%planet, Grid%shGrid, pCorot, Model%doGeoCorot)  ! [V]
        call calcGradIJ_cc(Model%planet%rp_m, Grid, isGCorner, pExB  , State%gradPotE_cc    , doLimO=.true. )  ! [V/m]
        call calcGradIJ_cc(Model%planet%rp_m, Grid, isGCorner, pCorot, State%gradPotCorot_cc, doLimO=.false.)  ! [V/m]

        ! GC drifts depend on grad(lambda * V^(-2/3))
        ! lambda is constant, so just need grad(V^(-2/3) )
        call calcGradVM_cc(Model%planet%rp_m, Model%planet%ri_m, Model%planet%magMoment, &
                            Grid, isGCorner, State%bvol, State%gradVM_cc, &
                            doSmoothO=Model%doSmoothGrads, doLimO=.true.)
        end associate
        
    end subroutine calcPotGrads_cc


    subroutine calcGradIJ_cc(Rp_m, Grid, isG, Q, gradQ, doLimO)
        !! Uses Green-Gauss theorem to get cell-averaged gradient, using corner-located values
        real(rp), intent(in) :: Rp_m
            !! Planet radius in meters
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: isG
            !! Mask for corners that are safe to use in calculating the gradient across the attached faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: Q
            !! Variable we are taking the gradient of across faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradQ
        logical, optional, intent(in) :: doLimO

        real(rp), dimension(:,:,:), allocatable :: gradQtmp
        integer :: i,j
        real(rp) :: qLow, qHigh 
        logical :: doLim
        
        if (present(doLimO)) then
            doLim = doLimO
        else
            doLim = .true.
        endif

        gradQ = 0.0

        associate(sh=>Grid%shGrid)

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,qLow,qHigh)
        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg
                if (all(isG(i:i+1, j:j+1))) then
                    ! Theta-direction gradient
                    !qLow  = Grid%lenFace(i  ,j,RAI_TH)/2.0 * (Q(i  ,j+1) + Q(i  ,j))
                    !qHigh = Grid%lenFace(i+1,j,RAI_TH)/2.0 * (Q(i+1,j+1) + Q(i+1,j))
                    !gradQ(i,j,RAI_TH) = (qHigh - qLow) / (Grid%areaCC(i,j) )
                    qLow  = 0.5 * (Q(i  ,j+1) + Q(i  ,j))
                    qHigh = 0.5 * (Q(i+1,j+1) + Q(i+1,j))
                    gradQ(i,j,RAI_TH) = (qHigh - qLow) / (Grid%lenFace(i,j,RAI_PH) )  ! [Q/Rp]

                    ! Phi-direction gradient
                    !qLow  = Grid%lenFace(i,j  ,RAI_PH)/2.0 * (Q(i+1,j  ) + Q(i, j  ))
                    !qHigh = Grid%lenFace(i,j+1,RAI_PH)/2.0 * (Q(i+1,j+1) + Q(i, j+1))
                    !gradQ(i,j,RAI_PH) = (qhigh - qLow) / (Grid%areaCC(i,j) )
                    qLow  = 0.5 * (Q(i+1,j  ) + Q(i, j  ))
                    qHigh = 0.5 * (Q(i+1,j+1) + Q(i, j+1))
                    gradQ(i,j,RAI_PH) = (qhigh - qLow) / (0.5*(Grid%lenFace(i,j,RAI_TH)+Grid%lenFace(i+1,j,RAI_TH)) )  ! [Q/Rp]
                endif
            enddo
        enddo

        if (doLim) then
            allocate(gradQtmp(sh%isg:sh%ieg,sh%jsg:sh%jeg, 2))
            gradQtmp = gradQ    
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j)
            do j=sh%jsg+1,sh%jeg-1
                do i=sh%isg+1,sh%ieg-1
                    gradQ(i,j,RAI_TH) = MCLim(gradQtmp(i-1:i+1,j      ,RAI_TH), isG(i-1:i+1,j      ))
                    gradQ(i,j,RAI_PH) = MCLim(gradQtmp(i      ,j-1:j+1,RAI_PH), isG(i      ,j-1:j+1))
                enddo
            enddo
            !deallocate(gradQtmp)
        endif

        call wrapJcc(Grid%shGrid, gradQ(:,:,RAI_TH))
        call wrapJcc(Grid%shGrid, gradQ(:,:,RAI_PH))

        gradQ = gradQ / Rp_m  ! [Q/m]

        end associate


        contains

        function MCLim(dq, isG_3) result(dqbar)
            real(rp), dimension(3), intent(in) :: dq
            logical, dimension(3), intent(in) :: isG_3
            real(rp) :: dqbar
            real(rp), dimension(2) :: dqmask
            real(rp) :: magdq
            real(rp) :: dvL, dvR, dvC
            
            dvC = dq(2)
            dvL = 0.5*(dq(1)+dq(2))
            dvR = 0.5*(dq(2)+dq(3))

            dqbar = 0.0

            if (all(isG_3)) then
                ! Standard MCLim
                !if (dq(1)*dq(3) <= 0.0) then
                if (dvL*dvR <= 0.0) then
                    dqbar = 0.0
                else
                    magdq = min(2*abs(dvL),2*abs(dvR),abs(dvC))
                    dqbar = sign(magdq,dvC)
                    !magdq = min(2*abs(dq(1)),2*abs(dq(3)),abs(dq(2)))
                    !!SIGN(A,B) returns the value of A with the sign of B
                    !dqbar = sign(magdq,dq(2))
                endif
            else
                dqmask(2) = dvC
                if (isG_3(1) .and. isG_3(2)) then
                    dqmask(1) = dvL
                else if (isG_3(3) .and. isG_3(2)) then
                    dqmask(1) = dvR
                endif
!
                if (dqmask(1)*dqmask(2) <= 0.0) then
                    dqbar = 0.0
                else
                    magdq = min(2*abs(dqmask(1)), abs(dqmask(2)))
                    dqbar = sign(magdq, dqmask(2))
                endif
                ! Trying one-sided lim
!                dqmask(2) = dq(2)
!                if (isG_3(1) .and. isG_3(2)) then
!                    dqmask(1) = dq(1)
!                else if (isG_3(3) .and. isG_3(2)) then
!                    dqmask(1) = dq(3)
!                endif
!
!                if (dqmask(1)*dqmask(2) <= 0.0) then
!                    dqbar = 0.0
!                else
!                    magdq = min(2*abs(dqmask(1)), abs(dqmask(2)))
!                    dqbar = sign(magdq, dqmask(2))
!                endif
            endif
        end function MCLim

    end subroutine calcGradIJ_cc


    subroutine calcGradVM_cc(Rp_m, RIon_m, B0, Grid, isGcorner, V, gradVM, doSmoothO, doLimO)
        !! Calculates the cell-centered gradient of the quantity vm = FTV^(-2/3)
        real(rp), intent(in) :: Rp_m
            !! Planet radius [m]]
        real(rp), intent(in) :: RIon_m
            !! Ionosphere radius [m]
        real(rp), intent(in) :: B0
            !! Planet's surface field strength [Gauss]
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: isGcorner
            !! Mask for corners that are safe to use in calculating the gradient across the attached faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: V
            !! Flux tube volume (FTV)
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradVM
            !! grad(FTV) we return [units(FTV)/m]
        logical, optional, intent(in) :: doSmoothO, doLimO

        integer :: i,j
        real(rp) :: dcl_dm, bVolcc
            !! d colat / d meters, used to convert dipole derivative w.r.t. colat to meters
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: V0, dV
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: dV0_dth
            !! V0 = dipole FTV, dV = V - V0
        logical :: doSmooth, doLim

        if (present(doSmoothO)) then
            doSmooth = doSmoothO
        else
            doSmooth = .true.
        endif
        if (present(doLimO)) then
            doLim = doLimO
        else
            doLim = .true.
        endif

        associate(sh=>Grid%shGrid)

        do i=sh%isg, sh%ieg+1
            ! DipFTV_colat takes a colat at 1 Rp, make sure we use that value instead of colat in ionosphere (shGrid%th)
            V0(i,:) = DipFTV_colat(Grid%thRp(i), B0)
        enddo

        dV = 0.0
        where(isGcorner)
            dV = V - V0
        end where


        ! Analytic gradient of dipole FTV at cell centers
        dV0_dth = 0.0
        dcl_dm = 1.0/RIon_m
        do i=sh%isg,sh%ieg
            dV0_dth(i,:) = DerivDipFTV(Grid%thcRp(i), B0) * dcl_dm
        enddo

        ! Gradient of perturbation
        if (doSmooth) then
            call smoothCornerVar(Grid%shGrid, isGcorner, dV)
        endif

        ! Building gradVM term by term
        ! grad(V^(-2/3)) = -2/3*V^(-5/3) * grad(V)
        gradVM = 0.0
        call calcGradIJ_cc(Rp_m, Grid, isGcorner, dV, gradVM, doLimO=doLim)
        gradVM(:,:,RAI_TH) = gradVM(:,:,RAI_TH) + dV0_dth
        
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,bVolcc)
        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg
                if(all(isGcorner(i:i+1,j:j+1))) then
                    !bVolcc = toCenter2D(dV(i:i+1,j:j+1)) + DipFTV_colat(Grid%thcRp(i), B0)  ! Will include smoothing of dV if enabled
                    bVolcc = toCenter2D(V(i:i+1,j:j+1))
                    gradVM(i,j,:) = (-2./3.)*bVolcc**(-5./3.)*gradVM(i,j,:)
                else
                    ! gradVM should be zero for this point coming out of calcGradIJ_cc, but set to dipole value just in case
                    gradVM(i,j,RAI_PH) = 0.0
                    gradVM(i,j,RAI_TH) = 0.0
                    !gradVM(i,j,RAI_TH) = (-2./3.)*DipFTV_colat(Grid%thcRp(i), B0)**(-5./3.)*dV0_dth(i,j)
                endif
            enddo
        enddo

        end associate

    end subroutine calcGradVM_cc

    subroutine smoothCornerVar(sh, isGc, V)
        type(ShellGrid_T), intent(in) :: sh
        logical , dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(in) :: isGc
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(inout) :: V
        
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1) :: Vsm
        real(rp), dimension(3,3) :: tmpV
        logical , dimension(3,3) :: tmpGood
        integer :: i,j

        Vsm = 0.0
        associate (isg=>sh%isg, ieg=>sh%ieg, jsg=>sh%jsg, jeg=>sh%jeg)
        ! Handle cells along grid extents first
        ! isg,jsg corner
        tmpV = 0.0
        tmpGood = .false.
        tmpV   (2:3,2:3) = V   (isg:isg+1,jsg:jsg+1)
        tmpGood(2:3,2:3) = isGc(isg:isg+1,jsg:jsg+1)
        Vsm(isg,jsg) = SmoothOperator33(tmpV, tmpGood)
        ! ieg,jsg corner
        tmpV = 0.0
        tmpGood = .false.
        tmpV   (1:2,2:3) = V   (ieg:ieg+1,jsg:jsg+1)
        tmpGood(1:2,2:3) = isGc(ieg:ieg+1,jsg:jsg+1)
        Vsm(ieg+1,jsg) = SmoothOperator33(tmpV, tmpGood)
        ! isg,jeg corner
        tmpV = 0.0
        tmpGood = .false.
        tmpV   (2:3,1:2) = V   (isg:isg+1,jeg:jeg+1)
        tmpGood(2:3,1:2) = isGc(isg:isg+1,jeg:jeg+1)
        Vsm(isg,jeg+1) = SmoothOperator33(tmpV, tmpGood)
        ! ieg,jeg corner
        tmpV = 0.0
        tmpGood = .false.
        tmpV   (1:2,1:2) = V   (ieg:ieg+1,jeg:jeg+1)
        tmpGood(1:2,1:2) = isGc(ieg:ieg+1,jeg:jeg+1)
        Vsm(ieg+1,jeg+1) = SmoothOperator33(tmpV, tmpGood)
        ! jsg, jeg edges
        do i=isg+1,ieg
            tmpV = 0.0
            tmpGood = .false.
            tmpV   (:,2:3) = V   (i-1:i+1,jsg:jsg+1)
            tmpGood(:,2:3) = isGc(i-1:i+1,jsg:jsg+1)
            Vsm(i,jsg) = SmoothOperator33(tmpV, tmpGood)

            tmpV = 0.0
            tmpGood = .false.
            tmpV   (:,1:2) = V   (i-1:i+1,jeg:jeg+1)
            tmpGood(:,1:2) = isGc(i-1:i+1,jeg:jeg+1)
            Vsm(i,jeg+1) = SmoothOperator33(tmpV, tmpGood)
        enddo
        ! isg,ieg edges
        do j=jsg+1,jeg
            tmpV = 0.0
            tmpGood = .false.
            tmpV   (2:3,:) = V   (isg:isg+1,j-1:j+1)
            tmpGood(2:3,:) = isGc(isg:isg+1,j-1:j+1)
            Vsm(isg,j) = SmoothOperator33(tmpV, tmpGood)

            tmpV = 0.0
            tmpGood = .false.
            tmpV   (1:2,:) = V   (ieg:ieg+1,j-1:j+1)
            tmpGood(1:2,:) = isGc(ieg:ieg+1,j-1:j+1)
            Vsm(ieg+1,j) = SmoothOperator33(tmpV, tmpGood)
        enddo
        ! Now everyone else
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j)
        do j=jsg+1,jeg
            do i=isg+1,ieg
                Vsm(i,j) = SmoothOperator33(V(i-1:i+1,j-1:j+1), isGc(i-1:i+1,j-1:j+1))
            enddo
        enddo

        call wrapJcorners(sh, Vsm)

        ! Write back to provided array
        V = Vsm
        end associate

    end subroutine smoothCornerVar

!------
! Velocity calculations
!------

    subroutine calcVelocity_cc(Model, Grid, State, k, Vtp, gradVMO)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: Vtp
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2), optional, intent(in) :: gradVMO
            !! Optional gradVM to use in place of State's gradVM, in case you wanna do testing on gradVM calculation
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradPot
        integer :: i,j
        
        ! Ionospheric ExB and gradient-curvature drifts
        if (present(gradVMO)) then
            gradPot = State%gradPotE_cc + Grid%alamc(k)*gradVMO
        else
            gradPot = State%gradPotE_cc + Grid%alamc(k)*State%gradVM_cc
        endif

        ! If doing our own corotation, assuming ExB pot didn't include it already. Add it ourselves
        if (Model%doOwnCorot) then
            gradPot = gradPot + State%gradPotCorot_cc
        endif

        Vtp(:,:,RAI_TH) =      gradPot(:,:,RAI_PH) / (Grid%Brcc*1.0e-9)  ! [m/s]
        Vtp(:,:,RAI_PH) = -1.0*gradPot(:,:,RAI_TH) / (Grid%Brcc*1.0e-9)  ! [m/s]

    end subroutine calcVelocity_cc


    subroutine reconVelocityLRs(Model, Grid, State, k, iVelL, iVelR)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: iVelL, iVelR

        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGCC
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1,2) :: tmpVelL, tmpVelR
        integer :: i,j
        
        iVelL = 0.0
        iVelR = 0.0
        isGCC = .false.
        ! The good mask here indicates if a cell can be used in reconstruction, which for us includes ACTIVE and BUFFER zones
        where (State%active .ne. RAIJUINACTIVE)
            isGCC = .true.
        end where

        ! ReconFaces is gonna take one cc value and reconstruct at cell faces
        ! But phi velocity at theta edge is meaningless for us
        ! So we make two temporary arrays, and then just pack the meaningful components into iVelL and iVelR at the end
        call ReconFaces(Model, Grid, isGCC, State%cVel(:,:,k,RAI_TH), tmpVelL, tmpVelR, doOMPO=.false.)
        iVelL(:,:,RAI_TH) = tmpVelL(:,:,RAI_TH)
        iVelR(:,:,RAI_TH) = tmpVelR(:,:,RAI_TH)
        call ReconFaces(Model, Grid, isGCC, State%cVel(:,:,k,RAI_PH), tmpVelL, tmpVelR, doOMPO=.false.)
        iVelL(:,:,RAI_PH) = tmpVelL(:,:,RAI_PH)
        iVelR(:,:,RAI_PH) = tmpVelR(:,:,RAI_PH)

    end subroutine reconVelocityLRs

!------
! time handling
!------

    function activeDt_LR(Model, Grid, State, k) result(dt)
        !! Calculates min dt needed to safele evolve active domain for given energy invariant channel k
        !! TODO: Consider dynamic CFL factor
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        
        integer :: i,j
        real(rp) :: dl, velij_th, velij_ph
        !real(rp), dimension(Grid%shGrid%is:Grid%shGrid%ie+1,Grid%shGrid%js:Grid%shGrid%je+1) :: dtArr
        real(rp), dimension(Grid%shGrid%isg+1:Grid%shGrid%ieg,Grid%shGrid%jsg+1:Grid%shGrid%jeg) :: dtArr
        real(rp) :: dt

        associate (sh => Grid%shGrid)

        dtArr = HUGE

        !!$OMP PARALLEL DO default(shared) &
        !!$OMP schedule(dynamic) &
        !!$OMP private(i,j, velij_th, velij_ph, dl)
        do j=sh%jsg+1,sh%jeg
            do i=sh%isg+1,sh%ieg  ! Excluding first and last face
                velij_th = TINY
                velij_ph = TINY
                ! NOTE: We are only checking faces bordering non-ghost cells because those are the only ones we use for evolution
                ! TODO: Strictly speaking, there are some faces that are included here that shouldn't be because they're never used to evolve anything
                !  so we should make the loop js:je, is:ie, and handle the last row and column afterwards

                ! We are responsible for face (i,j,TH) and (i,j,PH)
                ! Theta faces first

                ! This is a weird pattern. Basically, we can't cycle if the overall desired condition isn't met, because we are doing theta and phi directions in the same loop
                ! At the same time, we don't want to write one massive if condition with all options covered, or a bunch of nested if statements
                ! So instead, we break them up, such that if a single if condition is true it means we have failed the physical condition for us to calculate a valid timestep
                ! If all if conditions fail then the physical conditions are met and we can do our calculations in the else block
                if (Model%doActiveShell .and. ( .not. State%activeShells(i-1,k) .or. .not. State%activeShells(i,k)) ) then
                    ! In order for a theta-dir face to be usable, we need both sides to be active
                    continue
                else if ( State%active(i-1,j) .eq. RAIJUINACTIVE .or. State%active(i,j) .eq. RAIJUINACTIVE ) then
                    continue
                else
                    ! We are good lets calculate a dt for this face
                    velij_th = max(abs(State%iVelL(i,j,k,RAI_TH)), abs(State%iVelR(i,j,k,RAI_TH)), TINY)  ! [m/s]
                    !dtArr(i,j,RAI_TH) = ( Grid%delTh(i) * Model%planet%ri_m ) / velij  ! [s]
                endif
                
                ! Phi faces
                if (Model%doActiveShell .and. .not. State%activeShells(i,k)) then
                    ! In order for a phi-dir face to be usable, we just need this i shell to be active
                    continue
                else if (State%active(i,j-1) .eq. RAIJUINACTIVE  .or. State%active(i,j) .eq. RAIJUINACTIVE ) then 
                    continue
                else
                    velij_ph = max(abs(State%iVelL(i,j,k,RAI_PH)), abs(State%iVelR(i,j,k,RAI_PH)), TINY)
                    !dtArr(i,j,RAI_PH) = ( Grid%delPh(j) * Model%planet%ri_m ) / velij  ! [s]
                endif

                dl = abs(min(Grid%delTh(i), Grid%delPh(j)*sin(sh%thc(i))))
                dtArr(i,j) = (dl*Model%planet%ri_m) / sqrt(velij_th**2 + velij_ph**2)

            enddo
        enddo

        dt = Model%CFL*minval(dtArr)

        end associate
            
    end function activeDt_LR

end module raijuPreAdvancer