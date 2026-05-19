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
    use raijuMagHelpers

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
        ! If we should, convert some ion energy into new electrons
        if (Model%doIon2EleEnergyTransfer) then
            call ion2eleEnergyTransfer(Model, Grid, State)
        endif

        ! Handle edge cases that may effect the validity of information carried over from last coupling period
        call prepEtaLast(Grid%shGrid, State, State%isFirstCpl)

        call Tic("Calc velocities")
        call calcInterfaceVelocities(Model, Grid, State)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(k)
        do k=1,Grid%Nk
            ! Diagnostic cell-centered velocity from the face velocities that drive the solver.
            call calcVelocity_cc(Grid, State, k, State%cVel(:,:,k,:))
            ! Calc sub-time step. Each channel will do this on its own, but this way we can output the step sizes everyone is using.
            State%dtk(k) = activeDt_LR(Model, Grid, State, k)
        enddo
        call Toc("Calc velocities")

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
! Pre-advance cross-flavor, cross-lambda operations
!------

    subroutine ion2eleEnergyTransfer(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: i,j
            !! Fraction of total ion pressure we transver over tscl seconds

        call fillArray(State%P_transfer, 0.0_rp)
        associate(sh=>Grid%shGrid)
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j)
            do j=sh%jsg,sh%jeg
                do i=sh%isg,sh%ieg
                    call ion2eleEnergyTransfer_ij(Model, Grid, State, i, j)
                enddo
            enddo
        end associate


        contains

        subroutine ion2eleEnergyTransfer_ij(Model, Grid, State, i, j)
            type(raijuModel_T), intent(in) :: Model
            type(raijuGrid_T ), intent(in) :: Grid
            type(raijuState_T), intent(inout) :: State
            integer :: i,j

            integer :: k,s
            real(rp), dimension(Grid%nSpc) :: Pion_spc, Dion_spc, Pion_new_spc
            logical,  dimension(Grid%nSpc) :: isIonFlav
            real(rp) :: Pion_total, Dion_total, Pele_total, P_transfer, kT_MaxwellEle
            real(rp) :: Bmin_cc, tscl, dt_over_tscl, vm
            real(rp) :: etaMaxwell_k
            real(rp) :: fracPressTransfer = 0.05
            logical  :: doDipDriftPeriod = .true.
            real(rp) :: pressThresh = 0.01  ! [nPa] you must be at least this tall to get on the ride

            if (State%active(i,j) == RAIJUINACTIVE) then
                ! Note: we are applying this to both buffer and active cells, worth seeing if we should remove application in buffer
                return
            endif

            vm = State%bvol_cc(i,j)**(-2./3.)
            Bmin_cc = 0.25*sum( norm2(State%Bmin(i:i+1,j:j+1,:), dim=3) )
            Pion_spc = 0.0
            Dion_spc = 0.0
            P_transfer = 0.0
            isIonFlav = .false.

            ! First get total pressures of all the ions
            do s=1,Grid%nSpc
                associate(spc=>Grid%spc(s))
                if ((spc%spcType == RAIJUHPLUS .and. spc%flav .ne. F_PSPH ) .or. (spc%spcType==RAIJUOPLUS) ) then
                    isIonFlav(s) = .true.
                    Pion_spc(s) = SpcEta2Press(spc, State%eta(i,j,spc%kStart:spc%kEnd), State%bvol_cc(i,j))
                    Dion_spc(s) = SpcEta2Den  (spc, State%eta(i,j,spc%kStart:spc%kEnd), State%bvol_cc(i,j))
                else if (spc%flav == F_HOTE) then
                    ! This is the electron flavor, grab its pressure while we're here
                    Pele_total = SpcEta2Press(spc, State%eta(i,j,spc%kStart:spc%kEnd), State%bvol_cc(i,j))
                endif
                end associate
            enddo

            Pion_total = sum(Pion_spc)
            Dion_total = sum(Dion_spc)
            if (Pion_total < pressThresh .or. Dion_total < TINY) then
                ! No one here to create electrons, move on
                return
            endif
            ! Also get timescale, for now use 1 full GC drift period
            if (doDipDriftPeriod) then
                tscl = DipoleDriftPeriod(Dion_total, Pion_total, State%radMax_cc(i,j))
            else
                tscl = DriftPeriod(Dion_total, Pion_total, State%radMax_cc(i,j), Bmin_cc, State%rCurv(i,j), Model%planet%rp_m)
            endif
            dt_over_tscl = min(1.0_rp, State%dt/tscl)  ! Don't take bigger step than 1 
            if (State%dt/tscl > 1.0) then
                write(*,*)"dtCpl bigger than drift timescale for point i,j=",i,j
                write(*,*)"State%dt=",State%dt,"tscl=",tscl
            endif
            ! This is our very crude approach to calculating how much ion pressure/energy should go to electrons. We can put in data-informed calculation later
            P_transfer = fracPressTransfer*Pion_total  ! Transfer 5% of ion pressure over one drift timescale (we add timescale later)

            ! If Pele_total < Pion_total, then we're gonna do stuff. Otherwise there's enough electron pressure here and we will move on
            if (Pele_total >= P_transfer) then
                return
            endif

            ! Only transfer what we actually need to get Pele to P_transfer
            P_transfer = max(P_transfer - Pele_total, 0.0)  ! max just to make sure we can't transfer negative pressure, even though above check should catch it
            if (P_transfer <= 0.0) return  ! Once more for good measure. I really don't wanna touch etas using bad pressure values
            ! Finalize electron maxwellian we'll be adding
            kT_MaxwellEle = fracPressTransfer*DP2kT(Dion_total, Pion_total)  ! 5% of the ion's bulk temperature
            
            State%P_transfer(i,j) = P_transfer
            State%driftPeriod(i,j) = tscl

            ! Ok lets get to work
            do s=1,Grid%nSpc
                associate(spc=>Grid%spc(s))
                if (isIonFlav(s)) then
                    ! Update Pion_spc so its now the target pressure we want to adjust the ion flavors to
                    Pion_new_spc(s) = Pion_spc(s) - Pion_spc(s)/Pion_total * P_transfer * dt_over_tscl ! P_transfer is overall pressure transfer, Pion_spc/Pion_total is fraction this flavor contributs to that total P_transfer, dtCpl/tscl is how much we do this coupling/advance step
                    if (Pion_new_spc(s) < 0.0) then
                        write(*,*)"Pion_spc(s) target < 0 for s=",s
                        write(*,*)"Initial Pion_spc(s)=",Pion_new_spc(s)
                        write(*,*)"Pressure frac=",Pion_spc(s)/Pion_total
                        write(*,*)"P_transfer=",P_transfer
                        write(*,*)"dt_over_tscl=",dt_over_tscl
                    endif
                    call rescaleEtaDP(spc, State%eta(i,j,spc%kStart:spc%kEnd), Dion_spc(s), Pion_new_spc(s), State%bvol_cc(i,j), allowOneMomentRescaleO=.false.)
                else if (spc%flav == F_HOTE) then
                    ! We are at electron flavor, add in our Maxwellian
                    do k=spc%kStart,spc%kEnd
                        etaMaxwell_k = Maxwell2Eta(Model, Dion_total, kT_MaxwellEle, vm, abs(spc%alami(k)), abs(spc%alami(k+1)))
                        if (State%eta(i,j,k) < etaMaxwell_k) then
                            State%eta(i,j,k) = State%eta(i,j,k) + State%dt/tscl*etaMaxwell_k
                        endif
                    enddo
                endif
                end associate
            enddo
        end subroutine ion2eleEnergyTransfer_ij

    end subroutine ion2eleEnergyTransfer


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


    subroutine calcCornerBackgrounds(Model, Grid, State, pCorot, ftv0)
        !! Build unsmoothed corner-centered background fields for the new
        !! face-native velocity path.
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(inout) :: pCorot
            !! Corotation potential [V]
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(inout) :: ftv0
            !! Dipole flux tube volume [Rp/nT]

        integer :: i,j

        if (Model%doOwnCorot) then
            call potCorot(Model%planet, Grid%shGrid, pCorot, Model%doGeoCorot)
        else
            pCorot = State%pot_corot * 1.0e3_rp  ! [kV -> V]
        endif

        associate(sh=>Grid%shGrid)
        do j=sh%jsg,sh%jeg+1
            do i=sh%isg,sh%ieg+1
                ! DipFTV_colat expects colatitude mapped to 1 Rp, not the
                ! ionospheric shell colatitude.
                ftv0(i,j) = DipFTV_colat(Grid%thRp(i), Model%planet%magMoment)
            enddo
        enddo
        end associate

    end subroutine calcCornerBackgrounds


    subroutine calcCornerPerturbations(Model, Grid, State, pCorot, ftv0, dPotEB, logFtvRatio)
        !! Split total corner fields into perturbations plus analytic/background
        !! pieces. In coupled mode State%espot is already potential_total; in
        !! standalone mode RAIJU owns corotation and State%espot excludes it.
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: pCorot
            !! Corotation potential [V]
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: ftv0
            !! Dipole flux tube volume [Rp/nT]
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(inout) :: dPotEB
            !! Electric-potential perturbation [V]
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(inout) :: logFtvRatio
            !! Logarithmic flux-tube-volume perturbation, log(bvol/ftv0)

        dPotEB = State%espot * 1.0e3  ! [kV -> V]
        if (.not. Model%doOwnCorot) then
            dPotEB = dPotEB - pCorot
        endif

        logFtvRatio = 0.0_rp
        where (State%bvol > TINY .and. ftv0 > TINY)
            logFtvRatio = log(State%bvol/ftv0)
        end where

    end subroutine calcCornerPerturbations


    subroutine smoothCornerPerturbations(Grid, State, dPotEB, logFtvRatio)
        !! Smooth only the perturbation pieces. The analytic/background fields
        !! are deliberately not passed through this window.
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(inout) :: dPotEB
            !! Electric-potential perturbation [V]
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(inout) :: logFtvRatio
            !! Logarithmic flux-tube-volume perturbation, log(bvol/ftv0)

        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: isClosedCorner
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: isFtvCorner

        where (State%topo .eq. RAIJUCLOSED)
            isClosedCorner = .true.
        elsewhere
            isClosedCorner = .false.
        end where

        isFtvCorner = isClosedCorner .and. State%bvol > TINY

        ! Smoothing disabled for face-native velocity test runs.
        !call smoothCornerVar(Grid%shGrid, isClosedCorner, dPotEB)
        !call smoothCornerVar(Grid%shGrid, isFtvCorner, logFtvRatio)

    end subroutine smoothCornerPerturbations


    subroutine calcCornerPotPieces(Model, Grid, State)
        !! Build the two corner potential pieces used by the face-native
        !! velocity calculation without storing an Nk-sized effective potential.
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: pCorot
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: ftv0, ftv
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: dPotEB, logFtvRatio

        integer :: i,j

        call calcCornerBackgrounds(Model, Grid, State, pCorot, ftv0)
        call calcCornerPerturbations(Model, Grid, State, pCorot, ftv0, dPotEB, logFtvRatio)
        call smoothCornerPerturbations(Grid, State, dPotEB, logFtvRatio)

        State%pot_corot = pCorot * 1.0e-3_rp
        State%pot_eb = dPotEB + pCorot
        ftv = ftv0*exp(logFtvRatio)
        State%pot_ftv = 0.0_rp

        associate(sh=>Grid%shGrid)
        do j=sh%jsg,sh%jeg+1
            do i=sh%isg,sh%ieg+1
                if (State%topo(i,j) .ne. RAIJUCLOSED) cycle
                if (State%bvol(i,j) <= TINY) cycle
                if (ftv(i,j) <= TINY) cycle
                State%pot_ftv(i,j) = ftv(i,j)**(-2.0_rp/3.0_rp)
            enddo
        enddo

        if (sh%isPeriodic) then
            call wrapJcorners(sh, State%pot_eb)
            call wrapJcorners(sh, State%pot_ftv)
        endif
        end associate

    end subroutine calcCornerPotPieces


    subroutine smoothCornerVar(sh, isGc, V)
        type(ShellGrid_T), intent(in) :: sh
        logical , dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(in) :: isGc
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(inout) :: V
        
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1) :: Vsm
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1) :: Vwrk
        logical , dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1) :: isGwrk
        real(rp), dimension(3,3) :: tmpV
        logical , dimension(3,3) :: tmpGood
        integer :: i,j

        Vwrk = V
        isGwrk = isGc
        if (sh%isPeriodic) then
            call wrapJcorners(sh, Vwrk)
            call wrapJcorners(sh, isGwrk)
        endif

        Vsm = 0.0
        associate (isg=>sh%isg, ieg=>sh%ieg, jsg=>sh%jsg, jeg=>sh%jeg)
        ! Handle cells along grid extents first
        ! isg,jsg corner
        tmpV = 0.0
        tmpGood = .false.
        tmpV   (2:3,2:3) = Vwrk   (isg:isg+1,jsg:jsg+1)
        tmpGood(2:3,2:3) = isGwrk(isg:isg+1,jsg:jsg+1)
        Vsm(isg,jsg) = SmoothOperator33_goodCenter(tmpV, tmpGood)
        ! ieg,jsg corner
        tmpV = 0.0
        tmpGood = .false.
        tmpV   (1:2,2:3) = Vwrk   (ieg:ieg+1,jsg:jsg+1)
        tmpGood(1:2,2:3) = isGwrk(ieg:ieg+1,jsg:jsg+1)
        Vsm(ieg+1,jsg) = SmoothOperator33_goodCenter(tmpV, tmpGood)
        ! isg,jeg corner
        tmpV = 0.0
        tmpGood = .false.
        tmpV   (2:3,1:2) = Vwrk   (isg:isg+1,jeg:jeg+1)
        tmpGood(2:3,1:2) = isGwrk(isg:isg+1,jeg:jeg+1)
        Vsm(isg,jeg+1) = SmoothOperator33_goodCenter(tmpV, tmpGood)
        ! ieg,jeg corner
        tmpV = 0.0
        tmpGood = .false.
        tmpV   (1:2,1:2) = Vwrk   (ieg:ieg+1,jeg:jeg+1)
        tmpGood(1:2,1:2) = isGwrk(ieg:ieg+1,jeg:jeg+1)
        Vsm(ieg+1,jeg+1) = SmoothOperator33_goodCenter(tmpV, tmpGood)
        ! jsg, jeg edges
        do i=isg+1,ieg
            tmpV = 0.0
            tmpGood = .false.
            tmpV   (:,2:3) = Vwrk   (i-1:i+1,jsg:jsg+1)
            tmpGood(:,2:3) = isGwrk(i-1:i+1,jsg:jsg+1)
            Vsm(i,jsg) = SmoothOperator33_goodCenter(tmpV, tmpGood)

            tmpV = 0.0
            tmpGood = .false.
            tmpV   (:,1:2) = Vwrk   (i-1:i+1,jeg:jeg+1)
            tmpGood(:,1:2) = isGwrk(i-1:i+1,jeg:jeg+1)
            Vsm(i,jeg+1) = SmoothOperator33_goodCenter(tmpV, tmpGood)
        enddo
        ! isg,ieg edges
        do j=jsg+1,jeg
            tmpV = 0.0
            tmpGood = .false.
            tmpV   (2:3,:) = Vwrk   (isg:isg+1,j-1:j+1)
            tmpGood(2:3,:) = isGwrk(isg:isg+1,j-1:j+1)
            Vsm(isg,j) = SmoothOperator33_goodCenter(tmpV, tmpGood)

            tmpV = 0.0
            tmpGood = .false.
            tmpV   (1:2,:) = Vwrk   (ieg:ieg+1,j-1:j+1)
            tmpGood(1:2,:) = isGwrk(ieg:ieg+1,j-1:j+1)
            Vsm(ieg+1,j) = SmoothOperator33_goodCenter(tmpV, tmpGood)
        enddo
        ! Now everyone else
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j)
        do j=jsg+1,jeg
            do i=isg+1,ieg
                Vsm(i,j) = SmoothOperator33_goodCenter(Vwrk(i-1:i+1,j-1:j+1), isGwrk(i-1:i+1,j-1:j+1))
            enddo
        enddo

        if (sh%isPeriodic) then
            call wrapJcorners(sh, Vsm)
        endif

        ! Write back to provided array
        V = Vsm
        end associate

        contains

        function SmoothOperator33_goodCenter(Q, isG) result(Qs)
            real(rp), dimension(3,3), intent(in) :: Q
            logical , dimension(3,3), intent(in) :: isG
            real(rp) :: Qs

            if (isG(2,2)) then
                Qs = SmoothOperator33(Q, isG)
            else
                Qs = 0.0_rp
            endif

        end function SmoothOperator33_goodCenter

    end subroutine smoothCornerVar

!------
! Velocity calculations
!------

    subroutine calcInterfaceVelocities(Model, Grid, State)
        !! Calculate the single face-centered normal velocity used by the solver.
        !! The electric and FTV pieces are differenced at faces separately, then
        !! combined per lambda channel as v = vEB + lambda*vFTV.
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: vEB, vFTV
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: isFtvCorner
        integer :: k

        call calcCornerPotPieces(Model, Grid, State)
        isFtvCorner = State%topo == RAIJUCLOSED .and. State%bvol > TINY .and. State%pot_ftv > TINY
        call calcPotentialDriftFaces(Model, Grid, State, State%pot_eb , vEB )
        call calcPotentialDriftFaces(Model, Grid, State, State%pot_ftv, vFTV, isGoodCornerO=isFtvCorner)
        call applyLowLatFTVCutoff(Model, Grid, vFTV)

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(k)
        do k=1,Grid%Nk
            State%iVel(:,:,k,:) = vEB + Grid%alamc(k)*vFTV
        enddo

    end subroutine calcInterfaceVelocities


    subroutine applyLowLatFTVCutoff(Model, Grid, vFTV)
        !! Below the MHD inner boundary, RAIJU should keep only the electric
        !! drift piece.  The old cell-centered path did this by adding the
        !! lambda/FTV term only for cells i <= iLowLat_Evol.
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: vFTV

        associate(sh=>Grid%shGrid)
        if (Model%iLowLat_Evol < sh%ieg+1) then
            vFTV(Model%iLowLat_Evol+1:sh%ieg+1,:,RAI_TH) = 0.0_rp
        endif
        if (Model%iLowLat_Evol < sh%ieg) then
            vFTV(Model%iLowLat_Evol+1:sh%ieg+1,:,RAI_PH) = 0.0_rp
        endif
        end associate

    end subroutine applyLowLatFTVCutoff


    subroutine calcPotentialDriftFaces(Model, Grid, State, pot, vFace, isGoodCornerO)
        !! Convert one corner potential piece to its face-normal drift velocity
        !! coefficient. For pot_eb this is the full ExB/corotation drift. For
        !! pot_ftv this is the velocity per unit lambda.
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: pot
            !! Corner potential piece
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: vFace
            !! Face-normal drift velocity, or velocity coefficient for pot_ftv
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in), optional :: isGoodCornerO
            !! Optional validity mask for potential-specific corner requirements.

        integer :: i,j
        real(rp) :: dPot_ds, denom
        logical :: isGood1, isGood2

        vFace = 0.0_rp

        associate(sh=>Grid%shGrid)
        ! Theta-normal faces span the phi direction, so use the phi corner difference.
        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg+1
                if (present(isGoodCornerO)) then
                    isGood1 = isGoodCornerO(i,j)
                    isGood2 = isGoodCornerO(i,j+1)
                else
                    isGood1 = State%topo(i,j) .eq. RAIJUCLOSED
                    isGood2 = State%topo(i,j+1) .eq. RAIJUCLOSED
                endif
                if (.not. isGood1 .or. .not. isGood2) cycle
                denom = Grid%lenFace(i,j,RAI_TH)*Model%planet%rp_m
                if (abs(denom) <= TINY .or. abs(Grid%BrFace(i,j,RAI_TH)) <= TINY) cycle
                dPot_ds = (pot(i,j+1) - pot(i,j))/denom
                vFace(i,j,RAI_TH) = dPot_ds/(Grid%BrFace(i,j,RAI_TH)*1.0e-9_rp)
            enddo
        enddo

        ! Phi-normal faces span the theta direction, so use the theta corner difference.
        do j=sh%jsg,sh%jeg+1
            do i=sh%isg,sh%ieg
                if (present(isGoodCornerO)) then
                    isGood1 = isGoodCornerO(i,j)
                    isGood2 = isGoodCornerO(i+1,j)
                else
                    isGood1 = State%topo(i,j) .eq. RAIJUCLOSED
                    isGood2 = State%topo(i+1,j) .eq. RAIJUCLOSED
                endif
                if (.not. isGood1 .or. .not. isGood2) cycle
                denom = Grid%lenFace(i,j,RAI_PH)*Model%planet%rp_m
                if (abs(denom) <= TINY .or. abs(Grid%BrFace(i,j,RAI_PH)) <= TINY) cycle
                dPot_ds = (pot(i+1,j) - pot(i,j))/denom
                vFace(i,j,RAI_PH) = -dPot_ds/(Grid%BrFace(i,j,RAI_PH)*1.0e-9_rp)
            enddo
        enddo
        end associate

    end subroutine calcPotentialDriftFaces


    subroutine calcVelocity_cc(Grid, State, k, Vtp)
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: Vtp

        integer :: i,j

        Vtp = 0.0_rp

        ! Diagnostic cell-centered velocity from the face velocities that drive
        ! the solver.  A cell owns theta faces i/i+1 and phi faces j/j+1.
        associate(sh=>Grid%shGrid)
        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg
                Vtp(i,j,RAI_TH) = 0.5_rp*(State%iVel(i,j,k,RAI_TH) + State%iVel(i+1,j,k,RAI_TH))
                Vtp(i,j,RAI_PH) = 0.5_rp*(State%iVel(i,j,k,RAI_PH) + State%iVel(i,j+1,k,RAI_PH))
            enddo
        enddo
        end associate

    end subroutine calcVelocity_cc


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
                    velij_th = max(abs(State%iVel(i,j,k,RAI_TH)), TINY)  ! [m/s]
                    !dtArr(i,j,RAI_TH) = ( Grid%delTh(i) * Model%planet%ri_m ) / velij  ! [s]
                endif
                
                ! Phi faces
                if (Model%doActiveShell .and. .not. State%activeShells(i,k)) then
                    ! In order for a phi-dir face to be usable, we just need this i shell to be active
                    continue
                else if (State%active(i,j-1) .eq. RAIJUINACTIVE  .or. State%active(i,j) .eq. RAIJUINACTIVE ) then 
                    continue
                else
                    velij_ph = max(abs(State%iVel(i,j,k,RAI_PH)), TINY)
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
