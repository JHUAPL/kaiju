module raijuAdvancer
    use clocks
    use planethelper
    use kai2geo

    use raijudefs
    use raijutypes
    use raijugrids
    use raijuetautils
    use raijuBCs
    use raijuPreAdvancer
    use raijulosses
    use raijuRecon

    implicit none

    contains

!------
! Advance entry point
!------

    subroutine raijuAdvance(Model, Grid, State, dtCpl)
        !! Controls entirety of eta evolution over time dtCpl
        !! Assumes that any coupling setup has been completed
        !! Calculates velocities and dt, evolves all etas over 
        !! dtCpl, and accumulates certain quantities like precipitation info
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        real(rp), intent(in) :: dtCpl

        integer :: k

        State%dt = dtCpl

        call Tic("Pre-Advance",.true.)
        call raijuPreAdvance(Model, Grid, State)
        call Toc("Pre-Advance",.true.)
        State%isFirstCpl = .false.

        ! Step
        call Tic("AdvanceState",.true.)
        call AdvanceState(Model, Grid, State)
        call Toc("AdvanceState",.true.)

        ! etas back to moments
        call Tic("Moments Eval",.true.)
        call EvalMoments(Grid, State)
        call EvalMoments(Grid, State, doAvgO=.true.)
        call Toc("Moments Eval",.true.)

    end subroutine raijuAdvance


    subroutine AdvanceState(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: k

        ! Send each channel off on their own
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call AdvanceLambda(Model, Grid, State, k)
        enddo

        ! Wrap up, save new time and cpl step number
        ! Save current active domain to active_last
        State%active_last = State%active
        State%eta_last = State%eta  ! Also happens in AdvanceLambda, but we do it again here just to be sure
    end subroutine AdvanceState


    subroutine AdvanceLambda(Model, Grid, State, k)
        !! Advances a single lambda channel from State time t to t+dt
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k

        integer :: n, s , Nmax
            !! n = step counter
            !! s = species index
        integer :: i,j
            !! theta, phi loop indices
        real(rp) :: t, dt, odt, tEnd
            !! t, dt = Current time and delta-t for this channel
            !! odt = previous time step
            !! tEnd = Time we stop at
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGoodEvol
            !! Flag for whether we should be evolving certain cells
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGoodRecon
            !! Flag for whether certain cells should be used for reconstruction

        ! If we are plasmasphere channel and are not evolving plasmasphere, nothing to do
        if (.not. Model%doPsphEvol .and. Grid%spc(Grid%k2spc(k))%flav==F_PSPH) then
            return
        endif

        ! Initial settings
        s = Grid%k2spc(k)
        t = State%t
        tEnd = State%t + State%dt
        dt  = State%dtk(k)
        odt = State%dtk(k)

        ! NOTE: When we actually implement activeShells these will need to be updated within time loop
        if (Grid%spc(Grid%k2spc(k))%flav==F_PSPH) then
            call calcEvolDom_psph(Model, Grid, State, isGoodEvol)
        else
            where (State%active .eq. RAIJUACTIVE)
                isGoodEvol = .true.
            elsewhere
                isGoodEvol = .false.
            end where
        endif

        where (State%active .ne. RAIJUINACTIVE)
            isGoodRecon = .true.
        elsewhere
            isGoodRecon = .false.
        end where


        associate(sh=>Grid%shGrid, spc=>Grid%spc(s))

            ! Here we go!
            n = 0  ! counter
            do while ( tEnd-t > TINY)
                
                ! Calc new active shells
                !if (Model%doActiveShell) then
                !    State%activeShells(:,k) = setLambdaActiveShells(sh, spc, State%bVol, &
                !                                State%eta(:,:,spc%kStart:spc%kEnd), k, worthyFracO = Model%worthyFrac)
                !endif
                

                ! Calc next time step
                !dt = activeDt(Model, Grid, State, k)
                dt = activeDt_LR(Model, Grid, State, k)

                ! If needed, reduce dt to fit within remaining time
                ! This way, all channels will end exactly at tEnd
                if (t + dt > tEnd) then
                    dt = max(tEnd - t, TINY)  ! Make sure we never go back in time and advance at least a little bit so we will eventually end
                endif

                ! Advection
                call stepLambda(Model, Grid, State, k, dt, odt, isGoodEvol, isGoodRecon)
                
                ! Losses
                if (Model%doLosses) then
                    call applyLosses(Model, Grid, State, k, dt)
                endif

                ! Accumulate eta_avg
                ! Explicitly do loop here to help with speed
                do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                    do i=Grid%shGrid%isg,Grid%shGrid%ieg
                        ! Weighted by this step's dt to handle variable step sizes
                        State%eta_avg(i,j,k) = State%eta_avg(i,j,k) + dt*State%eta(i,j,k)
                    enddo
                enddo

                t = t + dt
                n = n+1
                odt = dt

                if (n/State%dt > Model%maxItersPerSec) then
                    write(*,*)"ERROR: Too many advance steps. Dying"
                    write(*,*)" k, dt, nSteps= ", k, dt, n
                    stop
                endif
            enddo

            ! Post-advance bookkeeping for this specific eta channel

            if (Model%doLosses) then
                ! Divide losses and precip fluxes by big dt to turn them into proper rates
                ! Quantities going into the ionosphere divided by 2 for per-hemisphere output
                State%precipNFlux(k)%data = State%precipNFlux(k)%data/State%dt/2.0_rp
                State%precipEFlux(k)%data = State%precipEFlux(k)%data/State%dt/2.0_rp
                State%CCHeatFlux(:,:,k)  = State%CCHeatFlux (:,:,k)/State%dt/2.0_rp
                State%dEta_dt(:,:,k)     = State%dEta_dt    (:,:,k)/State%dt
            endif
            State%eta_avg(:,:,k) = State%eta_avg(:,:,k)/State%dt

            State%nStepk(k) = State%nStepk(k) + n
            State%eta_last(:,:,k) = State%eta(:,:,k)

        end associate

    end subroutine AdvanceLambda


    subroutine stepLambda(Model, Grid, State, k, dt, odt, isGoodEvol, isGoodRecon)
        !! Advances a single lambda channel a single time step
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt
        real(rp), intent(in) :: odt
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isGoodEvol
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isGoodRecon

        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: Qflux
        integer :: i,j

        ! Calculate eta at half-step
        call calcEtaHalf(Grid, State, k, dt, odt)
        ! Save eta to eta_last
        State%eta_last(:,:,k) = State%eta(:,:,k)

        ! Calculate eta face fluxes
        call calcFluxes(Model, Grid, State, k, isGoodRecon, State%eta_half(:,:,k), Qflux)  ! [eta * Rp/s]

        !! This probably doesn't need to be at this scope. 
        !! Same with isG inside calcFluxes; can be pushed at least 1 level up

        ! Calc new eta
        do j=Grid%shGrid%js,Grid%shGrid%je
            do i=Grid%shGrid%is,Grid%shGrid%ie
                if (.not. isGoodEvol(i,j)) then
                    cycle
                endif

                if (Model%doHack_rcmEtaPush) then
                    ! Do eta conservation exactly 
                    ! (still finite volume for space but ignore different magnetic field between cells)
                    State%eta(i,j,k) = State%eta(i,j,k) + dt/Grid%areaCC(i,j) &
                                                    * ( Qflux(i  ,j  ,RAI_TH)*Grid%lenFace(i  ,j  ,RAI_TH)/Grid%BrFace(i  , j  , RAI_TH) &
                                                      - Qflux(i+1,j  ,RAI_TH)*Grid%lenFace(i+1,j  ,RAI_TH)/Grid%BrFace(i+1, j  , RAI_TH) &
                                                      + Qflux(i  ,j  ,RAI_PH)*Grid%lenFace(i  ,j  ,RAI_PH)/Grid%BrFace(i  , j  , RAI_PH) &
                                                      - Qflux(i  ,j+1,RAI_PH)*Grid%lenFace(i  ,j+1,RAI_PH)/Grid%BrFace(i  , j+1, RAI_PH) )
                else
                    State%eta(i,j,k) = State%eta(i,j,k) + dt/Grid%areaCC(i,j)/Grid%BrCC(i,j) &
                                                    * ( Qflux(i  ,j  ,RAI_TH)*Grid%lenFace(i  ,j  ,RAI_TH) &
                                                      - Qflux(i+1,j  ,RAI_TH)*Grid%lenFace(i+1,j  ,RAI_TH) &
                                                      + Qflux(i  ,j  ,RAI_PH)*Grid%lenFace(i  ,j  ,RAI_PH) &
                                                      - Qflux(i  ,j+1,RAI_PH)*Grid%lenFace(i  ,j+1,RAI_PH) )
                endif

                !if (State%eta(i,j,k) < -1.0*TINY) then
                !    write(*,*)"RAIJU ERROR: eta below zero!"
                !    write(*,*)"i,j,k=",i,j,k
                !    write(*,*)"eta=",State%eta(i,j,k)
                !    write(*,*)"eta_half=",State%eta_half(i,j,k)
                !    write(*,*)"eta_last=",State%eta_last(i,j,k)
                !    write(*,*)"Theta fluxes:",Qflux(i  ,j  ,RAI_TH),Qflux(i+1,j  ,RAI_TH)
                !    write(*,*)"Phi fluxes:",Qflux(i  ,j  ,RAI_PH),Qflux(i,j+1,RAI_PH)
                !    stop
                !endif

            enddo
        enddo

        call wrapJcc(Grid%shGrid, State%eta(:,:,k))

    end subroutine stepLambda


    subroutine calcEtaHalf(Grid, State, k, dt, odt)
        !! Estimate eta half step, currently using AB2
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt
        real(rp), intent(in) :: odt

        integer :: i,j
        real(rp) :: dQ

        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                dQ = (State%eta(i,j,k) - State%eta_last(i,j,k))/odt
                State%eta_half(i,j,k) = State%eta(i,j,k) + 0.5_rp*dt*dQ
            enddo
        enddo

    end subroutine calcEtaHalf


    subroutine calcEvolDom_psph(Model, Grid, State, isGoodEvol)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State

        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: isGoodEvol

        integer :: i,j
        real(rp) :: Req

        isGoodEvol = .false.  ! Can't use fillArray yet, no version for logicals

        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                Req = sqrt(State%xyzMincc(i,j,XDIR)**2 + State%xyzMincc(i,j,YDIR)**2)
                if (Req < Model%psphEvolRad) then
                    continue  ! Keep points within this radius set to no evol 
                else
                    if (State%active(i,j) == RAIJUACTIVE) then
                        isGoodEvol(i,j) = .true.
                    endif
                endif
            enddo
        enddo

    end subroutine

end module raijuAdvancer
