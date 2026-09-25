module raijuAdvancer
    use clocks
    use planethelper
    use earthhelper, only : GallagherXY
    use kai2geo

    use raijudefs
    use raijutypes
    use raijugrids
    use raijuetautils
    use raijuBCs
    use raijuPreAdvancer
    use raijuPostAdvancer
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
        call raijuAdvanceState(Model, Grid, State)
        call Toc("AdvanceState",.true.)

        call Tic("Post-Advance",.true.)
        call raijuPostAdvanceState(Model, Grid, State)
        call Toc("Post-Advance",.true.)


    end subroutine raijuAdvance


    subroutine raijuAdvanceState(Model, Grid, State)
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
    end subroutine raijuAdvanceState


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
                    call State%myRaiApp%WriteBlackBox()
                    stop
                endif
            enddo

            ! Post-advance bookkeeping for this specific eta channel

            if (Model%doLosses) then
                ! Divide losses and precip fluxes by big dt to turn them into proper rates
                ! Quantities going into the ionosphere divided by 2 for per-hemisphere output
                State%precipNFlux_Nk(k)%data = State%precipNFlux_Nk(k)%data/State%dt/2.0_rp
                State%precipEFlux_Nk(k)%data = State%precipEFlux_Nk(k)%data/State%dt/2.0_rp
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

        ! Low-lat theta ghosts feed eta_half and face reconstruction at the
        ! equatorward active edge, so impose their BC before flux calculation.
        call applyLowLatEtaBC(Grid, State, k)

        ! Calculate eta at half-step
        call calcEtaHalf(Grid, State, k, dt, odt)
        ! Save eta to eta_last
        State%eta_last(:,:,k) = State%eta(:,:,k)

        ! Calculate eta face fluxes
        call calcFluxes(Model, Grid, State, k, isGoodRecon, State%eta_half(:,:,k), Qflux)  ! [eta * Rp/s]

        ! Check to see if any fluxes need to be attenuated to avoid cell underflows
        call fixFluxes(Model, Grid, State%eta(:,:,k), Qflux, dt, isGoodEvol)
        
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

                !! TODO: Something better later
                State%eta(i,j,k) = max(State%eta(i,j,k),0.0_rp)
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

        ! Re-apply after updating active cells.
        call applyLowLatEtaBC(Grid, State, k)

    end subroutine stepLambda


    subroutine applyLowLatEtaBC(Grid, State, k)
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k

        integer :: i, j
        real(rp) :: nowKp, den, etaBnd, xeq, yeq

        associate(sh=>Grid%shGrid, spc=>Grid%spc(Grid%k2spc(k)))

        if (sh%ie+1 > sh%ieg) return

        if (spc%flav == F_PSPH) then
            ! Plasmasphere ghosts represent real lower-L shells, so reset each
            ! ghost row from Gallagher at its own mapped equatorial location.
            nowKp = State%KpTS%evalAt(State%t)
            do j=sh%js, sh%je
                do i=sh%ie+1, sh%ieg
                    xeq = State%xyzMincc(i,j,XDIR)
                    yeq = State%xyzMincc(i,j,YDIR)
                    den = GallagherXY(xeq, yeq, nowKp)
                    etaBnd = den*State%bvol_cc(i,j)*sclEta
                    State%eta     (i,j,k) = etaBnd
                    State%eta_last(i,j,k) = etaBnd
                enddo
            enddo
        else
            ! Hot/species channels use a zero-gradient low-lat eta boundary.
            do j=sh%js, sh%je
                State%eta     (sh%ie+1:sh%ieg,j,k) = State%eta     (sh%ie,j,k)
                State%eta_last(sh%ie+1:sh%ieg,j,k) = State%eta_last(sh%ie,j,k)
            enddo
        endif

        ! Keep periodic phi ghosts consistent before eta_half or face
        ! reconstruction sees the low-lat boundary values.
        call wrapJcc(sh, State%eta     (:,:,k))
        call wrapJcc(sh, State%eta_last(:,:,k))

        end associate

    end subroutine applyLowLatEtaBC


    subroutine fixFluxes(Model, Grid, eta, Qflux, dt, isGoodEvol)
        !! Positivity-preserving flux attenuation.  First compute how much
        !! conserved content each cell can afford to send out, then attenuate
        !! each shared face flux by the donor cell's alpha.
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: eta
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: Qflux
        real(rp), intent(in) :: dt
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isGoodEvol

        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: alpha
        real(rp) :: botConFlux
        real(rp) :: cellCon
        real(rp) :: conWeight
        real(rp) :: faceConFlux
        real(rp) :: leftConFlux
        real(rp) :: outCon
        real(rp) :: rightConFlux
        real(rp) :: topConFlux
        integer :: i,j

        alpha = 1.0_rp

        ! Cell pass: bound each cell's outgoing conserved-content transport.
        ! Here "conserved content" means the quantity that the update actually
        ! conserves: eta*areaCC*BrCC in the RAIJU form, or eta*areaCC in
        ! the RCM-like branch.  The limiter is deliberately pessimistic:
        ! it does not spend incoming flux before it arrives.  That keeps the
        ! rule local and easy to audit, at the price of adding diffusion when a
        ! cell is trying very hard to empty itself.  The numerics take a sip of
        ! tea, not a payday loan.
        do j=Grid%shGrid%js,Grid%shGrid%je
            do i=Grid%shGrid%is,Grid%shGrid%ie
                if (.not. isGoodEvol(i,j)) then
                    cycle
                endif

                if (Model%doHack_rcmEtaPush) then
                    ! Convert each face flux to the conserved-content transport
                    ! appearing in the doHack update.  Using this same measure
                    ! here and below is what keeps the limiter aligned with the
                    ! discrete update instead of with a prettier continuum story.
                    conWeight = Grid%areaCC(i,j)
                    leftConFlux  = Qflux(i  ,j  ,RAI_TH)*Grid%lenFace(i  ,j  ,RAI_TH)/Grid%BrFace(i  , j  , RAI_TH)
                    rightConFlux = Qflux(i+1,j  ,RAI_TH)*Grid%lenFace(i+1,j  ,RAI_TH)/Grid%BrFace(i+1, j  , RAI_TH)
                    botConFlux   = Qflux(i  ,j  ,RAI_PH)*Grid%lenFace(i  ,j  ,RAI_PH)/Grid%BrFace(i  , j  , RAI_PH)
                    topConFlux   = Qflux(i  ,j+1,RAI_PH)*Grid%lenFace(i  ,j+1,RAI_PH)/Grid%BrFace(i  , j+1, RAI_PH)
                else
                    ! Standard RAIJU content: eta is cell-centered, but the
                    ! finite-volume balance is for eta*area*Br.  Face transports
                    ! already carry BrFace from calcFluxes, so multiplying by
                    ! lenFace gives the signed content rate through each face.
                    conWeight = Grid%areaCC(i,j)*Grid%BrCC(i,j)
                    leftConFlux  = Qflux(i  ,j  ,RAI_TH)*Grid%lenFace(i  ,j  ,RAI_TH)
                    rightConFlux = Qflux(i+1,j  ,RAI_TH)*Grid%lenFace(i+1,j  ,RAI_TH)
                    botConFlux   = Qflux(i  ,j  ,RAI_PH)*Grid%lenFace(i  ,j  ,RAI_PH)
                    topConFlux   = Qflux(i  ,j+1,RAI_PH)*Grid%lenFace(i  ,j+1,RAI_PH)
                endif

                ! Match the signs in the eta update:
                !   + left  - right  + bottom  - top.
                ! A left/bottom face removes content when its signed transport
                ! is negative; a right/top face removes content when it is
                ! positive.  That is the "diode" part: only fluxes leaving this
                ! cell count against this cell's available content.
                cellCon = max(0.0_rp, eta(i,j)*conWeight)
                outCon = dt*(max(0.0_rp, -leftConFlux ) &
                            + max(0.0_rp,  rightConFlux) &
                            + max(0.0_rp, -botConFlux  ) &
                            + max(0.0_rp,  topConFlux  ))
                if (outCon > cellCon) then
                    alpha(i,j) = cellCon/outCon
                endif
            enddo
        enddo

        ! Face pass: scale each face once, using the donor cell's alpha.
        ! This is the conservation-critical part.  We never separately adjust
        ! the two cells sharing a face; we edit the one shared flux, so whatever
        ! one cell loses is still exactly what the neighbor gains.
        do j=Grid%shGrid%js,Grid%shGrid%je
            do i=Grid%shGrid%is,Grid%shGrid%ie+1
                if (Model%doHack_rcmEtaPush) then
                    faceConFlux = Qflux(i,j,RAI_TH)*Grid%lenFace(i,j,RAI_TH)/Grid%BrFace(i,j,RAI_TH)
                else
                    faceConFlux = Qflux(i,j,RAI_TH)*Grid%lenFace(i,j,RAI_TH)
                endif

                ! RAI_TH face i sits between cells i-1 and i.  Positive signed
                ! transport enters cell i from i-1, so i-1 is the donor; negative
                ! transport leaves cell i toward i-1, so i is the donor.
                if (faceConFlux > 0.0_rp) then
                    if (i > Grid%shGrid%is) then
                        if (isGoodEvol(i-1,j)) Qflux(i,j,RAI_TH) = alpha(i-1,j)*Qflux(i,j,RAI_TH)
                    endif
                else if (faceConFlux < 0.0_rp) then
                    if (i <= Grid%shGrid%ie) then
                        if (isGoodEvol(i,j)) Qflux(i,j,RAI_TH) = alpha(i,j)*Qflux(i,j,RAI_TH)
                    endif
                endif
            enddo
        enddo

        do j=Grid%shGrid%js,Grid%shGrid%je+1
            do i=Grid%shGrid%is,Grid%shGrid%ie
                if (Model%doHack_rcmEtaPush) then
                    faceConFlux = Qflux(i,j,RAI_PH)*Grid%lenFace(i,j,RAI_PH)/Grid%BrFace(i,j,RAI_PH)
                else
                    faceConFlux = Qflux(i,j,RAI_PH)*Grid%lenFace(i,j,RAI_PH)
                endif

                ! RAI_PH is the same story in the periodic longitude direction:
                ! face j is between j-1 and j.  The guard checks avoid asking
                ! the limiter to borrow alpha from outside the active loop; phi
                ! ghosts/wrapping are handled by the surrounding grid machinery.
                if (faceConFlux > 0.0_rp) then
                    if (j > Grid%shGrid%js) then
                        if (isGoodEvol(i,j-1)) Qflux(i,j,RAI_PH) = alpha(i,j-1)*Qflux(i,j,RAI_PH)
                    endif
                else if (faceConFlux < 0.0_rp) then
                    if (j <= Grid%shGrid%je) then
                        if (isGoodEvol(i,j)) Qflux(i,j,RAI_PH) = alpha(i,j)*Qflux(i,j,RAI_PH)
                    endif
                endif
            enddo
        enddo

        ! The periodic phi interface appears as both j=js and j=je+1 in the
        ! update stencil. Keep the duplicate face identical using the copy that
        ! was limited by the donor cell in the pass above.
        do i=Grid%shGrid%is,Grid%shGrid%ie
            if (Qflux(i,Grid%shGrid%je+1,RAI_PH) > 0.0_rp) then
                Qflux(i,Grid%shGrid%js,RAI_PH) = Qflux(i,Grid%shGrid%je+1,RAI_PH)
            else if (Qflux(i,Grid%shGrid%js,RAI_PH) < 0.0_rp) then
                Qflux(i,Grid%shGrid%je+1,RAI_PH) = Qflux(i,Grid%shGrid%js,RAI_PH)
            endif
        enddo

    end subroutine fixFluxes


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
                State%eta_half(i,j,k) = max(State%eta(i,j,k) + 0.5_rp*dt*dQ,0.0_rp)
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
                Req = State%radMax_cc(i,j)
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
