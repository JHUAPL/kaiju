module mhdgroup
    use gamtypes
    use clocks
    use gamutils
    use gridutils
    use stress
    use fields
    use ringav
    use multifluid

    implicit none
    logical, parameter, private :: doBorisGrav = .false.
    
    contains

    subroutine AdvanceMHD(Model,Grid,State,oState,ooState,Solver,dt)
        type(Model_T)    , intent(inout) :: Model
        type(Grid_T)     , intent(inout) :: Grid
        type(State_T)    , intent(inout) :: State,oState,ooState
        type(gamSolver_T), intent(inout) :: Solver
        real(rp), intent(in) :: dt

        integer :: n

        !Trap for case (at beginning of run) where oState has same time as state
        if (abs(oState%time-State%time)<TINY) then
            oState%time = State%time-Model%dt
        endif
        !Use predictor to create half state
        call Tic("Predictor")
        call Predictor(Model,Grid,State,oState,ooState,Solver%StateHf,0.5*dt)
        call Toc("Predictor")

        !Get electric field for MHD case
        if (Model%doMHD) then
            !NOTE, CalcElecField and CalcStress are independent of each other
            call Tic("E-Field")
            call CalcElecField(Model,Grid,State,Solver%StateHf,Solver%Vf,State%Efld)
            !Call user hack function if defined
            if (associated(Model%HackE)) then
                call Tic("HackE")
                call Model%HackE(Model,Grid,State)
                call Toc("HackE")
            endif
            call Toc("E-Field")
            
        endif
        
        !Get plasma stresses
        !Send dGasM in hydro/MHD, but unallocated w/ hydro
        !Calculate stresses for all species
        call Tic("Stress")
        call CalcStress(Model,Grid,State,Solver%StateHF,Solver%gFlx,Solver%mFlx,Solver%dGasH,Solver%dGasM)
        call Toc("Stress")

        !Finalize, apply stresses and save State->oState for next predictor step
        call Tic("Update")
        call applyUpdate(Model,Grid,State,oState,ooState,dt,Solver%dGasH,Solver%dGasM,State%Efld)
        call Toc("Update")

        !Apply gravity if necessary (before ring avg)
        if (Model%doGrav) then
            call Tic("Gravity")
            call applyGrav(Model,Grid,State,oState,ooState,dt)
            call Toc("Gravity")
        endif

        !Ring average if necessary
        !Ring average configuration and this grid contains part of singularity
        call Tic("Ring-Average")
        if ( Model%doRing .and. (Model%Ring%doS .or. Model%Ring%doE) ) then
            call RingAverage(Model,Grid,State)
        endif
        call Toc("Ring-Average")

    end subroutine AdvanceMHD
    
    !Updates plasma state using plasma deltas and electric field
    subroutine applyUpdate(Model,Grid,State,oState,ooState,dt,dGasH,dGasM,E)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State,oState,ooState
        real(rp), intent(in) :: dt
        real(rp), intent(in) :: dGasH(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NVAR,BLK:Model%nSpc)
        real(rp), optional, intent(in) :: dGasM(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), optional, intent(inout) :: E(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM)

        integer :: i,j,k,nv,s
        real(rp), dimension(NVAR,0:Model%nSpc) :: U,oU,dPg
        real(rp), dimension(NDIM) :: B,oB,dPm
        logical :: BorisUpdate
        real(rp) :: Ca
        !DIR$ ASSUME_ALIGNED dGasH: ALIGN
        !DIR$ ASSUME_ALIGNED dGasM: ALIGN
        !DIR$ ASSUME_ALIGNED E: ALIGN
    !--------------------
    !Copy current->old
        call Tic("Copy2Old")
        !Shift states backwards
        if (Model%doAB3) then
            call CopyState(Model,Grid,oState,ooState)
        endif

        call CopyState(Model,Grid,State,oState)
        State%time = State%time + dt

        call Toc("Copy2Old")

    !--------------------
    !Apply hydro stresses
        !State%Gas = State%Gas + dt*dGasH

        call Tic("Reynolds")
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,U,dPg)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=Grid%is,Grid%ie
                    U   = State%Gas(i,j,k,:,:)
                    dPg = dGasH    (i,j,k,:,:)
                    call CellReynolds(Model,U,dPg,Model%dt)
                    State%Gas(i,j,k,:,:) = U
                enddo
            enddo
        enddo
        call Toc("Reynolds")

        !If we're only doing hydro then we're done here
        if (.not. Model%doMHD) then
            return
        endif

    !--------------------
    !Update magnetic fields
        !Only here if we're doing MHD

        call Tic("CT-Update")
        !Now apply E-field update to calculate updated MAG FLUXES
        !Average pole field if doing ring average
        if ( Model%doRing .and. (Model%Ring%doS .or. Model%Ring%doE) ) then
            call PoleE(Model,Grid,E)
        endif

        call E2Flux(Model,Grid,State%magFlux,E,dt)
        !Calculate Bxyz's
        call bFlux2Fld(Model,Grid,State%magFlux,State%Bxyz)
        call Toc("CT-Update")
        
    !--------------------
    !Apply Maxwell/Boris stresses
        !If multifluid then always do Boris (use HUGE if no value specified)
        if (Model%doMultiF) then
            BorisUpdate = .true.
            if (Model%doBoris) then
                Ca = Model%Ca
            else
                Ca = HUGE
            endif
        else
            !Not multifluid
            if (Model%doBoris) then
                BorisUpdate = .true.
                Ca = Model%Ca
            else
                BorisUpdate = .false.
                Ca = 0.0
            endif
        endif

        call Tic("Maxwell")
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,U,oU,B,oB,dPg,dPm) 
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=Grid%is,Grid%ie
                    U  =  State%Gas(i,j,k,:,:)
                    oU = oState%Gas(i,j,k,:,:)
                    B  =  State%Bxyz(i,j,k,:)
                    oB = oState%Bxyz(i,j,k,:)
                    if (Model%doBackground) then
                        B  =  B + Grid%B0(i,j,k,:)
                        oB = oB + Grid%B0(i,j,k,:) 
                    endif
                    dPg = dGasH(i,j,k,:,:)
                    dPm = dGasM(i,j,k,:)
                    if (BorisUpdate) then
                        call CellBoris(Model,U,oU,B,oB,dPg,dPm,Ca,dt)
                    else
                        call CellMaxwell(Model,U,dPm,dt)
                    endif
                    State%Gas(i,j,k,:,:) = U

                enddo
            enddo
        enddo

        call Toc("Maxwell")

        contains

            !Copies iState => oState
            subroutine CopyState(Model,Grid,iState,oState)
                type(Model_T), intent(in)    :: Model
                type(Grid_T ), intent(in)    :: Grid
                type(State_T), intent(in)    :: iState
                type(State_T), intent(inout) :: oState
                
                integer :: i,j,k
                
            !This is the code we want, but it's slow so we can't do it this way
            !This is why we can't have nice things

                ! oState%Gas = iState%Gas
                ! if (Model%doMHD) then
                !     !Note, not bothering w/ Efld since that's only ever at half state
                !     oState%Bxyz    = iState%Bxyz
                !     oState%magFlux = iState%magFlux
                !     if (Model%doResistive) then
                !         oState%Deta = iState%Deta
                !     endif
                ! endif
                ! oState%time = iState%time


            !Let's goddamn do this thing w/ interleaved memory copies
                !Open one big-ass // block

                !$OMP PARALLEL default(shared) &
                !$OMP private(i,j,k)

                !$OMP DO collapse(2)
                do k=Grid%ksg,Grid%keg
                    do j=Grid%jsg,Grid%jeg
                        do i=Grid%isg,Grid%ieg
                            oState%Gas(i,j,k,:,:) = iState%Gas(i,j,k,:,:)
                            oState%Bxyz(i,j,k,:)  = iState%Bxyz(i,j,k,:)
                            if (Model%doResistive .and. Model%doMHD) then
                                oState%Deta(i,j,k,:) = iState%Deta(i,j,k,:)
                            endif
                        enddo
                    enddo
                enddo
                !$OMP END DO NOWAIT !Rely on barrier at end of parallel region
                !Oh, you're already done? Well start the next thing

                !$OMP DO collapse(2)
                do k=Grid%ksg,Grid%keg+1
                    do j=Grid%jsg,Grid%jeg+1
                        do i=Grid%isg,Grid%ieg+1
                            oState%magFlux(i,j,k,:) = iState%magFlux(i,j,k,:)
                        enddo
                    enddo
                enddo
                !$OMP END PARALLEL !You're done when I say you're done ... which is now
                oState%time = iState%time

            end subroutine CopyState

    end subroutine applyUpdate

    !Apply reynolds update to cell over all species
    !In multifluid case do accumulation of species->bulk
    subroutine CellReynolds(Model,U,dGasH,dt)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U(NVAR,BLK:Model%nSpc)
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(inout) :: dGasH
        real(rp), intent(in) :: dt

        integer :: s,n
        logical , dimension(Model%nSpc) :: isGood
        real(rp), dimension(NDIM) :: oB

        if (Model%doMultiF) then
            isGood = ( U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
            !Update individual species then accumulate
            do s=1,Model%nSpc
                if (isGood(s)) then
                    U(:,s) = U(:,s) + dt*dGasH(:,s)
                else
                    !Always apply update (even to cells that aren't above floor)
                    !Otherwise they'll never get above floor
                    !TODO: Test only updating "bad" fluids if delta-rho > 0, i.e. don't take from vacuum
                    U(DEN,s) = U(DEN,s) + dt*dGasH(DEN,s)
                endif
            enddo
            oB = 0.0
            call MultiF2Bulk(Model,U,oB)
            !Reset bulk delta to sum of component species
            do n=1,NVAR
                dGasH(n,BLK) = sum(dGasH(n,1:Model%nSpc),mask=isGood)
            enddo
        else
            !Just do bulk flow and get outta here
            U(:,BLK) = U(:,BLK) + dt*dGasH(:,BLK)
        endif

    end subroutine CellReynolds

    !Apply Maxwell update to cell (only used in single species case)
    subroutine CellMaxwell(Model,U,dGasM,dt)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U(NVAR,BLK:Model%nSpc)
        real(rp), dimension(NDIM), intent(in) :: dGasM
        real(rp), intent(in) :: dt

        real(rp) :: D,Mx,My,Mz,TotE,KinE,P

        D = max( U(DEN,BLK),dFloor )
        Mx = U(MOMX,BLK)
        My = U(MOMY,BLK)
        Mz = U(MOMZ,BLK)
        TotE = U(ENERGY,BLK)

        KinE = 0.5*(Mx**2.0 + My**2.0 + Mz**2.0)/D
        P = max((Model%gamma-1)*( TotE - KinE ),pFloor)

        !Perform momentum update
        Mx = Mx + dt*dGasM(XDIR)
        My = My + dt*dGasM(YDIR)
        Mz = Mz + dt*dGasM(ZDIR)
        KinE = 0.5*(Mx**2.0 + My**2.0 + Mz**2.0)/D
    
        !Update variables
        U(DEN   ,BLK) = D !Incorporate floor
        U(MOMX  ,BLK) = Mx
        U(MOMY  ,BLK) = My
        U(MOMZ  ,BLK) = Mz
        U(ENERGY,BLK) = KinE + P/(Model%gamma-1)

    end subroutine CellMaxwell


    !Apply Boris update to cell over all species
    !Assuming already have updated hydro w/ dGasH
    !U = hydro-updated state, do update in place
    !Assuming B2/B0 are *TOTAL* fields (ie includes background)
    subroutine CellBoris(Model,U2,U0,B2,B0,dGasH,dGasM,Ca,dt)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U2(NVAR,BLK:Model%nSpc)
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(in) :: U0,dGasH
        real(rp), dimension(NDIM), intent(in) :: B2, B0, dGasM
        real(rp), intent(in) :: Ca,dt

        integer :: s
        real(rp) :: D0 ,D2 ,P2 ,alpha ,pCon ,dRho
        real(rp) :: D0s,D2s,P2s,alphaS,pConS,dRhoS,Rs,Rsbar
        real(rp), dimension(NDIM) :: B1,bhat1,bhat2,dPMag
        real(rp), dimension(NDIM) :: V0,V2,V0s,Vtmp,Vperp
        real(rp), dimension(NDIM) :: dPGas,dPGasS,dMom
        real(rp), dimension(NVAR) :: W2

    !Start w/ bulk update even in multifluid case
        !Get D, pressure from hydro update
        call CellC2P(Model,U2(:,BLK),W2)
        D2 = max(W2(DEN),dFloor)
        P2 = max(W2(PRESSURE),pFloor)

        !Get old D/V
        D0 = U0(DEN      ,BLK)
        D0 = max(D0,dFloor)
        V0 = U0(MOMX:MOMZ,BLK)/D0

        !B0 = n, B1 = n+1/2, B2 = B
        B1 = 0.5*(B2+B0) !Half-step field
        !Get terms for updating
        alpha = dot_product(B1,B1)/(D2*Ca*Ca) !Ratio of magnetic/plasma mass
        dRho  = dt*dGasH(DEN      ,BLK)
        dPGas = dt*dGasH(MOMX:MOMZ,BLK)
        dPMag = dt*dGasM(XDIR:ZDIR)

        pCon = 1.0/(1.0+alpha) !Dg/(Db+Dg), ratio of plasma mass and total (w/ magnetic) mass
        bhat1 = normVec(B1)
        bhat2 = normVec(B2) !Final field direction

        if (alpha > TINY) then
            !Note, only using dPMag in perp 
            dMom = Vec2Para(dPGas,bhat1) + pCon*Vec2Perp(dPGas+dPMag,bhat1) &
                 + alpha*pCon*dRho*Vec2Perp(V0,bhat1) 
        else 
            !Null field region, just do Maxwell update
            dMom = dPGas + dPMag
        endif !alpha

        !Finish bulk fluid update using dMom, reset primitives then convert to conserved
        V2 = (D0*V0 + dMom)/D2
        W2(DEN      ) = D2
        W2(VELX:VELZ) = V2
        W2(PRESSURE ) = P2
        call CellP2C(Model,W2,U2(:,BLK))

        !Now handle multifluid update
        if (Model%doMultiF) then
            !Do species independent things
            Vperp = Vec2Perp(V2,bhat2)

            !Loop over species
            do s=1,Model%nSpc
                !Get D, pressure from hydro update
                call CellC2P(Model,U2(:,s),W2)
                D2s = max(W2(DEN),dFloor)
                P2s = max(W2(PRESSURE),pFloor)
                
                !Get old D/V
                D0s = U0(DEN      ,s)
                D0s = max(D0s,dFloor)
                V0s = U0(MOMX:MOMZ,s)/D0s

                if (D2s >= Spcs(s)%dVac) then
                    alphaS = alpha*0.5*(D2s+D0s)/D2s
                    pConS = 1.0/(1.0+alphaS)
                    Rs = D2s/D2
                    Rsbar = 0.5*( D2s/D2 + D0s/D0 )
                    dRhoS = D2s-D0s
                    dPGasS = dt*dGasH(MOMX:MOMZ,s)

                    !Note, dPMag depends on pCon (not pConS)
                    dMom = pCon*Rs*dPMag                          &
                         + alphaS*pConS*dRhoS*Vec2Perp(V0s,bhat1) &
                         + pConS*Rsbar*Vec2Perp(dPGas,bhat1)      &
                         + Vec2Para(dPGasS,bhat1)
                         
                    Vtmp = (D0s*V0s + dMom)/D2s
                    W2(DEN      ) = D2s
                    W2(VELX:VELZ) = Vperp + Vec2Para(Vtmp,bhat2)
                    W2(PRESSURE ) = P2s                    
                else 
                    !Not a good species
                    Vtmp = V2
                    W2(DEN      ) = D2s
                    W2(VELX:VELZ) = 0.0
                    W2(PRESSURE ) = pFloor
                endif
                
                call CellP2C(Model,W2,U2(:,s))
            enddo

            !Done calculating species conserved quantities, recalculate bulk
            call MultiF2Bulk(Model,U2,bhat2)

        endif !Multifluid

    end subroutine CellBoris

    !Extrapolate forward in time by pdt using current (State) and 1-2 old states (oState/ooState)
    subroutine Predictor(Model,Grid,State,oState,ooState,pState,pdt)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(in) :: State, oState, ooState
        type(State_T), intent(inout) :: pState
        real(rp), intent(in) :: pdt

        integer :: i,j,k
        real(rp) :: odt,oodt
        logical :: isCC
        real(rp), dimension(NVAR,BLK:Model%nSpc) :: pQ,Q,oQ,ooQ
        real(rp), dimension(NDIM) :: Bxyz

        !Do timing stuff
        odt = State%time-oState%time !Sep. between states
        if (Model%doAB3) then
            oodt = oState%time-ooState%time
        else
            oodt = 0.0
        endif

        if(Model%doMHD) then    
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k)
            do k=Grid%ksg,Grid%keg+1
                do j=Grid%jsg,Grid%jeg+1
                    do i=Grid%isg,Grid%ieg+1
                        !Do interface fluxes
                        if(Model%doAB3) then
                            pState%magFlux(i,j,k,:) = ExtrapAB3(pdt,State%magFlux(i,j,k,:),oState%magFlux(i,j,k,:),ooState%magFlux(i,j,k,:),odt,oodt)
                            if (Model%doResistive) then
                                pState%Deta(i,j,k,:) = ExtrapAB3(pdt,State%Deta(i,j,k,:),oState%Deta(i,j,k,:),ooState%Deta(i,j,k,:),odt,oodt)
                            endif !Resistive
                        else
                            !Do interface fluxes
                            pState%magFlux(i,j,k,:) = ExtrapAB2(pdt,State%magFlux(i,j,k,:),oState%magFlux(i,j,k,:),odt)
                            if (Model%doResistive) then
                                pState%Deta(i,j,k,:) = ExtrapAB2(pdt,State%Deta(i,j,k,:),oState%Deta(i,j,k,:),odt)
                            endif !Resistive
                        endif
                    enddo !I loop
                enddo
            enddo !K loop
        endif

        !Now finish by getting Bxyz's
        !bflux2fld will call ring/singularity fixes
        call bFlux2Fld(Model,Grid,pState%magFlux,pState%Bxyz)

        !Loop over cell centers and do predictor on each cell of fluid/s
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,ooQ,oQ,Q,pQ,Bxyz)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg
                    !Grab NVAR,BLK:nSpc fluid info from this cell across time
                    Q  =  State%Gas(i,j,k,:,:)
                    oQ = oState%Gas(i,j,k,:,:)
                    if (Model%doAB3) then
                        ooQ = oState%Gas(i,j,k,:,:)
                    else
                        ooQ = 0.0
                    endif
                    if (Model%doMHD) then
                        Bxyz = pState%Bxyz(i,j,k,:)
                    else
                        Bxyz = 0.0
                    endif

                    call PredictCell(Model,Bxyz,pdt,odt,oodt,Q,oQ,ooQ,pQ)
                    pState%Gas(i,j,k,:,:) = pQ
                enddo !i
            enddo !j
        enddo !k

        !TODO: Remove hackpredictor nonsense
        if (associated(Model%HackPredictor)) then
            call Model%HackPredictor(Model,Grid,pState)
        endif

    end subroutine Predictor

    !Predict single cell of fluid/s forward in time by pdt into conserved state pQ
    !using current/older conserved states Q,oQ,ooQ separated by odt and oodt in time
    subroutine PredictCell(Model,Bxyz,pdt,odt,oodt,Q,oQ,ooQ,pQ)
        type(Model_T), intent(in) :: Model
        real(rp)     , intent(in) :: pdt,odt,oodt,Bxyz(NDIM)
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(in)  :: Q,oQ,ooQ
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(out) :: pQ

        logical, dimension(1:Model%nSpc) :: isGoodX,isGood1,isGood2
        integer :: s,s0,sE
        real(rp), dimension(NVAR) :: pW,W,oW,ooW

        !Q = conserved, W = primitive
        pQ = 0.0
        
        !If not MF, just do BLK and get the hell outta here
        if (.not. Model%doMultiF) then
            !Start by converting states to primitive
            call CellC2P(Model,oQ(:,BLK),oW)
            call CellC2P(Model, Q(:,BLK), W)

            if (Model%doAB3) then
               call CellC2P(Model,ooQ(:,BLK),ooW) 
               pW = ExtrapAB3(pdt,W,oW,ooW,odt,oodt)
            else
                !AB2, so just do it already
                pW = ExtrapAB2(pdt,W,oW,odt)
            endif

            !Apply global clamps
            pW(DEN)      = max(pW(DEN)     ,dFloor)
            pW(PRESSURE) = max(pW(PRESSURE),pFloor)

            !Return to conserved for final predicted state
            call CellP2C(Model,pW,pQ(:,BLK))
            return
        endif

    !If still here, we're doing multifluid

        !Start by testing states
        !X = current state, 1 = oState, 2 = ooState
        isGoodX = (  Q(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        isGood1 = ( oQ(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        if (Model%doAB3) then
            isGood2 = (ooQ(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        else
            isGood2 = .false. !Treat ooState as bad always if doing AB2
        endif

        s0 = 1
        sE = Model%nSpc

        !Logic:
        !Either of X,1 is bad => most recent
        !Both of X,1 are good =>
        !    ifAB3 => Test 2
        !    ifAB2 => Do AB2
        do s=s0,sE
            if ( (.not. isGoodX(s0)) .or. (.not. isGood1(s0)) ) then
                !at least 1 state is bad, just use most recent
                pQ(:,s) = Q(:,s)
            else
                !Both current/oState are good, convert these states
                call CellC2P(Model,oQ(:,s),oW)
                call CellC2P(Model, Q(:,s), W)

                if (isGood2(s)) then
                    !We're doing AB3 and ooState is good, so do ab3
                    call CellC2P(Model,ooQ(:,s),ooW) 
                    pW = ExtrapAB3(pdt,W,oW,ooW,odt,oodt)
                else
                    !Just do AB2
                    pW = ExtrapAB2(pdt,W,oW,odt)
                endif

                !Now treat state
                !Apply global clamps
                pW(DEN)      = max(pW(DEN)     ,dFloor)
                pW(PRESSURE) = max(pW(PRESSURE),pFloor)

                !Return to conserved for final predicted state
                call CellP2C(Model,pW,pQ(:,s))
            endif 
        enddo
        
        !Finally, accumulate MF to bulk
        call MultiF2Bulk(Model,pQ,Bxyz)

    end subroutine PredictCell

    !Updates plasma state using static gravitational potential
    !Use time-averaged density from n,n+1 and update momentum
    subroutine applyGrav(Model,Gr,State,oState,ooState,dt)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State, oState, ooState
        real(rp), intent(in) :: dt

        integer :: i,j,k,s,s0,sE
        
        real(rp) :: Dhf,Dblk,IntE
        real(rp), dimension(NDIM) :: g,B,bMxyz,cMxyz
        real(rp), dimension(NDIM,NDIM) :: Lam,Laminv
        real(rp), dimension(NVAR) :: pW,pCon
        real(rp), dimension(0:Model%nSpc) :: RhoMin


        RhoMin(BLK) = 0.0
        if (Model%doMultiF) then
            !Don't do bulk
            s0 = 1
            sE = Model%nSpc
            RhoMin(1:Model%nSpc) = Spcs(:)%dVac
        else
            !Only do bulk
            s0 = BLK
            sE = BLK
        endif

        !Add grav forces
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(s,i,j,k,Dhf,Dblk,bMxyz,cMxyz) &
        !$OMP private(g,B,Lam,Laminv,pW,pCon,IntE)
        do s=s0,sE
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%is,Gr%ie
                        if ( State%Gas(i,j,k,DEN,s) >= RhoMin(s) ) then
                        !Get necessary state info
                            pCon = State%Gas(i,j,k,1:NVAR,s) !Conserved
                            call CellC2P(Model,pCon,pW) !Primitive
                            Dblk = State%Gas(i,j,k,DEN,BLK) !Bulk density
                            IntE = pW(PRESSURE)/(Model%gamma-1)

                            !Get average density (species)
                            Dhf = 0.5*( State%Gas(i,j,k,DEN,s) + oState%Gas(i,j,k,DEN,s) )
                            g = Gr%gxyz(i,j,k,XDIR:ZDIR)
                            cMxyz = pCon(MOMX:MOMZ) !Classical momentum

                            if (doBorisGrav .and. Model%doBoris) then

                            !Map classical momentum to boris momentum
                                !NOTE: transform is I if Boris isn't on
                                
                                B = State%Bxyz(i,j,k,:)
                                call Mom2Rel(Model,Dblk,B,Lam)
                                bMxyz = matmul(Lam,cMxyz)
                            !Apply gravitational force to boris momentum
                                
                                bMxyz = bMxyz + Model%dt*Dhf*g
                            !Map updated Boris momentum back to classical
                                call Rel2Mom(Model,Dblk,B,Laminv)
                                cMxyz = matmul(Laminv,bMxyz) !Classical momentum
                            !Set new primitive variables, convert and store
                                !Only V has changed
                                pW(VELX:VELZ) = cMxyz/max(pW(DEN),dFloor)
                                call CellP2C(Model,pW,pCon)
                                State%Gas(i,j,k,:,s) = pCon
                            else
                                !Regular update
                                cMxyz = cMxyz + Model%dt*Dhf*g
                                !Reset conserved state
                                pCon(MOMX:MOMZ) = cMxyz
                                pCon(ENERGY   ) = IntE + 0.5*dot_product(cMxyz,cMxyz)/pCon(DEN)
                                State%Gas(i,j,k,:,s) = pCon
                            endif !Classical vs. Boris momentum update

                        endif !RhoMin
                    enddo !i
                enddo !j
            enddo !k
        enddo !s

        if (Model%doMultiF) then
            call State2Bulk(Model,Gr,State)
        endif

    end subroutine applyGrav

!=====
!Adams Bashforth helper routines
    !Extrapolate forward in time by dt given Q(t) and oQ(t-odt)
    elemental function ExtrapAB2(dt,Q,oQ,odt) result(nQ)
        real(rp), intent(in) :: dt,Q,oQ,odt
        real(rp) :: nQ

        real(rp) :: dQ

        dQ = (Q-oQ)/odt
        nQ = Q + dt*dQ
    end function ExtrapAB2

    !Extrapolate forward in time by dt given:
    !Q(t),oQ(t-odt),ooQ(t-odt-oodt)
    elemental function ExtrapAB3(dt,Q,oQ,ooQ,odt,oodt) result(nQ)
        real(rp), intent(in) :: dt,Q,oQ,ooQ,odt,oodt
        real(rp) :: nQ

        real(rp) :: dQ,odQ,tScl
        
        dQ = (Q-oQ)/odt
        odQ = (oQ-ooQ)/oodt
        tScl = (dt + odt)/(odt + oodt)
        nQ = Q + dt*dQ + dt*tScl*(dQ - odQ)
    end function ExtrapAB3
    
end module mhdgroup
