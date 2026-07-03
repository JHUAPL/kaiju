module raijuColdStartHelper

    use raijudefs
    use raijutypes
    use imaghelper
    use earthhelper
    use arrayutil

    use raijuetautils
    use raijuloss_CX

    implicit none

    contains

!------
! ColdStarter type stuff
!------

    subroutine initRaijuColdStarter(Model, iXML, coldStarter, tEndO)
        type(raijuModel_T), intent(in) :: Model
        type(XML_Input_T), intent(in) :: iXML
        type(raijuColdStarter_T), intent(inout) :: coldStarter
        real(rp), intent(in), optional :: tEndO

        
        call iXML%Set_Val(coldStarter%doCX,'coldStarter/doCX',coldStarter%doCX)
        call iXML%Set_Val(coldStarter%doUpdate,'coldStarter/doUpdate',coldStarter%doUpdate)
        call iXML%Set_Val(coldStarter%evalCadence,'coldStarter/evalCadence',coldStarter%evalCadence)
        if (present(tEndO)) then
            call iXML%Set_Val(coldStarter%tEnd,'coldStarter/tEnd',tEndO)
        else
            call iXML%Set_Val(coldStarter%tEnd,'coldStarter/tEnd',coldStarter%evalCadence-TINY)  ! Don't do any updates as default
        endif

        coldStarter%doneFirstCS = .false.

    end subroutine initRaijuColdStarter




!------
! Worker routines
!------

    subroutine raijuGeoColdStart(Model, Grid, State, t0, dstModel, doAccumulateO)
        !! Cold start RAIJU assuming we are at Earth sometime around 21st century
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        real(rp), intent(in) :: t0
            !! Target time to pull SW values from
        real(rp), intent(in) :: dstModel
            !! Current dst of global model
        logical, optional, intent(in) :: doAccumulateO
            !! If true, keep State%eta and add coldstart to it.
            !! If false, replace State%eta with coldStart info
            !! Default: false

        logical :: doInitRC, doAccumulate
        integer :: i,j,k
        integer :: s, sIdx_p, sIdx_e
        real(rp) :: dstReal, dstTarget
        real(rp) :: dps_current, dps_preCX, dps_postCX, dps_rescale, dps_ele
        real(rp) :: etaScale
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, Grid%Nk) :: etaCS

        if (present(doAccumulateO)) then
            doAccumulate = doAccumulateO
        else
            doAccumulate = .false.
        endif

        call fillArray(etaCS, 0.0_rp)

        associate(cs=>State%coldStarter)

        where (State%active .ne. RAIJUINACTIVE)
            isGood = .true.
        elsewhere
            isGood = .false.
        endwhere

        sIdx_p = spcIdx(Grid, F_HOTP)
        sIdx_e = spcIdx(Grid, F_HOTE)

        ! Update Dst target
        dstReal = GetSWVal('symh', Model%tsF, t0)
        if (cs%doneFirstCS) then
            write(*,*)"Already coldstarted once, you shouldn't be here "
            return
        else
            ! On first try, we assume there is no existing ring current, and its our job to make up the entire difference
            dstTarget = dstReal - dstModel
        endif
        
        cs%lastEval = t0
        cs%lastTarget = dstTarget
        cs%doneFirstCS = .true.  ! Whether we do anything or not, we were at least called once
        
        ! Now decide if we need to add a starter ring current
        if (dstTarget >= 0) then  ! We've got nothing to contribute
            write(*,*)"RAIJU coldstart not adding starter ring current"
            return
        endif

        ! Init hot protons
        call raiColdStart_initHOTP(Model, Grid, State, t0, dstTarget, etaCS)
        !call raiColdStart_initHOTP_RCOnly(Model, Grid, State, t0, dstTarget, etaCS)
        dps_preCX  = spcEta2DPS(Model, Grid, State%bvol_cc, etaCS, Grid%spc(sIdx_p), isGood)
        ! Hit it with some charge exchange
        if (cs%doCX) then
            call raiColdStart_applyCX(Model, Grid, State, Grid%spc(sIdx_p), etaCS)
        endif
        dps_postCX = spcEta2DPS(Model, Grid, State%bvol_cc, etaCS, Grid%spc(sIdx_p), isGood)
        ! Use HOTP moments to set electrons
        call raiColdStart_initHOTE(Model, Grid, State, etaCS)
        dps_ele = spcEta2DPS(Model, Grid, State%bvol_cc, etaCS, Grid%spc(sIdx_e), isGood)
        dps_current = dps_postCX  ! Note: if using fudge we're gonna lose electrons immediately, don't include them in current dst for now

        if (dstTarget < 0) then
            etaScale = abs(dstTarget / dps_current)    
            etaCS(:,:,Grid%spc(sIdx_p)%kStart:Grid%spc(sIdx_p)%kEnd) = etaScale*etaCS(:,:,Grid%spc(sIdx_p)%kStart:Grid%spc(sIdx_p)%kEnd)
            dps_rescale = spcEta2DPS(Model, Grid, State%bvol_cc, etaCS, Grid%spc(sIdx_p), isGood)
        else
            dps_rescale = dps_current
        endif

        write(*,*)          "RAIJU Cold starting..."
        write(*,'(a,f7.2)') "  Real Dst             : ",dstReal
        write(*,'(a,f7.2)') "  Model Dst            : ",dstModel
        write(*,'(a,f7.2)') "  Target DPS-Dst       : ",dstTarget
        write(*,'(a,f7.2)') "  Hot proton pre-loss  : ",dps_preCX
        write(*,'(a,f7.2)') "            post-loss  : ",dps_postCX
        write(*,'(a,f7.2)') "         post-rescale  : ",dps_rescale
        write(*,'(a,f7.2)') "  Hot electron DPS-Dst : ",dps_ele

        end associate

        ! finally, put it into raiju state
        if(doAccumulate) then
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,k)
            do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                do i=Grid%shGrid%isg,Grid%shGrid%ieg
                    do k=1,Grid%Nk
                        if(etaCS(i,j,k) > State%eta(i,j,k)) then
                            State%eta(i,j,k) = etaCS(i,j,k)
                        endif
                    enddo
                enddo
            enddo
        else
            State%eta = etaCS
        endif

    end subroutine raijuGeoColdStart


    subroutine raiColdStart_initHOTP_RCOnly(Model, Grid, State, t0, dstTarget, etaCS)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), intent(in) :: t0
            !! Target time to pull SW values from
        real(rp), intent(in) :: dstTarget
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, Grid%Nk), intent(inout) :: etaCS

        real(rp) :: dstTarget_p
        logical :: isInTM03
        integer :: i,j,sIdx
        integer, dimension(2) :: ij_TM
        real(rp) :: vSW, dSW, dPS_emp, pPS_emp, ktPS_emp
        real(rp) :: x0_TM, y0_TM, T0_TM, Bvol0_TM, P0_ps, N0_ps
        real(rp) :: L, vm, P_rc, D_rc, kt_rc

        sIdx = spcIdx(Grid, F_HOTP)
        associate(sh=>Grid%shGrid, spc=>Grid%spc(sIdx))

        ! Set everything to zero to start
        etaCS(:,:,spc%kStart:spc%kEnd) = 0.0_rp

        ! Scale target Dst down to account for electrons contributing stuff later
        dstTarget_p = dstTarget / (1.0 + 1.0/Model%tiote)
        call SetQTRC(dstTarget_p,doVerbO=.false.) ! This sets a global QTRC_P0 inside earthhelper.F90


        ! Get reference TM value at -10 Re
        x0_TM = -10.0-TINY
        y0_TM = 0.0
        ! Empirical temperature
        call EvalTM03([x0_TM,y0_TM,0.0_rp],N0_ps,P0_ps,isInTM03)
        T0_TM = DP2kT(N0_ps, P0_ps)
        ! Model FTV
        ij_TM = minloc( sqrt( (State%xyzMincc(:,:,XDIR)-x0_TM)**2 + (State%xyzMincc(:,:,YDIR)**2) ) )
        Bvol0_TM = State%bvol_cc(ij_TM(IDIR), ij_TM(JDIR))

        ! Now set our initial density and pressure profile
        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ie  ! Note: Not setting low lat ghosts, we want them to be zero
                
                if (State%active(i,j) .eq. RAIJUINACTIVE) cycle

                L = norm2(State%xyzMincc(i,j,XDIR:YDIR))
                vm = State%bvol_cc(i,j)**(-2./3.)

                kt_rc = T0_TM*(Bvol0_TM/State%bvol_cc(i,j))**(2./3.)
                kt_rc = min(kt_rc, 4.0*T0_TM)  ! Limit cap. Not a big fan, but without cap we get stuff that's too energetic and won't go away (until FLC maybe)


                P_rc = P_QTRC(L)  ! From earthhelper.F90
                D_rc = PkT2Den(P_rc, kt_rc)
                

                ! Finally map it to HOTP etas
                call DkT2SpcEta(Model, spc, etaCS(i,j,spc%kStart:spc%kEnd), D_rc, kt_rc, vm)

            enddo
        enddo

        end associate

    end subroutine raiColdStart_initHOTP_RCOnly


    subroutine raiColdStart_initHOTP(Model, Grid, State, t0, dstTarget, etaCS)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), intent(in) :: t0
            !! Target time to pull SW values from
        real(rp), intent(in) :: dstTarget
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, Grid%Nk), intent(inout) :: etaCS

        real(rp) :: dstTarget_p
        logical :: isInTM03
        integer :: i,j,sIdx
        integer, dimension(2) :: ij_TM
        real(rp) :: vSW, dSW, dPS_emp, pPS_emp, ktPS_emp
        real(rp) :: x0_TM, y0_TM, T0_TM, Bvol0_TM, P0_ps, N0_ps
        real(rp) :: L, vm, P_rc, D_rc, kt_rc, P_final, D_final

        P_final = 0.0
        D_final = 0.0
        sIdx = spcIdx(Grid, F_HOTP)

        associate(sh=>Grid%shGrid, spc=>Grid%spc(sIdx))

        ! Initialize TM03
        call InitTM03(Model%tsF,t0)

        ! Scale target Dst down to account for electrons contributing stuff later
        if (dstTarget < 0) then
            dstTarget_p = dstTarget / (1.0 + 1.0/Model%tiote)
            call SetQTRC(dstTarget_p,doVerbO=.false.) ! This sets a global QTRC_P0 inside earthhelper.F90
        endif

        ! Get Borovsky statistical values
        vSW = abs(GetSWVal("Vx",Model%tsF,t0))
        dSW =     GetSWVal("D" ,Model%tsF,t0)
        dPS_emp  = 0.292*(dSW**0.49)
        ktPS_emp = -3.65 + 0.0190*vSW*1.0e-3 !m/s=>km/s
        ktPS_emp = max(ktPS_emp,TINY)
        pPS_emp = DkT2P(dPS_emp, ktPS_emp)

        ! Get reference TM value at -10 Re
        x0_TM = -10.0-TINY
        y0_TM = 0.0
        ! Empirical temperature
        call EvalTM03([x0_TM,y0_TM,0.0_rp],N0_ps,P0_ps,isInTM03)
        T0_TM = DP2kT(N0_ps, P0_ps)
        ! Model FTV
        ij_TM = minloc( sqrt( (State%xyzMincc(:,:,XDIR)-x0_TM)**2 + (State%xyzMincc(:,:,YDIR)**2) ) )
        Bvol0_TM = State%bvol_cc(ij_TM(IDIR), ij_TM(JDIR))
        if (.not. isInTM03) then
            write(*,*) "This should not happen w/ TM03, you should figure this out ..."
            stop
        endif

        ! Set everything to zero to start
        etaCS(:,:,spc%kStart:spc%kEnd) = 0.0_rp

        ! Now set our initial density and pressure profile
        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ie  ! Note: Not setting low lat ghosts, we want them to be zero
                
                if (State%active(i,j) .eq. RAIJUINACTIVE) cycle

                L = norm2(State%xyzMincc(i,j,XDIR:YDIR))
                vm = State%bvol_cc(i,j)**(-2./3.)

                kt_rc = T0_TM*(Bvol0_TM/State%bvol_cc(i,j))**(2./3.)
                kt_rc = min(kt_rc, 4.0*T0_TM)  ! Limit cap. Not a big fan, but without cap we get stuff that's too energetic and won't go away (until FLC maybe)

                call EvalTM03_SM(State%xyzMincc(i,j,:),N0_ps,P0_ps,isInTM03)

                if (dstTarget < 0) then
                    P_rc = P_QTRC(L)  ! From earthhelper.F90
                else
                    P_rc = 0.0
                endif

                if (.not. isInTM03) then
                    N0_ps = dPS_emp
                    P0_ps = pPS_emp
                endif

                if (P0_ps > P_rc) then
                    P_final = P0_ps
                    D_final = PkT2Den(P0_ps, kt_rc)
                else
                    P_final = P_rc
                    D_final = PkT2Den(P_rc, kt_rc)
                endif

                ! Finally map it to HOTP etas
                call DkT2SpcEta(Model, spc, etaCS(i,j,spc%kStart:spc%kEnd), D_final, kt_rc, vm)

            enddo
        enddo

        end associate

    end subroutine raiColdStart_initHOTP


    subroutine raiColdStart_applyCX(Model, Grid, State, spc, etaCS)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        type(raijuSpecies_T), intent(in) :: spc
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, Grid%Nk), intent(inout) :: etaCS

        integer :: i,j,k
        type(raiLoss_CX_T) :: lossCX
        type(XML_Input_T) :: nullXML  ! Needed for raiLoss inits, but CX doesn't actually need it, so make a dummy one
        real(rp) :: tCX = 12*3600  ! [s] Amount of time to apply CX for
        real(rp) :: tau

        call lossCX%doInit(Model, Grid, nullXML)
        
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,tau)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                do k = spc%kStart,spc%kEnd
                    tau = lossCX%calcTau(Model, Grid, State, i, j, k)
                    etaCS(i,j,k) = etaCS(i,j,k)*exp(-tCX/tau)
                enddo
            enddo
        enddo

    end subroutine raiColdStart_applyCX


    subroutine raiColdStart_initHOTE(Model, Grid, State, etaCS)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, Grid%Nk), intent(inout) :: etaCS

        integer :: sIdx_e, sIdx_p
        integer :: i,j
        real(rp) :: press_p, den_p
        real(rp) :: kt_p, kt_e, den, vm

        sIdx_p = spcIdx(Grid, F_HOTP)
        sIdx_e = spcIdx(Grid, F_HOTE)


        associate(spc_p=>Grid%spc(sIdx_p), spc_e=>Grid%spc(sIdx_e))
        ! Set everything to zero to start
        etaCS(:,:,spc_e%kStart:spc_e%kEnd) = 0.0_rp

        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,vm,den_p, press_p,kt_p,kt_e)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ie  ! Note: Not setting low lat ghosts, we want them to be zero
                if (State%active(i,j) .eq. RAIJUINACTIVE) cycle

                vm = State%bvol_cc(i,j)**(-2./3.)
                !den = State%Den(sIdx_p)%data(i,j)
                !kt_p = DP2kT(den, State%Press(sIdx_p)%data(i,j))

                den_p   = SpcEta2Den  (spc_p, etaCS(i,j,spc_p%kStart:spc_p%kEnd), State%bvol_cc(i,j))
                press_p = SpcEta2Press(spc_p, etaCS(i,j,spc_p%kStart:spc_p%kEnd), State%bvol_cc(i,j))
                kt_p = DP2kT(den_p, press_p)

                kt_e = kt_p / Model%tiote
                call DkT2SpcEta(Model, Grid%spc(sIdx_e), &
                                etaCS(i,j,Grid%spc(sIdx_e)%kStart:Grid%spc(sIdx_e)%kEnd), &
                                den_p, kt_e, vm)
            enddo
        enddo

        end associate

    end subroutine raiColdStart_initHOTE


    function GetSWVal(vID,fID,t) result(qSW)
        character(len=*), intent(in) :: vID,fID
        real(rp), intent(in) :: t
        real(rp) :: qSW

        type(TimeSeries_T) :: tsQ
        tsQ%wID = trim(fID)
        call tsQ%initTS(trim(vID),doLoudO=.false.)
        qSW = tsQ%evalAt(t)
    end function GetSWVal

end module raijuColdStartHelper