!remixCoupler_T procedures (mirrors raijuCplTypesSub.F90 structure).
!Data path: vApp%remixApp -> App%mixApp.

submodule (volttypes) remixCplTypesSub

    use mixtypes
    use mixdefs
    use mhd2mix_interface
    use mix2mhd_interface
    use imag2mix_interface
    use mix2shell_interface,  only: mapMixToShell
    use voltappHelper
    use kai2geo
    use mixmain,   only: updateF107  ! not host-associated via volttypes
    use shellutils,   only: wrapJ_SGV                      ! periodic phi ghost fill for the DK->mix map
    ! DK is now computed on the shell (shellDK in voltapp); the coupler only MAPS the
    ! full-sphere DK conductance to REMIX via mapDKToMix -> voltVarToMix (no dkmain/dkinit here).
    use planethelper, only: DP2Cs                          ! (D,P) -> sound speed for tubeShell forward feed
    use kdefs,        only: Mp_cgs, TINY, PI                ! #/cc -> mass density conversion; PI for conjugate geometry
    use gdefs,        only: BLK                             ! bulk fluid index for tubeShell bcN/bcP

    implicit none

    contains

    !InitModel — REMIX-internal init only; gApp-dependent map builds stay in voltapp.
    module subroutine remixCplInitModel(App, xml)
        class(remixCoupler_T), intent(inout) :: App
        type(XML_Input_T),     intent(inout) :: xml

        call App%mixApp%InitModel(xml)

    end subroutine remixCplInitModel


    !InitCouplers — REMIX-specific init behind the ionoCoupler contract: restart/IO,
    !F10.7, ionosphere radii. Called polymorphically (no select type).
    module subroutine remixCpl_InitCouplers(App, vApp, xml)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp
        type(XML_Input_T),     intent(inout) :: xml

        type(TimeSeries_T)    :: f107
        character(len=strLen) :: RunID, resID
        logical               :: isRestart
        real(rp)              :: maxF107
        integer               :: nRes

        isRestart = vApp%gApp%Model%isRestart
        RunID     = trim(vApp%gApp%Model%RunID)

        if (isRestart) then
            call App%ReadRestart(resId=RunID, nRes=vApp%IO%nRes)
        endif

        if (vApp%writeFiles) then
            call App%InitIO(xml)
        endif

        App%mixApp%ion%rad_iono_m   = vApp%planet%ri_m
        App%mixApp%ion%rad_planet_m = vApp%planet%rp_m

        ! Ensure remix and voltron restart numbers match
        if (isRestart .and. vApp%IO%nRes /= App%mixApp%ion(1)%P%nRes) then
            write(*,*) "Voltron and Remix disagree on restart number, you should sort that out."
            write(*,*) "Error code: A house divided cannot stand"
            write(*,*) "   Voltron nRes = ", vApp%IO%nRes
            write(*,*) "   Remix   nRes = ", App%mixApp%ion(1)%P%nRes
            stop
        endif

        if (App%mixApp%ion(NORTH)%P%doSWF107 .neqv. App%mixApp%ion(SOUTH)%P%doSWF107) then
            write(*,*) 'Something is wrong. doSWf107 is set differently for the two hemispheres.'
            write(*,*) 'Exiting ...'
            stop
        endif

        ! read f107 from the SW file and overwrite what's been read from .xml above (in init_mix)
        ! note, only checking for NORTH, because both hemispheres read the same xml file
        if (App%mixApp%ion(NORTH)%P%doSWF107) then
            f107%wID = vApp%tilt%wID
            call f107%initTS("f10.7", doLoudO=.false.)
            maxF107 = f107%getMax()
            call updateF107(App%mixApp%ion, maxF107)
            vApp%dkApp%params%f107 = maxF107
            if (vApp%isLoud) write(*,*) 'Using F10.7 = ', maxF107
        endif

        if (vApp%isLoud) write(*,*) 'Using MJD0  = ', vApp%gApp%Model%MJD0

        ! REMIX -> Gamera maps (needs gApp)
        call init_mix2MhdCoupler(vApp%gApp, App%mixApp)
        ! initialize additional coupled gamera data
        call vApp%gApp%InitMhdCoupler(vApp)
        if (isRestart) then
            call xml%Set_Val(resID, "/Kaiju/gamera/restart/resID", "msphere")
            call xml%Set_Val(nRes,  "/Kaiju/gamera/restart/nRes", -1)
            call vApp%gApp%ReadRestart(resID, nRes)
        endif

        ! App%mhd_Rin (equatorward-mask radius) is set in initializeFromGamera after this returns.

    end subroutine remixCpl_InitCouplers


    !toIono: pack Gamera + IMAG into REMIX
    module subroutine remixCpl_toIono(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        real(rp) :: curTilt

        ! Recalculate MJD for current time
        call MJDRecalc(vApp%MJD)

        ! IMAG precipitation → REMIX (only after spinup, only when deep)
        ! Legacy imag2mix path -- superseded by mapDKToMix + mapNonDKToMix
        !if (vApp%doDeep .and. (vApp%time > 0)) then
        !    call CoupleIMagToMix(vApp, App%mixApp)
        !endif

        ! voltron ShellGrid -> REMIX hemisphere grids
        call App%packInputs(vApp)

        ! Set tilt and time for conductance models
        call vApp%tilt%getValue(vApp%time, curTilt)
        App%mixApp%IO%tilt = curTilt
        App%mixApp%IO%time = vApp%time

        ! DK conductance/precip (computed on the shell by shellDK) -> both REMIX hemispheres
        call mapDKToMix(vApp%shGrid, vApp%dkApp%vars, App%mixApp)
        ! Non-DK IM fields -> REMIX.
        if (vApp%doDeep .and. (vApp%time > 0)) call mapNonDKToMix(vApp, App%mixApp)

        ! DK is always applied first. If enabled, imported GCM conductance then overwrites
        ! REMIX SIGMAP/SIGMAH in state space.
        if (vApp%doGCM) call mapGCMToMix(vApp%gcm, App%mixApp)

        ! REMIX will live exclusively on SM now.
        !! Geographic grid for the REMIX potential solve (non-GCM path).
        !if (.not. vApp%doGCM) call mix_transform_grid(App%mixApp%ion)

    end subroutine remixCpl_toIono


    !fromIono — map solved REMIX potential back to Voltron + Gamera.
    module subroutine remixCpl_fromIono(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        ! REMIX hemisphere grids -> voltron ShellGrid.
        ! (shell2gam, shell -> Gamera inner boundary, runs in the driver after this.)
        call mapMixToShell(App%mixApp, vApp)

        ! Solved potential onto the shell corners for DK's doET enhancement
        call mixVarToVoltron(App%mixApp, POT, vApp%shGrid, vApp%pot_iono)

        ! Store MJD/time now so I/O calls need no special setup
        App%mixApp%IO%mjd  = vApp%MJD
        App%mixApp%IO%time = vApp%time

    end subroutine remixCpl_fromIono


    !AdvanceModel — steady-state ionospheric potential solve.
    module subroutine remixCplAdvanceModel(App, dt)
        class(remixCoupler_T), intent(inout) :: App
        real(rp),              intent(in)    :: dt

        call App%mixApp%AdvanceModel(dt)

    end subroutine remixCplAdvanceModel


    !Remaining BaseApp_T procedures — delegate to mixApp.
    module subroutine remixCplInitIO(App, xml)
        class(remixCoupler_T), intent(inout) :: App
        type(XML_Input_T),     intent(inout) :: xml

        call App%mixApp%InitIO(xml)

    end subroutine remixCplInitIO

    module subroutine remixCplWriteRestart(App, nRes)
        class(remixCoupler_T), intent(inout) :: App
        integer,               intent(in)    :: nRes

        ! IO%mjd and IO%time were set in fromIono — no extra setup needed
        call App%mixApp%WriteRestart(nRes)

    end subroutine remixCplWriteRestart

    module subroutine remixCplReadRestart(App, resId, nRes)
        class(remixCoupler_T), intent(inout) :: App
        character(len=*),      intent(in)    :: resId
        integer,               intent(in)    :: nRes

        call App%mixApp%ReadRestart(resId, nRes)

    end subroutine remixCplReadRestart

    module subroutine remixCplWriteConsoleOutput(App)
        class(remixCoupler_T), intent(inout) :: App

        call App%mixApp%WriteConsoleOutput()

    end subroutine remixCplWriteConsoleOutput

    module subroutine remixCplWriteFileOutput(App, nStep)
        class(remixCoupler_T), intent(inout) :: App
        integer,               intent(in)    :: nStep

        call App%mixApp%WriteFileOutput(nStep)

    end subroutine remixCplWriteFileOutput

    module subroutine remixCplWriteSlimFileOutput(App, nStep)
        class(remixCoupler_T), intent(inout) :: App
        integer,               intent(in)    :: nStep

        call App%mixApp%WriteSlimFileOutput(nStep)

    end subroutine remixCplWriteSlimFileOutput

    module subroutine remixCplCleanup(App)
        class(remixCoupler_T), intent(inout) :: App

        call App%mixApp%Cleanup()

    end subroutine remixCplCleanup


    !packInputs — voltron ShellGrid -> REMIX. Forward feed from tubeShell (J11i/bcN/bcP).
    !tubeShell is empty during spinup (const_sigma then, so zeros are benign) and populated
    !post-spinup. Equatorward floor/mask is self-contained (App%mhd_Rin + per-hemi floors).
    module subroutine remix_packInputs(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        type(ShellGridVar_T)                 :: facV, denV, csV
        logical, dimension(:,:), allocatable :: mask
        real(rp) :: floorD(2), floorCs(2)
        logical  :: doDeepFeed
        integer  :: i, j, h, Nth, is_, ie_, js_, je_

        doDeepFeed = (vApp%doDeep .and. vApp%time >= 0)

        ! Build forward-feed ShellGridVars from rvs tubeShell (empty=0 during spinup).
        call initShellVar(vApp%shGrid, SHGR_CORNER, facV)
        call initShellVar(vApp%shGrid, SHGR_CORNER, denV)
        call initShellVar(vApp%shGrid, SHGR_CORNER, csV)

        facV%data = vApp%State%tubeShell%J11i%data                  ! uA/m^2 (J11i is FAC already)
        denV%data = vApp%State%tubeShell%bcN(BLK)%data * Mp_cgs      ! #/cc -> mass density g/cm^3
        do j = lbound(csV%data,2), ubound(csV%data,2)
            do i = lbound(csV%data,1), ubound(csV%data,1)
                csV%data(i,j) = DP2Cs(vApp%State%tubeShell%bcN(BLK)%data(i,j), &
                                      vApp%State%tubeShell%bcP(BLK)%data(i,j)) * 1.0e5_rp  ! km/s -> cm/s
            enddo
        enddo
        ! bcN/bcP = 0 marks no-sample (dipole-region) tubes, not vacuum: mask them so TSC
        ! renormalizes onto valid corners. FAC 0 is physical -> no mask.
        facV%mask = .true.
        denV%mask = (denV%data > 0._rp)
        csV%mask  = (vApp%State%tubeShell%bcN(BLK)%data > 0._rp) .and. &
                    (vApp%State%tubeShell%bcP(BLK)%data > 0._rp)

        ! voltron shell -> REMIX's two hemisphere grids (TSC + supersample).
        call mapShellVarToMix(vApp%shGrid, facV, App%mixApp%ion, FAC)
        call mapShellVarToMix(vApp%shGrid, denV, App%mixApp%ion, DENSITY)
        call mapShellVarToMix(vApp%shGrid, csV,  App%mixApp%ion, SOUND_SPEED)

        ! Per-hemisphere floor over valid (non-dipole) corners (deep feed only)
        is_ = vApp%shGrid%is; ie_ = vApp%shGrid%ie
        js_ = vApp%shGrid%js; je_ = vApp%shGrid%je
        floorD = 0.0_rp; floorCs = 0.0_rp
        if (doDeepFeed) then
            do h = NORTH, SOUTH
                floorD(h)  = 1.0e30_rp
                floorCs(h) = 1.0e30_rp
                do j = js_, je_+1
                    do i = is_, ie_+1
                        if (sqrt(App%mhd_Rin)*sin(vApp%shGrid%th(i)) > 1._rp) cycle  ! dipole region
                        if ((h == NORTH) .neqv. (cos(vApp%shGrid%th(i)) > 0._rp)) cycle ! wrong hemi
                        if (denV%data(i,j) <= 0._rp) cycle  ! skip empty (open-tube / polar-cap) corners
                        floorD(h)  = min(floorD(h),  denV%data(i,j))
                        floorCs(h) = min(floorCs(h), csV%data(i,j))
                    enddo
                enddo
                if (floorD(h)  >= 1.0e30_rp) floorD(h)  = 0.0_rp
                if (floorCs(h) >= 1.0e30_rp) floorCs(h) = 0.0_rp
            enddo
        endif

        ! Equatorward of the MHD inner boundary: FAC->0, D/Cs->per-hemi floor.
        do h = NORTH, SOUTH
            Nth = App%mixApp%ion(h)%G%Nt
            allocate(mask(App%mixApp%ion(h)%G%Np, App%mixApp%ion(h)%G%Nt))
            mask = (sqrt(App%mhd_Rin) * sin(App%mixApp%ion(h)%G%t) > 1._rp - TINY)
            where (mask) App%mixApp%ion(h)%St%Vars(:,:,FAC)         = 0._rp
            where (mask) App%mixApp%ion(h)%St%Vars(:,:,DENSITY)     = floorD(h)
            where (mask) App%mixApp%ion(h)%St%Vars(:,:,SOUND_SPEED) = floorCs(h)
            App%mixApp%ion(h)%St%Vars(:, Nth-1:Nth, FAC) = 0._rp
            deallocate(mask)
        enddo

    end subroutine remix_packInputs


    ! DK conductance/precip (computed on the shell by shellDK) -> both REMIX hemispheres.
    ! DK outputs are cell-centered (Np,Nt) (phi,theta) arrays; pack into a SHGR_CC
    ! ShellGridVar (theta,phi) and TSC-supersample onto REMIX.
    subroutine mapDKToMix(voltGrid, dkVars, mixApp)
        type(ShellGrid_T),   intent(in)    :: voltGrid
        type(dkVariables_T), intent(in)    :: dkVars
        class(mixApp_T),     intent(inout) :: mixApp
        type(ShellGridVar_T) :: sgv
        integer :: is_, ie_, js_, je_

        is_ = voltGrid%is; ie_ = voltGrid%ie
        js_ = voltGrid%js; je_ = voltGrid%je

        ! Fields needed for the REMIX potential solve
        call mapOneDK(dkVars%SigmaP%data,   SIGMAP)
        call mapOneDK(dkVars%SigmaH%data,   SIGMAH)
        call mapOneDK(dkVars%Avg_Eng%data,  AVG_ENG)
        call mapOneDK(dkVars%Num_Flux%data, NUM_FLUX)

        ! DK computed / IM-fed continuous fields (TSC supersample)
        call mapOneDK(dkVars%DeltaE%data,   DELTAE)    ! mono potential drop
        call mapOneDK(dkVars%Z_Eavg%data,   Z_EAVG)    ! MHD average energy
        call mapOneDK(dkVars%Z_Nflux%data,  Z_NFLUX)   ! MHD number flux
        call mapOneDK(dkVars%Im_Beta%data,  IM_BETA)   ! IM beta
        call mapOneDK(dkVars%Im_Eden%data,  IM_EDEN)   ! IM electron density
        call mapOneDK(dkVars%Im_Epre%data,  IM_EPRE)   ! IM electron pressure
        call mapOneDK(dkVars%Im_Enflx%data, IM_ENFLX)  ! IM number flux e-
        call mapOneDK(dkVars%Im_Eflux%data, IM_EFLUX)  ! IM energy flux e-
        call mapOneDK(dkVars%Im_Npsp%data,  IM_NPSP)   ! plasmasphere density
        call mapOneDK(dkVars%Im_Eavg%data,  IM_EAVG)   ! IM avg energy (derived in shellDK)
        call mapOneDK(dkVars%Im_Gtype%data, IM_GTYPE)  ! IM grid-type weight [0,1]

        ! CATEGORICAL field -- AUR_TYPE
        call mapOneDKcat(dkVars%Aur_Type%data, AUR_TYPE)

    contains

        subroutine mapOneDK(dkField, mixID)
            real(rp), dimension(:,:), intent(in) :: dkField   ! (Nt,Np) active, (theta,phi)
            integer,                  intent(in) :: mixID
            integer :: i, j
            call initShellVar(voltGrid, SHGR_CC, sgv)
            do i = is_, ie_          ! theta (active; voltron grid has no theta ghosts)
                do j = js_, je_      ! phi (active)
                    sgv%data(i,j) = dkField(i - is_ + 1, j - js_ + 1)
                enddo
            enddo
            call wrapJ_SGV(voltGrid, sgv)   ! periodic phi ghost fill (replaces padDKPhiPeriodic)
            sgv%mask = .true.               ! DK conductance/precip valid over the closed active region
            call mapShellVarToMix(voltGrid, sgv, mixApp%ion, mixID)
        end subroutine mapOneDK

        ! Categorical (nearest-neighbour) interpolation for integer fields 
        subroutine mapOneDKcat(dkField, mixID)
            real(rp), dimension(:,:), intent(in) :: dkField
            integer,                  intent(in) :: mixID
            integer :: i, j
            call initShellVar(voltGrid, SHGR_CC, sgv)
            do i = is_, ie_
                do j = js_, je_
                    sgv%data(i,j) = dkField(i - is_ + 1, j - js_ + 1)
                enddo
            enddo
            call wrapJ_SGV(voltGrid, sgv)
            sgv%mask = .true.
            call mapShellVarToMixNN(voltGrid, sgv, mixApp%ion, mixID)  ! nearest, no averaging
        end subroutine mapOneDKcat

    end subroutine mapDKToMix

    subroutine mapNonDKToMix(vApp, mixApp)
        class(voltApp_T), intent(inout) :: vApp
        class(mixApp_T),  intent(inout) :: mixApp

        type(ShellGridVar_T) :: sinLon, cosLon
        integer  :: i, j, jM
        real(rp) :: x1, x2, latc_v, lonc_v, s, c

        ! Protons: bypassed (no DK / no RAIJU ion source) -> deterministic zero.
        mixApp%ion(NORTH)%St%Vars(:,:,IM_IFLUX) = 0.0_rp
        mixApp%ion(NORTH)%St%Vars(:,:,IM_IAVG ) = 0.0_rp
        mixApp%ion(NORTH)%St%Vars(:,:,IM_INFLX) = 0.0_rp
        mixApp%ion(SOUTH)%St%Vars(:,:,IM_IFLUX) = 0.0_rp
        mixApp%ion(SOUTH)%St%Vars(:,:,IM_IAVG ) = 0.0_rp
        mixApp%ion(SOUTH)%St%Vars(:,:,IM_INFLX) = 0.0_rp

        ! Conjugate geometry from RAIJU/tube mapping. Interpolate lon via sin/cos (wrap-safe).
        sinLon = vApp%State%tubeShell%lonc; sinLon%data = sin(vApp%State%tubeShell%lonc%data)
        cosLon = vApp%State%tubeShell%lonc; cosLon%data = cos(vApp%State%tubeShell%lonc%data)

        associate(GN => mixApp%ion(NORTH)%G, GS => mixApp%ion(SOUTH)%G, &
                  ion => mixApp%ion, sh => vApp%shGrid, ts => vApp%State%tubeShell)
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,jM,x1,x2,latc_v,lonc_v,s,c)
            do j = 1, GN%Np
                do i = 1, GN%Nt
                    ! NH REMIX cell center
                    x1 = GN%t(j,i); x2 = GN%p(j,i)
                    call InterpShellVar_TSC_pnt(sh, ts%latc, x1, x2, latc_v)
                    call InterpShellVar_TSC_pnt(sh, sinLon,  x1, x2, s)
                    call InterpShellVar_TSC_pnt(sh, cosLon,  x1, x2, c)
                    lonc_v = atan2(s, c); if (lonc_v < 0.0_rp) lonc_v = lonc_v + 2.0_rp*PI
                    ion(NORTH)%St%Vars(j,i,IM_THCON) = PI/2.0_rp - latc_v  ! conjugate colat [0,pi]
                    ion(NORTH)%St%Vars(j,i,IM_PHCON) = lonc_v              ! conjugate lon   [0,2pi]

                    ! SH REMIX cell center (reflect theta about the equator)
                    jM = GS%Np + 1 - j
                    x1 = PI - GS%t(jM,i); x2 = GS%p(jM,i)
                    call InterpShellVar_TSC_pnt(sh, ts%latc, x1, x2, latc_v)
                    call InterpShellVar_TSC_pnt(sh, sinLon,  x1, x2, s)
                    call InterpShellVar_TSC_pnt(sh, cosLon,  x1, x2, c)
                    lonc_v = atan2(s, c); if (lonc_v < 0.0_rp) lonc_v = lonc_v + 2.0_rp*PI
                    ion(SOUTH)%St%Vars(j,i,IM_THCON) = PI/2.0_rp - latc_v
                    ion(SOUTH)%St%Vars(j,i,IM_PHCON) = lonc_v
                end do
            end do
        end associate

    end subroutine mapNonDKToMix

    ! Imported GCM conductance overwrite in REMIX state space.
    ! Uses each GCM  shell-grid as source and maps into both REMIX hemispheres.
    subroutine mapGCMToMix(gcm, mixApp)
        class(gcmCoupler_T),    intent(inout) :: gcm
        class(mixApp_T),intent(inout) :: mixApp

        if (gcm%GEO%gcm2mix_nvar > 0) call mapOneGCMGridToMix(gcm%GEO, mixApp)
        if (gcm%APEX%gcm2mix_nvar > 0) call mapOneGCMGridToMix(gcm%APEX, mixApp)
    end subroutine mapGCMToMix

    subroutine mapOneGCMGridToMix(gcmGrid, mixApp)
        type(gcm_grid_T), intent(inout) :: gcmGrid
        class(mixApp_T),  intent(inout) :: mixApp

        integer :: v, h, mixVar

        if (gcmGrid%gcm2mix_nvar == 0) return

        do v = 1, gcmGrid%gcm2mix_nvar

            select case (trim(gcmGrid%inlist(v)))
            case ("SIGMAP")
                mixVar = SIGMAP
            case ("SIGMAH")
                mixVar = SIGMAH
            case default
                write(*,*) "MIXCPL WARNING: Unknown GCM inlist variable: ", trim(gcmGrid%inlist(v))
                cycle
            end select

            call mapShellVarToMix(gcmGrid%grid, gcmGrid%gcmInput(v), mixApp%ion, mixVar)
            do h = NORTH, SOUTH
                select case (mixVar)
                case (SIGMAP)
                    mixApp%ion(h)%St%Vars(:,:,SIGMAP) = max(mixApp%ion(h)%P%pedmin, &
                                                             mixApp%ion(h)%St%Vars(:,:,SIGMAP))
                case (SIGMAH)
                    mixApp%ion(h)%St%Vars(:,:,SIGMAH) = max(mixApp%ion(h)%P%hallmin, &
                                                             mixApp%ion(h)%St%Vars(:,:,SIGMAH))
                end select
            enddo
            write(*,*) "Mapped GCM variable ", trim(gcmGrid%inlist(v)), " to REMIX variable ", mixVar
        end do

        write(*,*) "Val Check NORTH: ", maxval(mixApp%ion(NORTH)%St%Vars(:,:,SIGMAP)), minval(mixApp%ion(NORTH)%St%Vars(:,:,SIGMAP))
        write(*,*) "Val Check SOUTH: ", maxval(mixApp%ion(SOUTH)%St%Vars(:,:,SIGMAP)), minval(mixApp%ion(SOUTH)%St%Vars(:,:,SIGMAP))
    end subroutine mapOneGCMGridToMix
end submodule remixCplTypesSub
