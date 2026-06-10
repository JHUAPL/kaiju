!remixCoupler_T procedures (mirrors raijuCplTypesSub.F90 structure).
!Data path: vApp%remixApp -> App%mixApp.

submodule (volttypes) remixCplTypesSub

    use mixtypes
    use mixdefs
    use mhd2mix_interface
    use mix2mhd_interface
    use imag2mix_interface
    use mixinterfaceutils
    use shell2mix_interface,  only: mapShellToMix, init_shell2mix
    use mix2shell_interface,  only: mapMixToShell
    use voltappHelper
    use kai2geo
    use mixmain,   only: updateF107, mix_transform_grid  ! not host-associated via volttypes
    use dkmain,    only: getMixToDK, getDKToMix, DKdogcm  ! DragonKing conductance/precip flow
    use dkinit,    only: initDKmesh

    implicit none

    contains

    !InitModel — REMIX-internal init only; gApp-dependent map builds stay in voltapp.
    module subroutine remixCplInitModel(App, xml)
        class(remixCoupler_T), intent(inout) :: App
        type(XML_Input_T),     intent(inout) :: xml

        call App%mixApp%InitModel(xml)

    end subroutine remixCplInitModel


    !InitCouplers — REMIX-specific init behind the ionoCoupler contract: restart/IO,
    !F10.7, ionosphere radii, and the gApp-dependent map builds (mix2mhd, mhd2mix,
    !shell2mix). The driver calls this polymorphically (no select type), so a future
    !smiteCoupler_T supplies its own without touching voltapp. (gam2shell, the
    !model-agnostic Gamera->shell init, stays in the driver.)
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

        ! DragonKing init: build the DK grid from the REMIX ion grids, then init the model.
        call initDKmesh(App%dkApp, App%mixApp%ion)
        call App%dkApp%InitModel(xml)

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
            App%dkApp%params%f107 = maxF107
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

        ! (Direct mhd2mix path retired: its maps were consumed only by the removed
        !  convertGameraToRemix/mapGameraToRemix, so init_mhd2Mix is no longer needed.)

        ! Coupler staging layer: deep-copy voltron ShellGrid, init staging vars.
        App%shGr = vApp%shGrid
        call initShellVar(App%shGr, SHGR_CORNER, App%fac)
        call initShellVar(App%shGr, SHGR_CORNER, App%density)
        call initShellVar(App%shGr, SHGR_CORNER, App%soundSpd)
        call initShellVar(App%shGr, SHGR_CORNER, App%pot_iono)

        ! Coupler ShellGrid -> REMIX's two hemisphere grids.
        call init_shell2mix(App%shell2mix, App%shGr, App%mixApp)

    end subroutine remixCpl_InitCouplers


    !toIono — pack Gamera + IMAG into REMIX (was runRemix + convertGameraToRemix in DeepUpdate).
    module subroutine remixCpl_toIono(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        real(rp) :: curTilt

        ! Recalculate MJD for current time
        call MJDRecalc(vApp%MJD)

        ! IMAG precipitation → REMIX (only after spinup, only when deep)
        if (vApp%doDeep .and. (vApp%time > 0)) then
            call CoupleIMagToMix(vApp, App%mixApp)
        endif

        ! voltron ShellGrid -> REMIX hemisphere grids.
        ! (gam2shell, Gamera -> shell, runs in the driver before this.)
        call App%packInputs(vApp)

        ! Set tilt and time for conductance models
        call vApp%tilt%getValue(vApp%time, curTilt)
        App%mixApp%IO%tilt = curTilt
        App%mixApp%IO%time = vApp%time

        ! DragonKing conductance + precipitation (replaces REMIX-internal getConductance):
        ! mix St%Vars -> DK -> compute -> write SigmaP/SigmaH/precip back into St%Vars.
        call getMixToDK(App%mixApp, App%dkApp)
        App%dkApp%IO%tilt = curTilt
        App%dkApp%IO%time = vApp%time
        if (vApp%doGCM) call DKdogcm(App%dkApp, vApp%gcm)
        call App%dkApp%AdvanceModel(0.0_rp)
        call getDKToMix(App%mixApp, App%dkApp)

        ! Geographic grid for the REMIX potential solve (non-GCM path).
        if (.not. vApp%doGCM) call mix_transform_grid(App%mixApp%ion)

    end subroutine remixCpl_toIono


    !fromIono — map solved REMIX potential back to Voltron + Gamera.
    module subroutine remixCpl_fromIono(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        ! REMIX hemisphere grids -> voltron ShellGrid.
        ! (shell2gam, shell -> Gamera inner boundary, runs in the driver after this.)
        call mapMixToShell(App%mixApp, vApp)

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


    !packInputs — voltron ShellGrid -> coupler staging (direct copy) -> REMIX (shell2mix).
    module subroutine remix_packInputs(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        ! voltron ShellGrid -> coupler staging ShellGrid.
        App%fac%data      = vApp%gam2shell%fac%data
        App%fac%mask      = vApp%gam2shell%fac%mask
        App%density%data  = vApp%gam2shell%density%data
        App%density%mask  = vApp%gam2shell%density%mask
        App%soundSpd%data = vApp%gam2shell%soundSpd%data
        App%soundSpd%mask = vApp%gam2shell%soundSpd%mask

        ! coupler staging ShellGrid -> REMIX's two hemisphere grids.
        call mapShellToMix(App%shell2mix, App, vApp%gam2shell, App%mixApp)

    end subroutine remix_packInputs



end submodule remixCplTypesSub
