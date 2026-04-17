! Submodule implementing remixCoupler_T procedures.
! Mirrors raijuCplTypesSub.F90 in structure.
! All computation is identical to the pre-refactor code — only the
! data path changes: vApp%remixApp -> App%mixApp.

submodule (volttypes) remixCplTypesSub

    use mixtypes
    use mixdefs
    use mhd2mix_interface
    use mix2mhd_interface
    use imag2mix_interface
    use mixinterfaceutils
    use voltappHelper
    use voltapp, only: getConductance
    use kai2geo

    implicit none

    contains

    ! ----------------------------------------------------------------
    ! InitModel — REMIX-internal init only.
    ! The map-building calls (init_mhd2Mix, init_mix2MhdCoupler) that
    ! need gApp stay in voltapp.F90/initializeFromGamera for Step 1.
    ! ----------------------------------------------------------------
    module subroutine remixCplInitModel(App, xml)
        class(remixCoupler_T), intent(inout) :: App
        type(XML_Input_T),     intent(inout) :: xml

        call App%mixApp%InitModel(xml)

    end subroutine remixCplInitModel


    ! ----------------------------------------------------------------
    ! toIono — pack Gamera + IMAG data into REMIX.
    ! Contains the logic previously in runRemix() plus
    ! convertGameraToRemix() which was called just before it in
    ! DeepUpdate().
    ! ----------------------------------------------------------------
    module subroutine remixCpl_toIono(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        real(rp) :: curTilt

        ! Gamera → mixInput (unit conversions, shell J extraction)
        call convertGameraToRemix(vApp%mhd2mix, vApp%gApp, App%mixApp)

        ! Recalculate MJD for current time
        call MJDRecalc(vApp%MJD)

        ! IMAG precipitation → REMIX (only after spinup, only when deep)
        if (vApp%doDeep .and. (vApp%time > 0)) then
            call CoupleIMagToMix(vApp, App%mixApp)
        endif

        ! Interpolate mixInput onto REMIX grid
        call mapGameraToRemix(vApp%mhd2mix, App%mixApp)

        ! Set tilt and time for conductance models
        call vApp%tilt%getValue(vApp%time, curTilt)
        App%mixApp%IO%tilt = curTilt
        App%mixApp%IO%time = vApp%time

        ! Compute conductances (Pedersen, Hall, parallel)
        if (vApp%doGCM) then
            call getConductance(App%mixApp%ion, App%mixApp%IO, vApp%gcm)
        else
            call getConductance(App%mixApp%ion, App%mixApp%IO)
        endif

    end subroutine remixCpl_toIono


    ! ----------------------------------------------------------------
    ! fromIono — map solved REMIX potential back to Voltron + Gamera.
    ! Contains the logic previously in updateVoltPotential() and
    ! CouplePotentialToMhd().
    ! ----------------------------------------------------------------
    module subroutine remixCpl_fromIono(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        type(ShellGridVar_T) :: iono_var

        ! REMIX potential → Voltron ShellGrid
        call initShellVar(vApp%shGrid, SHGR_CORNER, iono_var)
        call calcCorotPotential(vApp%planet, vApp%shGrid, vApp%State%potential_corot, &
                                doGeoCorotO=vApp%doGeoCorot)
        call mixVarToVoltron(App%mixApp, POT, vApp%shGrid, iono_var)
        vApp%State%potential_total%data = iono_var%data + vApp%State%potential_corot%data
        vApp%State%potential_total%mask = .true.

        ! Store MJD/time now so I/O calls need no special setup
        App%mixApp%IO%mjd  = vApp%MJD
        App%mixApp%IO%time = vApp%time

        ! Voltron potential → Gamera inner boundary
        call mapRemixToGamera(vApp%gApp, App%mixApp)
        call convertRemixToGamera(vApp%gApp, App%mixApp)

    end subroutine remixCpl_fromIono


    ! ----------------------------------------------------------------
    ! AdvanceModel — steady-state ionospheric potential solve.
    ! ----------------------------------------------------------------
    module subroutine remixCplAdvanceModel(App, dt)
        class(remixCoupler_T), intent(inout) :: App
        real(rp),              intent(in)    :: dt

        call App%mixApp%AdvanceModel(dt)

    end subroutine remixCplAdvanceModel


    ! ----------------------------------------------------------------
    ! Remaining BaseApp_T procedures — delegate to mixApp.
    ! ----------------------------------------------------------------
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

end submodule remixCplTypesSub
