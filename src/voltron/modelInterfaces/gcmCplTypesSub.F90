submodule (volttypes) gcmCplTypesSub
    use shellInterp
    use mixdefs, only : NUM_FLUX, AVG_ENG
    use math
    use ioh5
    use mixio
    use kai2geo
    use clocks
    use gcmtypes

    implicit none

contains

    ! ----------------------------------------------------------------
    ! ionoCoupler_T deferred procedures
    ! ----------------------------------------------------------------

    ! GCM coupling is volt->GCM (export direction)
    ! interpolate Voltron shell fields onto each GCM grid.
    module subroutine gcmCpl_toIono(App, vApp)
        class(gcmCoupler_T), intent(inout) :: App
        class(voltApp_T),    intent(inout) :: vApp
        call gcmCpl_mapVoltToGCM(vApp%shGrid, vApp%State, vApp%dkApp%vars, App)
    end subroutine gcmCpl_toIono

    ! fromIono: GCM -> Volt
    ! This is currently performed inside of remix because we already loaded the packed the gcm array
    module subroutine gcmCpl_fromIono(App, vApp)
        class(gcmCoupler_T), intent(inout) :: App
        class(voltApp_T),    intent(inout) :: vApp
    end subroutine gcmCpl_fromIono

    ! No staging grid needed for GCM coupling at this stage.
    module subroutine gcmCpl_packInputs(App, vApp)
        class(gcmCoupler_T), intent(inout) :: App
        class(voltApp_T),    intent(inout) :: vApp
    end subroutine gcmCpl_packInputs

    module subroutine gcmCpl_InitCouplers(App, vApp, xml)
        class(gcmCoupler_T), intent(inout) :: App
        class(voltApp_T),    intent(inout) :: vApp
        type(XML_Input_T),   intent(inout) :: xml
    end subroutine gcmCpl_InitCouplers

    module subroutine gcmCpl_Coupling(App, vApp)
        class(gcmCoupler_T), intent(inout) :: App
        class(voltApp_T),    intent(inout) :: vApp
        ! No-op for base class; MPI version overrides
    end subroutine gcmCpl_Coupling

    ! ----------------------------------------------------------------
    ! BaseApp_T lifecycle stubs
    ! ----------------------------------------------------------------

    module subroutine gcmCplInitModel(App, xml)
        class(gcmCoupler_T), intent(inout) :: App
        type(XML_Input_T),   intent(inout) :: xml
    end subroutine gcmCplInitModel

    module subroutine gcmCplInitIO(App, xml)
        class(gcmCoupler_T), intent(inout) :: App
        type(XML_Input_T),   intent(inout) :: xml
    end subroutine gcmCplInitIO

    module subroutine gcmCplWriteRestart(App, nRes)
        class(gcmCoupler_T), intent(inout) :: App
        integer,             intent(in)    :: nRes
    end subroutine gcmCplWriteRestart

    module subroutine gcmCplReadRestart(App, resId, nRes)
        class(gcmCoupler_T), intent(inout) :: App
        character(len=*),    intent(in)    :: resId
        integer,             intent(in)    :: nRes
    end subroutine gcmCplReadRestart

    module subroutine gcmCplWriteConsoleOutput(App)
        class(gcmCoupler_T), intent(inout) :: App
    end subroutine gcmCplWriteConsoleOutput

    module subroutine gcmCplWriteFileOutput(App, nStep)
        class(gcmCoupler_T), intent(inout) :: App
        integer,             intent(in)    :: nStep
    end subroutine gcmCplWriteFileOutput

    module subroutine gcmCplWriteSlimFileOutput(App, nStep)
        class(gcmCoupler_T), intent(inout) :: App
        integer,             intent(in)    :: nStep
    end subroutine gcmCplWriteSlimFileOutput

    module subroutine gcmCplAdvanceModel(App, dt)
        class(gcmCoupler_T), intent(inout) :: App
        real(rp),            intent(in)    :: dt
    end subroutine gcmCplAdvanceModel

    module subroutine gcmCplCleanup(App)
        class(gcmCoupler_T), intent(inout) :: App
    end subroutine gcmCplCleanup

    ! ----------------------------------------------------------------
    ! Core implementation
    ! ----------------------------------------------------------------

    ! Interpolate Voltron shell fields onto each active GCM grid (GEO and APEX),
    ! dispatching on the string token in outlist(v).
    subroutine gcmCpl_mapVoltToGCM(shGrid, State, dkVars, gcm)
        type(ShellGrid_T),  intent(in)    :: shGrid
        type(voltState_t),  intent(in)    :: State
        type(dkVariables_T),intent(in)    :: dkVars
        class(gcmCoupler_T),        intent(inout) :: gcm

        integer :: v

        if (gcm%GEO%mix2gcm_nvar > 0) then
            do v = 1, gcm%GEO%mix2gcm_nvar
                call map_one(trim(gcm%GEO%outlist(v)), gcm%GEO%grid, gcm%GEO%gcmOutput(v))
            end do
        end if

        if (gcm%APEX%mix2gcm_nvar > 0) then
            do v = 1, gcm%APEX%mix2gcm_nvar
                call map_one(trim(gcm%APEX%outlist(v)), gcm%APEX%grid, gcm%APEX%gcmOutput(v))
            end do
        end if

    contains

        subroutine map_one(varName, dstGrid, dst)
            character(len=*),     intent(in)    :: varName
            type(ShellGrid_T),    intent(in)    :: dstGrid
            type(ShellGridVar_T), intent(inout) :: dst

            ! DK vars are now ShellGridVar_T on the shell grid (mask defaults .true., phi wrap
            ! handled inside the interp), so they map directly -- no repack, same as POT.
            select case (varName)
            case ("POT")
                call InterpShellVar_TSC_SG(shGrid, State%potential_base,   dstGrid, dst)
            case ("NUM_FLUX")
                call InterpShellVar_TSC_SG(shGrid, dkVars%Num_Flux, dstGrid, dst)
            case ("AVG_ENG")
                call InterpShellVar_TSC_SG(shGrid, dkVars%Avg_Eng, dstGrid, dst)
            case default
                dst%data = 0._rp
                write(*,*) "GCMCPL WARNING: Unknown outlist variable: ", varName
            end select
        end subroutine map_one

    end subroutine gcmCpl_mapVoltToGCM

end submodule gcmCplTypesSub
