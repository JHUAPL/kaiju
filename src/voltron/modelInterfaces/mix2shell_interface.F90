! Build the ionospheric potential on the voltron ShellGrid from REMIX's hemisphere grids:
! electrostatic part (mixVarToVoltron) + analytic corotation, kept separate so shell2gam
! can TSC-interp the smooth electrostatic part.
module mix2shell_interface
    use kdefs
    use mixdefs,           only: POT
    use volttypes
    use mixinterfaceutils, only: mixVarToVoltron
    use voltappHelper,     only: calcCorotPotential
    use shellGrid,         only: initShellVar, ShellGridVar_T, SHGR_CORNER

    implicit none
    private
    public :: mapMixToShell

    contains

    ! Total shell potential = electrostatic (REMIX -> shell) + corotation; corotation
    ! kept separate in potential_corot so shell2gam can subtract and re-add it sharply.
    subroutine mapMixToShell(mixApp, vApp)
        type(mixApp_T),   intent(inout) :: mixApp
        class(voltApp_T), intent(inout) :: vApp

        type(ShellGridVar_T) :: iono_var

        ! Corotation potential on the shell (kept separate for shell2gam).
        call initShellVar(vApp%shGrid, SHGR_CORNER, iono_var)
        call calcCorotPotential(vApp%planet, vApp%shGrid, vApp%State%potential_corot, &
                                doGeoCorotO=vApp%doGeoCorot)

        ! Electrostatic REMIX potential onto the voltron full-sphere shell.
        call mixVarToVoltron(mixApp, POT, vApp%shGrid, iono_var)

        ! Total = electrostatic + corotation.
        vApp%State%potential_total%data = iono_var%data + vApp%State%potential_corot%data
        vApp%State%potential_total%mask = .true.
    end subroutine mapMixToShell

end module mix2shell_interface
