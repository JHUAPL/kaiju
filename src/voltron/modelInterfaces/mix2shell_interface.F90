! REMIX hemisphere grids -> voltron ShellGrid electrostatic potential (potential_base).
! Corotation and the total are voltron's job (updateVoltPotential).
module mix2shell_interface
    use kdefs
    use mixdefs,           only: POT
    use volttypes
    use mixinterfaceutils, only: mixVarToVoltron
    use shellGrid,         only: initShellVar, ShellGridVar_T, SHGR_CORNER

    implicit none
    private
    public :: mapMixToShell

    contains

    subroutine mapMixToShell(mixApp, vApp)
        type(mixApp_T),   intent(inout) :: mixApp
        class(voltApp_T), intent(inout) :: vApp

        type(ShellGridVar_T) :: iono_var

        call initShellVar(vApp%shGrid, SHGR_CORNER, iono_var)
        call mixVarToVoltron(mixApp, POT, vApp%shGrid, iono_var)

        ! mask is zero-filled equatorward of remix's low-lat boundary
        vApp%State%potential_base%data = iono_var%data
        vApp%State%potential_base%mask = iono_var%mask
    end subroutine mapMixToShell

end module mix2shell_interface
