! Module containing base types from which other models will inherit

module basetypes

    use gdefs
    use ioclock
    use ioH5
    use xml_input

    implicit none

    type, abstract :: BaseApp_T
        contains
        ! the base App intentionally does NOT contain instances of the base model/grid/state
        ! this allows the specific sub-apps to include their specific instances of model/grid/state
        ! including the option to not use them, or have multiple copies of them

        ! functions which will be over-written by sub-apps
        procedure(InitModel_interface), deferred :: InitModel
        procedure(InitIO_interface), deferred :: InitIO
        procedure(WriteRestart_interface), deferred :: WriteRestart
        procedure(ReadRestart_interface), deferred :: ReadRestart
        procedure(WriteConsoleOutput_interface), deferred :: WriteConsoleOutput
        procedure(WriteFileOutput_interface), deferred :: WriteFileOutput
        procedure(WriteSlimFileOutput_interface), deferred :: WriteSlimFileOutput
        procedure(AdvanceModel_interface), deferred :: AdvanceModel
        procedure(Cleanup_interface), deferred :: Cleanup

    end type BaseApp_T

    type, abstract :: BaseOptions_T
        contains
        ! this base class contains nothing, but will be extended and passed into Init by other classes
    end type BaseOptions_T

    ! helper type that can contain a pointer to any type of BaseApp
    ! useful for storing arrays or other structures of BaseApps
    type AppHolder
        class(BaseApp_T), allocatable :: p
    end type AppHolder

    abstract interface
        subroutine InitModel_interface(App, Xml)
            import BaseApp_T ! I don't want to talk about it
            import XML_Input_T
            class(BaseApp_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine InitModel_interface

        subroutine InitIO_interface(App, Xml)
            import BaseApp_T
            import XML_Input_T
            class(BaseApp_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine InitIO_interface

        subroutine WriteRestart_interface(App,nRes)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
            integer, intent(in) :: nRes
        end subroutine WriteRestart_interface

        subroutine ReadRestart_interface(App,resId,nRes)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
            character(len=*), intent(in) :: resId
            integer, intent(in) :: nRes
        end subroutine ReadRestart_interface

        subroutine WriteConsoleOutput_interface(App)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
        end subroutine WriteConsoleOutput_interface

        subroutine WriteFileOutput_interface(App,nStep)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine WriteFileOutput_interface

        subroutine WriteSlimFileOutput_interface(App,nStep)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine WriteSlimFileOutput_interface

        subroutine AdvanceModel_interface(App, dt)
            import BaseApp_T
            import rp
            class(BaseApp_T), intent(inout) :: App
            real(rp), intent(in) :: dt
        end subroutine AdvanceModel_interface

        subroutine Cleanup_interface(App)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
        end subroutine
    end interface

end module

