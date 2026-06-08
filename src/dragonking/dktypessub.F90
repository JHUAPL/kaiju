submodule (dktypes) dktypessub
    use dkmain
    use dkinit
    use xml_input
    
    implicit none

    contains

    module subroutine dkInitModel(App, Xml)
        class(dkApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call initDK(App, Xml)        
    end subroutine dkInitModel


    module subroutine dkInitIO(App, Xml)
        class(dkApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

    end subroutine dkInitIO


    module subroutine dkWriteRestart(App, nRes)
        class(dkApp_T), intent(inout) :: App
        integer, intent(in) :: nRes

    end subroutine dkWriteRestart


    module subroutine dkReadRestart(App, resId, nRes)
        class(dkApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

    end subroutine dkReadRestart


    module subroutine dkWriteConsoleOutput(App)
        class(dkApp_T), intent(inout) :: App

    end subroutine dkWriteConsoleOutput


    module subroutine dkWriteFileOutput(App, nStep)
        class(dkApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

    end subroutine dkWriteFileOutput


    module subroutine dkWriteSlimFileOutput(App, nStep)
        class(dkApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

    end subroutine dkWriteSlimFileOutput

    module subroutine dkAdvanceModel(App, dt)
        class(dkApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt  ! no need for time step here

        call DKcompute(App)
        
    end subroutine dkAdvanceModel


    module subroutine dkCleanup(App)
        class(dkApp_T), intent(inout) :: App

    end subroutine dkCleanup

end submodule dktypessub