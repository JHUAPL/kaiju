submodule (dktypes) dktypessub
    use dkmain
    use dkinit
    use dkio
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

        character(len=strLen) :: runID
        logical :: isRestart
        integer :: nRes

        ! Same keys REMIX/Gamera use, so DK shares the run's RunID / restart cfg.
        call Xml%Set_Val(runID,     "/Kaiju/Gamera/sim/runid",     "msphere")
        call Xml%Set_Val(isRestart, "/Kaiju/gamera/restart/doRes", .false.)
        call Xml%Set_Val(nRes,      "/Kaiju/gamera/restart/nRes",  -1)

        call initDKIO(App, runID, isRestart, nRes)
    end subroutine dkInitIO


    module subroutine dkWriteRestart(App, nRes)
        class(dkApp_T), intent(inout) :: App
        integer, intent(in) :: nRes

        call writeDKRestart(App, nRes)
    end subroutine dkWriteRestart


    module subroutine dkReadRestart(App, resId, nRes)
        class(dkApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        call readDKRestart(App, resId, nRes)
    end subroutine dkReadRestart


    module subroutine dkWriteConsoleOutput(App)
        class(dkApp_T), intent(inout) :: App

    end subroutine dkWriteConsoleOutput


    module subroutine dkWriteFileOutput(App, nStep)
        class(dkApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call writeDKFileOutput(App, nStep)
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