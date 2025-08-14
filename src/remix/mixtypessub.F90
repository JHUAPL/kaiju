submodule (mixtypes) mixtypessub
    use mixio

    
    implicit none

    contains

    module subroutine mixInitModel(App, xml)
        class(mixApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call mixjuInit(App, Xml)
        ! Note: be careful about adding anything more here, because mixCpl_T's init doesn't call mixInitModel

    end subroutine mixInitModel


    module subroutine mixInitIO(App, Xml)
        class(mixApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        ! Initialize IO file
        call initMIXIO(App%ion, runID, isRestart, nRes) !doubt about the last 3

    end subroutine mixInitIO


    module subroutine mixWriteRestart(App, nRes)
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! synchronize restart output number
        App%ion%P%nRes = nRes

        call writeMIXRestart(App%ion,Step,mjd,time)

    end subroutine mixWriteRestart


    module subroutine mixReadRestart(App, resId, nRes)
        class(mixApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        type(ShellGrid_T) :: shRes
            !! ShellGrid saved to file. Just use to make sure we are using the same grid as was used previously

        ! synchronize restart output number
        App%State%IO%nRes = nRes
        ! Build restart filename
        App%Model%nResIn = nRes
        call genResInFname(App%Model, App%Model%ResF, runIdO=resId)
        ! Handle grid reading first
        call GenShellGridFromFile(shRes, RAI_SG_NAME, App%Model%ResF)
        if(.not. checkResGrid(App%Grid%shGrid, shRes)) then
            write(*,*)"RAIJU restart error: Grid generated from XML doesn't match that from restart file, that's not allowed"
            stop
        endif

        ! Now read State info
        call mixjuResInput(App%Model, App%Grid, App%State)

    end subroutine mixReadRestart


    module subroutine mixWriteConsoleOutput(App)
        class(mixApp_T), intent(inout) :: App

        call mixjuConsoleOut(App%Model, App%Grid, App%State)

    end subroutine mixWriteConsoleOutput


    module subroutine mixWriteFileOutput(App, nStep)
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        ! synchronize file output number
        App%State%IO%nOut = nStep

        call mixjuOutput(App%Model, App%Grid, App%State)

    end subroutine mixWriteFileOutput


    module subroutine mixWriteSlimFileOutput(App, nStep)
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call mixWriteFileOutput(App, nStep)

    end subroutine mixWriteSlimFileOutput


    module subroutine mixAdvanceModel(App, dt)
        class(mixApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt

        call mixjuAdvance(App%Model, App%Grid, App%State, dt)

    end subroutine mixAdvanceModel


    module subroutine mixCleanup(App)
        class(mixApp_T), intent(inout) :: App
        
        write(*,*) "RAIJU doing nothing for cleanup, idk what to do here yet"

    end subroutine mixCleanup

end submodule mixtypessub