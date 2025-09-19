submodule (mixtypes) mixtypessub
    use mixio

    
    implicit none

    contains

    module subroutine mixInitModel(App, Xml) !doubts with this one, since many things involved
        class(mixApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call mixInit(App, Xml)
        call init_mix(App%ion, [NORTH, SOUTH],Xml)
        ! call init_mix(I,hmsphrs,optFilename,RunID,isRestart,mixIOobj,nRes,optIO)
    end subroutine mixInitModel


    module subroutine mixInitIO(App, Xml) !doubts with this one too
        class(mixApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        ! Initialize IO file
        call initMIXIO(App%ion, runID, isRestart, App%ion%P%nRes)

    end subroutine mixInitIO


    module subroutine mixWriteRestart(App, nRes) !OK
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! synchronize restart output number
        App%ion%P%nRes = nRes

        call writeMIXRestart(App%ion,nRes,App%IO%mjd,App%IO%time)

    end subroutine mixWriteRestart


    module subroutine mixReadRestart(App, resId, nRes) !OK
        class(mixApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        ! synchronize restart number (needed?)
        App%ion%P%nRes = nRes

        call readMIXrestart(resId,nRes,App%ion)

    end subroutine mixReadRestart


    module subroutine mixWriteConsoleOutput(App) ! not sure about this one
        class(mixApp_T), intent(inout) :: App

    end subroutine mixWriteConsoleOutput


    module subroutine mixWriteFileOutput(App, nStep) !OK
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        ! synchronize file output number
        App%State%IO%nOut = nStep

        call writeMIX(App%ion,nStep,App%IO%mjd,App%IO%time)

    end subroutine mixWriteFileOutput


    module subroutine mixWriteSlimFileOutput(App, nStep) !OK
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call App%WriteFileOutput(App,nStep)  ! no separate slim output for remix

    end subroutine mixWriteSlimFileOutput


    module subroutine mixAdvanceModel(App, dt) !OK
        class(mixApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt  ! no need for time step here

        call run_mix(App%ion,App%IO%tilt)   ! ignoring doModelOpt and gcm (optional)

    end subroutine mixAdvanceModel


    module subroutine mixCleanup(App) !OK
        class(mixApp_T), intent(inout) :: App
        
        ! no cleanup for Remix

    end subroutine mixCleanup

end submodule mixtypessub