submodule (mixtypes) mixtypessub
    use mixio
    use mixmain
    use xml_input
    
    implicit none

    contains

    module subroutine mixInitModel(App, Xml)
        class(mixApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call init_mix(App%ion, [NORTH, SOUTH], Xml, App%IO)
        
    end subroutine mixInitModel


    module subroutine mixInitIO(App, Xml)
        class(mixApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        character(len=strLen) :: runID
        logical :: isRestart
        integer :: nRes

        call Xml%Set_Val(runID, "/Kaiju/Gamera/sim/runid","msphere")
        call Xml%Set_Val(isRestart, "/Kaiju/Gamera/restart/doRes",.false.)
        call Xml%Set_Val(nRes, "/Kaiju/Gamera/restart/nRes",-1)

        ! synchronize restart output number
        App%ion(:)%P%nRes = nRes

        ! Initialize IO file
        call initMIXIO(App, runID, isRestart, nRes)

    end subroutine mixInitIO


    module subroutine mixWriteRestart(App, nRes)
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! synchronize restart output number
        App%ion(:)%P%nRes = nRes

        call writeMIXRestart(App%ion,nRes,App%IO%mjd,App%IO%time)

    end subroutine mixWriteRestart


    module subroutine mixReadRestart(App, resId, nRes) ! OK!
        class(mixApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        ! For some reason, Remix uses the internal restart number
        ! And ignores the one provided by Voltron
        ! App%ion(:)%P%nRes != nRes

        call readMIXrestart(resId,App%ion(1)%P%nRes,App%ion)

    end subroutine mixReadRestart


    module subroutine mixWriteConsoleOutput(App)
        class(mixApp_T), intent(inout) :: App

    end subroutine mixWriteConsoleOutput


    module subroutine mixWriteFileOutput(App, nStep)
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call writeMIX(App%ion,nStep,App%IO%mjd,App%IO%time)

    end subroutine mixWriteFileOutput


    module subroutine mixWriteSlimFileOutput(App, nStep)
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call App%WriteFileOutput(nStep)  ! no separate slim output for remix

    end subroutine mixWriteSlimFileOutput

    module subroutine mixAdvanceModel(App, dt)
        class(mixApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt  ! no need for time step here
        integer :: h,NumH

        NumH = size(App%ion)

        do h=1,NumH
          call Tic("MIX-SOLVE")
          call run_solver(App%ion(h)%P,App%ion(h)%G,App%ion(h)%St,App%ion(h)%S)
          call Toc("MIX-SOLVE")

          call Tic("MIX-POT")
          call get_potential(App%ion(h))
          call Toc("MIX-POT")
        end do
        
    end subroutine mixAdvanceModel


    module subroutine mixCleanup(App)
        class(mixApp_T), intent(inout) :: App

    end subroutine mixCleanup

end submodule mixtypessub