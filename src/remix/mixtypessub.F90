submodule (mixtypes) mixtypessub
    use mixio
    use mixmain
    use xml_input
    
    implicit none

    contains

    module subroutine mixInitModel(App, Xml) !doubts with this one, since many things involved
        class(mixApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call init_mix(App%ion, [NORTH, SOUTH], Xml, App%IO)
        
    end subroutine mixInitModel


    module subroutine mixInitIO(App, Xml) ! OK!
        class(mixApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        character(len=strLen) :: runID   ! these two things are needed when we're coupled with Gamera
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


    module subroutine mixWriteRestart(App, nRes) ! OK!
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

        ! synchronize restart number
        App%ion(:)%P%nRes = nRes

        call readMIXrestart(resId,nRes,App%ion)

    end subroutine mixReadRestart


    module subroutine mixWriteConsoleOutput(App) ! OK! no console output
        class(mixApp_T), intent(inout) :: App

    end subroutine mixWriteConsoleOutput


    module subroutine mixWriteFileOutput(App, nStep) ! OK!
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call writeMIX(App%ion,nStep,App%IO%mjd,App%IO%time)

    end subroutine mixWriteFileOutput


    module subroutine mixWriteSlimFileOutput(App, nStep) ! OK!
        class(mixApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call App%WriteFileOutput(nStep)  ! no separate slim output for remix

    end subroutine mixWriteSlimFileOutput


    ! module subroutine mixAdvanceModel(App, dt) ! OK! ignoring gcm for now
    !     class(mixApp_T), intent(inout) :: App
    !     real(rp), intent(in) :: dt  ! no need for time step here

    !     if (App%IO%time<=0) then
    !         call run_mix(App%ion,App%IO%tilt,doModelOpt=.false.)   ! ignoring gcm
    !     else
    !         call run_mix(App%ion,App%IO%tilt,doModelOpt=.true.)   ! ignoring gcm
    !     end if

    ! end subroutine mixAdvanceModel

    module subroutine mixAdvanceModel(App, dt) ! OK! ignoring gcm for now
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


    module subroutine mixCleanup(App) ! OK! no cleanup
        class(mixApp_T), intent(inout) :: App

    end subroutine mixCleanup

end submodule mixtypessub