submodule (dktypes) dktypessub
    use dkaurora
    ! use dkmain
    use xml_input
    
    implicit none

    contains

    module subroutine dkInitModel(App, Xml)
        class(dkApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        ! call dragonking_init(App%aurora,App%Params,App%G)
        
    end subroutine dkInitModel


    module subroutine dkInitIO(App, Xml)
        class(dkApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        character(len=strLen) :: runID
        logical :: isRestart
        integer :: nRes

        ! call Xml%Set_Val(runID, "/Kaiju/Gamera/sim/runid","msphere")
        ! call Xml%Set_Val(isRestart, "/Kaiju/Gamera/restart/doRes",.false.)
        ! call Xml%Set_Val(nRes, "/Kaiju/Gamera/restart/nRes",-1)

        ! ! synchronize restart output number
        ! App%ion(:)%P%nRes = nRes

        ! ! Initialize IO file
        ! call initdkIO(App, runID, isRestart, nRes)

    end subroutine dkInitIO


    module subroutine dkWriteRestart(App, nRes)
        class(dkApp_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! ! synchronize restart output number
        ! App%ion(:)%P%nRes = nRes

        ! call writedkRestart(App%ion,nRes,App%IO%mjd,App%IO%time)

    end subroutine dkWriteRestart


    module subroutine dkReadRestart(App, resId, nRes) ! OK!
        class(dkApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        ! For some reason, Redk uses the internal restart number
        ! And ignores the one provided by Voltron
        ! App%ion(:)%P%nRes != nRes

        ! call readdkrestart(resId,App%ion(1)%P%nRes,App%ion)

    end subroutine dkReadRestart


    module subroutine dkWriteConsoleOutput(App)
        class(dkApp_T), intent(inout) :: App

    end subroutine dkWriteConsoleOutput


    module subroutine dkWriteFileOutput(App, nStep)
        class(dkApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        ! call writedk(App%ion,nStep,App%IO%mjd,App%IO%time)

    end subroutine dkWriteFileOutput


    module subroutine dkWriteSlimFileOutput(App, nStep)
        class(dkApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        ! call App%WriteFileOutput(nStep)  ! no separate slim output for redk

    end subroutine dkWriteSlimFileOutput

    module subroutine dkAdvanceModel(App, dt)
        class(dkApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt  ! no need for time step here
        ! integer :: h,NumH

        ! NumH = size(App%ion)

        ! do h=1,NumH
        !   call Tic("dk-SOLVE")
        !   call run_solver(App%ion(h)%P,App%ion(h)%G,App%ion(h)%St,App%ion(h)%S)
        !   call Toc("dk-SOLVE")

        !   call Tic("dk-POT")
        !   call get_potential(App%ion(h))
        !   call Toc("dk-POT")
        ! end do
        
    end subroutine dkAdvanceModel


    module subroutine dkCleanup(App)
        class(dkApp_T), intent(inout) :: App

    end subroutine dkCleanup

end submodule dktypessub