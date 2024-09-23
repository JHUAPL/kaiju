submodule (gamtypes) gamtypessub
    use gamapp

    implicit none

    contains

    ! procedures for gamApp_T, expanded on from src/base/types/gamtypes.F90
    module subroutine gamInitModel(App, Xml)
        class(gamApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call initGamera(App, Xml)

    end subroutine gamInitModel

    module subroutine gamInitIO(App, Xml)
        class(gamApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call initGamIO(App%Model, Xml)

    end subroutine gamInitIO

    module subroutine gamWriteRestart(App, nRes)
        class(gamApp_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! synchronize restart output number
        App%Model%IO%nRes = nRes

        call resOutput(App%Model, App%Grid, App%oState, App%State)

    end subroutine gamWriteRestart

    module subroutine gamReadRestart(App, resId, nRes)
        class(gamApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        character(len=strLen) :: inH5

        call getRestart(App%Model,App%Grid,resId,nRes,inH5)
        call readH5Restart(App%Model, App%Grid, App%State, App%oState,inH5)

    end subroutine gamReadRestart

    module subroutine gamWriteConsoleOutput(App)
        class(gamApp_T), intent(inout) :: App

        call consoleOutput(App%Model, App%Grid, App%State)

    end subroutine gamWriteConsoleOutput

    module subroutine gamWriteFileOutput(App, nStep)
        class(gamApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        ! synchronize file output number
        App%Model%IO%nOut = nStep

        call fOutput(App%Model, App%Grid, App%State)

    end subroutine gamWriteFileOutput

    module subroutine gamWriteSlimFileOutput(App, nStep)
        class(gamApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call App%WriteFileOutput(nStep)

    end subroutine gamWriteSlimFileOutput

    module subroutine gamAdvanceModel(App, dt)
        class(gamApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt

        real(rp) :: targetSimT 

        targetSimT = App%model%t+dt

        ! ensure gamera is stepped at least one time
        ! fortran doesn't allow checking at the end of a while loop
        call stepGamera(App)

        do while(App%model%t < targetSimT)
            call stepGamera(App)
        enddo

    end subroutine gamAdvanceModel

    module subroutine gamCleanup(App)
        class(gamApp_T), intent(inout) :: App
        ! nothing to cleanup in non-mpi gamera
    end subroutine gamCleanup
 
end submodule

