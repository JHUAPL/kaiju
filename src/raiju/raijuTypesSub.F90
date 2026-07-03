submodule (raijutypes) raijuTypesSub

    use raijuStarter
    use raijuGrids
    use raijuIO
    use raijuOut
    use raijuAdvancer

    implicit none

    contains

    module subroutine raiInitModel(App, xml)
        class(raijuApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call raijuInit(App, Xml)
        ! Note: be careful about adding anything more here, because raiCpl_T's init doesn't call raiInitModel

    end subroutine raiInitModel


    module subroutine raiInitIO(App, Xml)
        class(raijuApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        ! Init output file
        call raijuInitIO(App%Model, App%Grid, App%Model%writeGhosts)

    end subroutine raiInitIO


    module subroutine raiWriteRestart(App, nRes)
        class(raijuApp_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! synchronize restart output number
        App%State%IO%nRes = nRes

        call raijuResOutput(App%Model, App%Grid, App%State, App%opt)

    end subroutine raiWriteRestart


    module subroutine raiReadRestart(App, resId, nRes)
        class(raijuApp_T), intent(inout) :: App
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
        !call GenShellGridFromFile(shRes, RAI_SG_NAME, App%Model%ResF)
        !if(.not. checkResGrid(App%Grid%shGrid, shRes)) then
        !    write(*,*)"RAIJU restart error: Grid generated from XML doesn't match that from restart file, that's not allowed"
        !    stop
        !endif

        ! Now read State info
        call raijuResInput(App%Model, App%Grid, App%State)

    end subroutine raiReadRestart


    module subroutine raiWriteConsoleOutput(App)
        class(raijuApp_T), intent(inout) :: App

        call raijuConsoleOut(App%Model, App%Grid, App%State)

    end subroutine raiWriteConsoleOutput


    module subroutine raiWriteFileOutput(App, nStep)
        class(raijuApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        ! synchronize file output number
        App%State%IO%nOut = nStep

        call raijuOutput(App%Model, App%Grid, App%State)

    end subroutine raiWriteFileOutput


    module subroutine raiWriteSlimFileOutput(App, nStep)
        class(raijuApp_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call raiWriteFileOutput(App, nStep)

    end subroutine raiWriteSlimFileOutput


    module subroutine raiAdvanceModel(App, dt)
        class(raijuApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt

        call raijuAdvance(App%Model, App%Grid, App%State, dt)

    end subroutine raiAdvanceModel


    module subroutine raiCleanup(App)
        class(raijuApp_T), intent(inout) :: App
        
        write(*,*) "RAIJU doing nothing for cleanup, idk what to do here yet"

    end subroutine raiCleanup

end submodule raijuTypesSub