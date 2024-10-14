!Driver for Gamera (uncoupled)

program gamerax
    use gamapp
    use usergamic

    implicit none

    type(gamApp_T) :: gApp
    character(len=strLen) :: inpXML
    type(XML_Input_T) :: xmlInp
    real(rp) :: nextDT

    ! set options for gamera app
    gApp%gOptions%userInitFunc => initUser

    !call printConfigStamp()
    call initClocks()

    call getIDeckStr(inpXML)
    call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

    !Create XML reader
    write(*,*) 'Reading input deck from ', trim(inpXML)
    xmlInp = New_XML_Input(trim(inpXML),'Kaiju',.true.)

    call gApp%InitModel(xmlInp)
    call gApp%InitIO(xmlInp)

    do while (gApp%Model%t < gApp%Model%tFin)
        call Tic("Omega") !Start root timer

        !Step model
        nextDT = min(gApp%Model%tFin-gApp%Model%t, gApp%Model%IO%nextIOTime()-gApp%Model%t)
        call gApp%AdvanceModel(nextDT)

        !Output if necessary
        call Tic("IO")

        if (gApp%Model%IO%doConsole(gApp%Model%t)) then
            call gApp%WriteConsoleOutput()

            !Timing info
            if (gApp%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()

        elseif (gApp%Model%IO%doTimer(gApp%Model%t)) then
            if (gApp%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif

        if (gApp%Model%IO%doOutput(gApp%Model%t)) then
            call gApp%WriteFileOutput(gApp%Model%IO%nOut)
        endif

        if (gApp%Model%IO%doRestart(gApp%Model%t)) then
            call gApp%WriteRestart(gApp%Model%IO%nRes)
        endif

        call Toc("IO")
        call Toc("Omega")
    end do

   write(*,*) "Fin"

end program gamerax

