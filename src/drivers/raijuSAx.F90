program raijuSAx
    !! Stand-alone RAIJU driver

    use kdefs
    use xml_input
    use clocks
    use geopack, only : mjdRECALC

    ! RAIJU stuff
    use raijudefs
    use raijustarter
    use raijuBCs
    use raijuPreAdvancer, only : raijuPreAdvance
    use raijuAdvancer   , only : raijuAdvance

    implicit none

    integer :: testCode

    type(raijuApp_T   ) :: raiApp

    character(len=strLen) :: XMLStr, gStr
    type(XML_Input_T) :: inpXML    
    
    character(len=strLen) :: FLH5
    
    logical :: doChmpOut,doFLOut
    logical :: doPosFix = .false.    

    real(rp) :: mjd0



    call initClocks()
    call Tic("Omega")
                    
    ! Init xml
    call getIDeckStr(XMLStr)

    ! Personal xml flags
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/RAIJU",.true.)

    !! Break-out testing of specific pieces
    call inpXML%Set_Val(testCode,'testing/code',-1)
    select case(testCode)
        case (-1)  ! Nothing, just pass
            write(*,*) "No test code, doing standard SA driving"
        case default
            write(*,*) "Idk what this test code is, try again"
            stop
    end select

    call inpXML%Set_Val(raiApp%Model%doClockConsoleOut,'driver/doClockOut',.false.)

    call inpXML%Set_Val(doPosFix     ,'hax/doPosFix',doPosFix     )

    ! Init RAIJU
    call raiApp%InitModel(inpXML)
    call raiApp%InitIO(inpXML)
    raiApp%State%isFirstCpl = .false.
    if (raiApp%Model%isRestart) then
        call raiApp%ReadRestart(raiApp%Model%resID, raiApp%Model%nResIn)
        raiApp%State%isFirstCpl = .false.
    endif
    

    !> MJD at the start of the simulation (corresponds to sim t=0)
    call inpXML%Set_Val(mjd0,'prob/MJD0',51544.0)  ! default to 2000-01-01T00:00:00
    if (.not. raiApp%Model%isRestart) then
        raiApp%State%mjd = mjd0
    endif

    ! If we need geopack, make sure we init
    if (raiApp%Model%doGeoCorot) then
        call mjdRECALC(raiApp%State%mjd)
    endif

    ! Before we get started, run preAdvancer once just so that Step#0 has that information
    ! call raijuPreAdvance(raiApp%Model, raiApp%Grid, raiApp%State,isfirstCplO=isfirstCpl)
    !isfirstCpl = .false.

    ! Ready to loop
    do while ( raiApp%State%t < (raiApp%Model%tFin + 0.5) )
        
        call Tic("Output")
    ! Output if ready
        if (raiApp%State%IO%doRestart(raiApp%State%t)) then
            call raiApp%WriteRestart(raiApp%State%IO%nRes)
            !call raijuResOutput(raiApp%Model,raiApp%Grid,raiApp%State)
            !call raijuResInput(raiApp%Model,raiApp%Grid,raiApp%State)
        endif

        if (raiApp%State%IO%doOutput(raiApp%State%t)) then
            call raiApp%WriteFileOutput(raiApp%State%IO%nOut)
            !call raijuOutput(raiApp%Model,raiApp%Grid,raiApp%State)
        endif

        if (raiApp%State%IO%doConsole(raiApp%State%t)) then
            call raiApp%WriteConsoleOutput()
        endif
        call Toc("Output")

    ! Update potentials and magnetic field somehow if you want
    ! This used to be where one way driving voltron -> RAIJU happened

    ! Step RAIJU
        call Tic("RAIJU Advance")
        call raiApp%AdvanceModel(raiApp%Model%dt)
        call Toc("RAIJU Advance")

        !write(*,*)raiApp%State%t
        !write(*,*)"MJDs: "
        !write(*,*)raiApp%State%mjd

        ! Advance model times
        raiApp%State%t  = raiApp%State%t  + raiApp%Model%dt
        raiApp%State%ts = raiApp%State%ts + 1
        raiApp%State%mjd = T2MJD(raiApp%State%t,mjd0)

        ! Update geopack if we need to
        if (raiApp%Model%doGeoCorot) then
            call mjdRECALC(raiApp%State%mjd)
        endif

        if (raiApp%Model%doClockConsoleOut) then
            call printClocks()
            call cleanClocks()
        endif
    enddo


    call Toc("Omega")
    write(*,*)"Main Objective complete!"

    contains


end program raijuSAx