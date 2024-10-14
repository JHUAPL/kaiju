!Driver for Gamera with MPI decomposition

program gamera_mpix
    use gamapp_mpi
    use usergamic
    use mpidefs

    implicit none

    type(gamAppMpi_T) :: gAppMpi

    integer :: ierror, length
    integer :: required=MPI_THREAD_MULTIPLE
    integer :: provided
    character( len = MPI_MAX_ERROR_STRING) :: message
    character(len=strLen) :: inpXML
    type(XML_Input_T) :: xmlInp
    real(rp) :: nextDT

    ! initialize MPI
    !Set up MPI with or without thread support
#ifdef _OPENMP
    call MPI_INIT_THREAD(required,provided,ierror)
    if(provided < required) then
        print *,"Not support for MPI_THREAD_MULTIPLE, aborting!"
        call abort()
    end if
    print *,"MPI + OpenMP !!!"
#else
    print *," MPI without threading"
    call MPI_INIT(ierror)
#endif

    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    ! initialize mpi data type
    call setMpiReal()

    ! set options for gamera
    gAppMpi%gOptions%userInitFunc => initUser
    gAppMpi%gOptionsMpi%gamComm = MPI_COMM_WORLD

    !call printConfigStamp()
    call initClocks()

    call getIDeckStr(inpXML)
    call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

    !Create XML reader
    write(*,*) 'Reading input deck from ', trim(inpXML)
    xmlInp = New_XML_Input(trim(inpXML),'Kaiju',.true.)

    call gAppMpi%InitModel(xmlInp)
    call gAppMpi%InitIO(xmlInp)

    do while (gAppMpi%Model%t < gAppMpi%Model%tFin)
        call Tic("Omega") !Start root timer

        !Step model
        nextDT = min(gAppMpi%Model%tFin-gAppMpi%Model%t, gAppMpi%Model%IO%nextIOTime()-gAppMpi%Model%t)
        call gAppMpi%AdvanceModel(nextDT)

        !Output if necessary
        call Tic("IO")

        if (gAppMpi%Model%IO%doConsole(gAppMpi%Model%t)) then
            call gAppMpi%WriteConsoleOutput()

            !Timing info
            if (gAppMpi%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()

        elseif (gAppMpi%Model%IO%doTimer(gAppMpi%Model%t)) then
            if (gAppMpi%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif

        if (gAppMpi%Model%IO%doOutput(gAppMpi%Model%t)) then
            call gAppMpi%WriteFileOutput(gAppMpi%Model%IO%nOut)
        endif

        if (gAppMpi%Model%IO%doRestart(gAppMpi%Model%t)) then
            call gAppMpi%WriteRestart(gAppMpi%Model%IO%nRes)
        endif

        call Toc("IO")
        call Toc("Omega")
    end do

    call MPI_FINALIZE(ierror)
    write(*,*) "Fin Mpi"

end program gamera_mpix

