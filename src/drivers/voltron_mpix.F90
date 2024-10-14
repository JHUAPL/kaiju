!Driver for mpi decomposed Gamera coupled with Voltron (remix only)

program voltron_mpix
    use clocks
    use gamCouple_mpi_G2V
    use voltapp_mpi
    use couplingHelpers
    use output
    use voltio
    use uservoltic
    use mpi_f08
    use xml_input

    implicit none

    ! create allocatable objects for the possible app types
    ! only one will be used by each executable
    type(voltAppMpi_T), allocatable :: vApp
    type(gamCouplerMpi_gam_T), allocatable :: gApp

    procedure(StateIC_T), pointer :: userInitFunc => initUser
    integer :: ierror, length, provided, worldSize, worldRank, numHelpers, appId
    type(MPI_Comm) :: appComm, splittingComm
    integer :: required=MPI_THREAD_MULTIPLE
    character( len = MPI_MAX_ERROR_STRING) :: message
    character(len=strLen) :: inpXML, helpersBuf
    logical :: useHelpers, helperQuit
    integer(KIND=MPI_AN_MYADDR) :: tagMax
    logical :: tagSet
    type(XML_Input_T) :: xmlInp
    real(rp) :: nextDT
    integer :: divideSize,i

    ! initialize MPI
    !Set up MPI with or without thread support
#ifdef _OPENMP
    call MPI_INIT_THREAD(required,provided,ierror)
    if(provided < required) then
        print *,"Not support for MPI_THREAD_MULTIPLE, aborting!"
        call abort()
    end if
    !print *,"MPI + OpenMP !!!"
#else
    !print *," MPI without threading"
    call MPI_INIT(ierror)
#endif

    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    call MPI_Comm_Size(MPI_COMM_WORLD, worldSize, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    call MPI_Comm_Rank(MPI_COMM_WORLD, worldRank, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    ! initialize mpi data type
    call setMpiReal()

    call initClocks()

    call mpi_comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tagMax, tagSet, ierror)
    !print *, 'Tag Upper-Bound = ', tagMax

    call getIDeckStr(inpXML)
    call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

    !Create XML reader
    write(*,*) 'Reading input deck from ', trim(inpXML)
    xmlInp = New_XML_Input(trim(inpXML),'Kaiju',.true.)

    ! need to know how many voltron helpers there are
    call getIDeckStr(inpXML)
    call ReadXmlImmediate(trim(inpXML),'/Kaiju/Voltron/Helpers/useHelpers',helpersBuf,'F',.false.)
    read(helpersBuf,*) useHelpers
    call ReadXmlImmediate(trim(inpXML),'/Kaiju/Voltron/Helpers/numHelpers',helpersBuf,'0',.false.)
    read(helpersBuf,*) numHelpers
    if(.not. useHelpers) numHelpers = 0

    if(worldRank .gt. (worldSize-1-numHelpers)) then
        appId = helperId
    elseif(worldRank .eq. worldSize-1-numHelpers) then
        appId = voltId
    else
        appId = gamId
    endif
    ! every app splits into their own comunicator
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, appId, 0, appComm, ierror)

    if(appId == helperId) then
        allocate(vApp)
        vApp%vOptions%gamUserInitFunc => initUser
        allocate(vApp%vOptionsMpi%couplingPoolComm)
        vApp%vOptionsMpi%couplingPoolComm = MPI_COMM_WORLD
        call initVoltronHelper_mpi(vApp)

        ! do helper loop
        helperQuit = .false.
        do while(.not. helperQuit)
            call Tic("Omega", .true.)
            call helpVoltron(vApp, helperQuit)
            call Toc("Omega", .true.)
        end do
        call endVoltronWaits(vApp)
    elseif(appId == voltId) then

        allocate(vApp)

        splittingComm = MPI_COMM_WORLD

        ! *** CALL TGCM comms coupling first
        ! this should all be moved out of the driver and into TGCM coupling code later on
        call ReadXmlImmediate(trim(inpXML),'/Kaiju/Voltron/coupling/doGCM',helpersBuf,'F',.false.)
        read(helpersBuf,*) vApp%doGCM
        call voltronSplitWithApp(splittingComm, mageId, 0, vApp%mageCplComm)
        call MPI_Comm_Rank(vApp%mageCplComm, vApp%voltCplRank, ierror)
        call MPI_Comm_Size(vApp%mageCplComm, vApp%CplSize, ierror)
        if(vApp%mageCplComm /= MPI_COMM_NULL) then
            print *,'VOLTRON has created an mageCplComm with ', vApp%CplSize-1, ' other app(s)'
            print *,'VOLTRON using mageCplComm tag ', mageId
            write(*,*) "VOLTRON CPLCOMM: ",vApp%mageCplComm,vApp%voltCplRank

            ! Tell everyone who I am
            if (.not.allocated(vApp%IAm)) allocate(vApp%IAm(vApp%CplSize))
            vApp%IAm(vApp%voltCplRank+1) = voltId
            do i=1,vApp%CplSize
                call MPI_Bcast(vApp%IAm(i),1,MPI_INTEGER,i-1,vApp%mageCplComm,ierror)
            enddo

            ! IAm array starts at 1, MPI Ranks start at 0
            do i=1,vApp%CplSize
                ! Assign rank if match
                select case (vApp%IAm(i))
                    case (voltId)
                        if (i-1 .ne. vApp%voltCplRank) then
                            write(*,*) "I AM NOT MYSELF:",i-1,vApp%voltCplRank
                        endif
                    case (gamId)
                        write(*,*) "Gam not involved in MAGE2MAGE yet"
                    case (rcmId)
                        write(*,*) "RCM not involved in MAGE2MAGE yet"
                    case (hidraNId)
                        vApp%hidraNCplRank = i-1
                        write(*,*) "Volt coupling to hidraN"
                    case (hidraSId)
                        vApp%hidraSCplRank = i-1
                        write(*,*) "Volt coupling to hidraS"
                    case (hidraId)
                        vApp%hidraCplRank = i-1
                        write(*,*) "Volt coupling to hidra"
                     case (tgcmId)
                        vApp%gcmCplRank = i-1
                        write(*,*) "Volt coupling to TIEGCM"
                    case default
                        write(*,*) "Volt does not know about this Coupling ID: ", vApp%IAm(i)
                end select
            enddo
        endif
        if(vApp%CplSize == 1 .and. (.not. vApp%doGCM)) then
            write(*,*) "VOLTRON: We're not coupling to a GCM"
            call mpi_comm_free(vApp%mageCplComm, ierror)
            if(ierror /= MPI_Success) then
                call MPI_Error_string( ierror, message, length, ierror)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
            end if
            vApp%mageCplComm = MPI_COMM_NULL
        endif
        if (vApp%CplSize == 1 .and. vApp%doGCM) then
            write(*,*) "VOLTRON: Coupling to GCM Failed."
            call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
        endif
        ! *** end TIEGCM coupling code

        vApp%vOptions%gamUserInitFunc => initUser
        allocate(vApp%vOptionsMpi%couplingPoolComm)
        vApp%vOptionsMpi%couplingPoolComm = splittingComm
        call initVoltron_mpi(vApp)

        ! voltron run loop
        do while (vApp%time < vApp%tFin)
            !Start root timer
            call Tic("Omega", .true.)

            !Advance Voltron models one coupling step
            nextDT = min(vApp%tFin-vApp%time, vApp%IO%nextIOTime()-vApp%time)
            call Tic("StepVoltron")
            call stepVoltron_mpi(vApp, nextDT)
            call Toc("StepVoltron")

            !IO checks
            call Tic("IO", .true.)
            !Console output
            if (vApp%IO%doConsole(vApp%time)) then
                call consoleOutputV(vApp,vApp%gApp)
                if (vApp%IO%doTimerOut) then
                    call printClocks()
                endif
                call cleanClocks()
            elseif (vApp%IO%doTimer(vApp%time)) then
                if (vApp%IO%doTimerOut) then
                    call printClocks()
                endif
                call cleanClocks()
            endif

            !Restart output
            if (vApp%IO%doRestart(vApp%time)) then
                call resOutputV(vApp,vApp%gApp)
            endif
            !Data output
            if (vApp%IO%doOutput(vApp%time)) then
                call fOutputV(vApp,vApp%gApp)
            endif

            call Toc("IO", .true.)
            call Toc("Omega", .true.)
        enddo
        ! if using helpers, tell them to quit
        if(vApp%useHelpers) call vhReqHelperQuit(vApp)
        call endVoltronWaits(vApp)
    elseif(appId == gamId) then
        ! gamera rank
        allocate(gApp)
        gApp%gOptionsMpi%gamComm = appComm

        gApp%gOptionsCplMpiG%couplingPoolComm = MPI_COMM_WORLD
        gApp%gOptionsMpi%doIO = .false.
        gApp%gOptions%userInitFunc => initUser

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
    else
        write (*,*) "Unrecognized appId in voltron_mpi.x"
    endif

    call MPI_FINALIZE(ierror)
    write(*,*) "Fin Voltron Mpi"
    

end program voltron_mpix

