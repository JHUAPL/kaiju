!Driver for mpi decomposed Gamera coupled with Voltron (remix only)

program voltron_mpix
    use clocks
    use gamCouple_mpi_G2V
    use coupled_layout_mpi
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
    type(MPI_Comm) :: appComm, splittingComm, world_dup
    integer :: required=MPI_THREAD_MULTIPLE
    character( len = MPI_MAX_ERROR_STRING) :: message
    character(len=strLen) :: inpXML, helpersBuf
    logical :: useHelpers, helperQuit
    integer(KIND=MPI_AN_MYADDR) :: tagMax
    logical :: tagSet
    type(XML_Input_T) :: xmlInp
    type(coupledLayout_T) :: mpiLayout
    real(rp) :: nextDT
    integer :: divideSize,i
    logical :: doResetClocks = .false.

    ! initialize MPI
    ! This is in case we want to do duplicate communicator from the start
    world_dup = MPI_COMM_WORLD

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
        call mpi_Abort(world_dup, 1, ierror)
    end if

    call MPI_Comm_Size(world_dup, worldSize, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(world_dup, 1, ierror)
    end if
    call MPI_Comm_Rank(world_dup, worldRank, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(world_dup, 1, ierror)
    end if

    ! initialize mpi data type
    call setMpiReal()

    call initClocks()

    call mpi_comm_get_attr(world_dup, MPI_TAG_UB, tagMax, tagSet, ierror)
    !print *, 'Tag Upper-Bound = ', tagMax

    call getIDeckStr(inpXML)
    call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

    !Create XML reader
    write(*,*) 'Reading input deck from ', trim(inpXML)
    xmlInp = New_XML_Input(trim(inpXML),'Kaiju',.true.)

    ! are helpers being used
    call ReadXmlImmediate(trim(inpXML),'/Kaiju/Voltron/Helpers/useHelpers',helpersBuf,'F',.false.)
    read(helpersBuf,*) useHelpers

    call initVoltronCoupledLayout(xmlInp, worldSize, mpiLayout)

    if(mpiLayout%hasCustomLayout) then
        call getLayoutAppId(mpiLayout, worldRank, appId)
    else
        ! need to know how many voltron helpers there are
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
    endif
    ! every app splits into their own comunicator
    call MPI_COMM_SPLIT(world_dup, appId, worldRank, appComm, ierror)

    if(appId == helperId) then
        allocate(vApp)
        vApp%vOptions%gamUserInitFunc => initUser
        allocate(vApp%vOptionsMpi%couplingPoolComm)
        vApp%vOptionsMpi%couplingPoolComm = world_dup
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

        splittingComm = world_dup

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

            ! === TIEGCM ===
            call getCplRanks(vApp%IAm, vApp%cplSize, tgcmId, vApp%gcmCplRank, vApp%gcmNumRanks, vApp%gcmRankList)
            if (vApp%gcmNumRanks .ne. 0) then
                write(*,*) "V coupling to ", vApp%gcmNumRanks," ranks of T"
            endif

            ! === WACCM-X ===
            ! Check for waccm-x coupling only if we're not coupling tiegcm
            if (vApp%gcmNumRanks == 0) then
                call getCplRanks(vApp%IAm, vApp%cplSize, wcmxId, vApp%gcmCplRank, vApp%gcmNumRanks, vApp%gcmRankList)
                if (vApp%gcmNumRanks .ne. 0) then
                    write(*,*) "V coupling to ", vApp%gcmNumRanks," ranks of W"
                    vApp%gcm%isWACCMX = .true.
                    vApp%gcm%nWACCMX = vApp%gcmNumRanks
                endif
            endif
            
            ! === HIDRA-N ===
            call getCplRanks(vApp%IAm, vApp%cplSize, hidraNId, vApp%hidraNCplRank, vApp%hidraNNumRanks, vApp%hidraNRankList)
            if (vApp%hidraNNumRanks .ne. 0) then
                write(*,*) "V coupling to ", vApp%hidraNNumRanks," ranks of H-N"
            endif
            
            ! === HIDRA-S ===
            call getCplRanks(vApp%IAm, vApp%cplSize, hidraSId, vApp%hidraSCplRank, vApp%hidraSNumRanks, vApp%hidraSRankList)
            if (vApp%hidraSNumRanks .ne. 0) then
                write(*,*) "V coupling to ", vApp%hidraSNumRanks," ranks of H-S"
            endif
        endif
        if(vApp%CplSize == 1 ) then
            write(*,*) "VOLTRON: We're not coupling to a GCM"
            call mpi_comm_free(vApp%mageCplComm, ierror)
            if(ierror /= MPI_Success) then
                call MPI_Error_string( ierror, message, length, ierror)
                print *,message(1:length)
                call mpi_Abort(world_dup, 1, ierror)
            end if
            vApp%mageCplComm = MPI_COMM_NULL
        endif
        if (vApp%CplSize == 1 .and. vApp%doGCM) then
            write(*,*) "VOLTRON: Coupling to GCM Failed."
            call mpi_Abort(world_dup, 1, ierror)
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
                doResetClocks = .true.
            elseif (vApp%IO%doTimer(vApp%time)) then
                if (vApp%IO%doTimerOut) then
                    call printClocks()
                endif
                doResetClocks = .true.
            endif

            !Data output
            if (vApp%IO%doOutput(vApp%time)) then
                call fOutputV(vApp,vApp%gApp)
            endif

            !Restart output
            if (vApp%IO%doRestart(vApp%time)) then
                call resOutputV(vApp,vApp%gApp)
            endif
            
            !Reset clocks last so data is available for all output
            if (doResetClocks) then
                call cleanClocks()
                doResetClocks = .false.
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

        gApp%gOptionsCplMpiG%couplingPoolComm = world_dup
        gApp%gOptionsMpi%doIO = .false.
        gApp%gOptions%userInitFunc => initUser

        call gApp%InitModel(xmlInp)
        call gApp%InitIO(xmlInp)

        do while (gApp%Model%t < gApp%Model%tFin)
            call Tic("Omega", .true.) !Start root timer

            !Step model
            nextDT = min(gApp%Model%tFin-gApp%Model%t, gApp%Model%IO%nextIOTime()-gApp%Model%t)
            call gApp%AdvanceModel(nextDT)

            !Output if necessary
            call Tic("IO")

            if (gApp%Model%IO%doConsole(gApp%Model%t)) then
                call gApp%WriteConsoleOutput()

                !Timing info
                if (gApp%Model%IO%doTimerOut) call printClocks()
                doResetClocks = .true.

            elseif (gApp%Model%IO%doTimer(gApp%Model%t)) then
                if (gApp%Model%IO%doTimerOut) call printClocks()
                doResetClocks = .true.
            endif

            if (gApp%Model%IO%doOutput(gApp%Model%t)) then
                do while(gApp%Model%IO%doOutput(gApp%Model%t))
                    call gApp%WriteFileOutput(gApp%Model%IO%nOut)
                end do
            endif

            if (gApp%Model%IO%doRestart(gApp%Model%t)) then
                call gApp%WriteRestart(gApp%Model%IO%nRes)
            endif

            if (doResetClocks) then
                call cleanClocks()
                doResetClocks = .false.
            endif

            call Toc("IO")
            call Toc("Omega", .true.)
        end do
    else
        write (*,*) "Unrecognized appId in voltron_mpi.x"
    endif

    call MPI_FINALIZE(ierror)
    write(*,*) "Fin Voltron Mpi"
    

end program voltron_mpix
