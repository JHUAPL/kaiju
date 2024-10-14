! Collection of data and objects for the voltron middle man
! MPI version

module voltapp_mpi
    use volttypes_mpi
    use gamcouple_mpi_V2G
    use couplingHelpers
    use mpi_f08
    use ebsquish, only : SquishBlocksRemain, DoSquishBlock
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN, IEEE_QUIET_NAN
    use volthelpers_mpi
    use voltapp
    use gcm_mpi
    
    implicit none

    contains

    !end lingering MPI waits, if there are any
    subroutine endVoltronWaits(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        logical :: reqStat
        integer :: ierr

        call vApp%gApp%Cleanup()

        if(vApp%vHelpWin /= MPI_WIN_NULL) then
            call MPI_WIN_FREE(vApp%vHelpWin, ierr)
        endif

    end subroutine endVoltronWaits

    !Initialize Voltron (after Gamera has already been initialized)
    subroutine initVoltron_mpi(vApp, optFilename)
        type(voltAppMpi_T), intent(inout) :: vApp
        character(len=*), optional, intent(in) :: optFilename

        type(MPI_Comm) :: allVoltComm
        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        integer :: commSize, ierr, numCells, length, ic
        integer, allocatable, dimension(:) :: neighborRanks, inData, outData
        integer :: nHelpers, gamNRES, commId
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder, wasWeighted
        integer(KIND=MPI_BASE_MYADDR) :: winsize

        if(.not. allocated(vApp%gApp)) then
            ! mpi voltron uses mpi remotely coupled gamera
            allocate(gamCouplerMpi_volt_T :: vApp%gApp)
            ! set varible for polymorphic type
            SELECT type(cplApp=>vApp%gApp)
                TYPE IS (gamCouplerMpi_volt_T)
                    cplApp%gOptionsCplMpiV%couplingPoolComm => vApp%vOptionsMpi%couplingPoolComm
                CLASS DEFAULT
                    write(*,*) 'Gamera MPI coupler is wrong type'
                    stop
            END SELECT
        endif

        ! read helper XML options
        if(present(optFilename)) then
            inpXML = optFilename
        else
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck in initVoltron_mpi, exiting ...")
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Gamera',.true.)

        ! initialize F08 MPI objects
        vApp%vHelpComm = MPI_COMM_NULL
        vApp%vHelpWin = MPI_WIN_NULL

        call voltronSplitWithApp(vApp%vOptionsMpi%couplingPoolComm, helperId, 0, allVoltComm)

        ! get info about voltron-only mpi communicator
        call MPI_Comm_rank(allVoltComm, vApp%vHelpRank, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call xmlInp%Set_Val(vApp%useHelpers,"/Kaiju/Voltron/Helpers/useHelpers",.false.)
        call xmlInp%Set_Val(vApp%doSquishHelp,"/Kaiju/Voltron/Helpers/doSquishHelp",.true.)
        call xmlInp%Set_Val(vApp%squishLoadBalance,"/Kaiju/Voltron/Helpers/squishLoadBalance",.true.)
        call xmlInp%Set_Val(nHelpers,"/Kaiju/Voltron/Helpers/numHelpers",0)

        if(vApp%masterSquish .and. vApp%squishLoadBalance) then
            print *,"Dynamic load balancing of squish helpers is not supported if the"
            print *,"     voltron master rank is also calculating squish."
            print *,"     Please set either:"
            print *,"        /Kaiju/Voltron/Helpers/masterSquish       to False"
            print *,"     or"
            print *,"        /Kaiju/Voltron/Helpers/squishLoadBalance  to False"
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        call MPI_Comm_Size(allVoltComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.not. vApp%useHelpers) nHelpers = 0
        if(nHelpers .eq. 0) vApp%useHelpers = .false.
        if(nHelpers .ne. commSize-1) then
            print *,"The number of voltron helpers is not correct."
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        ! additional MPI options
        call xmlInp%Set_Val(vApp%doSerialMHD,"coupling/doSerial",.false.)

        ! now initialize basic voltron structures from gamera data
        if(vApp%gcmCplRank /= -1) then
          call init_gcm_mix_mpi(vApp%gcm,vApp%mageCplComm,vApp%gcmCplRank)
        endif
        if(present(optFilename)) then
            call initVoltron(vApp, optFilename)
        else
            call initVoltron(vApp)
        endif

        if ( ( vApp%doGCM ) .and. (vApp%gcmCplRank /= -1) ) then
            write(*,*) "Initializing GCM ..."
            call init_gcm_mpi(vApp%gcm,vApp%remixApp%ion,vApp%gApp%Model%isRestart)
        end if

        if(vApp%useHelpers) then
            if(vApp%doSquishHelp) then
                ! over-ride the number of squish blocks
                if(vApp%masterSquish) then
                    vApp%ebTrcApp%ebSquish%numSquishBlocks = nHelpers + 1
                else
                    vApp%ebTrcApp%ebSquish%numSquishBlocks = nHelpers
                endif

                if(vApp%ebTrcApp%ebSquish%numSquishBlocks /= size(vApp%ebTrcApp%ebSquish%blockStartIndices)) then
                    ! number of squish blocks changed
                    deallocate(vApp%ebTrcApp%ebSquish%blockStartIndices)
                    allocate(vApp%ebTrcApp%ebSquish%blockStartIndices(vApp%ebTrcApp%ebSquish%numSquishBlocks))
                    call LoadBalanceBlocks(vApp) ! redo the default load balancing for this number of blocks
                endif

            endif

            call MPI_Comm_Size(allVoltComm, commSize, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            ! add topology to the helper communicator to permit neighborhood operations
            reorder = .false. ! don't allow MPI to reorder the ranks, master must remain master
            allocate(neighborRanks(1:commSize-1))
            allocate(inData(1:commSize-1))
            allocate(outData(1:commSize-1))
            do ic=1,commSize-1
                neighborRanks(ic) = ic
            enddo
            ! currently all helpers only send and receive data from rank 0
            inData(:) = 1
            outData(:) = 1
            call mpi_dist_graph_create_adjacent(allVoltComm, &
                commSize-1,neighborRanks,inData, &
                commSize-1,neighborRanks,outData, &
                MPI_INFO_NULL, reorder, vApp%vHelpComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            deallocate(neighborRanks, inData, outData)

            ! send initial timing data to helper voltron ranks
            call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
            call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)

            ! allocate data for helpers to set status, and share with them
            allocate(vApp%vHelpIdle(nHelpers))
            vApp%vHelpIdle = 0
            call mpi_type_extent(MPI_INTEGER, length, ierr) ! size of integer
            winsize = nHelpers*length ! have to use this because it's an odd size integer
            call mpi_win_create(vApp%vHelpIdle, winsize, length, MPI_INFO_NULL, vApp%vHelpComm, vApp%vHelpWin, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
        endif

        if(.not. vApp%doSerialMHD) then
            !initial coupling
            call Tic("DeepUpdate")
            call DeepUpdate_mpi(vApp)
            call Toc("DeepUpdate")
        endif
        
    end subroutine initVoltron_mpi

    subroutine initVoltronHelper_mpi(vApp, optFilename)
        type(voltAppMpi_T), intent(inout) :: vApp
        character(len=*), optional, intent(in) :: optFilename

        type(MPI_Comm) :: allVoltComm
        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        integer :: commSize, ierr, numCells, length, ic
        integer :: nHelpers, gamNRES, commId
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder, wasWeighted
        integer(KIND=MPI_BASE_MYADDR) :: winsize

        ! replace the gamera option with one customized for squish helpers
        if(allocated(vApp%gApp)) deallocate(vApp%gApp)
        allocate(SHgamCoupler_T :: vApp%gApp)

        ! read helper XML options
        if(present(optFilename)) then
            inpXML = optFilename
        else
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck in initVoltron_mpi, exiting ...")
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Gamera',.true.)

        ! initialize F08 MPI objects
        vApp%vHelpComm = MPI_COMM_NULL
        vApp%vHelpWin = MPI_WIN_NULL

        call appWaitForVoltronSplit(vApp%vOptionsMpi%couplingPoolComm, helperId, 0, allVoltComm)

        ! get info about voltron-only mpi communicator
        call MPI_Comm_rank(allVoltComm, vApp%vHelpRank, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        vApp%isLoud = .false.
        vApp%writeFiles = .false.
        vApp%gApp%Model%isLoud = .false.

        if (.not. vApp%isLoud) call xmlInp%BeQuiet()

        call xmlInp%Set_Val(vApp%useHelpers,"/Kaiju/Voltron/Helpers/useHelpers",.false.)
        call xmlInp%Set_Val(vApp%doSquishHelp,"/Kaiju/Voltron/Helpers/doSquishHelp",.true.)
        call xmlInp%Set_Val(vApp%masterSquish,"/Kaiju/Voltron/Helpers/masterSquish",.false.)
        call xmlInp%Set_Val(vApp%squishLoadBalance,"/Kaiju/Voltron/Helpers/squishLoadBalance",.true.)
        call xmlInp%Set_Val(nHelpers,"/Kaiju/Voltron/Helpers/numHelpers",0)
        call MPI_Comm_Size(allVoltComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.not. vApp%useHelpers) then
            print *,"Voltron helpers were created, but the helping option is disabled."
            print *,"Please either turn the /Voltron/Helpers/useHelpers option on, or "
            print *,"  remove the unnecessary Voltron helper ranks."
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif
        if(nHelpers .ne. commSize-1) then
            print *,"The number of voltron helpers is not correct."
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        ! add topology to the helper communicator to permit neighborhood operations
        reorder = .false. ! don't allow MPI to reorder the ranks, master must remain master
        call mpi_dist_graph_create_adjacent(allVoltComm, &
            1,(/0/),(/1/), &
            1,(/0/),(/1/), &
            MPI_INFO_NULL, reorder, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! now initialize basic voltron structures from gamera data
        if(present(optFilename)) then
            call initVoltron(vApp, optFilename)
        else
            call initVoltron(vApp)
        endif

        ! get starting time values from master voltron rank
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)

        ! accept the mpi window on the voltron master rank, expose none of my own memory
        winsize = 0 ! have to use this because it's an odd size integer
        call mpi_win_create(MPI_BOTTOM, winsize, 1, MPI_INFO_NULL, vApp%vHelpComm, vApp%vHelpWin, ierr)
                    if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    !Step Voltron if necessary (currently just updating state variables)
    subroutine stepVoltron_mpi(vApp, dt)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: dt

        real(rp) :: stepEndTime

        stepEndTime = vApp%time + dt

        do while(stepEndTime .ge. vApp%DeepT)
            ! if Gamera is processing data, finish it
            ! if this is the first call, this will noop
            call vApp%gApp%FinishUpdateMhdData(vApp)

            vApp%time = vApp%DeepT
            vApp%MJD = T2MJD(vApp%time,vApp%gApp%Model%MJD0)

            ! update the next predicted coupling interval
            vApp%DeepT = vApp%DeepT + vApp%DeepDT

            if(.not. vApp%doSerialMHD) call vApp%gApp%StartUpdateMhdData(vApp)

            call Tic("DeepUpdate")
            call DeepUpdate_mpi(vApp)
            call Toc("DeepUpdate")

            if(vApp%doSerialMHD) call vApp%gApp%StartUpdateMhdData(vApp)

        enddo

        ! step end time is greater than, or equal to, the current DeepT
        ! advance to that partial deep step time
        vApp%time = stepEndTime
        vApp%MJD = T2MJD(vApp%time,vApp%gApp%Model%MJD0)

    end subroutine stepVoltron_mpi

!----------
!Deep coupling stuff (time coming from vApp%time, so in seconds)

    subroutine startDeep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        ! convert gamera data to mixInput
        call Tic("G2R")
        call convertGameraToRemix(vApp%mhd2mix, vApp%gApp, vApp%remixApp)
        call Toc("G2R")

        if (vApp%doGCM .and. vApp%time >=0 .and. vApp%gcmCplRank /= -1) then
            call Tic("GCM2MIX")
            call coupleGCM2MIX(vApp%gcm,vApp%remixApp%ion,vApp%MJD,vApp%time,vApp%mageCplComm,vApp%gcmCplRank)
            call Toc("GCM2MIX")
        end if

        ! run remix
        call Tic("ReMIX", .true.)
        call runRemix(vApp)
        call Toc("ReMIX", .true.)

        call Tic("R2G")
        call CouplePotentialToMhd(vApp)
        call Toc("R2G")

        ! only do imag after spinup
        if(vApp%doDeep .and. vApp%time >= 0) then
            call Tic("DeepUpdate", .true.)
            call PreDeep(vApp, vApp%gApp)
            call SquishStart(vApp)

            if(vApp%useHelpers .and. vApp%doSquishHelp) then
                call Tic("VoltHelpers", .true.)
                call vhReqStep(vApp)
                call vhReqSquishStart(vApp)
                call Toc("VoltHelpers", .true.)
            endif

            ! moving this to after voltron helpers are started
            call DoImag(vApp)

            vApp%deepProcessingInProgress = .true.
            call Toc("DeepUpdate", .true.)
        else
            vApp%gApp%Grid%Gas0 = 0
        endif

    end subroutine startDeep

    subroutine endDeep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        vApp%deepProcessingInProgress = .false.

        ! only do imag after spinup with deep enabled
        if(vApp%doDeep .and. vApp%time >= 0) then
            call Tic("DeepUpdate", .true.)

            do while(SquishBlocksRemain(vApp))
                call Tic("Squish",.true.)
                call DoSquishBlock(vApp)
                call Toc("Squish",.true.)
            enddo

            vApp%deepProcessingInProgress = .false.

            if(vApp%useHelpers .and. vApp%doSquishHelp) then
                call Tic("VoltHelpers", .true.)
                call vhReqSquishEnd(vApp)
                call Toc("VoltHelpers", .true.)
            endif

            call SquishEnd(vApp)
            call PostDeep(vApp, vApp%gApp)
            call Toc("DeepUpdate", .true.)
        endif

    end subroutine endDeep

    function deepInProgress(vApp)
        type(voltAppMpi_T), intent(in) :: vApp
        logical :: deepInProgress

        deepInProgress = vApp%deepProcessingInProgress

    end function deepInProgress

    subroutine doDeepBlock(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        if(.not. vApp%deepProcessingInProgress) return

        if(SquishBlocksRemain(vApp)) then
            call Tic("DeepUpdate")
            call Tic("Squish",.true.)
            call DoSquishBlock(vApp)
            call Toc("Squish",.true.)
            call Toc("DeepUpdate")
        endif

        if(.not. SquishBlocksRemain(vApp)) then
            if((.not. vApp%useHelpers) .or. &
               (vApp%useHelpers .and. .not. vApp%doSquishHelp) .or. &
               (vApp%useHelpers .and. vApp%doSquishHelp .and. allHelpersIdle(vApp)) ) then
                call endDeep(vApp)
            endif
        endif

    end subroutine doDeepBlock

    subroutine DeepUpdate_mpi(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        real(rp) :: tAdv
        integer :: ierr

        ! ensure deep update is complete
        do while(deepInProgress(vApp))
            call doDeepBlock(vApp)
        enddo

        ! setup squish operation but don't yet perform the computations
        call startDeep(vApp)
        call endDeep(vApp)

    end subroutine

    subroutine helpVoltron(vApp, helperQuit)
        type(voltAppMpi_T), intent(inout) :: vApp
        logical, intent(out) :: helperQuit ! should the helper quit

        integer :: ierr, helpType
        type(MPI_Request) :: helpReq

        helperQuit = .false. ! don't quit normally

        ! assumed to only be in this function if helpers are enabled

        ! async wait for command type
        call mpi_Ibcast(helpType, 1, MPI_INTEGER, 0, vApp%vHelpComm, helpReq, ierr)
        call mpi_wait(helpReq, MPI_STATUS_IGNORE, ierr)

        !write (*,*) 'Helper got request type: ', helpType

        ! not idle now
        call setIdleStatus(vApp, helpType)

        ! then call appropriate function to deal with command
        select case(helpType)
            CASE (VHSTEP)
                call vhHandleStep(vApp)
            CASE (VHSQUISHSTART)
                call vhHandleSquishStart(vApp)
            CASE (VHSQUISHEND)
                call vhHandleSquishEnd(vApp)
            case (VHQUIT)
                helperQuit = .true.
            CASE DEFAULT
                print *,"Unknown voltron helper request type received."
                stop
        end select

        ! idle now
        call setIdleStatus(vApp, 0)

    end subroutine

end module voltapp_mpi

