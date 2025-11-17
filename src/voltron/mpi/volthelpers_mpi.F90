! Collection of data and objects for additional voltron helper ranks

module volthelpers_mpi
    use volttypes_mpi
    use clocks
    use mpi_f08
    use loadBalance
    use ebsquish, only : SquishStart, GetSquishBds, SquishBlocksRemain, DoSquishBlock, GetAdjustedSquishStart
    use voltCplHelper, only : calcTubes
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN, IEEE_QUIET_NAN

    implicit none

    enum, bind(C)
        enumerator :: VHSTEP=1,VHQUIT,VHSQUISHSTART,VHSQUISHEND,VHTUBESTART,VHTUBEEND
    endenum

    type(MPI_Datatype), private :: tubeMpiType

    contains

    ! helpers idle status functions

    function allHelpersIdle(vApp)
        type(voltAppMpi_T), intent(in) :: vApp
        logical :: allHelpersIdle

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        if (.not. vApp%useHelpers) then
            ! shouldn't be here
            allHelpersIdle = .true.
            return
        endif

        call Tic("VoltHelpers", .true.)
        call Tic("VHWait")

        ! lock the data window to read status
        call mpi_win_lock(MPI_LOCK_SHARED, 0, 0, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        allHelpersIdle = all(vApp%vHelpIdle .le. 0)

        ! unlock the data window
        call mpi_win_unlock(0, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call Toc("VHWait")
        call Toc("VoltHelpers", .true.)

    end function

    subroutine setIdleStatus(vApp, newStatus)
        type(voltAppMpi_T), intent(inout) :: vApp
        integer, intent(in) :: newStatus

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message
        integer(KIND=MPI_BASE_MYADDR) :: disp

        ! lock the data window to set my status
        call mpi_win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        disp = vApp%vHelpRank-1 ! my slot in the target's array
        call mpi_put(newStatus, 1, MPI_INTEGER, 0, disp, 1, MPI_INTEGER, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! if desired, could try to use mpi_win_flush_local here to allow this rank to move on

        ! unlock the data window
        call mpi_win_unlock(0, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    ! one-time initialization routine to help set things up
    subroutine voltronAndHelpersInit(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        call createTubeMpiType()

    end subroutine

    ! chimp data update functions

    subroutine sendChimpStateData(ebState, vHelpComm)
        type(ebState_T), intent(inout) :: ebState
        type(MPI_Comm), intent(in) :: vHelpComm

        integer :: ierr, length
        integer :: sendCount
        integer, allocatable, save :: mpiCounts(:), mpiZeroes(:), mpiDispls(:)
        character( len = MPI_MAX_ERROR_STRING) :: message

        if(.not. allocated(mpiCounts)) then
            call mpi_comm_size(vHelpComm, length, ierr)
            allocate(mpiCounts(length-1))
            allocate(mpiZeroes(length-1))
            allocate(mpiDispls(length-1))
            mpiZeroes(:) = 0
            mpiCounts(:) = 0
            mpiDIspls(:) = 0
        endif

        ! Fields
        !  eb1
        sendCount = size(ebState%eb1%dB)
        mpiCounts(:) = sendCount
        call mpi_neighbor_alltoallv(ebState%eb1%dB, mpiCounts, mpiDispls, MPI_MYFLOAT, &
            ebState%eb1%dB, mpiZeroes, mpiDispls, MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sendCount = size(ebState%eb1%E)
        mpiCounts(:) = sendCount
        call mpi_neighbor_alltoallv(ebState%eb1%E, mpiCounts, mpiDispls, MPI_MYFLOAT, &
            ebState%eb1%E, mpiZeroes, mpiDispls, MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sendCount = size(ebState%eb1%W)
        mpiCounts(:) = sendCount
        call mpi_neighbor_alltoallv(ebState%eb1%W, mpiCounts, mpiDispls, MPI_MYFLOAT, &
            ebState%eb1%W, mpiZeroes, mpiDispls, MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb1%time, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        !  eb2
        sendCount = size(ebState%eb2%dB)
        mpiCounts(:) = sendCount
        call mpi_neighbor_alltoallv(ebState%eb2%dB, mpiCounts, mpiDispls, MPI_MYFLOAT, &
            ebState%eb2%dB, mpiZeroes, mpiDispls, MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sendCount = size(ebState%eb2%E)
        mpiCounts(:) = sendCount
        call mpi_neighbor_alltoallv(ebState%eb2%E, mpiCounts, mpiDispls, MPI_MYFLOAT, &
            ebState%eb2%E, mpiZeroes, mpiDispls, MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sendCount = size(ebState%eb2%W)
        mpiCounts(:) = sendCount
        call mpi_neighbor_alltoallv(ebState%eb2%W, mpiCounts, mpiDispls, MPI_MYFLOAT, &
            ebState%eb2%W, mpiZeroes, mpiDispls, MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb2%time, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine sendChimpStateData

    subroutine recvChimpStateData(ebState, vHelpComm)
        type(ebState_T), intent(inout) :: ebState
        type(MPI_Comm), intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message
        integer :: sendCount

        ! Fields
        !  eb1
        sendCount = size(ebState%eb1%dB)
        call mpi_neighbor_alltoallv(ebState%eb1%dB, (/0/), (/0/), MPI_MYFLOAT, &
            ebState%eb1%dB, (/sendCount/), (/0/), MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sendCount = size(ebState%eb1%E)
        call mpi_neighbor_alltoallv(ebState%eb1%E, (/0/), (/0/), MPI_MYFLOAT, &
            ebState%eb1%E, (/sendCount/), (/0/), MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sendCount = size(ebState%eb1%W)
        call mpi_neighbor_alltoallv(ebState%eb1%W, (/0/), (/0/), MPI_MYFLOAT, &
            ebState%eb1%W, (/sendCount/), (/0/), MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb1%time, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        !  eb2
        sendCount = size(ebState%eb2%dB)
        call mpi_neighbor_alltoallv(ebState%eb2%dB, (/0/), (/0/), MPI_MYFLOAT, &
            ebState%eb2%dB, (/sendCount/), (/0/), MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sendCount = size(ebState%eb2%E)
        call mpi_neighbor_alltoallv(ebState%eb2%E, (/0/), (/0/), MPI_MYFLOAT, &
            ebState%eb2%E, (/sendCount/), (/0/), MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sendCount = size(ebState%eb2%W)
        call mpi_neighbor_alltoallv(ebState%eb2%W, (/0/), (/0/), MPI_MYFLOAT, &
            ebState%eb2%W, (/sendCount/), (/0/), MPI_MYFLOAT, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb2%time, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine sendChimpUpdate(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call sendChimpStateData(vApp%ebTrcApp%ebState, vApp%vHelpComm)

        call mpi_bcast(vApp%iDeep, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%qkSquishStride, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%doQkSquish, 1, MPI_LOGICAL, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%rTrc, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%nTrc, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call mpi_bcast(vApp%ebTrcApp%ebSquish%Rinner, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine sendChimpUpdate

    subroutine recvChimpUpdate(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call recvChimpStateData(vApp%ebTrcApp%ebState, vApp%vHelpComm)

        call mpi_bcast(vApp%iDeep, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%qkSquishStride, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%doQkSquish, 1, MPI_LOGICAL, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%rTrc, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%nTrc, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if        
        call mpi_bcast(vApp%ebTrcApp%ebSquish%Rinner, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

! ----------------------------------------------------------

! Functions below for sending and receiving help requests

    subroutine vhRequestType(vApp, rType)
        type(voltAppMpi_T), intent(in) :: vApp
        integer, intent(in) :: rType

        integer :: ierr,wtf
        type(MPI_Request) :: helpReq

        wtf = rType
        ! async to match waiting helper nodes
        call mpi_Ibcast(wtf, 1, MPI_INTEGER, 0, vApp%vHelpComm, helpReq, ierr)
        call mpi_wait(helpReq, MPI_STATUS_IGNORE, ierr)

    end subroutine vhRequestType

    subroutine vhReqStep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call Tic("VHReqStep")
        call vhRequestType(vApp, VHSTEP)
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call Toc("VHReqStep")
    end subroutine

    subroutine vhHandleStep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
    end subroutine

    subroutine vhReqSquishStart(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr,b
        
        call Tic("VHReqSquishS")
        call vhRequestType(vApp, VHSQUISHSTART)

        ! send eb data
        call sendChimpUpdate(vApp)

        ! send which sections of squish data to solve
        if(vApp%masterSquish) then
             ! 1 block per volton rank, including master and helpers
            vApp%ebTrcApp%ebSquish%myNumBlocks = 1
        else
            ! 1 block per helper, none for master
            vApp%ebTrcApp%ebSquish%myNumBlocks = 0
        endif

        ! calculate block start indices from load balanced offsets
        do b=1,vApp%ebTrcApp%ebSquish%numSquishBlocks
            vApp%ebTrcApp%ebSquish%blockStartIndices(b) = &
                GetAdjustedSquishStart(vApp,vApp%squishLb%balStartInd(b))
        enddo

        call mpi_bcast(vApp%ebTrcApp%ebSquish%numSquishBlocks, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%ebTrcApp%ebSquish%blockStartIndices, vApp%ebTrcApp%ebSquish%numSquishBlocks, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        vApp%ebTrcApp%ebSquish%myFirstBlock = 1 ! master is always the first block, if it has one
        vApp%ebTrcApp%ebSquish%curSquishBlock = vApp%ebTrcApp%ebSquish%myFirstBlock
        call Toc("VHReqSquishS")

    end subroutine

    subroutine vhHandleSquishStart(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        integer :: clockStart, clockEnd, clockRate, clockMax

        call recvChimpUpdate(vApp)

        call mpi_bcast(vApp%ebTrcApp%ebSquish%numSquishBlocks, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(vApp%ebTrcApp%ebSquish%numSquishBlocks /= size(vApp%ebTrcApp%ebSquish%blockStartIndices)) then
            ! number of squish blocks changed or was wrong, reallocate
            deallocate(vApp%ebTrcApp%ebSquish%blockStartIndices)
            allocate(vApp%ebTrcApp%ebSquish%blockStartIndices(vApp%ebTrcApp%ebSquish%numSquishBlocks))
        endif
        call mpi_bcast(vApp%ebTrcApp%ebSquish%blockStartIndices, vApp%ebTrcApp%ebSquish%numSquishBlocks, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        vApp%ebTrcApp%ebSquish%myNumBlocks = 1 ! helpers always have 1 block
        if(vApp%masterSquish) then
            ! master voltron rank is helping, rank 1 works on block 2, etc...
            vApp%ebTrcApp%ebSquish%myFirstBlock = vApp%vHelpRank + 1
        else
            ! master voltron rank is not helping, rank 1 works on block 1, etc...
            vApp%ebTrcApp%ebSquish%myFirstBlock = vApp%vHelpRank
        endif
        
        call SquishStart(vApp)

        call system_clock(clockStart,clockRate,clockMax)
        do while(SquishBlocksRemain(vApp))
            call DoSquishBlock(vApp)
        end do
        call system_clock(count=clockEnd)
        if(clockEnd < clockStart) then ! wrap
            vApp%lastSquishTime = real(clockMax-clockStart+clockEnd, rp)/real(clockRate, rp)
        else
            vApp%lastSquishTime = real(clockEnd-clockStart, rp)/real(clockRate, rp)
        endif

    end subroutine

    subroutine vhReqSquishEnd(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr,length,i,firstBlock,ks,ke, oldSizes(4), newSizes(4), offsets(4)
        type(MPI_Datatype) :: newtype
        character( len = MPI_MAX_ERROR_STRING) :: message
        real(rp), dimension(:), allocatable, save :: helperTimes, helperLoads
        real(rp) :: targetTime

        real(rp), parameter :: hAlpha = 0.8 ! rate at which helpers adjust to changing squish loads

        call Tic("VHReqSquishE")
        call vhRequestType(vApp, VHSQUISHEND)

        oldSizes = shape(vApp%chmp2mhd%xyzSquish)
        newSizes = (/vApp%iDeep+2-vApp%ebTrcApp%ebState%ebGr%is, &
                     vApp%ebTrcApp%ebState%ebGr%je+2-vApp%ebTrcApp%ebState%ebGr%js, &
                     0, 2/)
        offsets = (/0, 0, 0, 0/)

        ! collect solved squish data
        if(vApp%masterSquish) then
            firstBlock = 2 ! master rank handles block 0
        else
            firstBlock = 1 ! collect all data
        endif
        do i = firstBlock,vApp%ebTrcApp%ebSquish%numSquishBlocks
            call GetSquishBds(vApp, ks, ke, i)
            offsets(3) = ks-1
            newSizes(3) = ke+1-ks

            call mpi_type_create_subarray(4, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, MPI_MYFLOAT, newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_type_commit(newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_recv(vApp%chmp2mhd%xyzSquish, 1, newtype, i-firstBlock+1, 78111, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_type_free(newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            ! isGood is the same size as xyzsquish, but only 3 dimensions
            call mpi_type_create_subarray(3, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, MPI_LOGICAL, newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_type_commit(newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_recv(vApp%chmp2mhd%isGood, 1, newtype, i-firstBlock+1, 78112, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_type_free(newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if 

        enddo

        if(vApp%squishLoadBalance) then
            call mpi_gather(MPI_IN_PLACE, 0, MPI_MYFLOAT, vAPp%squishLb%instantTimes, &
                    1, MPI_MYFLOAT, 0, vApp%vhelpComm, ierr)
            call updateLoads(vApp%squishLb)
        endif

        call Toc("VHReqSquishE")

    end subroutine

    subroutine vhHandleSquishEnd(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr,length,ks,ke, oldSizes(4), newSizes(4), offsets(4)
        type(MPI_Datatype) :: newtype
        character( len = MPI_MAX_ERROR_STRING) :: message

        oldSizes = shape(vApp%chmp2mhd%xyzSquish)
        newSizes = (/vApp%iDeep+2-vApp%ebTrcApp%ebState%ebGr%is, &
                     vApp%ebTrcApp%ebState%ebGr%je+2-vApp%ebTrcApp%ebState%ebGr%js, &
                     0, 2/)
        offsets = (/0, 0, 0, 0/)

        ! send solved squish data
        call GetSquishBds(vApp, ks, ke, vApp%ebTrcApp%ebSquish%myFirstBlock)
        offsets(3) = ks-1
        newSizes(3) = ke+1-ks

        call mpi_type_create_subarray(4, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, MPI_MYFLOAT, newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_type_commit(newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_send(vApp%chmp2mhd%xyzSquish, 1, newtype, 0, 78111, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_type_free(newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! isGood is the same size as xyzsquish, but only 3 dimensions
        call mpi_type_create_subarray(3, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, MPI_LOGICAL, newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_type_commit(newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_send(vApp%chmp2mhd%isGood, 1, newtype, 0, 78112, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_type_free(newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(vApp%squishLoadBalance) then
            call mpi_gather([vApp%lastSquishTime], 1, MPI_MYFLOAT, [vApp%lastSquishTime], 0, MPI_MYFLOAT, 0, vApp%vhelpComm, ierr)
        endif

    end subroutine

    subroutine vhReqHelperQuit(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call vhRequestType(vApp, VHQUIT)

    end subroutine

    subroutine vhReqTubeStart(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr, r
        call Tic("VHReqTubeS")
        call vhRequestType(vApp, VHTUBESTART)

        ! send eb data
        call sendChimpUpdate(vApp)

        ! load balance by controlling the size of each tube helper's shell grid
        do r=1,vApp%tubeLb%nL
            call mpi_send([vApp%shGrid%js+vApp%tubeLb%balStartInd(r)], 1, MPI_INTEGER, r, 88311, vApp%vHelpComm, ierr)
            call checkAndHandleMpiErrorCode(ierr)
            if (r == vApp%tubeLb%nL) then
                call mpi_send(vApp%shGrid%je, 1, MPI_INTEGER, r, 88312, vApp%vHelpComm, ierr)
            else
                ! minus 2 because the calculation loops to this value +1
                call mpi_send([vApp%shGrid%js+vApp%tubeLb%balStartInd(r+1)-2], 1, MPI_INTEGER, r, 88312, vApp%vHelpComm, ierr)
            endif
            call checkAndHandleMpiErrorCode(ierr)
        enddo

        call Toc("VHReqTubeS")

    end subroutine

    subroutine vhHandleTubeStart(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        integer :: clockStart, clockEnd, clockRate, clockMax

        call recvChimpUpdate(vApp)

        ! load balancing works by reducing how large each instance of chimp thinks the shell grid is
        call mpi_recv(vApp%shGrid%js, 1, MPI_INTEGER, 0, 88311, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
        call checkAndHandleMpiErrorCode(ierr)
        call mpi_recv(vApp%shGrid%je, 1, MPI_INTEGER, 0, 88312, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
        call checkAndHandleMpiErrorCode(ierr)

        call system_clock(clockStart,clockRate,clockMax)
        call calcTubes(vApp)
        call system_clock(count=clockEnd)
        if(clockEnd < clockStart) then ! wrap
            vApp%lastTubeTime = real(clockMax-clockStart+clockEnd, rp)/real(clockRate, rp)
        else
            vApp%lastTubeTime = real(clockEnd-clockStart, rp)/real(clockRate, rp)
        endif

    end subroutine

    subroutine vhReqTubeEnd(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr, r, js, je, oldSizes(2), newSizes(2), offsets(2)
        type(MPI_Datatype) :: newtype

        call Tic("VHReqTubeE")
        call vhRequestType(vApp, VHTUBEEND)

        oldSizes = shape(vApp%State%ijTubes)
        do r=1,vApp%tubeLb%nL
            js = vApp%shGrid%js+vApp%tubeLb%balStartInd(r)
            if (r == vApp%tubeLb%nL) then
                je = vApp%shGrid%je
            else
                je = vApp%shGrid%js+vApp%tubeLb%balStartInd(r+1)-2
            endif

            newSizes = (/vApp%shGrid%ie+2-vApp%shGrid%is, &
                        je+2-js/)
            offsets = (/0, js-1/)

            ! send solved tube data
            call mpi_type_create_subarray(2, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, tubeMpiType, newtype, ierr)
            call checkAndHandleMpiErrorCode(ierr)
            call mpi_type_commit(newtype, ierr)
            call checkAndHandleMpiErrorCode(ierr)
            call mpi_recv(vApp%State%ijTubes, 1, newtype, r, 88131, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
            call checkAndHandleMpiErrorCode(ierr)
            call mpi_type_free(newtype, ierr)
            call checkAndHandleMpiErrorCode(ierr)
        enddo

        if(vApp%tubeLoadBalance) then
            ! collect recent timing data and update load balancing
            call mpi_gather(MPI_IN_PLACE, 0, MPI_MYFLOAT, vApp%tubeLb%instantTimes, &
                        1, MPI_MYFLOAT, 0, vApp%vhelpComm, ierr)
            call checkAndHandleMpiErrorCode(ierr)
            call updateLoads(vApp%tubeLb)
        endif

        call Toc("VHReqTubeE")

    end subroutine

    subroutine vhHandleTubeEnd(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr, oldSizes(2), newSizes(2), offsets(2)
        type(MPI_Datatype) :: newtype

        oldSizes = shape(vApp%State%ijTubes)
        newSizes = (/vApp%shGrid%ie+2-vApp%shGrid%is, &
                     vApp%shGrid%je+2-vApp%shGrid%js/)
        offsets = (/0, vApp%shGrid%js-1/)

        ! send solved tube data
        call mpi_type_create_subarray(2, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, tubeMpiType, newtype, ierr)
        call checkAndHandleMpiErrorCode(ierr)
        call mpi_type_commit(newtype, ierr)
        call checkAndHandleMpiErrorCode(ierr)
        call mpi_send(vApp%State%ijTubes, 1, newtype, 0, 88131, vApp%vHelpComm, ierr)
        call checkAndHandleMpiErrorCode(ierr)
        call mpi_type_free(newtype, ierr)
        call checkAndHandleMpiErrorCode(ierr)

        if(vApp%tubeLoadBalance) then
            call mpi_gather([vApp%lastTubeTime], 1, MPI_MYFLOAT, [vApp%lastTubeTime], 0, MPI_MYFLOAT, 0, vApp%vhelpComm, ierr)
            call checkAndHandleMpiErrorCode(ierr)
        endif

    end subroutine

    subroutine createTubeMpiType()

        type(Tube_T) :: dummyTube
        integer :: ierr
        integer(kind=MPI_ADDRESS_KIND) :: baseLoc,endLoc,typeSize,expectedSize
        integer :: blockLengths(4)
        integer(kind=MPI_ADDRESS_KIND) :: blockDisps(4)
        type(MPI_Datatype) :: blockTypes(4)

        tubeMpiType = MPI_DATATYPE_NULL        

        ! quick sanity check that Tube_T hasn't change since this was written
        call MPI_GET_ADDRESS(dummyTube,baseLoc, ierr)
        call checkAndHandleMpiErrorCode(ierr)
        call MPI_GET_ADDRESS(dummyTube%nTrc,endLoc, ierr)
        call checkAndHandleMpiErrorCode(ierr)
        typeSize = MPI_Aint_diff(endLoc,baseLoc)
        expectedSize = 392 ! determined when code was written
        if(typeSize /= expectedSize) then
            write (*,*) 'Tube_T object has changed, helpers cannot be used. Exitting.'
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        ! now create the custom type
        ! this, again, assumes that Tube_T is largely unchanged from when the code was written
        blockLengths(1) = 5+NDIM
        blockLengths(2) = 1
        blockLengths(3) = 13+NDIM+4*(1+MAXTUBEFLUIDS)
        blockLengths(4) = 1
        blockDisps(1) = 0
        blockDisps(2) = blockLengths(1)*8
        blockDisps(3) = blockDisps(2) + blockLengths(2)*8
        blockDisps(4) = blockDisps(3) + blockLengths(3)*8
        blockTypes(1) = MPI_MYFLOAT
        blockTypes(2) = MPI_INTEGER
        blockTypes(3) = MPI_MYFLOAT
        blockTypes(4) = MPI_INTEGER
        call mpi_type_create_struct(4,blockLengths,blockDisps,blockTypes,tubeMpiType,ierr)
        call checkAndHandleMpiErrorCode(ierr)
        call mpi_type_commit(tubeMpiType, ierr)
        call checkAndHandleMpiErrorCode(ierr)

    end subroutine

end module

