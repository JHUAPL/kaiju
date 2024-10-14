module couplingHelpers
    use mpi_f08
    use mpidefs

    implicit none

    contains

    ! this helper routine wait its turn for Voltron to couple with this app
    ! participating in splitting until it is its turn
    subroutine appWaitForVoltronSplit(startingComm, appId, key, voltCoupledComm)
        type(MPI_Comm), intent(in) :: startingComm
        integer, intent(in) :: appId, key
        type(MPI_Comm), intent(inout) :: voltCoupledComm

        type(MPI_Comm) :: currentPool
        integer :: ierr, currentApp, voltRank

        currentApp = appId + 1 ! ensure the loop is entered
        currentPool = startingComm
        do while(currentApp .ne. appId)
            voltRank = -1
            call MPI_Allreduce(MPI_IN_PLACE, voltRank, 1, MPI_INTEGER, MPI_MAX, currentPool, ierr)
            call MPI_Bcast(currentApp, 1, MPI_INTEGER, voltRank, currentPool, ierr)
            if(currentApp == appId) then
                ! it's my turn, split with volt and then skip making the second comm
                call MPI_comm_split(currentPool, appId, key, voltCoupledComm, ierr)
                ! key is never used when making the exclusion pool, 0 is used to preserve order
                call MPI_comm_split(currentPool, MPI_UNDEFINED, 0, currentPool, ierr)
            else
                ! it's not my turn, don't split with volt, and then join the second comm
                call MPI_comm_split(currentPool, MPI_UNDEFINED, key, voltCoupledComm, ierr)
                ! key is never used when making the exclusion pool, 0 is used to preserve order
                call MPI_comm_split(currentPool, voltId, 0, currentPool, ierr)
            endif
        enddo

    end subroutine

    ! this helper routine allows voltron to couple with a specific app
    ! and create a smaller communicator containing the remaining apps
    subroutine voltronSplitWithApp(couplingPool, appId, key, coupledComm)
        type(MPI_Comm), intent(inout) :: couplingPool
        integer, intent(in) :: appId, key
        type(MPI_Comm), intent(inout) :: coupledComm

        integer :: ierr, myRank, appIdCpy
        appIdCpy = appId ! mpi_bcast doesn't interact well with intent(in)

        ! tell everyone I'm the broadcasting root
        ! broadcast which app I'm creating a communicator with, split with it, and then
        ! create a smaller pool that excludes that app
        call MPI_comm_rank(couplingPool, myRank, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, myRank, 1, MPI_INTEGER, MPI_MAX, couplingPool, ierr)
        call MPI_Bcast(appIdCpy, 1, MPI_INTEGER, myRank, couplingPool, ierr)
        call MPI_comm_split(couplingPool, appIdCpy, key, coupledComm, ierr)
        ! key is never used when making the exclusion pool, 0 is used to preserve order
        call MPI_comm_split(couplingPool, voltId, 0, couplingPool, ierr)

    end subroutine

end module

