module volttypes_mpi
    use volttypes
    use mpi_f08
    use mpidefs

    implicit none

    ! mpi voltron specific options
    type, extends(BaseOptions_T) :: VoltOptionsMpi_T
        type(MPI_Comm), pointer :: couplingPoolComm

        contains
    end type voltOptionsMpi_T

    type, extends(voltApp_T) :: voltAppMpi_T
        ! voltron to helpers comms variables
        type(MPI_Comm) :: vHelpComm
        integer :: vHelpRank
        real(rp) :: lastSquishTime
        type(MPI_Win) :: vHelpWin
        integer, dimension(:), allocatable :: vHelpIdle
        logical :: useHelpers = .false.
        logical :: doSquishHelp = .false., masterSquish = .false.
        logical :: squishLoadBalance = .true., deepProcessingInProgress=.false.

        ! coupling comms variables to be done on volt rank
        type(MPI_Comm) :: mageCplComm
        integer :: voltCplRank,CplSize
        integer :: gcmCplRank = -1, hidraNCplRank = -1,hidraSCplRank = -1,hidraCplRank = -1
        integer, dimension(:),allocatable :: IAm

        ! mpi voltron specific options
        type(VoltOptionsMpi_T) :: vOptionsMpi

    end type voltAppMpi_T

end module

