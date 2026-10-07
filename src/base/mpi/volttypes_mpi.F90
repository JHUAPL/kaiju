module volttypes_mpi
    use volttypes
    use coupled_layout_mpi
    use mpi_f08
    use mpidefs

    implicit none

    ! mpi voltron specific options
    type, extends(BaseOptions_T) :: VoltOptionsMpi_T
        type(MPI_Comm), pointer :: couplingPoolComm
        type(coupledLayout_T) :: mageLayout

        contains
    end type voltOptionsMpi_T

    type, extends(voltApp_T) :: voltAppMpi_T
        ! voltron to helpers comms variables
        type(MPI_Comm) :: vHelpComm
        integer :: vHelpRank
        type(MPI_Win) :: vHelpWin
        integer, dimension(:), allocatable :: vHelpIdle
        logical :: useHelpers = .false.

        real(rp) :: lastSquishTime
        type(loadBal_T) :: squishLb
        logical :: doSquishHelp = .false., masterSquish = .false.
        logical :: squishLoadBalance = .true., deepProcessingInProgress=.false.

        real(rp) :: lastTubeTime
        type(loadBal_T) :: tubeLb
        logical :: doTubeHelp = .false., tubeLoadbalance = .true.

        ! coupling comms variables to be done on volt rank
        type(MPI_Comm) :: mageCplComm
        integer :: voltCplRank,CplSize
        integer :: gcmCplRank = -1, hidraNCplRank = -1, hidraSCplRank = -1, hidraCplRank = -1
        integer :: gcmNumRanks, hidraNNumRanks, hidraSNumRanks
        integer, dimension(:), allocatable :: gcmRankList, hidraNRankList, hidraSRankList
        integer, dimension(:),allocatable :: IAm

        ! mpi voltron specific options
        type(VoltOptionsMpi_T) :: vOptionsMpi

    end type voltAppMpi_T

end module

