module gamtypes_mpi
    use gamtypes
    use mpidefs
    use mpi_f08

    implicit none

    type, extends(BaseOptions_T) :: gamOptionsMpi_T
        type(MPI_Comm) :: gamComm
        logical :: doIO = .true.

        contains
    end type gamOptionsMpi_T

    type, extends(gamApp_T) :: gamAppMpi_T
        type(MPI_Comm) :: gamMpiComm
        integer, dimension(:), allocatable :: sendRanks, recvRanks
        logical :: blockHalo = .false.

        ! Gas Data Transfer Variables
        integer, dimension(:), allocatable :: sendCountsGas
        type(MPI_Datatype), dimension(:), allocatable :: sendTypesGas
        integer, dimension(:), allocatable :: recvCountsGas
        type(MPI_Datatype), dimension(:), allocatable :: recvTypesGas
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: sendDisplsGas, recvDisplsGas

        ! Magnetic Flux Data Transfer Variables
        integer, dimension(:), allocatable :: sendCountsMagFlux
        type(MPI_Datatype), dimension(:), allocatable :: sendTypesMagFlux
        integer, dimension(:), allocatable :: recvCountsMagFlux
        type(MPI_Datatype), dimension(:), allocatable :: recvTypesMagFlux
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: sendDisplsMagFlux, recvDisplsMagFlux

        ! Bxyz Data Transfer Variables
        integer, dimension(:), allocatable :: sendCountsBxyz
        type(MPI_Datatype), dimension(:), allocatable :: sendTypesBxyz
        integer, dimension(:), allocatable :: recvCountsBxyz
        type(MPI_Datatype), dimension(:), allocatable :: recvTypesBxyz
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: sendDisplsBxyz, recvDisplsBxyz

        ! Debugging flags
        logical :: printMagFluxFaceError = .false.
        real(rp) :: faceError = 0.0_rp
        logical :: slowestRankPrints = .true.

        ! MPI-specific options
        type(gamOptionsMpi_T) :: gOptionsMpi

        contains

        ! only over-riding specific functions
        procedure :: InitModel => gamMpiInitModel
        !procedure :: InitIO => gamInitIO
        !procedure :: WriteRestart => gamWriteRestart
        !procedure :: ReadRestart => gamReadRestart
        procedure :: WriteConsoleOutput => gamMpiWriteConsoleOutput
        !procedure :: WriteFileOutput => gamWriteFileOutput
        !procedure :: WriteSlimFileOutput => gamWriteSlimFileOutput
        procedure :: AdvanceModel => gamMpiAdvanceModel

    end type gamAppMpi_T

    ! placeholders for app procedures
    interface
        module subroutine gamMpiInitModel(App, Xml)
            class(gamAppMpi_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine gamMpiInitModel

        module subroutine gamMpiWriteConsoleOutput(App)
            class(gamAppMpi_T), intent(inout) :: App
        end subroutine gamMpiWriteConsoleOutput

        module subroutine gamMpiAdvanceModel(App, dt)
            class(gamAppMpi_T), intent(inout) :: App
            real(rp), intent(in) :: dt
        end subroutine gamMpiAdvanceModel
    end interface

    contains

end module

