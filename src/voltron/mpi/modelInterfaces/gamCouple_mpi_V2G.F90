! implementation of coupling between voltron and mpi-decomposed gamera

module gamCouple_mpi_V2G
    use gamtypes_mpi
    use volttypes_mpi
    use couplingHelpers
    use uservoltic ! required to have IonInnerBC_T defined
    use init
    use gamCoupleHelper

    implicit none

    ! options for MPI voltron coupler
    type, extends(BaseOptions_T) :: gamOptionsCplMpiV_T
        type(MPI_Comm), pointer :: couplingPoolComm

        contains
    end type


    ! type used on the voltron side
    type, extends(gamCoupler_T) :: gamCouplerMpi_volt_T

        ! voltron to gamera comms variables
        type(MPI_Comm) :: couplingComm
        integer :: myRank
        logical :: doAsyncCoupling = .true.
        logical :: doDeep
        logical :: inProgress = .false.

        ! array of all zeroes to simplify various send/receive calls
        integer, dimension(:), allocatable :: zeroArrayCounts
        type(MPI_Datatype), dimension(:), allocatable :: zeroArrayTypes
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable ::  zeroArrayDispls

        ! list of gamera ranks to communicate with
        integer, dimension(:), allocatable :: sendRanks, recvRanks

        ! Remix COUPLING VARIABLES
        integer, dimension(:), allocatable :: sendGCountsIneijk
        type(MPI_Datatype), dimension(:), allocatable :: sendGTypesIneijk
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: sendGDisplsIneijk
        integer, dimension(:), allocatable :: sendGCountsInexyz
        type(MPI_Datatype), dimension(:), allocatable :: sendGTypesInexyz
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: sendGDisplsInexyz
        ! Remix ASYNCHRONOUS VARIABLES
        type(MPI_Request) :: ineijkGSendReq, inexyzGSendReq

        ! Imag COUPLING VARIABLES
        integer, dimension(:), allocatable :: recvGCountsGas
        type(MPI_Datatype), dimension(:), allocatable :: recvGTypesGas
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: recvGDisplsGas
        integer, dimension(:), allocatable :: recvGCountsBxyz
        type(MPI_Datatype), dimension(:), allocatable :: recvGTypesBxyz
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: recvGDisplsBxyz
        integer, dimension(:), allocatable :: sendGCountsGas0
        type(MPI_Datatype), dimension(:), allocatable :: sendGTypesGas0
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: sendGDisplsGas0
        logical :: imagProcessingInProgress = .false.

        ! coupler-specific options
        type(gamOptionsCplMpiV_T) :: gOptionsCplMpiV

        contains

        ! only over-riding specific functions
        procedure :: InitModel => gamCplMpiVInitModel
        procedure :: InitIO => gamCplMpiVInitIO
        procedure :: WriteRestart => gamCplMpiVWriteRestart
        procedure :: ReadRestart => gamCplMpiVReadRestart
        procedure :: WriteConsoleOutput => gamCplMpiVWriteConsoleOutput
        procedure :: WriteFileOutput => gamCplMpiVWriteFileOutput
        procedure :: WriteSlimFileOutput => gamCplMpiVWriteSlimFileOutput
        !procedure :: AdvanceModel => gamCplMpiVAdvanceModel
        procedure :: Cleanup => gamCplMpiVCleanup

        procedure :: InitMhdCoupler => gamCplMpiVInitMhdCoupler
        procedure :: StartUpdateMhdData => gamCplMpiVStartUpdateMhdData
        procedure :: FinishUpdateMhdData => gamCplMpiVFinishUpdateMhdData

    end type

    contains

    ! voltron to gamera functions

    subroutine gamCplMpiVInitModel(App, Xml)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        integer :: commSize, ierr, numCells, length, ic, numInNeighbors, numOutNeighbors
        type(MPI_Comm) :: voltComm
        integer :: gamNRES
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder, wasWeighted
        integer, allocatable, dimension(:) :: neighborRanks, inData, outData
        integer, allocatable, dimension(:) :: iRanks, jRanks, kRanks

        App%couplingComm = MPI_COMM_NULL
        App%ineijkGSendReq = MPI_REQUEST_NULL
        App%inexyzGSendReq = MPI_REQUEST_NULL

        ! only create gamera data structures so that it LOOKS like Gamera to other models
        App%Model%isLoud = .false.
        App%Grid%ijkShift(1:3) = 0
        call ReadCorners(App%Model,App%Grid,Xml,childGameraOpt=.true.)
        call SetRings(App%Model,App%Grid,Xml)
        call Corners2Grid(App%Model,App%Grid)
        call DefaultBCs(App%Model,App%Grid)
        call PrepState(App%Model,App%Grid,&
            App%oState,App%State,Xml,App%gOptions%userInitFunc)

        ! split couplingPoolComm into a communicator with only the non-helper voltron rank and Gamera ranks
        call voltronSplitWithApp(App%gOptionsCplMpiV%couplingPoolComm, gamId, 0, voltComm)

        call Xml%Set_Val(App%doAsyncCoupling,"/kaiju/voltron/coupling/doAsyncCoupling",.true.)
        call Xml%Set_Val(App%doDeep, "/kaiju/voltron/coupling/doDeep", .true.)

        ! create a new communicator using MPI Topology
        call MPI_Comm_Size(voltComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        allocate(neighborRanks(1:commSize-1))
        allocate(inData(1:commSize-1))
        allocate(outData(1:commSize-1))
        allocate(iRanks(1:commSize))
        allocate(jRanks(1:commSize))
        allocate(kRanks(1:commSize))

        allocate(App%zeroArrayCounts(1:commSize-1))
        allocate(App%zeroArrayTypes(1:commSize-1))
        allocate(App%zeroArrayDispls(1:commSize-1))
        App%zeroArrayCounts(:) = 0
        App%zeroArrayTypes(:) = MPI_INTEGER ! MPI_DATATYPE_NULL
        App%zeroArrayDispls(:) = 0

        ! doing a very very rough approximation of data transferred to help MPI reorder
        ! for deep updates, assume each rank sends data equal to its # physical cells

        call MPI_Comm_rank(voltComm, App%myRank, ierr)
        ! identify who is voltron
        call MPI_Allreduce(MPI_IN_PLACE, App%myRank, 1, MPI_INTEGER, MPI_MAX, voltComm, ierr)

        ! get i/j/k ranks from each Gamera mpi rank
	iRanks(:) = -1
        jRanks(:) = -1
        kRanks(:) = -1
        call mpi_gather(MPI_IN_PLACE, 1, MPI_INTEGER, iRanks, 1, MPI_INTEGER, App%myRank, voltComm, ierr)
        call mpi_gather(MPI_IN_PLACE, 1, MPI_INTEGER, jRanks, 1, MPI_INTEGER, App%myRank, voltComm, ierr)
        call mpi_gather(MPI_IN_PLACE, 1, MPI_INTEGER, kRanks, 1, MPI_INTEGER, App%myRank, voltComm, ierr)

        ! get the number of physical cells from rank 0
        call mpi_recv(numCells, 1, MPI_INTEGER, 0, 97500, voltComm, MPI_STATUS_IGNORE, ierr)

        do ic=1,commSize-1
            neighborRanks(ic) = ic-1
            inData(ic) = numCells
            outData(ic) = numCells
        enddo

        reorder = .true. ! allow MPI to reorder the ranks
         call mpi_dist_graph_create_adjacent(voltComm, &
            commSize-1,neighborRanks,inData, &
            commSize-1,neighborRanks,outData, &
            MPI_INFO_NULL, reorder, App%couplingComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Comm_rank(App%couplingComm, App%myRank, ierr)

        call mpi_dist_graph_neighbors_count(App%couplingComm,numInNeighbors,numOutNeighbors,wasWeighted,ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if (numInNeighbors /= numOutNeighbors) then
            print *,'Number of in edges and out edges did not match for voltron mpi'
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        allocate(App%sendRanks(numOutNeighbors))
        allocate(App%recvRanks(numInNeighbors))

        ! don't care about the weights, dump them into an existing array
        call mpi_dist_graph_neighbors(App%couplingComm, numInNeighbors, App%recvRanks, inData, &
                                      numOutNeighbors, App%sendRanks, inData, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! get i/j/k ranks again in case MPI ranks were reordered in the new communicator
	iRanks(:) = -1
        jRanks(:) = -1
        kRanks(:) = -1
        call mpi_gather(MPI_IN_PLACE, 1, MPI_INTEGER, iRanks, 1, MPI_INTEGER, App%myRank, App%couplingComm, ierr)
        call mpi_gather(MPI_IN_PLACE, 1, MPI_INTEGER, jRanks, 1, MPI_INTEGER, App%myRank, App%couplingComm, ierr)
        call mpi_gather(MPI_IN_PLACE, 1, MPI_INTEGER, kRanks, 1, MPI_INTEGER, App%myRank, App%couplingComm, ierr)

        ! create the MPI datatypes needed to transfer state data
        call createVoltDataTypes(App, iRanks, jRanks, kRanks)

        deallocate(neighborRanks, inData, outData, iRanks, jRanks, kRanks)

    end subroutine

    subroutine createVoltDataTypes(vApp, iRanks, jRanks, kRanks)
        class(gamCouplerMpi_volt_T), intent(inout) :: vApp
        integer, dimension(1:SIZE(vApp%recvRanks)+1), intent(in) :: iRanks, jRanks, kRanks

        integer :: ierr, NiRanks, NjRanks, NkRanks, NipT, NjpT, NkpT, dataSize
        integer :: r, rRank, recvDataOffset, sRank, sendDataOffset
        type(MPI_Datatype) :: recvDatatype
        type(MPI_Datatype) :: iPSI,iPSI1,Exyz2,Eijk2,Exyz3,Eijk3,Exyz4,Eijk4
        type(MPI_Datatype) :: iP,iPjP,iPjPkP,iPjPkP4Bxyz,iPjPkP4Gas,iPjPkP5Gas
        type(MPI_Datatype) :: iPG2,iPG2jPG2,iPG2jPG2kPG2,iPG2jPG2kPG24Gas,iPG2jPG2kPG25Gas

        associate(Grid=>vApp%Grid,Model=>vApp%Model)

        allocate(vApp%sendGCountsInexyz(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendGDisplsInexyz(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendGTypesInexyz(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendGCountsIneijk(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendGDisplsIneijk(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendGTypesIneijk(1:SIZE(vApp%sendRanks)))

        allocate(vApp%recvGCountsGas(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvGDisplsGas(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvGTypesGas(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvGCountsBxyz(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvGDisplsBxyz(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvGTypesBxyz(1:SIZE(vApp%recvRanks)))

        if(vApp%doDeep) then
            allocate(vApp%sendGCountsGas0(1:SIZE(vApp%sendRanks)))
            allocate(vApp%sendGDisplsGas0(1:SIZE(vApp%sendRanks)))
            allocate(vApp%sendGTypesGas0(1:SIZE(vApp%sendRanks)))
        endif

        ! counts are always 1 because we're sending a single (complicated) mpi datatype
        ! displacements are always 0 because the displacements are baked into each mpi datatype
        ! set all datatypes to null by default
        vApp%sendGCountsInexyz(:) = 1
        vApp%sendGCountsIneijk(:) = 1
        vApp%sendGDisplsInexyz(:) = 0
        vApp%sendGDisplsIneijk(:) = 0
        vApp%sendGTypesInexyz(:) = MPI_DATATYPE_NULL
        vApp%sendGTypesIneijk(:) = MPI_DATATYPE_NULL

        vApp%recvGCountsGas(:) = 1
        vApp%recvGCountsBxyz(:) = 1
        vApp%recvGDisplsGas(:) = 0
        vApp%recvGDisplsBxyz(:) = 0
        vApp%recvGTypesGas(:) = MPI_DATATYPE_NULL
        vApp%recvGTypesBxyz(:) = MPI_DATATYPE_NULL

        if(vApp%doDeep) then
            vApp%sendGCountsGas0(:) = 1
            vApp%sendGDisplsGas0(:) = 0
            vApp%sendGTypesGas0(:) = MPI_DATATYPE_NULL
        endif

        ! assemble the different datatypes
        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry
        NiRanks = maxval(iRanks)+1 ! gamera i decomposition
        NjRanks = maxval(jRanks)+1 ! gamera j decomposition
        NkRanks = maxval(kRanks)+1 ! gamera k decomposition
        NipT = Grid%Nip / NiRanks ! number of physical cells per gamera rank in the i direction
        NjpT = Grid%Njp / NjRanks ! number of physical cells per gamera rank in the j direction
        NkpT = Grid%Nkp / NkRanks ! number of physical cells per gamera rank in the k direction

        ! I dimension
        call mpi_type_contiguous(PsiSh, MPI_MYFLOAT, iPSI, ierr) ! PsiSh i
        call mpi_type_contiguous(PsiSh+1, MPI_MYFLOAT, iPSI1, ierr) ! PsiSh+1 i
        call mpi_type_contiguous(NipT, MPI_MYFLOAT, iP, ierr) ! physical i
        call mpi_type_contiguous(2*Model%nG+NipT, MPI_MYFLOAT, iPG2, ierr) ! physical + 2*ghosts i

         ! J dimension
        call mpi_type_hvector(2*Model%nG+NjpT, 1, PsiSh*dataSize, iPSI, Exyz2, ierr) ! PsiSh i - p+2g j
        call mpi_type_hvector(2*Model%nG+NjpT+1, 1, (PsiSh+1)*dataSize, iPSI1, Eijk2, ierr) ! PsiSh+1 i - p+2g+1 j
        call mpi_type_hvector(NjpT, 1, Grid%Ni*dataSize, iP, iPjP, ierr) ! physical i - physical j
        call mpi_type_hvector(2*Model%nG+NjpT, 1, Grid%Ni*dataSize, iPG2, iPG2jPG2, ierr) ! p+2g i - p+2g j

        ! K dimension - currently assume NO k decomposition
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Nj*PsiSh*dataSize, &
                              Exyz2, Exyz3, ierr) ! PsiSh i - p+2g j - p+2g k
        call mpi_type_hvector(2*Model%ng+NkpT+1, 1, (Grid%Nj+1)*(PsiSh+1)*dataSize,&
                              Eijk2, Eijk3, ierr) ! PsiSh+1 i - p+2g+1 j - p+2g+1 k
        call mpi_type_hvector(NkpT, 1, Grid%Ni*Grid%Nj*dataSize, iPjP, iPjPkP, ierr)
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Ni*Grid%Nj*dataSize, iPG2jPG2, iPG2jPG2kPG2, ierr)

        ! 4th dimension
        call mpi_type_hvector(NDIM, 1, PsiSh*Grid%Nj*Grid%Nk*dataSize, Exyz3, Exyz4, ierr)
        call mpi_type_hvector(NDIM, 1, (PsiSh+1)*(Grid%Nj+1)*(Grid%Nk+1)*dataSize, Eijk3, Eijk4, ierr)
        call mpi_type_hvector(NDIM, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPjPkP, iPjPkP4Bxyz, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPjPkP, iPjPkP4Gas, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPG2jPG2kPG2, iPG2jPG2kPG24Gas, ierr)

        ! 5th dimension
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,iPjPkP4Gas,iPjPkP5Gas,ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, &
                              iPG2jPG2kPG24Gas, iPG2jPG2kPG25Gas, ierr)

        ! figure out exactly what data needs to be sent to (and received from) each gamera rank
        ! create custom MPI datatypes to perform these transfers
        do r=1,SIZE(vApp%recvRanks)
            rRank = vApp%recvRanks(r)+1

            ! gas
            recvDataOffset = (Model%nG + kRanks(rRank)*NkpT)*Grid%Nj*Grid%Ni + &
                             (Model%nG + jRanks(rRank)*NjpT)*Grid%Ni + &
                             (Model%nG + iRanks(rRank)*NipT)
            call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, iPjPkP5Gas, &
                                   vApp%recvGTypesGas(r), ierr)
            ! Bxyz
            call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, iPjPkP4Bxyz, &
                                   vApp%recvGTypesBxyz(r), ierr)
        enddo

        do r=1,SIZE(vApp%sendRanks)
            sRank = vApp%sendRanks(r)+1

            if(iRanks(sRank) .gt. 0) then
                ! never send shallow data to any rank but minimum i
                vApp%sendGCountsInexyz(r) = 0
                vApp%sendGCountsIneijk(r) = 0
                ! set these types to non null because MPI complains
                vApp%sendGTypesInexyz(r) = MPI_INTEGER
                vApp%sendGTypesIneijk(r) = MPI_INTEGER
            else
                ! calculate the byte offset to the start of the data

                ! Inexyz
                sendDataOffset = kRanks(sRank)*NkpT*Grid%Nj*PsiSh + &
                                 jRanks(sRank)*NjpT*PsiSh

                call mpi_type_hindexed(1, (/1/), sendDataOffset*dataSize, Exyz4, &
                                       vApp%sendGTypesInexyz(r), ierr)

                !Ineijk
                sendDataOffset = kRanks(sRank)*NkpT*(Grid%Nj+1)*(PsiSh+1) + &
                                 jRanks(sRank)*NjpT*(PsiSh+1)

                call mpi_type_hindexed(1, (/1/), sendDataOffset*dataSize, Eijk4, &
                                       vApp%sendGTypesIneijk(r), ierr)
            endif

            if(vApp%doDeep) then
                ! gas0
                sendDataOffset = kRanks(sRank)*NkpT*Grid%Nj*Grid%Ni + &
                                 jRanks(sRank)*NjpT*Grid%Ni + &
                                 iRanks(sRank)*NipT
                call mpi_type_hindexed(1, (/1/), sendDataOffset*dataSize, iPG2jPG2kPG25Gas, &
                                       vApp%sendGTypesGas0(r), ierr)
            endif
        enddo

        do r=1,size(vApp%sendGTypesInexyz)
            call mpi_type_commit(vApp%sendGTypesInexyz(r), ierr)
            call mpi_type_commit(vApp%sendGTypesIneijk(r), ierr)

            call mpi_type_commit(vApp%recvGTypesGas(r), ierr)
            call mpi_type_commit(vApp%recvGTypesBxyz(r), ierr)

            if(vApp%doDeep) then
                call mpi_type_commit(vApp%sendGTypesGas0(r), ierr)
            endif
        enddo

        end associate

    end subroutine

    subroutine gamCplMpiVInitIO(App, Xml)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        ! initialize only my own IO, Gamera is in another process
        App%vh5File = trim(App%Model%RunID) // ".gamCpl.h5"

    end subroutine

    subroutine gamCplMpiVWriteRestart(App, Nres)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! write only my own restart data
        call writeGamCouplerRestart(App, nRes)

    end subroutine

    subroutine gamCplMpiVReadRestart(App, resId, nRes)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        ! read only my own restart data
        call readGamCouplerRestart(App, resId, nRes)

        ! we are never in progress after a restart
        App%inProgress = .false.

    end subroutine

    subroutine gamCplMpiVWriteConsoleOutput(App)
        class(gamCouplerMpi_volt_T), intent(inout) :: App

        real(rp) :: cpcp(2) = 0.0

        ! write only my own console info
        if(App%Model%isLoud) then
            call getCPCP(App%mixOutput,cpcp)
             write (*, '(a,2f8.3,a)')             '      CPCP = ' , cpcp(NORTH), cpcp(SOUTH), ' [kV, N/S]'
        endif

    end subroutine

    subroutine gamCplMpiVWriteFileOutput(App, nStep)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        integer, intent(in) :: nStep

        ! write ony my own file output
        call writeCouplerFileOutput(App, nStep)

    end subroutine

    subroutine gamCplMpiVWriteSlimFileOutput(App, nStep)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call gamCplMpiVWriteFileOutput(App, nStep)

    end subroutine

    subroutine gamCplMpiVCleanup(App)
        class(gamCouplerMpi_volt_T), intent(inout) :: App

        logical :: reqStat
        integer :: ierr

        call MPI_TEST(App%ineijkGSendReq,reqStat,MPI_STATUS_IGNORE,ierr)
        if(.not. reqStat) then
            ! async neighborhood ops don't support cancel
            call MPI_WAIT(App%ineijkGSendReq, MPI_STATUS_IGNORE, ierr)
        endif

        call MPI_TEST(App%inexyzGSendReq,reqStat,MPI_STATUS_IGNORE,ierr)
        if(.not. reqStat) then
            ! async neighborhood ops don't support cancel
            call MPI_WAIT(App%inexyzGSendReq, MPI_STATUS_IGNORE, ierr)
        endif
    end subroutine

    ! coupler routines

    subroutine gamCplMpiVInitMhdCoupler(App, voltApp)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: voltApp

        integer :: ierr

        ! call parent init function
        call gamInitMhdCoupler(App, voltApp)

        ! over-ride some of the initial voltron parameters on the gamera ranks
        ! local gamera has correct values from the above parent init function
        call mpi_bcast(App%Model%t, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%tFin, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%MJD0, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)

        ! send initial coupling time to Gamera
        call sendCplTimeMpi(App, voltApp%DeepT)

        ! perform initial coupling to collect data on voltron
        call recvGameraCplDataMpi(App)

        ! over-ride IO parameters on gamera during gamera's IO init subroutine
        call mpi_bcast(App%Model%IO%tRes, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%dtRes, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%nRes, 1, MPI_INTEGER, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%tOut, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%dtOut, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%nOut, 1, MPI_INTEGER, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%tCon, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%dtCon, 1, MPI_MYFLOAT, App%myRank, App%couplingComm, ierr)

    end subroutine

    subroutine gamCplMpiVStartUpdateMhdData(App, voltApp)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: voltApp

        if(App%inProgress) return

        App%inProgress = .true.

        call Tic("Coupling", .true.)
        call sendGameraCplDataMpi(App, voltApp%DeepT)
        call Toc("Coupling", .true.)

    end subroutine

    subroutine gamCplMpiVFinishUpdateMhdData(App, voltApp)
        class(gamCouplerMpi_volt_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: voltApp

        if(.not. App%inProgress) return

        call Tic("Coupling", .true.)
        call recvGameraCplDataMpi(App)
        call Toc("Coupling", .true.)

        App%inProgress = .false.

    end subroutine

    subroutine recvGameraCplDataMpi(gCplApp)
        class(gamCouplerMpi_volt_T), intent(inout) :: gCplApp

        integer :: ierr

        ! Receive Deep Gas Data
        call Tic("GameraSync", .true.)
        call mpi_neighbor_alltoallw(gCplApp%State%Gas, gCplApp%zeroArrayCounts, &
                                    gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                    gCplApp%State%Gas, gCplApp%recvGCountsGas, &
                                    gCplApp%recvGDisplsGas, gCplApp%recvGTypesGas, &
                                    gCplApp%couplingComm, ierr)
        call Toc("GameraSync", .true.)
        ! Receive Deep Bxyz Data
        call mpi_neighbor_alltoallw(gCplApp%State%Bxyz, gCplApp%zeroArrayCounts, &
                                    gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                    gCplApp%State%Bxyz, gCplApp%recvGCountsBxyz, &
                                    gCplApp%recvGDisplsBxyz, gCplApp%recvGTypesBxyz, &
                                    gCplApp%couplingComm, ierr)

    end subroutine

    subroutine sendGameraCplDataMpi(gCplApp, CouplingTargetT)
        class(gamCouplerMpi_volt_T), intent(inout) :: gCplApp
        real(rp), intent(in) :: CouplingTargetT

        call sendShallowCplDataMpi(gCplApp)
        if(gCplApp%doDeep) call sendDeepCplDataMpi(gCplApp)
        call sendCplTimeMpi(gCplApp, CouplingTargetT)

    end subroutine

    subroutine sendShallowCplDataMpi(gCplApp)
        class(gamCouplerMpi_volt_T), intent(inout) :: gCplApp

        integer :: ierr

        ! voltron updates inEijk and inExyz in the IonInnerBC_T
        ! find the remix BC to read data from
        SELECT type(iiBC=>gCplApp%Grid%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                if(gCplApp%doAsyncCoupling) then
                    ! asynchronous
                    call mpi_wait(gCplApp%ineijkGSendReq, MPI_STATUS_IGNORE, ierr)
                    call mpi_Ineighbor_alltoallw(iiBC%inEijk, gCplApp%sendGCountsIneijk, &
                                                 gCplApp%sendGDisplsIneijk, gCplApp%sendGTypesIneijk, &
                                                 0, gCplApp%zeroArrayCounts, &
                                                 gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                                 gCplApp%couplingComm, gCplApp%ineijkGSendReq, ierr)

                    call mpi_wait(gCplApp%inexyzGSendReq, MPI_STATUS_IGNORE, ierr)
                    call mpi_Ineighbor_alltoallw(iiBC%inExyz, gCplApp%sendGCountsInexyz, &
                                                 gCplApp%sendGDisplsInexyz, gCplApp%sendGTypesInexyz, &
                                                 0, gCplApp%zeroArrayCounts, &
                                                 gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                                 gCplApp%couplingComm, gCplApp%inexyzGSendReq, ierr)
                else
                    ! synchronous
                    call mpi_neighbor_alltoallw(iiBC%inEijk, gCplApp%sendGCountsIneijk, &
                                                gCplApp%sendGDisplsIneijk, gCplApp%sendGTypesIneijk, &
                                                0, gCplApp%zeroArrayCounts, &
                                                gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                                gCplApp%couplingComm, ierr)

                    call mpi_neighbor_alltoallw(iiBC%inExyz, gCplApp%sendGCountsInexyz, &
                                                gCplApp%sendGDisplsInexyz, gCplApp%sendGTypesInexyz, &
                                                0, gCplApp%zeroArrayCounts, &
                                                gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                                gCplApp%couplingComm, ierr)
                endif
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in sendGameraCplDataMpi'
                stop
        END SELECT

    end subroutine

    subroutine sendDeepCplDataMpi(gCplApp)
        class(gamCouplerMpi_volt_T), intent(inout) :: gCplApp

        integer :: ierr

        ! Send Deep Gas0 Data
        call mpi_neighbor_alltoallw(gCplApp%Grid%Gas0, gCplApp%sendGCountsGas0, &
                                    gCplApp%sendGDisplsGas0, gCplApp%sendGTypesGas0, &
                                    gCplApp%Grid%Gas0, gCplApp%zeroArrayCounts, &
                                    gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                    gCplApp%couplingComm, ierr)

    end subroutine

    subroutine sendCplTimeMpi(gCplApp, CouplingTargetT)
        class(gamCouplerMpi_volt_T), intent(inout) :: gCplApp
        real(rp), intent(in) :: CouplingTargetT

        integer :: ierr

        ! Send Target Time for next coupling
        call mpi_bcast(CouplingTargetT,1,MPI_MYFLOAT, gCplApp%myRank, gCplApp%couplingComm, ierr)

    end subroutine

end module

