! implementation of coupling between voltron and mpi-decomposed gamera

module gamCouple_mpi_G2V
    use gamtypes_mpi
    use volttypes_mpi
    use couplingHelpers
    use gamapp_mpi, only: CalcDT_mpi
    use uservoltic ! required to have IonInnerBC_T defined

    implicit none

    type, extends(BaseOptions_T) :: gamOptionsCplMpiG_T
        type(MPI_Comm) :: couplingPoolComm

        contains
    end type

    ! type used on the gamera side
    type, extends(gamAppMpi_T) :: gamCouplerMpi_gam_T

        type(MPI_Comm) :: couplingComm
        integer :: myRank, voltRank
        logical :: processingData = .false.

        real(rp) :: DeepT
        logical :: doDeep

        ! array of all zeroes to simplify various send/receive calls
        integer, dimension(1) :: zeroArrayCounts = (/ 0 /)
        type(MPI_Datatype), dimension(1) :: zeroArrayTypes
        integer(kind=MPI_AN_MYADDR), dimension(1) :: zeroArrayDispls = (/ 0 /)

        ! Remix COUPLING VARIABLES
        integer, dimension(1) :: recvVCountsIneijk
        type(MPI_Datatype), dimension(1) :: recvVTypesIneijk
        integer(kind=MPI_AN_MYADDR), dimension(1) :: recvVDisplsIneijk
        integer, dimension(1) :: recvVCountsInexyz
        type(MPI_Datatype), dimension(1) :: recvVTypesInexyz
        integer(kind=MPI_AN_MYADDR), dimension(1) :: recvVDisplsInexyz

        ! Imag COUPLING VARIABLES
        integer, dimension(1) :: sendVCountsGas
        type(MPI_Datatype), dimension(1) :: sendVTypesGas
        integer(kind=MPI_AN_MYADDR), dimension(1) :: sendVDisplsGas
        integer, dimension(1) :: sendVCountsBxyz
        type(MPI_Datatype), dimension(1) :: sendVTypesBxyz
        integer(kind=MPI_AN_MYADDR), dimension(1) :: sendVDisplsBxyz
        integer, dimension(1) :: recvVCountsGas0
        type(MPI_Datatype), dimension(1) :: recvVTypesGas0
        integer(kind=MPI_AN_MYADDR), dimension(1) :: recvVDisplsGas0

        ! coupler-specific options
        type(gamOptionsCplMpiG_T) :: gOptionsCplMpiG

        contains

        ! only over-riding specific functions
        procedure :: InitModel => gamCplMpiGInitModel
        procedure :: InitIO => gamCplMpiGInitIO
        !procedure :: WriteRestart => gamCplWriteRestart
        !procedure :: ReadRestart => gamCplReadRestart
        !procedure :: WriteConsoleOutput => gamCplWriteConsoleOutput
        !procedure :: WriteFileOutput => gamCplWriteFileOutput
        !procedure :: WriteSlimFileOutput => gamCplWriteSlimFileOutput
        procedure :: AdvanceModel => gamCplMpiGAdvanceModel

    end type

    contains

    ! gamera to voltron functions

    subroutine gamCplMpiGInitModel(App, Xml)
        class(gamCouplerMpi_gam_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        integer :: length, commSize, ierr, numCells, dataCount, numInNeighbors, numOutNeighbors
        type(MPI_Comm) :: voltComm
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder, wasWeighted
        integer, dimension(1) :: rankArray, weightArray
        real(rp) :: tIO

        ! first initialize base gamera mpi
        call gamMpiInitModel(App, Xml)

        ! initialize F08 MPI objects
        App%couplingComm = MPI_COMM_NULL
        App%zeroArraytypes = (/ MPI_DATATYPE_NULL /)

        ! split voltron helpers off of the communicator
        ! split couplingPoolComm into a communicator with only the non-helper voltron rank
        call appWaitForVoltronSplit(App%gOptionsCplMpiG%couplingPoolComm, gamId, 0, voltComm)

        call Xml%Set_Val(App%doDeep, "/kaiju/voltron/coupling/doDeep", .true.)

        ! create a new communicator using MPI Topology
        call MPI_Comm_Size(voltComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Comm_rank(voltComm, App%myRank, ierr)
        ! identify who is voltron
        App%voltRank = -1
        call MPI_Allreduce(MPI_IN_PLACE, App%voltRank, 1, MPI_INTEGER, MPI_MAX, voltComm, ierr)

        ! send my i/j/k ranks to the voltron rank
        call mpi_gather(App%Grid%Ri, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, App%voltRank, voltComm, ierr)
        call mpi_gather(App%Grid%Rj, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, App%voltRank, voltComm, ierr)
        call mpi_gather(App%Grid%Rk, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, App%voltRank, voltComm, ierr)

        numCells = App%Grid%Nip*App%Grid%Njp*App%Grid%Nkp
        ! rank 0 send the number of physical cells to voltron rank
        if(App%myRank == 0) then
            call mpi_send(numCells, 1, MPI_INTEGER, commSize-1, 97500, voltComm, ierr)
        endif

        ! each gamera rank only talks to the voltron rank on this new communicator
        ! appoximate amount of data transfer as num physical cells
        dataCount = numCells ! deep data

        reorder = .true. ! allow MPI to reorder the ranks
        call mpi_dist_graph_create_adjacent(voltComm, &
            1,(/commSize-1/),(/dataCount/), &
            1,(/commSize-1/),(/dataCount/), &
            MPI_INFO_NULL, reorder, App%couplingComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Comm_rank(App%couplingComm, App%myRank, ierr)
        call mpi_dist_graph_neighbors_count(App%couplingComm,numInNeighbors,numOutNeighbors,wasWeighted,ierr)
        if(numInNeighbors /= 1) then
            print *,'Number of in edges was not 1 for rank ', App%myRank
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif
        if(numOutNeighbors /= 1) then
            print *,'Number of out edges was not 1 for rank ', App%myRank
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        ! get the rank of voltron in the new communicator
        call mpi_dist_graph_neighbors(App%couplingComm, numInNeighbors, rankArray, weightArray, &
                                      numOutNeighbors, rankArray, weightArray, ierr)
        App%voltRank = rankArray(1)

        ! send i/j/k ranks again since my rank may have changed in the new communicator
        call mpi_gather(App%Grid%Ri, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, App%voltRank, App%couplingComm, ierr)
        call mpi_gather(App%Grid%Rj, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, App%voltRank, App%couplingComm, ierr)
        call mpi_gather(App%Grid%Rk, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, App%voltRank, App%couplingComm, ierr)

        ! create the MPI datatypes for communicating state data with voltron
        call createG2VDataTypes(App)

        ! over-ride some of the initial voltron parameters on the gamera ranks
        if(.not. App%Model%isRestart) then
            ! don't over-ride restart time
            call mpi_bcast(App%Model%t, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)
        endif
        call mpi_bcast(App%Model%tFin, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%MJD0, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)

        ! correct DT now that we now the actual tFin
        call CalcDT_mpi(App)

        ! receive the initial coupling time
        call recvCplTimeMpi(App)

        ! couple one time to update mhd data on voltron node
        call sendVoltronCplDataMpi(App)

        ! then over-ride some IO terms from voltron
        call mpi_bcast(App%Model%IO%tRes, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%dtRes, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%nRes, 1, MPI_INTEGER, App%voltRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%tOut, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%dtOut, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%nOut, 1, MPI_INTEGER, App%voltRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%tCon, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)
        call mpi_bcast(App%Model%IO%dtCon, 1, MPI_MYFLOAT, App%voltRank, App%couplingComm, ierr)

        if(.not. App%Model%isRestart) then
            ! re-write Gamera's first output with corrected time, save and restore initial output time
            tIO = App%Model%IO%tOut
            call App%WriteFileOutput(App%Model%IO%nOut)
            App%Model%IO%tOut = tIO
        else
            ! never processing when restarted
            App%processingData = .false.
        endif

    end subroutine

    subroutine gamCplMpiGInitIO(App, Xml)
        class(gamCouplerMpi_gam_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        real(rp) :: save_tRes, save_dtRes, save_tOut, save_dtOut, save_tCon, save_dtCon
        integer :: save_nRes, save_nOut
        integer :: ierr

        ! save and restore parameters over-written by voltron
        ! doing it this way allows all comms and setup with voltron to happen in the init function
        ! and prevents some issues
        save_tRes = App%Model%IO%tRes
        save_dtRes = App%Model%IO%dtRes
        save_nRes = App%Model%IO%nRes
        save_tOut = App%Model%IO%tOut
        save_dtOut = App%Model%IO%dtOut
        save_nOut = App%Model%IO%nOut
        save_tCon = App%Model%IO%tCon
        save_dtCon = App%Model%IO%dtCon
        
        ! initialize parent's IO
        call gamInitIO(App, Xml)

        ! restore saved IO options
        App%Model%IO%tRes = save_tRes
        App%Model%IO%dtRes = save_dtRes
        App%Model%IO%nRes = save_nRes
        App%Model%IO%tOut = save_tOut
        App%Model%IO%dtOut = save_dtOut
        App%Model%IO%nOut = save_nOut
        App%Model%IO%tCon = save_tCon
        App%Model%IO%dtCon = save_dtCon

    end subroutine

    subroutine gamCplMpiGAdvanceModel(App, dt)
        class(gamCouplerMpi_gam_T), intent(inout) :: App
        real(rp), intent(in) :: dt

        real(rp) :: targetSimT 

        targetSimT = App%Model%t+dt

        ! ensure the model is advanced at least one step
        if(App%Model%t >= targetSimT) call gamMpiAdvanceModel(App, 0.0_rp)

        ! may need to step around coupling intervals
        do while(App%Model%t < targetSimT)
            if(.not. App%processingData) then
                ! receive new data to process
                call recvVoltronCplDataMpi(App)
                App%processingData = .true.
            elseif(App%DeepT <= App%Model%t) then
                ! send results and get new data
                call sendVoltronCplDataMpi(App)
                call recvVoltronCplDataMpi(App)
                App%processingData = .true.
            else
                if(targetSimT < App%DeepT) then
                    ! advance to the current step target time
                    call gamMpiAdvanceModel(App, targetSimT-App%Model%t)
                else
                    ! advance to next coupling time
                    call gamMpiAdvanceModel(App, App%DeepT-App%Model%t)
                    ! send results and get new data
                    call sendVoltronCplDataMpi(App)
                    call recvVoltronCplDataMpi(App)
                    App%processingData = .true.
                endif
            endif
        end do

    end subroutine

    subroutine recvVoltronCplDataMpi(gCplApp)
        class(gamCouplerMpi_gam_T), intent(inout) :: gCplApp


        call recvShallowCplDataMpi(gCplApp)
        if(gCplApp%doDeep) call recvDeepCplDataMpi(gCplApp)
        call recvCplTimeMpi(gCplApp)

    end subroutine

    subroutine recvShallowCplDataMpi(gCplApp)
        class(gamCouplerMpi_gam_T), intent(inout) :: gCplApp

        integer :: ierr

        ! The data goes into inEijk and inExyz in the IonInnerBC_T
        ! find the remix BC to write data to
        if(gCplApp%Grid%hasLowerBC(IDIR)) then
            SELECT type(iiBC=>gCplApp%Grid%externalBCs(INI)%p)
                TYPE IS (IonInnerBC_T)

                    ! Recv inEijk Data
                    call Tic("VoltSync", .true.)
                    call mpi_neighbor_alltoallw(0, gCplApp%zeroArrayCounts, &
                                                gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                                iiBC%inEijk, gCplApp%recvVCountsIneijk, &
                                                gCplApp%recvVDisplsIneijk, gCplApp%recvVTypesIneijk, &
                                                gCplApp%couplingComm, ierr)
                    call Toc("VoltSync", .true.)

                    ! Recv inExyz Data
                    call mpi_neighbor_alltoallw(0, gcplApp%zeroArrayCounts, &
                                                gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                                iiBC%inExyz, gCplApp%recvVCountsInexyz, &
                                                gCplApp%recvVDisplsInexyz, gCplApp%recvVTypesInexyz, &
                                                gCplApp%couplingComm, ierr)
                CLASS DEFAULT
                    write(*,*) 'Could not find Ion Inner BC in Voltron MPI performShallowUpdate'
                    stop
            END SELECT
        else
            ! not a rank with remix BC, but still need to call mpi_neighbor_alltoallw
            ! Recv nothing step 1
            call Tic("VoltSync", .true.)
            call mpi_neighbor_alltoallw(0, gCplApp%zeroArrayCounts, &
                                        gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                        0, gCplApp%zeroArrayCounts, &
                                        gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                        gCplApp%couplingComm, ierr)
            call Toc("VoltSync", .true.)

            ! Recv nothing step 2
            call mpi_neighbor_alltoallw(0, gCplApp%zeroArrayCounts, &
                                        gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                        0, gCplApp%zeroArrayCounts, &
                                        gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                        gCplApp%couplingComm, ierr)
        endif

    end subroutine

    subroutine recvDeepCplDataMpi(gCplApp)
        class(gamCouplerMpi_gam_T), intent(inout) :: gCplApp

        integer :: ierr

        ! Receive Gas0 Data
        call mpi_neighbor_alltoallw(gCplApp%Grid%Gas0, gCplApp%zeroArrayCounts, &
                                    gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                    gCplApp%Grid%Gas0, gCplApp%recvVCountsGas0, &
                                    gCplApp%recvVDisplsGas0, gCplApp%recvVTypesGas0, &
                                    gCplApp%couplingComm, ierr)

    end subroutine

    subroutine recvCplTimeMpi(gCplApp)
        class(gamCouplerMpi_gam_T), intent(inout) :: gCplApp

        integer :: ierr

        ! Receive next Coupling Target Time
        call mpi_bcast(gCplApp%DeepT,1,MPI_MYFLOAT, gCplApp%voltRank, gCplApp%couplingComm, ierr)

        ! Convert to Gamera time units
        gCplApp%DeepT = gCplApp%DeepT / gCplApp%Model%Units%gT0

    end subroutine

    subroutine sendVoltronCplDataMpi(gCplApp)
         class(gamCouplerMpi_gam_T), intent(inout) :: gCplApp

        integer :: ierr

        ! send state data to voltron
        ! Send Gas Data
        call mpi_neighbor_alltoallw(gCplApp%State%Gas, gCplApp%sendVCountsGas, &
                                    gCplApp%sendVDisplsGas, gCplApp%sendVTypesGas, &
                                    gCplApp%State%Gas, gCplapp%zeroArrayCounts, &
                                    gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                    gCplApp%couplingComm, ierr)
        ! Send Bxyz Data
        call mpi_neighbor_alltoallw(gCplApp%State%Bxyz, gCplApp%sendVCountsBxyz, &
                                    gCplApp%sendVDisplsBxyz, gCplApp%sendVTypesBxyz, &
                                    gCplApp%State%Bxyz, gCplApp%zeroArrayCounts, &
                                    gCplApp%zeroArrayDispls, gCplApp%zeroArrayTypes, &
                                    gCplApp%couplingComm, ierr)

    end subroutine

    subroutine createG2VDataTypes(gCplApp)
        type(gamCouplerMpi_gam_T), intent(inout) :: gCplApp

        integer :: ierr, dataSize, sendDataOffset, recvDataOffset
        type(MPI_Datatype) :: iPSI, iPSI1, Eijk2, EIjk3, Eijk4, Exyz2, Exyz3, Exyz4
        type(MPI_Datatype) :: iP,iPjP,iPjPkP,iPjPkP4Gas,iPjPkP4Bxyz,iPjPkP5Gas
        type(MPI_Datatype) :: iPG2,iPG2jPG2,iPG2jPG2kPG2,iPG2jPG2kPG24Gas,iPG2jPG2kPG25Gas

        associate(Grid=>gCplApp%Grid,Model=>gCplApp%Model)

        ! no need to allocate, all arrays are fixed size (1) since we only talk to voltron

        ! counts always 1
        ! displacements always 0
        ! datatypes to null by default
        gCplApp%recvVCountsInexyz = 1
        gCplApp%recvVCountsIneijk = 1
        gCplApp%recvVDisplsInexyz = 0
        gCplApp%recvVDisplsIneijk = 0
        gCplApp%recvVTypesInexyz = MPI_DATATYPE_NULL
        gCplApp%recvVTypesIneijk = MPI_DATATYPE_NULL

        gCplApp%sendVCountsGas = 1
        gCplApp%sendVCountsBxyz = 1
        gCplApp%sendVDisplsGas = 0
        gCplApp%sendVDisplsBxyz = 0
        gCplApp%sendVTypesGas = MPI_DATATYPE_NULL
        gCplApp%sendVTypesBxyz = MPI_DATATYPE_NULL

        if(gCplApp%doDeep) then
            gCplApp%recvVCountsGas0 = 1
            gCplApp%recvVDisplsGas0 = 0
            gCplApp%recvVTypesGas0 = MPI_DATATYPE_NULL
        endif

        ! assemble datatypes
        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per value

        ! I dimension
        call mpi_type_contiguous(PsiSh, MPI_MYFLOAT, iPSI, ierr) ! PsiSh i
        call mpi_type_contiguous(PsiSh+1, MPI_MYFLOAT, iPSI1, ierr) ! PsiSh+1 i
        call mpi_type_contiguous(Grid%Nip, MPI_MYFLOAT, iP, ierr)
        call mpi_type_contiguous(Grid%Ni, MPI_MYFLOAT, iPG2, ierr)

        ! J dimension
        call mpi_type_hvector(Grid%Nj, 1, PsiSh*dataSize, iPSI, Exyz2, ierr) ! PsiSh i - p+2g j
        call mpi_type_hvector(Grid%Nj+1, 1, (PsiSh+1)*dataSize, iPSI1, Eijk2, ierr) ! PsiSh+1 i - p+2g+1 j
        call mpi_type_hvector(Grid%Njp, 1, Grid%Ni*datasize, iP, iPjP, ierr)
        call mpi_type_hvector(Grid%Nj, 1, Grid%Ni*datasize, iPG2, iPG2jPG2, ierr)

        ! K dimension
        call mpi_type_hvector(Grid%Nk, 1, Grid%Nj*PsiSh*dataSize, &
                              Exyz2, Exyz3, ierr) ! PsiSh i - p+2g j - p+2g k
        call mpi_type_hvector(Grid%Nk+1, 1, (Grid%Nj+1)*(PsiSh+1)*dataSize,&
                              Eijk2, Eijk3, ierr) ! PsiSh+1 i - p+2g+1 j - p+2g+1 k
        call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*datasize, iPjP, iPjPkP, ierr)
        call mpi_type_hvector(Grid%Nk, 1, Grid%Ni*Grid%Nj*datasize, iPG2jPG2, iPG2jPG2kPG2, ierr)

        ! 4th dimension
        call mpi_type_hvector(NDIM, 1, PsiSh*Grid%Nj*Grid%Nk*dataSize, Exyz3, Exyz4, ierr)
        call mpi_type_hvector(NDIM, 1, (PsiSh+1)*(Grid%Nj+1)*(Grid%Nk+1)*dataSize, Eijk3, Eijk4, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*datasize, iPjPkP, iPjPkP4Gas, ierr)
        call mpi_type_hvector(NDIM, 1, Grid%Ni*Grid%Nj*Grid%Nk*datasize, iPjPkP, iPjPkP4Bxyz, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*datasize, iPG2jPG2kPG2, iPG2jPG2kPG24Gas, ierr)

        ! 5th dimension
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,iPjPkP4Gas,iPjPkP5Gas,ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,iPG2jPG2kPG24Gas, &
                              iPG2jPG2kPG25Gas,ierr)

        ! create appropriate MPI Datatypes
        if(.not. Grid%hasLowerBC(IDIR)) then
            ! only gamera ranks with lower I boundary participate in shallow updates
            gCplApp%recvVCountsInexyz = 0
            gCplApp%recvVCountsIneijk = 0
            ! set these types to non null because MPI complains
            gCplApp%recvVTypesInexyz = MPI_INTEGER
            gCplApp%recvVTypesIneijk = MPI_INTEGER
        else
            !Inexyz
            recvDataOffset = 0
            call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, Exyz4, &
                                   gCplApp%recvVTypesInexyz(1), ierr)

            !Ineijk
            recvDataOffset = 0
            call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, Eijk4, &
                                   gCplApp%recvVTypesIneijk(1), ierr)
        endif

        ! Gas
        sendDataOffset = Model%nG*Grid%Nj*Grid%Ni + &
                         Model%nG*Grid%Ni + &
                         Model%nG
        call mpi_type_hindexed(1,(/1/),sendDataOffset*dataSize,iPjPkP5Gas,gCplApp%sendVTypesGas(1),ierr)

        ! Bxyz
        call mpi_type_hindexed(1,(/1/),sendDataOffset*dataSize,iPjPkP4Bxyz,gCplApp%sendVTypesBxyz(1),ierr)

        if(gCplApp%doDeep) then
            ! Gas0
            recvDataOffset = 0
            call mpi_type_hindexed(1,(/1/),recvDataOffset*dataSize,iPG2jPG2kPG25Gas, &
                                   gCplApp%recvVTypesGas0(1),ierr)
        endif

        ! commit new mpi datatypes
        call mpi_type_commit(gCplApp%recvVTypesIneijk(1), ierr)
        call mpi_type_commit(gCplApp%recvVTypesInexyz(1), ierr)
        call mpi_type_commit(gCplApp%sendVTypesGas(1), ierr)
        call mpi_type_commit(gCplApp%sendVTypesBxyz(1), ierr)
        if(gCplApp%doDeep) then
            call mpi_type_commit(gCplApp%recvVTypesGas0(1), ierr)
        endif

        end associate

    end subroutine

end module

