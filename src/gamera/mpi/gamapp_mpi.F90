! Main data objects and functions to perform a gamera simulation

module gamapp_mpi
    use gamtypes_mpi
    use step
    use init
    use mhdgroup
    use gamapp
    use bcs_mpi
    use mpidefs
    use mpi_f08

    implicit none

    contains

    subroutine initGamera_mpi(gamAppMpi, xmlInp)
        class(gamAppMpi_T), intent(inout) :: gamAppMpi
        type(XML_Input_T), intent(inout) :: xmlInp

        character(len=strLen) :: kaijuRoot
        logical :: doLoud
        integer :: rank,ierr

        gamAppMpi%Model%isLoud = .true.

        ! initialize F08 MPI objects
        gamAppMpi%gamMpiComm = MPI_COMM_NULL

        !Verbose for root Gamera rank only
        !NOTE: Doing this w/ direct MPI call b/c isTiled/etc doesn't get set until later
        call mpi_comm_rank(gamAppMpi%gOptionsMpi%gamComm, rank, ierr)
        if (rank == 0) then
            doLoud = .true.
        else
            doLoud = .false.
        endif

        call xmlInp%SetVerbose(doLoud)
        call xmlInp%SetRootStr('Kaiju/Gamera')

        ! try to verify that the XML file has "Kaiju" as a root element
        kaijuRoot = ""
        call xmlInp%Get_Key_Val("/Gamera/sim/H5Grid",kaijuRoot,.false.)
        if(len(trim(kaijuRoot)) /= 0) then
            write(*,*) "The input XML appears to be of an old style."
            write(*,*) "As of June 12th, 2021 it needs a root element of <Kaiju>."
            write(*,*) "Please modify your XML config file by adding this line at the top:"
            write(*,*) "<Kaiju>"
            write(*,*) "and this line at the bottom:"
            write(*,*) "</Kaiju>"
            write(*,*) "OR (preferred) convert your configuration to an INI file and use"
            write(*,*) " the XMLGenerator.py script to create conforming XML files."
            write(*,*) "Please refer to the python script or"
            write(*,*) " the [Generating XML Files] wiki page for additional info."
            stop
        endif

        call xmlInp%Set_Val(gamAppMpi%blockHalo,"coupling/blockHalo",.false.)

        ! read debug flags
        call xmlInp%Set_Val(writeGhosts,"debug/writeGhosts",.false.)
        call xmlInp%Set_Val(writeMagFlux,"debug/writeMagFlux",.false.)
        call xmlInp%Set_Val(gamAppMpi%printMagFluxFaceError,"debug/printMagFluxError",.false.)
        call xmlInp%Set_Val(gamAppMpi%slowestRankPrints,"debug/slowestRankPrints",.true.)

        !Initialize Grid/State/Model (Hatch Gamera)
        !Will enforce 1st BCs, caculate 1st timestep, set oldState
        call Hatch_mpi(gamAppMpi,xmlInp,gamAppMpi%gOptions%userInitFunc,gamAppMpi%gOptionsMpi%gamComm)
        call cleanClocks()

        if (gamAppMpi%gOptionsMpi%doIO) then
            if (.not. gamAppMpi%Model%isRestart) call fOutput(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)
            call consoleOutput(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)
        endif

    end subroutine initGamera_mpi

    subroutine Hatch_mpi(gamAppMpi,xmlInp,userInitFunc,gamComm,endTime)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        type(XML_Input_T), intent(inout) :: xmlInp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        type(MPI_Comm), intent(in) :: gamComm
        real(rp), optional, intent(in) :: endTime

        integer :: numNeighbors, ierr, length, commSize, rank, ic, jc, kc, rn
        integer :: targetRank, sendNumber, numInNeighbors, numOutNeighbors, listIndex
        integer, dimension(27) :: sourceRanks, sourceData
        logical :: reorder,periodicI,periodicJ,periodicK,wasWeighted
        character(len=strLen) :: message
        real(rp), dimension(:,:,:), allocatable :: tempX,tempY,tempZ
        real(rp) :: tmpDT
        character(len=strLen) :: inH5
        logical :: doH5g, allLoud
        integer, dimension(NDIM) :: dims
        integer, dimension(:), allocatable :: gamRestartNumbers

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        ! setup the distributed graph topology communicator for future MPI work
        call xmlInp%Set_Val(Grid%NumRi,'iPdir/N',1)
        call xmlInp%Set_Val(Grid%NumRj,'jPdir/N',1)
        call xmlInp%Set_Val(Grid%NumRk,'kPdir/N',1)
        call xmlInp%Set_Val(periodicI,'iPdir/bcPeriodic',.false.)
        call xmlInp%Set_Val(periodicJ,'jPdir/bcPeriodic',.false.)
        call xmlInp%Set_Val(periodicK,'kPdir/bcPeriodic',.false.)

        call xmlInp%Set_Val(doH5g ,"sim/doH5g" ,.false.)

        ! debugging flag to make all gamera ranks loud
        call xmlInp%Set_Val(allLoud, "output/allLoud", .false.)

        if(Grid%NumRi*Grid%NumRj*Grid%NumRk > 1) then
            Grid%isTiled = .true.
        endif

        Grid%ijkShift(:) = 0

        if(Grid%isTiled) then
            call mpi_comm_size(gamComm, commSize, ierr)
            if(commSize /= Grid%NumRi*Grid%NumRj*Grid%NumRk) then
                print *,'Expected ',Grid%NumRi*Grid%NumRj*Grid%NumRk,' ranks, but this job has ',commSize,' gamera ranks'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

            ! read number of physical cells from H5 file or xml
            if(doH5g .or. Model%isRestart) then
                ! Read basic grid size info from full grid file
                call xmlInp%Set_Val(inH5,"sim/H5Grid","gMesh.h5")

                dims = GridSizeH5(inH5)

                ! These assume Model%nG is 4
                Grid%Nip = dims(1) - 9
                Grid%Njp = dims(2) - 9
                Grid%Nkp = dims(3) - 9
            else
                call xmlInp%Set_Val(Grid%Nip,"idir/N",1)
                call xmlInp%Set_Val(Grid%Njp,"jdir/N",1)
                call xmlInp%Set_Val(Grid%Nkp,"kdir/N",1)
            endif

            if(modulo(Grid%Nip,Grid%NumRi) /= 0) then
                print *,'Number of cells in i dimension is ',Grid%Nip,' but it must be a multiple of the number of MPI ranks in that dimension, which was ',Grid%NumRi
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif
            if(modulo(Grid%Njp,Grid%NumRj) /= 0) then
                print *,'Number of cells in j dimension is ',Grid%Njp,' but it must be a multiple of the number of MPI ranks in that dimension, which was ',Grid%NumRj
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif
            if(modulo(Grid%Nkp,Grid%NumRk) /= 0) then
                print *,'Number of cells in k dimension is ',Grid%Nkp,' but it must be a multiple of the number of MPI ranks in that dimension, which was ',Grid%NumRk
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

            ! rank assignment. These are now 0 to N-1, matching the way that MPI counts it
            call mpi_comm_rank(gamComm, rank, ierr)
            Grid%Rk = modulo(rank, Grid%NumRk)
            rank = (rank-Grid%Rk)/Grid%NumRk
            Grid%Rj = modulo(rank, Grid%NumRj)
            rank = (rank-Grid%Rj)/Grid%NumRj
            Grid%Ri = rank

            ! adjust corner information to reflect this individual node's grid data
            Grid%Nip = Grid%Nip/Grid%NumRi
            Grid%Njp = Grid%Njp/Grid%NumRj
            Grid%Nkp = Grid%Nkp/Grid%NumRk

            ! assume Model%nG is 4
            Grid%Ni = Grid%Nip + 8
            Grid%Nj = Grid%Njp + 8
            Grid%Nk = Grid%Nkp + 8

            ! create MPI topology
            sourceRanks(:) = -1
            sourceData(:) = 0
            numNeighbors = 0
            call mpi_comm_rank(gamComm, rank, ierr)
            do ic=-1,1
                do jc=-1,1
                    do kc=-1,1
                        targetRank = 0

                        ! I rank offset
                        if((Grid%Ri+ic >= 0 .and. Grid%Ri+ic < Grid%NumRi) .or. periodicI) then
                            ! talking to a neighbor, or wrapping around periodic MPI boundaries
                            targetRank = targetRank + modulo(Grid%Ri+ic,Grid%NumRi)*Grid%NumRk*Grid%NumRj
                        else
                            ! talking to a non-periodic boundary, so this is my own I rank
                            targetRank = targetRank + Grid%Ri*Grid%NumRk*Grid%NumRj
                        endif

                        ! J rank offset
                        if((Grid%Rj+jc >= 0 .and. Grid%Rj+jc < Grid%NumRj) .or. periodicJ) then
                            ! talking to a neighbor, or wrapping around periodic MPI boundaries
                            targetRank = targetRank + modulo(Grid%Rj+jc,Grid%NumRj)*Grid%NumRk
                        else
                            ! talking to a non-periodic boundary, so this is my own J rank
                            targetRank = targetRank + Grid%Rj*Grid%NumRk
                        endif

                        ! K rank offset
                        if((Grid%Rk+kc >= 0 .and. Grid%Rk+kc < Grid%NumRk) .or. periodicK) then
                            ! talking to a neighbor, or wrapping around periodic MPI boundaries
                            targetRank = targetRank + modulo(Grid%Rk+kc,Grid%NumRk)
                        else
                            ! talking to a non-periodic boundary, so this is my own K rank
                            targetRank = targetRank + Grid%Rk
                        endif
                       
                        if(targetRank /= rank) then ! ensure I'm not talking to myself
#if defined __INTEL_COMPILER && __INTEL_COMPILER >= 1800
                        listIndex = findloc(sourceRanks, targetRank, 1)
#else
                        !Bypass as findloc does not work for gfortran<9
                        !Work-around code
                        listIndex = 0
                        do rn=1,size(sourceRanks)
                            if (sourceRanks(rn) .eq. targetRank) then
                                listIndex = rn
                                exit
                            endif
                        enddo
#endif
                            if(listIndex == 0) then ! this rank not in the list yet
                                numNeighbors = numNeighbors+1
                                listIndex = numNeighbors
                                sourceRanks(listIndex) = targetRank
                            endif
                            ! calculate the size of this region that will be transmitted
                            sendNumber = 1
                            if (ic == 0) then
                                sendNumber = sendNumber * Grid%Nip
                            else
                                sendNumber = sendNumber * Model%nG
                            endif
                            if (jc == 0) then
                                sendNumber = sendNumber * Grid%Njp
                            else
                                sendNumber = sendNumber * Model%nG
                            endif
                            if (kc == 0) then
                                sendNumber = sendNumber * Grid%Nkp
                            else
                                sendNumber = sendNumber * Model%nG
                            endif
                            ! add these cells to whatever we are already sending to that rank
                            sourceData(listIndex) = sourceData(listIndex) + sendNumber
                        endif
                    enddo
                enddo
            enddo

            reorder = .true. ! allow MPI to reorder the ranks
            call mpi_dist_graph_create_adjacent(gamComm, &
                numNeighbors,sourceRanks,sourceData, &
                numNeighbors,sourceRanks,sourceData, & ! comms symmetrical
                MPI_INFO_NULL, reorder, gamAppMpi%gamMpiComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            ! get our final ranks from the new MPI topology
            call mpi_comm_rank(gamAppMpi%gamMpiComm, rank, ierr)
            Grid%Rk = modulo(rank, Grid%NumRk)
            rank = (rank-Grid%Rk)/Grid%NumRk
            Grid%Rj = modulo(rank, Grid%NumRj)
            rank = (rank-Grid%Rj)/Grid%NumRj
            Grid%Ri = rank

            ! only rank 0 should be loud
            if((Grid%Ri==0 .and. Grid%Rj==0 .and. Grid%Rk==0) .or. allLoud) then
                Model%isLoud = .true.
            else
                Model%isLoud = .false.
            endif

            ! whether this rank has external BCs
            Grid%hasLowerBC(IDIR) = Grid%Ri == 0
            Grid%hasLowerBC(JDIR) = Grid%Rj == 0
            Grid%hasLowerBC(KDIR) = Grid%Rk == 0
            Grid%hasUpperBC(IDIR) = Grid%Ri == (Grid%NumRi-1)
            Grid%hasUpperBC(JDIR) = Grid%Rj == (Grid%NumRj-1)
            Grid%hasUpperBC(KDIR) = Grid%Rk == (Grid%NumRk-1)

        endif

        ! call appropriate subroutines to read corner info and mesh size data
        call ReadCorners(Model,Grid,xmlInp,endTime)

        if(Grid%isTiled) then
            ! now create the arrays that MPI will use to send and receive the data
            call mpi_dist_graph_neighbors_count(gamAppMpi%gamMpiComm,numInNeighbors,numOutNeighbors,wasWeighted,ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            if (numInNeighbors /= numOutNeighbors) then
                print *,'Number of in edges and out edges did not match for rank ', rank
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

            allocate(gamAppMpi%sendRanks(numOutNeighbors))
            allocate(gamAppMpi%recvRanks(numInNeighbors))
            ! don't care about the weights, dump them into an existing array
            call mpi_dist_graph_neighbors(gamAppMpi%gamMpiComm, numInNeighbors, gamAppMpi%recvRanks, sourceData, numOutNeighbors, gamAppMpi%sendRanks, sourceData, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            ! Create MPI datatypes
            call createDatatypes(gamAppMpi, periodicI, periodicJ, periodicK)

        endif

        ! call appropriate subroutines to calculate all appropriate grid data from the corner data
        call CalcGridInfo(Model,Grid,gamAppMpi%State,gamAppMpi%oState,gamAppMpi%Solver,xmlInp,userInitFunc)

        ! All Gamera ranks compare restart numbers to ensure they're the same
        if(Model%isRestart .and. Grid%isTiled) then
            if(Grid%Ri==0 .and. Grid%Rj==0 .and. Grid%Rk==0) then
                ! master rank receives data
                allocate(gamRestartNumbers(commSize))
                call mpi_gather(gamAppMpi%Model%IO%nRes, 1, MPI_INTEGER, gamRestartNumbers, 1, MPI_INTEGER, 0, gamAppMpi%gamMpiComm, ierr)
                if(.not. all(gamRestartNumbers .eq. minval(gamRestartNumbers))) then
                    write(*,*) "Gamera ranks did not all agree on restart numbers, you should sort that out."
                    write(*,*) "Error code: A house divided cannot stand"
                    write(*,*) "   Minimum Gamera nRes = ", minval(gamRestartNumbers)
                    write(*,*) "   Maximum Gamera nRes = ", maxval(gamRestartNumbers)
                    stop
                endif
                deallocate(gamRestartNumbers)
            else
                ! all other ranks only send data
                call mpi_gather(gamAppMpi%Model%IO%nRes, 1, MPI_INTEGER, 0, 0, MPI_INTEGER, 0, gamAppMpi%gamMpiComm, ierr)
            endif
        endif

        if(Grid%isTiled) then
            ! correct boundary conditions if necessary
            if(Grid%NumRi > 1) then
                ! MPI decomposed in I dimension

                ! check min I BC
                SELECT type(iiBC=>Grid%externalBCs(INI)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasLowerBC(IDIR) .and. .not. periodicI) then
                            print *, 'Min I BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasLowerBC(IDIR) .or. periodicI) then
                            ! print *, 'Over-writing min I BC to be an MPI BC'
                            deallocate(Grid%externalBCs(INI)%p)
                            allocate(Grid%externalBCs(INI)%p,source=mpiNullBc_T(INI))
                        endif
                END SELECT

                ! check max I BC
                SELECT type(iiBC=>Grid%externalBCs(OUTI)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasUpperBC(IDIR) .and. .not. periodicI) then
                            print *, 'Max I BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasUpperBC(IDIR) .or. periodicI) then
                            ! print *, 'Over-writing max I BC to be an MPI BC'
                            deallocate(Grid%externalBCs(OUTI)%p)
                            allocate(Grid%externalBCs(OUTI)%p,source=mpiNullBc_T(OUTI))
                        endif
                END SELECT

            endif

            if(Grid%NumRj > 1) then
                ! MPI decomposed in J dimension

                ! check min J BC
                SELECT type(iiBC=>Grid%externalBCs(INJ)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasLowerBC(JDIR) .and. .not. periodicJ) then
                            print *, 'Min J BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasLowerBC(JDIR) .or. periodicJ) then
                            ! print *, 'Over-writing min J BC to be an MPI BC'
                            deallocate(Grid%externalBCs(INJ)%p)
                            allocate(Grid%externalBCs(INJ)%p,source=mpiNullBc_T(INJ))
                        endif
                END SELECT

                ! check max J BC
                SELECT type(iiBC=>Grid%externalBCs(OUTJ)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasUpperBC(JDIR) .and. .not. periodicJ) then
                            print *, 'Max J BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasUpperBC(JDIR) .or. periodicJ) then
                            ! print *, 'Over-writing max J BC to be an MPI BC'
                            deallocate(Grid%externalBCs(OUTJ)%p)
                            allocate(Grid%externalBCs(OUTJ)%p,source=mpiNullBc_T(OUTJ))
                        endif
                END SELECT

            endif

            if(Grid%NumRk > 1) then
                ! MPI decomposed in K dimension

                ! check min K BC
                SELECT type(iiBC=>Grid%externalBCs(INK)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasLowerBC(KDIR) .and. .not. periodicK) then
                            print *, 'Min K BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasLowerBC(KDIR) .or. periodicK) then
                            ! print *, 'Over-writing min K BC to be an MPI BC'
                            deallocate(Grid%externalBCs(INK)%p)
                            allocate(Grid%externalBCs(INK)%p,source=mpiNullBc_T(INK))
                        endif
                END SELECT

                ! check max K BC
                SELECT type(iiBC=>Grid%externalBCs(OUTK)%p)
                    TYPE IS (mpiNullBc_T)
                        ! this BC is already an MPI BC
                        if(Grid%hasUpperBC(KDIR) .and. .not. periodicK) then
                            print *, 'Max K BC was set to be an MPI BC by the user, but this is not an MPI periodic case'
                            stop
                        endif
                    CLASS DEFAULT
                        ! this BC is something other than an MPI BC
                        if(.not. Grid%hasUpperBC(KDIR) .or. periodicK) then
                            ! print *, 'Over-writing max K BC to be an MPI BC'
                            deallocate(Grid%externalBCs(OUTK)%p)
                            allocate(Grid%externalBCs(OUTK)%p,source=mpiNullBc_T(OUTK))
                        endif
                END SELECT

            endif

            ! Make sure BCs are valid and setup for correct directions
            call ValidateBCs(gamAppMpi%Model, gamAppMpi%Grid)

            ! update all ghost values
            call updateMpiBCs(gamAppMpi, gamAppMpi%State)

            !Ensure all processes have the same starting timestep
            call CalcDT_mpi(gamAppMpi)

            if(Model%isRestart) then
                !Update the ghost cells in the old state
                call updateMpiBCs(gamAppMpi, gamAppMpi%oState)
            else
                !Manufacture an old state
                gamAppMpi%oState = gamAppMpi%State
                gamAppMpi%oState%time = gamAppMpi%State%time-Model%dt !Initial old state
            endif

        endif

        end associate
    end subroutine Hatch_mpi

    ! each gamera rank calculates their own DT, and then synchronizes them to use
    ! whoever has the lowest one
    subroutine CalcDT_mpi(gamAppMpi)
        class(gamAppMpi_T), intent(inout) :: gamAppMpi
        real(rp) :: tmpDT
        integer :: ierr

        call Tic("DT")
        gamAppMpi%Model%dt = CalcDT(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)

        if(gamAppMpi%gamMpiComm /= MPI_COMM_NULL) then
            call Tic("mpiDT")
            tmpDT = gamAppMpi%Model%dt
            call MPI_AllReduce(MPI_IN_PLACE, tmpDT, 1, MPI_MYFLOAT, MPI_MIN, gamAppMpi%gamMpiComm,ierr)
            gamAppMpi%Model%dt = tmpDT
            call Toc("mpiDT")
        endif
        call Toc("DT")

    end subroutine CalcDT_mpi

    ! this function checks logic to determine which rank should print debug timing
    ! returns true for the rank that should, and false for all others
    function debugPrintingRank(gamAppMpi) result(amPrintingRank)
        type(gamAppMpi_T), intent(in) :: gamAppMpi
        logical :: amPrintingRank
        real(rp) :: inData(2), outData(2)
        integer :: myRank, slowRank, ierr

        amPrintingRank = .false.

        if(.not. gamAppMpi%Grid%isTiled) then
            ! only a single rank
            amPrintingRank = .true.
        else if(gamAppMpi%slowestRankPrints) then
            ! only the slowest rank prints

            myRank = gamAppMpi%Grid%Ri*gamAppMpi%Grid%NumRj*gamAppMpi%Grid%NumRk + &
                     gamAppMpi%Grid%Rj*gamAppMpi%Grid%NumRk + &
                     gamAppMpi%Grid%Rk

            inData(1) = gamAppMpi%Model%kzcsMHD
            inData(2) = myRank
            call MPI_AllReduce(inData, outData, 1, MPI_2MYFLOAT, MPI_MINLOC, gamAppMpi%gamMpiComm, ierr)
            slowRank = outData(2) ! convert rank back to an integer
            if(myRank == slowRank) then
                amPrintingRank = .true.
            endif
        else
            ! only rank 0 prints
            amPrintingRank = gamAppMpi%Grid%Ri==0 .and. gamAppMpi%Grid%Rj==0 .and. gamAppMpi%Grid%Rk==0
        endif

    end function debugPrintingRank

    subroutine consoleOutput_mpi(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        real(rp) :: totalFaceError = 0.0_rp
        integer :: ierr

        !Always call console output (for certain performance metrics) but only be loud if root
        call consoleOutput(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)

        ! gamera mpi specific output
        if(gamAppMpi%printMagFluxFaceError) then
            if(gamAppMpi%Grid%isTiled) then
                call MPI_AllReduce(gamAppMpi%faceError, totalFaceError, 1, MPI_MYFLOAT, MPI_SUM, gamAppMpi%gamMpiComm, ierr)
            else
                totalFaceError = gamAppMpi%faceError
            endif
            if (gamAppMpi%Model%isLoud) then
                write(*,*) ANSICYAN
                write(*,*) 'GAMERA MPI'
                write(*,'(a,f8.3)') 'Total MF Face Error = ', totalFaceError
                write(*,'(a)',advance="no") ANSIRESET!, ''
            endif
        endif

    end subroutine consoleOutput_mpi

    subroutine updateMpiBCs(gamAppMpi, State)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        type(State_T), intent(inout) :: State

        integer :: i,ierr
        character(len=strLen) :: BCID

        !Enforce BCs
        call Tic("BCs")
        call EnforceBCs(gamAppMpi%Model,gamAppMpi%Grid,State)
        call Toc("BCs")

        !Track timing for all gamera ranks to finish physical BCs
        ! Only synchronize when timing
        if(gamAppMpi%Model%IO%doTimerOut .and. gamAppMpi%Grid%isTiled) then
            call Tic("Sync BCs")
            call MPI_BARRIER(gamAppMpi%gamMpiComm,ierr)
            call Toc("Sync BCs")
        endif

        !Update ghost cells
        call Tic("Halos")
        call HaloUpdate(gamAppMpi, State)
        call bFlux2Fld(gamAppMpi%Model, gamappMpi%Grid, State%magFlux, State%Bxyz) !Update Bxyz's
        call Toc("Halos")

        !Track timing for all gamera ranks to finish halo comms
        ! Only synchronize when timing
        if(gamAppMpi%Model%IO%doTimerOut .and. gamAppMpi%Grid%isTiled) then
            call Tic("Sync Halos")
            call MPI_BARRIER(gamAppMpi%gamMpiComm,ierr)
            call Toc("Sync Halos")
        endif

        ! Re-apply periodic BCs last
        do i=1,gamAppMpi%Grid%NumBC
            if(allocated(gamAppMpi%Grid%externalBCs(i)%p)) then
                SELECT type(bc=>gamAppMpi%Grid%externalBCs(i)%p)
                    TYPE IS (periodicInnerIBC_T)
                        write (BCID, '(A,I0)') "BC#", i
                        call Tic(BCID)
                        call bc%doBC(gamAppMpi%Model,gamAppMpi%Grid,State)
                        call Toc(BCID)
                    TYPE IS (periodicOuterIBC_T)
                        write (BCID, '(A,I0)') "BC#", i
                        call Tic(BCID)
                        call bc%doBC(gamAppMpi%Model,gamAppMpi%Grid,State)
                        call Toc(BCID)
                    TYPE IS (periodicInnerJBC_T)
                        write (BCID, '(A,I0)') "BC#", i
                        call Tic(BCID)
                        call bc%doBC(gamAppMpi%Model,gamAppMpi%Grid,State)
                        call Toc(BCID)
                    TYPE IS (periodicOuterJBC_T)
                        write (BCID, '(A,I0)') "BC#", i
                        call Tic(BCID)
                        call bc%doBC(gamAppMpi%Model,gamAppMpi%Grid,State)
                        call Toc(BCID)
                    TYPE IS (periodicInnerKBC_T)
                        write (BCID, '(A,I0)') "BC#", i
                        call Tic(BCID)
                        call bc%doBC(gamAppMpi%Model,gamAppMpi%Grid,State)
                        call Toc(BCID)
                    TYPE IS (periodicOuterKBC_T)
                        write (BCID, '(A,I0)') "BC#", i
                        call Tic(BCID)
                        call bc%doBC(gamAppMpi%Model,gamAppMpi%Grid,State)
                        call Toc(BCID)
                    CLASS DEFAULT
                        ! do nothing, not a periodic BC
                endselect
            endif
        enddo

        !Track timing for all gamera ranks to finish periodic BCs
        ! Only synchronize when timing
        if(gamAppMpi%Model%IO%doTimerOut .and. gamAppMpi%Grid%isTiled) then
            call Tic("Sync Periodics")
            call MPI_BARRIER(gamAppMpi%gamMpiComm,ierr)
            call Toc("Sync Periodics")
        endif

    end subroutine updateMpiBCs

    subroutine stepGamera_mpi(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi

        integer :: ierr,i
        real(rp) :: tmp

        !update the state variables to the next timestep
        call UpdateStateData(gamAppMpi)

        !Track timing for all gamera ranks to finish math
        ! Only synchronize when timing
        if(gamAppMpi%Model%IO%doTimerOut .and. gamAppMpi%Grid%isTiled) then
            call Tic("Sync Math")
            call MPI_BARRIER(gamAppMpi%gamMpiComm,ierr)
            call Toc("Sync Math")
        endif

        !Update BCs MPI style
        call updateMpiBCs(gamAppMpi, gamAppmpi%State)

        !Calculate new timestep
        call CalcDT_mpi(gamAppMpi)

    end subroutine stepGamera_mpi

    subroutine haloUpdate(gamAppMpi, State)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        type(State_T), intent(inout) :: State

        integer :: ierr, length
        type(MPI_Request) :: gasReq, mfReq, bReq
        character(len=strLen) :: message

        ! arrays for calculating mag flux face error, if applicable
        real(rp), dimension(:,:), allocatable :: maxI
        real(rp), dimension(:,:), allocatable :: maxJ
        real(rp), dimension(:,:), allocatable :: maxK

        if(gamAppMpi%gamMpiComm /= MPI_COMM_NULL) then
            ! just tell MPI to use the arrays we defined during initialization to send and receive data!

            ! Gas Cell Data
            if(gamAppMpi%blockHalo) then
                ! synchronous
                call mpi_neighbor_alltoallw(State%Gas, gamAppMpi%sendCountsGas, &
                                            gamAppMpi%sendDisplsGas, gamAppMpi%sendTypesGas, &
                                            State%Gas, gamAppMpi%recvCountsGas, &
                                            gamAppMpi%recvDisplsGas, gamAppMpi%recvTypesGas, &
                                            gamAppMpi%gamMpiComm, ierr)
            else
                !asynchronous
                call mpi_ineighbor_alltoallw(State%Gas, gamAppMpi%sendCountsGas, &
                                            gamAppMpi%sendDisplsGas, gamAppMpi%sendTypesGas, &
                                            State%Gas, gamAppMpi%recvCountsGas, &
                                            gamAppMpi%recvDisplsGas, gamAppMpi%recvTypesGas, &
                                            gamAppMpi%gamMpiComm, gasReq, ierr)
            endif

            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

            if(gamAppMpi%Model%doMHD) then
                if(gamAppMpi%printMagFluxFaceError) then
                    if(.not. allocated(maxI)) allocate(maxI(gamAppMpi%grid%js:gamAppMpi%grid%je,gamAppMpi%grid%ks:gamAppMpi%grid%ke))
                    if(.not. allocated(maxJ)) allocate(maxJ(gamAppMpi%grid%is:gamAppMpi%grid%ie,gamAppMpi%grid%ks:gamAppMpi%grid%ke))
                    if(.not. allocated(maxK)) allocate(maxK(gamAppMpi%grid%is:gamAppMpi%grid%ie,gamAppMpi%grid%js:gamAppMpi%grid%je))
                    maxI = State%magFlux(gamAppMpi%grid%ie+1,gamAppMpi%grid%js:gamAppMpi%grid%je,gamAppMpi%grid%ks:gamAppMpi%grid%ke,IDIR)
                    maxJ = State%magFlux(gamAppMpi%grid%is:gamAppMpi%grid%ie,gamAppMpi%grid%je+1,gamAppMpi%grid%ks:gamAppMpi%grid%ke,JDIR)
                    maxK = State%magFlux(gamAppMpi%grid%is:gamAppMpi%grid%ie,gamAppMpi%grid%js:gamAppMpi%grid%je,gamAppMpi%grid%ke+1,KDIR)
                endif

                ! Magnetic Face Flux Data and Cell-Centered Data
                if(gamAppMpi%blockHalo) then
                    ! synchronous
                    call mpi_neighbor_alltoallw(State%magFlux, gamAppMpi%sendCountsMagFlux, &
                                                gamAppMpi%sendDisplsMagFlux, gamAppMpi%sendTypesMagFlux, &
                                                State%magFlux, gamAppMpi%recvCountsMagFlux, &
                                                gamAppMpi%recvDisplsMagFlux, gamAppMpi%recvTypesMagFlux, &
                                                gamAppMpi%gamMpiComm, ierr)
                    call mpi_neighbor_alltoallw(State%Bxyz, gamAppMpi%sendCountsBxyz, &
                                                gamAppMpi%sendDisplsBxyz, gamAppMpi%sendTypesBxyz, &
                                                State%Bxyz, gamAppMpi%recvCountsBxyz, &
                                                gamAppMpi%recvDisplsBxyz, gamAppMpi%recvTypesBxyz, &
                                                gamAppMpi%gamMpiComm, ierr)
                else
                    ! asynchronous
                    call mpi_ineighbor_alltoallw(State%magFlux, gamAppMpi%sendCountsMagFlux, &
                                                gamAppMpi%sendDisplsMagFlux, gamAppMpi%sendTypesMagFlux, &
                                                State%magFlux, gamAppMpi%recvCountsMagFlux, &
                                                gamAppMpi%recvDisplsMagFlux, gamAppMpi%recvTypesMagFlux, &
                                                gamAppMpi%gamMpiComm, mfReq, ierr)
                    call mpi_ineighbor_alltoallw(State%Bxyz, gamAppMpi%sendCountsBxyz, &
                                                gamAppMpi%sendDisplsBxyz, gamAppMpi%sendTypesBxyz, &
                                                State%Bxyz, gamAppMpi%recvCountsBxyz, &
                                                gamAppMpi%recvDisplsBxyz, gamAppMpi%recvTypesBxyz, &
                                                gamAppMpi%gamMpiComm, bReq, ierr)
                endif

                if(ierr /= MPI_Success) then
                    call MPI_Error_string( ierr, message, length, ierr)
                    print *,message(1:length)
                    call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
                endif
                if(gamAppMpi%printMagFluxFaceError) then
                    gamAppMpi%faceError = &
                      sum(abs(State%magFlux(gamAppMpi%grid%ie+1,gamAppMpi%grid%js:gamAppMpi%grid%je,gamAppMpi%grid%ks:gamAppMpi%grid%ke,IDIR) - maxI)) + &
                      sum(abs(State%magFlux(gamAppMpi%grid%is:gamAppMpi%grid%ie,gamAppMpi%grid%je+1,gamAppMpi%grid%ks:gamAppMpi%grid%ke,JDIR) - maxJ)) + &
                      sum(abs(State%magFlux(gamAppMpi%grid%is:gamAppMpi%grid%ie,gamAppMpi%grid%js:gamAppMpi%grid%je,gamAppMpi%grid%ke+1,KDIR) - maxK))
                endif

            endif

            if(.not. gamAppMpi%blockHalo) then
                call mpi_wait(gasReq, MPI_STATUS_IGNORE, ierr)
                if(gamAppMpi%Model%doMHD) then
                    call mpi_wait(mfReq, MPI_STATUS_IGNORE, ierr)
                    call mpi_wait(bReq, MPI_STATUS_IGNORE, ierr)
                endif
            endif

        endif

    end subroutine haloUpdate

    subroutine createDatatypes(gamAppMpi, periodicI, periodicJ, periodicK)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        logical, intent(in) :: periodicI, periodicJ, periodicK

        integer :: iData,jData,kData,rankIndex,offset,dataSize,ierr
        type(MPI_Datatype) :: dType,dtGas4,dtGas5,dtBxyz4

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        allocate(gamAppMpi%sendCountsGas(1:SIZE(gamAppMpi%sendRanks)))
        allocate(gamAppMpi%sendDisplsGas(1:SIZE(gamAppMpi%sendRanks)))
        allocate(gamAppMpi%sendTypesGas(1:SIZE(gamAppMpi%sendRanks)))
        allocate(gamAppMpi%recvCountsGas(1:SIZE(gamAppMpi%recvRanks)))
        allocate(gamAppMpi%recvDisplsGas(1:SIZE(gamAppMpi%recvRanks)))
        allocate(gamAppMpi%recvTypesGas(1:SIZE(gamAppMpi%recvRanks)))

        ! counts are always 1 because we're sending a single (complicated) mpi datatype
        gamAppMpi%sendCountsGas(:) = 1
        gamAppMpi%recvCountsGas(:) = 1

        ! displacements are always 0 because the displacements are baked into each mpi datatype
        gamAppMpi%sendDisplsGas(:) = 0
        gamAppMpi%recvDisplsGas(:) = 0

        ! set all datatypes to null by default
        gamAppMpi%sendTypesGas(:) = MPI_DATATYPE_NULL
        gamAppMpi%recvTypesGas(:) = MPI_DATATYPE_NULL

        if(Model%doMHD) then
            ! do setup for mpi transfer of MHD data
            allocate(gamAppMpi%sendCountsMagFlux(1:SIZE(gamAppMpi%sendRanks)))
            allocate(gamAppMpi%sendDisplsMagFlux(1:SIZE(gamAppMpi%sendRanks)))
            allocate(gamAppMpi%sendTypesMagFlux(1:SIZE(gamAppMpi%sendRanks)))
            allocate(gamAppMpi%recvCountsMagFlux(1:SIZE(gamAppMpi%recvRanks)))
            allocate(gamAppMpi%recvDisplsMagFlux(1:SIZE(gamAppMpi%recvRanks)))
            allocate(gamAppMpi%recvTypesMagFlux(1:SIZE(gamAppMpi%recvRanks)))
            allocate(gamAppMpi%sendCountsBxyz(1:SIZE(gamAppMpi%sendRanks)))
            allocate(gamAppMpi%sendDisplsBxyz(1:SIZE(gamAppMpi%sendRanks)))
            allocate(gamAppMpi%sendTypesBxyz(1:SIZE(gamAppMpi%sendRanks)))
            allocate(gamAppMpi%recvCountsBxyz(1:SIZE(gamAppMpi%recvRanks)))
            allocate(gamAppMpi%recvDisplsBxyz(1:SIZE(gamAppMpi%recvRanks)))
            allocate(gamAppMpi%recvTypesBxyz(1:SIZE(gamAppMpi%recvRanks)))

            ! counts are always 1 because we're sending a single (complicated) mpi datatype
            gamAppMpi%sendCountsMagFlux(:) = 1
            gamAppMpi%recvCountsMagFlux(:) = 1
            gamAppMpi%sendCountsBxyz(:) = 1
            gamAppMpi%recvCountsBxyz(:) = 1

            ! displacements are always 0 because the displacements are baked into each mpi datatype
            gamAppMpi%sendDisplsMagFlux(:) = 0
            gamAppMpi%recvDisplsMagFlux(:) = 0
            gamAppMpi%sendDisplsBxyz(:) = 0
            gamAppMpi%recvDisplsBxyz(:) = 0

            ! set all datatypes to null by default
            gamAppMpi%sendTypesMagFlux(:) = MPI_DATATYPE_NULL
            gamAppMpi%recvTypesMagFlux(:) = MPI_DATATYPE_NULL
            gamAppMpi%sendTypesBxyz(:) = MPI_DATATYPE_NULL
            gamAppMpi%recvTypesBxyz(:) = MPI_DATATYPE_NULL
        endif

        ! figure out exactly what data needs to be sent to (and received from) each neighbor
        ! create custom MPI datatypes to perform these transfers
        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr)
        do iData = -1,1
            do jData=-1,1
                do kData=-1,1
                    do rankIndex = 1,SIZE(gamappMpi%recvRanks)
                        ! for each possibly adjacent rank

                        ! Start cell centered
                        ! Gas and Bxyz do faces and edges, not corners
                        call calcRecvDatatypeOffsetCC(gamAppMpi,gamAppMpi%recvRanks(rankIndex),iData,jData,kData,&
                                                      periodicI,periodicJ,periodicK,dType,offset,&
                                                      .true.,.true.,.false.)
                        if(dType /= MPI_DATATYPE_NULL) then
                            ! add 4th and 5th dimensions for cell centered Gas
                            call mpi_type_hvector(NVAR,1,Grid%Ni*Grid%Nj*Grid%Nk*dataSize,dType,dtGas4,ierr)
                            call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,&
                                                  dtGas4,dtGas5,ierr)
                            call appendDatatype(gamAppMpi%recvTypesGas(rankIndex),dtGas5,offset)

                            if(Model%doMHD) then
                                ! also add a 4th dimension for cell centered Bxyz
                                call mpi_type_hvector(NDIM,1,Grid%Ni*Grid%Nj*Grid%Nk*dataSize,dType,dtBxyz4,ierr)
                                call appendDatatype(gamAppMpi%recvTypesBxyz(rankIndex),dtBxyz4,offset)
                            endif
                        endif

                        ! Gas and Bxyz do faces and edges, not corners
                        call calcSendDatatypeOffsetCC(gamAppMpi,gamAppmpi%sendRanks(rankIndex),iData,jData,kData,&
                                                      periodicI,periodicJ,periodicK,dType,offset,&
                                                      .true.,.true.,.false.)
                        if(dType /= MPI_DATATYPE_NULL) then
                            ! add 4th and 5th dimensions for cell centered Gas
                            call mpi_type_hvector(NVAR,1,Grid%Ni*Grid%Nj*Grid%Nk*dataSize,dType,dtGas4,ierr)
                            call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,&
                                                  dtGas4,dtGas5,ierr)
                            call appendDatatype(gamAppMpi%sendTypesGas(rankIndex),dtGas5,offset)

                            if(Model%doMHD) then
                                ! also add a 4th dimension for cell centered Bxyz
                                call mpi_type_hvector(NDIM,1,Grid%Ni*Grid%Nj*Grid%Nk*dataSize,dType,dtBxyz4,ierr)
                                call appendDatatype(gamAppMpi%sendTypesBxyz(rankIndex),dtBxyz4,offset)
                            endif
                        endif
                        ! End cell centered

                        ! Start face centered
                        if(Model%doMHD) then
                            ! magFlux does faces and edges, not corners
                            call calcRecvDatatypeOffsetFC(gamAppMpi,gamAppMpi%recvRanks(rankIndex),iData,&
                                                          jData,kData,periodicI,periodicJ,periodicK,dType,offset,&
                                                          .true.,.true.,.false.)
                            if(dType /= MPI_DATATYPE_NULL) then
                                ! face centered datatype already has all 4 dimensions
                                call appendDatatype(gamAppMpi%recvTypesMagFlux(rankIndex),dType,offset)
                            endif

                            call calcSendDatatypeOffsetFC(gamAppMpi,gamAppmpi%sendRanks(rankIndex),iData,&
                                                    jData,kData,periodicI,periodicJ,periodicK,dType,offset,&
                                                    .true.,.true.,.false.)
                            if(dType /= MPI_DATATYPE_NULL) then
                                ! face centered datatype already has all 4 dimensions
                                call appendDatatype(gamAppMpi%sendTypesMagFlux(rankIndex),dType,offset)
                            endif
                        endif
                        ! End face centered
                    enddo
                enddo
            enddo
        enddo

        ! commit the created MPI datatypes
        do rankIndex=1,SIZE(gamAppMpi%sendRanks)
            if(gamAppMpi%sendTypesGas(rankIndex) == MPI_DATATYPE_NULL) then
                ! no comms between these ranks
                gamAppMpi%sendCountsGas(rankIndex) = 0
                gamAppMpi%sendDisplsGas(rankIndex) = 0
                gamAppMpi%sendTypesGas(rankIndex) = MPI_INTEGER ! NULL causes crash
            else
                ! there are comms between these ranks
                call mpi_type_commit(gamAppMpi%sendTypesGas(rankIndex), ierr)
            endif
            if(gamAppMpi%recvTypesGas(rankIndex) == MPI_DATATYPE_NULL) then
                gamAppMpi%recvCountsGas(rankIndex) = 0
                gamAppMpi%recvDisplsGas(rankIndex) = 0
                gamAppMpi%recvTypesGas(rankIndex) = MPI_INTEGER
            else
                call mpi_type_commit(gamAppMpi%recvTypesGas(rankIndex), ierr)
            endif
            if(Model%doMHD) then
                if(gamAppMpi%sendTypesMagFlux(rankIndex) == MPI_DATATYPE_NULL) then
                    gamAppMpi%sendCountsMagFlux(rankIndex) = 0
                    gamAppMpi%sendDisplsMagFlux(rankIndex) = 0
                    gamAppMpi%sendTypesMagFlux(rankIndex) = MPI_INTEGER
                else
                    call mpi_type_commit(gamAppMpi%sendTypesMagFlux(rankIndex), ierr)
                endif
                if(gamAppMpi%recvTypesMagFlux(rankIndex) == MPI_DATATYPE_NULL) then
                    gamAppMpi%recvCountsMagFlux(rankIndex) = 0
                    gamAppMpi%recvDisplsMagFlux(rankIndex) = 0
                    gamAppMpi%recvTypesMagFlux(rankIndex) = MPI_INTEGER
                else
                    call mpi_type_commit(gamAppMpi%recvTypesMagFlux(rankIndex), ierr)
                endif
                if(gamAppMpi%sendTypesBxyz(rankIndex) == MPI_DATATYPE_NULL) then
                    gamAppMpi%sendCountsBxyz(rankIndex) = 0
                    gamAppMpi%sendDisplsBxyz(rankIndex) = 0
                    gamAppMpi%sendTypesBxyz(rankIndex) = MPI_INTEGER
                else
                    call mpi_type_commit(gamAppMpi%sendTypesBxyz(rankIndex), ierr)
                endif
                if(gamAppMpi%recvTypesBxyz(rankIndex) == MPI_DATATYPE_NULL) then
                    gamAppMpi%recvCountsBxyz(rankIndex) = 0
                    gamAppMpi%recvDisplsBxyz(rankIndex) = 0
                    gamAppMpi%recvTypesBxyz(rankIndex) = MPI_INTEGER
                else
                    call mpi_type_commit(gamAppMpi%recvTypesBxyz(rankIndex), ierr)
                endif
            endif
        enddo

        end associate

    end subroutine createDatatypes

    subroutine calcRecvDatatypeOffsetCC(gamAppMpi,recvFromRank,iData,jData,kData,&
                                  periodicI,periodicJ,periodicK,dType,offset,&
                                  doFace,doEdge,doCorner)
        type(gamAppMpi_T), intent(in) :: gamAppMpi
        integer, intent(in) :: recvFromRank,iData,jData,kData
        logical, intent(in) :: periodicI, periodicJ, periodicK
        type(MPI_Datatype), intent(out) :: dType
        integer, intent(out) :: offset
        logical, intent(in) :: doFace, doEdge, doCorner

        integer :: tgtRank, calcOffset, ierr, dataSize

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        ! this function calculates the 3D offset and datatype for this rank to receive data in the
        !   iData,jData,kData position from recvFromRank (if that rank sends that data to this rank)
        if((Grid%Ri+iData >= 0 .and. Grid%Ri+iData < Grid%NumRi) .or. periodicI) then
            ! talking to a neighbor, or wrapping around periodic MPI boundaries
            tgtRank = modulo(Grid%Ri+iData,Grid%NumRi)*Grid%NumRk*Grid%NumRj
        else
            ! talking to a non-periodic boundary, so this is my own I rank
            tgtRank = Grid%Ri*Grid%NumRk*Grid%NumRj
        endif
        if((Grid%Rj+jData >= 0 .and. Grid%Rj+jData < Grid%NumRj) .or. periodicJ) then
            ! talking to a neighbor, or wrapping around periodic MPI boundaries
            tgtRank = tgtRank + modulo(Grid%Rj+jData,Grid%NumRj)*Grid%NumRk
        else
            ! talking to a non-periodic boundary, so this is my own J rank
            tgtRank = tgtRank + Grid%Rj*Grid%NumRk
        endif
        if((Grid%Rk+kData >= 0 .and. Grid%Rk+kData < Grid%NumRk) .or. periodicK) then
            ! talking to a neighbor, or wrapping around periodic MPI boundaries
            tgtRank = tgtRank + modulo(Grid%Rk+kData,Grid%NumRk)
        else
            ! talking to a non-periodic boundary, so this is my own K rank
            tgtRank = tgtRank + Grid%Rk
        endif

        if(tgtRank /= recvFromRank) then
            ! this rank does not receive this data from the specified other rank
            dType = MPI_DATATYPE_NULL
            offset = -1
            return
        endif

        ! assemble the datatype, and calculate the array offset
        call calcDatatypeCC(gamAppMpi,iData,jData,kData,dType,doFace,doEdge,doCorner)
        if(dType == MPI_DATATYPE_NULL) then
            offset = -1
            return
        endif

        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry
        SELECT CASE (iData)
            CASE (-1)
                calcOffset = 0
            CASE (0)
                calcOffset = Model%nG
            CASE (1)
                calcOffset = Model%nG+Grid%Nip
            CASE DEFAULT
                write (*,*) 'Unrecognized iData type in calcRecvDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (jData)
            CASE (-1)
                ! no change to calcOffset
            CASE (0)
                calcOffset = calcOffset + Model%nG*Grid%Ni
            CASE (1)
                calcOffset = calcOffset + (Model%nG + Grid%Njp)*Grid%Ni
            CASE DEFAULT
                write (*,*) 'Unrecognized jData type in calcRecvDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (kData)
            CASE (-1)
                ! no change to calcOffset
            CASE (0)
                calcOffset = calcOffset + Model%nG*Grid%Ni*Grid%Nj
            CASE (1)
                calcOffset = calcOffset + (Model%nG + Grid%Nkp)*Grid%Ni*Grid%Nj
            CASE DEFAULT
                write (*,*) 'Unrecognized kData type in calcRecvDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        ! return calculated results
        offset = calcOffset*dataSize

        end associate

    end subroutine calcRecvDatatypeOffsetCC

    subroutine calcSendDatatypeOffsetCC(gamAppMpi,sendToRank,iData,jData,kData,&
                                  periodicI,periodicJ,periodicK,dType,offset,&
                                  doFace,doEdge,doCorner)
        type(gamAppMpi_T), intent(in) :: gamAppMpi
        integer, intent(in) :: sendToRank,iData,jData,kData
        logical, intent(in) :: periodicI, periodicJ, periodicK
        type(MPI_Datatype), intent(out) :: dType
        integer, intent(out) :: offset
        logical, intent(in) :: doFace, doEdge, doCorner

        integer :: myRank, tempRank, sendToI, sendToJ, sendToK
        logical :: wrapI, wrapJ, wrapK
        integer :: tgtRank, calcOffset, ierr, dataSize

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        ! this function calculates the 3D offset and datatype for this rank to SEND data to the
        !   iData,jData,kData position in sendToRank (if that rank receives that data from this rank)

        ! calculate which rank the target rank receives the specified data from
        tempRank = sendToRank
        sendToK = modulo(tempRank, Grid%NumRk)
        tempRank = (tempRank-sendToK)/Grid%NumRk
        sendToJ = modulo(tempRank, Grid%NumRj)
        tempRank = (tempRank-sendToJ)/Grid%NumRj
        sendToI = tempRank
        if((sendToI+iData >= 0 .and. sendToI+iData < Grid%NumRi) .or. periodicI) then
            tgtRank = modulo(sendToI+iData,Grid%NumRi)*Grid%NumRk*Grid%NumRj
            wrapI = .true.
        else
            tgtRank = sendToI*Grid%NumRk*Grid%NumRj
            wrapI = .false.
        endif
        if((sendToJ+jData >= 0 .and. sendToJ+jData < Grid%NumRj) .or. periodicJ) then
            tgtRank = tgtRank + modulo(sendToJ+jData,Grid%NumRj)*Grid%NumRk
            wrapJ = .true.
        else
            tgtRank = tgtRank + sendToJ*Grid%NumRk
            wrapJ = .false.
        endif
        if((sendToK+kData >= 0 .and. sendToK+kData < Grid%NumRk) .or. periodicK) then
            tgtRank = tgtRank + modulo(sendToK+kData,Grid%NumRk)
            wrapK = .true.
        else
            tgtRank = tgtRank + sendToK
            wrapK = .false.
        endif

        call mpi_comm_rank(gamAppMpi%gamMpiComm, myRank, ierr)

        if(tgtRank /= myRank) then
            ! this rank does not send data to the specified position on the other rank
            dType = MPI_DATATYPE_NULL
            offset = -1
            return
        endif

        ! assemble the datatype, and calculate the array offset
        call calcDatatypeCC(gamAppMpi,iData,jData,kData,dType,doFace,doEdge,doCorner)
        if(dType == MPI_DATATYPE_NULL) then
            offset = -1
            return
        endif

        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry
        SELECT CASE (iData)
            CASE (-1)
                if(wrapI) then
                    ! sending upper physical data to lower receiver
                    calcOffset = Grid%Nip
                else
                    ! sending lower ghost data to lower receiver
                    calcOffset = 0
                endif
            CASE (0)
                ! this is always full physical data
                calcOffset = Model%nG
            CASE (1)
                if(wrapI) then
                    ! sending lower physical data to upper receiver
                    calcOffset = Model%nG
                else
                    ! sending upper ghosts to upper receiver
                    calcOffset = Model%nG + Grid%Nip
                endif
            CASE DEFAULT
                write (*,*) 'Unrecognized iData type in calcSendDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (jData)
            CASE (-1)
                if(wrapJ) then
                    calcOffset = calcOffset + Grid%Njp*Grid%Ni
                else
                    ! no change to calcOffset
                endif
            CASE (0)
                calcOffset = calcOffset + Model%nG*Grid%Ni
            CASE (1)
                if(wrapJ) then
                    calcOffset = calcOffset + Model%nG*Grid%Ni
                else
                    calcOffset = calcOffset + (Model%nG+Grid%Njp)*Grid%Ni
                endif
            CASE DEFAULT
                write (*,*) 'Unrecognized jData type in calcSendDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (kData)
            CASE (-1)
                if(wrapK) then
                    calcOffset = calcOffset + Grid%Nkp*Grid%Ni*Grid%Nj
                else
                    ! no change to calcOffset
                endif
            CASE (0)
                calcOffset = calcOffset + Model%nG*Grid%Ni*Grid%Nj
            CASE (1)
                if(wrapK) then
                    calcOffset = calcOffset + Model%nG*Grid%Ni*Grid%Nj
                else
                    calcOffset = calcOffset + (Model%nG+Grid%Nkp)*Grid%Ni*Grid%Nj
                endif
            CASE DEFAULT
                write (*,*) 'Unrecognized kData type in calcSendDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        ! return calculated results
        offset = calcOffset*dataSize

        end associate

    end subroutine calcSendDatatypeOffsetCC

    subroutine calcDatatypeCC(gamAppMpi,iData,jData,kData,dType,&
                              doFace,doEdge,doCorner)
        type(gamAppMpi_T), intent(in) :: gamAppMpi
        integer, intent(in) :: iData, jData, kData
        type(MPI_Datatype), intent(out) :: dType
        logical, intent(in) :: doFace, doEdge, doCorner

        type(MPI_Datatype) :: dType1D,dType2D,dType3D
        integer :: dataSize,ierr,dataSum

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        dataSum = abs(iData)+abs(jData)+abs(kData)
        if(dataSum == 1 .and. .not. doFace) then
            dType = MPI_DATATYPE_NULL
            return
        elseif(dataSum == 2 .and. .not. doEdge) then
            dType = MPI_DATATYPE_NULL
            return
        elseif(dataSum == 3 .and. .not. doCorner) then
            dType = MPI_DATATYPE_NULL
            return
        endif

        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry
        SELECT CASE (iData)
            CASE (-1)
                call mpi_type_contiguous(Model%nG, MPI_MYFLOAT, dType1D, ierr)
            CASE (0)
                call mpi_type_contiguous(Grid%Nip, MPI_MYFLOAT, dType1D, ierr)
            CASE (1)
                call mpi_type_contiguous(Model%nG, MPI_MYFLOAT, dType1D, ierr)
            CASE DEFAULT
                write (*,*) 'Unrecognized iData type in calcDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (jData)
            CASE (-1)
                call mpi_type_hvector(Model%nG, 1, Grid%Ni*dataSize, dType1D, dType2D, ierr)
            CASE (0)
                call mpi_type_hvector(Grid%Njp, 1, Grid%Ni*dataSize, dType1D, dType2D, ierr)
            CASE (1)
                call mpi_type_hvector(Model%nG, 1, Grid%Ni*dataSize, dType1D, dType2D, ierr)
            CASE DEFAULT
                write (*,*) 'Unrecognized jData type in calcDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (kData)
            CASE (-1)
                call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, dType2D, dType3D, ierr)
            CASE (0)
                call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*dataSize, dType2D, dType3D, ierr)
            CASE (1)
                call mpi_type_hvector(Model%nG, 1, Grid%Ni*Grid%Nj*dataSize, dType2D, dType3D, ierr)
            CASE DEFAULT
                write (*,*) 'Unrecognized kData type in calcDatatypeCC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        end associate

        dType = dType3D

    end subroutine calcDatatypeCC

    subroutine calcRecvDatatypeOffsetFC(gamAppMpi,recvFromRank,iData,jData,kData,&
                                  periodicI,periodicJ,periodicK,dType,offset,&
                                  doFace,doEdge,doCorner)
        type(gamAppMpi_T), intent(in) :: gamAppMpi
        integer, intent(in) :: recvFromRank,iData,jData,kData
        logical, intent(in) :: periodicI, periodicJ, periodicK
        type(MPI_Datatype), intent(out) :: dType
        integer, intent(out) :: offset
        logical, intent(in) :: doFace, doEdge, doCorner

        integer :: tgtRank, calcOffset, ierr, dataSize

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        ! this function calculates the 3D offset and datatype for this rank to receive data in the
        !   iData,jData,kData position from recvFromRank (if that rank sends that data to this rank)
        if((Grid%Ri+iData >= 0 .and. Grid%Ri+iData < Grid%NumRi) .or. periodicI) then
            ! talking to a neighbor, or wrapping around periodic MPI boundaries
            tgtRank = modulo(Grid%Ri+iData,Grid%NumRi)*Grid%NumRk*Grid%NumRj
        else
            ! talking to a non-periodic boundary, so this is my own I rank
            tgtRank = Grid%Ri*Grid%NumRk*Grid%NumRj
        endif
        if((Grid%Rj+jData >= 0 .and. Grid%Rj+jData < Grid%NumRj) .or. periodicJ) then
            ! talking to a neighbor, or wrapping around periodic MPI boundaries
            tgtRank = tgtRank + modulo(Grid%Rj+jData,Grid%NumRj)*Grid%NumRk
        else
            ! talking to a non-periodic boundary, so this is my own J rank
            tgtRank = tgtRank + Grid%Rj*Grid%NumRk
        endif
        if((Grid%Rk+kData >= 0 .and. Grid%Rk+kData < Grid%NumRk) .or. periodicK) then
            ! talking to a neighbor, or wrapping around periodic MPI boundaries
            tgtRank = tgtRank + modulo(Grid%Rk+kData,Grid%NumRk)
        else
            ! talking to a non-periodic boundary, so this is my own K rank
            tgtRank = tgtRank + Grid%Rk
        endif

        if(tgtRank /= recvFromRank) then
            ! this rank does not receive this data from the specified other rank
            dType = MPI_DATATYPE_NULL
            offset = -1
            return
        endif

        ! assemble the datatype, and calculate the array offset
        call calcDatatypeFC(gamAppMpi,iData,jData,kData,dType,doFace,doEdge,doCorner)
        if(dType == MPI_DATATYPE_NULL) then
            offset = -1
            return
        endif

        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry
        SELECT CASE (iData)
            CASE (-1)
                ! overall offset of the final struct from the start of the data array
                calcOffset = 0
            CASE (0)
                calcOffset = Model%nG
            CASE (1)
                calcOffset = Model%nG+Grid%Nip
            CASE DEFAULT
                write (*,*) 'Unrecognized iData type in calcRecvDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (jData)
            CASE (-1)
                ! no change to overall offset, calcOffset
            CASE (0)
                calcOffset = calcOffset + Model%nG*(Grid%Ni+1)
            CASE (1)
                calcOffset = calcOffset + (Model%nG+Grid%Njp)*(Grid%Ni+1)
            CASE DEFAULT
                write (*,*) 'Unrecognized jData type in calcRecvDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (kData)
            CASE (-1)
                ! no change to overall offset, calcOffset
            CASE (0)
                calcOffset = calcOffset + Model%nG*(Grid%Ni+1)*(Grid%Nj+1)
            CASE (1)
                calcOffset = calcOffset + (Model%nG+Grid%Nkp)*(Grid%Ni+1)*(Grid%Nj+1)
            CASE DEFAULT
                write (*,*) 'Unrecognized kData type in calcRecvDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        ! return calculated results
        offset = calcOffset*dataSize

        end associate

    end subroutine calcRecvDatatypeOffsetFC

    subroutine calcSendDatatypeOffsetFC(gamAppMpi,sendToRank,iData,jData,kData,&
                                  periodicI,periodicJ,periodicK,dType,offset,&
                                  doFace,doEdge,doCorner)
        type(gamAppMpi_T), intent(in) :: gamAppMpi
        integer, intent(in) :: sendToRank,iData,jData,kData
        logical, intent(in) :: periodicI, periodicJ, periodicK
        type(MPI_Datatype), intent(out) :: dType
        integer, intent(out) :: offset
        logical, intent(in) :: doFace, doEdge, doCorner

        integer :: myRank, tempRank, sendToI, sendToJ, sendToK
        logical :: wrapI, wrapJ, wrapK
        integer :: tgtRank, calcOffset, ierr, dataSize

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        ! this function calculates the 3D offset and datatype for this rank to SEND data to the
        !   iData,jData,kData position in sendToRank (if that rank receives that data from this rank)

        ! calculate which rank the target rank receives the specified data from
        tempRank = sendToRank
        sendToK = modulo(tempRank, Grid%NumRk)
        tempRank = (tempRank-sendToK)/Grid%NumRk
        sendToJ = modulo(tempRank, Grid%NumRj)
        tempRank = (tempRank-sendToJ)/Grid%NumRj
        sendToI = tempRank
        if((sendToI+iData >= 0 .and. sendToI+iData < Grid%NumRi) .or. periodicI) then
            tgtRank = modulo(sendToI+iData,Grid%NumRi)*Grid%NumRk*Grid%NumRj
            wrapI = .true.
        else
            tgtRank = sendToI*Grid%NumRk*Grid%NumRj
            wrapI = .false.
        endif
        if((sendToJ+jData >= 0 .and. sendToJ+jData < Grid%NumRj) .or. periodicJ) then
            tgtRank = tgtRank + modulo(sendToJ+jData,Grid%NumRj)*Grid%NumRk
            wrapJ = .true.
        else
            tgtRank = tgtRank + sendToJ*Grid%NumRk
            wrapJ = .false.
        endif
        if((sendToK+kData >= 0 .and. sendToK+kData < Grid%NumRk) .or. periodicK) then
            tgtRank = tgtRank + modulo(sendToK+kData,Grid%NumRk)
            wrapK = .true.
        else
            tgtRank = tgtRank + sendToK
            wrapK = .false.
        endif

        call mpi_comm_rank(gamAppMpi%gamMpiComm, myRank, ierr)

        if(tgtRank /= myRank) then
            ! this rank does not send data to the specified position on the other rank
            dType = MPI_DATATYPE_NULL
            offset = -1
            return
        endif

        ! assemble the datatype, and calculate the array offset
        call calcDatatypeFC(gamAppMpi,iData,jData,kData,dType,doFace,doEdge,doCorner)
        if(dType == MPI_DATATYPE_NULL) then
            offset = -1
            return
        endif

        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry
        SELECT CASE (iData)
            CASE (-1)
                if(wrapI) then
                    ! sending upper physical data to lower receiver
                    calcOffset = Grid%Nip
                else
                    ! sending lower ghost data to lower receiver
                    calcOffset = 0
                endif
            CASE (0)
                ! this is always full physical data
                calcOffset = Model%nG
            CASE (1)
                if(wrapI) then
                    ! sending lower physical data to upper receiver
                    calcOffset = Model%nG
                else
                    ! sending upper ghosts to upper receiver
                    calcOffset = Model%nG + Grid%Nip
                endif
            CASE DEFAULT
                write (*,*) 'Unrecognized iData type in calcSendDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (jData)
            CASE (-1)
                if(wrapJ) then
                    calcOffset = calcOffset + Grid%Njp*(Grid%Ni+1)
                else
                    ! no change to calcOffset
                endif
            CASE (0)
                calcOffset = calcOffset + Model%nG*(Grid%Ni+1)
            CASE (1)
                if(wrapJ) then
                    calcOffset = calcOffset + Model%nG*(Grid%Ni+1)
                else
                    calcOffset = calcOffset + (Model%nG+Grid%Njp)*(Grid%Ni+1)
                endif
            CASE DEFAULT
                write (*,*) 'Unrecognized jData type in calcSendDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (kData)
            CASE (-1)
                if(wrapK) then
                    calcOffset = calcOffset + Grid%Nkp*(Grid%Ni+1)*(Grid%Nj+1)
                else
                    ! no change to calcOffset
                endif
            CASE (0)
                calcOffset = calcOffset + Model%nG*(Grid%Ni+1)*(Grid%Nj+1)
            CASE (1)
                if(wrapK) then
                    calcOffset = calcOffset + Model%nG*(Grid%Ni+1)*(Grid%Nj+1)
                else
                    calcOffset = calcOffset + (Model%nG+Grid%Nkp)*(Grid%Ni+1)*(Grid%Nj+1)
                endif
            CASE DEFAULT
                write (*,*) 'Unrecognized kData type in calcSendDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        ! return calculated results
        offset = calcOffset*dataSize

        end associate

    end subroutine calcSendDatatypeOffsetFC

    subroutine calcDatatypeFC(gamAppMpi,iData,jData,kData,dType,doFace,doEdge,doCorner)
        type(gamAppMpi_T), intent(in) :: gamAppMpi
        integer, intent(in) :: iData, jData, kData
        type(MPI_Datatype), intent(out) :: dType
        logical, intent(in) :: doFace, doEdge, doCorner

        type(MPI_Datatype) :: dType1DI,dType1DJ,dType1DK,dType2DI,dType2DJ,dType2DK,dType3DI,dType3DJ,dType3DK
        integer :: offsetI, offsetJ, offsetK, dataSum, ierr, dataSize
        logical :: anyMaxDim,sendSharedFace
        integer :: countArray(3)
        integer(kind=MPI_BASE_MYADDR) :: offsetArray(3)

        associate(Grid=>gamAppMpi%Grid,Model=>gamAppMpi%Model)

        ! assemble the datatype, and calculate the array offset
        dataSum = abs(iData)+abs(jData)+abs(kData)
        anyMaxDim = iData==1 .or. jData==1 .or. kData==1
        ! due to face ownership "corners" share data with "edges", and
        !  "edges" share data with "faces", so special care must be taken
        if((dataSum == 2 .and. .not. doEdge .and. anyMaxDim .and. doFace) .or. &
           (dataSum == 3 .and. .not. doCorner .and. anyMaxDim .and. doEdge)) then
            sendSharedFace = .true.
        else
            sendSharedFace = .false.
        endif
        if(dataSum == 1 .and. .not. doFace) then
            dType = MPI_DATATYPE_NULL
            return
        elseif(dataSum == 2 .and. .not. doEdge .and. .not. sendSharedFace) then
            dType = MPI_DATATYPE_NULL
            return
        elseif(dataSum == 3 .and. .not. doCorner .and. .not. sendSharedFace) then
            dType = MPI_DATATYPE_NULL
            return
        endif

        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry
        SELECT CASE (iData)
            CASE (-1)
                ! receive lower I face and interior but not upper
                offsetI = 0
                call mpi_type_contiguous(Model%nG, MPI_MYFLOAT, dType1DI, ierr)
                call mpi_type_contiguous(Model%nG, MPI_MYFLOAT, dType1DJ, ierr)
                call mpi_type_contiguous(Model%nG, MPI_MYFLOAT, dType1DK, ierr)
            CASE (0)
                ! receive lower I face and interior but not upper
                offsetI = 0
                call mpi_type_contiguous(Grid%Nip, MPI_MYFLOAT, dType1DI, ierr)
                call mpi_type_contiguous(Grid%Nip, MPI_MYFLOAT, dType1DJ, ierr)
                call mpi_type_contiguous(Grid%Nip, MPI_MYFLOAT, dType1DK, ierr)
            CASE (1)
                ! overwrite lower, interior, and upper I faces
                offsetI = 0
                if(sendSharedFace) then
                    call mpi_type_contiguous(1, MPI_MYFLOAT, dType1DI, ierr)
                else
                    call mpi_type_contiguous(Model%nG+1, MPI_MYFLOAT, dType1DI, ierr)
                endif
                call mpi_type_contiguous(Model%nG,   MPI_MYFLOAT, dType1DJ, ierr)
                call mpi_type_contiguous(Model%nG,   MPI_MYFLOAT, dType1DK, ierr)
            CASE DEFAULT
                write (*,*) 'Unrecognized iData type in calcDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (jData)
            CASE (-1)
                offsetJ = 0
                call mpi_type_hvector(Model%nG, 1, dataSize*(Grid%Ni+1), dType1DI, dType2DI, ierr)
                call mpi_type_hvector(Model%nG, 1, dataSize*(Grid%Ni+1), dType1DJ, dType2DJ, ierr)
                call mpi_type_hvector(Model%nG, 1, dataSize*(Grid%Ni+1), dType1DK, dType2DK, ierr)
            CASE (0)
                offsetJ = 0
                call mpi_type_hvector(Grid%Njp, 1, dataSize*(Grid%Ni+1), dType1DI, dType2DI, ierr)
                call mpi_type_hvector(Grid%Njp, 1, dataSize*(Grid%Ni+1), dType1DJ, dType2DJ, ierr)
                call mpi_type_hvector(Grid%Njp, 1, dataSize*(Grid%Ni+1), dType1DK, dType2DK, ierr)
            CASE (1)
                offsetJ = 0
                call mpi_type_hvector(Model%nG,   1, dataSize*(Grid%Ni+1), dType1DI, dType2DI, ierr)
                if(sendSharedFace) then
                    call mpi_type_hvector(1, 1, dataSize*(Grid%Ni+1), dType1DJ, dType2DJ, ierr)
                else
                    call mpi_type_hvector(Model%nG+1, 1, dataSize*(Grid%Ni+1), dType1DJ, dType2DJ, ierr)
                endif
                call mpi_type_hvector(Model%nG,   1, dataSize*(Grid%Ni+1), dType1DK, dType2DK, ierr)
            CASE DEFAULT
                write (*,*) 'Unrecognized jData type in calcDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        SELECT CASE (kData)
            CASE (-1)
                offsetK = 0
                call mpi_type_hvector(Model%nG, 1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DI, dType3DI, ierr)
                call mpi_type_hvector(Model%nG, 1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DJ, dType3DJ, ierr)
                call mpi_type_hvector(Model%nG, 1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DK, dType3DK, ierr)
            CASE (0)
                offsetK = 0
                call mpi_type_hvector(Grid%Nkp, 1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DI, dType3DI, ierr)
                call mpi_type_hvector(Grid%Nkp, 1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DJ, dType3DJ, ierr)
                call mpi_type_hvector(Grid%Nkp, 1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DK, dType3DK, ierr)
            CASE (1)
                offsetK = 0
                call mpi_type_hvector(Model%nG,   1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DI, dType3DI, ierr)
                call mpi_type_hvector(Model%nG,   1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DJ, dType3DJ, ierr)
                if(sendSharedFace) then
                    call mpi_type_hvector(1, 1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DK, dType3DK, ierr)
                else
                    call mpi_type_hvector(Model%nG+1, 1, dataSize*(Grid%Ni+1)*(Grid%Nj+1), dType2DK, dType3DK, ierr)
                endif
            CASE DEFAULT
                write (*,*) 'Unrecognized kData type in calcDatatypeFC'
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        ENDSELECT

        ! move J offset into the 2nd index of the 4th dimension
        offsetJ = offsetJ + dataSize*(Grid%Ni+1)*(Grid%Nj+1)*(Grid%Nk+1)
        ! move K offset into the 3rd index of the 4th dimension
        offsetK = offsetK + 2*dataSize*(Grid%Ni+1)*(Grid%Nj+1)*(Grid%Nk+1)

        ! default counts and offsets
        countArray = (/1,1,1/)
        offsetArray = (/integer(kind=MPI_BASE_MYADDR):: offsetI, offsetJ, offsetK /)
        if(sendSharedFace) then
            ! this is only sending data for shared faces with an adjacent type
            if(iData /= 1) then
                countArray(1) = 0
                offsetArray(1) = 0
                dType3DI = MPI_INTEGER
            endif
            if(jData /= 1) then
                countArray(2) = 0
                offsetArray(2) = 0
                dType3DJ = MPI_INTEGER
            endif
            if(kData /= 1) then
                countArray(3) = 0
                offsetArray(3) = 0
                dType3DK = MPI_INTEGER
            endif
        endif
        call mpi_type_create_struct(3,countArray,offsetArray, &
                                    (/ dType3DI, dType3DJ, dType3DK /), dType, ierr)

        end associate

    end subroutine calcDatatypeFC

    subroutine appendDatatype(appendType,dType,offset)
        type(MPI_Datatype), intent(inout) :: appendType
        type(MPI_Datatype), intent(in) :: dType
        integer, intent(in) :: offset

        integer :: ierr
        integer(kind=MPI_BASE_MYADDR) :: tempOffsets(2)

        if(appendType == MPI_DATATYPE_NULL) then
            ! the root datatype is empty, just add the new type and offset
            call mpi_type_hindexed(1, (/ 1 /), (/ offset /), dType, appendType, ierr)
        else
            ! the root datatype already has defined structure, so merge with the new one into a struct
            ! need to use a temporary array so that the ints are of type MPI_BASE_MYADDR
            tempOffsets = (/ 0, offset /)
            call mpi_type_create_struct(2, (/ 1, 1 /), tempOffsets, (/ appendType, dType /), appendType, ierr)
        endif

    end subroutine appendDatatype

end module gamapp_mpi

