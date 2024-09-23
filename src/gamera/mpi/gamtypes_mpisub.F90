submodule (gamtypes_mpi) gamtypes_mpisub
    use gamapp_mpi
    use mpidefs
    use mpi_f08

    implicit none

    contains

    ! procedures for gamAppMpi_T
    module subroutine gamMpiInitModel(App, Xml)
        class(gamAppMpi_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call initGamera_mpi(App, Xml)

    end subroutine gamMpiInitModel

    module subroutine gamMpiWriteConsoleOutput(App)
        class(gamAppMpi_T), intent(inout) :: App

        call consoleOutput_mpi(App)

    end subroutine gamMpiWriteConsoleOutput

    module subroutine gamMpiAdvanceModel(App, dt)
        class(gamAppMpi_T), intent(inout) :: App
        real(rp), intent(in) :: dt

        real(rp) :: targetSimT 

        targetSimT = App%Model%t+dt

        ! ensure gamera is stepped at least one time
        ! fortran doesn't allow checking at the end of a while loop
        call stepGamera_mpi(App)

        do while(App%Model%t < targetSimT)
            call stepGamera_mpi(App)
        enddo

    end subroutine gamMpiAdvanceModel

end submodule

