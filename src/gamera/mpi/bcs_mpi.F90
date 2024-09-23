!MPI specific boundary conditions
module bcs_mpi
    use gamtypes
    use gambctypes
    use gamutils
    use math
    use gridutils
    use ringutils

    implicit none

    type, extends(multiBC_T) :: mpiNullBc_T
        contains
        procedure :: doBC => mpi_null_bc
    end type mpiNullBc_T

    contains

!--------------------------------------------    
    !do-nothing BC just so that we can tell that this BC is handled by the MPI halo update
    subroutine mpi_null_bc(bc,Model,Grid,State)
        class(mpiNullBc_T), intent(inout) :: bc
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(in) :: Grid
        class(State_T), intent(inout) :: State

        ! do nothing

    end subroutine mpi_null_bc

end module bcs_mpi

