!Routines for OpenMP assisted array utilities
module arrayutil

    use kdefs

    implicit none

    ! interface to fill arrays of all sizes, multiple types
    interface fillArray
        module procedure fillArray1D_R,fillArray1D_I,fillArray2D_R,fillArray2D_I,&
            fillArray3D_R,fillArray3D_I,fillArray4D_R,fillArray4D_I,fillArray5D_R,fillArray5D_I
    end interface

    contains


    subroutine fillArray1D_R(array, val)
        real(rp), dimension(:), intent(inout) :: array
        real(rp), intent(in) :: val

        integer :: lbds(1),ubds(1),i

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i)
        do i=lbds(1),ubds(1)
            array(i) = val
        enddo

#else
        array = val
#endif

    end subroutine fillArray1D_R

    subroutine fillArray1D_I(array, val)
        integer, dimension(:), intent(inout) :: array
        integer, intent(in) :: val

        integer :: lbds(1),ubds(1),i

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i)
        do i=lbds(1),ubds(1)
            array(i) = val
        enddo

#else
        array = val
#endif

    end subroutine fillArray1D_I

    subroutine fillArray2D_R(array, val)
        real(rp), dimension(:,:), intent(inout) :: array
        real(rp), intent(in) :: val

        integer :: lbds(2),ubds(2),i,j

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j)
        do j=lbds(2),ubds(2)
            do i=lbds(1),ubds(1)
                array(i,j) = val
            enddo
        enddo

#else
        array = val
#endif

    end subroutine fillArray2D_R

    subroutine fillArray2D_I(array, val)
        integer, dimension(:,:), intent(inout) :: array
        integer, intent(in) :: val

        integer :: lbds(2),ubds(2),i,j

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j)
        do j=lbds(2),ubds(2)
            do i=lbds(1),ubds(1)
                array(i,j) = val
            enddo
        enddo

#else
        array = val
#endif

    end subroutine fillArray2D_I

    subroutine fillArray3D_R(array, val)
        real(rp), dimension(:,:,:), intent(inout) :: array
        real(rp), intent(in) :: val

        integer :: lbds(3),ubds(3),i,j,k

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=lbds(3),ubds(3)
            do j=lbds(2),ubds(2)
                do i=lbds(1),ubds(1)
                    array(i,j,k) = val
                enddo
            enddo
        enddo

#else
        array = val
#endif

    end subroutine fillArray3D_R

    subroutine fillArray3D_I(array, val)
        integer, dimension(:,:,:), intent(inout) :: array
        integer, intent(in) :: val

        integer :: lbds(3),ubds(3),i,j,k

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=lbds(3),ubds(3)
            do j=lbds(2),ubds(2)
                do i=lbds(1),ubds(1)
                    array(i,j,k) = val
                enddo
            enddo
        enddo

#else
        array = val
#endif

    end subroutine fillArray3D_I

    subroutine fillArray4D_R(array, val)
        real(rp), dimension(:,:,:,:), intent(inout) :: array
        real(rp), intent(in) :: val

        integer :: lbds(4),ubds(4),i,j,k

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=lbds(3),ubds(3)
            do j=lbds(2),ubds(2)
                do i=lbds(1),ubds(1)
                    array(i,j,k,:) = val
                enddo
            enddo
        enddo

#else
        array = val
#endif

    end subroutine fillArray4D_R

    subroutine fillArray4D_I(array, val)
        integer, dimension(:,:,:,:), intent(inout) :: array
        integer, intent(in) :: val

        integer :: lbds(4),ubds(4),i,j,k

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=lbds(3),ubds(3)
            do j=lbds(2),ubds(2)
                do i=lbds(1),ubds(1)
                    array(i,j,k,:) = val
                enddo
            enddo
        enddo

#else
        array = val
#endif

    end subroutine fillArray4D_I

    subroutine fillArray5D_R(array, val)
        real(rp), dimension(:,:,:,:,:), intent(inout) :: array
        real(rp), intent(in) :: val

        integer :: lbds(5),ubds(5),i,j,k

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=lbds(3),ubds(3)
            do j=lbds(2),ubds(2)
                do i=lbds(1),ubds(1)
                    array(i,j,k,:,:) = val
                enddo
            enddo
        enddo

#else
        array = val
#endif

    end subroutine fillArray5D_R

    subroutine fillArray5D_I(array, val)
        integer, dimension(:,:,:,:,:), intent(inout) :: array
        integer, intent(in) :: val

        integer :: lbds(5),ubds(5),i,j,k

        lbds = lbound(array)
        ubds = ubound(array)

#ifdef _OPENMP
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=lbds(3),ubds(3)
            do j=lbds(2),ubds(2)
                do i=lbds(1),ubds(1)
                    array(i,j,k,:,:) = val
                enddo
            enddo
        enddo

#else
        array = val
#endif

    end subroutine fillArray5D_I


end module arrayutil

