module mpidefs

  use mpi_f08
  use kdefs, ONLY: sp,dp,rp

  implicit none

  ! datatype for the floating point used by kaiju, assigned below
  type(MPI_Datatype), public :: MPI_MYFLOAT
  type(MPI_Datatype), public :: MPI_2MYFLOAT

  ! datatype for address sizes. Adjustable because different mpi
  !   implementations don't always do it the same way
  !   address pointer for basic mpi functions, and
  !   address point for neighborhood mpi functions
#ifdef MPI_BASE_ADDR_SIZE
  integer, parameter :: MPI_BASE_MYADDR = MPI_BASE_ADDR_SIZE
#else
  integer, parameter :: MPI_BASE_MYADDR = MPI_ADDRESS_KIND ! this is the default
#endif

#ifdef MPI_AN_ADDR_SIZE
  integer, parameter :: MPI_AN_MYADDR   = MPI_AN_ADDR_SIZE
#else
  integer, parameter :: MPI_AN_MYADDR   = MPI_ADDRESS_KIND ! this is the default
#endif

!MPI Coupling Communicator Id
    integer, parameter :: voltId   = 116
    integer, parameter :: helperId = 64
    integer, parameter :: gamId    = 45
    integer, parameter :: rcmId    = 34
    integer, parameter :: tgcmId   = 57
    integer, parameter :: hidraId  = 40
    integer, parameter :: hidraNId = 54
    integer, parameter :: hidraSId = 59
    integer, parameter :: mageId   = 26

contains

  subroutine setMpiReal()

    ! make the mpi float datatype match the precision of the fortran float type being used
    MPI_MYFLOAT = MPI_DOUBLE_PRECISION
    MPI_2MYFLOAT = MPI_2DOUBLE_PRECISION
    if(rp ==  sp) then
        MPI_MYFLOAT = MPI_REAL
        MPI_2MYFLOAT = MPI_2REAL
    endif

  end subroutine setMpiReal

  ! helper functions to print info about custom MPI datatypes
  subroutine simplePrintDataType(datatype)
      type(MPI_Datatype), intent(in) :: datatype

      logical :: typeUsed
      typeUsed = printDataType(datatype)

  end subroutine

  recursive function printDataType(datatype) result(retVal)
      type(MPI_Datatype), intent(in) :: datatype
      logical :: retVal

      integer :: numInts, numAdds, numDTs, combiner, ierr, i
      integer, dimension(:), allocatable :: arrayInts
      type(MPI_Datatype), dimension(:), allocatable :: arrayDTs
      integer(kind=MPI_BASE_MYADDR), dimension(:), allocatable :: arrayAdds
      character(len=MPI_MAX_OBJECT_NAME) :: typeName
      call mpi_type_get_envelope(datatype, numInts, numAdds, numDTs, combiner, ierr)

      SELECT CASE(combiner)
          CASE (MPI_COMBINER_NAMED)
              call mpi_type_get_name(datatype, typeName, i, ierr)
              write (*,*) 'Datatype is named: ', trim(typeName)
              retVal = .false.
              RETURN
          CASE (MPI_COMBINER_STRUCT)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is struct containing ',arrayInts(1),' datatypes:'
              do i=1,arrayInts(1)
                  write (*,*) 'blocklength ',arrayInts(i+1),', displacement ',arrayAdds(i),' type:'
                  if (printDataType(arrayDTs(i))) then
                      call mpi_type_free(arrayDTs(i),ierr)
                  endif
              enddo
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE (MPI_COMBINER_HINDEXED)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is hindexed containing ',arrayInts(1),' blocks:'
              do i=1,arrayInts(1)
                  write (*,*) 'blocklength ',arrayInts(i+1),' displacement ',arrayAdds(i)
              enddo
              write (*,*) 'hindexed datatype:'
              if (printDataType(arrayDTs(1))) then
                  call mpi_type_free(arrayDTs(1), ierr)
              endif
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE (MPI_COMBINER_HVECTOR)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is hvector containing ',arrayInts(1),' blocks:'
              write (*,*) 'blocklength ',arrayInts(2),' displacement ',arrayAdds(1),' type:'
              if (printDataType(arrayDTs(1))) then
                  call mpi_type_free(arrayDTs(1), ierr)
              endif
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE (MPI_COMBINER_CONTIGUOUS)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is contiguous containing ',arrayInts(1),' repetitions of datatype:'
              if (printDataType(arrayDTs(1))) then
                  call mpi_type_free(arrayDTs(1), ierr)
              endif
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE (MPI_COMBINER_SUBARRAY)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is subarray containing ',arrayInts(1),' dimensions:'
              do i=1,arrayInts(1)
                  write (*,*) 'size ',arrayInts(1+i),' subsize ',arrayInts(arrayInts(1)+1+i),' start ',arrayInts(2*arrayInts(1)+1+i)
              enddo
              if(arrayInts(3*arrayInts(1)+2) .eq. MPI_ORDER_FORTRAN) then
                  write (*,*) 'data order is FORTRAN'
              else
                  write (*,*) 'data order is ',arrayInts(3*arrayInts(1)+2)
              endif
              write (*,*) 'subarray datatype:'
              if (printDataType(arrayDTs(1))) then
                  call mpi_type_free(arrayDTs(1), ierr)
              endif
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE DEFAULT
              write (*,*) 'Unknown combiner type in printDataType'
              retVal = .false.
              RETURN
      ENDSELECT

      retVal = .true.

  end function printDataType

end module mpidefs
