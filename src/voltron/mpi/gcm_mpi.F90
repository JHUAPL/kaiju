module gcm_mpi
  use gcminterp
  use mpi_f08
  use mpidefs

  implicit none

  interface exportgcm
    module procedure exportgcmmpi
  end interface
  interface importgcm
    module procedure importgcmmpi
  end interface

contains

    subroutine init_gcm_mpi(gcm,ion,isRestart)
      type(gcm_T),intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      logical, intent(in) :: isRestart

      write(*,*) "start init_gcm"

      if (.not. allocated(gcm%GEO%outlist)) allocate(gcm%GEO%outlist(gcm%GEO%mix2gcm_nvar))
      if (.not. allocated(gcm%APEX%outlist)) allocate(gcm%APEX%outlist(gcm%APEX%mix2gcm_nvar))
      gcm%APEX%outlist(1) = POT
      gcm%APEX%outlist(2) = NUM_FLUX
      gcm%GEO%outlist(1) = AVG_ENG
      gcm%GEO%outlist(2) = NUM_FLUX

      if ( .not. allocated(gcm%GEO%inlist)) allocate(gcm%GEO%inlist(gcm%GEO%gcm2mix_nvar))
      if ( .not. allocated(gcm%APEX%inlist)) allocate(gcm%APEX%inlist(gcm%APEX%gcm2mix_nvar))
      gcm%APEX%inlist(1) = SIGMAP
      gcm%APEX%inlist(2) = SIGMAH

      gcm%isRestart = isRestart

      call init_gcm_grid(gcm,ion)
      
    end subroutine init_gcm_mpi

    subroutine init_gcm_mix_mpi(gcm,gcmCplComm,gcmCplRank)
      type(gcm_T), intent(inout) :: gcm
      type(MPI_Comm) :: gcmCplComm
      integer, intent(in) :: gcmCplRank

      integer :: g,i

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Find grid dimension
 
      write(*,*) "MIXCPL Waiting for Grid"

      i = 0
      do g=1,2
        select case(g)
        case (1)
          call init_gcm_mpi_import(gcm%GEO,gcmCplComm,gcmCplRank,i)
          write(*,*) "MIXCPL LOAD GEO"
        case (2)
          call init_gcm_mpi_import(gcm%APEX,gcmCplComm,gcmCplRank,i)
          write(*,*) "MIXCPL LOAD APEX"
        end select

      end do
    end subroutine init_gcm_mix_mpi

    subroutine init_gcm_mpi_import(gcm,gcmCplComm,gcmCplRank,i)
        type(gcm_grid_T),intent(inout):: gcm
        type(MPI_Comm) :: gcmCplComm
       integer, intent(in) :: gcmCplRank

        integer :: i,nlon,nlat,Nh,ierr

        Nh = GCMhemispheres

        i = i + 1
        call mpi_recv(gcm%nlat, 1, MPI_INTEGER, gcmCplRank, (tgcmId+voltId)*100+i, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT NLAT: ",gcm%nlat

        i = i + 1
        call mpi_recv(gcm%nlon, 1, MPI_INTEGER, gcmCplRank, (tgcmId+voltId)*100+i, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT NLON: ",gcm%nlon


        if (.not.allocated(gcm%inlat)) allocate(gcm%inlat(gcm%nlat))
        if (.not.allocated(gcm%inlon)) allocate(gcm%inlon(gcm%nlon))

        i = i + 1
        call mpi_recv(gcm%inlat, gcm%nlat, MPI_DOUBLE_PRECISION, gcmCplRank, (tgcmId+voltId)*100+i, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT LAT: ",gcm%inlat
        i = i + 1
        call mpi_recv(gcm%inlon, gcm%nlon, MPI_DOUBLE_PRECISION, gcmCplRank, (tgcmId+voltId)*100+i, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT LON: ",gcm%inlon
        write(*,*) "MIXCPL GOT GRID INFO: ",gcm%nlat,gcm% nlon

    end subroutine init_gcm_mpi_import

    subroutine coupleGCM2MIX(gcm,ion,mjd,time,gcmCplComm,myRank)
        type(mixIon_T),dimension(:),intent(inout) :: ion
        type(gcm_T), intent(inout) :: gcm
        real(rp), intent(in) :: time, mjd
        type(MPI_Comm), optional :: gcmCplComm
        integer, optional, intent(in) :: myRank

        ! maybe create the SM and GEO list of points here
        ! since the destination grid does not need to be structured
        ! can do a simple loop over all grid points and transform
        ! transform all remix points to GEO
        ! transform all gcm points to SM

        !Skip MPI exchange on first restart
        if (present(gcmCplComm) .and. gcm%isRestart) then
            !gcm%isRestart = .False.
            return
        endif
        !Must MIX export first.  TIEGCM will also import first.
        call Tic("Export")
        if (present(gcmCplComm)) then
            call exportgcm(ion,gcm,mjd,time,gcmCplComm,myRank)
        else
            write(*,*) "Are we trying to Export to Couple GCM?"
        endif
        call Toc("Export")

        !Import gcm data
        call Tic("Import")
        if (present(gcmCplComm)) then
            call importgcm(gcm, gcmCplComm,myRank)
        else
            write(*,*) "Are we trying to Import to Couple GCM?"
        endif
        call process_gcmimport(gcm,ion)
        call Toc("Import")

        if (gcm%isRestart) gcm%isRestart=.false.
        ! We have a separate internal counter for coupling here.
        ! This may be used later on for WACCM-X coupling which is desync from remix coupling time
        ! TIEGCM coupling time is also 5s while WACCM-X will couple at 1 min default
        gcm%cplStep = gcm%cplStep + 1
        
    end subroutine coupleGCM2MIX


    subroutine importgcmmpi(gcm, gcmCplComm,gcmCplRank)
        type(gcm_T), intent(inout) :: gcm
        type(MPI_Comm) :: gcmCplComm
        integer, intent(in) :: gcmCplRank

        integer :: v

        call Tic("MpiExchange")
        write(*,*) "MIX: IMPORT GCM"

        if (gcmCplComm /= MPI_COMM_NULL) then
            call Tic("Passing")
            do v=1,2
            select case(v)
            case (1)
                !write(*,*) "Import GEO stuff"
                call import_gcm_per_grid(gcm%GEO,gcmCplComm,gcmCplRank)
            case (2)
                !write(*,*) "Import APEX stuff"
                call import_gcm_per_grid(gcm%APEX,gcmCplComm,gcmCplRank)
            end select

            enddo
            call Toc("Passing")

        endif
        call Toc("MpiExchange")

    end subroutine importgcmmpi

    subroutine import_gcm_per_grid(gcm,gcmCplComm,gcmCplRank)
        type(gcm_grid_T) :: gcm
        type(MPI_Comm) :: gcmCplComm
        integer, intent(in) :: gcmCplRank

        integer :: ierr,length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! Skip if not importing on this grid type
        if ( gcm%gcm2mix_nvar .eq. 0 ) return

        if (.not. allocated(gcm%invar2d)) allocate(gcm%invar2d(gcm%nlat,gcm%nlon,gcm%gcm2mix_nvar))

        call mpi_recv(gcm%invar2d, gcm%nlat*gcm%nlon*gcm%gcm2mix_nvar, MPI_DOUBLE_PRECISION, gcmCplRank, (tgcmId+voltId)*100, gcmCplComm, MPI_STATUS_IGNORE,ierr)

        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
    end subroutine import_gcm_per_grid

    subroutine exportgcmmpi(ion,gcm,mjd,time,gcmCplComm,gcmCplRank)
        type(gcm_T), intent(inout) :: gcm
        type(mixIon_T), dimension(:),intent(inout) :: ion
        real(rp), intent(in) :: time, mjd
        type(MPI_Comm) :: gcmCplComm
        integer, intent(in) :: gcmCplRank

        integer :: g,ierr,length
        real(rp),allocatable,dimension(:,:) :: auroralbc
        character( len = MPI_MAX_ERROR_STRING) :: message

        write(*,*) "MIX: EXPORT GCM"

        ! Prepare the export data
        call Tic("Transform")
        do g = 1,2
            select case (g)
            case (1)
                call transform_gcm_export(gcm%GEO,ion,g)
            case(2)
                call transform_gcm_export(gcm%APEX,ion,g)
            end select

        end do
        call Toc("Transform")

        !! Calculate the auroral boundary in APEX coordinates
        !call calculate_auroral_boundary(gcm%APEX,auroralbc)

        call Tic("Exchange")

        if (gcmCplComm /= MPI_COMM_NULL) then
        do g = 1,2
            select case (g)
            case (1)
            call export_gcm_per_grid(gcm%GEO,gcmCplComm,gcmCplRank)
            case(2)
            call export_gcm_per_grid(gcm%APEX,gcmCplComm,gcmCplRank)
            end select
        enddo

        !call mpi_send(auroralbc,gcm%APEX%nlon*2, MPI_DOUBLE_PRECISION, gcmCplRank,(tgcmId+voltId)*100, gcmCplComm, ierr)
        !    if(ierr /= MPI_Success) then
        !        call MPI_Error_string( ierr, message, length, ierr)
        !        print *,message(1:length)
        !        call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        !    end if
        endif
        
        call Toc("Exchange")
    end subroutine exportgcmmpi

    subroutine export_gcm_per_grid(gcm,gcmCplComm,gcmCplRank)
        type(gcm_grid_T),intent(inout) :: gcm
        type(MPI_Comm) :: gcmCplComm
        integer, intent(in) :: gcmCplRank

        integer :: ierr,length
        character( len = MPI_MAX_ERROR_STRING) :: message

        if (gcm%mix2gcm_nvar .eq. 0) return

        ! Send the coupling data
        !write(*,*) " MIXCPL: ", gcmCplRank,(tgcmId+voltId)*100,gcmCplComm,gcm%nlat,gcm%nlon

        call mpi_send(gcm%outvar2d, gcm%nlat*gcm%nlon*gcm%mix2gcm_nvar, MPI_DOUBLE_PRECISION, gcmCplRank, (tgcmId+voltId)*100, gcmCplComm, ierr)

        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine export_gcm_per_grid
end module gcm_mpi
