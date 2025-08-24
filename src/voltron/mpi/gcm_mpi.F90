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

      call init_gcm_grid(gcm,ion)

      write(*,*) "done init_gcm_grid"

      if (.not. allocated(gcm%GEO%outlist)) allocate(gcm%GEO%outlist(gcm%GEO%mix2gcm_nvar))
      if (.not. allocated(gcm%APEX%outlist)) allocate(gcm%APEX%outlist(gcm%APEX%mix2gcm_nvar))
      
      if (gcm%isWACCMX) then
        gcm%APEX%outlist(1) = POT
        gcm%APEX%outlist(2) = AVG_ENG
        gcm%APEX%outlist(3) = NUM_FLUX
      else ! is TIEGCM
        gcm%APEX%outlist(1) = POT
        gcm%APEX%outlist(2) = NUM_FLUX
        gcm%GEO%outlist(1) = AVG_ENG
        gcm%GEO%outlist(2) = NUM_FLUX
      endif

      if ( .not. allocated(gcm%GEO%inlist)) allocate(gcm%GEO%inlist(gcm%GEO%gcm2mix_nvar))
      if ( .not. allocated(gcm%APEX%inlist)) allocate(gcm%APEX%inlist(gcm%APEX%gcm2mix_nvar))
      gcm%APEX%inlist(1) = SIGMAP
      gcm%APEX%inlist(2) = SIGMAH

      gcm%isRestart = isRestart

      write(*,*) "done init_gcm_mpi"
      
    end subroutine init_gcm_mpi

    subroutine init_gcm_mix_mpi(gcm,gcmCplComm,gcmCplRank)
      type(gcm_T), intent(inout) :: gcm
      type(MPI_Comm) :: gcmCplComm
      integer, intent(in) :: gcmCplRank

      integer :: g,icolor

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call init_gcm_vars(gcm)
      
      if (gcm%isWACCMX) then 
        gcm%gcmId = wcmxId
      else 
        gcm%gcmId = tgcmId 
      endif

      !Find grid dimension
 
      write(*,*) "MIXCPL Waiting for Grid"

      ! Only obtain the grid information for grids we will couple with
      icolor = 0
      do g=1,2
        select case(g)
        case (1)
          if (gcm%GEO%mix2gcm_nvar > 0 .or. gcm%GEO%gcm2mix_nvar > 0) then
            call init_gcm_mpi_import(gcm%GEO,gcmCplComm,gcmCplRank,gcm%gcmId,icolor)
            write(*,*) "MIXCPL LOAD GEO"
          endif
        case (2)
          if (gcm%APEX%mix2gcm_nvar > 0 .or. gcm%APEX%gcm2mix_nvar > 0) then
            call init_gcm_mpi_import(gcm%APEX,gcmCplComm,gcmCplRank,gcm%gcmId,icolor)
            write(*,*) "MIXCPL LOAD APEX"
          endif
        end select
      end do
    end subroutine init_gcm_mix_mpi

    subroutine init_gcm_mpi_import(gcm,gcmCplComm,gcmCplRank,gcmId,icolor)
        type(gcm_grid_T),intent(inout):: gcm
        type(MPI_Comm), intent(in) :: gcmCplComm
        integer, intent(in) :: gcmCplRank,gcmId
        integer, intent(inout) :: icolor

        integer :: nlon,nlat,Nh,ierr

        Nh = GCMhemispheres

        icolor = icolor + 1
        call mpi_recv(gcm%nlat, 1, MPI_INTEGER, gcmCplRank, (gcmId+voltId)*100+icolor, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT NLAT: ",gcm%nlat

        icolor = icolor + 1
        call mpi_recv(gcm%nlon, 1, MPI_INTEGER, gcmCplRank, (gcmId+voltId)*100+icolor, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT NLON: ",gcm%nlon


        if (.not.allocated(gcm%inlat)) allocate(gcm%inlat(gcm%nlat))
        if (.not.allocated(gcm%inlon)) allocate(gcm%inlon(gcm%nlon))

        icolor = icolor + 1
        call mpi_recv(gcm%inlat, gcm%nlat, MPI_DOUBLE_PRECISION, gcmCplRank, (gcmId+voltId)*100+icolor, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT LAT: ",gcm%inlat
        icolor = icolor + 1
        call mpi_recv(gcm%inlon, gcm%nlon, MPI_DOUBLE_PRECISION, gcmCplRank, (gcmId+voltId)*100+icolor, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT LON: ",gcm%inlon
        write(*,*) "MIXCPL GOT GRID INFO: ",gcm%nlat,gcm% nlon

    end subroutine init_gcm_mpi_import

    subroutine coupleGCM2MIX(gcm,ion,mjd,time,gcmCplComm,gcmRank,myRank)
        type(mixIon_T),dimension(:),intent(inout) :: ion
        type(gcm_T), intent(inout) :: gcm
        real(rp), intent(in) :: time, mjd
        type(MPI_Comm), optional :: gcmCplComm
        integer, optional, intent(in) :: gcmRank,myRank

        ! maybe create the SM and GEO list of points here
        ! since the destination grid does not need to be structured
        ! can do a simple loop over all grid points and transform
        ! transform all remix points to GEO
        ! transform all gcm points to SM

        !Must MIX export first.  TIEGCM will also import first.
        call Tic("Export")
        if (present(gcmCplComm)) then
            call exportgcm(ion,gcm,mjd,time,gcmCplComm,gcmRank,myRank)
        else
            write(*,*) "Are we trying to Export to Couple GCM?"
        endif
        call Toc("Export")

        !Import gcm data
        call Tic("Import")
        if (present(gcmCplComm)) then
            call importgcm(gcm, gcmCplComm,gcmRank)
        else
            write(*,*) "Are we trying to Import to Couple GCM?"
        endif
        call process_gcmimport(gcm,ion)
        call Toc("Import")

        if (gcm%isRestart) gcm%isRestart=.false.
        ! We have a separate internal counter for coupling here.
        ! This may be used later on for WACCM-X coupling which is desync from remix coupling time
        ! TIEGCM coupling time is also 5s while WACCM-X will couple at 1 min default
        !gcm%cplStep = gcm%cplStep + 1
        
    end subroutine coupleGCM2MIX


    subroutine importgcmmpi(gcm, gcmCplComm,gcmCplRank)
        type(gcm_T), intent(inout) :: gcm
        type(MPI_Comm) :: gcmCplComm
        integer, intent(in) :: gcmCplRank

        integer :: v

        call Tic("MpiExchange")
        write(*,*) "MIX: IMPORT GCM"

        call Tic("Passing")
        do v=1,2
          select case(v)
          case (1)
            if (gcm%GEO%gcm2mix_nvar > 0) then
            !write(*,*) "Import GEO stuff"
              call import_gcm_per_grid(gcm%GEO,gcmCplComm,gcmCplRank,gcm%gcmId)
            end if
          case (2)
            if (gcm%APEX%gcm2mix_nvar > 0) then
            !write(*,*) "Import APEX stuff"
              call import_gcm_per_grid(gcm%APEX,gcmCplComm,gcmCplRank,gcm%gcmId)
            end if
          end select

        enddo
        call Toc("Passing")

        call Toc("MpiExchange")

    end subroutine importgcmmpi

    subroutine import_gcm_per_grid(gcm,gcmCplComm,gcmCplRank,gcmId)
        type(gcm_grid_T) :: gcm
        type(MPI_Comm) :: gcmCplComm
        integer, intent(in) :: gcmCplRank,gcmId

        integer :: ierr,length
        character( len = MPI_MAX_ERROR_STRING) :: message

        if (.not. allocated(gcm%invar2d)) allocate(gcm%invar2d(gcm%nlat,gcm%nlon,gcm%gcm2mix_nvar))

        call mpi_recv(gcm%invar2d, gcm%nlat*gcm%nlon*gcm%gcm2mix_nvar, MPI_DOUBLE_PRECISION, gcmCplRank, (gcmId+voltId)*100, gcmCplComm, MPI_STATUS_IGNORE,ierr)

        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
    end subroutine import_gcm_per_grid

    subroutine exportgcmmpi(ion,gcm,mjd,time,gcmCplComm,gcmCplRank,myCplRank)
        type(gcm_T), intent(inout) :: gcm
        type(mixIon_T), dimension(:),intent(inout) :: ion
        real(rp), intent(in) :: time, mjd
        type(MPI_Comm) :: gcmCplComm
        integer, intent(in) :: gcmCplRank,myCplRank

        integer :: g,ierr,length
        real(rp),allocatable,dimension(:,:) :: auroralbc
        character( len = MPI_MAX_ERROR_STRING) :: message

        write(*,*) "MIX: EXPORT GCM"

        ! Prepare the export data
        call Tic("Transform")
        do g = 1,2
            select case (g)
            case (1)
              if (gcm%GEO%mix2gcm_nvar > 0 .or. gcm%GEO%gcm2mix_nvar > 0) then
                call transform_gcm_export(gcm%GEO,ion,g)
              endif
            case(2)
              if (gcm%APEX%mix2gcm_nvar > 0 .or. gcm%APEX%gcm2mix_nvar > 0) then
                call transform_gcm_export(gcm%APEX,ion,g)
              endif
            end select

        end do
        call Toc("Transform")

        !! Calculate the auroral boundary in APEX coordinates
        !call calculate_auroral_boundary(gcm%APEX,auroralbc)

        call Tic("Exchange")

        do g = 1,2
            select case (g)
            case (1)
              if (gcm%GEO%mix2gcm_nvar > 0) then ! Skip if not passing
                call export_gcm_per_grid(gcm%GEO,gcmCplComm,gcmCplRank,myCplRank,gcm%gcmId,gcm%isWACCMX,gcm%nWACCMX)
              end if
            case(2)
              if (gcm%APEX%mix2gcm_nvar > 0) then ! Skip if not passing
                call export_gcm_per_grid(gcm%APEX,gcmCplComm,gcmCplRank,myCplRank,gcm%gcmId,gcm%isWACCMX,gcm%nWACCMX)
              end if
            end select
        enddo
        
        call Toc("Exchange")
    end subroutine exportgcmmpi

    subroutine export_gcm_per_grid(gcm,gcmCplComm,gcmCplRank,myCplRank,gcmId,isWACCMX,nWACCMX)
        type(gcm_grid_T),intent(inout) :: gcm
        type(MPI_Comm) :: gcmCplComm
        integer, intent(in) :: gcmCplRank, myCplRank, gcmId
        logical, intent(in) :: isWACCMX
        integer, intent(in) :: nWACCMX

        integer :: ierr,length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! Send the coupling data
        !write(*,*) " MIXCPL: ", gcmCplRank,(tgcmId+voltId)*100,gcmCplComm,gcm%nlat,gcm%nlon
        ! If WACCMX then each WACCMX rank will receive the data

        ! If this is slow, we may need to consider switching to isends

        if (isWACCMX) then
            write(*,*) " MIXCPL export1: ",minval(gcm%outvar2d(:,:,1)),maxval(gcm%outvar2d(:,:,1))
            write(*,*) " MIXCPL export2: ",minval(gcm%outvar2d(:,:,2)),maxval(gcm%outvar2d(:,:,2))
            write(*,*) " MIXCPL export3: ",minval(gcm%outvar2d(:,:,3)),maxval(gcm%outvar2d(:,:,3))
            call mpi_bcast(gcm%outvar2d, gcm%nlat*gcm%nlon*gcm%mix2gcm_nvar, MPI_DOUBLE_PRECISION, myCplRank, gcmCplComm, ierr)
        else
            call mpi_send(gcm%outvar2d, gcm%nlat*gcm%nlon*gcm%mix2gcm_nvar, MPI_DOUBLE_PRECISION, gcmCplRank, (gcmId+voltId)*100, gcmCplComm, ierr)
        endif

        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine export_gcm_per_grid

    subroutine getCplRanks(IAm, cplSize, appId, cplRank, nRanks, rankList)
        integer, dimension(:), intent(inout) :: IAm
        integer, intent(in) :: cplSize, appId
        integer, dimension(:), allocatable, intent(out) :: rankList
        integer, intent(out) :: nRanks, cplRank
        integer :: i, tempList(cplSize)  ! Temporary storage for ranks
        nRanks = 0
        do i=1, CplSize
            if (IAm(i) == appId*100 .or. IAm(i) == appId) then
                nRanks = nRanks + 1
                tempList(nRanks) = i-1
                IAm(i) = appId
            elseif (iAM(i) == appId*100 + 1) then
                ! is master rank
                nRanks = nRanks + 1
                tempList(nRanks) = i-1
                cplRank = i-1
                IAm(i) = appId
            endif
        enddo
        if (nRanks > 0) then
            allocate(rankList(nRanks))
            rankList = tempList(1:nRanks)
        endif
        if (nRanks == 1) cplRank = rankList(1)
    end subroutine getCplRanks

end module gcm_mpi
