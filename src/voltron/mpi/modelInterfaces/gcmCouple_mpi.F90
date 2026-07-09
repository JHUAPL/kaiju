module gcmCouple_mpi
  use volttypes_mpi
  use volttypes
  use mpi_f08
  use mpidefs
  use clocks

  implicit none

  ! ================================================================
  ! gcmCoupler_mpi_T — MPI extension of gcmCoupler_T
  ! Holds MPI-specific data and implements MPI-based initialization
  ! ================================================================
  type, extends(gcmCoupler_T) :: gcmCoupler_mpi_T

      ! MPI coupling communicator and rank information
      type(MPI_Comm) :: couplingComm = MPI_COMM_NULL
      integer :: gcmCplRank = -1
      integer :: myCplRank = -1   !! My rank in the coupling communicator

      contains

      ! MPI-aware initialization of GCM coupling
      !procedure :: toIono       => gcmCpl_toIono
      !procedure :: fromIono     => gcmCpl_fromIono
      !procedure :: packInputs   => gcmCpl_packInputs
      procedure :: InitCouplers  => gcmCpl_InitCouplers_mpi
      procedure :: Coupling      => gcmCpl_Coupling_mpi

      !procedure :: InitModel           => gcmCplInitModel
      !procedure :: InitIO              => gcmCplInitIO
      !procedure :: WriteRestart        => gcmCplWriteRestart
      !procedure :: ReadRestart         => gcmCplReadRestart
      !procedure :: WriteConsoleOutput  => gcmCplWriteConsoleOutput
      !procedure :: WriteFileOutput     => gcmCplWriteFileOutput
      !procedure :: WriteSlimFileOutput => gcmCplWriteSlimFileOutput
      !procedure :: AdvanceModel        => gcmCplAdvanceModel
      !procedure :: Cleanup             => gcmCplCleanup

  end type gcmCoupler_mpi_T

contains

    ! ================================================================
    ! gcmCpl_initCoupler_mpi — Initialize GCM coupler with MPI support
    ! 
    ! Converts this%gcmCoupler_mpi_T to a temporary gcm_T for compatibility
    ! with existing gcm_mpi initialization routines, then copies results back.
    ! Calls init_gcm_mix_mpi to:
    ! - Initialize grid information from GCM ranks
    ! - Set up MPI communicator for data exchange
    ! ================================================================
    subroutine gcmCpl_InitCouplers_mpi(App, vApp, xml)
        class(gcmCoupler_mpi_T), intent(inout) :: App
        class(voltApp_T),     intent(inout) :: vApp
        type(XML_Input_T),       intent(inout) :: xml

        call init_gcm_mix_mpi(App)

        call init_gcm(App,vApp%gApp%Model%isRestart)

        ! Type guard to access MPI-specific members
        select type(vApp)
        type is (voltAppMpi_T)
            App%couplingComm = vApp%mageCplComm
            App%myCplRank = vApp%voltCplRank
            App%gcmCplRank = vApp%gcmCplRank
        class default
            ! Non-MPI call or error handling
            write(*,*) "ERROR: gcmCpl_InitCouplers_mpi requires voltAppMpi_T"
        end select

    end subroutine gcmCpl_InitCouplers_mpi

    subroutine gcmCpl_Coupling_mpi(App, vApp)
        class(gcmCoupler_mpi_T), intent(inout) :: App
        class(voltApp_T),     intent(inout) :: vApp

        call coupleGCM2MIX(App)

    end subroutine gcmCpl_Coupling_mpi

    
    subroutine init_gcm_mix_mpi(gcm)
      class(gcmCoupler_mpi_T), intent(inout) :: gcm

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
            call init_gcm_mpi_import(gcm%GEO,gcm%couplingComm,gcm%gcmCplRank,gcm%gcmId,icolor)
            write(*,*) "MIXCPL LOAD GEO"
          endif
        case (2)
          if (gcm%APEX%mix2gcm_nvar > 0 .or. gcm%APEX%gcm2mix_nvar > 0) then
            call init_gcm_mpi_import(gcm%APEX,gcm%couplingComm,gcm%gcmCplRank,gcm%gcmId,icolor)
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

        integer :: nlon,nlat,ierr

        icolor = icolor + 1
        write(*,*) "MIXCPL Import 1:", gcmCplRank, (gcmId+voltId)*100+icolor
        call mpi_recv(gcm%nlat, 1, MPI_INTEGER, gcmCplRank, (gcmId+voltId)*100+icolor, gcmCplComm, MPI_STATUS_IGNORE,ierr)
        write(*,*) "MIXCPL GOT NLAT: ",gcm%nlat

        icolor = icolor + 1
        write(*,*) "MIXCPL Import 1:", gcmCplRank, (gcmId+voltId)*100+icolor
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
        write(*,*) "MIXCPL GOT GRID INFO: ",gcm%nlat,gcm%nlon

    end subroutine init_gcm_mpi_import

    subroutine coupleGCM2MIX(gcm)
        class(gcmCoupler_mpi_T), intent(inout) :: gcm

        !Must MIX export first.  TIEGCM will also import first.
        call Tic("Export")
        if (gcm%couplingComm /= MPI_COMM_NULL) then
            call exportgcmmpi(gcm)
        else
            write(*,*) "FAILURE: Are we trying to Export to Couple GCM?"
        endif
        call Toc("Export")

        !Import gcm data
        call Tic("Import")
        if (gcm%couplingComm /= MPI_COMM_NULL) then
            call importgcmmpi(gcm)
        else
            write(*,*) "FAILURE: Are we trying to Import to Couple GCM?"
        endif
        call Toc("Import")

        if (gcm%isRestart) gcm%isRestart=.false.
        ! We have a separate internal counter for coupling here.
        ! This may be used later on for WACCM-X coupling which is desync from remix coupling time
        ! TIEGCM coupling time is also 5s while WACCM-X will couple at 1 min default
        !gcm%cplStep = gcm%cplStep + 1
        
    end subroutine coupleGCM2MIX


    subroutine importgcmmpi(gcm)
        class(gcmCoupler_mpi_T), intent(inout) :: gcm

        integer :: v

        call Tic("MpiExchange")
        write(*,*) "MIX: IMPORT GCM"

        call Tic("Passing")
        do v=1,2
          select case(v)
          case (1)
            if (gcm%GEO%gcm2mix_nvar > 0) then
            !write(*,*) "Import GEO stuff"
              call import_gcm_per_grid(gcm%GEO,gcm%couplingComm,gcm%gcmCplRank,gcm%gcmId)
            end if
          case (2)
            if (gcm%APEX%gcm2mix_nvar > 0) then
            !write(*,*) "Import APEX stuff"
              call import_gcm_per_grid(gcm%APEX,gcm%couplingComm,gcm%gcmCplRank,gcm%gcmId)
            end if
          end select

        enddo
        call Toc("Passing")

        if (gcm%GEO%gcm2mix_nvar > 0) then
            call process_gcm_import(gcm%GEO)
        end if
        if (gcm%APEX%gcm2mix_nvar > 0) then
            call process_gcm_import(gcm%APEX)
        end if

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

    subroutine exportgcmmpi(gcm)
        class(gcmCoupler_mpi_T), intent(inout) :: gcm

        integer :: g,ierr,length
        real(rp),allocatable,dimension(:,:) :: auroralbc
        character( len = MPI_MAX_ERROR_STRING) :: message

        write(*,*) "MIX: EXPORT GCM"

        ! Prepare the export data
        call Tic("Transform",.true.)
        do g = 1,2
            select case (g)
            case (1)
              if (gcm%GEO%mix2gcm_nvar > 0 .or. gcm%GEO%gcm2mix_nvar > 0) then
                call load_gcm_export(gcm%GEO)
              endif
            case(2)
              if (gcm%APEX%mix2gcm_nvar > 0 .or. gcm%APEX%gcm2mix_nvar > 0) then
                call load_gcm_export(gcm%APEX)
              endif
            end select

        end do
        call Toc("Transform",.true.)

        !! Calculate the auroral boundary in APEX coordinates
        !call calculate_auroral_boundary(gcm%APEX,auroralbc)

        call Tic("Exchange",.true.)

        do g = 1,2
            select case (g)
            case (1)
              if (gcm%GEO%mix2gcm_nvar > 0) then ! Skip if not passing
                call export_gcm_per_grid(gcm%GEO,gcm%couplingComm,gcm%gcmCplRank,gcm%myCplRank,gcm%gcmId,gcm%isWACCMX,gcm%nWACCMX)
              end if
            case(2)
              if (gcm%APEX%mix2gcm_nvar > 0) then ! Skip if not passing
                call export_gcm_per_grid(gcm%APEX,gcm%couplingComm,gcm%gcmCplRank,gcm%myCplRank,gcm%gcmId,gcm%isWACCMX,gcm%nWACCMX)
              end if
            end select
        enddo
        
        call Toc("Exchange",.true.)
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

    
    ! ================================================================
    ! GCM initialization routines (from former tgcm.F90 / gcminterp)
    ! ================================================================

    subroutine init_gcm(gcm, isRestart)
        class(gcmCoupler_T), intent(inout) :: gcm
        logical, intent(in) :: isRestart

        write(*,*) "start init_gcm"

        call init_gcm_grid(gcm)

        write(*,*) "done init_gcm_grid"

        if (.not. allocated(gcm%GEO%outlist)) allocate(gcm%GEO%outlist(gcm%GEO%mix2gcm_nvar))
        if (.not. allocated(gcm%APEX%outlist)) allocate(gcm%APEX%outlist(gcm%APEX%mix2gcm_nvar))

        if (gcm%isWACCMX) then
            gcm%APEX%outlist(1) = "POT"
            gcm%APEX%outlist(2) = "AVG_ENG"
            gcm%APEX%outlist(3) = "NUM_FLUX"
        else ! is TIEGCM
            gcm%APEX%outlist(1) = "POT"
            gcm%APEX%outlist(2) = "NUM_FLUX"
            gcm%GEO%outlist(1) = "AVG_ENG"
            gcm%GEO%outlist(2) = "NUM_FLUX"
        endif

        if ( .not. allocated(gcm%GEO%inlist)) allocate(gcm%GEO%inlist(gcm%GEO%gcm2mix_nvar))
        if ( .not. allocated(gcm%APEX%inlist)) allocate(gcm%APEX%inlist(gcm%APEX%gcm2mix_nvar))
        gcm%APEX%inlist(1) = "SIGMAP"
        gcm%APEX%inlist(2) = "SIGMAH"

        gcm%isRestart = isRestart

        write(*,*) "done init_gcm"

    end subroutine init_gcm

    subroutine init_gcm_grid(gcm)
        class(gcmCoupler_T), intent(inout) :: gcm

        if (gcm%GEO%mix2gcm_nvar > 0) call create_gcm_grid(gcm%GEO, "GEO")
        if (gcm%APEX%mix2gcm_nvar > 0) call create_gcm_grid(gcm%APEX, "APEX")

    end subroutine init_gcm_grid

    subroutine init_gcm_vars(gcm)
        class(gcmCoupler_T), intent(inout) :: gcm

        if (gcm%isWACCMX) then
            write(*,*) "MIXCPL Coupling to WACCMX: ", gcm%isWACCMX
            gcm%GEO%gcm2mix_nvar = 0 ! Neutral Winds (eventually)
            gcm%GEO%mix2gcm_nvar = 0 ! Nothing passed on GEO

            gcm%APEX%gcm2mix_nvar = 2 ! SigmaP, sigmaH
            gcm%APEX%mix2gcm_nvar = 3 ! pot, eng, flux
        else ! is TIEGCM
            write(*,*) "MIXCPL Coupling to TIEGCM: ", (.not. gcm%isWACCMX)
            gcm%GEO%gcm2mix_nvar = 0 ! Neutral Winds (eventually)
            gcm%GEO%mix2gcm_nvar = 2 ! eng, flux

            gcm%APEX%gcm2mix_nvar = 2 ! SigmaP, sigmaH
            gcm%APEX%mix2gcm_nvar = 2 ! pot, flux (for auroral boundary)
        endif

    end subroutine init_gcm_vars

    subroutine create_gcm_grid(gcm, coord)
        type(gcm_grid_T), intent(inout) :: gcm
        character(len=*), intent(in) :: coord

        integer :: nGhosts(4)
        integer :: Np, Nt, Nh, i, j, h, v

        ! nlon includes periodic ghost

        if (.not.allocated(gcm%gclat)) allocate(gcm%gclat(gcm%nlat))
        if (.not.allocated(gcm%glon)) allocate(gcm%glon(gcm%nlon))

        ! Ascending Lat starts in southern hemisphere
        ! Ascending Clat starts in northern hemisphere
        gcm%gclat = 90. - gcm%inlat
        gcm%gclat = gcm%gclat(gcm%nlat:1:-1) ! reverse order. Ascending.

        ! Check first and last point of clat to make sure it's 0 and 180
        if (.not. (gcm%gclat(1) < TINY .and. gcm%gclat(gcm%nlat) > 180.-TINY)) &
            write(*,*) 'CLAT RANGE ERROR: ', gcm%gclat(1), gcm%gclat(gcm%nlat)

        ! find the point where longitude is 0, so we know how much we need to shift by
        gcm%lonshift = minloc(abs(gcm%inlon(:gcm%nlon-1)), 1) - 1

        gcm%glon(:gcm%nlon-1) = modulo(cshift(gcm%inlon(1:gcm%nlon-1), gcm%lonshift)+360., 360.)

        if (abs(gcm%glon(1)) > 0.01) write(*,*) "WARNING: GCM lonshift is not near 0.0 degrees: ", gcm%glon(1)

        gcm%glon(1) = 0.0_rp ! Clamp to 0.0 just in case
        gcm%glon(gcm%nlon) = 360.0 ! Clamp final value exactly to 360.0

        write(*,*) "GCM Grid GLON: ", gcm%glon
        write(*,*) "GCM Grid GCLAT: ", gcm%gclat

        gcm%glon = gcm%glon * deg2rad
        gcm%gclat = gcm%gclat * deg2rad

        nGhosts = (/0, 0, 0, 0/)
        call GenShellGrid(gcm%grid, gcm%gclat, gcm%glon, "GCM", nGhosts=nGhosts, coord=coord)

        if (.not. allocated(gcm%gcmInput)) allocate(gcm%gcmInput(gcm%gcm2mix_nvar))
        do i = 1, gcm%gcm2mix_nvar
            call initShellVar(gcm%grid, SHGR_CORNER, gcm%gcmInput(i))
        end do

        if (.not. allocated(gcm%gcmOutput)) allocate(gcm%gcmOutput(gcm%mix2gcm_nvar))
        do i = 1, gcm%mix2gcm_nvar
            call initShellVar(gcm%grid, SHGR_CORNER, gcm%gcmOutput(i))
        end do

        write(*,*) "Created GCM grid with nlat,nlon,nhlat: ", gcm%nlat, gcm%nlon
        write(*,*) "GCM Grid Check: ", shape(gcm%gcmOutput(1)%data), gcm%grid%Nt, gcm%grid%Np
        write(*,*) "GCM Lower bounds: ", lbound(gcm%gcmOutput(1)%data)
        write(*,*) "GCM Upper bounds: ", ubound(gcm%gcmOutput(1)%data)

    end subroutine create_gcm_grid

    subroutine load_gcm_export(gcm)
        type(gcm_grid_T), intent(inout) :: gcm

        integer :: v

        if (gcm%mix2gcm_nvar .eq. 0) return

        if (.not. allocated(gcm%outvar2d)) allocate(gcm%outvar2d(gcm%nlat, gcm%nlon, gcm%mix2gcm_nvar))

        gcm%outvar2d = 0._rp
        do v = 1, gcm%mix2gcm_nvar
            ! Invert latitude to map colat -> lat
            ! undo the lonshift by applying it backwards
            gcm%outvar2d(gcm%nlat:1:-1, :gcm%nlon-1, v) = cshift(gcm%gcmOutput(v)%data(:, :gcm%nlon-1), -1*gcm%lonshift, dim=2)
            gcm%outvar2d(:, gcm%nlon, v) = gcm%outvar2d(:, 1, v)
        end do

    end subroutine load_gcm_export

    subroutine process_gcm_import(gcm)
        type(gcm_grid_T), intent(inout) :: gcm

        integer :: v

        ! gcmInput lives in colat space on the Shell Grid
        ! invar2d is the direct lat space from GCM exchange
        do v = 1, gcm%gcm2mix_nvar
            gcm%gcmInput(v)%data = 0._rp

            gcm%gcmInput(v)%data(gcm%nlat:1:-1, :gcm%nlon-1) = cshift(gcm%invar2d(:, :gcm%nlon-1, v), gcm%lonshift, dim=2)
            gcm%gcmInput(v)%data(:, gcm%nlon) = gcm%gcmInput(v)%data(:, 1)
        end do

    end subroutine process_gcm_import

end module gcmCouple_mpi