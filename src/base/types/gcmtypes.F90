module gcmtypes
  use shellGrid
  use volttypes, only : ionoCoupler_T

  
  implicit none
  
  !integer,parameter :: geo2mix_nvar = 0 ! Neutral Winds(eventually)
  !integer,parameter :: apx2mix_nvar = 2 ! SigmaP, SigmaH, Neutral Winds(eventually)
  !integer :: gcm2mix_nvar = geo2mix_nvar + apx2mix_nvar ! SigmaP, SigmaH, Neutral Winds(eventually)
  !integer,parameter :: mix2apx_nvar = 1 ! pot
  !integer,parameter :: mix2geo_nvar = 2 ! eng, flux
  !integer :: mix2gcm_nvar = mix2apx_nvar + mix2geo_nvar ! pot, eng, flux

  integer,parameter :: GCMSIGMAP = 1,GCMSIGMAH = 2, GCMNORTH=NORTH, GCMSOUTH=SOUTH
  integer,parameter :: GCMhemispheres=2

  
  type gcm_grid_T
      type(ShellGrid_T) :: grid
      type(ShellGridVar_T), allocatable, dimension(:) :: gcmInput,gcmOutput
      real(rp),dimension(:),allocatable :: inlat, inlon
      real(rp),dimension(:),allocatable :: glon,gclat
      integer :: nlat,nlon,nlev,nhlat,lonshift
      character(len=32), dimension(:), allocatable :: outlist,inlist
      integer, dimension(:), allocatable :: t2N, t2S
      real(rp), dimension(:,:,:), allocatable :: invar2d, outvar2d
      integer :: mix2gcm_nvar
      !real(rp), dimension(:), allocatable :: lev,lon,clat,lat
      integer :: gcm2mix_nvar

  end type gcm_grid_T

  ! ================================================================
  ! gcm_T — Legacy GCM data container (used internally by tgcm.F90)
  ! ================================================================
  type :: gcm_T
      type(gcm_grid_T) :: GEO, APEX
      logical :: isRestart = .false.
      logical :: isWACCMX = .false.
      integer :: gcmId = -1
      integer :: nWACCMX = 0
      real(rp) :: altmax = 140.0_rp
      real(rp) :: ionalt = 0.0_rp
  end type gcm_T

  ! =========================================================
  ! gcmCoupler_T — GCM coupler (TIEGCM / WACCM-X)
  ! Extends ionoCoupler_T so it participates in the standard
  ! Voltron iono cycle and can live in vApp%ionoApp.
  ! Replaces direct use of gcm_T in the coupling layer.
  ! =========================================================
  type, extends(ionoCoupler_T) :: gcmCoupler_T

      ! Grid data for GEO (geographic) and APEX coordinate systems
      type(gcm_grid_T) :: GEO, APEX

      ! Control/state flags
      logical :: isRestart = .false.
      logical :: isWACCMX = .false.
      integer :: gcmId = -1
      integer :: nWACCMX = 0

      ! Physical parameters
      real(rp) :: altmax = 140.0_rp  !! [km] altitude max for APEX
      real(rp) :: ionalt = 0.0_rp    !! [km] ionosphere altitude

      contains

      procedure :: toIono       => gcmCpl_toIono
      procedure :: fromIono     => gcmCpl_fromIono
      procedure :: packInputs   => gcmCpl_packInputs
      procedure :: InitCouplers => gcmCpl_InitCouplers

      procedure :: InitModel           => gcmCplInitModel
      procedure :: InitIO              => gcmCplInitIO
      procedure :: WriteRestart        => gcmCplWriteRestart
      procedure :: ReadRestart         => gcmCplReadRestart
      procedure :: WriteConsoleOutput  => gcmCplWriteConsoleOutput
      procedure :: WriteFileOutput     => gcmCplWriteFileOutput
      procedure :: WriteSlimFileOutput => gcmCplWriteSlimFileOutput
      procedure :: AdvanceModel        => gcmCplAdvanceModel
      procedure :: Cleanup             => gcmCplCleanup

  end type gcmCoupler_T

end module gcmtypes

! ============================================================================
! gcminterp — Backward-compatible module interface to GCM routines
! These procedures are now implemented in volttypes submodule (gcmCplTypesSub)
! This wrapper module provides access for existing code that imports gcminterp.
! ============================================================================
module gcminterp
    use volttypes, only : init_gcm, init_gcm_grid, init_gcm_vars, create_gcm_grid, &
                         transform_gcm_export, process_gcm_import

    implicit none

    ! Re-export the submodule procedures
    public :: init_gcm, init_gcm_grid, init_gcm_vars, create_gcm_grid, &
              transform_gcm_export, process_gcm_import

end module gcminterp
