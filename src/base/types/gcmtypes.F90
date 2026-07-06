module gcmtypes
  use shellGrid

  
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

  type gcm_T
      type(gcm_grid_T) :: GEO,APEX
      logical :: isRestart, isWACCMX = .false.
      integer :: gcmId, nWACCMX
      real(rp) :: altmax = 140. !used for apex, setting it randomly for now, km
      real(rp) :: ionalt
  end type gcm_T


  
end module gcmtypes
