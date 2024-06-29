module gcmtypes
  use mixtypes
  
  implicit none
  
  !integer,parameter :: geo2mix_nvar = 0 ! Neutral Winds(eventually)
  !integer,parameter :: apx2mix_nvar = 2 ! SigmaP, SigmaH, Neutral Winds(eventually)
  !integer :: gcm2mix_nvar = geo2mix_nvar + apx2mix_nvar ! SigmaP, SigmaH, Neutral Winds(eventually)
  !integer,parameter :: mix2apx_nvar = 1 ! pot
  !integer,parameter :: mix2geo_nvar = 2 ! eng, flux
  !integer :: mix2gcm_nvar = mix2apx_nvar + mix2geo_nvar ! pot, eng, flux

  integer,parameter :: GCMSIGMAP = 1,GCMSIGMAH = 2, GCMNORTH=NORTH, GCMSOUTH=SOUTH
  integer,parameter :: GCMhemispheres=2

  
  type var_T
      real(rp), dimension(:,:),allocatable :: var
  end type var_T

  type gcm_grid_T
      type(Map_T), allocatable, dimension(:) :: r2tMaps, t2rMaps
      type(mixGrid_T), allocatable, dimension(:) :: SM,MyGrid
      real(rp),dimension(:),allocatable :: inlat, inlon
      real(rp),dimension(:,:,:),allocatable :: glon,gclat
      type(var_T), allocatable, dimension(:,:) :: gcmInput,gcmOutput
      type(var_T), allocatable, dimension(:,:) :: mixInput,mixOutput
      integer :: nlat,nlon,nlev,nhlat,lonshift
      integer, dimension(:), allocatable :: outlist,inlist
      integer, dimension(:), allocatable :: t2N, t2S
      real(rp), dimension(:,:,:), allocatable :: invar2d, outvar2d
      integer :: mix2gcm_nvar
      !real(rp), dimension(:), allocatable :: lev,lon,clat,lat
      integer :: gcm2mix_nvar
      character (len=strLen) :: Coord
  end type gcm_grid_T

  type gcm_T
      !type(mixGrid_T) :: rGEOGrid,tGEOGrid,tSMGrid ! GEO&SM Grid for tgcm, remix GEO grid
      !type(mixGrid_T) :: mixGrid,gcmGrid
      type(gcm_grid_T) :: GEO,APEX
      !type(var_T), dimension(GCMhemispheres,geo2mix_nvar) :: geoInput
      !type(var_T), dimension(GCMhemispheres,apx2mix_nvar) :: apxInput
      !type(var_T), dimension(GCMhemispheres,gcm2mix_nvar) :: mixInput
      !type(var_T), dimension(GCMhemispheres,mix2gcm_nvar) :: geoOutput,apxOutput,mixOutput
      !real(rp), dimension(:,:,:,:), allocatable :: mixInput, mixOutput
      !real(rp), dimension(:,:,:,:), allocatable :: gcmInput, gcmOutput
      !integer :: nlat,nlon,nlev,nhlat,ntime,order,nhemi,lonshift
      !real(rp),dimension(:,:),allocatable :: gx,gy
      integer :: nhemi
      !integer :: gcm2mix_nvar, mix2gcm_nvar
      integer :: cplStep = 1
      !character(len=strlen) :: mix2gcmH5,gcm2mixH5,mix2gcmLock,gcm2mixLock
      !character(len=strlen) :: mix2gcmH5 = "mix4gcm.h5"
      !character(len=strlen) :: mix2gcmLock = "mix2gcm.txt"
      !character(len=strlen) :: gcm2mixH5 = "gcm4mix.h5"
      !character(len=strlen) :: gcm2mixLock = "gcm2mix.txt"
      logical :: isRestart
      real(rp) :: altmax = 140. !used for apex, setting it randomly for now, km
      real(rp) :: ionalt
  end type gcm_T


  
end module gcmtypes
