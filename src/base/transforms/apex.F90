module apex
  use kdefs
  use dates
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! tiegcmlicense.txt.
!
! 2018 May: A.D. Richmond modifications to add full base vectors f3,g1,g2,g3;
!   to update comments, to unify specification of constants, and to add
!   subroutine apex_setup that simplifies initialization for many purposes.
! 2013 April: B. Foster (NCAR/HAO)
!   This is a refactored version of the legacy apex code, originally written
!   by Art Richmond and Roy Barnes in the 1993-2000 timeframe, using
!   field-line tracing routines developed in 1973 by Wally Clark at NOAA.
!   This version is written in free-format fortran90. Subroutines and
!   module data may be use-associated from this module.
!
! Reference: Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
!   Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995,
!   https://doi.org/10.5636/jgg.47.191
!
! A typical calling sequence for a code calling this module is as follows:
!
! subroutine apex_setup(date,altmax):
!   Set up 3D arrays over the Earth for geomagnetic-field epoch date,
!    up to an altitude of at least altmax,
!    for interpolation use when calling apex_mall, apex_q2g, or apex_m2g.
!    (This must be called before apex_mall, apex_q2g, or apex_m2g)
!   The resolution parameter nvert is set to 40 here, which results in an array
!    of lat,lon dimensions 121,201, adequate for most purposes.
!    The nalt dimension depends on nvert and altmax; for nvert=40 and altmax
!    between 910.2 and 1124.3, nalt=7 .  For very large altmax the maximum
!    nalt is nvert+1.
!
! subroutine apex_mall:
!   Calculate modified Apex coordinates, basis vectors described in
!    Richmond [1995], and other magnetic field parameters
!    (often called from lat,lon,alt nested loop)
!
! subroutine apex_q2g or apex_m2g:
!   Convert from quasi-dipole or modified magnetic-apex coordinates
!     to geodetic coordinates.
!
! Subroutines called by those above can be invoked directly to do coordinate
!   conversions; to get geomagnetic-field quantities and grid-point locations;
!   and to calculate the subsolar location, cosine of the solar zenith angle,
!   and magnetic local time.
!
  implicit none
  real(rp),parameter :: re = 6371.2,            & ! Reference radius of IGRF
                    eps = 1.e-5,            & ! Small number
                    fltnvrs = 298.257223563   ! Inverse flatness of geoid
  real(rp),allocatable,save :: &
    xarray(:,:,:), & ! cos(quasi-dipole latitude)*cos(apex longitude)
                     !   minus dipole value
    yarray(:,:,:), & ! cos(quasi-dipole latitude)*sin(apex longitude)
                     !   minus dipole value
    zarray(:,:,:), & ! sin(quasi-dipole latitude)
                     !   minus dipole value
    varray(:,:,:)    ! magnetic potential, normalized by radially varying
                     !   dipole value at S pole,
                     !   minus dipole value
!
  integer :: nglat,nglon,ngalt ! numbers of lat,lon,alt grid points,
!  set in apex_mka to nlat,nlon,nalt as determined by ggrid.
!
! This grid (geolat,geolon,geoalt is equivalent to gdlat,gdlon,gdalt,
! as passed to apex_mka.
  real(rp),allocatable,save :: geolat(:), geolon(:), geoalt(:)
  integer,parameter :: nmax=13 ! maximum degree of IGRF coefficients
  integer,parameter :: ncoef = nmax*nmax + 2*nmax + 1 ! 196
  real(rp),dimension(ncoef) :: &
    gb, & ! Coefficients for magnetic field calculation
    gv    ! Coefficients for magnetic potential calculation
!
 real(rp), parameter :: &
    dtr=0.0174532925199432957692369076847, & ! degrees to radians
    rtd=57.2957795130823208767981548147      ! radians to degrees
  real(rp) :: &
    pola ! Limiting latitude magnitude (deg); when the geographic
!          latitude is poleward of pola, then xarray,yarray,zarray,varray
!          are forced to be constant for all longitudes at each altitude.
  real(rp),parameter ::        & ! Formerly common /APXCON/
    req  = 6378.137,       & ! Equatorial earth radius
    precise = 7.6e-11,     & ! Precision factor
    glatlim = 89.9,        & ! Limit above which gradients are recalculated
    xmiss = -32767.,       & ! Set missing values to this
    req2 = req*req,        &
    req2e2 = req2*(2. - 1./fltnvrs)/fltnvrs
!
  real(rp) ::  & ! Formerly /APXDIPL/ and /DIPOLE/
    colat, & ! Geocentric colatitude of geomagnetic dipole north pole (deg)
    elon,  & ! East longitude of geomagnetic dipole north pole (deg)
    vp,    & ! Magnitude, in T.m, of dipole component of magnetic
             ! potential at geomagnetic pole and geocentric radius re
    ctp,stp  ! cosine and sine of colat
!
  real(rp) ::  & ! Formerly /FLDCOMD/
    bx,    & ! X comp. of field vector at the current tracing point (Gauss)
    by,    & ! Y comp. of field vector at the current tracing point (Gauss)
    bz,    & ! Z comp. of field vector at the current tracing point (Gauss)
    bb       ! Magnitude of field vector at the current tracing point (Gauss)
  real(rp) ::      & ! Formerly /APXIN/
    yapx(3,3)    ! Matrix of cartesian coordinates (loaded columnwise)
!
! /ITRA/ was only in subs linapx and itrace, so it could be removed from module
!  data if those routines are properly restructured to pass the /ITRA/
!  variables.
!
  integer ::   & ! Formerly /ITRA/
    nstp         ! Step count. Incremented in sub linapx.
  real(rp)    ::   &
    y(3),      & ! Array containing current tracing point cartesian coordinates.
    yp(3),     & ! Array containing previous tracing point cartesian coordinates.
    sgn,       & ! Determines direction of trace. Set in subprogram linapx
    ds           ! Step size (Km) Computed in subprogram linapx.
  real(rp) ::         & ! limits beyond which east-west gradients are computed
    glatmn,glatmx   ! differently to avoid potential underflow (apex_mka)
contains

!-----------------------------------------------------------------------

    subroutine init_apex(mjd,altmax)
      real(rp), intent(in) :: mjd,altmax
      integer :: iyr,idoy,imon,iday,ihr,imin
      real(rp) :: rsec
      real(rp) :: mjdin
      real(rp) :: yrfrac

      mjdin = mjd ! mjdin of mjd2ut is for inout intent.
      call mjd2ut(mjdin,iyr,idoy,imon,iday,ihr,imin,rsec)

      ! Apex setup requires fraction of year
      if (mod(iyr,4) .eq. 0) then
        yrfrac = iyr + idoy/366. + ihr/24./366. + imin/60./24./366. + rsec/60./60./24./366.
      else
        yrfrac = iyr + idoy/365. + ihr/24./365. + imin/60./24./365. + rsec/60./60./24./365.
      end if
        !write (*,'('' set up apex for date='',f7.2,'' altmax='',f10.1)') date,altmax
      call apex_setup(yrfrac,altmax)

    end subroutine init_apex

!-----------------------------------------------------------------------
subroutine apex_setup(date,altmax)
! Set up 3D arrays over the Earth for geomagnetic-field epoch date,
!  up to an altitude of at least altmax,
!  for use when calling apex_mall, apex_q2g, or apex_m2g.
  real(rp), intent(in) ::  &
               date,   & ! Date in year and fraction (e.g., 2018.5)
               altmax    ! Maximum altitude that will be used (km).
!
  integer,parameter :: nvert = 120 ! Resolution parameter, corresponding
!   to the maximum number of vertical grid increments when altmax=infinity.
!   Points are spaced uniformly in 1/r between the Earth's surface and an
!   altitude that is at least altmax (or a very large value if altmax=infinity).
  integer,parameter ::            &
               mxlat = 3*nvert+1, & ! number of latitude grid points
               mxlon = 5*nvert+1, & ! number of longitude grid points
               mxalt =   nvert+1    ! maximum number of height grid points
! Local
  integer :: nlat,nlon,nalt ! Numbers of lat, lon, alt grid points
  real(rp) :: gplat(mxlat),gplon(mxlon),gpalt(mxalt) ! grid latitudes (deg),
!   longitudes (deg), and altitudes (km), respectively
  integer :: ier ! Error flag, non-zero if there is a problem.
  nlat = 0
  nlon = 0
  nalt = 0
  call ggrid(nvert,-90._rp,90._rp,-180._rp,180._rp,0._rp,altmax, &
                 gplat,gplon,gpalt,mxlat,mxlon,mxalt,nlat,nlon,nalt)
  call apex_mka(date,gplat,gplon,gpalt,nlat,nlon,nalt,ier)
    if (ier /= 0) then
      write(6,"('>>> apex_setup error from apex_mka')")
      stop 'apex_setup apex_mka'
    endif
!
!
end subroutine apex_setup
!-----------------------------------------------------------------------
subroutine ggrid(nvert,glatmin,glatmax,glonmin,glonmax,altmin,altmax, &
                 gplat,gplon,gpalt,mxlat,mxlon,mxalt,nlat,nlon,nalt)
!
! Given desired range of geographic latitude, longitude and altitude,
! choose an appropriate grid that can be used in subsequent calls to
! subs apex_mka, apex_mall, apex_q2g.
!
! Input args:
  integer,intent(in) :: nvert,mxlat,mxlon,mxalt
  real(rp),intent(in) :: glatmin,glatmax,glonmin,glonmax,altmin,altmax
!
! Output args:
  integer,intent(out) :: nlat,nlon,nalt
  real(rp),intent(out) :: gplat(mxlat),gplon(mxlon),gpalt(mxalt)
!
! Local:
  real(rp) :: dlon,dlat,diht,dnv,glonmaxx,x
  integer :: nlatmin,nlatmax,nlonmin,nlonmax,naltmin,naltmax
  integer :: i,j,k,kk
!
! Check inputs:
  if (glatmin > glatmax) then
    write(6,"('>>> ggrid: glatmin=',f9.2,' must be <= glatmax=',f9.2)") &
      glatmin,glatmax
    stop 'ggrid'
  endif
  if (glonmin > glonmax) then
    write(6,"('>>> ggrid: glonmin=',f9.2,' must be <= glonmax=',f9.2)") &
      glonmin,glonmax
    stop 'ggrid'
  endif
  if (altmin > altmax) then
    write(6,"('>>> ggrid: altmin=',f9.2,' must be <= altmax=',f9.2)") &
      altmin,altmax
    stop 'ggrid'
  endif
!
! Init outputs:
  nlat = 0 ; nlon = 0 ; nalt = 0
  gplat = 0. ; gplon = 0. ; gpalt = 0.
!
  dnv = float(nvert)
  dlon = 360. / (5.*dnv)
  dlat = 180. / (3.*dnv)
  diht = 1.   / dnv
  nlatmin = max(int((glatmin+90.)/dlat),0)
  nlatmax = min(int((glatmax+90.)/dlat+1.),3*nvert)
  nlonmin = max(int((glonmin+180.)/dlon),0)

  glonmaxx = min(glonmax,glonmin+360.)
  nlonmax = min(int((glonmaxx+180.)/dlon+1.),10*nvert)

  x = re/(re+altmax)/diht-eps
  naltmin = max(x,1.)
  naltmin = min(naltmin,nvert-1)
  x = re/(re+altmin)/diht+eps
  i = x + 1.
  naltmax = min(i,nvert)
  nlat = nlatmax - nlatmin + 1
  nlon = nlonmax - nlonmin + 1
  nlon = min(nlon,5*nvert+1)
  nalt = naltmax - naltmin + 1
  do j=1,nlat
    gplat(j) = dlat*float(nlatmin+j-1) - 90.
  enddo
  gplat(nlat) = min(gplat(nlat),90.)
  do i=1,nlon
    gplon(i) = dlon*float(nlonmin+i-1) - 180.
  enddo
  do k=1,nalt
    kk = naltmax - k +1
    gpalt(k) = re*(float(nvert-kk) - eps) / (float(kk)+eps)
  enddo
  if (gplon(nlon-1) >= glonmax) nlon = nlon-1
  gpalt(1) = max(gpalt(1),0.)
! write(6,"('ggrid: nlat=',i4,' gplat=',/,(6f9.2))") nlat,gplat
! write(6,"('ggrid: nlon=',i4,' gplon=',/,(6f9.2))") nlon,gplon
! write(6,"('ggrid: nalt=',i4,' gpalt=',/,(6f9.2))") nalt,gpalt
end subroutine ggrid
!-----------------------------------------------------------------------
subroutine apex_mka(date,gplat,gplon,gpalt,nlat,nlon,nalt,ier)
!
! Given a 3d lat,lon,altitude grid, calculate xarray,yarray,zarray,varray in
! module data above. These arrays are used later for calculating quantities
! involving gradients of Apex coordinates, such as base vectors in the
! Modified-Apex and Quasi-Dipole systems.
!
! This defines module 3d data xarray,yarray,zarray,varray
!
! Input args:
  real(rp),intent(in) :: date              ! year and fraction
  integer,intent(in) :: nlat,nlon,nalt ! dimensions of 3d grid
  real(rp),intent(inout) :: gplat(nlat),gplon(nlon),gpalt(nalt)
!
! Output args:
  integer,intent(out) :: ier
!
! Local:
  integer :: i,j,k,kpol,istat
  real(rp) :: reqore,rqorm1,cp,ct,st,sp,stmcpm,stmspm,ctm
  real(rp) :: aht,alat,phia,bmag,xmag,ymag,zdown,vmp ! apex_sub output
  real(rp) :: vnor,rp,reqam1,a,slp,clp,phiar
!
! pola:
!   Limiting latitude magnitude (deg); when the geographic latitude is
!    poleward of pola, then xarray,yarray,zarray,varray are forced to
!    be constant for all longitudes at each altitude.
!
!   This makes pola = 89.9995 for precise = 7.6e-11:
  pola = 90.-sqrt(precise)*rtd    ! Pole angle (deg)
  ier = 0
!
! Allocate 3d xarray,yarray,zarray,varray arrays:
! Deallocate if necessary in case this routine is called more than once
! (e.g., called by each of coupled models with different grids)
!
  if (allocated(xarray)) deallocate(xarray)
  allocate(xarray(nlat,nlon,nalt),stat=istat)
  if (istat /= 0) stop 'allocate xarray'
  xarray = 0.
  if (allocated(yarray)) deallocate(yarray)
  allocate(yarray(nlat,nlon,nalt),stat=istat)
  if (istat /= 0) stop 'allocate yarray'
  yarray = 0.
  if (allocated(zarray)) deallocate(zarray)
  allocate(zarray(nlat,nlon,nalt),stat=istat)
  if (istat /= 0) stop 'allocate zarray'
  zarray = 0.
  if (allocated(varray)) deallocate(varray)
  allocate(varray(nlat,nlon,nalt),stat=istat)
  if (istat /= 0) stop 'allocate varray'
  varray = 0.
!
! Set geographic grid in module data for later reference:
! (these also are not deallocated by this module)
!
  nglon=nlon ; nglat=nlat ; ngalt=nalt
  if (allocated(geolat)) deallocate(geolat)
  if (allocated(geolon)) deallocate(geolon)
  if (allocated(geoalt)) deallocate(geoalt)
  allocate(geolat(nglat),stat=istat)
  allocate(geolon(nglon),stat=istat)
  allocate(geoalt(ngalt),stat=istat)
  geolat(:) = gplat(:)
  geolon(:) = gplon(:)
  geoalt(:) = gpalt(:)
!
! Set coefficients gb,gv (module data) for requested year:
!
  call cofrm(date)
! 20180531 ADR: cofrm also calculates dipole quantities colat,elon,vp,ctp,stp
! write(6,"('apex_mka after cofrm: ncoef=',i4,' gb=',/,(6f12.3))") ncoef,gb
! write(6,"('apex_mka after cofrm: ncoef=',i4,' gv=',/,(6f12.3))") ncoef,gv
  reqore = req/re
  rqorm1 = reqore-1.
  do j=1,nlat
    ct = sin(gplat(j)*dtr)
    st = cos(gplat(j)*dtr)
    kpol = 0
    if (abs(gplat(j)) > pola) kpol = 1
    do i=1,nlon
      cp = cos((gplon(i)-elon)*dtr)
      sp = sin((gplon(i)-elon)*dtr)
!
!  ctm   is pseudodipole component of z
! -ctm   is pseudodipole component of v
!  stmcpm is pseudodipole component of x
!  stmspm is pseudodipole component of y
!
      ctm = ctp*ct + stp*st*cp
      stmcpm = st*ctp*cp - ct*stp
      stmspm = st*sp
      do k=1,nalt
        call apex_sub(gplat(j),gplon(i),gpalt(k),&
          aht,alat,phia,bmag,xmag,ymag,zdown,vmp)
        vnor = vmp/vp
        rp = 1. + gpalt(k)/re
        varray(j,i,k) = vnor*rp*rp + ctm
        reqam1 = req*(aht-1.)
        slp = sqrt(max(reqam1-gpalt(k),0.)/(reqam1+re))
!
! Reverse sign of slp in southern magnetic hemisphere
!
        if (zdown.lt.0.) slp = -slp
        clp = sqrt (rp/(reqore*aht-rqorm1))
        phiar = phia*dtr
        xarray(j,i,k) = clp*cos (phiar) - stmcpm
        yarray(j,i,k) = clp*sin (phiar) - stmspm
        zarray(j,i,k) = slp - ctm
      enddo ! k=1,nalt
      if (kpol==1.and.i > 1) then
        xarray(j,i,:) = xarray(j,1,:)
        yarray(j,i,:) = yarray(j,1,:)
        zarray(j,i,:) = zarray(j,1,:)
        varray(j,i,:) = varray(j,1,:)
        cycle
      endif
    enddo ! i=1,nlon
  enddo ! j=1,nlat
!
! Establish for this grid polar latitude limits beyond which east-west
! gradients are computed differently to avoid potential underflow
! (glatmx,glatmn are in module data, glatlim is parameter constant)
!
  glatmx = max( glatlim,gplat(nlat-2))
  glatmn = min(-glatlim,gplat(2))
end subroutine apex_mka
!-----------------------------------------------------------------------
subroutine apex_mall(glat,glon,alt,hr, b,bhat,bmag,si,alon,xlatm,vmp,W,&
  D,be3,sim,d1,d2,d3,e1,e2,e3,xlatqd,F,f1,f2,f3,g1,g2,g3,ier)
!
! Compute Modified Apex coordinates, quasi-dipole coordinates,
! base vectors and other parameters by interpolation from
! precalculated arrays. Subroutine apex_setup or apex_mka must be called
! before calling this subroutine.
!
! Args:
  real(rp),intent(in)  :: & ! Input
    glat             ,& ! Geographic (geodetic) latitude (deg)
    glon             ,& ! Geographic (geodetic) longitude (deg)
    alt              ,& ! Altitude (km)
    hr                  ! Reference altitude (km)
  real(rp),intent(out) :: & ! Output
    b(3)             ,& ! Magnetic field components (east, north, up), in nT
    bhat(3)          ,& ! components (east, north, up) of unit vector along
                        ! geomagnetic field direction
    bmag             ,& ! Magnitude of magnetic field (nT)
    si               ,& ! sine of magnetic inclination
    alon             ,& ! Apex longitude = modified apex longitude =
                        ! quasi-dipole longitude (deg)
    xlatm            ,& ! Modified Apex latitude (deg)
    vmp              ,& ! Magnetic potential (T.m)
    W                ,& ! W of Richmond [1995] reference above, in km**2 /nT
                        !  (i.e., 10**15 m**2 /T)
    D                ,& ! D of Richmond [1995]
    be3              ,& ! B_e3 of Richmond [1995] (= Bmag/D), in nT
    sim              ,& ! sin(I_m) described in Richmond [1995]
    xlatqd           ,& ! Quasi-dipole latitude (deg)
    F                   ! F described in Richmond [1995] for quasi-dipole
                        !  coordinates
!
   real(rp),dimension(3),intent(out) :: d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3 ! Components of base vectors
   integer,intent(out) :: ier ! error return
!
! Local:
  real(rp) :: glonloc,cth,sth,glatx,clm,r3_2
  real(rp) :: fx,fy,fz,fv
  real(rp) :: dfxdth,dfydth,dfzdth,dfvdth, &
          dfxdln,dfydln,dfzdln,dfvdln, &
          dfxdh ,dfydh ,dfzdh ,dfvdh
  real(rp),dimension(3) :: gradx,grady,gradz,gradv, grclm,clmgrp,rgrlp
  real(rp) ::                        & ! dummies for polar calls to intrp
    fxdum,fydum,fzdum,fvdum,     &
    dmxdth,dmydth,dmzdth,dmvdth, &
    dmxdh,dmydh,dmzdh,dmvdh
!
! Init:
!
  ier = 0
  glonloc = glon
  call intrp (glat,glonloc,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
             fx,fy,fz,fv,                                               &
             dfxdth,dfydth,dfzdth,dfvdth,                               &
             dfxdln,dfydln,dfzdln,dfvdln,                               &
             dfxdh ,dfydh ,dfzdh ,dfvdh, ier)
  if (ier /= 0) then
    call setmiss(xmiss,xlatm,alon,vmp,b,bmag,be3,sim,si,F,D,W, &
      bhat,d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3)
    write(6,"('apex_mall called setmiss: glat,glon,alt=',3f12.3)") &
      glat,glon,alt
    return
  endif
  call adpl(glat,glonloc,cth,sth,fx,fy,fz,fv, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
  call gradxyzv(alt,cth,sth, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)
!
! If the point is very close to either the North or South
! geographic pole, recompute the east-west gradients after
! stepping a small distance from the pole.
!
  if (glat > glatmx .or. glat < glatmn) then
    glatx = glatmx
    if (glat < 0.) glatx = glatmn
    call intrp (glatx,glonloc,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
               fxdum,fydum,fzdum,fvdum,                                    &
               dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,           &
               dfvdln,dmxdh,dmydh,dmzdh,dmvdh, ier)
    call adpl(glatx,glonloc,cth,sth,fxdum,fydum,fzdum,fvdum, &
      dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,dfvdln)
    call grapxyzv(alt,cth,sth,dfxdln,dfydln,dfzdln,dfvdln, &
      gradx,grady,gradz,gradv)
  endif
  call gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
    xlatm,alon,vmp,grclm,clmgrp,xlatqd,rgrlp,b,clm,r3_2)
  call basevec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
     bmag,sim,si,F,D,W,bhat,d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3)
  be3 = bmag/D
  ier = 0
end subroutine apex_mall
!-----------------------------------------------------------------------
subroutine apex_q2g(qdlat_in,qdlon,alt,gdlat,gdlon,ier)
!
! Convert from quasi-dipole to geodetic coordinates. This subroutine
! (input magnetic, output geodetic) is the functional inverse of
! subroutine apex_mall (input geodetic, output magnetic). Sub apex_setup
!  or apex_mka must be called before this routine.
!
! 20160811 ADR If qdlat_in is +/- 90, move it a little off the pole so
!  that when gdlat, gdlon are used in apxmall the directions of the base
!  vectors vary appropriately with longitude.
!
! Args:
  real(rp),intent(in) ::  & ! inputs
    qdlat_in,         & ! quasi-dipole latitude (deg)
    qdlon,            & ! quasi-dipole longitude (deg)
    alt                 ! altitude (km)
  real(rp),intent(out) :: & ! outputs
    gdlat,            & ! geodetic latitude (deg)
    gdlon               ! geodetic longitude (deg)
  integer,intent(out) :: ier ! error return
!
! Local:
  real(rp) :: x0,y0,z0,xnorm,xdif,ydif,zdif,dist2,hgrd2e,hgrd2n,hgrd2,&
    angdist,distlon,glatx,cal,sal,coslm,slm,cad,sad,slp,clm2,slm2,&
    sad2,cal2,clp2,clp,dylon,qdlat
  real(rp) :: ylat,ylon ! first guess output by gm2gc, input to intrp
  integer :: iter
  integer,parameter :: niter=20
  real(rp) ::                        & ! output of sub intrp
    fx,fy,fz,fv,                 & ! interpolated values of x,y,z,v
    dfxdth,dfydth,dfzdth,dfvdth, & ! derivatives of x,y,z,v wrt colatitude
    dfxdln,dfydln,dfzdln,dfvdln, & ! derivatives of x,y,z,v wrt longitude
    dfxdh ,dfydh ,dfzdh ,dfvdh     ! derivatives of x,y,z,v wrt altitude
  real(rp) ::                        & ! dummies for polar calls to intrp
    fxdum,fydum,fzdum,fvdum,     &
    dmxdth,dmydth,dmzdth,dmvdth, &
    dmxdh,dmydh,dmzdh,dmvdh
  real(rp) :: cth,sth  ! output of adpl
  character(len=5) :: edge
  ier = 0 ; gdlat = 0. ; gdlon = 0.
! Keep qdlat away from poles
  qdlat = max(qdlat_in,-90.+sqrt(precise)*rtd)
  qdlat = min(qdlat   , 90.-sqrt(precise)*rtd)
!
! Determine quasi-cartesian coordinates on a unit sphere of the
! desired magnetic lat,lon in quasi-dipole coordinates.
!
  x0 = cos (qdlat*dtr) * cos (qdlon*dtr)
  y0 = cos (qdlat*dtr) * sin (qdlon*dtr)
  z0 = sin (qdlat*dtr)
!
! Initial guess:  use centered dipole, convert to geocentric coords
!
  call gm2gc (qdlat,qdlon,ylat,ylon)
!
! Iterate until (angular distance)**2 (units: radians) is within
! precise of location (qdlat,qdlon) on a unit sphere.
! (precise is a parameter in module data)
!
  do iter=1,niter
!
! geolat,lon,alt and nglat,lon,alt are in module data (set by apex_mka)
!
    call intrp (ylat,ylon,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
               fx,fy,fz,fv,                                            &
               dfxdth,dfydth,dfzdth,dfvdth,                            &
               dfxdln,dfydln,dfzdln,dfvdln,                            &
               dfxdh ,dfydh ,dfzdh ,dfvdh, ier)
    if (ier /= 0) then
      write(6,"('>>> apex_q2g error from intrp')")
      stop 'qpex_q2g intrp'
    endif
!
!  Add-back of pseudodipole component to x,y,z,v and their derivatives.
!
    call adpl(ylat,ylon,cth,sth,fx,fy,fz,fv, &
      dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
    distlon = cos(ylat*dtr)
    if (ylat > glatmx .or. ylat < glatmn) then ! glatmx,glatmn are module data
      glatx = glatmx
      if (ylat.lt.0.) glatx = glatmn
      distlon = cos (glatx*dtr)
      call intrp (glatx,ylon,alt, geolat,geolon,geoalt,nglat,nglon,ngalt, &
                 fxdum,fydum,fzdum,fvdum,                                 &
                 dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,        &
                 dfvdln,dmxdh,dmydh,dmzdh,dmvdh, ier)
      if (ier /= 0) then
        write(6,"('>>> apex_q2g error from polar intrp')")
        stop 'qpex_q2g intrp'
      endif
      call adpl(glatx,ylon,cth,sth,fxdum,fydum,fzdum,fvdum, &
        dmxdth,dmydth,dmzdth,dmvdth,dfxdln,dfydln,dfzdln,dfvdln)
    endif
!
! At this point, fx,fy,fz are approximate quasi-cartesian
! coordinates on a unit sphere for the quasi-dipole coordinates
! corresponding to the geodetic coordinates ylat, ylon.
! Normalize the vector length of (fx,fy,fz) to unity using xnorm
! so that the resultant vector can be directly compared with the
! target vector (x0,y0,z0).
!
    xnorm = sqrt(fx*fx + fy*fy + fz*fz)
    xdif = fx/xnorm - x0
    ydif = fy/xnorm - y0
    zdif = fz/xnorm - z0
!
! dist2 = square of distance between normalized (fx,fy,fz) and x0,y0,z0.
!
    dist2 = xdif*xdif + ydif*ydif + zdif*zdif
    if (dist2 <= precise) then
      ier = 0
      gdlat = ylat
      gdlon = ylon
      return
    endif
!
! hgrd2* = one-half of east or north gradient of dist2 on unit sphere.
!
    hgrd2e =  (xdif*dfxdln + ydif*dfydln + zdif*dfzdln)/distlon
    hgrd2n = -(xdif*dfxdth + ydif*dfydth + zdif*dfzdth)
    hgrd2  = sqrt(hgrd2e*hgrd2e + hgrd2n*hgrd2n)
!
! angdist = magnitude of angular distance to be moved for new guess
!           of ylat, ylon.
!
    angdist = dist2/hgrd2
!
! Following spherical trigonometry moves ylat,ylon to new location,
! in direction of grad(dist2), by amount angdist.
!
    cal = -hgrd2n/hgrd2
    sal = -hgrd2e/hgrd2
    coslm = cos(ylat*dtr)
    slm = sin(ylat*dtr)
    cad = cos(angdist)
    sad = sin(angdist)
    slp = slm*cad + coslm*sad*cal
    clm2 = coslm*coslm
    slm2 = slm*slm
    sad2 = sad*sad
    cal2 = cal*caL
    clp2 = clm2 + slm2*sad2 - 2.*slm*cad*coslm*sad*cal -clm2*sad2*cal2
    clp = sqrt (max(0.,clp2))
    ylat = atan2(slp,clp)*rtd
!
! Restrict latitude iterations to stay within the interpolation grid
! limits, but let intrp find any longitude exceedence.  This is only
! an issue when the interpolation grid does not cover the entire
! magnetic pole region.
!
    ylat = min(ylat,geolat(nglat))
    ylat = max(ylat,geolat(1))
    dylon = atan2 (sad*sal,cad*coslm-sad*slm*cal)*rtd
    ylon = ylon + dylon
    if (ylon > geolon(nglon)) ylon = ylon - 360.
    if (ylon < geolon(1))     ylon = ylon + 360.
  enddo ! iter=1,niter
  write(6,"('>>> apex_q2g: ',i3,' iterations only reduced the angular')") niter
  write(6,"('              difference to ',f10.5,' degrees, where test criterion')") &
    sqrt(dist2)*rtd
  write(6,"('              is ',f10.5,' degrees.')") sqrt(precise)*rtd
  edge = '     '
  if (ylat == geolat(nglat)) edge = 'north'
  if (ylat == geolat(1))     edge = 'south'
  if (edge /= '     ') then
    write(6,"('Coordinates are on the ',a,' edge of the interpolation grid ')") edge
    write(6,"('and latitude is constrained to stay within grid limits when iterating.')")
  endif
  ier = 1
end subroutine apex_q2g
!-----------------------------------------------------------------------
subroutine gradxyzv(alt,cth,sth, &
    dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh,gradx,grady,gradz,gradv)
!
! Calculates east,north,up components of gradients of x,y,z,v in
! geodetic coordinates.  All gradients are in inverse km.
! 940803 A. D. Richmond; updated 180510
!
! Args:
  real(rp),intent(in) :: alt,cth,sth
  real(rp),dimension(3),intent(out) :: gradx,grady,gradz,gradv
  real(rp),intent(in) ::             &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh,dfydh,dfzdh,dfvdh
!
! Local:
  real(rp) :: d,d2,rho,dddthod,drhodth,dzetdth,ddisdth
!
  d2 = req2 - req2e2*cth*cth
  d = sqrt(d2)
  rho = sth*(alt + req2/d)
  dddthod = req2e2*cth*sth/d2
  drhodth = alt*cth + (req2/d)*(cth-sth*dddthod)
  dzetdth =-alt*sth - ((req2-req2e2)/d)*(sth+cth*dddthod)
  ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)
  gradx(1) = dfxdln/rho
  grady(1) = dfydln/rho
  gradz(1) = dfzdln/rho
  gradv(1) = dfvdln/rho
  gradx(2) = -dfxdth/ddisdth
  grady(2) = -dfydth/ddisdth
  gradz(2) = -dfzdth/ddisdth
  gradv(2) = -dfvdth/ddisdth
  gradx(3) = dfxdh
  grady(3) = dfydh
  gradz(3) = dfzdh
  gradv(3) = dfvdh
end subroutine gradxyzv
!-----------------------------------------------------------------------
subroutine grapxyzv(alt,cth,sth, &
  dfxdln,dfydln,dfzdln,dfvdln,gradx,grady,gradz,gradv)
!
! Calculates east component of gradient near pole.
!
! Args:
  real(rp),intent(in) :: alt,cth,sth
  real(rp),intent(in) :: dfxdln,dfydln,dfzdln,dfvdln
  real(rp),dimension(3),intent(out) :: gradx,grady,gradz,gradv
!
! Local:
  real(rp) :: d,d2,rho,dddthod,drhodth,dzetdth,ddisdth
!
  d2 = req2 - req2e2*cth*cth
  d = sqrt(d2)
  rho = sth*(alt + req2/d)
  dddthod = req2e2*cth*sth/d2
  drhodth = alt*cth + (req2/d)*(cth-sth*dddthod)
  dzetdth =-alt*sth - ((req2-req2e2)/d)*(sth+cth*dddthod)
  ddisdth = sqrt(drhodth*drhodth + dzetdth*dzetdth)
  gradx(1) = dfxdln/rho
  grady(1) = dfydln/rho
  gradz(1) = dfzdln/rho
  gradv(1) = dfvdln/rho
end subroutine grapxyzv
!-----------------------------------------------------------------------
subroutine gradlpv(hr,alt,fx,fy,fz,fv,gradx,grady,gradz,gradv, &
  xlatm,xlonm,vmp,grclm,clmgrp,qdlat,rgrlp,b,clm,r3_2)
!
! Uses gradients of x,y,z,v to compute geomagnetic field and
! gradients of apex latitude, longitude.
!
! Args:
  real(rp),intent(in) :: & ! scalar inputs
    hr,              & ! reference altitude (km)
    alt,             & ! altitude (km)
    fx,fy,fz,fv        ! interpolated values of x,y,z,v, plus
                       ! pseudodipole component
  real(rp),dimension(3),intent(in) :: & ! 3-component inputs
    gradx,grady,gradz,gradv ! interpolated gradients of x,y,z,v,
                            ! including pseudodipole components (east,north,up)
!
! Local:
  integer :: i
  real(rp) :: rr,r,rn,sqrror,cpm,spm,bo,rn2,x2py2,xnorm,xlp,slp,clp,grclp
  real(rp),intent(out) :: & ! scalar outputs
    xlatm,  & !  modified apex latitude (lambda_m), degrees
    xlonm,  & !  apex longitude (phi_a), degrees
    vmp,    & !  magnetic potential, in T.m.
    qdlat,  & !  quasi-dipole latitude, degrees
    clm,    & !  cos(lambda_m)
    r3_2      !  ((re + alt)/(re + hr))**(3/2)
  real(rp),dimension(3),intent(out) :: & ! 3-component outputs
    grclm,   & ! grad(cos(lambda_m)), in km-1
    clmgrp,  & ! cos(lambda_m)*grad(phi_a), in km-1
    rgrlp,   & ! (re + alt)*grad(lambda')
    b          ! magnetic field, in nT
  xlatm=0. ; xlonm=0. ; vmp=0 ; grclm=0. ; clmgrp=0. ; rgrlp = 0. ; b=0.
  clm=0. ; r3_2=0. ; qdlat=0.
  rr = re + hr
  r  = re + alt
  rn = r/re
  sqrror = sqrt(rr/r)
  r3_2 = 1./sqrror/sqrror/sqrror
  xlonm = atan2(fy,fx)
  cpm = cos(xlonm)
  spm = sin(xlonm)
  xlonm = rtd*xlonm ! output
  bo = vp*1.e6 ! vp is module data; 1.e6 converts T to nT and km-1 to m-1
  rn2 = rn*rn
  vmp = vp*fv/rn2   ! output
  b(1) = -bo*gradv(1)/rn2
  b(2) = -bo*gradv(2)/rn2
  b(3) = -bo*(gradv(3)-2.*fv/r)/rn2
  x2py2 = fx*fx + fy*fy
  xnorm = sqrt(x2py2 + fz*fz)
  xlp = atan2(fz,sqrt(x2py2))
  slp = sin(xlp)
  clp = cos(xlp)
  qdlat = xlp*rtd   ! output
  clm = sqrror*clp  ! output
  if (clm > 1.) then
    write(6,"('>>> gradlpv: hr=',f12.3,' alt=',f12.3)") hr,alt
    write(6,"('    Point lies below field line that peaks at reference height.')")
    stop 'gradlpv'
  endif
  xlatm = rtd*acos(clm)
!
!  If southern magnetic hemisphere, reverse sign of xlatm
!
  if (slp < 0.) xlatm = -xlatm
  do i=1,3
    grclp = cpm*gradx(i) + spm*grady(i)
    rgrlp(i) = r*(clp*gradz(i) - slp*grclp)
    grclm(i) = sqrror*grclp
    clmgrp(i) = sqrror*(cpm*grady(i)-spm*gradx(i))
  enddo
  grclm(3) = grclm(3) - sqrror*clp/(2.*r)
end subroutine gradlpv
!-----------------------------------------------------------------------
subroutine basevec(hr,xlatm,grclm,clmgrp,rgrlp,b,clm,r3_2, &
     bmag,sim,si,F,D,W,bhat,d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3)
!
! Computes base vectors and other parameters for apex coordinates.
! Vector components:  east, north, up
!
! Args:
  real(rp),intent(in) :: & ! scalar inputs
    hr,      & ! reference altitude
    xlatm,   & ! modified apex latitude (deg)
    clm,     & ! cos(lambda_m)
    r3_2       ! ((re + altitude)/(re + hr))**(3/2)
  real(rp),dimension(3),intent(in) :: & ! 3-component inputs
    grclm,   & ! grad(cos(lambda_m)), in km-1
    clmgrp,  & ! cos(lambda_m)*grad(phi_a), in km-1
    rgrlp,   & ! (re + altitude)*grad(lambda')
    b          !  magnetic field, in nT
  real(rp),intent(out) :: & ! scalar output
    bmag,    & ! magnitude of magnetic field, in nT
    sim,     & ! sin(I_m) of Richmond reference
    si,      & ! sin(I)
    F,       & ! F of Richmond reference
    D,       & ! D of Richmond reference
    W          ! W of Richmond reference
  real(rp),dimension(3),intent(out) :: & ! 3-component outputs
    bhat,             & ! unit vector along geomagnetic field direction
    d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3   ! base vectors of Richmond reference
!
! Local:
  integer :: i
  real(rp) :: rr,simoslm,d1db,d2db
  rr = re + hr
  simoslm = 2./sqrt(4. - 3.*clm*clm)
  sim = simoslm*sin(xlatm*dtr)
  bmag = sqrt(b(1)*b(1) + b(2)*b(2) + b(3)*b(3))
  d1db = 0.
  d2db = 0.
  do i=1,3
    bhat(i) = b(i)/bmag
    d1(i) = rr*clmgrp(i)
    d1db = d1db + d1(i)*bhat(i)
    d2(i) = rr*simoslm*grclm(i)
    d2db = d2db + d2(i)*bhat(i)
  enddo
!
! Ensure that d1,d2 are exactly perpendicular to B:
!
  do i=1,3
    d1(i) = d1(i) - d1db*bhat(i)
    d2(i) = d2(i) - d2db*bhat(i)
  enddo
  e3(1) = d1(2)*d2(3) - d1(3)*d2(2)
  e3(2) = d1(3)*d2(1) - d1(1)*d2(3)
  e3(3) = d1(1)*d2(2) - d1(2)*d2(1)
  D = bhat(1)*e3(1) + bhat(2)*e3(2) + bhat(3)*e3(3)
  do i=1,3
    d3(i) = bhat(i)/D
    e3(i) = bhat(i)*D ! Ensure that e3 lies along bhat.
  enddo
  e1(1) = d2(2)*d3(3) - d2(3)*d3(2)
  e1(2) = d2(3)*d3(1) - d2(1)*d3(3)
  e1(3) = d2(1)*d3(2) - d2(2)*d3(1)
  e2(1) = d3(2)*d1(3) - d3(3)*d1(2)
  e2(2) = d3(3)*d1(1) - d3(1)*d1(3)
  e2(3) = d3(1)*d1(2) - d3(2)*d1(1)
  W = rr*rr*clm*abs(sim)/(bmag*D)
  si = -bhat(3)
  f1(1) =  rgrlp(2)
  f1(2) = -rgrlp(1)
  f2(1) = -d1(2)*r3_2
  f2(2) =  d1(1)*r3_2
  F = f1(1)*f2(2) - f1(2)*f2(1)
! Added output 2015 February 17
  f1(3) = 0.
  f2(3) = 0.
  g1(1) = r3_2*d1(1)/F
  g1(2) = r3_2*d1(2)/F
  g1(3) = r3_2*d1(3)/F
  g2(1) = rgrlp(1)/F
  g2(2) = rgrlp(2)/F
  g2(3) = rgrlp(3)/F
  g3(1) = 0.
  g3(2) = 0.
  g3(3) = F
  f3(1) = g1(2)*g2(3) - g1(3)*g2(2)
  f3(2) = g1(3)*g2(1) - g1(1)*g2(3)
  f3(3) = g1(1)*g2(2) - g1(2)*g2(1)
end subroutine basevec
!-----------------------------------------------------------------------
subroutine apex_sub(dlat,dlon,alt,aht,alat,alon,bmag,xmag,ymag,zmag,vmp)
!
!          Calculate apex radius, latitude, longitude; and magnetic field and
!          scalar magnetic potential.
!
!          INPUTS:
!            dlat = Geodetic latitude in degrees
!            dlon = Geodetic longitude in degrees
!            alt = Altitude in km
!
!          RETURNS:
!            aht  = (Apex height + req)/req, where req = equatorial Earth radius
!                   aht is analogous to the L value in invariant coordinates.
!            alat = Apex latitude in degrees (negative in S magnetic hemisphere)
!            alon = Apex longitude (geomagnetic longitude of apex) in degrees
!            bmag = geomagnetic field magnitude (nT)
!            xmag = geomagnetic field component (nT): north
!            ymag = geomagnetic field component (nT): east
!            zmag = geomagnetic field component (nT): downward
!            vmp  = geomagnetic potential (T-m)
!
! Args:
  real(rp),intent(inout) :: dlat,dlon,alt
  real(rp),intent(out) :: aht,alat,alon,bmag,xmag,ymag,zmag,vmp
!
! Local:
  real(rp) :: x,y,z,xre,yre,zre
  integer :: iflag
  vmp = 0.
!
! Last 7 args of linapx are output:
!
  call linapx(dlat,dlon,alt, aht,alat,alon,xmag,ymag,zmag,bmag)
  xmag = xmag*1.e5
  ymag = ymag*1.e5
  zmag = zmag*1.e5
  bmag = bmag*1.e5
  call gd2cart (dlat,dlon,alt,x,y,z)
  iflag = 3
  xre = x/re ; yre = y/re ; zre = z/re
  call feldg(iflag,xre,yre,zre,bx,by,bz,vmp)
end subroutine apex_sub
!-----------------------------------------------------------------------
subroutine linapx(gdlat,glon,alt,aht,alat,alon,xmag,ymag,zmag,fmag)
!
! Transform geographic coordinates to Apex coordinates.
! Based on code originally written in 1973 by Wally Clark at NOAA
! Modified by A.D. Richmond, 1994, to return meaningful values when the
!  maximum number of steps has been reached, by treating the field line
!  above that high altitude as dipolar.
! Rewritten in Fortran 90, eliminating common blocks, by B. Foster, 2013
! Reference: Stassinopoulos E.G., Mead Gilbert D., X-841-72-17 (1971),
!  NASA GSFC, Greenbelt, Maryland
!
! Input Args:
!
  real(rp),intent(inout) :: & ! These may be changed by convrt, depending on iflag
    gdlat,              & ! latitude of starting point (deg)
    glon,               & ! longitude of starting point (deg)
    alt                   ! height of starting point (km)
!
! Output Args:
!
  real(rp),intent(out) :: &
    aht,              & ! (Apex height+req)/req, where req is equatorial earth radius
    alat,             & ! Apex latitude (deg)
    alon,             & ! Apex longitude (deg)
    xmag,             & ! North component of magnetic field at starting point
    ymag,             & ! East component of magnetic field at starting point
    zmag,             & ! Down component of magnetic field at starting point
    fmag                ! Magnetic field magnitude at starting point
!
! Trace the geomagnetic field line from the given location to find the
! apex of the field line.  Before starting iterations to trace along
! the field line: (1) Establish a step size (ds, arc length in km)
! based on the geomagnetic dipole latitude; (2) determine the step
! direction from the sign of the vertical component of the geomagnetic
! field; and (3) convert to geocentric cartesian coordinates.  Each
! iteration increments a step count (nstp) and calls itrace to move
! along the the field line until reaching the iteration count limit
! (maxs) or passing the apex (iapx=2) and then calling fndapx to
! determine the apex location from the last three step locations
! (yapx); however, if reaching the iteration limit, apex coordinates
! are calculated by dipapx which assumes a simplified dipole field.
!
! Local:
!
  real(rp) :: gclat,r,singml,cgml2,rho,xlat,xlon,ht
  real(rp) :: bnrth,beast,bdown,babs,y1,y2,y3
  integer :: iflag,iapx
  integer,parameter :: maxs = 200
!
! Determine step size as a function of geomagnetic dipole
! coordinates of the starting point
!
  iflag = 2 ! gclat,r are returned
  call convrt(iflag,gdlat,alt,gclat,r,'linapx')
  singml = ctp*sin(gclat*dtr) + stp*cos(gclat*dtr)*cos((glon-elon)*dtr)
  cgml2 = max(0.25,1.-singml*singml)
  ds = .06*r/cgml2 - 370. ! ds is in module data
  yapx = 0. ! init (module data)
!
! Convert from geodetic to earth centered cartesian coordinates:
!
  call gd2cart(gdlat,glon,alt,y(1),y(2),y(3))
  nstp = 0
!
! Get magnetic field components to determine the direction for tracing field line:
!
  iflag = 1
  call feldg(iflag,gdlat,glon,alt,xmag,ymag,zmag,fmag)
  sgn = sign(1.,-zmag)
!
! Use cartesian coordinates to get magnetic field components
! (from which gradients steer the tracing)
!
100 continue
  iflag = 2 ! module data bx,by,bz,bb are returned
  y1 = y(1)/re ; y2 = y(2)/re ; y3 = y(3)/re
  call feldg(iflag,y1,y2,y3,bx,by,bz,bb)
  nstp = nstp + 1
!
! Quit if too many steps.
!
  if (nstp >= maxs) then
    rho = sqrt(y(1)*y(1) + y(2)*y(2))
    iflag = 3 ! xlat and ht are returned
    call convrt(iflag,xlat,ht,rho,y(3),'linapx')
    xlon = rtd*atan2(y(2),y(1))
    iflag = 1
    call feldg(iflag,xlat,xlon,ht,bnrth,beast,bdown,babs)
    call dipapx(xlat,xlon,ht,bnrth,beast,bdown,aht,alon)
    alat = -sgn*rtd*acos(sqrt(1./aht))
    return
  endif
!
! Find next point using Adams algorithm after 7 points
!
  call itrace(iapx)
  if (iapx == 1) goto 100
!
! Maximum radius just passed.  Find apex coords
!
  call fndapx(alt,zmag,aht,alat,alon)
end subroutine linapx
!-----------------------------------------------------------------------
subroutine convrt(iflag,gdlat,alt,x1,x2,caller)
!
! Convert space point from geodetic to geocentric or vice versa.
!
! iflag = 1: Convert from geodetic to cylindrical
!   Input:  gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!   Output: x1    = Distance from Earth's rotation axis (km)
!           x2    = Distance above (north of) Earth's equatorial plane (km)
!
! iflag = 2: Convert from geodetic to geocentric spherical
!   Input:  gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!   Output: x1    = Geocentric latitude (deg)
!           x2    = Geocentric distance (km)
!
! iflag = 3: Convert from cylindrical to geodetic
!   Input:  x1    = Distance from Earth's rotation axis (km)
!           x2    = Distance from Earth's equatorial plane (km)
!   Output: gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!
! iflag = 4: Convert from geocentric spherical to geodetic
!   Input:  x1    = Geocentric latitude (deg)
!           x2    = Geocentric distance (km)
!   Output: gdlat = Geodetic latitude (deg)
!           alt   = Altitude above reference ellipsoid (km)
!
! Args:
  integer,intent(in) :: iflag
  real(rp),intent(inout) :: gdlat,alt
  real(rp),intent(inout) :: x1,x2
  character(len=*),intent(in) :: caller
!
! Local:
  real(rp) :: sinlat,coslat,d,z,rho,rkm,scl,gclat,ri,a2,a4,a6,a8,&
    ccl,s2cl,c2cl,s4cl,c4cl,s8cl,s6cl,dltcl,sgl
  real(rp),parameter ::                                         &
    e2=(2.-1./fltnvrs)/fltnvrs                            , &
    e4=e2*e2, e6=e4*e2, e8=e4*e4                          , &
    ome2req = (1.-e2)*req                                 , &
    a21 =     (512.*e2 + 128.*e4 + 60.*e6 + 35.*e8)/1024. , &
    a22 =     (                        e6 +     e8)/  32. , &
    a23 = -3.*(                     4.*e6 +  3.*e8)/ 256. , &
    a41 =    -(           64.*e4 + 48.*e6 + 35.*e8)/1024. , &
    a42 =     (            4.*e4 +  2.*e6 +     e8)/  16. , &
    a43 =                                   15.*e8 / 256. , &
    a44 =                                      -e8 /  16. , &
    a61 =  3.*(                     4.*e6 +  5.*e8)/1024. , &
    a62 = -3.*(                        e6 +     e8)/  32. , &
    a63 = 35.*(                     4.*e6 +  3.*e8)/ 768. , &
    a81 =                                   -5.*e8 /2048. , &
    a82 =                                   64.*e8 /2048. , &
    a83 =                                 -252.*e8 /2048. , &
    a84 =                                  320.*e8 /2048.
  if (iflag < 3) then ! geodetic to geocentric
!
! Compute rho,z
    sinlat = sin(gdlat*dtr)
    coslat = sqrt(1.-sinlat*sinlat)
    d      = sqrt(1.-e2*sinlat*sinlat)
    z      = (alt+ome2req/d)*sinlat
    rho    = (alt+req/d)*coslat
    x1 = rho
    x2 = z
    if (iflag == 1) return
!
! Compute gclat,rkm
    rkm   = sqrt(z*z+rho*rho)
    gclat = rtd*atan2(z,rho)
    x1 = gclat
    x2 = rkm
    return    ! iflag == 2
  endif ! iflag < 3
!
! Geocentric to geodetic:
  if (iflag == 3) then
    rho = x1
    z   = x2
    rkm = sqrt(z*z+rho*rho)
    scl = z/rkm
    gclat = asin(scl)*rtd
  elseif (iflag == 4) then
    gclat = x1
    rkm = x2
    scl = sin(gclat*dtr)
  else
    return
  endif
!
! iflag == 3 or 4:
!
  ri = req/rkm
  a2 = ri*(a21+ri*(a22+ri* a23))
  a4 = ri*(a41+ri*(a42+ri*(a43+ri*a44)))
  a6 = ri*(a61+ri*(a62+ri* a63))
  a8 = ri*(a81+ri*(a82+ri*(a83+ri*a84)))
  ccl = sqrt(1.-scl*scl)
  s2cl = 2.*scl*ccL
  c2cl = 2.*ccl*ccl-1.
  s4cl = 2.*s2cl*c2cl
  c4cl = 2.*c2cl*c2cl-1.
  s8cl = 2.*s4cl*c4cl
  s6cl = s2cl*c4cl+c2cl*s4cl
  dltcl = s2cl*a2+s4cl*a4+s6cl*a6+s8cl*a8
  gdlat = dltcl*rtd+gclat
  sgl = sin(gdlat*dtr)
  alt = rkm*cos(dltcl)-req*sqrt(1.-e2*sgl*sgl)
end subroutine convrt
!-----------------------------------------------------------------------
subroutine gd2cart(gdlat,glon,alt,x,y,z)
!
! Arg:
  real(rp),intent(inout) :: gdlat,alt,z
  real(rp),intent(in) :: glon
  real(rp),intent(out) :: x,y
!
! Local:
  real(rp) :: ang,rho
  integer :: iflag
  iflag = 1 ! Convert from geodetic to cylindrical (rho,z are output)
  call convrt(iflag,gdlat,alt,rho,z,'gd2cart')
  ang = glon*dtr
  x = rho*cos(ang)
  y = rho*sin(ang)
end subroutine gd2cart
!-----------------------------------------------------------------------
subroutine feldg(iflag,glat,glon,alt,bnrth,beast,bdown,babs)
!
! Compute the DGRF/IGRF field components at the point glat,glon,alt.
! cofrm must be called to establish coefficients for correct date
! prior to calling FELDG.
!
! iflag = 1:
!   Inputs:
!     glat = Latitude of point (deg)
!     glon = Longitude of point (deg)
!     alt  = Height of point (km)
!   Outputs:
!     bnrth = North component of field vector (Gauss)
!     beast = East component of field vector (Gauss)
!     bdown = Downward component of field vector (Gauss)
!     babs  = Magnitude of field vector (Gauss)
!
! iflag = 2:
!   Inputs:
!     glat = x coordinate (in units of earth radii 6371.2 km)
!     glon = y coordinate (in units of earth radii 6371.2 km)
!     alt  = z coordinate (in units of earth radii 6371.2 km)
!   Outputs:
!     bnrth = x component of field vector (Gauss)
!     beast = y component of field vector (Gauss)
!     bdown = z component of field vector (Gauss)
!     babs  = Magnitude of field vector (Gauss)
!
! iflag = 3:
!   Inputs:
!     glat = x coordinate (in units of earth radii 6371.2 km)
!     glon = y coordinate (in units of earth radii 6371.2 km)
!     alt  = z coordinate (in units of earth radii 6371.2 km)
!   Outputs:
!     bnrth = Dummy variable
!     beast = Dummy variable
!     babs  = magnetic potential in T.m
!
! Args:
  integer,intent(in) :: iflag
  real(rp),intent(in)    :: glon
  real(rp),intent(inout) :: glat,alt
  real(rp),intent(out)   :: bnrth,beast,bdown,babs
!
! Local:
  integer :: i,is,ihmax,last,imax,mk,k,ih,m,il,ihm,ilm
  real(rp) :: rlat,ct,st,rlon,cp,sp,xxx,yyy,zzz,rq,f,x,y,z
  real(rp) :: xi(3),h(ncoef),g(ncoef)
  real(rp) :: s,t,bxxx,byyy,bzzz,brho
  if (iflag == 1) then
    is   = 1
    rlat = glat*dtr
    ct   = sin(rlat)
    st   = cos(rlat)
    rlon = glon*dtr
    cp   = cos(rlon)
    sp   = sin(rlon)
    call gd2cart(glat,glon,alt,xxx,yyy,zzz)
    xxx = xxx/re
    yyy = yyy/re
    zzz = zzz/re
  else
    is  = 2
    xxx = glat
    yyy = glon
    zzz = alt
  endif
  rq    = 1./(xxx**2+yyy**2+zzz**2)
  xi(1) = xxx*rq
  xi(2) = yyy*rq
  xi(3) = zzz*rq
  ihmax = nmax*nmax+1
  last  = ihmax+nmax+nmax
  imax  = nmax+nmax-1
!
  if (iflag /= 3) then
    do i=1,last
      g(i) = gb(i) ! gb is module data from cofrm
    enddo
  else
    do i=1,last
      g(i) = gv(i) ! gv is module data from cofrm
    enddo
  endif
  do i=ihmax,last
    h(i) = g(i)
  enddo
  mk = 3
  if (imax == 1) mk = 1
  do k=1,mk,2
    i = imax
    ih = ihmax
100 continue
    il = ih-i
    f = 2./float(i-k+2)
    x = xi(1)*f
    y = xi(2)*f
    z = xi(3)*(f+f)
    i = i-2
    if (i < 1) then
      h(il) = g(il) + z*h(ih) + 2.*(x*h(ih+1)+y*h(ih+2))
    elseif (i == 1) then
      h(il+2) = g(il+2) + z*h(ih+2) + x*h(ih+4) - y*(h(ih+3)+h(ih))
      h(il+1) = g(il+1) + z*h(ih+1) + y*h(ih+4) + x*(h(ih+3)-h(ih))
      h(il)   = g(il)   + z*h(ih)   + 2.*(x*h(ih+1)+y*h(ih+2))
    else
      do m=3,i,2
        ihm = ih+m
        ilm = il+m
        h(ilm+1) = g(ilm+1)+ z*h(ihm+1) + x*(h(ihm+3)-h(ihm-1))- &
                   y*(h(ihm+2)+h(ihm-2))
        h(ilm)   = g(ilm)  + z*h(ihm)   + x*(h(ihm+2)-h(ihm-2))+ &
                   y*(h(ihm+3)+h(ihm-1))
      enddo
      h(il+2) = g(il+2) + z*h(ih+2) + x*h(ih+4) - y*(h(ih+3)+h(ih))
      h(il+1) = g(il+1) + z*h(ih+1) + y*h(ih+4) + x*(h(ih+3)-h(ih))
      h(il)   = g(il)   + z*h(ih)   + 2.*(x*h(ih+1)+y*h(ih+2))
    endif
    ih = il
    if (i >= k) goto 100
  enddo ! k=1,mk,2
  s = .5*h(1)+2.*(h(2)*xi(3)+h(3)*xi(1)+h(4)*xi(2))
  t = (rq+rq)*sqrt(rq)
  bxxx = t*(h(3)-s*xxx)
  byyy = t*(h(4)-s*yyy)
  bzzz = t*(h(2)-s*zzz)
  babs = sqrt(bxxx**2+byyy**2+bzzz**2)
  if (is .eq. 1) then            ! (convert back to geodetic)
    beast = byyy*cp-bxxx*sp
    brho  = byyy*sp+bxxx*cp
    bnrth = bzzz*st-brho*ct
    bdown = -bzzz*ct-brho*st
  elseif (is .eq. 2) then        ! (leave in earth centered cartesian)
    bnrth = bxxx
    beast = byyy
    bdown = bzzz
  endif
!
! Magnetic potential computation (V) makes use of the fact that the
! calculation of V is identical to that for r*Br, if coefficients
! in the latter calculation have been divided by (n+1) (coefficients
! gv).  Factor .1 converts km to m and gauss to tesla.
!
  if (iflag == 3) babs = (bxxx*xxx + byyy*yyy + bzzz*zzz)*re*.1
end subroutine feldg
!-----------------------------------------------------------------------
subroutine dipapx(gdlat,gdlon,alt,bnorth,beast,bdown,A,alon)
!
! Compute A, alon from local magnetic field using dipole and spherical approx.
! 940501 A. D. Richmond
!
! Input:
!   gdlat  = geodetic latitude, degrees
!   gdlon  = geodetic longitude, degrees
!   alt    = altitude, km
!   bnorth = geodetic northward magnetic field component (any units)
!   beast  = eastward magnetic field component
!   bdown  = geodetic downward magnetic field component
! Output:
!   A      = apex radius, 1 + h_A/R_eq
!   alon   = apex longitude, degrees
!
! Algorithm:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let G be point at gdlat,gdlon.
!   Let E be point on sphere below apex of dipolar field line passing through G.
!   Let td be dipole colatitude of point G, found by applying dipole formula
!     for dip angle to actual dip angle.
!   Let b be Pi plus local declination angle.  b is in the direction
!     from G to E.
!   Let tg be colatitude of G.
!   Let ang be longitude angle from GM to G.
!   Let te be colatitude of E.
!   Let tp be colatitude of GM.
!   Let a be longitude angle from G to E.
!   Let apang = a + ang
!   Let PA be geomagnetic longitude, i.e., Pi minus angle measured
!     counterclockwise from arc GM-E to arc GM-GP.
!   Let TF be arc length between GM and E.
!   Then, using notation c=cos, s=sin, cot=cot, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: stfcpa,
!     stfspa are sin(TF) times cos(PA), sin(PA), respectively.
!
  real(rp),intent(in)  :: gdlat,gdlon,alt,bnorth,beast,bdown
  real(rp),intent(out) :: a,alon
!
! Local:
  real(rp) :: bhor,std,ctd,sb,cb,ctg,stg,ang,sang,cang,cte,ste,sa,ca, &
          cottd,capang,sapang,stfcpa,stfspa,ha,r
  bhor = sqrt(bnorth*bnorth + beast*beast)
  if (bhor == 0.) then
! To avoid dividing by zero, alon and a are set here when bhor=0.
    alon = 0.
    a = 1.e34
    return
  endif
  cottd = bdown*.5/bhor
  std = 1./sqrt(1.+cottd*cottd)
  ctd = cottd*std
  sb = -beast/bhor
  cb = -bnorth/bhor
  ctg = sin(gdlat*dtr)
  stg = cos(gdlat*dtr)
  ang = (gdlon-elon)*dtr
  sang = sin(ang)
  cang = cos(ang)
  cte = ctg*std + stg*ctd*cb
  ste = sqrt(1. - cte*cte)
  sa = sb*ctd/ste
  ca = (std*stg - ctd*ctg*cb)/ste
  capang = ca*cang - sa*sang
  sapang = ca*sang + sa*cang
  stfcpa = ste*ctp*capang - cte*stp
  stfspa = sapang*ste
  alon = atan2(stfspa,stfcpa)*rtd
  r = alt + re
  ha = alt + r*cottd*cottd
  A = 1. + ha/req
end subroutine dipapx
!-----------------------------------------------------------------------
subroutine itrace(iapx)
! Based on code originally written in 1973 by Wally Clark at NOAA
! Rewritten in Fortran 90 by B. Foster, 2013, eliminating common blocks.
  save
!
! Uses 4-point ADAMS formula after initialization.
! First 7 iterations advance point by 3 steps.
!
! y(3), yp(3), yapx(3,3), sgn and nstp are in module data
! yploc(3,4) is local
!
! Arg:
  integer,intent(out) :: iapx
!
! Local:
  integer :: i,j
  real(rp) :: yploc(3,4) ! local yp (i.e., not module data yp)
  real(rp) :: term,d2,d6,d12,d24,rc,rp
  iapx = 1
!
! Field line is defined by the following differential equations
! in cartesian coordinates.
! (yapx,yp,y are module data)
!
  yploc(1,4) = sgn*bx/bb
  yploc(2,4) = sgn*by/bb
  yploc(3,4) = sgn*bz/bb
  if (nstp > 7) then
    do i=1,3
      yapx(i,1) = yapx(i,2)
      yapx(i,2) = y(i)
      yp(i) = y(i)
      term = 55.*yploc(i,4)-59.*yploc(i,3)+37.*yploc(i,2)-9.*yploc(i,1)
      y(i) = yp(i) + d24*term
      yapx(i,3) = y(i)
      do j=1,3
        yploc(i,j) = yploc(i,j+1)
      enddo
    enddo
    rc = rdus ( y(1),  y(2),  y(3))
    rp = rdus (yp(1), yp(2), yp(3))
    if (rc < rp) iapx=2
    return
  endif
  do i=1,3
    select case (nstp)
      case (1)
        d2  = ds/2.
        d6  = ds/6.
        d12 = ds/12.
        d24 = ds/24.
        yploc(i,1)= yploc(i,4)
        yp(i)     = y(i)
        yapx(i,1) = y(i)
        y(i) = yp(i) + ds*yploc(i,1)
      case (2)
        yploc(i,2) = yploc(i,4)
        y(i) = yp(i) + d2*(yploc(i,2)+yploc(i,1))
      case (3)
        y(i) = yp(i) + d6*(2.*yploc(i,4)+yploc(i,2)+3.*yploc(i,1))
      case (4)
        yploc(i,2) = yploc(i,4)
        yapx(i,2)  = y(i)
        yp(i)      = y(i)
        y(i)       = yp(i) + d2*(3.*yploc(i,2)-yploc(i,1))
      case (5)
        y(i) = yp(i) + d12*(5.*yploc(i,4)+8.*yploc(i,2)-yploc(i,1))
      case (6)
        yploc(i,3) = yploc(i,4)
        yp(i)      = y(i)
        yapx(i,3)  = y(i)
        y(i)       = yp(i) + d12*(23.*yploc(i,3)-16.*yploc(i,2)+5.*yploc(i,1))
      case (7)
        yapx(i,1) = yapx(i,2)
        yapx(i,2) = yapx(i,3)
        y(i) = yp(i) + d24*(9.*yploc(i,4)+19.*yploc(i,3)-5.*yploc(i,2)+yploc(i,1))
        yapx(i,3) = y(i)
      case default
        write(6,"('>>> itrace: unresolved case nstp=',i4)") nstp
        stop 'itrace'
    end select
  enddo
!
! Signal if apex passed:
!
  if (nstp == 6 .or. nstp == 7) then
    rc = rdus( yapx(1,3), yapx(2,3), yapx(3,3))
    rp = rdus( yapx(1,2), yapx(2,2), yapx(3,2))
    if (rc < rp) iapx=2
  endif
end subroutine itrace
!-----------------------------------------------------------------------
real function rdus(d,e,f)
  real(rp),intent(in) :: d,e,f
  rdus = sqrt(d**2 + e**2 + f**2)
end function rdus
!-----------------------------------------------------------------------
subroutine fndapx(alt,zmag,A,alat,alon)
!
! Find apex coords once tracing has signalled that the apex has been passed.
! Based on code written by Wally Clark, NOAA, 1973.
! Modified by A.D. Richmond, 1994.
! Modified by R. Barnes, 2004, to repair problem noted by Dan Weimer.
!  where the apex location produced by fint may still have a non-zero
!  vertical magnetic field component.
! Rewritten in Fortran 90 by B. Foster, 2013, eliminating common blocks.
!
! Args:
  real(rp),intent(in) ::  &
                 alt, & ! altitude in km
                 zmag   ! Z component of magnetic field (sign only is needed)
  real(rp),intent(out) :: &
                 A,   & ! apex radius
                 alat,& ! apex latitude
                 alon   ! apex longitude
!
! Local:
  integer :: i,iflag_convrt, iflag_feldg
  real(rp) :: z(3),ht(3),yloc(3),gdlt,gdln,x,ydum,f,rho,xinter,rasq,xlon,ang,&
    cang,sang,r,cte,ste,stfcpa,stfspa
!
! Get geodetic field components.
!
  iflag_feldg = 1
  iflag_convrt = 3
  do i=1,3
    rho = sqrt(yapx(1,i)**2+yapx(2,i)**2)
    call convrt(iflag_convrt,gdlt,ht(i),rho,yapx(3,i),'fndapx')
    gdln = rtd*atan2(yapx(2,i),yapx(1,i))
    call feldg(iflag_feldg,gdlt,gdln,ht(i),x,ydum,z(i),f)
  enddo
!
! Find cartesian coordinates at dip equator by interpolation
!
  do i=1,3
    call fint(z(1),z(2),z(3),yapx(i,1),yapx(i,2),yapx(i,3),0._rp,yloc(i))
  enddo
!
! Find apex height by interpolation
!
  call fint(z(1),z(2),z(3),ht(1),ht(2),ht(3),0._rp,xinter)
!
! Ensure that xinter is not less than original starting altitude:
  xinter = max(alt,xinter)
  A = (req+xinter)/req
!
! Find apex coordinates, giving alat the sign of dip at starting point.
! alon is the value of the geomagnetic longitude at the apex.
!
  if (A < 1.) then
    write(6,"('>>> fndapx: A=',e12.4,' < 1.')") A
    stop 'fndapx'
  endif
  rasq = rtd*acos(sqrt(1./A))
  alat = sign(rasq,zmag)
!
! Algorithm for alon:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude colat, east longitude elon).
!   Let xlon be longitude of apex.
!   Let te be colatitude of apex.
!   Let ang be longitude angle from GM to apex.
!   Let tp be colatitude of GM.
!   Let TF be arc length between GM and apex.
!   Let PA = be geomagnetic longitude in radians, i.e., Pi minus angle measured
!     counterclockwise from arc GM-apex to arc GM-GP.
!   Then, using notation c=cos, s=sin, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: stfcpa,
!     stfspa are sin(TF) times cos(PA), sin(PA), respectively.
!
  xlon = atan2(yloc(2),yloc(1))
  ang  = xlon-elon*dtr
  cang = cos(ang)
  sang = sin(ang)
  r    = sqrt(yloc(1)**2+yloc(2)**2+yloc(3)**2)
  cte  = yloc(3)/r
  ste  = sqrt(1.-cte*cte)
  stfcpa = ste*ctp*cang - cte*stp
  stfspa = sang*ste
  alon = atan2(stfspa,stfcpa)*rtd
end subroutine fndapx
!-----------------------------------------------------------------------
subroutine fint(a1,a2,a3,a4,a5,a6,a7,result)
!
! Second degree interpolation
!
! Args:
  real(rp),intent(in) :: a1,a2,a3,a4,a5,a6,a7
  real(rp),intent(out) :: result
  result = ((a2-a3)*(a7-a2)*(a7-a3)*a4-(a1-a3)*(a7-a1)*(a7-a3)*a5+ &
    (a1-a2)*(a7-a1)*(a7-a2)*a6)/((a1-a2)*(a1-a3)*(a2-a3))
end subroutine fint
!-----------------------------------------------------------------------
subroutine gm2gc(gmlat,gmlon,gclat,gclon)
!
! Args:
  real(rp),intent(in)  :: gmlat,gmlon ! geomagnetic dipole lat,lon
  real(rp),intent(out) :: gclat,gclon ! geocentric lat,lon
!
! Local:
  real(rp) :: ylat,ylon,stm,ctm,ctc
  stm = cos(gmlat*dtr)
  ctm = sin(gmlat*dtr)
  ctc = ctp*ctm - stp*stm*cos(gmlon*dtr) ! ctp,stp are module data
  ctc = min(ctc,1.)
  ctc = max(ctc,-1.)
  gclat = asin(ctc)*rtd
  gclon = atan2(stp*stm*sin(gmlon*dtr),ctm-ctp*ctc)
!
! elon is in module data, and was set by cofrm (called from apex_mka)
!
  gclon = gclon*rtd + elon
  if (gclon < -180.) gclon = gclon + 360.
end subroutine gm2gc
!-----------------------------------------------------------------------
subroutine intrp(glat,glon,alt, gplat,gplon,gpalt, nlat,nlon,nalt, &
                 fx,fy,fz,fv,                                      &
                 dfxdth,dfydth,dfzdth,dfvdth,                      &
                 dfxdln,dfydln,dfzdln,dfvdln,                      &
                 dfxdh ,dfydh ,dfzdh ,dfvdh, ier)
!
! Args:
!
  real(rp),intent(in)    :: glat,glon,alt
  integer,intent(in) :: nlat,nlon,nalt
  real(rp),intent(in)    :: gplat(nlat),gplon(nlon),gpalt(nalt)
  real(rp),intent(out)   ::          &
    fx,fy,fz,fv,                 &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln, &
    dfxdh ,dfydh ,dfzdh ,dfvdh
  integer,intent(out) :: ier
!
! Local:
!
  integer :: i,j,k,i0,j0,k0
  real(rp) :: glonloc,xi,dlon,yj,dlat,hti,diht,zk,fac,omfac
  real(rp) :: dfxdn,dfxde,dfxdd, &
          dfydn,dfyde,dfydd, &
          dfzdn,dfzde,dfzdd, &
          dfvdn,dfvde,dfvdd, &
          dmf,dmdfdn,dmdfde,dmdfdd
  ier = 0
  glonloc = glon
  if (glonloc < gplon(1))    glonloc = glonloc + 360.
  if (glonloc > gplon(nlon)) glonloc = glonloc - 360.
!
  i0 = 0
  do i=1,nlat-1
    if (glat >= gplat(i).and.glat <= gplat(i+1)) then
      i0 = i
      dlat = gplat(i+1)-gplat(i)
      xi = (glat - gplat(i)) / dlat
      exit
    endif
  enddo
  if (i0==0) then
    write(6,"('>>> intrp: could not bracket glat=',f9.3,' in gplat=',/,(6f9.2))") &
      glat,gplat
    ier = 1
    return
  endif
  j0 = 0
  do j=1,nlon-1
    if (glon >= gplon(j).and.glon <= gplon(j+1)) then
      j0 = j
      dlon = gplon(j+1)-gplon(j)
      yj = (glon - gplon(j)) / dlon
      exit
    endif
  enddo
  if (j0==0) then
    write(6,"('>>> intrp: could not bracket glon=',f9.3,' in gplon=',/,(6f9.2))") &
      glon,gplon
    ier = 1
    return
  endif
  k0 = 0
  do k=1,nalt-1
    if (alt >= gpalt(k).and.alt <= gpalt(k+1)) then
      k0 = k
      hti = re/(re+alt)
      diht = re/(re+gpalt(k+1)) - re/(re+gpalt(k))
      zk = (hti - re/(re+gpalt(k))) / diht
      exit
    endif
  enddo
  if (k0==0) then
    write(6,"('>>> intrp: could not bracket alt=',f12.3,' in gpalt=',/,(6f12.2))") &
      alt,gpalt
    ier = 1
    return
  endif
  call trilin(xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fx,dfxdn,dfxde,dfxdd)
  dfxdth = -dfxdn*rtd/dlat
  dfxdln =  dfxde*rtd/dlon
  dfxdh  = -hti*hti*dfxdd/(re*diht)
  call trilin(yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fy,dfydn,dfyde,dfydd)
  dfydth = -dfydn*rtd/dlat
  dfydln =  dfyde*rtd/dlon
  dfydh  = -hti*hti*dfydd/(re*diht)
  call trilin(zarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fz,dfzdn,dfzde,dfzdd)
  dfzdth = -dfzdn*rtd/dlat
  dfzdln =  dfzde*rtd/dlon
  dfzdh  = -hti*hti*dfzdd/(re*diht)
  call trilin(varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,fv,dfvdn,dfvde,dfvdd)
  dfvdth = -dfvdn*rtd/dlat
  dfvdln =  dfvde*rtd/dlon
  dfvdh  = -hti*hti*dfvdd/(re*diht)
  if (nlat < 3) return
!
! Improve calculation of longitudinal derivatives near poles
!
  if (glat < dlat-90.) then
    fac = .5*xi
    omfac = 1. - fac
    xi = xi - 1.
    i0 = i0 + 1
    call trilin (xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfxdln = dfxdln*omfac + fac*dmdfde*rtd/dlon
    call trilin (yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfydln = dfydln*omfac + fac*dmdfde*rtd/dlon
    call trilin (varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfvdln = dfvdln*omfac + fac*dmdfde*rtd/dlon
  endif
  if (glat > 90.-dlat) then
    fac = .5*(1.-xi)
    omfac = 1. - fac
    xi = xi + 1.
    i0 = i0 - 1
    call trilin (xarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfxdln = dfxdln*omfac + fac*dmdfde*rtd/dlon
    call trilin (yarray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfydln = dfydln*omfac + fac*dmdfde*rtd/dlon
    call trilin (varray(i0:i0+1,j0:j0+1,k0:k0+1),nlat,nlon,xi,yj,zk,dmf,dmdfdn,dmdfde,dmdfdd)
    dfvdln = dfvdln*omfac + fac*dmdfde*rtd/dlon
  endif
end subroutine intrp
!-----------------------------------------------------------------------
subroutine trilin(u,nlat,nlon,xi,yj,zk,fu,dfudx,dfudy,dfudz)
!
! Args:
  integer,intent(in) :: &
    nlat,               & ! first dimension of u from calling routine
    nlon                  ! second dimension of u from calling routine
  real(rp),intent(in)    :: &
    u(1:2,1:2,1:2),     & ! u(1,1,1) is address of lower corner of interpolation box
    xi,  & ! fractional distance across box in x direction
    yj,  & ! fractional distance across box in y direction
    zk     ! fractional distance across box in z direction
  real(rp),intent(out)   :: &
    fu,                 & ! interpolated value of u
    dfudx,              & ! interpolated derivative of u with respect to i (x direction)
    dfudy,              & ! interpolated derivative of u with respect to j (y direction)
    dfudz                 ! interpolated derivative of u with respect to k (z direction)
!
! Local:
  real(rp) :: omxi,omyj,omzk
! write(6,"('Enter trilin: xi,yj,zk=',3e12.4)") xi,yj,zk
! write(6,"('Enter trilin: u(1,1,1),u(1,2,1),u(1,1,2),u(1,2,2)=',4e12.4)") &
!                          u(1,1,1),u(1,2,1),u(1,1,2),u(1,2,2)
! write(6,"('Enter trilin: u(2,1,1),u(2,2,1),u(2,1,2),u(2,2,2)=',4e12.4)") &
!                          u(2,1,1),u(2,2,1),u(2,1,2),u(2,2,2)
  omxi = 1. - xi
  omyj = 1. - yj
  omzk = 1. - zk
  fu = u(1,1,1)*omxi*omyj*omzk &
     + u(2,1,1)*xi*omyj*omzk   &
     + u(1,2,1)*omxi*yj*omzk   &
     + u(1,1,2)*omxi*omyj*zk   &
     + u(2,2,1)*xi*yj*omzk     &
     + u(2,1,2)*xi*omyj*zk     &
     + u(1,2,2)*omxi*yj*zk     &
     + u(2,2,2)*xi*yj*zk
  dfudx = (u(2,1,1)-u(1,1,1))*omyj*omzk &
        + (u(2,2,1)-u(1,2,1))*yj*omzk   &
        + (u(2,1,2)-u(1,1,2))*omyj*zk   &
        + (u(2,2,2)-u(1,2,2))*yj*zk
  dfudy = (u(1,2,1)-u(1,1,1))*omxi*omzk &
        + (u(2,2,1)-u(2,1,1))*xi*omzk   &
        + (u(1,2,2)-u(1,1,2))*omxi*zk   &
        + (u(2,2,2)-u(2,1,2))*xi*zk
  dfudz = (u(1,1,2)-u(1,1,1))*omxi*omyj &
        + (u(2,1,2)-u(2,1,1))*xi*omyj   &
        + (u(1,2,2)-u(1,2,1))*omxi*yj   &
        + (u(2,2,2)-u(2,2,1))*xi*yj
end subroutine trilin
!-----------------------------------------------------------------------
subroutine adpl(glat,glon,cth,sth,fx,fy,fz,fv, &
                dfxdth,dfydth,dfzdth,dfvdth,dfxdln,dfydln,dfzdln,dfvdln)
!
!  Add-back of pseudodipole component to x,y,z,v and their derivatives.
!
! Args:
  real(rp),intent(in)      :: glat,glon
  real(rp),intent(out)     :: cth,sth
  real(rp),intent(inout)   ::        &
     fx,fy,fz,fv,                &
    dfxdth,dfydth,dfzdth,dfvdth, &
    dfxdln,dfydln,dfzdln,dfvdln
!
! Local:
  real(rp) :: cph,sph,ctm
  cph = cos((glon-elon)*dtr)
  sph = sin((glon-elon)*dtr)
  cth = sin(glat*dtr)
  sth = cos(glat*dtr)
  ctm = ctp*cth + stp*sth*cph
  fx = fx + sth*ctp*cph - cth*stp
  fy = fy + sth*sph
  fz = fz + ctm
  fv = fv - ctm
  dfxdth = dfxdth + ctp*cth*cph + stp*sth
  dfydth = dfydth + cth*sph
  dfzdth = dfzdth - ctp*sth + stp*cth*cph
  dfvdth = dfvdth + ctp*sth - stp*cth*cph
  dfxdln = dfxdln - ctp*sth*sph
  dfydln = dfydln + sth*cph
  dfzdln = dfzdln - stp*sth*sph
  dfvdln = dfvdln + stp*sth*sph
end subroutine adpl
!-----------------------------------------------------------------------
subroutine setmiss(xmiss,xlatm,alon,vmp,b,bmag,be3,sim,si,F,D,W, &
  bhat,d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3)
!
! Args:
  real(rp),intent(in)  :: xmiss
  real(rp),intent(out) :: xlatm,alon,vmp,bmag,be3,sim,si,F,D,W
  real(rp),dimension(3),intent(out) :: bhat,d1,d2,d3,e1,e2,e3,b,f1,f2,f3,g1,g2,g3
  xlatm = xmiss
  alon  = xmiss
  vmp   = xmiss
  bmag  = xmiss
  be3   = xmiss
  sim   = xmiss
  si    = xmiss
  F     = xmiss
  D     = xmiss
  W     = xmiss
  bhat  = xmiss
  d1    = xmiss
  d2    = xmiss
  d3    = xmiss
  e1    = xmiss
  e2    = xmiss
  e3    = xmiss
  b     = xmiss
  f1    = xmiss
  f2    = xmiss
  f3    = xmiss
  g1    = xmiss
  g2    = xmiss
  g3    = xmiss
end subroutine setmiss
!-----------------------------------------------------------------------
subroutine cofrm(date)
  implicit none
!
! Input arg:
  real(rp),intent(in) :: date
!
! Local:
  integer :: m,n,i,l,ll,lm,nmx,nc,kmx,k,mm,nn
  real(rp) :: t,one,tc,r,f,f0, gpl
  integer,parameter :: n1=120, n2=195, isv=0
  integer,parameter :: &
    ncn1=19, & ! number of coefficients dimensioned n1
    ncn2=7     ! number of coefficients dimensioned n2: increase with each IGRF update
  integer,parameter :: ngh = n1*ncn1 + n2*ncn2 + 1 ! not sure why the extra +1
  real(rp),save :: g1(n1,ncn1), g2(n2,ncn2), gh(ngh)
  real(rp),parameter :: alt = 0.
  if (date < 1900. .or. date > 2030.) then
    write(6,"('>>> cofrm: date=',f8.2,' Date must be >= 1900 and <= 2030')") date
    stop 'cofrm'
  endif
  if (date > 2025.) then
    write(6,"('>>> WARNING cofrm:')")
    write(6,"(/,'   This version of IGRF is intended for use up to ')")
    write(6,"('     2025. Values for ',f9.3,' will be computed but')") date
    write(6,"('     may be of reduced accuracy.',/)")
  endif
  g1(:,1) = (/                                               & ! 1900
   -31543.,-2298., 5922., -677., 2905.,-1061.,  924., 1121., &
     1022.,-1469., -330., 1256.,    3.,  572.,  523.,  876., &
      628.,  195.,  660.,  -69., -361., -210.,  134.,  -75., &
     -184.,  328., -210.,  264.,   53.,    5.,  -33.,  -86., &
     -124.,  -16.,    3.,   63.,   61.,   -9.,  -11.,   83., &
     -217.,    2.,  -58.,  -35.,   59.,   36.,  -90.,  -69., &
       70.,  -55.,  -45.,    0.,  -13.,   34.,  -10.,  -41., &
       -1.,  -21.,   28.,   18.,  -12.,    6.,  -22.,   11., &
        8.,    8.,   -4.,  -14.,   -9.,    7.,    1.,  -13., &
        2.,    5.,   -9.,   16.,    5.,   -5.,    8.,  -18., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,   -1., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6.  /)
  g1(:,2) = (/                                               & ! 1905
   -31464.,-2298., 5909., -728., 2928.,-1086., 1041., 1065., &
     1037.,-1494., -357., 1239.,   34.,  635.,  480.,  880., &
      643.,  203.,  653.,  -77., -380., -201.,  146.,  -65., &
     -192.,  328., -193.,  259.,   56.,   -1.,  -32.,  -93., &
     -125.,  -26.,   11.,   62.,   60.,   -7.,  -11.,   86., &
     -221.,    4.,  -57.,  -32.,   57.,   32.,  -92.,  -67., &
       70.,  -54.,  -46.,    0.,  -14.,   33.,  -11.,  -41., &
        0.,  -20.,   28.,   18.,  -12.,    6.,  -22.,   11., &
        8.,    8.,   -4.,  -15.,   -9.,    7.,    1.,  -13., &
        2.,    5.,   -8.,   16.,    5.,   -5.,    8.,  -18., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6.  /)
  g1(:,3) = (/                                               & ! 1910
   -31354.,-2297., 5898., -769., 2948.,-1128., 1176., 1000., &
     1058.,-1524., -389., 1223.,   62.,  705.,  425.,  884., &
      660.,  211.,  644.,  -90., -400., -189.,  160.,  -55., &
     -201.,  327., -172.,  253.,   57.,   -9.,  -33., -102., &
     -126.,  -38.,   21.,   62.,   58.,   -5.,  -11.,   89., &
     -224.,    5.,  -54.,  -29.,   54.,   28.,  -95.,  -65., &
       71.,  -54.,  -47.,    1.,  -14.,   32.,  -12.,  -40., &
        1.,  -19.,   28.,   18.,  -13.,    6.,  -22.,   11., &
        8.,    8.,   -4.,  -15.,   -9.,    6.,    1.,  -13., &
        2.,    5.,   -8.,   16.,    5.,   -5.,    8.,  -18., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    2.,    4.,    2.,    0.,    0.,   -6.  /)
  g1(:,4) = (/                                               & ! 1915
   -31212.,-2306., 5875., -802., 2956.,-1191., 1309.,  917., &
     1084.,-1559., -421., 1212.,   84.,  778.,  360.,  887., &
      678.,  218.,  631., -109., -416., -173.,  178.,  -51., &
     -211.,  327., -148.,  245.,   58.,  -16.,  -34., -111., &
     -126.,  -51.,   32.,   61.,   57.,   -2.,  -10.,   93., &
     -228.,    8.,  -51.,  -26.,   49.,   23.,  -98.,  -62., &
       72.,  -54.,  -48.,    2.,  -14.,   31.,  -12.,  -38., &
        2.,  -18.,   28.,   19.,  -15.,    6.,  -22.,   11., &
        8.,    8.,   -4.,  -15.,   -9.,    6.,    2.,  -13., &
        3.,    5.,   -8.,   16.,    6.,   -5.,    8.,  -18., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    8.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    1.,    4.,    2.,    0.,    0.,   -6.  /)
  g1(:,5) = (/                                               & ! 1920
   -31060.,-2317., 5845., -839., 2959.,-1259., 1407.,  823., &
     1111.,-1600., -445., 1205.,  103.,  839.,  293.,  889., &
      695.,  220.,  616., -134., -424., -153.,  199.,  -57., &
     -221.,  326., -122.,  236.,   58.,  -23.,  -38., -119., &
     -125.,  -62.,   43.,   61.,   55.,    0.,  -10.,   96., &
     -233.,   11.,  -46.,  -22.,   44.,   18., -101.,  -57., &
       73.,  -54.,  -49.,    2.,  -14.,   29.,  -13.,  -37., &
        4.,  -16.,   28.,   19.,  -16.,    6.,  -22.,   11., &
        7.,    8.,   -3.,  -15.,   -9.,    6.,    2.,  -14., &
        4.,    5.,   -7.,   17.,    6.,   -5.,    8.,  -19., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    9.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6.  /)
  g1(:,6) = (/                                               & ! 1925
   -30926.,-2318., 5817., -893., 2969.,-1334., 1471.,  728., &
     1140.,-1645., -462., 1202.,  119.,  881.,  229.,  891., &
      711.,  216.,  601., -163., -426., -130.,  217.,  -70., &
     -230.,  326.,  -96.,  226.,   58.,  -28.,  -44., -125., &
     -122.,  -69.,   51.,   61.,   54.,    3.,   -9.,   99., &
     -238.,   14.,  -40.,  -18.,   39.,   13., -103.,  -52., &
       73.,  -54.,  -50.,    3.,  -14.,   27.,  -14.,  -35., &
        5.,  -14.,   29.,   19.,  -17.,    6.,  -21.,   11., &
        7.,    8.,   -3.,  -15.,   -9.,    6.,    2.,  -14., &
        4.,    5.,   -7.,   17.,    7.,   -5.,    8.,  -19., &
        8.,   10.,  -20.,    1.,   14.,  -11.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    9.,    2.,   10.,    0., &
       -2.,   -1.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6.  /)
  g1(:,7) = (/                                               & ! 1930
   -30805.,-2316., 5808., -951., 2980.,-1424., 1517.,  644., &
     1172.,-1692., -480., 1205.,  133.,  907.,  166.,  896., &
      727.,  205.,  584., -195., -422., -109.,  234.,  -90., &
     -237.,  327.,  -72.,  218.,   60.,  -32.,  -53., -131., &
     -118.,  -74.,   58.,   60.,   53.,    4.,   -9.,  102., &
     -242.,   19.,  -32.,  -16.,   32.,    8., -104.,  -46., &
       74.,  -54.,  -51.,    4.,  -15.,   25.,  -14.,  -34., &
        6.,  -12.,   29.,   18.,  -18.,    6.,  -20.,   11., &
        7.,    8.,   -3.,  -15.,   -9.,    5.,    2.,  -14., &
        5.,    5.,   -6.,   18.,    8.,   -5.,    8.,  -19., &
        8.,   10.,  -20.,    1.,   14.,  -12.,    5.,   12., &
       -3.,    1.,   -2.,   -2.,    9.,    3.,   10.,    0., &
       -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -2.,    1.,    4.,    3.,    0.,    0.,   -6.  /)
  g1(:,8) = (/                                               & ! 1935
   -30715.,-2306., 5812.,-1018., 2984.,-1520., 1550.,  586., &
     1206.,-1740., -494., 1215.,  146.,  918.,  101.,  903., &
      744.,  188.,  565., -226., -415.,  -90.,  249., -114., &
     -241.,  329.,  -51.,  211.,   64.,  -33.,  -64., -136., &
     -115.,  -76.,   64.,   59.,   53.,    4.,   -8.,  104., &
     -246.,   25.,  -25.,  -15.,   25.,    4., -106.,  -40., &
       74.,  -53.,  -52.,    4.,  -17.,   23.,  -14.,  -33., &
        7.,  -11.,   29.,   18.,  -19.,    6.,  -19.,   11., &
        7.,    8.,   -3.,  -15.,   -9.,    5.,    1.,  -15., &
        6.,    5.,   -6.,   18.,    8.,   -5.,    7.,  -19., &
        8.,   10.,  -20.,    1.,   15.,  -12.,    5.,   11., &
       -3.,    1.,   -3.,   -2.,    9.,    3.,   11.,    0., &
       -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6.  /)
  g1(:,9) = (/                                               & ! 1940
   -30654.,-2292., 5821.,-1106., 2981.,-1614., 1566.,  528., &
     1240.,-1790., -499., 1232.,  163.,  916.,   43.,  914., &
      762.,  169.,  550., -252., -405.,  -72.,  265., -141., &
     -241.,  334.,  -33.,  208.,   71.,  -33.,  -75., -141., &
     -113.,  -76.,   69.,   57.,   54.,    4.,   -7.,  105., &
     -249.,   33.,  -18.,  -15.,   18.,    0., -107.,  -33., &
       74.,  -53.,  -52.,    4.,  -18.,   20.,  -14.,  -31., &
        7.,   -9.,   29.,   17.,  -20.,    5.,  -19.,   11., &
        7.,    8.,   -3.,  -14.,  -10.,    5.,    1.,  -15., &
        6.,    5.,   -5.,   19.,    9.,   -5.,    7.,  -19., &
        8.,   10.,  -21.,    1.,   15.,  -12.,    5.,   11., &
       -3.,    1.,   -3.,   -2.,    9.,    3.,   11.,    1., &
       -2.,   -2.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    6.,   -4.,    4.,    0., &
        0.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6.  /)
  g1(:,10) = (/                                              & ! 1945
   -30594.,-2285., 5810.,-1244., 2990.,-1702., 1578.,  477., &
     1282.,-1834., -499., 1255.,  186.,  913.,  -11.,  944., &
      776.,  144.,  544., -276., -421.,  -55.,  304., -178., &
     -253.,  346.,  -12.,  194.,   95.,  -20.,  -67., -142., &
     -119.,  -82.,   82.,   59.,   57.,    6.,    6.,  100., &
     -246.,   16.,  -25.,   -9.,   21.,  -16., -104.,  -39., &
       70.,  -40.,  -45.,    0.,  -18.,    0.,    2.,  -29., &
        6.,  -10.,   28.,   15.,  -17.,   29.,  -22.,   13., &
        7.,   12.,   -8.,  -21.,   -5.,  -12.,    9.,   -7., &
        7.,    2.,  -10.,   18.,    7.,    3.,    2.,  -11., &
        5.,  -21.,  -27.,    1.,   17.,  -11.,   29.,    3., &
       -9.,   16.,    4.,   -3.,    9.,   -4.,    6.,   -3., &
        1.,   -4.,    8.,   -3.,   11.,    5.,    1.,    1., &
        2.,  -20.,   -5.,   -1.,   -1.,   -6.,    8.,    6., &
       -1.,   -4.,   -3.,   -2.,    5.,    0.,   -2.,   -2.  /)
  g1(:,11) = (/                                              & ! 1950
   -30554.,-2250., 5815.,-1341., 2998.,-1810., 1576.,  381., &
     1297.,-1889., -476., 1274.,  206.,  896.,  -46.,  954., &
      792.,  136.,  528., -278., -408.,  -37.,  303., -210., &
     -240.,  349.,    3.,  211.,  103.,  -20.,  -87., -147., &
     -122.,  -76.,   80.,   54.,   57.,   -1.,    4.,   99., &
     -247.,   33.,  -16.,  -12.,   12.,  -12., -105.,  -30., &
       65.,  -55.,  -35.,    2.,  -17.,    1.,    0.,  -40., &
       10.,   -7.,   36.,    5.,  -18.,   19.,  -16.,   22., &
       15.,    5.,   -4.,  -22.,   -1.,    0.,   11.,  -21., &
       15.,   -8.,  -13.,   17.,    5.,   -4.,   -1.,  -17., &
        3.,   -7.,  -24.,   -1.,   19.,  -25.,   12.,   10., &
        2.,    5.,    2.,   -5.,    8.,   -2.,    8.,    3., &
      -11.,    8.,   -7.,   -8.,    4.,   13.,   -1.,   -2., &
       13.,  -10.,   -4.,    2.,    4.,   -3.,   12.,    6., &
        3.,   -3.,    2.,    6.,   10.,   11.,    3.,    8.  /)
  g1(:,12) = (/                                              & ! 1955
   -30500.,-2215., 5820.,-1440., 3003.,-1898., 1581.,  291., &
     1302.,-1944., -462., 1288.,  216.,  882.,  -83.,  958., &
      796.,  133.,  510., -274., -397.,  -23.,  290., -230., &
     -229.,  360.,   15.,  230.,  110.,  -23.,  -98., -152., &
     -121.,  -69.,   78.,   47.,   57.,   -9.,    3.,   96., &
     -247.,   48.,   -8.,  -16.,    7.,  -12., -107.,  -24., &
       65.,  -56.,  -50.,    2.,  -24.,   10.,   -4.,  -32., &
        8.,  -11.,   28.,    9.,  -20.,   18.,  -18.,   11., &
        9.,   10.,   -6.,  -15.,  -14.,    5.,    6.,  -23., &
       10.,    3.,   -7.,   23.,    6.,   -4.,    9.,  -13., &
        4.,    9.,  -11.,   -4.,   12.,   -5.,    7.,    2., &
        6.,    4.,   -2.,    1.,   10.,    2.,    7.,    2., &
       -6.,    5.,    5.,   -3.,   -5.,   -4.,   -1.,    0., &
        2.,   -8.,   -3.,   -2.,    7.,   -4.,    4.,    1., &
       -2.,   -3.,    6.,    7.,   -2.,   -1.,    0.,   -3.  /)
  g1(:,13) = (/                                              & ! 1960
   -30421.,-2169., 5791.,-1555., 3002.,-1967., 1590.,  206., &
     1302.,-1992., -414., 1289.,  224.,  878., -130.,  957., &
      800.,  135.,  504., -278., -394.,    3.,  269., -255., &
     -222.,  362.,   16.,  242.,  125.,  -26., -117., -156., &
     -114.,  -63.,   81.,   46.,   58.,  -10.,    1.,   99., &
     -237.,   60.,   -1.,  -20.,   -2.,  -11., -113.,  -17., &
       67.,  -56.,  -55.,    5.,  -28.,   15.,   -6.,  -32., &
        7.,   -7.,   23.,   17.,  -18.,    8.,  -17.,   15., &
        6.,   11.,   -4.,  -14.,  -11.,    7.,    2.,  -18., &
       10.,    4.,   -5.,   23.,   10.,    1.,    8.,  -20., &
        4.,    6.,  -18.,    0.,   12.,   -9.,    2.,    1., &
        0.,    4.,   -3.,   -1.,    9.,   -2.,    8.,    3., &
        0.,   -1.,    5.,    1.,   -3.,    4.,    4.,    1., &
        0.,    0.,   -1.,    2.,    4.,   -5.,    6.,    1., &
        1.,   -1.,   -1.,    6.,    2.,    0.,    0.,   -7.  /)
  g1(:,14) = (/                                              & ! 1965
   -30334.,-2119., 5776.,-1662., 2997.,-2016., 1594.,  114., &
     1297.,-2038., -404., 1292.,  240.,  856., -165.,  957., &
      804.,  148.,  479., -269., -390.,   13.,  252., -269., &
     -219.,  358.,   19.,  254.,  128.,  -31., -126., -157., &
      -97.,  -62.,   81.,   45.,   61.,  -11.,    8.,  100., &
     -228.,   68.,    4.,  -32.,    1.,   -8., -111.,   -7., &
       75.,  -57.,  -61.,    4.,  -27.,   13.,   -2.,  -26., &
        6.,   -6.,   26.,   13.,  -23.,    1.,  -12.,   13., &
        5.,    7.,   -4.,  -12.,  -14.,    9.,    0.,  -16., &
        8.,    4.,   -1.,   24.,   11.,   -3.,    4.,  -17., &
        8.,   10.,  -22.,    2.,   15.,  -13.,    7.,   10., &
       -4.,   -1.,   -5.,   -1.,   10.,    5.,   10.,    1., &
       -4.,   -2.,    1.,   -2.,   -3.,    2.,    2.,    1., &
       -5.,    2.,   -2.,    6.,    4.,   -4.,    4.,    0., &
        0.,   -2.,    2.,    3.,    2.,    0.,    0.,   -6.  /)
  g1(:,15) = (/                                              & ! 1970
   -30220.,-2068., 5737.,-1781., 3000.,-2047., 1611.,   25., &
     1287.,-2091., -366., 1278.,  251.,  838., -196.,  952., &
      800.,  167.,  461., -266., -395.,   26.,  234., -279., &
     -216.,  359.,   26.,  262.,  139.,  -42., -139., -160., &
      -91.,  -56.,   83.,   43.,   64.,  -12.,   15.,  100., &
     -212.,   72.,    2.,  -37.,    3.,   -6., -112.,    1., &
       72.,  -57.,  -70.,    1.,  -27.,   14.,   -4.,  -22., &
        8.,   -2.,   23.,   13.,  -23.,   -2.,  -11.,   14., &
        6.,    7.,   -2.,  -15.,  -13.,    6.,   -3.,  -17., &
        5.,    6.,    0.,   21.,   11.,   -6.,    3.,  -16., &
        8.,   10.,  -21.,    2.,   16.,  -12.,    6.,   10., &
       -4.,   -1.,   -5.,    0.,   10.,    3.,   11.,    1., &
       -2.,   -1.,    1.,   -3.,   -3.,    1.,    2.,    1., &
       -5.,    3.,   -1.,    4.,    6.,   -4.,    4.,    0., &
        1.,   -1.,    0.,    3.,    3.,    1.,   -1.,   -4.  /)
  g1(:,16) = (/                                              & ! 1975
   -30100.,-2013., 5675.,-1902., 3010.,-2067., 1632.,  -68., &
     1276.,-2144., -333., 1260.,  262.,  830., -223.,  946., &
      791.,  191.,  438., -265., -405.,   39.,  216., -288., &
     -218.,  356.,   31.,  264.,  148.,  -59., -152., -159., &
      -83.,  -49.,   88.,   45.,   66.,  -13.,   28.,   99., &
     -198.,   75.,    1.,  -41.,    6.,   -4., -111.,   11., &
       71.,  -56.,  -77.,    1.,  -26.,   16.,   -5.,  -14., &
       10.,    0.,   22.,   12.,  -23.,   -5.,  -12.,   14., &
        6.,    6.,   -1.,  -16.,  -12.,    4.,   -8.,  -19., &
        4.,    6.,    0.,   18.,   10.,  -10.,    1.,  -17., &
        7.,   10.,  -21.,    2.,   16.,  -12.,    7.,   10., &
       -4.,   -1.,   -5.,   -1.,   10.,    4.,   11.,    1., &
       -3.,   -2.,    1.,   -3.,   -3.,    1.,    2.,    1., &
       -5.,    3.,   -2.,    4.,    5.,   -4.,    4.,   -1., &
        1.,   -1.,    0.,    3.,    3.,    1.,   -1.,   -5.  /)
  g1(:,17) = (/                                              & ! 1980
   -29992.,-1956., 5604.,-1997., 3027.,-2129., 1663., -200., &
     1281.,-2180., -336., 1251.,  271.,  833., -252.,  938., &
      782.,  212.,  398., -257., -419.,   53.,  199., -297., &
     -218.,  357.,   46.,  261.,  150.,  -74., -151., -162., &
      -78.,  -48.,   92.,   48.,   66.,  -15.,   42.,   93., &
     -192.,   71.,    4.,  -43.,   14.,   -2., -108.,   17., &
       72.,  -59.,  -82.,    2.,  -27.,   21.,   -5.,  -12., &
       16.,    1.,   18.,   11.,  -23.,   -2.,  -10.,   18., &
        6.,    7.,    0.,  -18.,  -11.,    4.,   -7.,  -22., &
        4.,    9.,    3.,   16.,    6.,  -13.,   -1.,  -15., &
        5.,   10.,  -21.,    1.,   16.,  -12.,    9.,    9., &
       -5.,   -3.,   -6.,   -1.,    9.,    7.,   10.,    2., &
       -6.,   -5.,    2.,   -4.,   -4.,    1.,    2.,    0., &
       -5.,    3.,   -2.,    6.,    5.,   -4.,    3.,    0., &
        1.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6.  /)
  g1(:,18) = (/                                              & ! 1985
   -29873.,-1905., 5500.,-2072., 3044.,-2197., 1687., -306., &
     1296.,-2208., -310., 1247.,  284.,  829., -297.,  936., &
      780.,  232.,  361., -249., -424.,   69.,  170., -297., &
     -214.,  355.,   47.,  253.,  150.,  -93., -154., -164., &
      -75.,  -46.,   95.,   53.,   65.,  -16.,   51.,   88., &
     -185.,   69.,    4.,  -48.,   16.,   -1., -102.,   21., &
       74.,  -62.,  -83.,    3.,  -27.,   24.,   -2.,   -6., &
       20.,    4.,   17.,   10.,  -23.,    0.,   -7.,   21., &
        6.,    8.,    0.,  -19.,  -11.,    5.,   -9.,  -23., &
        4.,   11.,    4.,   14.,    4.,  -15.,   -4.,  -11., &
        5.,   10.,  -21.,    1.,   15.,  -12.,    9.,    9., &
       -6.,   -3.,   -6.,   -1.,    9.,    7.,    9.,    1., &
       -7.,   -5.,    2.,   -4.,   -4.,    1.,    3.,    0., &
       -5.,    3.,   -2.,    6.,    5.,   -4.,    3.,    0., &
        1.,   -1.,    2.,    4.,    3.,    0.,    0.,   -6.  /)
  g1(:,19) = (/                                              & ! 1990
   -29775.,-1848., 5406.,-2131., 3059.,-2279., 1686., -373., &
     1314.,-2239., -284., 1248.,  293.,  802., -352.,  939., &
      780.,  247.,  325., -240., -423.,   84.,  141., -299., &
     -214.,  353.,   46.,  245.,  154., -109., -153., -165., &
      -69.,  -36.,   97.,   61.,   65.,  -16.,   59.,   82., &
     -178.,   69.,    3.,  -52.,   18.,    1.,  -96.,   24., &
       77.,  -64.,  -80.,    2.,  -26.,   26.,    0.,   -1., &
       21.,    5.,   17.,    9.,  -23.,    0.,   -4.,   23., &
        5.,   10.,   -1.,  -19.,  -10.,    6.,  -12.,  -22., &
        3.,   12.,    4.,   12.,    2.,  -16.,   -6.,  -10., &
        4.,    9.,  -20.,    1.,   15.,  -12.,   11.,    9., &
       -7.,   -4.,   -7.,   -2.,    9.,    7.,    8.,    1., &
       -7.,   -6.,    2.,   -3.,   -4.,    2.,    2.,    1., &
       -5.,    3.,   -2.,    6.,    4.,   -4.,    3.,    0., &
        1.,   -2.,    3.,    3.,    3.,   -1.,    0.,   -6.  /)
  g2(:,1) = (/                                               & ! 1995
   -29692.,-1784., 5306.,-2200., 3070.,-2366., 1681., -413., &
     1335.,-2267., -262., 1249.,  302.,  759., -427.,  940., &
      780.,  262.,  290., -236., -418.,   97.,  122., -306., &
     -214.,  352.,   46.,  235.,  165., -118., -143., -166., &
      -55.,  -17.,  107.,   68.,   67.,  -17.,   68.,   72., &
     -170.,   67.,   -1.,  -58.,   19.,    1.,  -93.,   36., &
       77.,  -72.,  -69.,    1.,  -25.,   28.,    4.,    5., &
       24.,    4.,   17.,    8.,  -24.,   -2.,   -6.,   25., &
        6.,   11.,   -6.,  -21.,   -9.,    8.,  -14.,  -23., &
        9.,   15.,    6.,   11.,   -5.,  -16.,   -7.,   -4., &
        4.,    9.,  -20.,    3.,   15.,  -10.,   12.,    8., &
       -6.,   -8.,   -8.,   -1.,    8.,   10.,    5.,   -2., &
       -8.,   -8.,    3.,   -3.,   -6.,    1.,    2.,    0., &
       -4.,    4.,   -1.,    5.,    4.,   -5.,    2.,   -1., &
        2.,   -2.,    5.,    1.,    1.,   -2.,    0.,   -7., &
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  &
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  &
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  &
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. /)
  g2(:,2) = (/                                               & ! 2000
   -29619.4,-1728.2, 5186.1,-2267.7, 3068.4,-2481.6, 1670.9, &
     -458.0, 1339.6,-2288.0, -227.6, 1252.1,  293.4,  714.5, &
     -491.1,  932.3,  786.8,  272.6,  250.0, -231.9, -403.0, &
      119.8,  111.3, -303.8, -218.8,  351.4,   43.8,  222.3, &
      171.9, -130.4, -133.1, -168.6,  -39.3,  -12.9,  106.3, &
       72.3,   68.2,  -17.4,   74.2,   63.7, -160.9,   65.1, &
       -5.9,  -61.2,   16.9,    0.7,  -90.4,   43.8,   79.0, &
      -74.0,  -64.6,    0.0,  -24.2,   33.3,    6.2,    9.1, &
       24.0,    6.9,   14.8,    7.3,  -25.4,   -1.2,   -5.8, &
       24.4,    6.6,   11.9,   -9.2,  -21.5,   -7.9,    8.5, &
      -16.6,  -21.5,    9.1,   15.5,    7.0,    8.9,   -7.9, &
      -14.9,   -7.0,   -2.1,    5.0,    9.4,  -19.7,    3.0, &
       13.4,   -8.4,   12.5,    6.3,   -6.2,   -8.9,   -8.4, &
       -1.5,    8.4,    9.3,    3.8,   -4.3,   -8.2,   -8.2, &
        4.8,   -2.6,   -6.0,    1.7,    1.7,    0.0,   -3.1, &
        4.0,   -0.5,    4.9,    3.7,   -5.9,    1.0,   -1.2, &
        2.0,   -2.9,    4.2,    0.2,    0.3,   -2.2,   -1.1, &
       -7.4,    2.7,   -1.7,    0.1,   -1.9,    1.3,    1.5, &
       -0.9,   -0.1,   -2.6,    0.1,    0.9,   -0.7,   -0.7, &
        0.7,   -2.8,    1.7,   -0.9,    0.1,   -1.2,    1.2, &
       -1.9,    4.0,   -0.9,   -2.2,   -0.3,   -0.4,    0.2, &
        0.3,    0.9,    2.5,   -0.2,   -2.6,    0.9,    0.7, &
       -0.5,    0.3,    0.3,    0.0,   -0.3,    0.0,   -0.4, &
        0.3,   -0.1,   -0.9,   -0.2,   -0.4,   -0.4,    0.8, &
       -0.2,   -0.9,   -0.9,    0.3,    0.2,    0.1,    1.8, &
       -0.4,   -0.4,    1.3,   -1.0,   -0.4,   -0.1,    0.7, &
        0.7,   -0.4,    0.3,    0.3,    0.6,   -0.1,    0.3, &
        0.4,   -0.2,    0.0,   -0.5,    0.1,   -0.9 /)
  g2(:,3) = (/                                               & ! 2005
   -29554.63,-1669.05, 5077.99,-2337.24, 3047.69,-2594.50,   &
    1657.76, -515.43, 1336.30,-2305.83, -198.86, 1246.39,    &
     269.72,  672.51, -524.72,  920.55,  797.96,  282.07,    &
     210.65, -225.23, -379.86,  145.15,  100.00, -305.36,    &
    -227.00,  354.41,   42.72,  208.95,  180.25, -136.54,    &
    -123.45, -168.05,  -19.57,  -13.55,  103.85,   73.60,    &
      69.56,  -20.33,   76.74,   54.75, -151.34,   63.63,    &
     -14.58,  -63.53,   14.58,    0.24,  -86.36,   50.94,    &
      79.88,  -74.46,  -61.14,   -1.65,  -22.57,   38.73,    &
       6.82,   12.30,   25.35,    9.37,   10.93,    5.42,    &
     -26.32,    1.94,   -4.64,   24.80,    7.62,   11.20,    &
     -11.73,  -20.88,   -6.88,    9.83,  -18.11,  -19.71,    &
      10.17,   16.22,    9.36,    7.61,  -11.25,  -12.76,    &
      -4.87,   -0.06,    5.58,    9.76,  -20.11,    3.58,    &
      12.69,   -6.94,   12.67,    5.01,   -6.72,  -10.76,    &
      -8.16,   -1.25,    8.10,    8.76,    2.92,   -6.66,    &
      -7.73,   -9.22,    6.01,   -2.17,   -6.12,    2.19,    &
       1.42,    0.10,   -2.35,    4.46,   -0.15,    4.76,    &
       3.06,   -6.58,    0.29,   -1.01,    2.06,   -3.47,    &
       3.77,   -0.86,   -0.21,   -2.31,   -2.09,   -7.93,    &
       2.95,   -1.60,    0.26,   -1.88,    1.44,    1.44,    &
      -0.77,   -0.31,   -2.27,    0.29,    0.90,   -0.79,    &
      -0.58,    0.53,   -2.69,    1.80,   -1.08,    0.16,    &
      -1.58,    0.96,   -1.90,    3.99,   -1.39,   -2.15,    &
      -0.29,   -0.55,    0.21,    0.23,    0.89,    2.38,    &
      -0.38,   -2.63,    0.96,    0.61,   -0.30,    0.40,    &
       0.46,    0.01,   -0.35,    0.02,   -0.36,    0.28,    &
       0.08,   -0.87,   -0.49,   -0.34,   -0.08,    0.88,    &
      -0.16,   -0.88,   -0.76,    0.30,    0.33,    0.28,    &
       1.72,   -0.43,   -0.54,    1.18,   -1.07,   -0.37,    &
      -0.04,    0.75,    0.63,   -0.26,    0.21,    0.35,    &
       0.53,   -0.05,    0.38,    0.41,   -0.22,   -0.10,    &
      -0.57,   -0.18,   -0.82 /)
  g2(:,4) = (/                                                & ! 2010
     -29496.57,-1586.42, 4944.26,-2396.06, 3026.34,-2708.54,  &
       1668.17, -575.73, 1339.85,-2326.54, -160.40, 1232.10,  &
    	251.75,  633.73, -537.03,  912.66,  808.97,  286.48,  &
    	166.58, -211.03, -356.83,  164.46,   89.40, -309.72,  &
       -230.87,  357.29,   44.58,  200.26,  189.01, -141.05,  &
       -118.06, -163.17,   -0.01,   -8.03,  101.04,   72.78,  &
    	 68.69,  -20.90,   75.92,   44.18, -141.40,   61.54,  &
    	-22.83,  -66.26,   13.10,    3.02,  -78.09,   55.40,  &
    	 80.44,  -75.00,  -57.80,   -4.55,  -21.20,   45.24,  &
    	  6.54,   14.00,   24.96,   10.46,    7.03,    1.64,  &
    	-27.61,    4.92,   -3.28,   24.41,    8.21,   10.84,  &
    	-14.50,  -20.03,   -5.59,   11.83,  -19.34,  -17.41,  &
    	 11.61,   16.71,   10.85,    6.96,  -14.05,  -10.74,  &
    	 -3.54,    1.64,    5.50,    9.45,  -20.54,    3.45,  &
    	 11.51,   -5.27,   12.75,    3.13,   -7.14,  -12.38,  &
    	 -7.42,   -0.76,    7.97,    8.43,    2.14,   -8.42,  &
    	 -6.08,  -10.08,    7.01,   -1.94,   -6.24,    2.73,  &
    	  0.89,   -0.10,   -1.07,    4.71,   -0.16,    4.44,  &
    	  2.45,   -7.22,   -0.33,   -0.96,    2.13,   -3.95,  &
    	  3.09,   -1.99,   -1.03,   -1.97,   -2.80,   -8.31,  &
    	  3.05,   -1.48,    0.13,   -2.03,    1.67,    1.65,  &
    	 -0.66,   -0.51,   -1.76,    0.54,    0.85,   -0.79,  &
    	 -0.39,    0.37,   -2.51,    1.79,   -1.27,    0.12,  &
    	 -2.11,    0.75,   -1.94,    3.75,   -1.86,   -2.12,  &
    	 -0.21,   -0.87,    0.30,    0.27,    1.04,    2.13,  &
    	 -0.63,   -2.49,    0.95,    0.49,   -0.11,    0.59,  &
    	  0.52,    0.00,   -0.39,    0.13,   -0.37,    0.27,  &
    	  0.21,   -0.86,   -0.77,   -0.23,    0.04,    0.87,  &
    	 -0.09,   -0.89,   -0.87,    0.31,    0.30,    0.42,  &
    	  1.66,   -0.45,   -0.59,    1.08,   -1.14,   -0.31,  &
    	 -0.07,    0.78,    0.54,   -0.18,    0.10,    0.38,  &
    	  0.49,    0.02,    0.44,    0.42,   -0.25,   -0.26,  &
    	 -0.53,   -0.26,   -0.79/)

   g2(:,5) = (/                                               & ! 2015
     -29441.46,-1501.77, 4795.99,-2445.88, 3012.20,-2845.41,  &
       1676.35, -642.17, 1350.33,-2352.26, -115.29, 1225.85,  &
     	245.04,  581.69, -538.70,  907.42,  813.68,  283.54,  &
     	120.49, -188.43, -334.85,  180.95,   70.38, -329.23,  &
       -232.91,  360.14,   46.98,  192.35,  196.98, -140.94,  &
       -119.14, -157.40,   15.98,    4.30,  100.12,   69.55,  &
     	 67.57,  -20.61,   72.79,   33.30, -129.85,   58.74,  &
     	-28.93,  -66.64,   13.14,    7.35,  -70.85,   62.41,  &
     	 81.29,  -75.99,  -54.27,   -6.79,  -19.53,   51.82,  &
     	  5.59,   15.07,   24.45,    9.32,    3.27,   -2.88,  &
     	-27.50,    6.61,   -2.32,   23.98,    8.89,   10.04,  &
     	-16.78,  -18.26,   -3.16,   13.18,  -20.56,  -14.60,  &
     	 13.33,   16.16,   11.76,    5.69,  -15.98,   -9.10,  &
     	 -2.02,    2.26,    5.33,    8.83,  -21.77,    3.02,  &
     	 10.76,   -3.22,   11.74,    0.67,   -6.74,  -13.20,  &
     	 -6.88,   -0.10,    7.79,    8.68,    1.04,   -9.06,  &
     	 -3.89,  -10.54,    8.44,   -2.01,   -6.26,    3.28,  &
     	  0.17,   -0.40,    0.55,    4.55,   -0.55,    4.40,  &
     	  1.70,   -7.92,   -0.67,   -0.61,    2.13,   -4.16,  &
     	  2.33,   -2.85,   -1.80,   -1.12,   -3.59,   -8.72,  &
     	  3.00,   -1.40,    0.00,   -2.30,    2.11,    2.08,  &
     	 -0.60,   -0.79,   -1.05,    0.58,    0.76,   -0.70,  &
     	 -0.20,    0.14,   -2.12,    1.70,   -1.44,   -0.22,  &
     	 -2.57,    0.44,   -2.01,    3.49,   -2.34,   -2.09,  &
     	 -0.16,   -1.08,    0.46,    0.37,    1.23,    1.75,  &
     	 -0.89,   -2.19,    0.85,    0.27,    0.10,    0.72,  &
     	  0.54,   -0.09,   -0.37,    0.29,   -0.43,    0.23,  &
     	  0.22,   -0.89,   -0.94,   -0.16,   -0.03,    0.72,  &
     	 -0.02,   -0.92,   -0.88,    0.42,    0.49,    0.63,  &
     	  1.56,   -0.42,   -0.50,    0.96,   -1.24,   -0.19,  &
     	 -0.10,    0.81,    0.42,   -0.13,   -0.04,    0.38,  &
     	  0.48,    0.08,    0.48,    0.46,   -0.30,   -0.35,  &
     	 -0.43,   -0.36,   -0.71/)

   g2(:,6) = (/                                               & ! 2020
      -29404.8, -1450.9,  4652.5, -2499.6,  2982.0, -2991.6,  &
        1676.9,  -734.6,  1363.2, -2381.2,   -82.1,  1236.2,  &
         241.9,   525.7,  -543.4,   903.0,   809.5,   281.9,  &
          86.3,  -158.4,  -309.4,   199.7,    48.0,  -349.7,  &
        -234.3,   363.2,    47.7,   187.8,   208.3,  -140.6,  &
        -121.2,  -151.2,    32.3,    13.5,    98.8,    66.0,  &
          65.5,   -19.1,    72.9,    25.1,  -121.5,    52.8,  &
         -36.2,   -64.5,    13.5,     8.9,   -64.7,    68.1,  &
          80.6,   -76.7,   -51.5,    -8.2,   -16.9,    56.5,  &
           2.2,    15.8,    23.5,     6.4,    -2.2,    -7.2,  &
         -27.2,     9.8,    -1.8,    23.7,     9.7,	8.4,  &
         -17.5,   -15.3,    -0.5,    12.8,   -21.1,   -11.7,  &
          15.3,    14.9,    13.7,     3.6,   -16.5,    -6.9,  &
          -0.3,     2.8,     5.0,     8.4,   -23.4,	2.9,  &
          11.0,    -1.5,     9.8,    -1.1,    -5.1,   -13.2,  &
          -6.3,     1.1,     7.8,     8.8,     0.4,    -9.3,  &
          -1.4,   -11.9,     9.6,    -1.9,    -6.2,	3.4,  &
          -0.1,    -0.2,     1.7,     3.5,    -0.9,	4.8,  &
           0.7,    -8.6,    -0.9,    -0.1,     1.9,    -4.3,  &
           1.4,    -3.4,    -2.4,    -0.1,    -3.8,    -8.8,  &
           3.0,    -1.4,     0.0,    -2.5,     2.5,	2.3,  &
          -0.6,    -0.9,    -0.4,     0.3,     0.6,    -0.7,  &
          -0.2,    -0.1,    -1.7,     1.4,    -1.6,    -0.6,  &
          -3.0,     0.2,    -2.0,     3.1,    -2.5,    -2.0,  &
          -0.1,    -1.2,     0.4,     0.5,     1.3,	1.4,  &
          -1.1,    -1.8,     0.7,     0.1,     0.3,	0.8,  &
           0.5,    -0.2,    -0.3,     0.6,    -0.5,	0.2,  &
           0.1,    -0.9,    -1.1,    -0.0,    -0.3,	0.5,  &
           0.1,    -0.9,    -0.9,     0.5,     0.6,	0.7,  &
           1.4,    -0.3,    -0.4,     0.8,    -1.3,    -0.0,  &
          -0.1,     0.8,     0.3,    -0.0,    -0.1,	0.4,  &
           0.5,     0.1,     0.5,     0.5,    -0.4,    -0.5,  &
          -0.3,    -0.4,    -0.6/)
!
    g2(1:80,7) = (/                                           & ! 2022
        5.7,     7.4,   -25.9,   -11.0,    -7.0,   -30.2,     &
       -2.1,   -22.4,	  2.2,    -5.9,     6.0,     3.1,     &
       -1.1,   -12.0,	  0.5,    -1.2,    -1.5,    -0.1,     &
       -5.9,	 6.5,	  5.2,     3.5,    -5.1,    -5.0,     &
       -0.3,	 0.5,	 -0.0,    -0.6,     2.5,     0.2,     &
       -0.6,	 1.3,	  3.0,     0.9,     0.3,    -0.5,     &
       -0.3,	 0.0,	  0.4,    -1.6,     1.3,    -1.3,     &
       -1.4,	 0.8,	 -0.0,     0.0,     0.9,     1.0,     &
       -0.1,	-0.2,	  0.6,    -0.0,     0.6,     0.7,     &
       -0.8,	 0.1,	 -0.2,    -0.5,    -1.1,    -0.8,     &
     	0.1,	 0.8,	  0.3,    -0.0,     0.1,    -0.2,     &
       -0.1,	 0.6,	  0.4,    -0.2,    -0.1,     0.5,     &
     	0.4,	-0.3,	  0.3,    -0.4,    -0.1,     0.5,     &
     	0.4,	 0.0/)
    g2(81:n2,7) = 0.0
!
! Set gh from g1,g2:
!
  do n=1,ncn1
    i = (n-1)*n1
    gh(i+1:i+n1) = g1(:,n)
!   write(6,"('cofrm: n=',i3,' i+1:i+n1=',i4,':',i4)") n,i+1,i+n1
  enddo
  do n=1,ncn2
    i = n1*ncn1 + (n-1)*n2
    gh(i+1:i+n2) = g2(:,n)
!   write(6,"('cofrm: n=',i3,' i+1:i+n2=',i4,':',i4)") n,i+1,i+n2
  enddo
  gh(ngh) = 0. ! not sure why gh is dimensioned with the extra element, so set it to 0.

  if (date < 2020.) then
    t   = 0.2*(date - 1900.)
    ll  = t
    one = ll
    t   = t - one
    if (date < 1995.) then
      nmx   = 10
      nc    = nmx*(nmx+2)
      ll    = nc*ll
      kmx   = (nmx+1)*(nmx+2)/2
    else
      nmx   = 13
      nc    = nmx*(nmx+2)
      ll    = 0.2*(date - 1995.0)
      ll    = 120*19 + nc*ll
      kmx   = (nmx+1)*(nmx+2)/2
    endif
    tc    = 1.0 - t
    if (isv.eq.1) then
      tc = -0.2
      t = 0.2
    endif
  else ! date >= 2020
    t     = date - 2020.0
    tc    = 1.0
    if (isv.eq.1) then
      t = 1.0
      tc = 0.0
    end if
    ll    = 3255
    nmx   = 13
    nc    = nmx*(nmx+2)
    kmx   = (nmx+1)*(nmx+2)/2
  endif ! date < 2020
  r = alt
  l = 1
  m = 1
  n = 0
!
! Set outputs gb(ncoef) and gv(ncoef)
! These are module data above.
!
  gb(1) = 0.
  gv(1) = 0.
  f0 = -1.e-5
  do k=2,kmx
    if (n < m) then
      m = 0
      n = n+1
    endif ! n < m
    lm = ll + l
    if (m == 0) f0 = f0 * float(n)/2.
    if (m == 0) f  = f0 / sqrt(2.0)
    nn = n+1
    mm = 1

    if (m /= 0) then
      f = f / sqrt(float(n-m+1) / float(n+m) )
      gb(l+1)  = (tc*gh(lm) + t*gh(lm+nc))* f
    else
      gb(l+1)  = (tc*gh(lm) + t*gh(lm+nc))* f0
    endif
    gv(l+1) = gb(l+1)/float(nn)
    if (m /= 0) then
      gb(l+2)  = (tc*gh(lm+1) + t*gh(lm+nc+1))*f
      gv(l+2) = gb(l+2)/float(nn)
      l = l+2
    else
      l = l+1
    endif
    m = m+1
  enddo
! write(6,"('cofrm: ncoef=',i4,' gb=',/,(6f12.3))") ncoef,gb
! write(6,"('cofrm: ncoef=',i4,' gv=',/,(6f12.3))") ncoef,gv
!
! Compute geocentric colatitude (colat) and longitude (elon) of the
!  north pole of earth centered dipole, and cosine (ctp) and sine (stp)
!  of colat, for use in other routines.
!
  gpl = sqrt( gb(2  )**2+ gb(3  )**2+ gb(4  )**2)
  ctp = gb(2  )/gpl
  stp = sqrt(1. - ctp*ctp)
  colat = (acos(ctp))*rtd
  elon = atan2( gb(4  ), gb(3  ))*rtd
!
! Compute magnitude of magnetic potential at pole, radius Re.
!      .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through f0).
!
  vp = .2*gpl*re
end subroutine cofrm
!-----------------------------------------------------------------------
subroutine subsol(iyr,iday,ihr,imn,sec,sbsllat,sbsllon)
!
! Find subsolar geographic latitude and longitude given the
! date and time (Universal Time).
!
! This is based on formulas in Astronomical Almanac for the
! year 1996, p.  C24. (U.S.  Government Printing Office,
! 1994).  According to the Almanac, results are good to at
! least 0.01 degree latitude and 0.025 degree longitude
! between years 1950 and 2050.  Accuracy for other years has
! not been tested although the algorithm has been designed to
! accept input dates from 1601 to 2100.  Every day is assumed
! to have exactly 86400 seconds; thus leap seconds that
! sometimes occur on June 30 and December 31 are ignored:
! their effect is below the accuracy threshold of the algorithm.
!
! 961026 A. D. Richmond, NCAR
!
! Input Args:
  integer,intent(in) :: &
    iyr,   & ! Year (e.g., 1994). IYR must be in the range: 1601 to 2100.
    iday,  & ! Day number of year (e.g., IDAY = 32 for Feb 1)
    ihr,   & ! Hour of day    (e.g., 13 for 13:49)
    imn      ! Minute of hour (e.g., 49 for 13:49)
  real(rp),intent(in) :: sec ! Second and fraction after the hour/minute.
!
! Output Args:
  real(rp),intent(out) :: &
    sbsllat, & ! geographic latitude of subsolar point (degrees)
    sbsllon    ! geographic longitude of subsolar point (-180 to +180)
!
! Local:
  integer,parameter :: minyear=1601, maxyear = 2100

  real(rp) :: yr,l0,g0,ut,df,lf,gf,l,g,grad,n,epsilon,epsrad,alpha,delta,&
    etdeg,aptime,lambda,lamrad,sinlam
  integer :: nleap,ncent,nrot
  sbsllat=0. ; sbsllon=0.
  yr = iyr-2000
!
! nleap (final) = number of leap days from (2000 January 1) to (IYR January 1)
!                 (negative if iyr is before 1997)
  nleap = (iyr-1601)/4
  nleap = nleap - 99
  if (iyr <= 1900) then
    if (iyr < minyear) then
      write(6,*) 'subsolr invalid before ',minyear,': input year = ',iyr
      stop 'subsolr'
    endif
    ncent = (iyr-minyear)/100
    ncent = 3 - ncent
    nleap = nleap + ncent
  endif
  if (iyr > maxyear) then
    write(6,*) 'subsolr invalid after ',maxyear,':  input year = ',iyr
    stop 'subsolr'
  endif
!
! L0 = Mean longitude of Sun at 12 UT on January 1 of IYR:
!     L0 = 280.461 + .9856474*(365*(YR-NLEAP) + 366*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = 280.461 + .9856474*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = (280.461 - 360.) + (.9856474*365 - 360.)*(YR-4*NLEAP)
!          + (.9856474*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
!
      l0 = -79.549 + (-.238699*(yr-4*nleap) + 3.08514e-2*nleap)
!
! G0 = Mean anomaly at 12 UT on January 1 of IYR:
!     G0 = 357.528 + .9856003*(365*(YR-NLEAP) + 366*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = 357.528 + .9856003*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = (357.528 - 360.) + (.9856003*365 - 360.)*(YR-4*NLEAP)
!          + (.9856003*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
!
      g0 = -2.472 + (-.2558905*(yr-4*nleap) - 3.79617e-2*nleap)
!
! Universal time in seconds:
      ut = float(ihr*3600 + imn*60) + sec
!
! Days (including fraction) since 12 UT on January 1 of IYR:
      df = (ut/86400. - 1.5) + iday
!
! Addition to Mean longitude of Sun since January 1 of IYR:
      lf = .9856474*df
!
! Addition to Mean anomaly since January 1 of IYR:
      gf = .9856003*df
!
! Mean longitude of Sun:
      l = l0 + lf
!
! Mean anomaly:
      g = g0 + gf
      grad = g*dtr
!
! Ecliptic longitude:
      lambda = l + 1.915*sin(grad) + .020*sin(2.*grad)
      lamrad = lambda*dtr
      sinlam = sin(lamrad)
!
! Days (including fraction) since 12 UT on January 1 of 2000:
      n = df + 365.*yr + float(nleap)
!
! Obliquity of ecliptic:
      epsilon = 23.439 - 4.e-7*n
      epsrad = epsilon*dtr
!
! Right ascension:
      alpha = atan2(cos(epsrad)*sinlam,cos(lamrad))*rtd
!
! Declination:
      delta = asin(sin(epsrad)*sinlam)*rtd
!
! Subsolar latitude (output argument):
      sbsllat = delta
!
! Equation of time (degrees):
      etdeg = l - alpha
      nrot = nint(etdeg/360.)
      etdeg = etdeg - float(360*nrot)
!
! Apparent time (degrees):
! Earth rotates one degree every 240 s.
      aptime = ut/240. + etdeg
!
! Subsolar longitude (output argument):
      sbsllon = 180. - aptime
      nrot = nint(sbsllon/360.)
      sbsllon = sbsllon - float(360*nrot)
end subroutine subsol
!-----------------------------------------------------------------------
subroutine solgmlon(gclat,gclon,colat,elon,mlon)
!
! Compute geomagnetic longitude of the point with geocentric spherical
!  latitude and longitude of gclat and gclon, respectively.
! 940719 A. D. Richmond, NCAR
!
! Algorithm:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude colat, east longitude elon).
!   Let gclon be longitude of point P.
!   Let TE be geocentric colatitude of point P ( = pi/2 - gclat*dtr).
!   Let ang be longitude angle from GM to P.
!   Let TP be colatitude of GM.
!   Let TF be arc length between GM and P.
!   Let PA = be geomagnetic longitude in radians, i.e., Pi minus angle measured
!     counterclockwise from arc GM-P to arc GM-GP.
!   Then, using notation c=cos, s=sin, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: stfcpa,
!     stfspa are sin(TF) times cos(PA), sin(PA), respectively.
!
! Input Args:
  real(rp),intent(in)  :: &
               gclat, & ! geocentric latitude, deg.
               gclon, & ! geocentric longitude, deg.
               colat, & ! geocentric colatitude of dipole north pole, deg.
               elon     ! geocentric east longitude of dipole north pole, deg.
!
! Output Arg: Geomagnetic dipole longitude of the point (deg, -180. to 180.)
  real(rp),intent(out) :: mlon
!
! Local:
  real(rp) :: ctp,stp,ang,cang,sang,cte,ste,stfcpa,stfspa
  ctp = cos(colat*dtr)
  stp = sqrt(1. - ctp*ctp)
  ang = (gclon-elon)*dtr
  cang = cos(ang)
  sang = sin(ang)
  cte = sin(gclat*dtr)
  ste = sqrt(1.-cte*cte)
  stfcpa = ste*ctp*cang - cte*stp
  stfspa = sang*ste
  mlon = atan2(stfspa,stfcpa)*rtd
end subroutine solgmlon
!-----------------------------------------------------------------------
subroutine mlt2alon(xmlt,sbsllat,sbsllon,clatp,polon,alonx)
!
! To go from mlt to alon
!   Inputs:
!    xmlt    = magnetic local time for the apex longitude alonx (hours,
!		 0. to 24.)
!    sbsllat = geographic latitude of subsolar point (degrees)
!    sbsllon = geographic longitude of subsolar point (degrees)
!    clatp   = Geocentric colatitude of geomagnetic dipole north pole (deg)
!    polon   = East longitude of geomagnetic dipole north pole (deg)
!
!   Output:
!    alonx   = apex magnetic longitude of the point (deg, -180. to 180.)
!
  real(rp), intent(in) :: xmlt,sbsllat,sbsllon,clatp,polon
  real(rp), intent(out):: alonx
! Local
  real(rp) :: smlon
!
! 2018 May 17 ADR: Technically, geocentric instead of geographic (geodetic)
!   latitude should be used in call to solgmlon, but the error caused is
!   likely to be small.  Correcting for this error would involve
!   modifying the algorithm to be height-dependent, since the difference
!   between geocentric and geographic latitudes is height-dependent.
!
  call solgmlon(sbsllat,sbsllon,clatp,polon,smlon)
  alonx = 15.*(xmlt-12)+smlon
  if(alonx.gt.180) alonx = alonx-360.
  if(alonx.le.-180) alonx = alonx+360.
end subroutine mlt2alon
!-----------------------------------------------------------------------
subroutine magloctm(alon,sbsllat,sbsllon,clatp,polon,mlt)
!  Computes magnetic local time from magnetic longitude, subsolar coordinates,
!   and geomagnetic pole coordinates.
!  950302 A. D. Richmond, NCAR
!  Algorithm:  mlt is calculated from the difference of the apex longitude,
!   alon, and the geomagnetic dipole longitude of the subsolar point.
!
!   Inputs:
!    alon    = apex magnetic longitude of the point (deg)
!    sbsllat = geographic latitude of subsolar point (degrees)
!    sbsllon = geographic longitude of subsolar point (degrees)
!    clatp   = Geocentric colatitude of geomagnetic dipole north pole (deg)
!    polon   = East longitude of geomagnetic dipole north pole (deg)
!
!   Output:
!    mlt (real(rp)) = magnetic local time for the apex longitude alon (hours)
!
 real(rp), intent(in) :: alon,sbsllat,sbsllon,clatp,polon
 real(rp), intent(out):: mlt
! Local
  real(rp) :: smlon

!
! 2018 May 17 ADR: Technically, geocentric instead of geographic (geodetic)
!   latitude should be used in call to solgmlon, but the error caused is
!   likely to be small.  Correcting for this error would involve
!   modifying the algorithm to be height-dependent, since the difference
!   between geocentric and geographic latitudes is height-dependent.
!
  call solgmlon(sbsllat,sbsllon,clatp,polon,smlon)
  mlt = (alon-smlon)/15.+12.
  if(mlt.ge.24.) mlt = mlt-24
  if(mlt.lt. 0.) mlt = mlt + 24.
!
end subroutine magloctm
!-----------------------------------------------------------------------
subroutine cossza(glat,glon,sbsllat,sbsllon,csza)
!  Compute cosine of solar zenith angle from geographic coordinates of
!  point in question and of subsolar point.  Note that the zenith angle
!  is the angle between the vertical vector at the point in question
!  and the vertical vector at the subsolar point, and these vectors are
!  described by geographic, not geocentric, coordinates.
!
!  INPUTS:
!    glat    = geographic latitude of point (degrees)
!    glon    = geographic longitude of point (degrees)
!    sbsllat = geographic latitude of subsolar point (degrees)
!    sbsllon = geographic longitude of subsolar point (degrees)
!  OUTPUT:
!    csza    = cosine of solar zenith angle
!
 real(rp), intent(in) :: glat,glon,sbsllat,sbsllon
 real(rp), intent(out):: csza
!Local
 real(rp) :: ct,st,cs,ss,ang
!
 ct = sin(glat*dtr)
 st = sqrt(1-ct*ct)
 cs = sin(sbsllat*dtr)
 ss = sqrt(1-cs*cs)
 ang = (glon-sbsllon)*dtr
 csza = ct*cs + st*ss*cos(ang)
!
end subroutine cossza
!-----------------------------------------------------------------------
subroutine apex_m2g(xlatm,alon,alt,hr,gdlat,gdlon)
!  Given Modified Apex coordinates determine geographic
!   coordinates (gdlat, gdlon).
  real(rp), intent(in) ::  &
             xlatm,  & ! modified magnetic apex latitude
             alon,   & ! modified magnetic apex longitude
             alt,    & ! altitude
             hr        ! reference height for modified magnetic apex coordinates
  real(rp), intent(out) :: &
             gdlat,  & ! geographic (geodetic) latitude
             gdlon     ! geographic (geodetic) longitude
! Local
  integer :: ier
  real(rp) :: qdlat, qdlon ! Quasi-Dipole latitude, longitude
!
! equation 6.2 (Richmond, 1995)
  qdlat = acos(sqrt((re+alt)/(re+hr))*cos(xlatm*dtr))*rtd
  qdlat= sign(qdlat,xlatm)
  qdlon = alon
!
! convert from quasi-dipole to geodetic coordinates
  call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ier)
  if (ier /= 0) then
    write(6,"('>>> apex_q2g error ')")
    stop 'apex_m2g apex_q2g'
  endif
!
end subroutine apex_m2g
!-----------------------------------------------------------------------
subroutine gc2gm(gclat,gclon,gmlat,gmlon)
!
! Args:
  real(rp),intent(in)  :: gclat,gclon ! geocentric lat,lon (deg)
  real(rp),intent(out) :: gmlat,gmlon ! geomagnetic dipole lat,lon (deg)
!
! Local:
  real(rp) :: stc,ctc,ctm
  stc = cos(gclat*dtr)
  ctc = sin(gclat*dtr)
  ctm = ctp*ctc + stp*stc*cos((gclon-elon)*dtr)
!
! elon is in module data, and was set by cofrm (called from apex_mka)
! ctp,stp are in module data
! dtr, rtd are in module data
!
  ctm = min(ctm,1.)
  ctm = max(ctm,-1.)
  gmlat = asin(ctm)*rtd
  gmlon = atan2(stp*stc*sin((gclon-elon)*dtr),ctp*ctm-ctc)*rtd
!
!
end subroutine gc2gm
!-----------------------------------------------------------------------
end module apex
