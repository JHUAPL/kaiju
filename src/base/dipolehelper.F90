!Various routines to help w/ dipole coordinate geometry
!Using conventions from Lotko+04
!Dipole coordinates: phi,nu,mu w/ phi = spherical geometry azimuthal, mu = along field
!Spherical coordinates: RTP (P=phi)
!Dipole coordinates: PHI,DNU,DMU=PNM


module dipolehelper
	use kdefs
	use math
	
	implicit none

	contains

! Helpful dipole and FTV functions

    !Dipole field from moment
    function MagsphereDipole(xyz,M0) result(Bd)
        real(rp), intent(in) :: xyz(NDIM), M0
        real(rp) :: Bd(NDIM)

        real(rp) :: rad
        real(rp), dimension(NDIM) :: m

        rad = norm2(xyz)
        m = [0.0_rp,0.0_rp,M0]
        Bd = 3*dot_product(m,xyz)*xyz/rad**5.0 - m/rad**3.0

    end function MagsphereDipole

    !Take point xyz0 and push along dipole to point with radius r
    function DipoleShift(xyz0,r) result(xyz)
        real(rp), intent(in) :: xyz0(NDIM), r
        real(rp), dimension(NDIM) :: xyz

        real(rp) :: L,mlat,mlon
        !Find L of this point
        L = DipoleL(xyz0)

        !Avoid bad values if L<r, push as far as possible
        L = max(L,r)

        !Use r = L*cos^2(latitude)
        mlat = abs(acos(sqrt(r/L)))
        mlon = atan2(xyz0(YDIR),xyz0(XDIR)) !No change in longitude
        if (mlon<0) mlon = mlon+2*PI

        if (xyz0(ZDIR)<0) then
            mlat = -abs(mlat)
        endif
        
        !Get cartesian coordinates
        xyz(XDIR) = r*cos(mlat)*cos(mlon)
        xyz(YDIR) = r*cos(mlat)*sin(mlon)
        xyz(ZDIR) = r*sin(mlat)

    end function DipoleShift

    !Calculate dipole L shell for point
    function DipoleL(r) result(Leq)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: Leq

        real(rp) :: z,rad,lat
        z = r(ZDIR)

        rad = norm2(r)
        lat = abs( asin(z/rad) )
        Leq = rad/( cos(lat)*cos(lat) )
    end function DipoleL

    !Calculate invariant latitude (RADIANS) for x,y,z vector (in Rx)
    function InvLatitude(r) result(invlat)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: invlat

        real(rp) :: z,rad,lat,Leq

        z = r(ZDIR)
        rad = norm2(r)

        lat = abs( asin(z/rad) )
        Leq = rad/( cos(lat)*cos(lat) )
        invlat = abs(acos(sqrt(1.0/Leq)))

    end function InvLatitude

    !Same as above, but get signed invariant latitude
    function sInvLatitude(r) result(sinvlat)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: sinvlat

        sinvlat = InvLatitude(r)
        if (r(ZDIR) < 0) sinvlat = -sinvlat
    end function sInvLatitude
    
    !Get mirror ratio at R for a given invariant latitude
    function MirrorRatio(invlat,rad) result(Rm)
        real(rp), intent(in) :: invlat,rad
        real(rp) :: Rm
        real(rp) :: mlat,mlon,M0
        real(rp), dimension(NDIM) :: xyzIon,xyzMir
        M0 = 1.0

        mlat = invlat
        mlon = 0.0

        xyzIon(XDIR) = cos(mlat)*cos(mlon)
        xyzIon(YDIR) = cos(mlat)*sin(mlon)
        xyzIon(ZDIR) = sin(mlat)

        xyzMir = DipoleShift(xyzIon,rad)
        Rm = norm2(MagsphereDipole(xyzIon,M0))/norm2(MagsphereDipole(xyzMir,M0))

    end function MirrorRatio

    function DipColat2L(colat) result(L)
        real(rp), intent(in) :: colat
        real(rp) :: L
        L = 1.0/sin(colat)**2
    end function DipColat2L

    !Calculate FTV of dipole, Rx/nT
    !M0g is optional mag moment in Gauss, otherwise use Earth
    function DipFTV_L(L,M0gO) result(V)
        real(rp), intent(in) :: L
        real(rp), intent(in), optional :: M0gO
        real(rp) :: V
        real(rp) :: M0g,colat
        if (present(M0gO)) then
            M0g = M0gO
        else
            M0g = EarthM0g
        endif
        
        !Now get colat
        colat = asin( sqrt(1.0/L) )
        V = DipFTV_colat(colat,M0g)
    end function DipFTV_L

    !Same units as above but w/ colat as input
    function DipFTV_colat(colat,M0gO) result(V)
        real(rp), intent(in) :: colat
        real(rp), intent(in), optional :: M0gO
        real(rp) :: V
        real(rp) :: M0g,M0,cSum,S8

        if (present(M0gO)) then
            M0g = M0gO
        else
            M0g = EarthM0g
        endif
        M0 = abs(M0g*G2nT) !Convert to nano-tesa
        cSum =  35.0     *cos(1.0*colat) -      7.0 *cos(3.0*colat) &
              +(7.0/5.0) *cos(5.0*colat) - (1.0/7.0)*cos(7.0*colat)

        cSum = cSum/64.0
        S8 = sin(colat)**8.0
        V = 2*cSum/S8/M0
    end function DipFTV_colat

    !Dipole volume between L1,L2 (m3)
    !See Gkioulidou 2016
    function MagDV(L1,L2)
        real(rp), intent(in) :: L1,L2
        real(rp) :: MagDV

        real(rp) :: L21,Re3,a
        L21 = L2**3.0 - L1**3.0
        Re3 = REarth**3.0 !m^3
        a = 64.0*PI/105.0
        MagDV = a*Re3*L21
    end function MagDV

!---
!Dipole coordinates: phi,nu,mu w/ phi = spherical geometry azimuthal, mu = along field

	!Convert vector Axyz at position xyz -> dipole components (phi,nu,mu)
    function xyz2dip(xyz,Axyz) result(Adip)
        real(rp), dimension(NDIM), intent(in) :: xyz,Axyz
        real(rp), dimension(NDIM) :: Adip   ! Ap,An,Am

        real(rp) :: r,theta,Ct,St,delta
        real(rp), dimension(NDIM) :: Artp

        r = norm2(xyz)
        Adip = 0.0
        if (r < TINY) return

        theta = acos(xyz(ZDIR)/r)
        Ct = cos(theta)
        St = sin(theta)
        delta = sqrt(1.0 + 3.0*Ct*Ct)

        Artp = xyz2rtp(xyz,Axyz)
        Adip(1) = Artp(3)                                                 ! Aphi
        Adip(2) =  -(St/delta)      *Artp(1) + (2.0_rp*Ct/delta)*Artp(2)  ! Anu
        Adip(3) =  (2.0_rp*Ct/delta)*Artp(1) + (St/delta)       *Artp(2)  ! Amu        
    end function xyz2dip

    !Convert dipole components (phi,nu,mu) at position xyz -> Cartesian vector Axyz
    function dip2xyz(xyz,Adip) result(Axyz)
        real(rp), dimension(NDIM), intent(in) :: xyz,Adip
        real(rp), dimension(NDIM) :: Axyz     ! Ax,Ay,Az

        real(rp) :: r,theta,Ct,St,delta
        real(rp), dimension(NDIM) :: Artp     ! Ar,At,Ap

        r = norm2(xyz)
        Axyz = 0.0
        if (r < TINY) return

        theta = acos(xyz(ZDIR)/r)
        Ct = cos(theta)
        St = sin(theta)
        delta = sqrt(1.0 + 3.0*Ct*Ct)

        Artp(1) = -(St/delta)       *Adip(2) + (2.0_rp*Ct/delta)*Adip(3)   ! Ar
        Artp(2) =  (2.0_rp*Ct/delta)*Adip(2) + (St/delta)       *Adip(3)   ! At
        Artp(3) = Adip(1)                                                  ! Ap

        Axyz = rtp2xyz(xyz,Artp)
    end function dip2xyz

    !Do the Lotko+04 mapping of a magnetic field BxyzP @ xyzP => xyzG
    !Using PNM
    function LotkoMap(BxyzP, xyzP, xyzG) result(BxyzG)
        real(rp), dimension(NDIM), intent(in) :: BxyzP, xyzP, xyzG
        real(rp), dimension(NDIM) :: BxyzG

        real(rp), dimension(NDIM) :: BdipP, BdipG
        real(rp) :: hphiP, hnuP, hmuP
        real(rp) :: hphiG, hnuG, hmuG

        !--- Convert to dipole (PHI, DNU, DMU)
        BdipP = xyz2dip(xyzP, BxyzP)

        !--- Metric factors
        call DipoleMetrics(xyzP, hphiP, hnuP, hmuP)
        call DipoleMetrics(xyzG, hphiG, hnuG, hmuG)

        !--- Lotko invariants (perpendicular components)
        if (hphiG > TINY) then
            BdipG(1) = BdipP(1) * (hphiP / hphiG)   ! PHI
        else
            BdipG(1) = 0.0_rp
        endif

        if (hnuG > TINY) then
            BdipG(2) = BdipP(2) * (hnuP / hnuG)     ! DNU
        else
            BdipG(2) = 0.0_rp
        endif

        !--- Kill parallel component (DMU)
        BdipG(3) = 0.0_rp

        !--- Back to Cartesian
        BxyzG = dip2xyz(xyzG, BdipG)

    end function LotkoMap

    subroutine DipoleMetrics(xyz, hphi, hnu, hmu)
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(out) :: hphi, hnu, hmu

        real(rp) :: r, Ct, St, delta

        r = norm2(xyz)
        if (r < TINY) then
            hphi = 0.0_rp
            hnu  = 0.0_rp
            hmu  = 0.0_rp
            return
        endif

        Ct = xyz(ZDIR)/r
        Ct = max(-1.0_rp, min(1.0_rp, Ct))
        St = sqrt(max(0.0_rp, 1.0_rp - Ct*Ct))

        delta = sqrt(1.0_rp + 3.0_rp*Ct*Ct)

        hphi = r * St

        if (St > TINY) then
            hnu = r*r / (St * delta)
        else
            hnu = 0.0_rp
        endif

        hmu = hphi * hnu

    end subroutine DipoleMetrics

end module dipolehelper
