
module planethelper
    use kdefs
    use cmidefs
    use gamtypes
    use helpertypes
    use strings, only : toUpper
    implicit none
    
    contains

    subroutine getPlanetParams(planet, xmlInp, pStrO, doLoudO)
        type(planet_T), intent(inout) :: planet
        type(XML_Input_T), intent(in) :: xmlInp
        character(len=*), intent(in), optional :: pStrO
        logical, intent(in), optional :: doLoudO

        character(len=strLen) :: pID !Planet ID string
        logical :: doCorot, doLoud
        real(rp) :: RIon
        real(rp) :: period  ! Rotation period in seconds

        if (present(pStrO)) then
            pID = pStrO
        else
            call xmlInp%Set_Val(pID,"/Kaiju/Gamera/prob/planet","Earth")
        endif

        if (present(doLoudO)) then
            doLoud = doLoudO
        else
            doLoud = .false.
        endif

        planet%name = trim(toUpper(pID))

        select case (trim(toUpper(pID)))

        case("Earth","earth","EARTH")
            planet%rp_m = REarth
            planet%ri_m = RionE*1.e6  ! Defined in kdefs in 10000km
            planet%grav = 9.807
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", EarthM0g)
            call xmlInp%Set_val(period, "/Kaiju/Gamera/prob/period", 86400.0)
            planet%psiCorot = CorotPotential(planet%rp_m, period, planet%magMoment)
            if (doLoud) write(*,*)" Calculated corotation potential [kV]: ", planet%psiCorot
            planet%doGrav = .true.
        case("Saturn","saturn","SATURN")
            planet%rp_m = RSaturnXE*REarth
            planet%ri_m = 1.01*RSaturnXE*REarth
            planet%grav = 10.43
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", SaturnM0g)
            planet%psiCorot = -137.83*92 !kV
            planet%doGrav = .true.
        case("Jupiter", "jupiter", "JUPITER")
            planet%rp_m = RJupiterXE*REarth
            planet%ri_m = 1.01*RJupiterXE*REarth
            planet%grav = 24.79
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", JupiterM0g)
            planet%psiCorot = -2.5*1702.9*92.0 !kV
            planet%doGrav = .true.
        case("Mercury","mercury","MERCURY")
            planet%rp_m = RMercuryXE*REarth
            planet%ri_m = 1.05*RMercuryXE*REarth
            planet%grav = 0.0
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", MercuryM0g)
            planet%psiCorot = 0.0
            planet%doGrav = .false.
         case("Neptune","NEPTUNE")
            planet%rp_m = RNeptuneXE*REarth
            planet%ri_m = 1.01*RNeptuneXE*REarth
            planet%grav = 11.15
            call xmlInp%Set_val(planet%magMoment, "/Kaiju/Gamera/prob/M0", NeptuneM0g)
            planet%psiCorot = 10.024*92.0
            planet%doGrav = .false.
        case("Other","other","OTHER") ! Defaults to Earth values
            call xmlInp%Set_Val(planet%rp_m,"/Kaiju/Gamera/prob/x0",REarth)          ! [m]
            call xmlInp%Set_Val(Rion,"/Kaiju/Gamera/prob/Rion",RionE*1.e6/REarth)    ! [Rp]
            planet%ri_m =RIon*planet%rp_m                             ! [m]
            call xmlInp%Set_Val(planet%grav,"/Kaiju/Gamera/prob/G0",9.807)           ! [m/s2] 
            call xmlInp%Set_Val(planet%magMoment,"/Kaiju/Gamera/prob/M0",EarthM0g)   ! [gauss]
            call xmlInp%Set_Val(planet%psiCorot,"/Kaiju/Gamera/prob/Psi0",EarthPsi0) ! [kV]
            call xmlInp%Set_Val(planet%doGrav,"/Kaiju/Gamera/prob/doGrav",.true.)
        end select

        call xmlInp%Set_Val(doCorot,"/Kaiju/Gamera/prob/doCorot",.true.)
        if (.not. doCorot) then
            !Zero out corotation potential
            planet%psiCorot = 0.0
        endif

    end subroutine

    !Use planet params to calculate Gamera's normalization values
    !Placed here so that Chimp can get them easily as well
    subroutine getGamPlanetNorms(planet, gv0, gT0, gB0, gP0, M0, GM0)
        type(planet_T), intent(in) :: planet
        real(rp), intent(out), optional :: gv0, gT0, gB0, gP0, M0, GM0 !Return whichever values the caller wants

        real(rp) :: my_gv0, my_gB0  ! For local use

        my_gv0 = defV0
        my_gB0 = defB0

        if(present(gv0)) then 
            gv0 = my_gv0
        endif
        if(present(gT0)) then
            gT0 = planet%rp_m/my_gv0 !Set time scaling
        endif
        if(present(gB0)) then
            gB0 = my_gB0
        endif
        if(present(gP0)) then 
            gP0 = defP0 
        endif
        if(present(M0)) then
            M0  = -planet%magMoment*1.0e+5/my_gB0 !Magnetic moment
        endif
        if(present(GM0)) then
            GM0 = planet%grav*planet%rp_m/(my_gv0*my_gv0)
        endif
    end subroutine

    subroutine printPlanetParams(planet)
        type(planet_T), intent(in) :: planet

        write(*,*) "-------------"
        write(*,*) "Planet params"
        write(*,*) "Name             : ",trim(planet%name)
        write(*,*) "Planet Radius [m]: ",planet%rp_m
        write(*,*) "Iono Radius   [m]: ",planet%ri_m
        write(*,*) "Gravity    [m/s2]: ",planet%grav
        write(*,*) "Use gravity      : ",planet%doGrav
        write(*,*) "Mag. Moment   [G]: ",planet%magMoment
        write(*,*) "Corot Pot.   [kV]: ",planet%psiCorot
        write(*,*) "-------------"
    end subroutine

    subroutine writePlanetParams(planet, doIOp, h5File, gStrO)
        !! Write planet parameters as attributes to a "Planet" group
        !! Include a gStrO if you want it to be nested within some other group for some reason
        type(planet_T), intent(in) :: planet
        logical, intent(in) :: doIOp
        character(len=strLen) :: h5File
        character(len=strLen), optional :: gStrO

        type(IOVAR_T), dimension(8) :: IOVars

        call clearIO(IOVars)
        ! All attributes
        call AddOutVar(IOVars, "Name", planet%name)
        call AddOutVar(IOVars, "Rad_surface"   , planet%rp_m, uStr="m")
        call AddOutVar(IOVars, "Rad_ionosphere", planet%ri_m, uStr="m")
        call AddOutVar(IOVars, "Gravity", planet%grav, uStr="m/s^2")
        call AddOutVar(IOVars, "Mag Moment", planet%magMoment, uStr="Gauss")
        call AddOutVar(IOVars, "Psi Corot", planet%psiCorot, uStr="kV")
        if (planet%doGrav) then
            call AddOutVar(IOVars, "doGrav", "True")
        else
            call AddOutVar(IOVars, "doGrav", "False")
        endif

        if(present(gStrO)) then
            call WriteVars(IOVars, doIOp, h5File, gStrO, "Planet")
        else
            call WriteVars(IOVars, doIOp, h5File, "Planet")
        endif
    end subroutine writePlanetParams


    function CorotPotential(Rp_m, period, Bmag) result(cPot)
        !! Calculates corotation potential, assuming dipole and rotational axes are aligned
        real(rp), intent(in) :: Rp_m  ! Planetary radius in meters
        real(rp), intent(in) :: period  ! Rotation period in seconds
        real(rp), intent(in) :: Bmag  ! Magnetic moment in Gauss
        real(rp) :: cPot  ! Corotation potential in kV

        cPot = Rp_m**2.0 * (2.0*PI/period) * (Bmag*G2T) * 1.0e-3 ! [kV]
    end function CorotPotential

! Helpful thermo conversion functions

    !Turn pressure [nPa] and temperature [keV] to density [#/cc]
    function PkT2Den(P,kT) result(D)
        real(rp), intent(in) :: P,kT
        real(rp) :: D

        D = 6.25*P/max(kT,TINY)

    end function PkT2Den

    !Turn density [#/cc] and temperature [keV] to pressure [nPa]
    function DkT2P(D,kT) result(P)
        real(rp), intent(in) :: D,kT
        real(rp) :: P

        P = max(kT,TINY)*D/6.25
    end function DkT2P

    !Turn density [#/cc] and pressure [nPa] to temperature [keV]
    function DP2kT(D,P) result(kT)
        real(rp), intent(in) :: D,P
        real(rp) :: kT
        kT = 6.25*P/max(D,TINY)
    end function DP2kT

    !Turn density [#/cc] and velocity [km/s] to dynamic pressure
    function PV2PDyn(D,V) result(PDyn)
        real(rp), intent(in) :: D,V
        real(rp) :: PDyn

        PDyn = (1.94e-6)*(D*V*V) !nPa
    end function PV2PDyn

    !Turn density [#/cc] and pressure [nPa] to sound speed [km/s]
    !Assuming protons/gamma=5/3
    function DP2Cs(D,P,GamO) result(Cs)
        real(rp), intent(in) :: D,P
        real(rp), intent(in), optional :: GamO
        real(rp) :: Cs
        real(rp) :: TiEV,Gam
        
        if (present(GamO)) then
            Gam = GamO
        else
            Gam = 5.0/3
        endif
        !From NRL plasma formulary,
        !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
        TiEV = (1.0e+3)*DP2kT(D,P) !Temp in eV
        Cs = 9.79*sqrt( Gam*TieV )
    end function DP2Cs
    
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

    !Derivative wrt colat of dipole FTV
    function DerivDipFTV(colat,M0gO) result(dVdcol)
        real(rp), intent(in) :: colat
        real(rp), intent(in), optional :: M0gO
        real(rp) :: dVdcol
        real(rp) :: M0g,M0,cSum,dSum,S8
        real(rp) :: wtfgnu

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

        !Do cotan w/ tan b/c of gnu not having cotan by default
        if (abs(colat)>TINY) then
            wtfgnu = 1/tan(colat) ! = cotan(colat)
        else
            dVdcol = 0.0
            return
        endif

        !Deriv of csum wrt colat
        dSum = (-35.0*sin(1.0*colat) + 21.0*sin(3.0*colat) &
                - 7.0*sin(5.0*colat) +  1.0*sin(7.0*colat) )/64.0
        dVdcol = (2.0/S8/M0)*( -8.0*wtfgnu*cSum + dSum )

    end function DerivDipFTV

end module planethelper
