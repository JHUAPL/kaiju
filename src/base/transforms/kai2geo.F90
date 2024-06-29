!Various kaiju wrappers around geopack things
module kai2geo
    use kdefs
    use geopack
    use planethelper
    use math
    use apex

    implicit none

    logical, parameter, private :: doSpyro = .false. !Whether to sneak spyro

    contains
    
    !Returns SM-theta,phi of rotation axis
    subroutine GeoAxisTP(axTheta,axPhi)
        real(rp), intent(out) :: axTheta,axPhi

        real(rp), dimension(NDIM) :: xyzGEO,xyzSM

        if (doSpyro) then
            !Sneak the spyro assumption in by setting geo axis to Z-hat
            axPhi = 0.0
            axTheta = 0.0
            return
        endif

        xyzGEO = [0.0,0.0,1.0] !Coordinates of axis in geo
        call GEO2SM(xyzGEO,xyzSM)
        axPhi = atan2(xyzSM(YDIR),xyzSM(XDIR))
        axTheta = acos( xyzSM(ZDIR)/norm2(xyzSM) )

    end subroutine GeoAxisTP

    !Return geocorotation potential [kV] from ionospheric SM colat/lon [radians]
    subroutine geocorotation_from_SMTP(smclatIn,smlon,crpot)
        real(rp),intent(in) :: smclatIn,smlon
        real(rp),intent(out) :: crpot
        
        real(rp) :: r,smlat,xyz(NDIM)

        r =  (RionE*1.0e+6)/REarth !Radius of ionosphere in Re
        smlat = pi/2-smclatIn

        ! convert to XYZ
        xyz(XDIR) = r*cos(smlat)*cos(smlon)
        xyz(YDIR) = r*cos(smlat)*sin(smlon)
        xyz(ZDIR) = r*sin(smlat)

        call geocorotation_from_SMXYZ(xyz,crpot)

    end subroutine geocorotation_from_SMTP

    subroutine geocorotation_from_GEOTP(geoclat,geolon,crpot)
        real(rp), intent(in)  :: geoclat,geolon
        real(rp), intent(out) :: crpot

        real(rp) :: smclat,smlon

        call geoTP2smTP(geoclat,geolon,smclat,smlon)
        call geocorotation_from_SMTP(smclat,smlon,crpot)
        
    end subroutine geocorotation_from_GEOTP

    subroutine geocorotation_from_SMXYZ(xyz,crpot)
        real(rp), intent(in)  :: xyz(NDIM)
        real(rp), intent(out) :: crpot

        real(rp) :: r,smlon,smclat,axT,axP,A0
        !Using Eqn 14 from Hones & Bergeson 1965
        r = norm2(xyz)
        smlon  = atan2(xyz(YDIR),xyz(XDIR))
        smclat = acos (xyz(ZDIR)/norm2(xyz))
        !Get theta/phi of axis
        call GeoAxisTP(axT,axP)

        !NOTE: Negative sign gets absorbed into order of A0 x (...)
        A0 = EarthPsi0/r
        crpot = A0*( sin(smclat)*cos(smclat)*cos(smlon-axP)*sin(axT) - cos(axT)*sin(smclat)*sin(smclat) )

    end subroutine geocorotation_from_SMXYZ

    subroutine smTP2geoTP(incolat,inlon,outcolat,outlon)
        real(rp), intent(in)  :: incolat,inlon
        real(rp), intent(out) :: outcolat,outlon

        real(rp) :: r,tempt,tempp
        real(rp),dimension(NDIM) :: indim,outdim

        ! Calculate radius of ionosphere in km
        r =  (RionE*1.0e+6)/REarth

        ! temporarily convert to lat
        tempt = pi/2-incolat
      
        ! use temp values everywhere for consistency
        tempp = inlon

        ! convert to XYZ    
        indim(XDIR) = r*cos(tempt)*cos(tempp)
        indim(YDIR) = r*cos(tempt)*sin(tempp)
        indim(ZDIR) = r*sin(tempt)

        ! Do transform
        call SM2GEO(indim,outdim)
        ! Convert back to colatitude
        ! I'm going to be lazy here and not even care to convert to funky colat output
        outcolat = acos(outdim(ZDIR)/r) ! convert back to standard colatitude
        outlon = modulo((atan2(outdim(YDIR),outdim(XDIR))+2*pi),(2*pi)) 
    end subroutine smTP2geoTP

    subroutine geoTP2smTP(geocolat,geolon,smcolat,smlon)
        real(rp), intent(in)  :: geocolat,geolon
        real(rp), intent(out) :: smcolat,smlon

        real(rp) :: r,geolat
        real(rp), dimension(NDIM) :: geoXYZ,smXYZ
        r =  (RionE*1.0e+6)/REarth !Radius of ionosphere in Re

        geolat = PI/2 - geocolat
        !Convert to xyz
        geoXYZ(XDIR) = r*cos(geolat)*cos(geolon)
        geoXYZ(YDIR) = r*cos(geolat)*sin(geolon)
        geoXYZ(ZDIR) = r*sin(geolat)

        !Do transform
        call GEO2SM(geoXYZ,smXYZ)

        !Convert back to TP
        smcolat = acos(smXYZ(ZDIR)/r)
        smlon   = katan2(smXYZ(YDIR),smXYZ(XDIR))

    end subroutine geoTP2smTP

    !Take the signed invariant latitude and calculated the mirror ratio between R=1 and R=rad
    !Using dipole mapping across the gap
    function IGRFMirrorRatio(mlat,mlon,rad) result(RM)
        real(rp), intent(in) :: mlat,mlon,rad
        real(rp) :: RM

        real(rp) :: rIon
        real(rp), dimension(NDIM) :: xyzIon,xyzMir

        rIon = 1.0
        xyzIon(XDIR) = rIon*cos(mlat)*cos(mlon)
        xyzIon(YDIR) = rIon*cos(mlat)*sin(mlon)
        xyzIon(ZDIR) = rIon*sin(mlat)

        xyzMir = DipoleShift(xyzIon,rad)
        RM = norm2( IGRF_SM(xyzIon) )/norm2( IGRF_SM (xyzMir) )
    end function IGRFMirrorRatio

    subroutine transform_grid(G1,G2,transform,hemisphere,ym1,ym2)
      use gcmtypes
      use mixdefs
      type(mixGrid_T), intent(in)  :: G1
      type(mixGrid_T), intent(out) :: G2
      integer, intent(in) :: transform,hemisphere
      real(rp) :: r,tempt,tempp,ymod1,ymod2
      real(rp),dimension(NDIM) :: indim
      real(rp),dimension(:,:),allocatable :: tempz
      integer,optional :: ym1, ym2
      integer :: i,j
      
      if (present(ym1)) then; ymod1 = ym1; else; ymod1 = 1; endif
      if (present(ym2)) then; ymod2 = ym2; else; ymod2 = 1; endif
      ! Calculate radius of ionosphere in km
      r =  (RionE*1.0e+6)/REarth
      ! Copy grid info from G1 to G2
      G2%Nt = G1%Nt; G2%Np = G1%Np
      if (.not.allocated(G2%x)) allocate(G2%x(G2%Np,G2%Nt))
      if (.not.allocated(G2%y)) allocate(G2%y(G2%Np,G2%Nt))
      if (.not.allocated(G2%t)) allocate(G2%t(G2%Np,G2%Nt))
      if (.not.allocated(G2%p)) allocate(G2%p(G2%Np,G2%Nt))
      if (.not.allocated(tempz)) allocate(tempz(G2%Np,G2%Nt))
      ! interpolant
      ! True size (Np,Nt-1) (note, we include the cell between Np and 1 for periodic; 
      ! otherwise, it's unused and interpolant set to 0. below)      
      if (.not.allocated(G2%Interpolant)) allocate(G2%Interpolant(G2%Np,G2%Nt,4,4)) 
      ! Reset interpolant because we're changing grid
      G2%Interpolant = 0.0
      
      do j = 1,G2%Np
        do i=1,G2%Nt
          if (hemisphere == NORTH) then
            tempt = pi/2.-G1%t(j,i) ! temporarily convert to lat
          elseif (hemisphere == SOUTH) then
            ! Check if it's a funky mapped colat
            if (G1%t(j,i) < pi/2) then
              tempt = G1%t(j,i)-pi/2
            else
              tempt = pi/2-G1%t(j,i)
            endif
          endif
          tempp = G1%p(j,i)*ymod1       ! use temp values everywhere for consistency
          indim(XDIR) = r*cos(tempt)*cos(tempp)
          indim(YDIR) = r*cos(tempt)*sin(tempp)
          indim(ZDIR) = r*sin(tempt)
          select case (transform)
            case (iGEOtoSM)
              call GEO2SM(indim(XDIR),indim(YDIR),indim(ZDIR),G2%x(j,i),G2%y(j,i),tempz(j,i))
            case (iSMtoGEO)
              call SM2GEO(indim(XDIR),indim(YDIR),indim(ZDIR),G2%x(j,i),G2%y(j,i),tempz(j,i))
            !case (iMAGtoSM) ! This is used for WACCMX. Disabling for now.
            !  call MAG2SM(indim(XDIR),indim(YDIR),indim(ZDIR),G2%x(j,i),G2%y(j,i),tempz(j,i))
            !case (iSMtoMAG)
            !  call SM2MAG(indim(XDIR),indim(YDIR),indim(ZDIR),G2%x(j,i),G2%y(j,i),tempz(j,i))
            case (iSMtoAPEX)
              call SM2APEX(indim(XDIR),indim(YDIR),indim(ZDIR),G2%x(j,i),G2%y(j,i),tempz(j,i))
            case (iAPEXtoSM)
              call APEX2SM(indim(XDIR),indim(YDIR),indim(ZDIR),G2%x(j,i),G2%y(j,i),tempz(j,i))
            case default
              ! No transform
              G2%x(j,i) = indim(XDIR)
              G2%y(j,i) = indim(YDIR)
              tempz(j,i) = indim(ZDIR)
          end select
        end do
      end do
      
      G2%y = G2%y*ymod2
      ! Convert back to colatitude
      if (hemisphere == NORTH) then
        G2%t = acos(tempz/r) ! convert back to standard colatitude
      else if (hemisphere == SOUTH) then
        G2%t = pi-acos(tempz/r)   ! convert to funky remix colatitude
      endif
      G2%p = modulo((atan2(G2%y,G2%x)+2*pi),(2*pi)) 
      if (allocated(tempz)) deallocate(tempz)
    end subroutine transform_grid

    subroutine transform_TP(incolat,inlon,outcolat,outlon,transform,hemisphere)
      use geopack
      use mixdefs
      real(rp), intent(in) :: incolat,inlon
      real(rp), intent(out) :: outcolat,outlon
      integer, intent(in) :: transform
      integer, intent(in),optional :: hemisphere
      real(rp) :: r,tempt,tempp
      real(rp),dimension(NDIM) :: indim,outdim
      ! Calculate radius of ionosphere in km
      r =  (RionE*1.0e+6)/REarth
      
      ! temporarily convert to lat
      if (present(hemisphere) .and. hemisphere == SOUTH .and. incolat < pi/2) then
        ! Check if it's a funky mapped colat (i.e. remix's south)
        tempt = incolat-pi/2
      else
        tempt = pi/2-incolat
      end if
      ! use temp values everywhere for consistency
      tempp = inlon
      ! convert to XYZ    
      indim(XDIR) = r*cos(tempt)*cos(tempp)
      indim(YDIR) = r*cos(tempt)*sin(tempp)
      indim(ZDIR) = r*sin(tempt)
      ! Do transform
      select case (transform)
        case (iGEOtoSM)
          call GEO2SM(indim(XDIR),indim(YDIR),indim(ZDIR),outdim(XDIR),outdim(YDIR),outdim(ZDIR))
        case (iSMtoGEO)
          call SM2GEO(indim(XDIR),indim(YDIR),indim(ZDIR),outdim(XDIR),outdim(YDIR),outdim(ZDIR))
        case default
          write(*,*) "This transform has not yet been implemented: ",transform
      end select
      
      ! Convert back to colatitude
      ! I'm going to be lazy here and not even care to convert to funky colat output
      outcolat = acos(outdim(ZDIR)/r) ! convert back to standard colatitude
      outlon = modulo((atan2(outdim(YDIR),outdim(XDIR))+2*pi),(2*pi)) 
    end subroutine transform_TP

    subroutine SM2APEX(XSM,YSM,ZSM,XAPX,YAPX,ZAPX)
      real(rp), intent(in) :: XSM,YSM,ZSM
      real(rp), intent(out) :: XAPX,YAPX,ZAPX

      real(rp) :: XGEO,YGEO,ZGEO
      real(rp) :: ionalt, r, tempt,tempp ,qdlon ,malat ,qdlat

      ! Non-scalar arguments returned by APXMALL:
      real(rp) :: b0,si,vmp,w,d,be3,sim,f, &
          b(3),bhat(3), &
          d1(3),d2(3),d3(3), &
          e1(3),e2(3),e3(3), &
          f1(3),f2(3),f3(3), &
          g1(3),g2(3),g3(3)
      integer :: ist

      ! Calculate radius of ionosphere in km
      r =  (RionE*1.0e+6)/REarth
      ! Calculate altitude of ionosphere in km
      ionalt = (RIonE*1.e+3 - REarth*1.e-3)

      call SM2GEO(XSM,YSM,ZSM,XGEO,YGEO,ZGEO)

      tempt = asin(ZGEO/r)*rad2deg
      tempp = atan2(YGEO,XGEO) *rad2deg

      call apex_mall (tempt,tempp,ionalt,ionalt, &
             b,bhat,b0,si,qdlon, &
             malat,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, &
             qdlat,f,f1,f2,f3,g1,g2,g3,ist)
      if (ist /= 0) stop 'apex_mall error'

      tempt = malat *deg2rad
      tempp = qdlon *deg2rad

      XAPX = r*cos(tempt)*cos(tempp)
      YAPX = r*cos(tempt)*sin(tempp)
      ZAPX = r*sin(tempt)

    end subroutine SM2APEX

    subroutine APEX2SM(XAPX,YAPX,ZAPX,XSM,YSM,ZSM)
      real(rp), intent(in) :: XAPX,YAPX,ZAPX
      real(rp), intent(out) :: XSM,YSM,ZSM

      real(rp) :: r, ionalt, tempt, tempp
      real(rp) :: XGEO,YGEO,ZGEO,geolat,geolon

      integer :: ist

      ! Calculate radius of ionosphere in km
      r =  (RionE*1.0e+6)/REarth
      ! Calculate altitude of ionosphere in km
      ionalt = (RIonE*1.e+3 - REarth*1.e-3)

      tempt = asin(ZAPX/r)*rad2deg
      tempp = atan2(YAPX,XAPX) *rad2deg

      call apex_q2g(tempt,tempp,ionalt,geolat,geolon,ist)
      if (ist /= 0) stop 'apex_q2g error'

      geolon = geolon*deg2rad
      geolat = geolat*deg2rad
      XGEO = r*cos(geolat)*cos(geolon)
      YGEO = r*cos(geolat)*sin(geolon)
      ZGEO = r*sin(geolat)

      call GEO2SM(XGEO,YGEO,ZGEO,XSM,YSM,ZSM)

    end subroutine APEX2SM

end module kai2geo
