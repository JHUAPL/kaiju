!Various kaiju wrappers around geopack things
module kai2geo
    use kdefs
    use geopack
    use planethelper
    use math
    use apex

    implicit none

    logical, parameter, private :: doSpyro = .false. !Whether to sneak spyro

    enum, bind(C)
      enumerator :: iSMtoGEO=1,iGEOtoSM,iSMtoGSM,iGSMtoSM,iGEOtoGSM,iGSMtoGEO,iAPEXtoSM,iSMtoAPEX,iAPEXtoGEO,iGEOtoAPEX
    end enum

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

    subroutine transform_TP(incolat,inlon,outcolat,outlon,transform,hemisphere)
      use geopack
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
        case (0)
          ! No transform
          outdim(XDIR) = indim(XDIR)
          outdim(YDIR) = indim(YDIR)
          outdim(ZDIR) = indim(ZDIR)
        case (iSMtoAPEX)
          call SM2APEX(indim(XDIR),indim(YDIR),indim(ZDIR),outdim(XDIR),outdim(YDIR),outdim(ZDIR))
        case (iAPEXtoSM)
          call APEX2SM(indim(XDIR),indim(YDIR),indim(ZDIR),outdim(XDIR),outdim(YDIR),outdim(ZDIR))
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

    subroutine GEO2APEX(XGEO,YGEO,ZGEO,XAPX,YAPX,ZAPX)
      real(rp), intent(in) :: XGEO,YGEO,ZGEO
      real(rp), intent(out) :: XAPX,YAPX,ZAPX

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

    end subroutine GEO2APEX

    subroutine APEX2GEO(XAPX,YAPX,ZAPX,XGEO,YGEO,ZGEO)
      real(rp), intent(in) :: XAPX,YAPX,ZAPX
      real(rp), intent(out) :: XGEO,YGEO,ZGEO

      real(rp) :: r, ionalt, tempt, tempp
      real(rp) :: geolat,geolon

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

    end subroutine APEX2GEO

    ! Given a coordinate string original and destination
    ! This will return the appropriate mapping variable to use in the transform routines
    subroutine get_coord_mapping(sgCoord, sgDest, imap)
        character(len=*), intent(in) :: sgCoord, sgDest

        integer, intent(out) :: imap

        ! Check if coordinate systems match
        if (trim(sgCoord) == trim(sgDest)) then
            imap = 0  ! No mapping needed (placeholder value)
        else
            select case (trim(sgCoord))
            case ('GEO')
                select case (trim(sgDest))
                case ('SM')
                    imap = iSMtoGEO
                case ('APEX')
                    imap = iAPEXtoGEO
                case default
                    imap = 999
                    write(*,*) "Unsupported coordinate system Dest: ", trim(sgDest)
                end select
            case ('SM')
                select case (trim(sgDest))
                case ('GEO')
                    imap = iGEOtoSM
                case ('APEX')
                    imap = iAPEXtoSM
                case default
                    imap = 999
                    write(*,*) "Unsupported coordinate system Dest: ", trim(sgDest)
                end select
            case ('APEX')
                select case (trim(sgDest))
                case ('GEO')
                    imap = iGEOtoAPEX
                case ('SM')
                    imap = iSMtoAPEX
                case default
                    imap = 999
                    write(*,*) "Unsupported coordinate system Dest: ", trim(sgDest)
                end select
            case default
                imap = 999
                write(*,*) "Unsupported coordinate system Source: ", trim(sgCoord)
            end select
        end if

      end subroutine get_coord_mapping

end module kai2geo
