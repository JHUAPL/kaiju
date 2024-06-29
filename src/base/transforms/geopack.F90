!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ !
!
!          ##########################################################################
!          #                                                                        #
!          #                          GEOPACK-2008_dp                               #
!          #                     (MAIN SET OF FORTRAN CODES)                        #
!          #               (double-precision version of 03/21/08)                   #
!          #                (IGRF coefficients updated 01/01/20)                    #
!          ##########################################################################
!
!
!  This collection of subroutines is a result of several upgrades of the original package
!  written by N. A. Tsyganenko in 1978-1979.
!
!  PREFATORY NOTE TO THE VERSION OF FEBRUARY 4, 2008:
!
!  To avoid inappropriate use of obsolete subroutines from earlier versions, a suffix 08 was
!  added to the name of each subroutine in this release.
!
!  A possibility has been added in this version to calculate vector components in the
!  "Geocentric Solar Wind" (GSW) coordinate system, which, to our knowledge, was first
!  introduced by Hones et al., Planet. Space Sci., v.34, p.889, 1986 (aka GSWM, see Appendix,
!  Tsyganenko et al., JGRA, v.103(A4), p.6827, 1998). The GSW system is analogous to the
!  standard GSM, except that its X-axis is antiparallel to the currently observed solar wind
!  flow vector, rather than aligned with the Earth-Sun line. The orientation of axes in the
!  GSW system can be uniquely defined by specifying three components (VGSEX,VGSEY,VGSEZ) of
!  the solar wind velocity, and in the case of a strictly radial anti-sunward flow (VGSEY=
!  VGSEZ=0) the GSW system becomes identical to the standard GSM, which fact was used here
!  to minimize the number of subroutines in the package. To that end, instead of the special
!  case of the GSM coordinates, this version uses a more general GSW system, and three more
!  input parameters are added in the subroutine RECALC_08, the observed values (VGSEX,VGSEY,
!  VGSEZ) of the solar wind velocity. Invoking RECALC_08 with VGSEY=VGSEZ=0 restores the
!  standard (sunward) orientation of the X axis, which allows one to easily convert vectors
!  between GSW and GSM, as well as to/from other standard and commonly used systems. For more
!  details, see the documentation file GEOPACK-2008.DOC.
!
!  Another modification allows users to have more control over the procedure of field line
!  mapping using the subroutine TRACE_08. To that end, three new input parameters were added
!  in that subroutine, allowing one to set (i) an upper limit, DSMAX, on the automatically
!  adjusted step size, (ii) a permissible step error, ERR, and (iii) maximal length, LMAX,
!  of arrays where field line point coordinates are stored. Minor changes were also made in
!  the tracing subroutine, to make it more compact and easier for understanding, and to
!  prevent the algorithm from making uncontrollable large number of multiple loops in some
!  cases with plasmoid-like field structures.
!
!  One more subroutine, named GEODGEO_08, was added to the package, allowing one to convert
!  geodetic coordinates of a point in space (altitude above the Earth''s WGS84 ellipsoid and
!  geodetic latitude) to geocentric radial distance and colatitude, and vice versa.
!
!  For a complete list of modifications made earlier in previous versions, see the
!  documentation file GEOPACK-2008.DOC.
!
!----------------------------------------------------------------------------------
!
module geopack

      use kdefs
      use dates
      implicit none
      ! COMMON /GEOPACK1/
      real(rp), private :: ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI, &
         SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33, &
         E11,E21,E31,E12,E22,E32,E13,E23,E33
      ! COMMON /GEOPACK2/
      real(rp), private, dimension(105) :: G, H, REC

      !Overloader for conversion routines
      interface SM2GEO
        module procedure SM2GEO_xyz,SM2GEO_vec
      end interface

      interface GEO2SM
        module procedure GEO2SM_xyz,GEO2SM_vec
      end interface

      interface SM2GSE
        module procedure SM2GSE_xyz,SM2GSE_vec
      end interface

      interface GSE2SM
        module procedure GSE2SM_xyz,GSE2SM_vec
      end interface

      contains

!*******************************************************************
!
      SUBROUTINE SUN_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
!
!  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
!  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
!
!-------  INPUT PARAMETERS:
!  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
!    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
!
!-------  OUTPUT PARAMETERS:
!  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
!  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
!  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
!  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
!
!  LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
!
!     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
!
      implicit none
      integer, intent(in) :: IYEAR,IDAY,IHOUR,MIN,ISEC
      real(rp), intent(out) :: GST,SLONG,SRASN,SDEC
      real(rp) :: RAD = 57.295779513D0
      real(rp) :: FDAY,DJ,T,VL,G,OBLIQ,SOB,SLP,SIND,COSD,SC
!<<<>>>>
!      IMPLICIT REAL*8 (A-H,O-Z)
!      DATA RAD/57.295779513D0/

      IF(IYEAR<1901.OR.IYEAR>2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.D0
      VL=DMOD(279.696678D0+0.9856473354D0*DJ,360.D0)
      GST=DMOD(279.690983D0+.9856473354D0*DJ+360.D0*FDAY+180.D0,360.D0)/RAD
      G=DMOD(358.475845D0+0.985600267D0*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946D0-0.004789D0*T)*DSIN(G)+0.020094D0*DSIN(2.D0*G))/RAD
      IF(SLONG>6.2831853D0) SLONG=SLONG-6.283185307D0
      IF (SLONG<0.D0) SLONG=SLONG+6.283185307D0
      OBLIQ=(23.45229D0-0.0130125D0*T)/RAD
      SOB=DSIN(OBLIQ)
      SLP=SLONG-9.924D-5

!   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION DUE TO
!   EARTH''S ORBITAL MOTION

      SIND=SOB*DSIN(SLP)
      COSD=DSQRT(1.D0-SIND**2)
      SC=SIND/COSD
      SDEC=DATAN(SC)
      SRASN=3.141592654D0-DATAN2(DCOS(OBLIQ)/SOB*SC,-DCOS(SLP)/COSD)
      RETURN
      END SUBROUTINE SUN_08

!=====================================================================================

!     Wraps RECALC_08 with {VGSEX,VGSEY,VGSEZ}={-400.0,0.0,0.0} for GSM
!     & SM transformations

      SUBROUTINE mjdRECALC(mjd)
          implicit none
          real(rp), intent(in) :: mjd
          integer :: iyr,idoy,imon,iday,ihr,imin
          real(rp) :: rsec
          real(rp) :: mjdin

          mjdin = mjd ! mjdin of mjd2ut is for inout intent.
          call mjd2ut(mjdin,iyr,idoy,imon,iday,ihr,imin,rsec)
          call RECALC(iyr,idoy,ihr,imin,nint(rsec))
      END SUBROUTINE mjdRECALC

      SUBROUTINE RECALC(IYEAR,IDAY,IHOUR,MIN,ISEC)
          implicit none
          integer,intent(in) :: IYEAR,IDAY,IHOUR,MIN,ISEC
          call RECALC_08(IYEAR, IDAY, IHOUR, MIN, ISEC, -400.D0, 0.D0, 0.D0)
      END

!=====================================================================================

      SUBROUTINE RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VGSEX,VGSEY,VGSEZ)

!  1. PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
!     SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.

!  2. PREPARES COEFFICIENTS USED IN THE CALCULATION OF THE MAIN GEOMAGNETIC FIELD
!      (IGRF MODEL)

!  THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
!  IGRF_GEO_08, IGRF_GSW_08, DIP_08, GEOMAG_08, GEOGSW_08, MAGSW_08, SMGSW_08, GSWGSE_08,
!  GEIGEO_08, TRACE_08, STEP_08, RHAND_08.

!  THERE IS NO NEED TO REPEATEDLY INVOKE RECALC_08, IF MULTIPLE CALCULATIONS ARE MADE
!    FOR THE SAME DATE/TIME AND SOLAR WIND FLOW DIRECTION.

!-----INPUT PARAMETERS:

!     IYEAR   -  YEAR NUMBER (FOUR DIGITS)
!     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
!     IHOUR -  HOUR OF DAY (00 TO 23)
!     MIN   -  MINUTE OF HOUR (00 TO 59)
!     ISE!  -  SECONDS OF MINUTE (00 TO 59)
!     VGSEX,VGSEY,VGSEZ - GSE (GEOCENTRIC SOLAR-ECLIPTIC) COMPONENTS OF THE OBSERVED
!                              SOLAR WIND FLOW VELOCITY (IN KM/S)

!  IMPORTANT: IF ONLY QUESTIONABLE INFORMATION (OR NO INFORMATION AT ALL) IS AVAILABLE
!             ON THE SOLAR WIND SPEED, OR, IF THE STANDARD GSM AND/OR SM COORDINATES ARE
!             INTENDED TO BE USED, THEN SET VGSEX=-400.0 AND VGSEY=VGSEZ=0. IN THIS CASE,
!             THE GSW COORDINATE SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM.

!             IF ONLY SCALAR SPEED V OF THE SOLAR WIND IS KNOWN, THEN SETTING
!             VGSEX=-V, VGSEY=29.78, VGSEZ=0.0 WILL TAKE INTO ACCOUNT THE ~4 degs
!             ABERRATION OF THE MAGNETOSPHERE DUE TO EARTH''S ORBITAL MOTION

!             IF ALL THREE GSE COMPONENTS OF THE SOLAR WIND VELOCITY ARE AVAILABLE,
!             PLEASE NOTE THAT IN SOME SOLAR WIND DATABASES THE ABERRATION EFFECT
!             HAS ALREADY BEEN TAKEN INTO ACCOUNT BY SUBTRACTING 29.78 KM/S FROM VYGSE;
!             IN THAT CASE, THE UNABERRATED (OBSERVED) VYGSE VALUES SHOULD BE RESTORED
!             BY ADDING BACK THE 29.78 KM/S CORRECTION. WHETHER OR NOT TO DO THAT, MUST
!             BE EITHER VERIFIED WITH THE DATA ORIGINATOR OR DETERMINED BY AVERAGING
!             VGSEY OVER A SUFFICIENTLY LONG TIME INTERVAL.

!-----OUTPUT PARAMETERS:  NONE (ALL OUTPUT QUANTITIES ARE PLACED
!                         INTO THE COMMON BLOCKS /GEOPACK1/ AND /GEOPACK2/)

!    OTHER SUBROUTINES CALLED BY THIS ONE: SUN_08

!    AUTHOR:  N.A. TSYGANENKO

!    ORIGINALLY WRITTEN:    DEC.1, 1991
!
!    MOST RECENT REVISION: JANUARY 01, 2020:

!     The table of IGRF coefficients was extended to include those for the epoch 2020 (igrf-12)
!         (for details, see https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt)
!-----------------------------------------------------------------------------------
!    EARLIER REVISIONS:

!    REVISION OF NOVEMBER 30, 2010:

!     The table of IGRF coefficients was extended to include those for the epoch 2010
!         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)

!    REVISION OF NOVEMBER 15, 2007: ADDED THE POSSIBILITY TO TAKE INTO ACCOUNT THE OBSERVED
!     DEFLECTION OF THE SOLAR WIND FLOW FROM STRICTLY RADIAL DIRECTION. TO THAT END, THREE
!     GSE COMPONENTS OF THE SOLAR WIND VELOCITY WERE ADDED TO THE INPUT PARAMETERS.

!    CORRECTION OF MAY 9, 2006:  INTERPOLATION OF THE COEFFICIENTS (BETWEEN
!     LABELS 50 AND 105) IS NOW MADE THROUGH THE LAST ELEMENT OF THE ARRAYS
!     G(105)  AND H(105) (PREVIOUSLY MADE ONLY THROUGH N=66, WHICH IN SOME
!     CASES CAUSED RUNTIME ERRORS)

!    REVISION OF MAY 3, 2005:
!     The table of IGRF coefficients was extended to include those for the epoch 2005
!       the maximal order of spherical harmonics was also increased up to 13
!         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)

!    REVISION OF APRIL 3, 2003:
!    The code now includes preparation of the model coefficients for the subroutines
!    IGRF_08 and GEOMAG_08. This eliminates the need for the SAVE statements, used
!    in the old versions, making the codes easier and more compiler-independent.

!      IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
      integer, intent(in) :: IYEAR,IDAY,IHOUR,MIN,ISEC
      real(rp), intent(in) :: VGSEX,VGSEY,VGSEZ

!      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI, &
!         SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33, &
!         E11,E21,E31,E12,E22,E32,E13,E23,E33

!  THE COMMON BLOCK /GEOPACK1/ CONTAINS ELEMENTS OF THE ROTATION MATRICES AND OTHER
!   PARAMETERS RELATED TO THE COORDINATE TRANSFORMATIONS PERFORMED BY THIS PACKAGE

!      COMMON /GEOPACK2/ G(105),H(105),REC(105)

!  THE COMMON BLOCK /GEOPACK2/ CONTAINS COEFFICIENTS OF THE IGRF FIELD MODEL, CALCULATED
!    FOR A GIVEN YEAR AND DAY FROM THEIR STANDARD EPOCH VALUES. THE ARRAY REC CONTAINS
!    COEFFICIENTS USED IN THE RECURSION RELATIONS FOR LEGENDRE ASSOCIATE POLYNOMIALS.

      real(rp), dimension(105) :: G65,H65,G70,H70,G75,H75,G80,H80,G85,H85, &
          G90,H90,G95,H95,G00,H00,G05,H05,G10,H10,G15,H15,G20,H20
      real(rp), dimension(45) :: DG15,DH15,DG20,DH20
!<<<>>>
!      DIMENSION G65(105),H65(105),G70(105),H70(105),G75(105),H75(105),
!     + G80(105),H80(105),G85(105),H85(105),G90(105),H90(105),G95(105),
!     + H95(105),G00(105),H00(105),G05(105),H05(105),G10(105),H10(105),
!     + G15(105),H15(105),DG15(45),DH15(45),
!     + G20(105),H20(105),DG20(45),DH20(45)

      integer :: IY,N,N2,MN,M,MNN
      real(rp) :: DT,F2,F1,S,P,AA,G_10,G_11,H_11,SQ,SQQ,SQR
      real(rp) :: GST,SLONG,SRASN,SDEC,S1,S2,S3,DJ,T,OBLIQ
      real(rp) :: DZ1,DZ2,DZ3,DY1,DY2,DY3,V,DX1,DX2,DX3,DIP1,DIP2,DIP3
      real(rp) :: X1,X2,X3,Y1,Y2,Y3,Y,Z1,Z2,Z3,EXMAGX,EXMAGY,EXMAGZ,EYMAGX,EYMAGY

      G65(1:66) = [0.D0,-30334.D0,-2119.D0,-1662.D0,2997.D0,1594.D0,1297.D0, & 
          -2038.D0,1292.D0,856.D0,957.D0,804.D0,479.D0,-390.D0,252.D0, &
          -219.D0,358.D0,254.D0,-31.D0,-157.D0,-62.D0,45.D0,61.D0,8.D0, &
          -228.D0,4.D0,1.D0,-111.D0,75.D0,-57.D0,4.D0,13.D0,-26.D0,-6.D0, &
          13.D0,1.D0,13.D0,5.D0,-4.D0,-14.D0,0.D0,8.D0,-1.D0,11.D0,4.D0, &
          8.D0,10.D0,2.D0,-13.D0,10.D0,-1.D0,-1.D0,5.D0,1.D0,-2.D0,-2.D0, &
          -3.D0,2.D0,-5.D0,-2.D0,4.D0,4.D0,0.D0,2.D0,2.D0,0.D0]
      G65(67:105) = 0.D0

      H65(1:66) = [0.D0,0.D0,5776.D0,0.D0,-2016.D0,114.D0,0.D0,-404.D0, &
          240.D0,-165.D0,0.D0,148.D0,-269.D0,13.D0,-269.D0,0.D0,19.D0, &
          128.D0,-126.D0,-97.D0,81.D0,0.D0,-11.D0,100.D0,68.D0,-32.D0,-8.D0, &
          -7.D0,0.D0,-61.D0,-27.D0,-2.D0,6.D0,26.D0,-23.D0,-12.D0,0.D0,7.D0, &
          -12.D0,9.D0,-16.D0,4.D0,24.D0,-3.D0,-17.D0,0.D0,-22.D0,15.D0,7.D0, &
          -4.D0,-5.D0,10.D0,10.D0,-4.D0,1.D0,0.D0,2.D0,1.D0,2.D0,6.D0,-4.D0, &
          0.D0,-2.D0,3.D0,0.D0,-6.D0]
      H65(67:105) = 0.D0

      G70(1:66) = [0.D0,-30220.D0,-2068.D0,-1781.D0,3000.D0,1611.D0,1287.D0, &
         -2091.D0,1278.D0,838.D0,952.D0,800.D0,461.D0,-395.D0,234.D0, &
          -216.D0,359.D0,262.D0,-42.D0,-160.D0,-56.D0,43.D0,64.D0,15.D0, &
          -212.D0,2.D0,3.D0,-112.D0,72.D0,-57.D0,1.D0,14.D0,-22.D0,-2.D0, &
          13.D0,-2.D0,14.D0,6.D0,-2.D0,-13.D0,-3.D0,5.D0,0.D0,11.D0,3.D0, &
          8.D0,10.D0,2.D0,-12.D0,10.D0,-1.D0,0.D0,3.D0,1.D0,-1.D0,-3.D0, &
          -3.D0,2.D0,-5.D0,-1.D0,6.D0,4.D0,1.D0,0.D0,3.D0,-1.D0]
      G70(67:105) = 0.D0

      H70(1:66) = [0.D0,0.D0,5737.D0,0.D0,-2047.D0,25.D0,0.D0,-366.D0, &
          251.D0,-196.D0,0.D0,167.D0,-266.D0,26.D0,-279.D0,0.D0,26.D0, &
          139.D0,-139.D0,-91.D0,83.D0,0.D0,-12.D0,100.D0,72.D0,-37.D0,-6.D0, &
          1.D0,0.D0,-70.D0,-27.D0,-4.D0,8.D0,23.D0,-23.D0,-11.D0,0.D0,7.D0, &
          -15.D0,6.D0,-17.D0,6.D0,21.D0,-6.D0,-16.D0,0.D0,-21.D0,16.D0,6.D0, &
          -4.D0,-5.D0,10.D0,11.D0,-2.D0,1.D0,0.D0,1.D0,1.D0,3.D0,4.D0,-4.D0, &
          0.D0,-1.D0,3.D0,1.D0,-4.D0]
      H70(67:105) = 0.D0

      G75(1:66) = [0.D0,-30100.D0,-2013.D0,-1902.D0,3010.D0,1632.D0,1276.D0, &
          -2144.D0,1260.D0,830.D0,946.D0,791.D0,438.D0,-405.D0,216.D0, &
          -218.D0,356.D0,264.D0,-59.D0,-159.D0,-49.D0,45.D0,66.D0,28.D0, &
          -198.D0,1.D0,6.D0,-111.D0,71.D0,-56.D0,1.D0,16.D0,-14.D0,0.D0, &
          12.D0,-5.D0,14.D0,6.D0,-1.D0,-12.D0,-8.D0,4.D0,0.D0,10.D0,1.D0, &
          7.D0,10.D0,2.D0,-12.D0,10.D0,-1.D0,-1.D0,4.D0,1.D0,-2.D0,-3.D0, &
          -3.D0,2.D0,-5.D0,-2.D0,5.D0,4.D0,1.D0,0.D0,3.D0,-1.D0]
     G75(67:105) = 0.D0

      H75(1:66) = [0.D0,0.D0,5675.D0,0.D0,-2067.D0,-68.D0,0.D0,-333.D0, &
          262.D0,-223.D0,0.D0,191.D0,-265.D0,39.D0,-288.D0,0.D0,31.D0, &
          148.D0,-152.D0,-83.D0,88.D0,0.D0,-13.D0,99.D0,75.D0,-41.D0,-4.D0, &
          11.D0,0.D0,-77.D0,-26.D0,-5.D0,10.D0,22.D0,-23.D0,-12.D0,0.D0, &
          6.D0,-16.D0,4.D0,-19.D0,6.D0,18.D0,-10.D0,-17.D0,0.D0,-21.D0, &
          16.D0,7.D0,-4.D0,-5.D0,10.D0,11.D0,-3.D0,1.D0,0.D0,1.D0,1.D0,3.D0, &
          4.D0,-4.D0,-1.D0,-1.D0,3.D0,1.D0,-5.D0]
      H75(67:105) = 0.D0

      G80(1:66) = [0.D0,-29992.D0,-1956.D0,-1997.D0,3027.D0,1663.D0,1281.D0, &
          -2180.D0,1251.D0,833.D0,938.D0,782.D0,398.D0,-419.D0,199.D0, &
          -218.D0,357.D0,261.D0,-74.D0,-162.D0,-48.D0,48.D0,66.D0,42.D0, &
          -192.D0,4.D0,14.D0,-108.D0,72.D0,-59.D0,2.D0,21.D0,-12.D0,1.D0, &
          11.D0,-2.D0,18.D0,6.D0,0.D0,-11.D0,-7.D0,4.D0,3.D0,6.D0,-1.D0, &
          5.D0,10.D0,1.D0,-12.D0,9.D0,-3.D0,-1.D0,7.D0,2.D0,-5.D0,-4.D0, &
          -4.D0,2.D0,-5.D0,-2.D0,5.D0,3.D0,1.D0,2.D0,3.D0,0.D0]
      G80(67:105) = 0.D0

      H80(1:66) = [ 0.D0,0.D0,5604.D0,0.D0,-2129.D0,-200.D0,0.D0,-336.D0, &
          271.D0,-252.D0,0.D0,212.D0,-257.D0,53.D0,-297.D0,0.D0,46.D0, &
          150.D0,-151.D0,-78.D0,92.D0,0.D0,-15.D0,93.D0,71.D0,-43.D0,-2.D0, &
          17.D0,0.D0,-82.D0,-27.D0,-5.D0,16.D0,18.D0,-23.D0,-10.D0,0.D0, &
          7.D0,-18.D0,4.D0,-22.D0,9.D0,16.D0,-13.D0,-15.D0,0.D0,-21.D0, &
          16.D0,9.D0,-5.D0,-6.D0,9.D0,10.D0,-6.D0,2.D0,0.D0,1.D0,0.D0,3.D0, &
          6.D0,-4.D0,0.D0,-1.D0,4.D0,0.D0,-6.D0]
      H80(67:105) = 0.D0

      G85(1:66) = [0.D0,-29873.D0,-1905.D0,-2072.D0,3044.D0,1687.D0,1296.D0, &
          -2208.D0,1247.D0,829.D0,936.D0,780.D0,361.D0,-424.D0,170.D0, &
          -214.D0,355.D0,253.D0,-93.D0,-164.D0,-46.D0,53.D0,65.D0,51.D0, &
          -185.D0,4.D0,16.D0,-102.D0,74.D0,-62.D0,3.D0,24.D0,-6.D0,4.D0, &
          10.D0,0.D0,21.D0,6.D0,0.D0,-11.D0,-9.D0,4.D0,4.D0,4.D0,-4.D0,5.D0, &
          10.D0,1.D0,-12.D0,9.D0,-3.D0,-1.D0,7.D0,1.D0,-5.D0,-4.D0,-4.D0, &
          3.D0,-5.D0,-2.D0,5.D0,3.D0,1.D0,2.D0,3.D0,0.D0]
      G85(67:105) = 0.D0

      H85(1:66) = [ 0.D0,0.D0,5500.D0,0.D0,-2197.D0,-306.D0,0.D0,-310.D0, &
          284.D0,-297.D0,0.D0,232.D0,-249.D0,69.D0,-297.D0,0.D0,47.D0, &
          150.D0,-154.D0,-75.D0,95.D0,0.D0,-16.D0,88.D0,69.D0,-48.D0,-1.D0, &
          21.D0,0.D0,-83.D0,-27.D0,-2.D0,20.D0,17.D0,-23.D0,-7.D0,0.D0,8.D0, &
          -19.D0,5.D0,-23.D0,11.D0,14.D0,-15.D0,-11.D0,0.D0,-21.D0,15.D0, &
          9.D0,-6.D0,-6.D0,9.D0,9.D0,-7.D0,2.D0,0.D0,1.D0,0.D0,3.D0,6.D0, &
          -4.D0,0.D0,-1.D0,4.D0,0.D0,-6.D0]
      H85(67:105) = 0.D0

      G90(1:66) = [0.D0,-29775.D0,-1848.D0,-2131.D0,3059.D0,1686.D0,1314.D0, &
               -2239.D0,  1248.D0,  802.D0,  939.D0, 780.D0, 325.D0,-423.D0, &
                 141.D0,  -214.D0,  353.D0,  245.D0,-109.D0,-165.D0, -36.D0, &
                  61.D0,    65.D0,   59.D0, -178.D0,   3.D0,  18.D0, -96.D0, &
                  77.D0,   -64.D0,    2.D0,   26.D0,  -1.D0,   5.D0,   9.D0, &
                   0.D0,    23.D0,    5.D0,   -1.D0, -10.D0, -12.D0,   3.D0, &
                   4.D0,     2.D0,   -6.D0,    4.D0,   9.D0,   1.D0, -12.D0, &
                   9.D0,    -4.D0,   -2.D0,    7.D0,   1.D0,  -6.D0,  -3.D0, &
                  -4.D0,     2.D0,   -5.D0,   -2.D0,   4.D0,   3.D0,   1.D0, &
                   3.D0,     3.D0,    0.D0]
      G90(67:105) = 0.D0

      H90(1:66) = [0.D0,  0.D0,5406.D0,   0.D0,-2279.D0,-373.D0,  0.D0, &
                -284.D0,293.D0,-352.D0,   0.D0,  247.D0,-240.D0, 84.D0, &
                -299.D0,  0.D0,  46.D0, 154.D0, -153.D0, -69.D0, 97.D0, &
                   0.D0,-16.D0,  82.D0,  69.D0,  -52.D0,   1.D0, 24.D0, &
                   0.D0,-80.D0, -26.D0,   0.D0,   21.D0,  17.D0,-23.D0, &
                  -4.D0,  0.D0,  10.D0, -19.D0,    6.D0, -22.D0, 12.D0, &
                  12.D0,-16.D0, -10.D0,   0.D0,  -20.D0,  15.D0, 11.D0, &
                  -7.D0, -7.D0,   9.D0,   8.D0,   -7.D0,   2.D0,  0.D0, &
                   2.D0,  1.D0,   3.D0,   6.D0,   -4.D0,   0.D0, -2.D0, &
                   3.D0, -1.D0,  -6.D0]
      H90(67:105) = 0.D0

      G95(1:66) = [0.D0,-29692.D0,-1784.D0,-2200.D0,3070.D0,1681.D0,1335.D0, &
               -2267.D0,  1249.D0,  759.D0,  940.D0, 780.D0, 290.D0,-418.D0, &
                 122.D0,  -214.D0,  352.D0,  235.D0,-118.D0,-166.D0, -17.D0, &
                  68.D0,    67.D0,   68.D0, -170.D0,  -1.D0,  19.D0, -93.D0, &
                  77.D0,   -72.D0,    1.D0,   28.D0,   5.D0,   4.D0,   8.D0, &
                  -2.D0,    25.D0,    6.D0,   -6.D0,  -9.D0, -14.D0,   9.D0, &
                   6.D0,    -5.D0,   -7.D0,    4.D0,   9.D0,   3.D0, -10.D0, &
                   8.D0,    -8.D0,   -1.D0,   10.D0,  -2.D0,  -8.D0,  -3.D0, &
                  -6.D0,     2.D0,   -4.D0,   -1.D0,   4.D0,   2.D0,   2.D0, &
                   5.D0,     1.D0,    0.D0]
      G95(67:105) = 0.D0

      H95(1:66) = [0.D0,  0.D0,5306.D0,  0.D0,-2366.D0,-413.D0,  0.D0, &
                -262.D0,302.D0,-427.D0,  0.D0,  262.D0,-236.D0, 97.D0, &
                -306.D0,  0.D0,  46.D0,165.D0, -143.D0, -55.D0,107.D0, &
                   0.D0,-17.D0,  72.D0, 67.D0,  -58.D0,   1.D0, 36.D0, &
                   0.D0,-69.D0, -25.D0,  4.D0,   24.D0,  17.D0,-24.D0, &
                  -6.D0,  0.D0,  11.D0,-21.D0,    8.D0, -23.D0, 15.D0, &
                  11.D0,-16.D0,  -4.D0,  0.D0,  -20.D0,  15.D0, 12.D0, &
                  -6.D0, -8.D0,   8.D0,  5.D0,   -8.D0,   3.D0,  0.D0, &
                   1.D0,  0.D0,   4.D0,  5.D0,   -5.D0,  -1.D0, -2.D0, &
                   1.D0, -2.D0,  -7.D0]
      H95(67:105) = 0.D0

      G00 = [0.D0,-29619.4D0,-1728.2D0,-2267.7D0,3068.4D0,1670.9D0, &
               1339.6D0,  -2288.D0, 1252.1D0,  714.5D0, 932.3D0, 786.8D0, &
                 250.D0,   -403.D0,  111.3D0, -218.8D0, 351.4D0, 222.3D0, &
               -130.4D0,  -168.6D0,  -12.9D0,   72.3D0,  68.2D0,  74.2D0, &
               -160.9D0,    -5.9D0,   16.9D0,  -90.4D0,  79.0D0, -74.0D0, &
                   0.D0,    33.3D0,    9.1D0,    6.9D0,   7.3D0,  -1.2D0, &
                 24.4D0,     6.6D0,   -9.2D0,   -7.9D0, -16.6D0,   9.1D0, &
                  7.0D0,    -7.9D0,    -7.D0,     5.D0,   9.4D0,    3.D0, &
                - 8.4D0,     6.3D0,   -8.9D0,   -1.5D0,   9.3D0,  -4.3D0, &
                 -8.2D0,    -2.6D0,    -6.D0,    1.7D0,  -3.1D0,  -0.5D0, &
                  3.7D0,      1.D0,     2.D0,    4.2D0,   0.3D0,  -1.1D0, &
                  2.7D0,    -1.7D0,   -1.9D0,    1.5D0,  -0.1D0,   0.1D0, &
                 -0.7D0,     0.7D0,    1.7D0,    0.1D0,   1.2D0,   4.0D0, &
                 -2.2D0,    -0.3D0,    0.2D0,    0.9D0,  -0.2D0,   0.9D0, &
                 -0.5D0,     0.3D0,   -0.3D0,   -0.4D0,  -0.1D0,  -0.2D0, &
                 -0.4D0,    -0.2D0,   -0.9D0,    0.3D0,   0.1D0,  -0.4D0, &
                  1.3D0,    -0.4D0,    0.7D0,   -0.4D0,   0.3D0,  -0.1D0, &
                  0.4D0,      0.D0,    0.1D0]

      H00 = [0.D0,   0.D0,5186.1D0,   0.D0,-2481.6D0,-458.0D0,   0.D0, &
               -227.6D0,293.4D0,-491.1D0,   0.D0,  272.6D0,-231.9D0,119.8D0, &
               -303.8D0,   0.D0,  43.8D0,171.9D0, -133.1D0, -39.3D0,106.3D0, &
                   0.D0,-17.4D0,  63.7D0, 65.1D0,  -61.2D0,   0.7D0, 43.8D0, &
                   0.D0,-64.6D0, -24.2D0,  6.2D0,    24.D0,  14.8D0,-25.4D0, &
                 -5.8D0,  0.0D0,  11.9D0,-21.5D0,    8.5D0, -21.5D0, 15.5D0, &
                  8.9D0,-14.9D0,  -2.1D0,  0.0D0,  -19.7D0,  13.4D0, 12.5D0, &
                 -6.2D0, -8.4D0,   8.4D0,  3.8D0,   -8.2D0,   4.8D0,  0.0D0, &
                  1.7D0,  0.0D0,   4.0D0,  4.9D0,   -5.9D0,  -1.2D0, -2.9D0, &
                  0.2D0, -2.2D0,  -7.4D0,  0.0D0,    0.1D0,   1.3D0, -0.9D0, &
                 -2.6D0,  0.9D0,  -0.7D0, -2.8D0,   -0.9D0,  -1.2D0, -1.9D0, &
                 -0.9D0,  0.0D0,  -0.4D0,  0.3D0,    2.5D0,  -2.6D0,  0.7D0, &
                  0.3D0,  0.0D0,   0.0D0,  0.3D0,   -0.9D0,  -0.4D0,  0.8D0, &
                  0.0D0, -0.9D0,   0.2D0,  1.8D0,   -0.4D0,  -1.0D0, -0.1D0, &
                  0.7D0,  0.3D0,   0.6D0,  0.3D0,   -0.2D0,  -0.5D0, -0.9D0]

      G05 = [0.D0,  -29554.6D0,-1669.0D0,-2337.2D0,3047.7D0,1657.8D0, &
               1336.3D0,   -2305.8D0, 1246.4D0,  672.5D0, 920.6D0, 798.0D0, &
                210.7D0,    -379.9D0,  100.0D0, -227.0D0, 354.4D0, 208.9D0, &
               -136.5D0,    -168.1D0,  -13.6D0,   73.6D0,  69.6D0,  76.7D0, &
               -151.3D0,     -14.6D0,   14.6D0,  -86.4D0,  79.9D0, -74.5D0, &
                 -1.7D0,      38.7D0,   12.3D0,    9.4D0,   5.4D0,   1.9D0, &
                 24.8D0,       7.6D0,  -11.7D0,   -6.9D0, -18.1D0,  10.2D0, &
                  9.4D0,     -11.3D0,   -4.9D0,    5.6D0,   9.8D0,   3.6D0, &
                 -6.9D0,       5.0D0,  -10.8D0,   -1.3D0,   8.8D0,  -6.7D0, &
                 -9.2D0,      -2.2D0,   -6.1D0,    1.4D0,  -2.4D0,  -0.2D0, &
                  3.1D0,       0.3D0,    2.1D0,    3.8D0,  -0.2D0,  -2.1D0, &
                  2.9D0,      -1.6D0,   -1.9D0,    1.4D0,  -0.3D0,   0.3D0, &
                 -0.8D0,       0.5D0,    1.8D0,    0.2D0,   1.0D0,   4.0D0, &
                 -2.2D0,      -0.3D0,    0.2D0,    0.9D0,  -0.4D0,   1.0D0, &
                 -0.3D0,       0.5D0,   -0.4D0,   -0.4D0,   0.1D0,  -0.5D0, &
                 -0.1D0,      -0.2D0,   -0.9D0,    0.3D0,   0.3D0,  -0.4D0, &
                  1.2D0,      -0.4D0,    0.8D0,   -0.3D0,   0.4D0,  -0.1D0, &
                  0.4D0,      -0.1D0,   -0.2D0]

      H05 = [0.D0,  0.0D0,5078.0D0,  0.0D0,-2594.5D0,-515.4D0,  0.0D0, &
               -198.9D0,269.7D0,-524.7D0,  0.0D0,  282.1D0,-225.2D0,145.2D0, &
               -305.4D0,  0.0D0,  42.7D0,180.3D0, -123.5D0, -19.6D0,103.9D0, &
                  0.0D0,-20.3D0,  54.8D0, 63.6D0,  -63.5D0,   0.2D0, 50.9D0, &
                  0.0D0,-61.1D0, -22.6D0,  6.8D0,   25.4D0,  10.9D0,-26.3D0, &
                 -4.6D0,  0.0D0,  11.2D0,-20.9D0,    9.8D0, -19.7D0, 16.2D0, &
                  7.6D0,-12.8D0,  -0.1D0,  0.0D0,  -20.1D0,  12.7D0, 12.7D0, &
                 -6.7D0, -8.2D0,   8.1D0,  2.9D0,   -7.7D0,   6.0D0,  0.0D0, &
                  2.2D0,  0.1D0,   4.5D0,  4.8D0,   -6.7D0,  -1.0D0, -3.5D0, &
                 -0.9D0, -2.3D0,  -7.9D0,  0.0D0,    0.3D0,   1.4D0, -0.8D0, &
                 -2.3D0,  0.9D0,  -0.6D0, -2.7D0,   -1.1D0,  -1.6D0, -1.9D0, &
                 -1.4D0,  0.0D0,  -0.6D0,  0.2D0,    2.4D0,  -2.6D0,  0.6D0, &
                  0.4D0,  0.0D0,   0.0D0,  0.3D0,   -0.9D0,  -0.3D0,  0.9D0, &
                  0.0D0, -0.8D0,   0.3D0,  1.7D0,   -0.5D0,  -1.1D0,  0.0D0, &
                  0.6D0,  0.2D0,   0.5D0,  0.4D0,   -0.2D0,  -0.6D0, -0.9D0]

      G10 = [0.00D0,-29496.57D0,-1586.42D0,-2396.06D0,3026.34D0, &
                1668.17D0,  1339.85D0,-2326.54D0, 1232.10D0, 633.73D0, &
                 912.66D0,   808.97D0,  166.58D0, -356.83D0,  89.40D0, &
                -230.87D0,   357.29D0,  200.26D0, -141.05D0,-163.17D0, &
                  -8.03D0,    72.78D0,   68.69D0,   75.92D0,-141.40D0, &
                 -22.83D0,    13.10D0,  -78.09D0,   80.44D0, -75.00D0, &
                  -4.55D0,    45.24D0,   14.00D0,   10.46D0,   1.64D0, &
                   4.92D0,    24.41D0,    8.21D0,  -14.50D0,  -5.59D0, &
                 -19.34D0,    11.61D0,   10.85D0,  -14.05D0,  -3.54D0, &
                   5.50D0,     9.45D0,    3.45D0,   -5.27D0,   3.13D0, &
                 -12.38D0,    -0.76D0,    8.43D0,   -8.42D0, -10.08D0, &
                  -1.94D0,    -6.24D0,    0.89D0,   -1.07D0,  -0.16D0, &
                   2.45D0,    -0.33D0,    2.13D0,    3.09D0,  -1.03D0, &
                  -2.80D0,     3.05D0,   -1.48D0,   -2.03D0,   1.65D0, &
                  -0.51D0,     0.54D0,   -0.79D0,    0.37D0,   1.79D0, &
                   0.12D0,     0.75D0,    3.75D0,   -2.12D0,  -0.21D0, &
                   0.30D0,     1.04D0,   -0.63D0,    0.95D0,  -0.11D0, &
                   0.52D0,    -0.39D0,   -0.37D0,    0.21D0,  -0.77D0, &
                   0.04D0,    -0.09D0,   -0.89D0,    0.31D0,   0.42D0, &
                  -0.45D0,     1.08D0,   -0.31D0,    0.78D0,  -0.18D0, &
                   0.38D0,     0.02D0,    0.42D0,   -0.26D0,  -0.26D0]

      H10 = [0.00D0,  0.00D0,4944.26D0,   0.00D0,-2708.54D0, &
               -575.73D0,   0.00D0,-160.40D0, 251.75D0, -537.03D0, 0.00D0, &
                286.48D0,-211.03D0, 164.46D0,-309.72D0,    0.00D0,44.58D0, &
                189.01D0,-118.06D0,  -0.01D0, 101.04D0,    0.00D0,-20.90D0, &
                 44.18D0,  61.54D0, -66.26D0,   3.02D0,   55.40D0,  0.00D0, &
                -57.80D0, -21.20D0,   6.54D0,  24.96D0,    7.03D0,-27.61D0, &
                 -3.28D0,   0.00D0,  10.84D0, -20.03D0,   11.83D0,-17.41D0, &
                 16.71D0,   6.96D0, -10.74D0,   1.64D0,    0.00D0,-20.54D0, &
                 11.51D0,  12.75D0,  -7.14D0,  -7.42D0,    7.97D0,  2.14D0, &
                 -6.08D0,   7.01D0,   0.00D0,   2.73D0,   -0.10D0,  4.71D0, &
                  4.44D0,  -7.22D0,  -0.96D0,  -3.95D0,   -1.99D0, -1.97D0, &
                 -8.31D0,   0.00D0,   0.13D0,   1.67D0,   -0.66D0, -1.76D0, &
                  0.85D0,  -0.39D0,  -2.51D0,  -1.27D0,   -2.11D0, -1.94D0, &
                 -1.86D0,   0.00D0,  -0.87D0,   0.27D0,    2.13D0, -2.49D0, &
                  0.49D0,   0.59D0,   0.00D0,   0.13D0,    0.27D0, -0.86D0, &
                 -0.23D0,   0.87D0,   0.00D0,  -0.87D0,    0.30D0,  1.66D0, &
                 -0.59D0,  -1.14D0,  -0.07D0,   0.54D0,    0.10D0,  0.49D0, &
                  0.44D0,  -0.25D0,  -0.53D0,  -0.79D0]

      G15 = [0.00D0,-29441.46D0,-1501.77D0,-2445.88D0,3012.20D0, &
           1676.35D0,1350.33D0,-2352.26D0,1225.85D0, 581.69D0,907.42D0, &
            813.68D0, 120.49D0,-334.85D0,   70.38D0,-232.91D0,360.14D0, &
            192.35D0,-140.94D0,-157.40D0,    4.30D0,  69.55D0, 67.57D0, &
             72.79D0,-129.85D0, -28.93D0,   13.14D0, -70.85D0, 81.29D0, &
            -75.99D0,  -6.79D0,  51.82D0,   15.07D0,   9.32D0, -2.88D0, &
              6.61D0,  23.98D0,   8.89D0,  -16.78D0,  -3.16D0,-20.56D0, &
             13.33D0,  11.76D0, -15.98D0,   -2.02D0,   5.33D0,  8.83D0, &
              3.02D0,  -3.22D0,   0.67D0,  -13.20D0,  -0.10D0,  8.68D0, &
             -9.06D0, -10.54D0,  -2.01D0,   -6.26D0,   0.17D0,  0.55D0, &
             -0.55D0,   1.70D0,  -0.67D0,    2.13D0,   2.33D0, -1.80D0, &
             -3.59D0,   3.00D0,  -1.40D0,   -2.30D0,   2.08D0, -0.79D0, &
              0.58D0,  -0.70D0,   0.14D0,    1.70D0,  -0.22D0,  0.44D0, &
              3.49D0,  -2.09D0,  -0.16D0,    0.46D0,   1.23D0, -0.89D0, &
              0.85D0,   0.10D0,   0.54D0,   -0.37D0,  -0.43D0,  0.22D0, &
             -0.94D0,  -0.03D0,  -0.02D0,   -0.92D0,   0.42D0,  0.63D0, &
             -0.42D0,   0.96D0,  -0.19D0,    0.81D0,  -0.13D0,  0.38D0, &
              0.08D0,   0.46D0,  -0.35D0,  -0.36D0]

      H15 = [0.00D0,0.00D0,4795.99D0,0.00D0,-2845.41D0,-642.17D0,0.D0, &
          -115.29D0,245.04D0,-538.7D0,   0.0D0, 283.54D0,-188.43D0,180.95D0, &
          -329.23D0,   0.0D0, 46.98D0,196.98D0,-119.14D0,  15.98D0,100.12D0, &
             0.00D0,-20.61D0, 33.30D0, 58.74D0, -66.64D0,   7.35D0, 62.41D0, &
             0.00D0,-54.27D0,-19.53D0,  5.59D0,  24.45D0,   3.27D0,-27.50D0, &
            -2.32D0,  0.00D0, 10.04D0,-18.26D0,  13.18D0, -14.60D0, 16.16D0, &
             5.69D0, -9.10D0,  2.26D0,  0.00D0, -21.77D0,  10.76D0, 11.74D0, &
            -6.74D0, -6.88D0,  7.79D0,  1.04D0,  -3.89D0,   8.44D0,  0.00D0, &
             3.28D0, -0.40D0,  4.55D0,  4.40D0,  -7.92D0,  -0.61D0, -4.16D0, &
            -2.85D0, -1.12D0, -8.72D0,  0.00D0,   0.00D0,   2.11D0, -0.60D0, &
            -1.05D0,  0.76D0, -0.20D0, -2.12D0,  -1.44D0,  -2.57D0, -2.01D0, &
            -2.34D0,  0.00D0, -1.08D0,  0.37D0,   1.75D0,  -2.19D0,  0.27D0, &
             0.72D0, -0.09D0,  0.29D0,  0.23D0,  -0.89D0,  -0.16D0,  0.72D0, &
             0.00D0, -0.88D0,  0.49D0,  1.56D0,  -0.50D0,  -1.24D0, -0.10D0, &
             0.42D0, -0.04D0,  0.48D0,  0.48D0,  -0.30D0,  -0.43D0, -0.71D0]

      G20 = [0.0D0,-29404.8D0,-1450.9D0,-2499.6D0,2982.0D0,1677.0D0, &
           1363.2D0,-2381.2D0,1236.2D0, 525.7D0, 903.0D0, 809.5D0,  86.3D0, &
           -309.4D0,   48.0D0,-234.3D0, 363.2D0, 187.8D0,-140.7D0,-151.2D0, &
             13.5D0,   66.0D0,  65.5D0,  72.9D0,-121.5D0, -36.2D0,  13.5D0, &
            -64.7D0,   80.6D0, -76.7D0,  -8.2D0,  56.5D0,  15.8D0,   6.4D0, &
             -7.2D0,    9.8D0,  23.7D0,   9.7D0, -17.6D0,  -0.5D0, -21.1D0, &
             15.3D0,   13.7D0, -16.5D0,  -0.3D0,   5.0D0,   8.4D0,   2.9D0, &
             -1.5D0,   -1.1D0, -13.2D0,   1.1D0,   8.8D0,  -9.3D0, -11.9D0, &
             -1.9D0,   -6.2D0,  -0.1D0,   1.7D0,  -0.9D0,   0.7D0,  -0.9D0, &
              1.9D0,    1.4D0,  -2.4D0,  -3.8D0,   3.0D0,  -1.4D0,  -2.5D0, &
              2.3D0,   -0.9D0,   0.3D0,  -0.7D0,  -0.1D0,   1.4D0,  -0.6D0, &
              0.2D0,    3.1D0,  -2.0D0,  -0.1D0,   0.5D0,   1.3D0,  -1.2D0, &
              0.7D0,    0.3D0,   0.5D0,  -0.3D0,  -0.5D0,   0.1D0,  -1.1D0, &
             -0.3D0,    0.1D0,  -0.9D0,   0.5D0,   0.7D0,  -0.3D0,   0.8D0, &
              0.0D0,    0.8D0,   0.0D0,   0.4D0,   0.1D0,   0.5D0,  -0.5D0, &
             -0.4D0]

      H20 = [0.0D0,0.0D0,4652.5D0,0.0D0,-2991.6D0,-734.6D0,   0.0D0, &
           -82.1D0,241.9D0,-543.4D0,   0.0D0,   281.9D0, -158.4D0,  199.7D0, &
          -349.7D0,  0.0D0,  47.7D0, 208.3D0,  -121.2D0,   32.3D0,   98.9D0, &
             0.0D0,-19.1D0,  25.1D0,  52.8D0,   -64.5D0,    8.9D0,   68.1D0, &
             0.0D0,-51.5D0, -16.9D0,   2.2D0,    23.5D0,   -2.2D0,  -27.2D0, &
            -1.8D0,  0.0D0,   8.4D0, -15.3D0,    12.8D0,  -11.7D0,   14.9D0, &
             3.6D0, -6.9D0,   2.8D0,   0.0D0,   -23.4D0,   11.0D0,    9.8D0, &
            -5.1D0, -6.3D0,   7.8D0,   0.4D0,    -1.4D0,    9.6D0,    0.0D0, &
             3.4D0, -0.2D0,   3.6D0,   4.8D0,    -8.6D0,   -0.1D0,   -4.3D0, &
            -3.4D0, -0.1D0,  -8.8D0,   0.0D0,     0.0D0,    2.5D0,   -0.6D0, &
            -0.4D0,  0.6D0,  -0.2D0,  -1.7D0,    -1.6D0,   -3.0D0,   -2.0D0, &
            -2.6D0,  0.0D0,  -1.2D0,   0.5D0,     1.4D0,   -1.8D0,    0.1D0, &
             0.8D0, -0.2D0,   0.6D0,   0.2D0,    -0.9D0,    0.0D0,    0.5D0, &
             0.0D0, -0.9D0,   0.6D0,   1.4D0,    -0.4D0,   -1.3D0,   -0.1D0, &
             0.3D0, -0.1D0,   0.5D0,   0.5D0,    -0.4D0,   -0.4D0,   -0.6D0]

      DG20 = [0.0D0, 5.7D0,  7.4D0,-11.0D0,-7.0D0,  -2.1D0,     2.2D0, &
                   -5.9D0, 3.1D0,-12.0D0, -1.2D0,-1.6D0,  -5.9D0,     5.2D0, &
                   -5.1D0,-0.3D0,  0.5D0, -0.6D0, 0.2D0,   1.3D0,     0.9D0, &
                   -0.5D0,-0.3D0,  0.4D0,  1.3D0,-1.4D0,   0.0D0,     0.9D0, &
                   -0.1D0,-0.2D0,  0.0D0,  0.7D0, 0.1D0,  -0.5D0,    -0.8D0, &
                    0.8D0, 0.0D0,  0.1D0, -0.1D0, 0.4D0,  -0.1D0,     0.4D0, &
                    0.3D0,-0.1D0,  0.4D0]

      DH20 = [0.0D0, 0.0D0,-25.9D0, 0.0D0,-30.2D0, -22.4D0,   0.0D0, &
                    6.0D0,-1.1D0,  0.5D0, 0.0D0, -0.1D0,   6.5D0,   3.6D0, &
                   -5.0D0, 0.0D0,  0.0D0, 2.5D0, -0.6D0,   3.0D0,   0.3D0, &
                    0.0D0, 0.0D0, -1.6D0,-1.3D0,  0.8D0,   0.0D0,   1.0D0, &
                    0.0D0, 0.6D0,  0.6D0,-0.8D0, -0.2D0,  -1.1D0,   0.1D0, &
                    0.3D0, 0.0D0, -0.2D0, 0.6D0, -0.2D0,   0.5D0,  -0.3D0, &
                   -0.4D0, 0.5D0,  0.0D0]

      IY=IYEAR

!  WE ARE RESTRICTED BY THE INTERVAL 1965-2025, FOR WHICH EITHER THE IGRF/DGRF COEFFICIENTS OR SECULAR VELOCITIES
!    ARE KNOWN; IF IYEAR IS OUTSIDE THIS INTERVAL, THEN THE SUBROUTINE USES THE
!      NEAREST LIMITING VALUE AND PRINTS A WARNING:

      IF(IY<1965) THEN
       IY=1965
       WRITE(*,*) ' **** RECALC_08 WARNS: YEAR IS OUT OF INTERVAL 1965-2025: IYEAR=',IYEAR
       WRITE(*,*) '      CALCULATIONS WILL BE DONE FOR IYEAR=',IY
!       WRITE (*,10) IYEAR,IY
      ENDIF

      IF(IY>2025) THEN
       IY=2025
       WRITE(*,*) ' **** RECALC_08 WARNS: YEAR IS OUT OF INTERVAL 1965-2025: IYEAR=',IYEAR
       WRITE(*,*) '      CALCULATIONS WILL BE DONE FOR IYEAR=',IY
!       WRITE (*,10) IYEAR,IY
      ENDIF

!  CALCULATE THE ARRAY REC, CONTAINING COEFFICIENTS FOR THE RECURSION RELATIONS,
!  USED IN THE IGRF SUBROUTINE FOR CALCULATING THE ASSOCIATE LEGENDRE POLYNOMIALS
!  AND THEIR DERIVATIVES:

      DO N=1,14
         N2=2*N-1
         N2=N2*(N2-2)
         DO M=1,N
            MN=N*(N-1)/2+M
            REC(MN)=DFLOAT((N-M)*(N+M-2))/DFLOAT(N2)
         ENDDO
      ENDDO
!<<<>>>
!      DO 20 N=1,14
!         N2=2*N-1
!         N2=N2*(N2-2)
!         DO 20 M=1,N
!            MN=N*(N-1)/2+M
!20    REC(MN)=DFLOAT((N-M)*(N+M-2))/DFLOAT(N2)

!      IF (IY<1970) GOTO 50          !INTERPOLATE BETWEEN 1965 - 1970
!      IF (IY<1975) GOTO 60          !INTERPOLATE BETWEEN 1970 - 1975
!      IF (IY<1980) GOTO 70          !INTERPOLATE BETWEEN 1975 - 1980
!      IF (IY<1985) GOTO 80          !INTERPOLATE BETWEEN 1980 - 1985
!      IF (IY<1990) GOTO 90          !INTERPOLATE BETWEEN 1985 - 1990
!      IF (IY<1995) GOTO 100         !INTERPOLATE BETWEEN 1990 - 1995
!      IF (IY<2000) GOTO 110         !INTERPOLATE BETWEEN 1995 - 2000
!      IF (IY<2005) GOTO 120         !INTERPOLATE BETWEEN 2000 - 2005
!      IF (IY<2010) GOTO 130         !INTERPOLATE BETWEEN 2005 - 2010
!      IF (IY<2015) GOTO 140         !INTERPOLATE BETWEEN 2010 - 2015
!      IF (IY<2020) GOTO 150         !INTERPOLATE BETWEEN 2015 - 2020

!       INTERPOLATE BETWEEEN 1965 - 1970:
      IF (IY<1970) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1965)/5.D0
          F1=1.D0-F2
          DO N=1,105
             G(N)=G65(N)*F1+G70(N)*F2
             H(N)=H65(N)*F1+H70(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!50    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1965)/5.D0
!      F1=1.D0-F2
!      DO 55 N=1,105
!         G(N)=G65(N)*F1+G70(N)*F2
!55       H(N)=H65(N)*F1+H70(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 1970 - 1975:
      ELSEIF (IY<1975) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1970)/5.D0
          F1=1.D0-F2
          DO N=1,105
             G(N)=G70(N)*F1+G75(N)*F2
             H(N)=H70(N)*F1+H75(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!60    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1970)/5.D0
!      F1=1.D0-F2
!      DO 65 N=1,105
!         G(N)=G70(N)*F1+G75(N)*F2
!65       H(N)=H70(N)*F1+H75(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 1975 - 1980:
      ELSEIF (IY<1980) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1975)/5.D0
          F1=1.D0-F2
          DO N=1,105
             G(N)=G75(N)*F1+G80(N)*F2
             H(N)=H75(N)*F1+H80(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!70    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1975)/5.D0
!      F1=1.D0-F2
!      DO 75 N=1,105
!         G(N)=G75(N)*F1+G80(N)*F2
!75       H(N)=H75(N)*F1+H80(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 1980 - 1985:
      ELSEIF (IY<1985) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1980)/5.D0
          F1=1.D0-F2
          DO N=1,105
             G(N)=G80(N)*F1+G85(N)*F2
             H(N)=H80(N)*F1+H85(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!80    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1980)/5.D0
!      F1=1.D0-F2
!      DO 85 N=1,105
!         G(N)=G80(N)*F1+G85(N)*F2
!85       H(N)=H80(N)*F1+H85(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 1985 - 1990:
      ELSEIF (IY<1990) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1985)/5.D0
          F1=1.D0-F2
          DO N=1,105
             G(N)=G85(N)*F1+G90(N)*F2
             H(N)=H85(N)*F1+H90(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!90    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1985)/5.D0
!      F1=1.D0-F2
!      DO 95 N=1,105
!         G(N)=G85(N)*F1+G90(N)*F2
!95       H(N)=H85(N)*F1+H90(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 1990 - 1995:
      ELSEIF (IY<1995) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1990)/5.D0
          F1=1.D0-F2
          DO N=1,105
             G(N)=G90(N)*F1+G95(N)*F2
             H(N)=H90(N)*F1+H95(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!100   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1990)/5.D0
!      F1=1.D0-F2
!      DO 105 N=1,105
!         G(N)=G90(N)*F1+G95(N)*F2
!105      H(N)=H90(N)*F1+H95(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 1995 - 2000:
      ELSEIF (IY<2000) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1995)/5.D0
          F1=1.D0-F2
          DO N=1,105   !  THE 2000 COEFFICIENTS (G00) GO THROUGH THE ORDER 13, NOT 10
             G(N)=G95(N)*F1+G00(N)*F2
             H(N)=H95(N)*F1+H00(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!110   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1995)/5.D0
!      F1=1.D0-F2
!      DO 115 N=1,105   !  THE 2000 COEFFICIENTS (G00) GO THROUGH THE ORDER 13, NOT 10
!         G(N)=G95(N)*F1+G00(N)*F2
!115      H(N)=H95(N)*F1+H00(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 2000 - 2005:
      ELSEIF (IY<2005) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-2000)/5.D0
          F1=1.D0-F2
          DO N=1,105
             G(N)=G00(N)*F1+G05(N)*F2
             H(N)=H00(N)*F1+H05(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!120   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-2000)/5.D0
!      F1=1.D0-F2
!      DO 125 N=1,105
!         G(N)=G00(N)*F1+G05(N)*F2
!125      H(N)=H00(N)*F1+H05(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 2005 - 2010:
      ELSEIF (IY<2010) THEN
          F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-2005)/5.
          F1=1.-F2
          DO N=1,105
             G(N)=G05(N)*F1+G10(N)*F2
             H(N)=H05(N)*F1+H10(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!130   F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-2005)/5.
!      F1=1.-F2
!      DO 135 N=1,105
!         G(N)=G05(N)*F1+G10(N)*F2
!135      H(N)=H05(N)*F1+H10(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 2010 - 2015:
      ELSEIF (IY<2015) THEN
          F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25-2010)/5.
          F1=1.-F2
          DO N=1,105
             G(N)=G10(N)*F1+G15(N)*F2
             H(N)=H10(N)*F1+H15(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!140   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25-2010)/5.
!      F1=1.-F2
!      DO 145 N=1,105
!         G(N)=G10(N)*F1+G15(N)*F2
!145      H(N)=H10(N)*F1+H15(N)*F2
!      GOTO 300

!       INTERPOLATE BETWEEN 2015 - 2020:
      ELSEIF (IY<2020) THEN
          F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-2015)/5.
          F1=1.-F2
          DO N=1,105
             G(N)=G15(N)*F1+G20(N)*F2
             H(N)=H15(N)*F1+H20(N)*F2
          ENDDO
!      GOTO 300
!<<<>>>
!150   F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-2015)/5.
!      F1=1.-F2
!      DO 155 N=1,105
!         G(N)=G15(N)*F1+G20(N)*F2
!155      H(N)=H15(N)*F1+H20(N)*F2
!      GOTO 300

      ELSE!IF (IY>=2020) THEN
!       EXTRAPOLATE BEYOND 2020:
          DT=DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-2020.D0
          DO N=1,105
             G(N)=G20(N)
             H(N)=H20(N)
             IF (N>45) CYCLE
             G(N)=G(N)+DG20(N)*DT
             H(N)=H(N)+DH20(N)*DT
          ENDDO
!      GOTO 300
!<<<>>>
!      DO 40 N=1,105
!         G(N)=G20(N)
!         H(N)=H20(N)
!         IF (N>45) GOTO 40
!         G(N)=G(N)+DG20(N)*DT
!         H(N)=H(N)+DH20(N)*DT
!40    CONTINUE
!      GOTO 300

      ENDIF

!   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
!   THEM BY SCHMIDT NORMALIZATION FACTORS:

          S=1.D0
    !300   S=1.D0
          DO N=2,14
             MN=N*(N-1)/2+1
             S=S*DFLOAT(2*N-3)/DFLOAT(N-1)
             G(MN)=G(MN)*S
             H(MN)=H(MN)*S
             P=S
             DO M=2,N
                AA=1.D0
                IF (M.EQ.2) AA=2.D0
                P=P*DSQRT(AA*DFLOAT(N-M+1)/DFLOAT(N+M-2))
                MNN=MN+M-1
                G(MNN)=G(MNN)*P
                H(MNN)=H(MNN)*P
    !250         H(MNN)=H(MNN)*P
             ENDDO
          ENDDO

           G_10=-G(2)
           G_11= G(3)
           H_11= H(3)
!
!  NOW CALCULATE GEO COMPONENTS OF THE UNIT VECTOR EzMAG, PARALLEL TO GEODIPOLE AXIS:
!   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
!         ST0 * CL0                ST0 * SL0                CT0
!
          SQ=G_11**2+H_11**2
          SQQ=DSQRT(SQ)
          SQR=DSQRT(G_10**2+SQ)
          SL0=-H_11/SQQ
          CL0=-G_11/SQQ
          ST0=SQQ/SQR
          CT0=G_10/SQR
          STCL=ST0*CL0
          STSL=ST0*SL0
          CTSL=CT0*SL0
          CTCL=CT0*CL0
!
!  NOW CALCULATE GEI COMPONENTS (S1,S2,S3) OF THE UNIT VECTOR S = EX_GSE
!    POINTING FROM THE EARTH''S CENTER TO SUN
!
          CALL SUN_08 (IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
    !
          S1=DCOS(SRASN)*DCOS(SDEC)
          S2=DSIN(SRASN)*DCOS(SDEC)
          S3=DSIN(SDEC)
!
!  NOW CALCULATE GEI COMPONENTS (DZ1,DZ2,DZ3) OF THE UNIT VECTOR EZGSE
!  POINTING NORTHWARD AND ORTHOGONAL TO THE ECLIPTIC PLANE, AS
!  (0,-SIN(OBLIQ),COS(OBLIQ)). FOR THE EPOCH 1978, OBLIQ = 23.44214 DEGS.
!  HERE WE USE A MORE ACCURATE TIME-DEPENDENT VALUE, DETERMINED AS:
!
          DJ=DFLOAT(365*(IY-1900)+(IY-1901)/4 +IDAY) &
             -0.5+DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
          T=DJ/36525.D0
          OBLIQ=(23.45229D0-0.0130125D0*T)/57.2957795D0
          DZ1=0.D0
          DZ2=-DSIN(OBLIQ)
          DZ3=DCOS(OBLIQ)
!
!  NOW OBTAIN GEI COMPONENTS OF THE UNIT VECTOR EYGSE=(DY1,DY2,DY3),
!  COMPLETING THE RIGHT-HANDED SYSTEM. THEY CAN BE FOUND FROM THE VECTOR
!  PRODUCT EZGSE x EXGSE = (DZ1,DZ2,DZ3) x (S1,S2,S3):

          DY1=DZ2*S3-DZ3*S2
          DY2=DZ3*S1-DZ1*S3
          DY3=DZ1*S2-DZ2*S1

!  NOW CALCULATE GEI COMPONENTS OF THE UNIT VECTOR X = EXGSW, DIRECTED ANTIPARALLEL
!  TO THE OBSERVED SOLAR WIND FLOW. FIRST, CALCULATE ITS COMPONENTS IN GSE:

          V=DSQRT(VGSEX**2+VGSEY**2+VGSEZ**2)
          DX1=-VGSEX/V
          DX2=-VGSEY/V
          DX3=-VGSEZ/V

!  THEN IN GEI:

          X1=DX1*S1+DX2*DY1+DX3*DZ1
          X2=DX1*S2+DX2*DY2+DX3*DZ2
          X3=DX1*S3+DX2*DY3+DX3*DZ3

!  NOW CALCULATE GEI COMPONENTS (DIP1,DIP2,DIP3) OF THE UNIT VECTOR DIP = EZ_SM = EZ_MAG,
!   ALIGNED WITH THE GEODIPOLE AND POINTING NORTHWARD FROM ECLIPTIC PLANE:

          CGST=DCOS(GST)
          SGST=DSIN(GST)

          DIP1=STCL*CGST-STSL*SGST
          DIP2=STCL*SGST+STSL*CGST
          DIP3=CT0

!  THIS ALLOWS US TO CALCULATE GEI COMPONENTS OF THE UNIT VECTOR Y = EYGSW
!   BY TAKING THE VECTOR PRODUCT DIP x X AND NORMALIZING IT TO UNIT LENGTH:

          Y1=DIP2*X3-DIP3*X2
          Y2=DIP3*X1-DIP1*X3
          Y3=DIP1*X2-DIP2*X1
          Y=DSQRT(Y1*Y1+Y2*Y2+Y3*Y3)
          Y1=Y1/Y
          Y2=Y2/Y
          Y3=Y3/Y

!   AND GEI COMPONENTS OF THE UNIT VECTOR Z = EZGSW = EXGSW x EYGSW = X x Y:

          Z1=X2*Y3-X3*Y2
          Z2=X3*Y1-X1*Y3
          Z3=X1*Y2-X2*Y1

!   ELEMENTS OF THE MATRIX GSE TO GSW ARE THE SCALAR PRODUCTS:

!  E11=(EXGSE,EXGSW)  E12=(EXGSE,EYGSW)  E13=(EXGSE,EZGSW)
!  E21=(EYGSE,EXGSW)  E22=(EYGSE,EYGSW)  E23=(EYGSE,EZGSW)
!  E31=(EZGSE,EXGSW)  E32=(EZGSE,EYGSW)  E33=(EZGSE,EZGSW)

          E11= S1*X1 +S2*X2 +S3*X3
          E12= S1*Y1 +S2*Y2 +S3*Y3
          E13= S1*Z1 +S2*Z2 +S3*Z3
          E21=DY1*X1+DY2*X2+DY3*X3
          E22=DY1*Y1+DY2*Y2+DY3*Y3
          E23=DY1*Z1+DY2*Z2+DY3*Z3
          E31=DZ1*X1+DZ2*X2+DZ3*X3
          E32=DZ1*Y1+DZ2*Y2+DZ3*Y3
          E33=DZ1*Z1+DZ2*Z2+DZ3*Z3

!   GEODIPOLE TILT ANGLE IN THE GSW SYSTEM: PSI=ARCSIN(DIP,EXGSW)

          SPS=DIP1*X1+DIP2*X2+DIP3*X3
          CPS=DSQRT(1.D0-SPS**2)
          PSI=DASIN(SPS)

!   ELEMENTS OF THE MATRIX GEO TO GSW ARE THE SCALAR PRODUCTS:

!   A11=(EXGEO,EXGSW), A12=(EYGEO,EXGSW), A13=(EZGEO,EXGSW),
!   A21=(EXGEO,EYGSW), A22=(EYGEO,EYGSW), A23=(EZGEO,EYGSW),
!   A31=(EXGEO,EZGSW), A32=(EYGEO,EZGSW), A33=(EZGEO,EZGSW),

!   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:

!  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
!  EXGSW=(X1,X2,X3),  EYGSW=(Y1,Y2,Y3),   EZGSW=(Z1,Z2,Z3)
!                                                           AND  THEREFORE:

          A11=X1*CGST+X2*SGST
          A12=-X1*SGST+X2*CGST
          A13=X3
          A21=Y1*CGST+Y2*SGST
          A22=-Y1*SGST+Y2*CGST
          A23=Y3
          A31=Z1*CGST+Z2*SGST
          A32=-Z1*SGST+Z2*CGST
          A33=Z3

!  NOW CALCULATE ELEMENTS OF THE MATRIX MAG TO SM (ONE ROTATION ABOUT THE GEODIPOLE AXIS);
!   THEY ARE FOUND AS THE SCALAR PRODUCTS: CFI=GM22=(EYSM,EYMAG)=(EYGSW,EYMAG),
!                                          SFI=GM23=(EYSM,EXMAG)=(EYGSW,EXMAG),
!    DERIVED AS FOLLOWS:

! IN GEO, THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
!  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THEIR COMPONENTS ARE:
!  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
!            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
!            -ST0
!  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
!            -SL0*SIN(GST)+CL0*COS(GST)
!             0
!  NOW, NOTE THAT GEI COMPONENTS OF EYSM=EYGSW WERE FOUND ABOVE AS Y1, Y2, AND Y3,
!  AND WE ONLY HAVE TO COMBINE THESE QUANTITIES INTO SCALAR PRODUCTS:

          EXMAGX=CT0*(CL0*CGST-SL0*SGST)
          EXMAGY=CT0*(CL0*SGST+SL0*CGST)
          EXMAGZ=-ST0
          EYMAGX=-(SL0*CGST+CL0*SGST)
          EYMAGY=-(SL0*SGST-CL0*CGST)
          CFI=Y1*EYMAGX+Y2*EYMAGY
          SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ

! 10   FORMAT(//1X, &
!         '**** RECALC_08 WARNS: YEAR IS OUT OF INTERVAL 1965-2025: IYEAR=', &
!         I4,/,6X,'CALCULATIONS WILL BE DONE FOR IYEAR=',I4,/)
          RETURN
      END SUBROUTINE RECALC_08

!=====================================================================================

       SUBROUTINE SMGSW_08 (XSM,YSM,ZSM,XGSW,YGSW,ZGSW,J)

!  CONVERTS SOLAR MAGNETIC (SM) TO GEOCENTRIC SOLAR-WIND (GSW) COORDINATES OR VICE VERSA.

!                  J>0                 J<0
!-----INPUT: J,XSM,YSM,ZSM        J,XGSW,YGSW,ZGSW
!----OUTPUT:  XGSW,YGSW,ZGSW       XSM,YSM,ZSM

!  ATTENTION:  SUBROUTINE RECALC_08 MUST BE INVOKED BEFORE SMGSW_08 IN THREE CASES:
!     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
!     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
!     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED

!    IMPORTANT NOTE:

!        A NON-STANDARD DEFINITION IS IMPLIED HERE FOR THE SOLAR MAGNETIC (SM) COORDINATE
!        SYSTEM:  IT IS ASSUMED THAT THE XSM AXIS LIES IN THE PLANE DEFINED BY THE
!        GEODIPOLE AXIS AND THE OBSERVED VECTOR OF THE SOLAR WIND FLOW (RATHER THAN
!        THE EARTH-SUN LINE).  IN ORDER TO CONVERT MAG COORDINATES TO AND FROM THE
!        STANDARD SM COORDINATES, INVOKE RECALC_08 WITH VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0

!     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)

!     AUTHOR:  N. A. TSYGANENKO

          implicit none
          real(rp), intent(inout) :: XSM,YSM,ZSM,XGSW,YGSW,ZGSW
          integer, intent(in) :: J
!<<<>>>
!      IMPLICIT REAL*8 (A-H,O-Z)
!      COMMON /GEOPACK1/ A(10),SPS,CPS,B(22)

          IF (J>0) THEN
           XGSW=XSM*CPS+ZSM*SPS
           YGSW=YSM
           ZGSW=ZSM*CPS-XSM*SPS
          ELSE
           XSM=XGSW*CPS-ZGSW*SPS
           YSM=YGSW
           ZSM=XGSW*SPS+ZGSW*CPS
          ENDIF

          RETURN
      END SUBROUTINE SMGSW_08


!==========================================================================================

      SUBROUTINE GSWGSE_08 (XGSW,YGSW,ZGSW,XGSE,YGSE,ZGSE,J)
         implicit none
         real(rp), intent(inout) :: XGSW,YGSW,ZGSW,XGSE,YGSE,ZGSE
         integer , intent(in)    :: J

!
!  THIS SUBROUTINE TRANSFORMS COMPONENTS OF ANY VECTOR BETWEEN THE STANDARD GSE
!  COORDINATE SYSTEM AND THE GEOCENTRIC SOLAR-WIND (GSW, aka GSWM), DEFINED AS FOLLOWS
!  (HONES ET AL., PLANET.SPACE SCI., V.34, P.889, 1986; TSYGANENKO ET AL., JGRA,
!  V.103(A4), P.6827, 1998):
!
!  IN THE GSW SYSTEM, X AXIS IS ANTIPARALLEL TO THE OBSERVED DIRECTION OF THE SOLAR WIND FLOW.
!  TWO OTHER AXES, Y AND Z, ARE DEFINED IN THE SAME WAY AS FOR THE STANDARD GSM, THAT IS,
!  Z AXIS ORTHOGONAL TO X AXIS, POINTS NORTHWARD, AND LIES IN THE PLANE DEFINED BY THE X-
!  AND GEODIPOLE AXIS. THE Y AXIS COMPLETES THE RIGHT-HANDED SYSTEM.
!
!  THE GSW SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM IN THE CASE OF
!   A STRICTLY RADIAL SOLAR WIND FLOW.
!
!  AUTHOR:  N. A. TSYGANENKO
!  ADDED TO 2008 VERSION OF GEOPACK: JAN 27, 2008.
!
!                    J>0                       J<0
!-----INPUT:   J,XGSW,YGSW,ZGSW          J,XGSE,YGSE,ZGSE
!-----OUTPUT:    XGSE,YGSE,ZGSE            XGSW,YGSW,ZGSW
!
!  IMPORTANT THINGS TO REMEMBER:
!
!   (1) BEFORE CALLING GSWGSE_08, BE SURE TO INVOKE SUBROUTINE RECALC_08, IN ORDER
!       TO DEFINE ALL NECESSARY ELEMENTS OF TRANSFORMATION MATRICES
!
!   (2) IN THE ABSENCE OF INFORMATION ON THE SOLAR WIND DIRECTION, E.G., WITH ONLY SCALAR
!       SPEED V KNOWN, THIS SUBROUTINE CAN BE USED TO CONVERT VECTORS TO ABERRATED
!       COORDINATE SYSTEM, TAKING INTO ACCOUNT EARTH'S ORBITAL SPEED OF 29 KM/S.
!       TO DO THAT, SPECIFY THE LAST 3 PARAMETERS IN RECALC_08 AS FOLLOWS:
!       VGSEX=-V, VGSEY=29.0, VGSEZ=0.0.
!
!       IT SHOULD ALSO BE KEPT IN MIND THAT IN SOME SOLAR WIND DATABASES THE ABERRATION
!       EFFECT HAS ALREADY BEEN TAKEN INTO ACCOUNT BY SUBTRACTING 29 KM/S FROM VYGSE;
!       IN THAT CASE, THE ORIGINAL VYGSE VALUES SHOULD BE RESTORED BY ADDING BACK THE
!       29 KM/S CORRECTION. WHETHER OR NOT TO DO THAT, MUST BE VERIFIED WITH THE DATA
!       ORIGINATOR (OR CAN BE DETERMINED BY CALCULATING THE AVERAGE VGSEY OVER
!       A SUFFICIENTLY LONG TIME INTERVAL)
!
!   (3) IF NO INFORMATION IS AVAILABLE ON THE SOLAR WIND SPEED, THEN SET VGSEX=-400.0
!       AND  VGSEY=VGSEZ=0. IN THAT CASE, THE GSW COORDINATE SYSTEM BECOMES
!       IDENTICAL TO THE STANDARD ONE.
!
!     COMMON /GEOPACK1/ AAA(25),E11,E21,E31,E12,E22,E32,E13,E23,E33
!
!  DIRECT TRANSFORMATION:
!
      IF (J.GT.0) THEN
        XGSE=XGSW*E11+YGSW*E12+ZGSW*E13
        YGSE=XGSW*E21+YGSW*E22+ZGSW*E23
        ZGSE=XGSW*E31+YGSW*E32+ZGSW*E33
      ENDIF
!
!   INVERSE TRANSFORMATION: CARRIED OUT USING THE TRANSPOSED MATRIX:
!
      IF (J.LT.0) THEN
        XGSW=XGSE*E11+YGSE*E21+ZGSE*E31
        YGSW=XGSE*E12+YGSE*E22+ZGSE*E32
        ZGSW=XGSE*E13+YGSE*E23+ZGSE*E33
      ENDIF

      RETURN
      END

!==========================================================================================

      SUBROUTINE GEOGSW_08 (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW,J)

! CONVERTS GEOGRAPHIC (GEO) TO GEOCENTRIC SOLAR-WIND (GSW) COORDINATES OR VICE VERSA.

!                   J>0                   J<0
!----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSW,YGSW,ZGSW
!---- OUTPUT:    XGSW,YGSW,ZGSW      XGEO,YGEO,ZGEO

!  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEOGSW_08 IN THREE CASES:
!     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES, OR
!     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  HAVE CHANGED, AND/OR
!     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED

!  NOTE: THIS SUBROUTINE CONVERTS GEO VECTORS TO AND FROM THE SOLAR-WIND GSW COORDINATE
!        SYSTEM, TAKING INTO ACCOUNT POSSIBLE DEFLECTIONS OF THE SOLAR WIND DIRECTION FROM
!        STRICTLY RADIAL.  BEFORE CONVERTING TO/FROM STANDARD GSM COORDINATES, INVOKE RECALC_08
!        WITH VGSEX=-400.0 and VGSEY=0.0, VGSEZ=0.0

!     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)

!     AUTHOR:  N. A. TSYGANENKO

          implicit none
          real(rp), intent(inout) :: XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW
          integer, intent(in) :: J
!<<<>>>
!      IMPLICIT REAL*8 (A-H,O-Z)
!      COMMON /GEOPACK1/ AA(16),A11,A21,A31,A12,A22,A32,A13,A23,A33,B(9)

          IF (J>0) THEN
           XGSW=A11*XGEO+A12*YGEO+A13*ZGEO
           YGSW=A21*XGEO+A22*YGEO+A23*ZGEO
           ZGSW=A31*XGEO+A32*YGEO+A33*ZGEO
          ELSE
           XGEO=A11*XGSW+A21*YGSW+A31*ZGSW
           YGEO=A12*XGSW+A22*YGSW+A32*ZGSW
           ZGEO=A13*XGSW+A23*YGSW+A33*ZGSW
          ENDIF

          RETURN
      END SUBROUTINE GEOGSW_08

!==========================================================================================
      SUBROUTINE GEO2GSW (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW)
          implicit none
          real(rp), intent(in) :: XGEO,YGEO,ZGEO
          real(rp), intent(out) :: XGSW,YGSW,ZGSW
          XGSW=A11*XGEO+A12*YGEO+A13*ZGEO
          YGSW=A21*XGEO+A22*YGEO+A23*ZGEO
          ZGSW=A31*XGEO+A32*YGEO+A33*ZGEO
      END SUBROUTINE GEO2GSW

      SUBROUTINE GSW2GEO (XGSW,YGSW,ZGSW,XGEO,YGEO,ZGEO)
          implicit none
          real(rp), intent(in) :: XGSW,YGSW,ZGSW
          real(rp), intent(out) :: XGEO,YGEO,ZGEO
          XGEO=A11*XGSW+A21*YGSW+A31*ZGSW
          YGEO=A12*XGSW+A22*YGSW+A32*ZGSW
          ZGEO=A13*XGSW+A23*YGSW+A33*ZGSW
      END SUBROUTINE GSW2GEO

      SUBROUTINE SM2GSW (XSM,YSM,ZSM,XGSW,YGSW,ZGSW)
          implicit none
          real(rp), intent(in) :: XSM,YSM,ZSM
          real(rp), intent(out) :: XGSW,YGSW,ZGSW
          XGSW=XSM*CPS+ZSM*SPS
          YGSW=YSM
          ZGSW=ZSM*CPS-XSM*SPS
      END SUBROUTINE SM2GSW

      SUBROUTINE GSW2SM (XGSW,YGSW,ZGSW,XSM,YSM,ZSM)
          implicit none
          real(rp), intent(in) :: XGSW,YGSW,ZGSW
          real(rp), intent(out) :: XSM,YSM,ZSM
          XSM=XGSW*CPS-ZGSW*SPS
          YSM=YGSW
          ZSM=XGSW*SPS+ZGSW*CPS
      END SUBROUTINE GSW2SM

      SUBROUTINE GEO2SM_xyz (XGEO,YGEO,ZGEO,XSM,YSM,ZSM)
          implicit none
          real(rp), intent(in) :: XGEO,YGEO,ZGEO
          real(rp), intent(out) :: XSM,YSM,ZSM
          real(rp) :: XGSW,YGSW,ZGSW
          call GEO2GSW (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW)
          call GSW2SM (XGSW,YGSW,ZGSW,XSM,YSM,ZSM)
      END SUBROUTINE GEO2SM_xyz

      SUBROUTINE GEO2SM_vec(xyzGEO,xyzSM)
        implicit none
        real(rp), intent(in)  :: xyzGEO(NDIM)
        real(rp), intent(out) :: xyzSM (NDIM)

        call GEO2SM_xyz(xyzGEO(XDIR),xyzGEO(YDIR),xyzGEO(ZDIR),xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR))
      END SUBROUTINE GEO2SM_vec

      SUBROUTINE SM2GEO_xyz (XSM,YSM,ZSM,XGEO,YGEO,ZGEO)
          implicit none
          real(rp), intent(in) :: XSM,YSM,ZSM
          real(rp), intent(out) :: XGEO,YGEO,ZGEO
          real(rp) :: XGSW,YGSW,ZGSW
          call SM2GSW (XSM,YSM,ZSM,XGSW,YGSW,ZGSW)
          call GSW2GEO (XGSW,YGSW,ZGSW,XGEO,YGEO,ZGEO)
      END SUBROUTINE SM2GEO_xyz

      SUBROUTINE SM2GEO_vec(xyzSM,xyzGEO)
        implicit none
        real(rp), intent(in ) :: xyzSM (NDIM)
        real(rp), intent(out) :: xyzGEO(NDIM)
        call SM2GEO_xyz(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGEO(XDIR),xyzGEO(YDIR),xyzGEO(ZDIR))
      END SUBROUTINE SM2GEO_vec

      !SM=>GSE
      SUBROUTINE SM2GSE_xyz (XSM,YSM,ZSM,XGSE,YGSE,ZGSE)
          implicit none
          real(rp), intent(in) :: XSM,YSM,ZSM
          real(rp), intent(out) :: XGSE,YGSE,ZGSE

          real(rp) :: XGSW,YGSW,ZGSW

          call SM2GSW(XSM,YSM,ZSM,XGSW,YGSW,ZGSW)
          XGSE=XGSW*E11+YGSW*E12+ZGSW*E13
          YGSE=XGSW*E21+YGSW*E22+ZGSW*E23
          ZGSE=XGSW*E31+YGSW*E32+ZGSW*E33
      END SUBROUTINE SM2GSE_xyz
      
      SUBROUTINE SM2GSE_vec(xyzSM,xyzGSE)
        implicit none
        real(rp), intent(in ) :: xyzSM (NDIM)
        real(rp), intent(out) :: xyzGSE(NDIM)
        call SM2GSE_xyz(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGSE(XDIR),xyzGSE(YDIR),xyzGSE(ZDIR))
      END SUBROUTINE SM2GSE_vec

      !GSE=>SM
      SUBROUTINE GSE2SM_xyz (XGSE,YGSE,ZGSE,XSM,YSM,ZSM)
          implicit none
          real(rp), intent(in ) :: XGSE,YGSE,ZGSE
          real(rp), intent(out) :: XSM,YSM,ZSM

          real(rp) :: XGSW,YGSW,ZGSW
          XGSW=XGSE*E11+YGSE*E21+ZGSE*E31
          YGSW=XGSE*E12+YGSE*E22+ZGSE*E32
          ZGSW=XGSE*E13+YGSE*E23+ZGSE*E33

          call GSW2SM(XGSW,YGSW,ZGSW,XSM,YSM,ZSM)
      END SUBROUTINE GSE2SM_xyz

      SUBROUTINE GSE2SM_vec(xyzGSE,xyzSM)
        implicit none
        real(rp), intent(in ) :: xyzGSE(NDIM)
        real(rp), intent(out) :: xyzSM (NDIM)
        
        call GSE2SM_xyz(xyzGSE(XDIR),xyzGSE(YDIR),xyzGSE(ZDIR),xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR))
      END SUBROUTINE GSE2SM_vec

    !--------------------
    !IGRF stuff

      !Takes xyz (SM) and returns IGRF Bxyz (SM) in nT
      FUNCTION IGRF_SM(xyzSM) result(BxyzIGRF_SM)
        real(rp), intent(in ) :: xyzSM(NDIM)
        real(rp) :: BxyzIGRF_SM(NDIM)

        real(rp) :: XGSW,YGSW,ZGSW,HXGSW,HYGSW,HZGSW
        call SM2GSW(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),XGSW,YGSW,ZGSW)
        call IGRF_GSW_08(XGSW,YGSW,ZGSW,HXGSW,HYGSW,HZGSW)
        call GSW2SM(HXGSW,HYGSW,HZGSW,BxyzIGRF_SM(XDIR),BxyzIGRF_SM(YDIR),BxyzIGRF_SM(ZDIR))
        
      END FUNCTION IGRF_SM
      
      SUBROUTINE IGRF_GSW_08 (XGSWi,YGSWi,ZGSWi,HXGSW,HYGSW,HZGSW)
          !Using XGSWi to maintain intent(in) and call geopack routines
          implicit none
          real(rp), intent(in ) :: XGSWi,YGSWi,ZGSWi
          real(rp), intent(out) :: HXGSW,HYGSW,HZGSW

!
!  CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE GEOCENTRIC SOLAR-WIND
!  (GSW) COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL COEFFICIENTS
!  (e.g., https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt, revised 01 January, 2020)
!
!  THE GSW SYSTEM IS ESSENTIALLY SIMILAR TO THE STANDARD GSM (THE TWO SYSTEMS BECOME IDENTICAL
!  TO EACH OTHER IN THE CASE OF STRICTLY ANTI-SUNWARD SOLAR WIND FLOW). FOR A DETAILED
!  DEFINITION, SEE INTRODUCTORY COMMENTS FOR THE SUBROUTINE GSWGSE_08 .
!
!  BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR, IF THE DATE/TIME (IYEAR,IDAY,IHOUR,MIN,ISEC),
!  OR THE SOLAR WIND VELOCITY COMPONENTS (VGSEX,VGSEY,VGSEZ) HAVE CHANGED, THE MODEL COEFFICIENTS
!  AND GEO-GSW ROTATION MATRIX ELEMENTS SHOULD BE UPDATED BY CALLING THE SUBROUTINE RECALC_08.
!
!-----INPUT PARAMETERS:
!
!     XGSW,YGSW,ZGSW - CARTESIAN GEOCENTRIC SOLAR-WIND COORDINATES (IN UNITS RE=6371.2 KM)
!
!-----OUTPUT PARAMETERS:
!
!     HXGSW,HYGSW,HZGSW - CARTESIAN GEOCENTRIC SOLAR-WIND COMPONENTS OF THE MAIN GEOMAGNETIC
!                           FIELD IN NANOTESLA
!
!     LAST MODIFICATION:  FEB 07, 2008.
!     THIS VERSION OF THE CODE ACCEPTS DATES FROM 1965 THROUGH 2025.
!
!     AUTHOR: N. A. TSYGANENKO
!
!
      
        real(rp), dimension(14) :: A,B

        real(rp) :: XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW
        real(rp) :: RHO,RHO2,R,S,CF,SF,PP,P,IRP3
        real(rp) :: C,D,BBR,BBT,BBF,X,Y,Q,Z,W
        real(rp) :: BI,P2,D2,AN,E,HH,QQ,XK,DP,PM,BR,BT,BF
        real(rp) :: HE,HXGEO,HYGEO,HZGEO
        integer  :: N,NM,K,M,MM,MN
        XGSW = XGSWi
        YGSW = YGSWi
        ZGSW = ZGSWi

        CALL GEOGSW_08 (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW,-1)
        RHO2=XGEO**2+YGEO**2
        R=SQRT(RHO2+ZGEO**2)
        C=ZGEO/R
        RHO=SQRT(RHO2)
        S=RHO/R
        IF (S.LT.1.E-5) THEN
          CF=1.
          SF=0.
        ELSE
          CF=XGEO/RHO
          SF=YGEO/RHO
        ENDIF

        PP=1./R
        P=PP
!
!  IN THIS VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
!    HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
!      ON THE VALUE OF THE RADIAL DISTANCE R:
!
        IRP3=R+2
        NM=3+30/IRP3
        IF (NM.GT.13) NM=13

        K=NM+1
        DO N=1,K
           P=P*PP
           A(N)=P
          B(N)=P*N
        ENDDO

        P=1.
        D=0.
        BBR=0.
        BBT=0.
        BBF=0.

      DO M=1,K
         IF(M.EQ.1) THEN
            X=0.
            Y=1.
         ELSE
            MM=M-1
            W=X
            X=W*CF+Y*SF
            Y=Y*CF-W*SF
         ENDIF
         Q=P
         Z=D
         BI=0.
         P2=0.
         D2=0.
         DO N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X                                                                                  
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.NE.1) THEN ! GOTO 180
               QQ=Q
               IF(S.LT.1.E-5) QQ=Z
               BI=BI+AN*(E*X-HH*Y)*QQ
            ENDIF
            XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
            Q=PM
         ENDDO
         D=S*D+C*P
         P=S*P
         IF(M.NE.1) THEN !GOTO 200
            BI=BI*MM
            BBF=BBF+BI
         ENDIF
      ENDDO

      BR=BBR
      BT=BBT
      IF(S.GE.1.E-5) THEN !GOTO 210
         BF=BBF/S
      ELSE
         IF(C.LT.0.) BBF=-BBF
         BF=BBF
      ENDIF

      HE=BR*S+BT*C
      HXGEO=HE*CF-BF*SF
      HYGEO=HE*SF+BF*CF
      HZGEO=BR*C-BT*S

      CALL GEOGSW_08 (HXGEO,HYGEO,HZGEO,HXGSW,HYGSW,HZGSW,1)

      RETURN
    END SUBROUTINE IGRF_GSW_08 
end module geopack
