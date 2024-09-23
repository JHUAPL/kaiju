!
    MODULE Rcm_mod_subs
    use kdefs, ONLY : PI,Mp_cgs,Me_cgs,EarthM0g,eCharge,kev2erg
    use conversion_module, ONLY : almdel
    use rice_housekeeping_module, ONLY: use_plasmasphere, EWMTauInput
    use constants, ONLY: nt, radius_earth_m
    use rcmdefs
    use rcm_precision
    use clocks
    use math

    IMPLICIT NONE
    SAVE
!
!
    INTEGER, PARAMETER :: LUN = 11 !, LUN_2 = 12, LUN_3 = 13

!      Define a number of universal useful constants and parameters:
!      Part 1 is machine-dependent parameters and they should not be changed
!      under any circumstances.
!      Part 2 is physical constants; these may require editing, in which case
!      all the code must be recompiled.
!
    REAL (RPREC), PARAMETER ::     &
!                 Part 1: machine-specific and mathematical parameters
                               zero         = 0.0_rprec, &
                               one          = 1.0_rprec, &
                               two          = 2.0_rprec, &
                               three        = 3.0_rprec, &
                               four         = 4.0_rprec, &
                               five         = 5.0_rprec, &
                               six          = 6.0_rprec, &
                               eight        = 8.0_rprec, &
                               half         = 0.5_rprec, &
                               qtr          = 0.25_rprec,&
                               machine_eps1 = EPSILON (1.0_rprec), &
                               machine_eps2 = machine_eps1*10_rprec, &
                               ! machine_tiny = TINY (one),&
                               ! machine_huge = HUGE (one),&
                               pi_two       = two * pi, &
                               pi_by_two    = pi / two, &
                               rtd          = 180.0_rprec/pi, &
                               dtr          = pi/180.0_rprec, &
                               rth          = 12.0_rprec / pi,&
                               htr          = one / rth      ,&
!
!                 Part 2: physical constants          ! EDITING ALLOWED HERE
                               xmass(RCMNUMFLAV) = [Me_cgs*1.0e-3,Mp_cgs*1.0e-3], &
                               !xmass (2)    = (/ 9.1E-31_rprec, &
                               !                  1.67E-27_rprec /), &
                               !besu         = 3.0584E+4_rprec, &
                               besu         = EarthM0g*1.0e+5, & !Use consistent moment, G => nT
                               signbe       = one, &
                               romeca       = zero, &
                               !charge_e     = 1.6E-19_rprec, &
                               charge_e     = eCharge, & !Take from kdefs
                               sgn (ksize)  = one, &
!                 Part 3: conversion constants 
                               ev2erg       = kev2erg*1.0e-3, & ! conversion from eV to erg
                               m2cm         = 100., & ! conversion from meter to centimeter
                               nT2T         = nt, &  ! conversion from nT to T  
                               dfactor      = nt/radius_earth_m ! conversion for density
                  INTEGER (iprec) :: ie_el = 1, ie_hd = 2 ! coding for e and proton
!
!
!   Potential solver GMRESM tolerance:
    REAL (rprec) :: tol_gmres
    logical :: doRCMVerbose = .FALSE.
    logical :: doRCMVerboseH5 = .FALSE.
!
!
!   This is a definition of the label structure, for I/O:
    TYPE :: label_def
       INTEGER (iprec)   :: intg (20)
       REAL (rprec)      :: real (20)
       CHARACTER(LEN=80) :: char
    END TYPE label_def
    TYPE (label_def) :: label
!
!
!   Define an ellipse:
    TYPE :: ellipse_def
       REAL(rprec) :: aa, bb, xx, yy
    END TYPE ellipse_def
!
    TYPE (ellipse_def) :: boundary (2)
!
!
!   Grid info:
    REAL (rprec) :: dlam, dpsi, Ri, Re, &
                    alpha (isize, jsize), &
                    beta  (isize, jsize), &
                    colat (isize, jsize), &
                    aloct (isize, jsize), &
                    bir   (isize, jsize), &
                    sini  (isize, jsize), &
                    vcorot(isize, jsize), &
                    fac   (isize, jsize)
    INTEGER (iprec) :: i1, i2, iint, j1, j2, jint, imin, imin_j(jsize), ibnd_type
!
!
    LOGICAL ::  L_move_plasma_grid = .TRUE.
    LOGICAL ::  L_doOMPClaw        = .TRUE.
    LOGICAL ::  L_doOMPprecip      = .FALSE.
    LOGICAL ::  doVAvgInit         = .TRUE. !Whether we need to initialize v_avg

    LOGICAL ::  doNoBndFlow = .FALSE.
    integer :: nNBFL = 2
!
!
!   Plasma on grid:
    REAL (rprec) :: alamc (kcsize), etac (kcsize), fudgec (kcsize), &
                    eeta (isize,jsize,kcsize), eeta_cutoff, cmax, &
                    eeta_avg (isize,jsize,kcsize), deleeta(isize,jsize,kcsize), &
                    lossratep(isize,jsize,kcsize), lossmodel(isize,jsize,kcsize), Dpp(isize,jsize), &
                    last_veff(isize,jsize,ksize), &
                    last_ocbDist(isize,jsize)

    INTEGER (iprec) :: ikflavc (kcsize), i_advect, i_eta_bc, i_birk
    LOGICAL :: L_dktime
    INTEGER (iprec), PARAMETER :: irdk=18, inrgdk=13, isodk=2, iondk=2
    REAL (rprec) :: dktime (irdk, inrgdk, isodk, iondk), sunspot_number
    REAL (rprec) :: dtAvg_v,dtMHD
    LOGICAL :: advChannel(kcsize) = .true. !Which channels to advance

    logical :: kill_fudge
!
!
!   Magnetic field:
    REAL (rprec) :: xmin (isize,jsize), ymin (isize,jsize), zmin (isize,jsize), &
                    bmin (isize,jsize), vm (isize,jsize), &
                    rmin (isize,jsize), pmin(isize,jsize),&
                    x1 (isize,jsize), x2 (isize,jsize), &
                    y1 (isize,jsize), y2 (isize,jsize), &
                    b1 (isize,jsize), b2 (isize,jsize), &
                    vm1(isize,jsize), vm2(isize,jsize), &
                    bndloc (jsize),radcurv(isize,jsize),losscone(isize,jsize)
    INTEGER (iprec), ALLOCATABLE :: ibtime (:)
    REAL    (rprec) :: fstoff, fclps, fdst, fmeb, ftilt
    INTEGER (iprec) :: itype_bf
         ! itype_bf = 1 -- read time sequence from input files and interpolate in time
         ! itype_bf = 2 -- read single B-field configuration from input files           
         ! itype_bf = 3 -- expect that an external program will assign B-field arrays
         !                 including pmin and rmin, so do nothing.
!
!
!   Ionospheric quantities:
    REAL (rprec) :: qtped (isize,jsize), pedpsi (isize,jsize), &
                    qtplam(isize,jsize), pedlam (isize,jsize), &
                    qthall(isize,jsize), hall   (isize,jsize), &
                    ss (jsize), &
                    pwe (isize,jsize), pwn (isize,jsize), &
                    hwe (isize,jsize), hwn (isize,jsize), &
                    sw  (jsize), &
                    eflux (isize,jsize,iesize), eavg (isize,jsize,iesize), nflux (isize,jsize,iesize), &
                    efluxk (isize,jsize,kcsize,iesize), eavgk (isize,jsize,kcsize,iesize), nfluxk (isize,jsize,kcsize,iesize)
    INTEGER (iprec) :: icond, nsmthi, nsmthj, iwind
    LOGICAL :: ifloor, icorrect
!
!
!   Magnetospheric quantities:
    REAL (rprec) :: v (isize,jsize), vpar (isize,jsize), vbnd (jsize), &
                    birk (isize,jsize), pvgamma (isize,jsize,iesize), &
                    pressrcm (isize,jsize), &
                    v_avg (isize,jsize), birk_avg (isize,jsize), &
                    densrcm(isize,jsize),denspsph(isize,jsize)
    INTEGER (iprec) :: ipcp_type, ipot
!
!   Input PCP drop and its current value:
    INTEGER (iprec), ALLOCATABLE :: ivtime (:)
    REAL    (rprec), ALLOCATABLE :: vinput (:), vinput_phase(:)
    REAL    (rprec)              :: vdrop,      vdrop_phase
!
!
 
    INCLUDE 'rcmdir.h'

!    Logical :: IsCoupledExternally = .false.  ! flag to determine if RCM is standalone or not

    ! Variables for internal RCM timing:
    INTEGER(iprec) :: timer_start(10) = 0, timer_stop(10) = 0, count_rate
    REAL (rprec) :: timer_values (10)=0.0_rprec

!
    INTERFACE Gntrp
       MODULE PROCEDURE Gntrp_2d_ang
    END INTERFACE
 
    INTERFACE Interp
      MODULE PROCEDURE Interp_1d, Interp_2d, Interp_2d_of3d
    END INTERFACE

    INTERFACE Bjmod
       MODULE PROCEDURE Bjmod_int, Bjmod_real
    END INTERFACE

    INTERFACE Circle
       MODULE PROCEDURE Circle_1d, Circle_2d
    END INTERFACE
!
!
!
    CONTAINS
!
!
!
!      SUBROUTINE Comput (jtime, dt )
      SUBROUTINE Comput (dtCpl)
      IMPLICIT NONE
!      INTEGER (iprec), INTENT (IN) :: jtime
!      REAL (rprec),    INTENT (IN) :: dt
!
      REAL (rprec), INTENT(IN)  :: dtCpl
      INTEGER (iprec) :: j
      REAL (rprec)  ::  a(3), b(3), dx(3), dy(3), deqdt
!
!

    CALL Tic("GET_JBIRK")
    CALL Get_jbirk
    if (doRCMVerbose) then
      write(6,*)'RCM: finish getting jbirk' 
    endif
    CALL Toc("GET_JBIRK")

    CALL Tic("PRECIP")
    CALL kdiffPrecip(dtCpl)
    !CALL diffusePrecip (dtCpl)
    !Call diffusePrecipMaxwellian ()
    if (doRCMVerbose) then
      write(6,*)'RCM: finish getting diffuse precipitation'
    endif
    CALL Toc("PRECIP")

!      IF (ibnd_type == 4) THEN
!        DO NOTHING, VBND IS ALREADY SET
!         DO j = 1, jsize
!            v (1:imin_j(j)-1,j) = vbnd (j)
!         END DO
!      ELSE
!         STOP 'COMPUT: ibnd_type not implemented'
!      END IF
!
!
      RETURN 
      END SUBROUTINE Comput                         
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_jbirk ( )
      IMPLICIT NONE
!__________________________________________________________________________
!                                                                       
!  Program written by: r.w. spiro        
!  last update:
!     04-05-88          
!     01-29-96 frt                 - added ain,min_j arr
!  Algorithm by: r.a. wolf                                              
!                                                                       
!  This subroutine computes birk(i,j) given inner edge
!  locations 
!  modified 04-05-88 to include effects of gradients in eta.            
!  see raw document re including eeta in computation of jbirk           
!    dated feb 6, 1988.                                                 
!  birk is current density (2 hemispheres) in units of
!  microamp/m**2    
!
!  birk(i,j) here is the field-aligned current density per
!  unit of ionospheric area, so that it already includes
!  the factor sin(I); this is J_parallel*sin(I) in the RHS
!  of the Vasyliunas equation.
!
!  Issues with non-integer boundary (Stanislav's notes):
!  for BIRK from inner edge segments, this is not an issue
!  (except that if a segment is entirely outside the bndry,
!  then we don't compute its contribution); of course, we 
!  have to care about this somewhere else where motion of 
!  test particles is computed. For BIRK from gradients of 
!  EETA,  
!  removed edges 3/19 frt
!                                                                       
!______________________________________________________________________________
!
      REAL (rprec), PARAMETER :: cf1 = one / pi_two, &
                                 cf2 =  - (three/four)*( (two / pi) - half)
!
      INTEGER (iprec) :: i, j, k, kc, klbeg, klend, kl, klnext, &
                 ibmin, ibmax, jbmin, jbmax, jb1, jb2,  &
                 ig, jj, jindex, ib1, ib2
      REAL (rprec) :: detadi(isize,jsize), detadj(isize,jsize), &
                      dvmdi(isize,jsize), dvmdj(isize,jsize), dbirk, &
                      vmkl, vmnext, sum, b1, b2, x, y, el, umax, umin, ss, &
                      z, dg1, dg2, dg3, qmin, qmax, qn, qx,  &
                      denom, a1, a2, bjm, range, bim, gkl (5000)
      LOGICAL, dimension(1:isize,1:jsize) :: isOpen

!                                                                       
!
      birk (:,:) = zero
!
!
!     Compute J_parallel due to continuous channel:
!             
      !Replacing gradient w/ slope-limited gradient used in advection       
      isOpen = (vm < 0)
      call Grad_IJ(vm,isOpen,dvmdi,dvmdj)

      do kc=1,kcsize
        !Using new gradient, Grad_IJ is internally threaded so don't call it from OMP loop
        call Grad_IJ(eeta(:,:,kc),isOpen,detadi,detadj)

        !NOTE: Not great to have OMP inside k loop but easier than writing an unthreaded Grad_IJ
        !$OMP PARALLEL DO &
        !$OMP schedule(dynamic) &
        !$OMP DEFAULT (NONE) &
        !$OMP PRIVATE(i,j,dbirk) &
        !$OMP SHARED(kc,j1,j2,i2,alamc,dlam,dpsi,Ri) &
        !$OMP SHARED(alpha,beta,detadi,detadj,dvmdj,dvmdi,eeta,birk,isOpen)
        DO  j = j1, j2
          !Calculate Vasyliunas FAC wherever possible, even using MHD buffer cells  
          DO i = 1,i2
            if (isOpen(i,j)) CYCLE

            dbirk  = charge_e * signbe * ABS(alamc(kc)) * &
                     (detadj(i,j) * dvmdi(i,j) - detadi(i,j)*dvmdj(i,j)) / &
                     (alpha(i,j)*beta(i,j)*dlam*dpsi*Ri**2)
            birk (i, j) = birk (i, j) + dbirk    
          ENDDO !i
        ENDDO !j

      enddo

      CALL Circle (birk)

      END SUBROUTINE Get_jbirk
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_jbirk2 ( )
      IMPLICIT NONE
!__________________________________________________________________________
!                                                                       
!  Program written by: r.w. spiro        
!  last update:
!     04-05-88          
!     01-29-96 frt                 - added ain,min_j arr
!     11-14-02 frt uses a new version on the grid-based scheme
!  Algorithm by: r.a. wolf                                              
!                                                                       
!  This subroutine computes birk(i,j) given inner edge
!  locations 
!  modified 04-05-88 to include effects of gradients in eta.            
!  see raw document re including eeta in computation of jbirk           
!    dated feb 6, 1988.                                                 
!  birk is current density (2 hemispheres) in units of
!  microamp/m**2    
!
!  birk(i,j) here is the field-aligned current density per
!  unit of ionospheric area, so that it already includes
!  the factor sin(I); this is J_parallel*sin(I) in the RHS
!  of the Vasyliunas equation.
!
!  Issues with non-integer boundary (Stanislav's notes):
!  for BIRK from inner edge segments, this is not an issue
!  (except that if a segment is entirely outside the bndry,
!  then we don't compute its contribution); of course, we 
!  have to care about this somewhere else where motion of 
!  test particles is computed. For BIRK from gradients of 
!  EETA,  
! 
!  This version uses a new way to compute jbirk from the cts channel - frt
!                                                                       
!______________________________________________________________________________
!
      REAL (rprec), PARAMETER :: cf1 = one / pi_two, &
                                 cf2 =  - (three/four)*( (two / pi) - half)
!
      INTEGER (iprec) :: i, j, k, kc, klbeg, klend, kl, klnext, &
                 ibmin, ibmax, jbmin, jbmax, jb1, jb2,  &
                 ig, jj, jindex, ib1, ib2
      REAL (rprec) :: detadi(isize,jsize), detadj(isize,jsize), &
                      dvmdi(isize,jsize), dvmdj(isize,jsize), dbirk, &
                      vmkl, vmnext, sum, b1, b2, x, y, el, umax, umin, ss, &
                      z, dg1, dg2, dg3, qmin, qmax, qn, qx,  &
                      denom, a1, a2, bjm, range, bim, gkl (5000)
      REAL (rprec) :: eeta2(isize,jsize,ksize),vm2(isize,jsize)
!                                                                       
!
      birk (:,:) = zero
!
!
!     Compute J_parallel due to continuous channel:
!                                                                       
! define new temporary work arrays eeta2 and vm2
      eeta2 = eeta
      vm2   = vm
      DO j=1,jsize
       DO i=1,imin_j(j)-1
        DO k=1,kcsize
         eeta2(i,j,k) = eeta(imin_j(j),j,k)
        END DO
         vm2(i,j) = vm(imin_j(j),j)
       END DO
      END DO

!    1          2            3
! i-1,j+1------i,j+1------i+1,j+1
!    |          |            |  
!    |          |            |  
!    |          |            |  
!    4          |            5
! i-1,j--------i,j--------i+1,j
!    |          |            |  
!    |          |            |  
!    |          |            |  
!    6          7            8
! i-1,j-1------i,j-1------+1,j-1 
! the basic equation is
! dbirk  is proportional to 
! (eta_1*(vm_2-vm_1)+eta_2*(vm_3-vm_1)+eta_3*(vm_5-vm_2)+eta_4*(vm_1-vm_6)
!  eta_5*(vm_8-vm_3)+eta_6*(vm_4-vm_7)+eta_7*(vm_6-vm_8)+eta_8*(vm_7-vm_5))/8
      DO kc = 1, kcsize
!
         DO  j = j1, j2 
           DO  i = imin_j(j)+1, i2

            dbirk =(eeta2(i-1,j+1,kc)*(vm2(i  ,j+1)-vm2(i-1,j  )) + &
                    eeta2(i  ,j+1,kc)*(vm2(i+1,j+1)-vm2(i-1,j+1)) + &
                    eeta2(i+1,j+1,kc)*(vm2(i+1,j  )-vm2(i  ,j+1)) + &
                    eeta2(i-1,j  ,kc)*(vm2(i-1,j+1)-vm2(i-1,j-1)) + &
                    eeta2(i+1,j  ,kc)*(vm2(i+1,j-1)-vm2(i+1,j+1)) + &
                    eeta2(i-1,j-1,kc)*(vm2(i-1,j  )-vm2(i  ,j-1)) + &
                    eeta2(i  ,j-1,kc)*(vm2(i-1,j-1)-vm2(i+1,j-1)) + &
                    eeta2(i+1,j-1,kc)*(vm2(i  ,j-1)-vm2(i+1,j  ))) /8.

               dbirk  = charge_e * signbe * ABS(alamc(kc)) * dbirk/ &
                        (alpha(i,j)*beta(i,j)*dlam*dpsi*Ri**2)
               birk (i, j) = birk (i, j) + dbirk 
           END DO 
         END DO 
      END DO 
!                                                                       
!print*,'zeroing birk'
      CALL Circle (birk)
!
      RETURN 
      END SUBROUTINE Get_jbirk2
!
!
!*************************************************************************
!
!
      SUBROUTINE diffusePrecip (dtCpl)
      IMPLICIT NONE

!--------------------------------------------------------------------------
! sbao 05/2021
! This subroutine calculates diffuse electron precipitation using deleeta
! The equation of the differential flux is adapted from M. Gkioulidou et al (doi:10.1029/2012JA018032)
      REAL (rprec), INTENT(IN)  :: dtCpl
      INTEGER (iprec) :: i, j, ie, iedim_local, kc, klow
      REAL (rprec)    :: en, delEn, Jk, sum1 (iesize), sum2 (iesize)
      LOGICAL, dimension(1:isize,1:jsize) :: isOpen
      REAL (rprec) :: JkConst 


      !Try to do calculation everywhere possible including MHD buffer region
      isOpen = (vm < 0)
 
      !Set lowest RC channel
      if (use_plasmasphere) then
         klow = 2
      else
         klow = 1
      endif

      iedim_local = 2 ! # of species, electron and proton
!
      eavg  (:,:,:) = zero
      eflux (:,:,:) = zero

      loop_j: DO j = j1, j2
      !loop_i: DO i = imin_j(j), isize
      loop_i: DO i = 1, isize
            if (isOpen(i,j)) CYCLE
!           Now for each grid point, consider all species
!           present at that grid point, and compute sum1 and
!           sum2 for positive and negative particles separately:

!           For each grid point, clear sum1 and sum2:
!
            sum1 (1:iedim_local) = zero
            sum2 (1:iedim_local) = zero
!
            GRID_BASED: DO kc = klow, kcsize
              IF (alamc (kc) < zero) THEN
                 ie = 1  ! electron
              ELSE
                 ie = 2  ! proton 
              END IF
              en = ABS(alamc(kc))*vm(i,j) ! channel energy in eV
              delEn = ABS(almdel(kc))*vm(i,j) ! channel width in eV 
              JkConst = 1./(SQRT(8.*xmass(ie))*pi)*SQRT(charge_e)*nt/m2cm**2/radius_earth_m ! Constant for Jk
              !Jk = JkConst*SQRT(ABS(alamc(kc)))* deleeta(i,j,kc)/dtCpl*vm(i,j)/almdel(kc)   ! differential energy flux in 1/(eV cm^2 s sr)
              Jk = JkConst*SQRT(ABS(alamc(kc)))* 1./3.*eeta(i,j,kc)/dtCpl*vm(i,j)/almdel(kc)
              sum1(ie) = sum1(ie) + en*Jk*delEn !  in eV/(cm^2 s sr)
              sum2(ie) = sum2(ie) + Jk*delEn ! in 1/(cm^2 s sr)
            END DO GRID_BASED
              
            DO ie = 1, iedim_local
!                                                                       
               IF (sum2 (ie) > 10.*machine_tiny) THEN  ! zero  sbao 07/2019
!
!                compute thermal electron current, field-aligned
!                potential drop, electron energy flux,
!                and average electron energy at (i,j):          
!
                  eflux(i,j,ie) = ev2erg*pi*sum1(ie) ! energy flux in erg/(cm^2 s), pi comes from the vel. space integral 
                  eavg(i,j,ie) = sum1(ie)/sum2(ie)  ! averge energy in eV
                 
               ELSE
!                 we want eflux=0 and eavg=0 for no precipitation.
                  eflux (i, j, ie) = zero
                  eavg  (i, j, ie) = zero
!
               END IF
           
            END DO

      END DO loop_i
      END DO loop_j

      CALL Circle (eflux (:, :, ie_el))
      CALL Circle (eavg  (:, :, ie_el))
      CALL Circle (eflux (:, :, ie_hd))
      CALL Circle (eavg  (:, :, ie_hd))

      END SUBROUTINE diffusePrecip

      ! K: A brute force diffuse precipitation.
      ! Particles lost through scattering should precipitate
      subroutine kdiffPrecip(dtCpl)
        IMPLICIT NONE
        REAL (rprec), INTENT(IN)  :: dtCpl
        LOGICAL, dimension(1:isize,1:jsize) :: isOpen
        real(rprec), dimension(RCMNUMFLAV) :: nflx,eflx
        integer(iprec) :: klow,i,j,k,ie
        real(rprec) :: dn
        real(rprec) :: kin
        !Try to do calculation everywhere possible including MHD buffer region
        isOpen = (vm < 0)
        !Set lowest RC channel
        if (use_plasmasphere) then
            klow = 2
        else
            klow = 1
        endif
        eavg  (:,:,:) = 0.0
        eflux (:,:,:) = 0.0
        nflux (:,:,:) = 0.0
        do j=1,jsize
            do i=1,isize
                if (isOpen(i,j)) CYCLE
                nflx = 0.0
                eflx = 0.0
                do k=klow,kcsize
                    ! only accumulate diff for keV below 30 keV.
                    kin = abs(alamc(k)*vm(i,j))*1.0e-3 !Energy [keV]
                    IF (alamc (k) < -TINY .and. kin<=30.0) THEN
                        ie = RCMELECTRON
                    else if (alamc (k) > +TINY) then
                        ie = RCMPROTON
                    else
                        cycle
                    endif
                    !Now accumulate, 0.5 for single hemisphere, see detailed deriviation of the diffuse precipitation on kaiju wiki    
                    dn = 0.5*deleeta(i,j,k)*abs(bir(i,j)/sini(i,j))*nT2T/(m2cm**2)/dtCpl ! #/cm^2/s
                    nflx(ie) = nflx(ie) + dn !Num flux, #/cm2/s
                    eflx(ie) = eflx(ie) + dn*ABS(alamc(k))*vm(i,j) !Energy flux, eV/cm^2/s
                enddo
                eflux(i,j,:) = eflx*ev2erg  ! energy flux in erg/(cm^2 s)
                eavg (i,j,:) = eflx/nflx    ! Average energy in eV
                nflux(i,j,:) = nflx         ! Num flux in #/cm^2/s
               
                DO ie = 1, RCMNUMFLAV
                      IF (nflx(ie) > TINY) THEN 
                          ! The ratio of eavg is only meaningful when nflx is meaningful.
                          ! Note in REMIX, a higher floor will be applied toward the final precipitation.
                          eavg (i,j,ie) = eflx(ie)/nflx(ie) ! Average energy in eV
                      ELSE
!                         we want eflux=0 and eavg=0 for no precipitation.
                          eflux (i, j, ie) = zero
                          eavg  (i, j, ie) = zero
                          nflux (i, j, ie) = zero
!
                      END IF
                END DO

            enddo
        enddo
        
        CALL Circle (nflux (:, :, ie_el))
        CALL Circle (eflux (:, :, ie_el))
        CALL Circle (eavg  (:, :, ie_el))
        CALL Circle (nflux (:, :, ie_hd))
        CALL Circle (eflux (:, :, ie_hd))
        CALL Circle (eavg  (:, :, ie_hd))

        where(eflux<0.01 .or. eavg<0.01)
        !Do both or neither
          eavg  = 0.0
          eflux = 0.0
          nflux = 0.0
        end where
      end subroutine kdiffPrecip

 
      SUBROUTINE diffusePrecipChannel ()
      IMPLICIT NONE

!--------------------------------------------------------------------------
! sbao 05/2021
! This subroutine calculates diffuse electron precipitation using deleeta for each energy channel
! The equation of the differential flux is adapted M. Gkioulidou et al (doi:10.1029/2012JA018032)


      END SUBROUTINE diffusePrecipChannel


      SUBROUTINE diffusePrecipMaxwellian ()
      IMPLICIT NONE

!--------------------------------------------------------------------------
! sbao 01/2021
! This subroutine calculates diffuse electron precipitation using Maxwellian distribuiton, adapted from Get_vparallel

      INTEGER (iprec) :: i, j, ie, iedim_local, kc
      REAL (rprec)    :: en, ekt, therm, sum1 (iesize), sum2 (iesize)
      LOGICAL, dimension(1:isize,1:jsize) :: isOpen

      !Try to do calculation everywhere possible including MHD buffer region
      isOpen = (vm < 0)

      iedim_local = 2
!
      vpar  (:,:)   = zero
      eavg  (:,:,:) = zero
      eflux (:,:,:) = zero


      !$OMP PARALLEL DO if (L_doOMPprecip) &
      !$OMP schedule(dynamic) &
      !$OMP DEFAULT (NONE) &
      !$OMP PRIVATE(i,j,kc,ie,sum1,sum2) &
      !$OMP PRIVATE(en,ekt,therm) &
      !$OMP SHARED(j1,j2,iedim_local,imin_j,alamc,eeta) &
      !$OMP SHARED(vpar,vm,fudgec,birk,eflux,eavg,isOpen) 

      loop_j: DO j = j1, j2
      !loop_i: DO i = imin_j(j), isize
      loop_i: DO i = 1, isize
            if (isOpen(i,j)) CYCLE
!
!           For each grid point, clear sum1 and sum2:
!
            sum1 (1:iedim_local) = zero
            sum2 (1:iedim_local) = zero
!
!
!           Now for each grid point, consider all species
!           present at that grid point, and compute sum1 and
!           sum2 for positive and negative particles separately:
!
            GRID_BASED: DO kc = 1, kcsize
            ! IF ( ABS(alamc(kc))*vm(i,j) > 500.0_rprec) THEN
              IF (alamc (kc) < zero) THEN
                 ie = 1 
              ELSE 
                 ie = 2 
!              STOP 'BALGN4: ie is 2'
              END IF
              sum1(ie) = sum1(ie) + eeta(i,j,kc)*fudgec(kc)
              sum2(ie) = sum2(ie) + eeta(i,j,kc)*fudgec(kc)*ABS(alamc(kc))
             !END IF
            END DO GRID_BASED 
!
!           For positive and negative particles separately,
!           compute precipitating number flux, average energy,
!           and parallel potential drop:
!
            DO ie = 1, iedim_local 
!                                                                       
               IF (sum1 (ie) > 10.*machine_tiny) THEN  ! zero  sbao 07/2019
!
!                compute thermal electron current, field-aligned
!                potential drop, electron energy flux,
!                and average electron energy at (i,j):          
!
                  en    = sum1 (ie) * vm (i, j)**1.5 / 6.38E+21
                  ekt   = (two/three) * sum2 (ie) * vm (i,j) / sum1 (ie)
                  therm = 0.02675 * en * SQRT(ekt*xmass(1)/xmass(ie))
!
                  IF (therm < 1.E-30) therm = zero

                  eflux(i,j,ie) = 0.002 * therm * ekt 
                  eavg(i,j,ie) = two*ekt
                  ! sbao 6/19 detect Nan 
                  if (ISNAN(eflux(i,j,ie)))then
                       eflux(i,j,ie) = 0.0
                       eavg(i,j,ie) = 0.0
                  end if

               ELSE 
!                                                                       
!                 Case fudge=0: we want eflux=0 and eavg=0 for no precipitation.
!
                  eflux (i, j, ie) = zero
                  eavg  (i, j, ie) = zero
!
               END IF 
               ! corrections to eavg at eflux(i,j) == 0.0     sbao 07/2019 
               ! == does not work well with real number, use lt threshold instead. ldong 04/2020
               if ( (eflux(i,j,ie) .lt. 0.01) .or. (eavg(i,j,ie) .lt. 0.01) ) then
                  !Do both or neither
                  eavg(i,j,ie) = 0.0
                  eflux(i,j,ie) = 0.0
               endif

               ! IF (eflux(i,j,ie) .lt. 0.01) eavg(i,j,ie) = 0.0
               ! IF (eavg(i,j,ie) .lt. 0.01) eflux(i,j,ie) = 0.0
!                                                                       
            END DO
!
      END DO loop_i
      END DO loop_j 
!                                                                       
!
      CALL Circle (eflux (:, :, ie_el))
      CALL Circle (eavg  (:, :, ie_el))
      CALL Circle (eflux (:, :, ie_hd))
      CALL Circle (eavg  (:, :, ie_hd))
!
      RETURN
      END SUBROUTINE diffusePrecipMaxwellian
!
!
!==============================================================================
!
!
      SUBROUTINE Get_vparallel ()
      IMPLICIT NONE
!______________________________________________________________________________
!  last update: 
!     05-05-87       by:rws                              
!     02-10-96          frt - added arrays ain,min_j     
!                                                                       
!  Birk is sum of current densities into both hemispheres.              
!  (micro amp/m**2).  Before activating parallel potential drop        
!  we need to check if birk is being used correctly in
!  this routine.   
!
!  Stanislav: VPAR is computed inside IE loop (for both
!             negative and positive particles), and will
!             be the one for the largest IE value. Which
!             is nonsense.
!  Stanislav: this subroutine needs grid-based formulation
!             of plasma (EETA). Before it was done by
!             computing EETA for electrons from the inner
!             edges of electrons, then it was changed to
!             use directly grid-based population. In the
!             latter case, array PVEC returned by this
!             routine is the electron pressure (without
!             the factor of 2/3) and is the same as what
!             routine PV returns as array PVGAM. If the
!             electrons are on the grid only, as in my case,
!             then we call PV in rcm main program to compute
!             the ion pressure, and we use PVEC from this
!             routine for the electron pressure. (04/20/99)
!  Stanislav, may 18,99: make all loops over electrons only,
!             by using iedim_local and setting it to 1.
!______________________________________________________________________________
!
      INTEGER (iprec) :: i, j, ie, iedim_local, kc
      REAL (rprec)    :: en, ekt, therm, sum1 (iesize), sum2 (iesize)
!                                                                       
!
!                                  
      iedim_local = 2
!
      vpar  (:,:)   = zero
      eavg  (:,:,:) = zero
      eflux (:,:,:) = zero


      loop_j: DO j = j1, j2
      loop_i: DO i = imin_j(j), isize
!
!           For each grid point, clear sum1 and sum2:
!
            sum1 (1:iedim_local) = zero
            sum2 (1:iedim_local) = zero
!
!
!           Now for each grid point, consider all species
!           present at that grid point, and compute sum1 and
!           sum2 for positive and negative particles separately:
!
            GRID_BASED: DO kc = 1, kcsize
             IF ( ABS(alamc(kc))*vm(i,j) > 500.0_rprec) THEN
               IF (alamc (kc) < zero) THEN
                  ie = 1 
               ELSE 
                  ie = 2 
!                 STOP 'BALGN4: ie is 2'
               END IF
               sum1(ie) = sum1(ie) + eeta(i,j,kc)*fudgec(kc)
               sum2(ie) = sum2(ie) + eeta(i,j,kc)*fudgec(kc)*ABS(alamc(kc))
             END IF
            END DO GRID_BASED 
!
!           For positive and negative particles separately,
!           compute precipitating number flux, average energy,
!           and parallel potential drop:
!
            DO ie = 1, iedim_local 
!                                                                       
               IF (sum1 (ie) > 10.*machine_tiny) THEN  ! zero  sbao 07/2019
!
!                compute thermal electron current, field-aligned
!                potential drop, electron energy flux,
!                and average electron energy at (i,j):          
!
                  en    = sum1 (ie) * vm (i, j)**1.5 / 6.38E+21
                  ekt   = (two/three) * sum2 (ie) * vm (i,j) / sum1 (ie)
                  therm = 0.02675 * en * SQRT(ekt*xmass(1)/xmass(ie))
!
                  IF (therm < 1.E-30) THEN 
                     therm      = zero
                     vpar (i,j) = zero
                  ELSE 
                     IF (- birk (i, j) / therm > one) THEN
                        vpar (i,j) = ekt * (- birk (i,j) / therm - one)
                     ELSE 
                        vpar (i,j) = one
                     END IF
                     vpar(i,j) = MIN (vpar (i, j), 10000.0_rprec)
                  END IF 
!
!    !!!!!!!      ALERT: VPAR(I,J) IS SET TO 0 !!!!!!!!!!!!!!!!!
!
                  vpar (i, j) = zero
!
!
                  eflux(i,j,ie) = 0.002 * therm * &
                                 ( ekt + vpar(i,j) + half*vpar(i,j)**2/ekt)
                  eavg(i,j,ie) = two*(ekt+vpar(i,j)+half*vpar(i,j)**2 /ekt) / &
                                 (one + vpar (i, j) / ekt)
                  ! sbao 6/19 detect Nan 
                  if (ISNAN(eflux(i,j,ie)))then
                       if (.not. doQuietRCM) write(*,*)'eflux,i,j,therm,ekt,vpar,sum1,sum2,vm',eflux(i,j,ie),i,j,therm,ekt,vpar(i,j),sum1(ie),sum2(ie),vm(i,j)
                       eflux(i,j,ie) = 0.0
                       eavg(i,j,ie) = 0.0
                  end if

               ELSE 
!                                                                       
!                 Case fudge=0: we want eflux=0 and eavg=0 for no precipitation.
!
                  eflux (i, j, ie) = zero
                  eavg  (i, j, ie) = zero
!
               END IF 
               ! corrections to eavg at eflux(i,j) == 0.0     sbao 07/2019 
               ! == does not work well with real number, use lt threshold instead. ldong 04/2020
               IF (eflux(i,j,ie) .lt. 0.01) eavg(i,j,ie) = 0.0
               IF (eavg(i,j,ie) .lt. 0.01) eflux(i,j,ie) = 0.0

!                                                                       
            END DO
!
      END DO loop_i
      END DO loop_j 
!                                                                       
!
!
      CALL Circle (vpar)
      CALL Circle (eflux (:, :, ie_el))
      CALL Circle (eavg  (:, :, ie_el))
      CALL Circle (eflux (:, :, ie_hd))
      CALL Circle (eavg  (:, :, ie_hd))
!
      RETURN
      END SUBROUTINE Get_vparallel
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Floor_for_eflux ()
      IMPLICIT NONE
      INTEGER (iprec) :: ivalue_max, i, j
      REAL (rprec)    :: eflux_max
      DO j = 1, jsize
         eflux_max = eflux (isize, j, ie_el)
         ivalue_max = isize
         DO i = isize-1, imin_j(j), -1
            IF (eflux(i,j,ie_el) > eflux(i+1,j,ie_el)) THEN
               eflux_max  = eflux(i,j,ie_el)
               ivalue_max = i
            END IF
         END DO
         DO i = imin_j(j), ivalue_max - 1
            eflux(i,j,ie_el) = MAX (half*eflux_max, eflux(i,j,ie_el))
         END DO
      END DO
      !ion precipitation
      DO j = 1, jsize
         eflux_max = eflux (isize, j, ie_hd)
         ivalue_max = isize
         DO i = isize-1, imin_j(j), -1
            IF (eflux(i,j,ie_hd) > eflux(i+1,j,ie_hd)) THEN
               eflux_max  = eflux(i,j,ie_hd)
               ivalue_max = i
            END IF
         END DO
         DO i = imin_j(j), ivalue_max - 1
            eflux(i,j,ie_hd) = MAX (half*eflux_max, eflux(i,j,ie_hd))
         END DO
      END DO
      RETURN
      END SUBROUTINE Floor_for_eflux
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!

    SUBROUTINE Move_plasma_grid (dt, i_start, i_stop, j_start, j_stop,ie_ask)
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: dt
    INTEGER (iprec), INTENT (IN) :: ie_ask, i_start, i_stop, j_start, j_stop
!_____________________________________________________________________________
!   Subroutine to advance eta distribution for a time step
!   by a lifetime-based algorithm (raw doc dated 5/12/87)
!                                                                       
!   Last update: 05-11-88
!                01-29-96 ain ,min_j and calls to bndy added - frt
!   rws          06-05-97 etamov changed to reflect new use of
!                        eeta array in rcm697 version
!
!   CALLED FROM: TSTEP1
!_____________________________________________________________________________
!
!
    REAL (rprec) :: eeta2 (isize,jsize), veff  (isize,jsize), &
                    dvefdi(isize,jsize), dvefdj(isize,jsize), &
                    didt, djdt, biold, bjold, rate, mass_factor, a1, a2, fi, fj
    INTEGER (iprec) :: i, j, kc, ie, i_1, i_2, j_1, j_2
    LOGICAL :: pt_1_1, pt_1_2, pt_2_1, pt_2_2
!
	real :: v_1_1, v_1_2, v_2_1, v_2_2
!
    DO kc = 1, kcsize
!
       IF (alamc(kc) < zero) THEN
          ie = 1
       ELSE
          ie = 2 ! but must change if O+ is added
       END IF
       IF (ie /= ie_ask) CYCLE
       mass_factor = SQRT (xmass(1) / xmass(ie))
       veff = v + vcorot - vpar  + alamc(kc)*vm

       last_veff(:,:,kc) = veff  ! Storing for output
!
       dvefdi = Deriv_i (veff, imin_j)
       dvefdj = Deriv_j (veff, imin_j, j1, j2, 1.0E+25_rprec)
       WHERE (ABS(dvefdj) > 1.0E+24)
          dvefdj = 0.0
          dvefdi = 0.0
       END WHERE
!
!
       eeta2 (:,:) = eeta (:,:,kc)
       DO j = j1, j2
!      DO i = imin_j(j), i2
       DO i = i_start, i_stop
          IF (i <= imin_j(j)) CYCLE
          didt   =   dvefdj (i,j) / fac (i,j)
          djdt   = - dvefdi (i,j) / fac (i,j)
          biold  = REAL(i,rprec) - didt * dt 
          bjold  = Bjmod (REAL(j,rprec) - djdt * dt, jwrap, jsize )
          rate   = RatefnFDG (fudgec(kc), alamc(kc), sini (i,j), bir (i,j), &
                           vm (i,j), mass_factor)
!          IF (biold > Bndy(bndloc,bjold)) THEN 
          IF (biold > bndloc(j)) THEN 
!            Particle came from within the modeling region, find ETA_old by interp:
             i_1 = INT (biold)
             i_2 = i_1 + 1
             j_1 = INT (bjold)
             j_2 = j_1 + 1
!
!   (i_1,j_1) x--------------x (i_1,j_2)
!             |              |
!             |              |
!             |              |
!             |              |
!   (i_2,j_1) x--------------x (i_2,j_2)
!
             pt_1_1 = i_1 >= imin_j(j_1)
             pt_1_2 = i_1 >= imin_j(j_2)
             pt_2_1 = i_2 >= imin_j(j_1)
             pt_2_2 = i_2 >= imin_j(j_2)
!
             v_1_1 = eeta2 (i_1,j_1)
             v_1_2 = eeta2 (i_1,j_2)
             v_2_1 = eeta2 (i_2,j_1)
             v_2_2 = eeta2 (i_2,j_2)
!
             IF ((.NOT.pt_2_1) .OR. (.NOT. pt_2_2)) THEN
                 STOP 'ONE OF I_2 POINTS OUT OF BNDY'
             ELSE IF ((.NOT.pt_1_1) .AND. (.NOT.pt_1_2)) THEN
                 STOP 'BOTH I_1 POINTS OUT OF BNDY'
             ELSE IF (pt_1_1) THEN ! get 1,2 by interp.
                 
             ELSE ! get 1,1 by interp.
                 
             END IF
!            IF (i_1 < imin_j(j_1) .OR. i_2 < imin_j(j_1) .OR. &
!                i_1 < imin_j(j_2) .OR. i_2 < imin_j(j_2)) THEN
!              STOP 'PNT IS OUTSIDE MODELING REGION'
!            END IF
             fi = REAL(i_1)
             fj = REAL(j_1)
             a1 = (1.0-(biold-fi))*eeta2(i_1,j_1) + (biold-fi)*eeta2(i_2,j_1)
             a2 = (1.0-(biold-fi))*eeta2(i_1,j_2) + (biold-fi)*eeta2(i_2,j_2)
             eeta (i,j,kc) = ((1.0 - (bjold-fj)) * a1 + (bjold-fj) * a2)*EXP(-rate*dt)
          ELSE
!            Particle came from outside the modeling region, find ETA_old from b.c.:
             STOP 'PARTICLE IS FLOWING IN, NOT IMPLEMENTED'
          END IF
       END DO
       END DO
!
       CALL Circle (eeta(:,:,kc))
!
    END DO
!      
       ! refill the plasmasphere  04012020 sbao       
       CALL Plasmasphere_Refilling_Model(eeta(:,:,1), rmin, aloct, vm, dt)
       CALL Circle (eeta(:,:,1))

    RETURN

    END SUBROUTINE Move_plasma_grid
!
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Circle_2d (r)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN OUT) :: r (:,:)
!
      INTEGER (iprec) :: jlast, i, j, imax, jmax
!
      imax = SIZE (r, DIM = 1)
      jmax = SIZE (r, DIM = 2)
      jlast = jmax - jwrap 
      
      DO i = 1, imax 
        DO  j = 1, jwrap - 1
          r (i, j) = r (i, jlast + j)
        END DO
        r (i, jmax) = r (i, jwrap) 
      END DO
      RETURN 
      END SUBROUTINE Circle_2d
!
!
!
!
!
      SUBROUTINE Circle_1d (r)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN OUT) :: r (:)
!
      INTEGER (iprec) :: jlast, j, jmax
!
      jmax = SIZE (r, DIM = 1)
      jlast = jmax - jwrap

      DO  j = 1, jwrap - 1
        r (j) = r (jlast + j)
      END DO
      r (jmax) = r (jwrap)
      RETURN
      END SUBROUTINE Circle_1d
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    FUNCTION Eta_lambda_vgamma ( kbeg, kend, kcbeg, kcend, gamma)  RESULT (out)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: kbeg, kend, kcbeg, kcend
    REAL (rprec), INTENT (IN) :: gamma
    REAL (rprec) :: out (isize,jsize,iesize)
!______________________________________________________________________________
!
!   Subroutine computes quantity ETA*ABS(ALAM)*V**GAMMA at each grid point
!   for electrons and ions separately. KBEG, KEND, KCBEG, KCEND can be used
!   to restrict species to be included in the sum (use 1, ksize, 1, kcsize for
!   no restrictions). GAMMA is an input parameter:
!   ** if GAMMA = 0, then the computed sum is the adiabatic parameter PVGAMMA
!   ** if GAMMA = -5/3, compute energy density (or pressure without the 2/3 factor)
!   ** if GAMMA = -2/3, compute total energy of particles
!______________________________________________________________________________
!
    INTEGER (iprec) :: nbi, k, m, mbeg, mend, ipmax, ncount, i,j,n, ie,kc
    REAL (rprec) :: q, bimax, bicrss (100), charge
!
!
    out = zero
!
    DO ie = 1, iesize
!
        IF (ie == 1) THEN
           charge = - one
        ELSE
           charge = + one
        END IF
!
!       I. Compute sum for plasma on inner edges, electrons and ions separately :

!
!   II. Compute the sum for grid_based electrons or ions:
!
        DO j = 1, jsize
        DO i = 1, isize
!           IF (REAL(i,rprec) < Bndy(bndloc, REAL(j,rprec)) ) CYCLE
           IF (REAL(i,rprec) < bndloc(j) ) CYCLE
           DO kc = kcbeg, kcend
              q = alamc(kc) / charge
              IF (q > zero) THEN
                 out (i,j,ie) = out (i,j,ie) + ABS (alamc(kc) * eeta(i,j,kc))
              END IF
           END DO
        END DO
        END DO
!
    END DO
!
!
    DO ie = 1, iesize
       CALL Circle (out (:,:,ie))
    END DO
!
!
    DO ie = 1, iesize
       out (:,:,ie) = out(:,:,ie) * (vm (:,:)**((-three/two)*gamma))
    END DO
!
    RETURN
    END FUNCTION Eta_lambda_vgamma
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    FUNCTION Get_imin_for_grid (j)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: j
    INTEGER (iprec) :: Get_imin_for_grid
!
    Get_imin_for_grid = MIN (imin_j(Bjmod(j,jwrap,jsize)), &
                             imin_j(Bjmod(j+1,jwrap,jsize)), &
                             imin_j(Bjmod(j-1,jwrap,jsize)) )
    RETURN
    END FUNCTION Get_imin_for_grid
!
!
!
!
     FUNCTION Lt_from_aloct (phi)
     IMPLICIT NONE
     REAL (rprec), INTENT (IN) :: phi
     REAL (rprec) :: Lt_from_aloct
!
!    Convert an RCM phi angle (aloct in ionosphere) to MLT
!    Output (result) is:   0.0  <=  result < 24.00
!
     IF (phi < zero .OR. phi > pi_two) THEN
        WRITE (*,*) 'IN LT_FROM_ALOCT, PHI IS OUT OF BOUNDS'
        STOP
     ELSE
        Lt_from_aloct = MODULO((phi-pi)*RTH, 24.0_rprec)
        Lt_from_aloct = MODULO(Lt_from_aloct, 24.0_rprec)
     END IF
     RETURN
     END FUNCTION Lt_from_aloct

!
!
     SUBROUTINE Read_plasma_H5
        use ioh5
        use files
        implicit none
        logical :: doSP
        type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
        integer :: nvari, tauDim, Nk, Nm,Nl,Ne
        integer :: dims(4) ! update when add higher dimensions

        write(*,*) "Read rcmconfig.h5...."

        doSP = .false.
        call ClearIO(IOVars) !Reset IO chain
        call AddInVar(IOVars,"alamc") 
        call AddInVar(IOVars,"ikflavc") 
        call AddInVar(IOVars,"fudgec")
        call ReadVars(IOVars,doSP,RCMGAMConfig)

        !Store data for energy channels
        alamc(:) = IOVars(FindIO(IOVars, "alamc"))%data
        ikflavc(:) = IOVars(FindIO(IOVars, "ikflavc"))%data
        fudgec(:) = IOVars(FindIO(IOVars, "fudgec"))%data

        ! reset to make sure species if ikflav ==1 alamc is set to negative, for electrons
        where(ikflavc==1)alamc = -abs(alamc)

        ! Check if rcmconfig.h5 is up-to-date 
        if (ioExist(RCMGAMConfig,"Tau1i")) then
            write(*,*) "An old rcmconfig.h5 is used. Please make a new one using genRCM.py."
            stop
        endif

        !Store data for wave models
        !Dimension check: only compatible with tau(MLT,L,Kp,Ek)
        if (ioExist(RCMGAMConfig,"Taui")) then
            EWMTauInput%useWM = .true.
            !Chorus wave
            call AddInVar(IOVars,"Kpi")  
            call AddInVar(IOVars,"MLTi") 
            call AddInVar(IOVars,"Li") 
            call AddInVar(IOVars,"Eki") 
            call AddInVar(IOVars,"Taui") 
            call ReadVars(IOVars,doSP,RCMGAMConfig)
            tauDim = IOVars(FindIO(IOVars, "Taui"))%Nr
            if (tauDim /= 4) then
                write(*,*) "tauDim:",tauDim
                write(*,*) 'Currently only support tau model files in the form tau(Kp,MLT,L,Ek)'
                write(*,*)"tau:",IOVars(FindIO(IOVars, "Taui"))%dims
                stop
            endif

            dims = IOVars(FindIO(IOVars, "Taui"))%dims(1:tauDim)
            Nk   = IOVars(FindIO(IOVars, "Kpi"))%N
            Nm   = IOVars(FindIO(IOVars, "MLTi"))%N
            Nl   = IOVars(FindIO(IOVars, "Li"))%N
            Ne  = IOVars(FindIO(IOVars, "Eki"))%N
            if (Nk /=  dims(1) .or. Nm /= dims(2) .or. Nl /= dims(3) .or. Ne /= dims(4)) then
                write(*,*) "dims:",dims,"Nk:",Nk,"Nm:",Nm,"Nl:",Nl,"Ne:",Ne
                write(*,*) 'Dimensions of tau arrays are not compatible'
                stop
            endif

           !Store arrays
           EWMTauInput%ChorusTauInput%Nk = Nk
           EWMTauInput%ChorusTauInput%Nm = Nm
           EWMTauInput%ChorusTauInput%Nl = Nl
           EWMTauInput%ChorusTauInput%Ne = Ne


           if(allocated(EWMTauInput%ChorusTauInput%Kpi)) deallocate(EWMTauInput%ChorusTauInput%Kpi)
           if(allocated(EWMTauInput%ChorusTauInput%MLTi)) deallocate(EWMTauInput%ChorusTauInput%MLTi)
           if(allocated(EWMTauInput%ChorusTauInput%Li)) deallocate(EWMTauInput%ChorusTauInput%Li)
           if(allocated(EWMTauInput%ChorusTauInput%Eki)) deallocate(EWMTauInput%ChorusTauInput%Eki)
           if(allocated(EWMTauInput%ChorusTauInput%taui)) deallocate(EWMTauInput%ChorusTauInput%taui)

           allocate(EWMTauInput%ChorusTauInput%Kpi(Nk))
           allocate(EWMTauInput%ChorusTauInput%MLTi(Nm))
           allocate(EWMTauInput%ChorusTauInput%Li(Nl))
           allocate(EWMTauInput%ChorusTauInput%Eki(Ne))
           allocate(EWMTauInput%ChorusTauInput%taui(Nk,Nm,Nl,Ne))

           call IOArray1DFill(IOVars,"Kpi",EWMTauInput%ChorusTauInput%Kpi)
           call IOArray1DFill(IOVars,"MLTi",EWMTauInput%ChorusTauInput%MLTi)
           call IOArray1DFill(IOVars,"Li", EWMTauInput%ChorusTauInput%Li)
           call IOArray1DFill(IOVars,"Eki",EWMTauInput%ChorusTauInput%Eki)
           call IOArray4DFill(IOVars,"Taui",EWMTauInput%ChorusTauInput%taui)

           call ClearIO(IOVars)

           !Array order check: array is in acsending order
           !Chorus
           if(EWMTauInput%ChorusTauInput%Kpi(1) > EWMTauInput%ChorusTauInput%Kpi(Nk)) then
              write(*,*) "Kp: ",EWMTauInput%ChorusTauInput%Kpi
              write(*,*) "reorder wave model so Kp is in ascending order"
              stop
           end if

           if(EWMTauInput%ChorusTauInput%Li(1) > EWMTauInput%ChorusTauInput%Li(Nl)) then
              write(*,*) "L: ",EWMTauInput%ChorusTauInput%Li
              write(*,*) "reorder wave model so L shell is in ascending order"
              stop
           end if

           if(EWMTauInput%ChorusTauInput%MLTi(1) > EWMTauInput%ChorusTauInput%MLTi(Nm)) then
              write(*,*) "MLT: ",EWMTauInput%ChorusTauInput%MLTi
              write(*,*) "reorder wave model so MLT is in ascending order"
              stop
           end if

           if(EWMTauInput%ChorusTauInput%Eki(1) > EWMTauInput%ChorusTauInput%Eki(Ne)) then
              write(*,*) "Ek: ",EWMTauInput%ChorusTauInput%Eki
              write(*,*) "reorder wave model so Ek is in ascending order"
              stop
           end if        
 
        endif

        
      END SUBROUTINE Read_plasma_H5
!
!
!
      FUNCTION Dipole_Bfield (theta, phi, arg)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: theta, phi
      INTEGER (iprec), INTENT (IN) :: arg
      REAL (rprec) :: Dipole_Bfield
!
!       THETA and PHI are in radians, THETA measured from n.pole down,
!       PHI measured from noon to dusk to midnight etc.
!       Distances RMIN, XMIN, and YMIN are in units of RE
!       Since in the RCM, BESU is in [nT], and the factor of RE is ommited
!       from the formula, VM has units of (RE/nT)**(-2/3)
!
      REAL (rprec) :: rmin, xmin, ymin, bmin, vm
!
      rmin = one / SIN(theta)**2
      xmin = rmin * COS (phi)
      ymin = rmin * SIN (phi)
      bmin = besu * (one/rmin**3)
      vm   = (32.0_rprec/35.0_rprec* rmin**4 / Besu *  &
              SQRT(one-one/rmin)* &
              (one+half/rmin+three/eight/rmin**2+five/eight/two/rmin**3) &
              ) ** (-two/three)
!
      IF (arg == 1) THEN
         Dipole_Bfield = xmin
      ELSE IF (arg == 2) THEN
         Dipole_Bfield = ymin
      ELSE IF (arg == 3) THEN
         Dipole_Bfield = bmin
      ELSE IF (arg == 4) THEN
         Dipole_Bfield = vm
      ELSE
         STOP 'ILLEGAL ARGUMENT FOR DIPOLE_BFIELD'
      END IF
!
      RETURN
      END FUNCTION Dipole_Bfield
!
!
!
!
!
      FUNCTION Deriv_i (a, imin_j) RESULT (d_di)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: imin_j(:)
      REAL (rprec), INTENT (IN) :: a(:,:)
      REAL (rprec) :: d_di (SIZE(a,DIM=1), SIZE(a,DIM=2))
!
!_________________________________________________________________________
!     The idea is to use central differences for a second-order accuracy
!     of the first-order derivative, if we can. If one of the neighboring
!     points is outside the boundary, use forward or back difference
!     (first-order accuracy) if we can. If both points are outside, set
!     derivative to zero.
!_________________________________________________________________________
!
      INTEGER (iprec) :: i, j, j_size, i_size
!
      i_size = SIZE (a, DIM=1)
      j_size = SIZE (a, DIM=2)
!
      DO j = 1, j_size
!
         d_di (1:imin_j(j)-1,j) = 0.0_rprec
!
         i = imin_j(j)
         d_di (i,j) = -1.5_rprec*a(i,j) + 2.0_rprec*a(i+1,j) - 0.5_rprec*a(i+2,j)
!
         DO i = imin_j(j)+1, i_size - 1
            d_di (i,j) = 0.5_rprec*(a(i+1,j) - a(i-1,j))
         END DO
!
         i = i_size
         d_di (i,j) = +1.5_rprec*a(i,j) - 2.0_rprec*a(i-1,j) + 0.5_rprec*a(i-2,j)
!
      END DO
!
      CALL Circle (d_di)
!
      RETURN
      END FUNCTION Deriv_i
!
!
!
      FUNCTION Deriv_j (a, imin_j, j1, j2, error_value) RESULT (d_dj)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: a (:,:), error_value
      INTEGER (iprec), INTENT (IN) :: imin_j(:), j1, j2
      REAL (rprec) :: d_dj (SIZE(a,DIM=1), SIZE(a,DIM=2))
!___________________________________________________________________________
!
!     Take derivative of A with respect to J. The idea is to use
!     grid points that are inside the modeling region. Stanislav 10/27/2000.
!     3/6/2001: modified to be the same as Frank's version in rcm296.
!     Notice that pt_jpp is defined differently (>) than pt_jp(>=), and ditto
!     for pt_jmm and pt_jm. Took this from Frank's version.
!___________________________________________________________________________
!
      INTEGER (iprec) :: i_bnd, i, j, jm, jmm, jp, jpp, j_size
      LOGICAL :: pt_jp, pt_jpp, pt_jm, pt_jmm
!
!
      i_bnd = MAXVAL (imin_j)
      j_size = SIZE (imin_j)
!
      DO j = j1, j2
!
         jp = j + 1
         jm = j - 1
!
         d_dj (1:imin_j(j)-1,j) = 0.0_rprec
!
         DO i = imin_j(j), i_bnd-1
!
            pt_jp  = (i >= imin_j(jp))
            pt_jm  = (i >= imin_j(jm))
!
            IF (pt_jp .AND. pt_jm ) THEN
               d_dj (i,j) = half*(a (i,jp) - a (i,jm))
            ELSE IF (pt_jp ) THEN !j-1 is outside
               d_dj (i,j) = a(i,jp) - a(i,j)
            ELSE IF (pt_jm ) THEN  ! j+1 is outside
               d_dj (i,j) = (a (i,j) - a (i,jm))
            ELSE
               d_dj (i,j) = error_value
            END IF 
!
         END DO
!
         d_dj (i_bnd:,j) = 0.5_rprec*( a(i_bnd:,jp) - a(i_bnd:,jm))
      END DO
!
      CALL Circle (d_dj)
!
      RETURN
      END FUNCTION Deriv_j
!
!
!
      SUBROUTINE Rcm (itimei_in, itimef_in, nstep_in, icontrol, stropt, nslcopt, iXML)
!---------------------------------------------
! notes
!     icontrol controls behaviour of the RCM
!       0:  initialize RCM size params and go back
!       (RCMINIT) 1:  initialize RCM grid, energy channels, quit:
!       2:  read in inputs, quit:
!       3: Set initial conditions on plasma (edges and grid-based)
!       4:  run RCM from itimei to itimef with time step idt, quit:
!       ICONWRITERESTART = 31337: write a restart record to RCM
!       ICONWRITEOUTPUT = ICONWRITERESTART + 1: write an HDF5 output
!       ICONRESTART = ICONWRITERESTART - 1:  Read HDF5 restart


      USE xml_input
      USE rcm_timing_module
      IMPLICIT NONE
!
      type(XML_Input_T), intent(in), optional :: iXML
      REAL (rprec), INTENT (IN) :: itimei_in, itimef_in    !, &   
                                      !idt_in, idt1_in,& 
                                      !idt2_in, 
      INTEGER (iprec), INTENT (IN) :: nstep_in
      INTEGER (iprec), INTENT (IN) :: icontrol
      character(len=*), intent(in), optional :: stropt
      integer(iprec)  , intent(in), optional :: nslcopt
      CHARACTER(LEN=8) :: real_date
      CHARACTER (LEN=8) :: time_char
      CHARACTER(LEN=10) ::real_time
      CHARACTER(LEN=80) :: ST='', PS='', HD='', string_null=''
      LOGICAL :: FD, logical_flag
!                                                                       
!
      REAL (rprec), SAVE :: itimei, itimef      !, idt, idt1, idt2 
!      REAL (rprec), SAVE :: itout1, itout2,  itcln,  i_time
      INTEGER (iprec), SAVE :: idebug, k, kc, n, nstep
      INTEGER (iprec) :: i_avg, i_step
      REAL (rprec) :: dt,wAvg
      REAL (rprec), PARAMETER :: tinyT = 1e-6     !10.0*machine_tiny 

      CALL SYSTEM_CLOCK (timer_start(1), count_rate)

      itimef = itimef_in
      itimei = itimei_in
      nstep = nstep_in
     
      if (doRCMVerbose) then
        write(6,*)'RCM: itimei= ',itimei,' seconds',' itimef= ',itimef,' seconds ','nStep= ',nstep
      endif

      IF (icontrol == ICONWRITERESTART) then  ! write a restart record to RCM
         call WriteRCMH5(stropt,nslcopt,isRestart=.true.)
         return
      ENDIF
      IF (icontrol == ICONWRITEOUTPUT) then  ! write an HDF5 output
        call WriteRCMH5(stropt,nslcopt,isRestart=.false.)
        return
      ENDIF
      IF (icontrol == ICONRESTART) then
        !Read HDF5 restart
        call ReadRCMRestart(stropt,nslcopt)
        return
      ENDIF

   IF (icontrol == 0) then  ! initialize RCM size params and go back:


      !  Grid limit parameters
      i1   = 2          ! will reset anyway
      i2   = isize - 1
      j1   = jwrap
      j2   = jsize - 1
      iint = 1 
      jint = 1      

      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN

   END IF

   IF (icontrol == 1) then   ! initialize RCM grid, energy channels, quit:
   !   CALL Read_grid ()
      CALL Read_plasma_H5()
      
      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN
   END IF


   IF (icontrol == 2) then ! read in inputs, quit:
      
      if(present(iXML)) then
        call RCM_Params_XML(iXML)
      else
        call RCM_Params_XML()
      endif
      CALL Read_dktime_H5(L_dktime)

      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN

   END IF


   IF (icontrol == 3) then   !-->  Set initial conditions on plasma (grid-based):
!
      ! Open file for formatted output and do initial print out :
      CALL Date_and_time (real_date, real_time)

      i1 = imin + 1
        
      IF (itimei > tinyT) imin_j = CEILING (bndloc)


      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN

   END IF



   IF (icontrol == 4) then  ! run RCM from itimei to itimef with time step idt, quit:
      call Tic("Main_Loop")
      
      !Do check
      if (.not. any(advChannel)) then
        write(*,*) "RCM has no good channels, dying"
        stop
      endif

      CALL SYSTEM_CLOCK (timer_start(2), count_rate)

      !NOTE: v_avg behaves differently than birk_avg
      !v_avg is the running average over numerous rcm couplings, birk_avg is over just the current one

      if (doVAvgInit) then
        v_avg = v !Initialize v_avg = v at first time
        doVAvgInit = .FALSE.
      else
        !Figure out weighting for exponential moving average (EMA)
        !Want weighting such that ~95% of the weight comes from the last dtAvg seconds
        dt = (itimef-itimei) !Full RCM step
        wAvg = 1.0 - exp(-3*dt/max(dtAvg_v,dt))
        v_avg = wAvg*v + (1-wAvg)*v_avg
      endif

      birk_avg = zero
      eeta_avg = zero
      i_avg    = 0

      deleeta = 0.0
      lossratep = 0.0
      lossmodel= -1.0

!*******************  main time loop  *************************
!
      !Save most recent coupling dt
      dtMHD = itimef-itimei
      dt = (itimef - itimei)/REAL(nstep) 

      if (doRCMVerbose) then
         write(6,*)'RCM: substep length = ',dt,' seconds'
      endif
!
      !Q: Does this have any point since it gets recalculated in move-plasma?
      !fac = 1.0E-3_rprec * bir * alpha * beta * dlam * dpsi * ri**2 * signbe
!
      birk_avg = birk_avg + birk
 
      IF (nstep < 1) STOP 'Number of substep in RCM should be at least 1'
      
      !NOTE: Pushing averaging into move-plasma to avoid multiple recalculations of static arrays
      call Move_plasma_grid_MHD(dt,nstep)
      
      CALL Comput (itimef-itimei)
      if (doRCMVerbose) then
         write(6,*)'RCM: : finishing Comput'
      endif
      birk_avg = (birk_avg + birk)/2.     ! brik_avg takes two data points at itimei and itimef

      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      CALL SYSTEM_CLOCK (timer_stop(2), count_rate)      
      timer_values (2) = (timer_stop (2) - timer_start (2))/count_rate

      call Toc("Main_Loop")

      RETURN

   END IF

   
   WRITE (*,*) ' RCM was called with an invalid value of Icontrol, aborting ...'
   STOP

      CONTAINS
!
        !HDF5 Restart reader
        subroutine ReadRCMRestart(runid,nStp)
          use ioh5
          use files
          implicit none
          character(len=*), intent(in) :: runid
          integer(iprec), intent(in) :: nStp
          logical :: doSP !Do single precision
          character(len=strLen) :: H5File
          type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
          integer(iprec) :: nvar,nres,Ni,Nj,Nk

        !Prepare for reading
          doSP = .false. !Restarts are always double precision
          nres = nStp-1 !nStp holds number for *NEXT* restart output
          if (nres == -1) then
            !Use sym link
            H5File = trim(runid) // ".RCM.Res.XXXXX.h5"
          else
            !Use actual #
            write (H5File, '(A,A,I0.5,A)') trim(runid), ".RCM.Res.", nres, ".h5"
          endif
          write(*,*) 'Restarting RCM with file, ', trim(H5File)
          call ClearIO(IOVars) !Reset IO chain

        !List variables to read
          !Scalars (need to specify integers), order doesn't matter

          call AddInVar(IOVars,"itimei",vTypeO=IOINT )
          call AddInVar(IOVars,"isize" ,vTypeO=IOINT )
          call AddInVar(IOVars,"jsize" ,vTypeO=IOINT )
          call AddInVar(IOVars,"ksize" ,vTypeO=IOINT )
          call AddInVar(IOVars,"cmax"  ,vTypeO=IOREAL)
          call AddInVar(IOVars,"fmeb"  ,vTypeO=IOREAL)
          call AddInVar(IOVars,"fstoff",vTypeO=IOREAL)
          call AddInVar(IOVars,"fdst"  ,vTypeO=IOREAL)
          call AddInVar(IOVars,"fclps" ,vTypeO=IOREAL)
          call AddInVar(IOVars,"vdrop" ,vTypeO=IOREAL)
          call AddInVar(IOVars,"kp"    ,vTypeO=IOREAL)
          call AddInVar(IOVars,"i_avg" ,vTypeO=IOREAL)

          !Arrays
          call AddInVar(IOVars,"rcmxmin"  )
          call AddInVar(IOVars,"rcmymin"  )
          call AddInVar(IOVars,"rcmzmin"  )
          call AddInVar(IOVars,"rcmvm"    )
          call AddInVar(IOVars,"rcmbmin"  )
          call AddInVar(IOVars,"rcmbndloc")

          call AddInVar(IOVars,"rcmetac"   )
          call AddInVar(IOVars,"rcmeeta"   )
          call AddInVar(IOVars,"rcmeetaavg")
          call AddInVar(IOVars,"rcmlosspre")
          call AddInVar(IOVars,"rcmlossmod")
          !call AddInVar(IOVars,"rcmDpp")


          call AddInVar(IOVars,"rcmpedlam" )
          call AddInVar(IOVars,"rcmpedpsi" )
          call AddInVar(IOVars,"rcmhall"   )
          call AddInVar(IOVars,"rcmeavg"   )
          call AddInVar(IOVars,"rcmeflux"  )
          call AddInVar(IOVars,"rcmnflux"  )
          call AddInVar(IOVars,"rcmbirk"   )
          call AddInVar(IOVars,"rcmbirkavg")

          call AddInVar(IOVars,"rcmv")
          call AddInVar(IOVars,"rcmvavg")

          call AddInVar(IOVars,"alpha")
          call AddInVar(IOVars,"aloct")
          call AddInVar(IOVars,"colat")
          call AddInVar(IOVars,"beta")
          call AddInVar(IOVars,"bir")
          call AddInVar(IOVars,"sini")

          !Extra stuff (not in write arrays)
          call AddInVar(IOVars,"alamc")
          
        !Now do actual reading
          call ReadVars(IOVars,doSP,H5File)
        !Do some testing to make sure sizes match
          nvar = FindIO(IOVars,"rcmeeta",doFailO=.true.)
          Ni = IOVars(nvar)%dims(1)
          Nj = IOVars(nvar)%dims(2)
          Nk = IOVars(nvar)%dims(3)
          if ( (isize /= Ni) .or. (jsize /= Nj) .or. (ksize /= Nk) ) then
            write(*,*) 'RCM Restart Mismatch!'
            write(*,*) 'Input size: ',Ni,Nj,Nk
            write(*,*) 'RCM   size: ',isize,jsize,ksize
            stop
          endif

        !Parse data and put it where it goes, need to do each variable
          !Scalars
          itimei = GetIOInt(IOVars,"itimei")
          
          cmax   = GetIOReal(IOVars,"cmax")
          fmeb   = GetIOReal(IOVars,"fmeb")
          fstoff = GetIOReal(IOVars,"fstoff")
          fdst   = GetIOReal(IOVars,"fdst")
          fclps  = GetIOReal(IOVars,"fclps")
          vdrop  = GetIOReal(IOVars,"vdrop")
          i_avg  = GetIOReal(IOVars,"i_avg")

          !Pull 2D arrays
          call IOArray2DFill(IOVars,"rcmxmin",xmin)
          call IOArray2DFill(IOVars,"rcmymin",ymin)
          call IOArray2DFill(IOVars,"rcmzmin",zmin)
          call IOArray2DFill(IOVars,"rcmbmin",bmin)
          
          call IOArray2DFill(IOVars,"rcmv",v)
          call IOArray2DFill(IOVars,"rcmvavg",v_avg)
          doVAvgInit = .FALSE. !Don't need to start fresh v_avg

          call IOArray2DFill(IOVars,"rcmvm",vm)

          call IOArray2DFill(IOVars,"rcmbirk",birk)
          call IOArray2DFill(IOVars,"rcmbirkavg",birk_avg)

          call IOArray2DFill(IOVars,"rcmhall",hall)
          call IOArray2DFill(IOVars,"rcmpedlam",pedlam)
          call IOArray2DFill(IOVars,"rcmpedpsi",pedpsi)

          !Disable replacing geometry (calculated fresh anyways)
          !Better for when using upscaled restart
          !call IOArray2DFill(IOVars,"alpha",alpha)
          !call IOArray2DFill(IOVars,"aloct",aloct)
          !call IOArray2DFill(IOVars,"colat",colat)
          !call IOArray2DFill(IOVars,"beta",beta)
          !call IOArray2DFill(IOVars,"bir",bir)
          !call IOArray2DFill(IOVars,"sini",sini)
          !call IOArray2DFill(IOVars,"rcmDpp",Dpp)


          !Pull 1D arrays
          call IOArray1DFill(IOVars,"rcmetac",etac)
          call IOArray1DFill(IOVars,"rcmbndloc",bndloc)
          call IOArray1DFill(IOVars,"alamc",alamc)

          !Pull 3D arrays
          call IOArray3DFill(IOVars,"rcmeavg",eavg)
          call IOArray3DFill(IOVars,"rcmeeta",eeta)
          call IOArray3DFill(IOVars,"rcmeetaavg",eeta_avg)
          call IOArray3DFill(IOVars,"rcmeflux",eflux)
          if(ioExist(H5File,"rcmnflux")) then
            call IOArray3DFill(IOVars,"rcmnflux",nflux)
          endif
          if(ioExist(H5File,"rcmlosspre")) then
            call IOArray3DFill(IOVars,"rcmlosspre",lossratep)
          endif
          if(ioExist(H5File,"rcmlossmod")) then
            call IOArray3DFill(IOVars,"rcmlossmod",lossmodel)
          endif
          
        end subroutine ReadRCMRestart

        !HDF5 output routine
        !isRestart = Whether we're writing restart dump or regular output slice
        subroutine WriteRCMH5(runid,nStp,isRestart)
          use ioh5
          use files
          implicit none
          character(len=*), intent(in) :: runid
          integer(iprec), intent(in) :: nStp
          logical, intent(in) :: isRestart

          type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
          logical :: doSP !Do single precision output
          character(len=strLen) :: H5File,gStr,lnResF

        !Prepare for output
          !Reset IO chain
          call ClearIO(IOVars)
          !Distinguish output slices vs restarts
          if (isRestart) then
            doSP = .false. !Double precision restarts
            write (H5File, '(A,A,I0.5,A)') trim(runid), ".RCM.Res.", nStp, ".h5"
          else
            !Regular output
            doSP = .true.
            H5File = trim(runid) // ".rcm.h5"
            write (gStr, '(A,I0)') "Step#", nStp
          endif

        !Attributes
          call AddOutVar(IOVars,"time",1.0_rp*itimei)
          call AddOutVar(IOVars,"itimei",itimei)
          call AddOutVar(IOVars,"isize" ,isize )
          call AddOutVar(IOVars,"jsize" ,jsize )
          call AddOutVar(IOVars,"ksize" ,ksize )
          call AddOutVar(IOVars,"cmax"  ,cmax  )
          call AddOutVar(IOVars,"fmeb"  ,fmeb  )
          call AddOutVar(IOVars,"fstoff",fstoff)
          call AddOutVar(IOVars,"fdst"  ,fdst  )
          call AddOutVar(IOVars,"fclps" ,fclps )
          call AddOutVar(IOVars,"vdrop" ,vdrop )
          call AddOutVar(IOVars,"i_avg" ,i_avg )
          call AddOutVar(IOVars,"dtCpl" ,dtMHD )

        !Arrays
          call AddOutVar(IOVars,"rcmxmin",xmin)
          call AddOutVar(IOVars,"rcmymin",ymin)
          call AddOutVar(IOVars,"rcmzmin",zmin)
          call AddOutVar(IOVars,"rcmvm"  ,vm  )
          call AddOutVar(IOVars,"rcmbmin",bmin)
          call AddOutVar(IOVars,"rcmbndloc",bndloc)

          call AddOutVar(IOVars,"rcmetac"   ,etac)
          call AddOutVar(IOVars,"rcmeeta"   ,eeta)
          call AddOutVar(IOVars,"rcmeetaavg",eeta_avg)
          call AddOutVar(IOVars,"rcmdeleeta",deleeta)
          call AddOutVar(IOVars,"rcmlosspre",lossratep)
          call AddOutVar(IOVars,"rcmDpp",Dpp)
          call AddOutVar(IOVars,"rcmlossmod",lossmodel)

          call AddOutVar(IOVars,"rcmpedlam" ,pedlam  )
          call AddOutVar(IOVars,"rcmpedpsi" ,pedpsi  )
          call AddOutVar(IOVars,"rcmhall"   ,hall    )
          call AddOutVar(IOVars,"rcmeavg"   ,eavg    )
          call AddOutVar(IOVars,"rcmeflux"  ,eflux   )
          call AddOutVar(IOVars,"rcmnflux"  ,nflux   )
          call AddOutVar(IOVars,"rcmbirk"   ,birk    )
          call AddOutVar(IOVars,"rcmbirkavg",birk_avg)

          call AddOutVar(IOVars,"rcmv",v)
          call AddOutVar(IOVars,"rcmvavg",v_avg)

        !Extra stuff not in write_array
          call AddOutVar(IOVars,"alamc",alamc)
          call AddOutVar(IOVars,"aloct",aloct)
          call AddOutVar(IOVars,"colat",colat)
          call AddOutVar(IOVars,"alpha",alpha)
          call AddOutVar(IOVars,"beta" ,beta )
          call AddOutVar(IOVars,"bir"  ,bir  )
          call AddOutVar(IOVars,"sini" ,sini )

          if (doRCMVerboseH5) then
            !Good place to store useful but large 3D outputs
            call AddOutVar(IOVars,"rcmveff",last_veff,uStr="Volts")
            call AddOutVar(IOVars,"rcmocbDist",last_ocbDist)
          endif
          
        !Done staging output, now let er rip
          if (isRestart) then
            call AddOutVar(IOVars,"nRes",nStp)
            call CheckAndKill(H5File) !Always overwrite restarts
            call WriteVars(IOVars,doSP,H5File)
            !Create link to latest restart
            write (lnResF, '(A,A,A,A)') trim(runid), ".RCM.Res.", "XXXXX", ".h5"
            call MapSymLink(H5File,lnResF)
          else
            call WriteVars(IOVars,doSP,H5File,gStr)
          endif
        end subroutine WriteRCMH5

        subroutine RCM_Params_XML(iXML)
          use xml_input
          use strings

          type(XML_Input_T), intent(in), optional :: iXML
          character(len=strLen) :: inpXML
          type(XML_Input_T) :: xmlInp

          if(present(iXML)) then
            call iXML%GetFileStr(inpXML)
          else
            !Find input deck filename
            call getIDeckStr(inpXML)
          endif
          
          !Create new XML reader w/ RCM as root
          xmlInp = New_XML_Input(trim(inpXML),'Kaiju/RCM',.true.)

          call xmlInp%Set_Val(label%char,"sim/runid","MAGE sim")

          !Output
          call xmlInp%Set_Val(idebug,"output/idebug",1) ! 6.  0 <=> do disk printout
          call xmlInp%Set_Val(doRCMVerbose,"output/doDebug",doRCMVerbose)
          call xmlInp%Set_Val(doRCMVerboseH5,"output/doDebugH5",doRCMVerboseH5)
          !eflux
          call xmlInp%Set_Val(ifloor,"eflux/ifloor",.true.) ! 18. if true, install a floor for EFLUX
          call xmlInp%Set_Val(icorrect,"eflux/icorrect",.true.) ! 19. if true, make lat. correction to EFLUX

          !Grid
          call xmlInp%Set_Val(imin,"grid/imin",1)
          call xmlInp%Set_Val(ibnd_type,"grid/ibnd_type",4) ! 14.  type of bndy (1-eq.p, 2-iono)
          call xmlInp%Set_Val(ipcp_type,"grid/ipcp_type",13) ! 14.  type of bndy (1-eq.p, 2-iono)
          call xmlInp%Set_Val(nsmthi,"grid/nsmthi",0) ! 15.  How much to smooth cond in I
          call xmlInp%Set_Val(nsmthj,"grid/nsmthj",0) ! 16.  How much to smooth cond in J
          call xmlInp%Set_Val(L_move_plasma_grid,"grid/L_move_plasma_grid",.true.)

          !Catch-all params
          call xmlInp%Set_Val(ipot,"params/ipot",-1) !  6.  which potential solver to use
          call xmlInp%Set_Val(iwind,"params/iwind",0) !  9.  0 is no neutral winds
          call xmlInp%Set_Val(icond,"params/icond",3) ! 1 is active conductances, 2 is Hardy with kp, 3 is input

          call xmlInp%Set_Val(cmax,"params/cmax",3.0) ! in rcm_mod_balgn
          call xmlInp%Set_Val(eeta_cutoff,"params/eeta_cutoff",0.05) ! as a fraction

          !Charge exchange
          call xmlInp%Set_Val(kill_fudge,"chargex/kill_fudge",.false.) ! .true. means no loss
          if (kill_fudge) then
            fudgec = 0.0
          endif
          call xmlInp%Set_Val(L_dktime,"chargex/L_dktime",.true.)
          call xmlInp%Set_Val(sunspot_number,"chargex/sunspot_number",96.0)

          !Clawpack options
          call xmlInp%Set_Val(L_doOMPClaw,"clawpack/doOMPClaw",L_doOMPClaw)

          !Averaging timescale for plasmasphere
          call xmlInp%Set_Val(dtAvg_v,"plasmasphere/tAvg",60.0)

          call xmlInp%Set_Val(doNoBndFlow,"experimental/doNoBndFlow",.false.)
          call xmlInp%Set_Val(nNBFL,"experimental/NBFLayers",nNBFL)

          !Some values just setting
          tol_gmres = 1.0e-5
          itype_bf = 3 ! 1 is interpolate for HV, 2--MHD code, 3--receive through module
          i_advect = 3  ! 1-interpolate, 2rCLAWPACK/inter, 3-CLAWPACK
          i_eta_bc = 2! 1-time-dep. from file, 2-constant for run
          i_birk  = 1 ! birk calculation 1=default 3 = new
        end subroutine RCM_Params_XML

 !
      END SUBROUTINE Rcm
!
!

    !K: HDF5 version of lifetime reader
    SUBROUTINE Read_dktime_H5(L_dktime)
      use ioh5
      use files
      IMPLICIT NONE
      LOGICAL, INTENT (IN) :: L_dktime
      logical :: doSP
      type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
      
      !real(rprec) :: dktime2(irdk,inrgdk,isodk, iondk)
      Call CheckFileOrDie(RCMGAMConfig,"RCM-Config H5 file does not exist.")

      if (L_dktime) then
        !Read from HDF5
        doSP = .false.
        call ClearIO(IOVars) !Reset IO chain
        call AddInVar(IOVars,"dktable")
        call ReadVars(IOVars,doSP,RCMGAMConfig)

        dktime = reshape(IOVars(1)%data,[irdk,inrgdk,isodk, iondk])

      endif
    END SUBROUTINE Read_dktime_H5



      FUNCTION Cexrat (isp,enrg,rloc,ssn,dktime,irdk,inrgdk,isoldk, &
                       iondk)
      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: isp, irdk, inrgdk, isoldk, iondk
      REAL(rprec), INTENT (IN) :: enrg, rloc, ssn, dktime (irdk,inrgdk,isoldk,iondk)
      REAL(rprec) :: Cexrat
!
!-------------------------------------------------------------------------
!  copyright rice university, 1993
!
!  version 1.00                                 05.09.90
!          2.00                                 02.04.90
!                                       msm delivery version
!          2.10                                 06.11.93
!               error output routed to unit 9
!
!  programmer: r. w. spiro
!
!  purpose:  function subprogram to return charge exchange loss rate
!          (sec**(-1)) for ions of species isp, energy enrg (ev) at
!          l=rloc (re) for sunspot number ssn.  this routine is based
!          on a table generated by james bishop of u. of michigan.
!
!  calling parameters
!        isp       species identifier
!                    isp=2 for h+ ions
!                    isp=3 for o+ ions
!        enrg      energy in ev
!        rloc      radial location (re)
!        ssn       sunspot number
!        dktime    table of ion decay times
!        irdk      radial dimension of dktime array
!        inrgdk    energy dimension of dktime array
!        isoldk    sunspot number dimension of dktime array
!        iondk     number of ion species in dktime array
!--------------------------------------------------------------------------------
!
      INTEGER(iprec), PARAMETER ::irsiz=18,inrgsz=13,isolsz=2,ionsiz=2
      REAL(rprec) ::  elgvec(inrgsz), rvec(irsiz),ssnvec(2), &
               enrglg, br, bnrg, ssnuse, bssn, decayt
      INTEGER(iprec) :: ispndx, ir, inrg
!
      DATA elgvec /2.50,2.75,3.00,3.25,3.50,3.75,4.00, &
                   4.25,4.50,4.75,5.00,5.25,5.50/
!
      DATA rvec /1.50,2.00,2.50,3.00, &
                 3.50,4.00,4.50,5.00, &
                 5.50,6.00,6.50,7.00, &
                 7.50,8.00,8.50,9.00, &
                 9.50,10.00/
!
      DATA ssnvec /0.0,100./
!
!
      IF (irsiz /= irdk .OR. inrgsz /= inrgdk .OR. &
          ionsiz /= iondk .OR. isolsz /= isoldk) THEN
         write(*,*) 'dimension error in function cexrat'
         write(*,*) 'irdk,inrgdk,iondk,isoldk',irdk,inrgdk,iondk,isoldk
         write(*,*) 'irsiz,inrgsz,ionsiz,isolsz',irsiz,inrgsz,ionsiz,isolsz
         write(*,*) 'stopping program in cexrat'
         STOP
      END IF
!
      enrglg = LOG10(enrg) !  work with log10 of particle energy
      ispndx=isp-1
!
      if_1: IF (rloc <= rvec(1)) THEN !  find br for interpolation
         br=1.0
      ELSE IF (rloc > rvec(irdk)) THEN
         br=irdk
      ELSE
         do_1: DO ir=1,irdk-1
            IF (rloc <= rvec(ir+1)) THEN
               br=ir+(rloc-rvec(ir))/(rvec(ir+1)-rvec(ir))
               EXIT do_1
            END IF
         END DO do_1
      END IF if_1
!
      if_2: IF (enrglg.le.elgvec(1)) THEN !  find bnrg for interpolation
         bnrg = 1.0
      ELSE IF (enrglg > elgvec(inrgdk)) THEN
         bnrg = inrgdk
      ELSE
         do_2: DO inrg=1,inrgdk-1
            IF (enrglg <= elgvec(inrg+1)) THEN
               bnrg=inrg+(enrglg-elgvec(inrg))/(elgvec(inrg+1)-elgvec(inrg))
               EXIT do_2
            END IF
         END DO do_2
      END IF if_2
!
!**********  change 9/30/91  *****************************************
!  if ssn.gt.ssnvec(2), then use ssnvec(2) for ssn
      ssnuse=ssn
      IF (ssnuse > ssnvec(2)) ssnuse=ssnvec(2)
!
!*********  end change  9/30/91  ************************************
!
!  find bssn for interpolation
      bssn=1.0+(ssnuse-ssnvec(1))/(ssnvec(2)-ssnvec(1))
!
!  decayt is decay time in seconds
    
      decayt = G3ntrp (dktime(1_iprec,1_iprec,1_iprec,ispndx),irdk,inrgdk,isoldk,br,bnrg,bssn)
!
      IF (ABS(decayt) < 1.0E-20) THEN
         write(*,*) 'decayt is less than 1.e-20 sec in cexrat'
         write(*,*) 'decayt=',decayt,'  br=',br,'  bnrg=',bnrg,'bssn=',bssn
         write(*,*) 'isp=',isp,'  enrg=',enrg,' rloc=',rloc,' ssn=',ssn
         write(*,*) 'ssnuse=',ssnuse
      END IF
!
!  to get charge exchange rate (sec**9-1)) cexrat, invert decayt
!
      cexrat=1.0/decayt
      RETURN
      END FUNCTION Cexrat
!
!
      FUNCTION G3ntrp (a,imax,jmax,kmax,bi,bj,bk)
      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: imax, jmax, kmax
      REAL(rprec), INTENT (IN) :: a(imax,jmax,kmax), bi, bj, bk
      REAL(rprec) :: G3ntrp
!
!---------------------------------------------------------------------------
!  copyright Rice University, 1993
!
!  VERSION 1.00                            DATE: 01.11.88
!          1.01A                                 02.02.89
!          2.00  MSM DELIVERY VERSION            01.28.93
!
!  PURPOSE: FUNCTION SUBPROGRAM TO PERFORM A GENERAL 3-D LINEAR
!           INTERPOLATION OF ARRAY A(I,J,K) AT PT(BV(1),BV(2),BV(3))
!
!  INPUT:
!       A          3-D ARRAY TO BE INTERPOLATED
!       IMAX       I DIMENSION OF ARRAY A
!       JMAX       J DIMENSION OF ARRAY A
!       KMAX       K DIMENSION OF ARRAY A
!       BI         FLOATING POINT VALUE TO INTERPOLATE IN I DIMENSION
!       BJ         FLOATING POINT VALUE TO INTERPOLATE IN J DIMENSION
!       BK         FLOATING POINT VALUE TO INTERPOLATE IN K DIMENSION
!
!
!  OUTPUT:
!       G3NTRP     INTERPOLATED VALUES OF ARRAY A
!----------------------------------------------------------------------
!
!
      INTEGER(iprec) ::  ndx(3),ndim(3), kstop, jstop, L, i, j, k
      REAL(rprec) ::  BV(3),COEF(3,2), fndx
!
      NDIM(1)=IMAX
      NDIM(2)=JMAX
      NDIM(3)=KMAX
      BV(1)=BI
      BV(2)=BJ
      BV(3)=BK
      DO L=1,3
         NDX(L)=BV(L)
         IF(NDX(L).LT.1) NDX(L)=1
         IF(NDX(L).GT.NDIM(L)-1) NDX(L)=NDIM(L)-1
         IF(NDX(L).LE.0) NDX(L)=1
         FNDX=REAL(NDX(L))
         COEF(L,1)=1.-BV(L)+FNDX
         COEF(L,2)=BV(L)-FNDX
      END DO
!
      G3NTRP=0.
      kstop = MIN(KMAX,2)
      jstop = MIN(JMAX,2)
      DO I=1,2
      DO J=1,jstop
      DO K=1,kstop
         G3ntrp=G3ntrp+ &
          coef(1,i)*coef(2,j)*coef(3,k)*a(ndx(1)+i-1,ndx(2)+j-1,ndx(3)+k-1)
      END DO
      END DO
      END DO
!
      RETURN
      END FUNCTION G3ntrp
!
!=========================================================================


!=========================================================================
!
!Advance eeta by dt nstep times, dtcpl=dt x nstep

SUBROUTINE Move_plasma_grid_MHD (dt,nstep)
    use rice_housekeeping_module, ONLY : LowLatMHD,doNewCX,ELOSSMETHOD,doFLCLoss,dp_on,doPPRefill,doSmoothDDV,staticR,NowKp
    use math, ONLY : SmoothOpTSC,SmoothOperator33
    use lossutils, ONLY : CXKaiju,FLCRat
    use planethelper, ONLY : DipFTV_colat,DerivDipFTV
    use constants, ONLY : nt,radius_earth_m
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: dt
    INTEGER (iprec), INTENT(IN) :: nstep

    !Clawpack-sized grids
    REAL (rprec), dimension(-1:isize+2,-1:jsize-1) :: etaC
    !Clawpack x Nk sized grids
    REAL (rprec), dimension(-1:isize+2,-1:jsize-1,1:kcsize) :: didt,djdt,rateC
    !RCM-sized grids
    REAL (rprec), dimension( 1:isize  , 1:jsize  ) :: rate,dvedi,dvedj,vv,dvvdi,dvvdj,dvmdi,dvmdj
    REAL (rprec), dimension( 1:isize  , 1:jsize  ) :: vv_avg,dvvdi_avg,dvvdj_avg

    REAL (rprec), dimension( 1:isize  , 1:jsize  ) :: ftv,dftvi,dftvj

    LOGICAL, dimension(1:isize,1:jsize) :: isOpen
    INTEGER, dimension(1:isize,1:jsize) :: ocbDist
    INTEGER (iprec) :: iOCB_j(1:jsize)
    REAL (rprec) :: mass_factor,r_dist,lossCX,lossFLC
    REAL (rprec), dimension(2) :: lossFT
    REAL (rprec), save :: xlower,xupper,ylower,yupper, T1,T2 !Does this need save?
    INTEGER (iprec) :: i, j, kc, ie, iL,jL,iR,jR,iMHD,n
    INTEGER (iprec) :: CLAWiter, joff

    REAL (rprec) :: T1k,T2k !Local loop variables b/c clawpack alters input
    LOGICAL, save :: FirstTime=.true.

    call Tic("Move_Plasma_Init")
    if (jwrap /= 3) then
        write(*,*) 'Somebody should rewrite this code to not assume that jwrap=3'
    stop
    endif

    !Doing silly thing to find i of MHD's lowlat BC
    do i=1,isize
        if (0.5*PI-colat(i,jwrap) <= LowLatMHD) exit
    enddo
    iMHD = i !low-lat boundary for MHD on RCM grid

    !---
    !Do prep work
    where (eeta<0)
        eeta = 0.0
    endwhere

    eeta_avg = 0.0
    joff=jwrap-1

    if (FirstTime) then
        T1=0.
        FirstTime = .false.
    else
        T1=T2
    end if

    T2=T1+dt

    xlower = 1
    xupper = isize
    ylower = 0.0
    yupper = jsize-3
  
    !---
    !Get OCB
    isOpen = (vm < 0)
    do j=1,jsize
        if (any(isOpen(:,j))) then
            !Some open cells on this column
            do i=isize,1,-1
                if (isOpen(i,j)) exit
            enddo
            iOCB_j(j) = i
        else
            !No open cells here
            iOCB_j(j) = 0
        endif
    enddo !j loop

    if (doNoBndFlow) then
        ocbDist = OCBMap(isOpen,nNBFL)
    endif

    !---
    !Calculate node-centered IJ gradients for use inside loop (instead of redoing for each channel)
    !veff = v + vcorot - vpar + vm*alamc(k) = vv + vm*alamc(k)

    !Do array-sized prep work
    !$OMP PARALLEL WORKSHARE if (L_doOMPClaw)
    fac = 1.0E-3*signbe*bir*alpha*beta*dlam*dpsi*ri**2
    vv     = v     + vcorot - vpar !Current potential
    vv_avg = v_avg + vcorot - vpar !Time-averaged potential for plasmasphere

    where (.not. isOpen)
        !Using ftv directly w/ possible intermediate smoothing
        ftv = vm**(-3.0/2)
    elsewhere
        ftv = 0.0
    endwhere
    !$OMP END PARALLEL WORKSHARE

    !Get IJ gradients of potential
    call Grad_IJ(vv    ,isOpen,dvvdi    ,dvvdj    )
    call Grad_IJ(vv_avg,isOpen,dvvdi_avg,dvvdj_avg)

    !Zero out velocities below staticR if necessary
    if (staticR > TINY) then
        where ( rmin <= staticR )
            dvvdi_avg = 0.0
            dvvdj_avg = 0.0
        endwhere
    endif

    !Now get energy-dep. portion, grad_ij vm
    call FTVGrad(ftv,isOpen,dftvi,dftvj)

    !$OMP PARALLEL WORKSHARE if (L_doOMPClaw)
    dvmdi = (-2.0/3.0)*(ftv**(-5.0/3.0))*dftvi
    dvmdj = (-2.0/3.0)*(ftv**(-5.0/3.0))*dftvj
    !Calculate plasmasphere density forall i,j once 
    Dpp = (1.0e-6)*eeta(:,:,1)*dfactor*vm**1.5 !Convert eta to #/cc  
    !$OMP END PARALLEL WORKSHARE

    call Toc("Move_Plasma_Init")
    
!---
!Now calculate things that won't change over the substepping
    call Tic("Move_Plasma_preAdv")
    !ie, di/dj-dt, lossratep/rate
    rate = 0.0
    lossratep =  0.0
    lossmodel = -1.0
    didt = 0.0
    djdt = 0.0

    !$OMP PARALLEL DO if (L_doOMPClaw) &
    !$OMP schedule(dynamic) &
    !$OMP DEFAULT(SHARED) &
    !$OMP private(i,j,kc,ie,iL,jL,iR,jR,rate,dvedi,dvedj) &
    !$OMP private(mass_factor,r_dist,lossCX,lossFLC,lossFT)
    DO kc = kcsize,1,-1
        !Skip boring channels
        IF (.not. advChannel(kc)) CYCLE
        IF (MAXVAL(eeta(:,:,kc)) < machine_tiny) THEN
            eeta(:,:,kc) = 0.0
            advChannel(kc) = .false.
            CYCLE
        ENDIF

        !If oxygen is to be added, must change this!
        IF (alamc(kc) <= 0.0) THEN
            ie = RCMELECTRON
        ELSE
            ie = RCMPROTON
        END IF

        mass_factor = SQRT (xmass(1)/xmass(ie))

  !---
  !Get "interface" velocities on clawpack grid, |-1:isize+2,-1:jsize-1|
        !Start by calculating dvedi,dvedj = grad_ij (veff) = grad_ij (vv) + alamc(k)*grad_ij vm
        if ( (abs(alamc(kc))<TINY) .and. (dtAvg_v>TINY) ) then
            !Do plasmasphere effective potential, uses averaged potential and no energy dep. portion
            dvedi = dvvdi_avg
            dvedj = dvvdj_avg
        else
            !Any other RC channel
            dvedi = dvvdi + alamc(kc)*dvmdi
            dvedj = dvvdj + alamc(kc)*dvmdj
        endif
        !Now loop over clawpack grid interfaces and calculate velocities
        didt(:,:,kc) = 0.0
        djdt(:,:,kc) = 0.0
    
        do j=1,jsize-1 !clawpack jdim
            do i=isize,2,-1

                !I interface

                !Clawpack i,j I-interface is betwen RCM nodes i-1,j+jwrap-1 and i,j+wrap-1
                ! i.e., i,j:I => i-1,j+joff / i,j+joff
                iL = i-1; jL = WrapJ(j+joff)
                iR = i  ; jR = WrapJ(j+joff)

                didt(i,j,kc) = CalcInterface(isOpen(iL,jL),dvedj(iL,jL),fac(iL,jL), &
                                             isOpen(iR,jR),dvedj(iR,jR),fac(iR,jR) )

                !J interface
                !Clawpack i,j J-interface is between RCM nodes i,j+joff-1 and i,j+joff
                iL = i; jL = WrapJ(j+joff-1)
                iR = i; jR = WrapJ(j+joff  )

                !Note extra - in dvedi part of call
                djdt(i,j,kc) = CalcInterface(isOpen(iL,jL),-dvedi(iL,jL),fac(iL,jL), &
                                             isOpen(iR,jR),-dvedi(iR,jR),fac(iR,jR) )

            enddo
        enddo !j loop

        !Freeze flow too close to MHD inner boundary
        didt(iMHD-1:,:,kc) = 0.0 
        djdt(iMHD+1:,:,kc) = 0.0

        !Freeze flow into the domain, only move stuff around from MHD buffer
        didt(1:2,:,kc) = 0.0
        djdt(1  ,:,kc) = 0.0

        !Halt inflow from buffer cells right next to ocb
        if (doNoBndFlow) then
            call NoBoundaryFlow(ocbDist,didt(:,:,kc),djdt(:,:,kc),nNBFL)
        endif

        call PadClaw(didt(:,:,kc))
        call PadClaw(djdt(:,:,kc))

        !---
        !Calculate loss terms on clawpack grid
        !Start w/ loss term on RCM grid
        do j=1,jsize
            do i=1,isize
                !Do some init
                lossCX  = 0.0
                lossFLC = 0.0
                lossFT  = 0.0
                lossratep(i,j,kc) = 0.0
                lossmodel(i,j,kc) = -1.0 ! -1: undefined; 0: C05; 1: chorus; 2: hiss; 3: C+H; 4: strong diffusion; 5: fudge; 10: ion FLC.
                rate(i,j) = 0.0

                if (isOpen(i,j)) then
                    rate(i,j) = -TINY !Set negative value to signal clawpack source term
                    cycle
                endif

                !Calculate losses and keep track of total losses/precip losses
                if ( (ie == RCMELECTRON) .and. (kc /= 1) ) then
                    !Do electron losses
                    lossFT = Ratefn(xmin(i,j),ymin(i,j),alamc(kc),vm(i,j),bmin(i,j),losscone(i,j),Dpp(i,j),NowKp,fudgec(kc),sini(i,j),bir(i,j),mass_factor,ELOSSMETHOD)
                    lossratep(i,j,kc) = lossratep(i,j,kc) + lossFT(1)
                    lossmodel(i,j,kc) = lossFT(2)
                    rate(i,j) = rate(i,j) + lossFT(1)
                endif

                if (ie == RCMPROTON) then
                !Do ion losses
                    r_dist = sqrt(xmin(i,j)**2+ymin(i,j)**2)
                    if ( L_dktime ) then
                        lossCX = CXKaiju(ie,abs(alamc(kc))*vm(i,j),r_dist)
                    endif
                    if (doFLCLoss) then
                        lossFLC = FLCRat(ie,alamc(kc),vm(i,j),bmin(i,j),radcurv(i,j),losscone(i,j))
                    endif
                    lossratep(i,j,kc) = lossratep(i,j,kc) + lossFLC
                    lossmodel(i,j,kc) = 10.0
                    rate(i,j) = rate(i,j) + lossFLC + lossCX
                endif

            enddo !i loop
        enddo !j loop

        !Have loss on RCM grid, now get claw grid
        call rcm2claw(rate,rateC(:,:,kc))
    ENDDO !kc loop
    call Toc("Move_Plasma_preAdv")

!Done static (per coupling) things, now substep and advect
    call Tic("Move_Plasma_Adv")
    !---
    !Main channel loop
    !NOTE: T1k/T2k need to be private b/c they're altered by claw2ez

    !$OMP PARALLEL DO if (L_doOMPClaw) &
    !$OMP schedule(dynamic) &
    !$OMP DEFAULT(SHARED) &
    !$OMP private(i,j,kc,n,T1k,T2k,CLAWiter,etaC)
    DO kc = kcsize,1,-1
        !Skip boring channels
        IF (.not. advChannel(kc)) CYCLE
        IF (MAXVAL(eeta(:,:,kc)) < machine_tiny) THEN
            eeta(:,:,kc) = 0.0
            CYCLE
        ENDIF

        eeta_avg(:,:,kc) = 0.0
        if ( (kc==1) .and. (.not. dp_on)) then  ! We are plasmasphere and we don't want to evolve it
            ! Just set eeta_avg to whatever eta is there, and leave
            !! Note: Regions that have ever been outside of active domain will never have psph again, based on current implementation (2023-08-24)
            eeta_avg(:,:,kc) = eeta(:,:,kc)
            cycle
        else  ! We are hot channel or we are evolving plasmasphere channel
              ! Add current weighted eeta as first contribution
            eeta_avg(:,:,kc) = eeta(:,:,kc)/(nstep+1)
        endif
       
       !Sub-step nstep times
        do n=1,nstep
            !---
            !Tally precipitation losses
            deleeta(:,:,kc) = deleeta(:,:,kc) + eeta(:,:,kc)*(1.0-exp(-lossratep(:,:,kc)*dt))
            call circle(deleeta(:,:,kc))

        !---
        !Do clawpack call
            call rcm2claw(eeta(:,:,kc),etaC)
            !Call clawpack, always as first time
            !Need local copies b/c clawpack alters T1/T2
            T1k = T1
            T2k = T2
            call claw2ez(.true.,T1k,T2k,xlower,xupper,ylower,yupper, &
                         CLAWiter,2,isize-1+1,jsize-3,etaC,didt(:,:,kc),djdt(:,:,kc),rateC(:,:,kc))
        !---
        !Unpack and finish up
            !Copy out
            do j=j1,j2 !jwrap,jsize-1
                do i=1,isize-1
                    if (isOpen(i,j)) then
                        eeta(i,j,kc) = 0.0
                    else
                        eeta(i,j,kc) = max(etaC(i,j-joff),0.0)
                    endif
                enddo
            enddo
            eeta(:,jsize,kc) = eeta(:,jwrap,kc)
            call circle(eeta(:,:,kc))

            if ( (kc==1) .and. dp_on .and. doPPRefill) then
                !refill the plasmasphere  04012020 sbao
                !K: Added kc==1 check 8/11/20
                call Kaiju_Plasmasphere_Refill(eeta(:,:,1), xmin,ymin, aloct, vm, imin_j,dt)
                call circle(eeta(:,:,kc)) !Probably don't need to re-circle
            endif

            eeta_avg(:,:,kc) = eeta_avg(:,:,kc) + eeta(:,:,kc)/(nstep+1)
        enddo !substep loop

    enddo !Main kc loop

    call Toc("Move_Plasma_Adv")

  contains

    !Calculate RCM-node centered gradient of FTV
    subroutine FTVGrad(ftv,isOpen,dftvdi,dftvdj)
      REAL (rprec), dimension(1:isize,1:jsize), intent(IN)  :: ftv
      REAL (rprec), dimension(1:isize,1:jsize), intent(OUT) :: dftvdi,dftvdj
      LOGICAL     , dimension(1:isize,1:jsize), intent(IN)  :: isOpen

      REAL (rprec), dimension(1:isize,1:jsize) :: V0,dV
      REAL (rprec), dimension(1:isize,1:jsize) :: dV0i,dV0j,ddVi,ddVj

      INTEGER (iprec) :: i
      REAL (rprec) :: cl,dcldi,dv0dcl

      !Calculate dipole FTV
      do i=1,isize
        cl = colat(i,jwrap)
        V0(i,:) = DipFTV_colat(cl)
      enddo

      !Now decompose the two contributions
      dV = 0.0
      where (.not. isOpen)
        dV = ftv - V0
      endwhere

    !Take gradients of each
      !Grad of dipole, analytic
      dV0i = 0.0
      dV0j = 0.0
      do i=2,isize-1
        dcldi = 0.5*(colat(i+1,jwrap)-colat(i-1,jwrap))
        cl = colat(i,jwrap)
        dv0dcl = DerivDipFTV(cl)
        dV0i(i,:) = dv0dcl*dcldi
      enddo
      
      !Grad of perturbation
      call Grad_IJ(dV,isOpen,ddVi,ddVj,doLimO=.true. )

      !Possibly smooth grad of perturbation
      if (doSmoothDDV) then
        call Smooth_IJ(ddVi,isOpen)
        call Smooth_IJ(ddVj,isOpen)
      endif

      !Recombine pieces
      dftvdi = dV0i + ddVi
      dftvdj = dV0j + ddVj

      !Old calculation, just do raw gradient
      !call Grad_IJ(ftv,isOpen,dftvdi,dftvdj)

    end subroutine FTVGrad


    !Do smoothing window on RCM grid quantity
    subroutine Smooth_IJ(Q,isOpen)
      REAL (rprec), dimension(1:isize,1:jsize), intent(INOUT)  :: Q
      LOGICAL     , dimension(1:isize,1:jsize), intent(IN)  :: isOpen
      REAL (rprec), dimension(1:isize,1:jsize) :: Qs
      REAL (rprec), dimension(3,3) :: Q33
      LOGICAL     , dimension(3,3) :: G33

      INTEGER (iprec) :: i,j

      Qs = Q

      !$OMP PARALLEL DO if (L_doOMPClaw) &
      !$OMP DEFAULT (SHARED) &
      !$OMP PRIVATE(i,j,Q33,G33)
      do j=j1,j2 !jwrap,jsize-1
        do i=2,isize-1
          Q33(:,:) = Q(i-1:i+1,j-1:j+1)
          G33(:,:) = .not. isOpen(i-1:i+1,j-1:j+1) !Only smooth w/ good cells
          Qs(i,j) = SmoothOperator33(Q33,G33)
        enddo
      enddo
      
      Qs(:,jsize) = Qs(:,jwrap)
      call circle(Qs)
      Q = Qs !Save back smoothed array

    end subroutine Smooth_IJ

    !Copy variable from rcm to clawpack grid
    subroutine rcm2claw(qR,qC)
      REAL (rprec), dimension( 1:isize  , 1:jsize  ), intent(IN)  :: qR
      REAL (rprec), dimension(-1:isize+2,-1:jsize-1), intent(OUT) :: qC

      !Center patch
      qC(1:isize,1:jsize-jwrap) = qR(1:isize,jwrap:jsize-1)
      call PadClaw(qC)

    end subroutine rcm2claw

    !Fill padding of a clawpack grid
    subroutine PadClaw(qC)
      REAL (rprec), dimension(-1:isize+2,-1:jsize-1), intent(INOUT) :: qC
      INTEGER (iprec) :: i, j
      !Pole
      do i=-1,0
        qC(i,j1-joff:j2-joff) = qC(1,j1-joff:j2-joff)
      enddo

      !Equator
      do i=isize+1,isize+2
        qC(i,j1-joff:j2-joff) = qC(isize,j1-joff:j2-joff)
      enddo

      !Periodic
      qC(:,-1:0)                    = qC(:,jsize-4:jsize-3)
      qC(:,jsize-joff:jsize-joff+1) = qC(:,1:2)
    end subroutine PadClaw

    subroutine NoBoundaryFlow(ocbDist,didt,djdt,nL0)
        !Zero flow coming from the direction of the open/closed boundary or cells nL distance from ocb
        !Flow towards the ocb is left alone
        integer, dimension(1:isize,1:jsize), intent(in) :: ocbDist
        real (rp), dimension(-1:isize+2,-1:jsize-1), intent(inout) :: didt,djdt
        integer, intent(in), optional :: nL0  ! Number of i layers inward of o/c boundary to act upon
        
        integer :: nL
        integer :: i,j, jW, jWm1, jWp1

        if(present(nL0)) then
            nL = nL0
        else
            nL = 2
        endif

        !Gonna ignore odd dvdti,j indicies and assume PadClaw is gonna work things out later

        !First, populate open/closed boundary dist map
        ! !! Only strictly a distance for values <= nL, if greater then we don't care
        
        !Now kill some velocities
        !$OMP PARALLEL DO if (L_doOMPClaw) &
        !$OMP schedule(dynamic) &
        !$OMP DEFAULT(SHARED) &
        !$OMP private(i,j,jW,jWm1,jWp1)
        do j=1,jsize-1 !clawpack jdim
            do i=2,isize
                jW = WrapJ(j + joff)
                jWm1 = WrapJ(j + joff - 1)
                jWp1 = WrapJ(j + joff + 1)

                if (ocbDist(i,jW) .eq. 0) then
                    cycle
                else if (ocbDist(i,jW) .le. nL ) then
                    ! Zap
                    if ( (ocbDist(i-1,jW) .lt. ocbDist(i,jW)) .and. (didt(i-1,j) .gt. 0) ) then ! lower i
                        didt(i-1,j) = 0
                    else if ( (ocbDist(i+1,jW) .lt. ocbDist(i,jW)) .and. (didt(i+1,j) .lt. 0) ) then ! upper i
                        didt(i+1,j) = 0
                    else if ( (ocbDist(i,jWm1) .lt. ocbDist(i,jW)) .and. (djdt(i,j-1) .gt. 0) ) then ! lower j
                        djdt(i,j-1) = 0
                    else if ( (ocbDist(i,jWp1) .lt. ocbDist(i,jW)) .and. (djdt(i,j+1) .lt. 0) ) then ! upper j
                        djdt(i,j+1) = 0
                    endif
                !else ! If we aren't next to any layers lower than the one we're currently on, we're done
                !    exit
                endif
            enddo !i loop
        enddo !j loop

    end subroutine NoBoundaryFlow

    function CalcInterface(isOpL,dvL,facL,isOpR,dvR,facR) result(dxdt)
      LOGICAL, intent(IN) :: isOpL,isOpR
      REAL (rprec), intent(IN) :: dvL,dvR,facL,facR
      REAL (rprec) :: dxdt

      REAL (rprec) :: dvAvg,fAvg
      if (isOpL .and. isOpR) then
        !Both sides are bad, no flow
        dxdt = 0.0
        
      else if ( (.not. isOpL) .and. (.not. isOpR) ) then
        !Both sides are good, do basic averaging
        dvAvg = 0.5*( dvL +  dvR)
        fAvg  = 0.5*(facL + facR)
        dxdt = dvAvg/fAvg
      else if (isOpL) then
        !Left = Open / Right = Closed
        !Use right velocity if it's leftward
        dxdt = min(0.0,dvR/facR)
      else if (isOpR) then
        !Left = Closed / Right = Open
        !Use left velocity if it's rightward
        dxdt = max(0.0,dvL/facL)
      endif
    end function CalcInterface


    !Wrap around large indices, jwrap<=>jsize seem to be repeated points on axis
    !So jsize+1 => jwrap+1, i.e. j=>j-jsize+jwrap
    function WrapJ(j) result(jp)
      INTEGER (iprec), intent(IN) :: j
      INTEGER (iprec) :: jp

      if (j>jsize) then
        jp = j-jsize+jwrap
      else
        jp = j
      endif      
    end function WrapJ

    function OCBMap(isOp, nL0) result (ocbDist)
        !Calculate and label cells within nL cells (including diagonal direction) from the ocb
        logical, dimension(1:isize,1:jsize), intent(in) :: isOp
        integer, intent(in), optional :: nL0  ! Number of i layers inward of o/c boundary to act upon
        
        integer :: nL
        integer, dimension(1:isize,1:jsize) :: ocbDist  ! Cell's distance from closest open field cell
        integer :: i,j,iL

        if(present(nL0)) then
            nL = nL0
        else
            nL = 2
        endif

        where (isOp)
            ocbDist = 0
        elsewhere
            ocbDist = nL + 1
        end where

        do iL=1,nL
            !write(*,*) "il=",iL
           do j=2,jsize-1 !clawpack jdim
              do i=2,isize-1
                    if (isOp(i,j)) then
                        cycle
                    else if (any(ocbDist(i-1:i+1,j-1:j+1) .eq. iL-1) .and. (ocbDist(i,j) .eq. nL+1)) then  ! ignore value of current cell
                        !write(*,*) "i,j=",i,j
                        !write(*,*) ocbDist(i-1:i+1,j-1:j+1)
                        !write(*,*) any(ocbDist(i-1:i+1,j-1:j+1) .eq. iL-1)
                        !write(*,*) (ocbDist(i,j) .ne. iL-1)
                        ocbDist(i,j) = iL
                    !!Find a better way to determine when we can quit
                    !else ! If we aren't next to any layers lower than the one we're currently on, we're done
                    !    exit  !This doesn't do anything except try to save time
                            !! Might not want to do this, in case you get some odd open-closed-open lobes
                    endif
                enddo !i loop
            enddo !j loop
        enddo !layer loop

        if (doRCMVerboseH5) then
            last_ocbDist = ocbDist
        endif

        if (.false.) then
            do j=1,jsize
                write(*,*) "j=",j
                write(*,*) ocbDist(:,j)
            enddo
        endif
    end function OCBMap

END SUBROUTINE Move_plasma_grid_MHD

!Calculate RCM-node centered gradient of veff
subroutine Grad_IJ(veff,isOpen,dvedi,dvedj,doLimO)
  REAL (rprec), dimension(1:isize,1:jsize), intent(IN)  :: veff
  REAL (rprec), dimension(1:isize,1:jsize), intent(OUT) :: dvedi,dvedj
  LOGICAL     , dimension(1:isize,1:jsize), intent(IN)  :: isOpen
  LOGICAL, intent(in), optional :: doLimO

  INTEGER (iprec) :: i,j

  LOGICAL :: isOp(3),doLim
  REAL (rprec) :: Q(3)

  if (present(doLimO)) then
    doLim = doLimO
  else
    doLim = .true.
  endif

  dvedi = 0.0
  dvedj = 0.0

  !$OMP PARALLEL DO if (L_doOMPClaw) &
  !$OMP DEFAULT (NONE) &
  !$OMP PRIVATE(i,j,isOp,Q) &
  !$OMP SHARED(dvedi,dvedj,veff,isOpen,doLim)
  do j=2,jsize-1
    do i=2,isize-1
      !Do I deriv
      Q          = veff  (i-1:i+1,j)
      isOp       = isOpen(i-1:i+1,j)
      dvedi(i,j) = Deriv_IJ(Q,isOp,doLim)

      !Do J deriv
      Q          = veff  (i,j-1:j+1)
      isOp       = isOpen(i,j-1:j+1)
      dvedj(i,j) = Deriv_IJ(Q,isOp,doLim)
    enddo
  enddo

!Lower/Upper boundary
  dvedi(1,:) = dvedi(2,:)
  dvedj(1,:) = dvedj(2,:)

  dvedi(isize,:) = dvedi(isize-1,:)
  dvedj(isize,:) = dvedj(isize-1,:)
!Periodic
  dvedi(:,jsize) = dvedi(:,jwrap  )
  dvedi(:,1    ) = dvedi(:,jsize-2)

  dvedj(:,jsize) = dvedj(:,jwrap  )
  dvedj(:,1    ) = dvedj(:,jsize-2)
end subroutine Grad_IJ

!Take derivative from 3-point stencil if possible
function Deriv_IJ(Q,isOp,doLim) result(dvdx)
  LOGICAL     , intent(IN) :: isOp(-1:+1)
  REAL (rprec), intent(IN) ::    Q(-1:+1)
  LOGICAL     , intent(IN)  :: doLim

  LOGICAL, parameter :: doSuperBee = .false. !Use superbee (instead of minmod/MC)
  REAL (rprec) :: dvdx
  REAL (rprec) :: dvL,dvR,dvC
  dvdx = 0.0

  if (isOp(0)) return

  !If still here then central point is closed
  if (isOp(-1) .and. isOp(+1)) then
    !Nothing to work with
    return
  endif
  !Have at least two points
  if ((.not. isOp(-1)) .and. (.not. isOp(+1))) then
    !Both sides closed
    
    !dvdx = 0.5*(Q(+1)-Q(-1)) !Straight up centered derivative
    dvL =       Q( 0) - Q(-1)
    dvR =       Q(+1) - Q( 0)
    dvC = 0.5*( Q(+1) - Q(-1) )

    if (doLim) then
      !Do slope limiter, either minmod or superbee
      if (doSuperBee) then
        !Superbee slope-lim on gradient
        dvdx = qkmaxmod( qkminmod(dvR,2*dvL),qkminmod(2*dvR,dvL) )
      else
        dvdx = MCLim(dvL,dvR,dvC)
        !dvdx = qkminmod(dvL,dvR) !Just minmod lim
      endif
    else
      !Take straight centered difference
      dvdx = dvC
    endif

  else if (.not. isOp(-1)) then
    !-1 is closed, do backward difference
    dvdx = Q(0)-Q(-1)
  else
    !+1 is closed, do forward difference
    dvdx = Q(+1)-Q(0)
  endif

  contains
    function MCLim(dqL,dqR,dqC) result(dqbar)
      REAL (rprec), intent(in) :: dqL,dqR,dqC
      REAL (rprec) :: dqbar
      REAL (rprec) :: magdq

      if (dqL*dqR <= 0) then
        !Sign flip, clamp
        dqbar = 0.0
      else
        !Consistent sense, use MC limiter
        magdq = min(2*abs(dqL),2*abs(dqR),abs(dqC))
        !SIGN(A,B) returns the value of A with the sign of B
        dqbar = sign(magdq,dqC)
      endif
    end function MCLim

    !Quick and lazy minmod limiter
    function qkminmod(a,b) result(c)
      REAL (rprec), intent(in) :: a,b
      REAL (rprec) :: c

      if (a*b > 0) then
        !Pick min modulus
        if (abs(a) < abs(b)) then
          c = a
        else
          c = b
        endif !No sign flip
      else
        c = 0.0
      endif
    end function qkminmod

    !Quick and laxy maxmod limiter
    function qkmaxmod(a,b) result(c)
      REAL (rprec), intent(in) :: a,b
      REAL (rprec) :: c

      if (a*b > 0) then
        !Pick max modulus
        if (abs(a) < abs(b)) then
          c = b
        else
          c = a
        endif !No sign flip
      else
        c = 0.0
      endif
    end function qkmaxmod

end function Deriv_IJ

!Adapted by K: from S. Bao's adaptation of Colby Lemon's code, 09/20

SUBROUTINE Kaiju_Plasmasphere_Refill(eeta0,xmin,ymin,aloct,vm,imin_j,idt)
  use rice_housekeeping_module, ONLY : NowKp
  use earthhelper, ONLY : GallagherXY
  use rcmdefs, ONLY : DenPP0

  implicit none

  REAL (rprec), intent(inout), dimension(isize,jsize) :: eeta0
  REAL (rprec), intent(in), dimension(isize,jsize) :: xmin,ymin, aloct, vm
  REAL (rprec), intent(in)  :: idt
  INTEGER (iprec), intent(in), dimension(jsize) :: imin_j

  integer :: i,j
  REAL (rprec) , parameter :: day2s = 24.0*60.0*60,s2day=1.0/day2s
  REAL (rprec) :: dppT,dpsph,eta2cc,tau,etaT,deta,dndt
  REAL (rprec) :: dpp0,rad,maxX

  dpp0 = 10*DenPP0 !Use 10x the plasmasphere cutoff density to decide on refilling
  maxX = 2.0 !Max over-filling relative to target, i.e. don't go above maxX x den-target


  do j=1,jsize
    do i=1,isize
      if (vm(i,j) <= 0) cycle
      if (i < imin_j(j)+1) cycle !Don't refill outside active domain

      rad = sqrt( xmin(i,j)**2.0 + ymin(i,j)**2.0 )

      !Closed field line, calculate Berbue+ 2005 density (#/cc)
      !Or use Gallagher on nightside w/ NowKp (current Kp)
      !dppT = 10.0**(-0.66*rad + 4.89) !Target refilled density [#/cc]
      dppT = GallagherXY(xmin(i,j),ymin(i,j),NowKp)
      
      eta2cc = (1.0e-6)*dfactor*vm(i,j)**1.5 !Convert eta to #/cc
      dpsph = eta2cc*eeta0(i,j) !Current plasmasphere density [#/cc]

      !Check for other outs before doing anything
      if (dppT  <  dpp0) cycle !Target too low
      !if (dpsph <  dpp0) cycle !Current density too low to bother w/
      if (dpsph >= maxX*dppT) cycle !Too much already there

      etaT = dppT/eta2cc !Target eta for refilling

      !Now calculate refilling
      dndt = 10.0**(3.48-0.331*rad) !cm^-3/day, Denton+ 2012 eqn 1
      !dndt = (cos(aloct(i,j))+1)*dndt !Bias refilling towards dayside

      deta = (idt*s2day)*dndt/eta2cc !Change in eta over idt
      !deta = min(deta,etaT-eeta0(i,j)) !Don't overfill

      eeta0(i,j) = eeta0(i,j) + deta

    enddo
  enddo

END SUBROUTINE Kaiju_Plasmasphere_Refill

!Adapted by S.Bao from Colby Lemon's original code. 04012020 sbao
SUBROUTINE Plasmasphere_Refilling_Model(eeta0, rmin, aloct, vm, idt)

      implicit none
      REAL (rprec), intent(inout), dimension(isize,jsize) :: eeta0
      REAL (rprec), intent(in), dimension(isize,jsize) :: rmin, aloct, vm
      REAL (rprec) :: den_increase(isize, jsize), ftv(isize,jsize)
      REAL (rprec) :: idt
      REAL (rprec) , parameter :: m_per_Re = 6380.e3
      REAL (rprec) , parameter :: nT_per_T = 1.e9
      REAL (rprec) , parameter :: cm_per_m = 1.e2
      where (vm > 0)
        ftv = vm**(-3.0/2.0)
      elsewhere
        ftv = 0.0  ! open field lines, most likely. Set ftv to zero because we want eeta to be zero there
      end where

      den_increase = (idt/1000.0/(24*60*60)) * 10**(3.01 - 0.322*rmin) * ftv * (m_per_Re * nT_per_T * cm_per_m**3)   ! ple/cc
      where (aloct < pi/2 .OR. aloct > 3*pi/2)  ! If we are on the dayside
        eeta0 = eeta0 + 1.8 * den_increase
      elsewhere
        eeta0 = eeta0 + 0.2 * den_increase
      end where


    ! Keep eeta0 between 0 and two times the Berube et al. 2005 density.
      eeta0 = min(eeta0, 2*10**(4.56 - 0.51*rmin) * ftv * (m_per_Re*nT_per_T*cm_per_m**3))
      eeta0 = max(eeta0, 0.0)

END SUBROUTINE

FUNCTION RatefnFDG (fudgx, alamx, sinix, birx, vmx, xmfact)
  IMPLICIT NONE
  REAL (rprec), INTENT (IN) :: fudgx,alamx,sinix,birx,vmx,xmfact
  REAL (rprec)              :: RatefnFDG
  !                                                                       
  !   Function subprogram to compute precipitation rate
  !   Last update:  04-04-88
  !
  RatefnFDG = 0.0466_rprec*fudgx*SQRT(ABS(alamx))*(sinix/birx)*vmx**2
  RatefnFDG = xmfact * RatefnFDG
  RETURN
END FUNCTION RatefnFDG

FUNCTION Ratefn (xx,yy,alamx,vmx,beqx,losscx,nex,kpx,fudgxO,sinixO,birxO,xmfactO,ELOSSMETHOD)

 use lossutils, ONLY : RatefnC_tau_s
 IMPLICIT NONE
 INTEGER (iprec), INTENT (IN) :: ELOSSMETHOD
 REAL (rprec), INTENT (IN) :: xx,yy,alamx,vmx,beqx,losscx,nex,kpx
 REAL (rprec), INTENT (IN), OPTIONAL :: fudgxO,sinixO,birxO,xmfactO
 REAL (rprec) :: fudgx,sinix,birx,xmfact
 REAL (rprec), dimension(2) :: Ratefn
 REAL (rprec) :: L, MLT, K, tau
  if (present(fudgxO)) then
    fudgx = fudgxO
  else
    fudgx = 0.0
  endif
  if (present(sinixO)) then
    sinix = sinixO
  else
    sinix = 0.0
  endif
  if (present(birxO)) then
    birx = birxO
  else
    birx = 0.0
  endif
  if (present(xmfactO)) then
    xmfact = xmfactO
  else
    xmfact = 0.0
  endif

 Ratefn = [1.D-10,-1.D0] ! default rate is 1e-10/s, type is -1.
 select case (ELOSSMETHOD)
         case (ELOSS_FDG)
            Ratefn(1)= RatefnFDG(fudgx, alamx, sinix, birx, vmx, xmfact) !1/s
            Ratefn(2)= -2.0 
         case (ELOSS_SS)
            tau = RatefnC_tau_s(alamx,vmx,beqx,losscx)
            Ratefn(1) = 1.D0/tau !/s
            Ratefn(2) = -1.0
         case (ELOSS_WM)
            if (EWMTauInput%useWM) then
                if (kpx > 10.0) then
                    write(*,*) "Kp = ", kpx, ", invalid Kp input for the wave models. Please use 'FDG' or 'SS' in the electron loss model instead."
                    stop
                endif
                Ratefn = RatefnWM(xx,yy,alamx,vmx,nex,kpx,beqx,losscx)
            else
                write(*,*) "Wave models are missing in rcmconfig.h5"
                stop
            endif
         case default
            write(*,*) "The electron loss rate model type entered is not supported."
            stop
 end select

END FUNCTION Ratefn

FUNCTION RatefnWM(xx,yy,alamx,vmx,nex,kpx,beqx,losscx)
        use lossutils, ONLY : LossR8_PSHEET, LossR8_IMAG, Lo, Li, SSCATTER
        IMPLICIT NONE

        REAL (rprec), INTENT (IN) :: xx,yy,alamx,vmx,nex,kpx,beqx,losscx
        REAL (rprec), dimension(2) :: RatefnWM,rIMAG
        REAL (rprec) :: L,rPS,wIMAG

        L = sqrt(xx**2+yy**2)

        !Wave type: Hiss: 1.0; Chorus: 2.0 (Inner Mag,L<7); strong scattering: 3.0 (Plasmasheet,L>8); Blending zone, 7<L<8
        if (L > Lo) then
            !Only plasma sheet
            rPS = LossR8_PSHEET(alamx,vmx,beqx,losscx)
            RatefnWM(1) = rPS
            RatefnWM(2) = SSCATTER !wave type number for strong scattering
        else if (L < Li) then
            !IMAG only
            rIMAG = LossR8_IMAG(xx,yy,alamx,vmx,nex,kpx)
            RatefnWM(1) = rIMAG(1)
            RatefnWM(2) = rIMAG(2)

        else
            !Middle
            rIMAG = LossR8_IMAG  (xx,yy,alamx,vmx,nex,kpx)
            rPS   = LossR8_PSHEET(alamx,vmx,beqx,losscx)
            wIMAG = RampDown(L,Li,Lo-Li) !Ramp from 1 to 0
            RatefnWM(1) = wIMAG*rIMAG(1) + (1-wIMAG)*rPS
            RatefnWM(2) = SSCATTER-wIMAG
        endif

END FUNCTION RatefnWM

!-------------------------------------
    FUNCTION Gntrp_2d_ang (array, bi, bj, ikind)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: ikind
    REAL (rprec), INTENT (IN)    :: array (:,:), bi, bj
    REAL (rprec)                 :: Gntrp_2d_ang
!
    INTEGER (iprec) :: ii, jj, ni, nj
    REAL (rprec)    :: fi,fj,a1,a2,v_1, v_2
!
!
!   Prepare indices for interpolation:
!
    ni = SIZE (array, 1)
    nj = SIZE (array, 2)
    IF (bj < 1.0 .OR. bj > nj) write(*,*) 'BJ out of range in gntrp_2d_ang'
    IF (bi < 1.0 .OR. bi > ni) write(*,*) 'BI out of range in gntrp_2d_ang'
!
    ii = INT (bi)
    fi = REAL (ii,rprec)
    jj = INT (bj)
    fj = REAL (jj,rprec)
    IF (ii == ni) ii = ii - 1
    IF (jj == nj) jj = jwrap   ! periodicity in J
!                                                                       
    v_1 = array (ii,jj)
    v_2 = array (ii+1,jj)
!   IF (ii < imin_j(jj))   v_1 = array(imin_j(jj),jj)
!   IF (ii+1 < imin_j(jj)) v_2 = array(imin_j(jj),jj)
    a1 = (1.0-(bi-fi))*v_1 + (bi-fi)*v_2
!
    v_1 = array (ii  , jj+1)
    v_2 = array (ii+1, jj+1)
!   IF (ii   < imin_j(jj+1)) v_1 = array (imin_j(jj+1),jj+1)
!   IF (ii+1 < imin_j(jj+1)) v_2 = array (imin_j(jj+1),jj+1)
    a2 = (1.0-(bi-fi))*v_1 + (bi-fi)*v_2
!                                                                       
    IF (ikind == 1) THEN
       IF (jj+1 == nj) a2 = a2 + pi_two
       IF (jj  == jwrap)      a1 = zero
       IF (jj+1 == jwrap)      a2 = a2 + pi_two      ! sts, feb 22
    END IF
!
    Gntrp_2d_ang = (1.0 - (bj-fj)) * a1 + (bj-fj) * a2
!
    RETURN
    END FUNCTION Gntrp_2d_ang
!
!
!
!
!
    FUNCTION Intrp_2d_grid (array, bi, bbj, jwrap, imin_j)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: jwrap, imin_j(:)
    REAL (rprec), INTENT (IN)    :: array (:,:), bi, bbj
    REAL (rprec)                 :: Intrp_2d_grid
!
    INTEGER (iprec) :: ii, jj, imax_array, jmax_array
    REAL (rprec)    :: bj, ca, cb, cc, cd, di, dj
!
!
    imax_array = SIZE (array, 1)
    jmax_array = SIZE (array, 2)
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    bj  = Bjmod ( bbj,jwrap,jmax_array)
    jj  = INT (bj)
!
!   (i,j)---------------------(i,j+1)
!   | A          |               B  |
!   |            |                  |
!   |            |di                |
!   |   dj       |       1-dj       |
!   |---------(bi,bj)---------------|
!   |            |                  |
!   |            |                  |
!   |            |1-di              |
!   |            |                  |
!   | C          |              D   |
!   (i+1,j)-----------------(i+1,j+1)
!
    di = bi - ii
    dj = bj - jj
    IF (ii < imin_j(jj)) THEN
       ca = 0.0
    ELSE
       ca = (one-di)*(one-dj) 
    END IF
    IF (ii+1 < imin_j(jj)) THEN
       cc = 0.0
    ELSE
       cc = di * (one-dj)
    END IF
    IF (ii < imin_j(jj+1)) THEN
       cb = 0.0
    ELSE
       cb = (one-di)*dj
    END IF
    IF (ii+1 < imin_j(jj+1)) THEN
       cd = 0.0
    ELSE
       cd = di*dj
    END IF
!                                                                       
    Intrp_2d_grid = ca*array(ii,jj)+cb*array(ii,jj+1)+ &
                      cc*array(ii+1,jj)+cd*array(ii+1,jj+1)
    Intrp_2d_grid = Intrp_2d_grid / (ca+cb+cc+cd)
!
    RETURN
    END FUNCTION Intrp_2d_grid
    FUNCTION Interp_1d (array, bi )
    IMPLICIT NONE
    REAL (rprec), INTENT (IN)    :: array (:), bi
    REAL (rprec)                 :: Interp_1d
!                                                                       
!   This function subprogram interpolates 1-dim. ARRAY to return
!   the value of ARRAY at the non-integer point BI.
!
!   Stanislav: if Bi < 1, then the array is extrapolated
!              linearly based on the values A(1,:) and A(2,:).
!              If Bi > imax, then array is linearly
!              extrapolated on the values A(imax-1,:) and A(imax,:).
!
    INTEGER (iprec) :: ii, imax_array
    REAL (rprec)    :: fi
!
    imax_array = SIZE (array)
!
    IF (bi < one .OR. bi > imax_array) STOP 'OUT OF BOUNDS IN INTERP_1D'
!
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
    Interp_1d = (one - (bi-fi) ) * array (ii) + (bi-fi) * array (ii+1)
    RETURN
    END FUNCTION Interp_1d



    FUNCTION Interp_2d (array, bi, bj)
    IMPLICIT NONE
    REAL (rprec), INTENT (IN)    :: array (:,:), bi, bj
    REAL (rprec)                 :: Interp_2d
!
    INTEGER (iprec) :: ii, jn, jj, jp1, imax_array, jmax_array
    REAL (rprec)    :: fi,fj,a1,a2
!
!
!   Prepare indices for interpolation:
!
    imax_array = SIZE (array, 1)
    jmax_array = SIZE (array, 2)
!
!
    IF (bi < one .OR. bi > imax_array .OR. bj < 1 .OR. bj > jmax_array) &
        STOP 'OUT OF BOUNDS IN INTERP_2D'
!
!
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
!   Decide which interpolation to perform and proceed:
!
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp. of 2-d array
!
      Interp_2d = (one-(bi-fi)) * array(ii,jn) + (bi-fi)*array(ii+1,jn)
!
    ELSE    !        2-d interpolation of 2-d array:
!
!         If jwrap <= bj < jmax, then jwrap-1 <= INT(bj) <= jmax-1
!         and jwrap <= INT(bj)+1 <= jmax
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (one-(bi-fi))*array(ii,jj)  + (bi-fi)*array(ii+1,jj)
       a2 = (one-(bi-fi))*array(ii,jp1) + (bi-fi)*array(ii+1,jp1)
!                                                                       
       Interp_2d = (one - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Interp_2d




    FUNCTION Interp_2d_of3d (array, bi, bj, index_3)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: index_3
    REAL (rprec), INTENT (IN)    :: array (:,:,:), bi, bj
    REAL (rprec)                 :: Interp_2d_of3d
!                                                                       
!   This is the same as Gntrp_2d but for a 3-dim array, see comments for
!   Gntrp_2d. A separate function is needed since if Gntrp_2d were used,
!   then we would need to pass array sections (the other option is Fortran
!   77 style of passing an offset array, but that should be avoided for
!   compiler checking and parallelization reasons).
!
    INTEGER (iprec) :: ii, jn, jj, jp1, imax_array, jmax_array, kmax_array
    REAL (rprec)    :: fi,fj,a1,a2
!
!
!   Prepare indices for interpolation:
!
    imax_array = SIZE (array, 1)
    jmax_array = SIZE (array, 2)
    kmax_array = SIZE (array, DIM = 3)
!
!
    IF (bi < one .OR. bi > imax_array .OR. bj < one .OR. bj > jmax_array) THEN
       WRITE (*,*) 'OUT OF BOUNDS IN INTERP_2D_OF_3D'
       STOP
    END IF
!
!
    IF (index_3 > kmax_array .OR. index_3 < 1) STOP 'INTRP_2D_OF3D: index_3 OUT OF RANGE'
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
!   Decide which interpolation to perform and proceed:
!
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp. of 2-d array
!
      Interp_2d_of3d = (one-(bi-fi)) * array(ii,jn,index_3) + &
                  (bi-fi)*array(ii+1,jn,index_3)
!
    ELSE    !        2-d interpolation of 2-d array:
!
!         If jwrap <= bj < jmax, then jwrap-1 <= INT(bj) <= jmax-1
!         and jwrap <= INT(bj)+1 <= jmax
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (one-(bi-fi))*array(ii,jj,index_3)  + (bi-fi)*array(ii+1,jj,index_3)
       a2 = (one-(bi-fi))*array(ii,jp1,index_3) + (bi-fi)*array(ii+1,jp1,index_3)
!                                                                       
       Interp_2d_of3d = (one - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Interp_2d_of3d
!
!   7. A few routines needed for circulariation and normalization in J:
!
    FUNCTION Bjmod_real (bj, jwrap, jsize)
    IMPLICIT NONE
    REAL (rprec),    INTENT (IN) :: bj
    INTEGER(iprec),  INTENT (IN) :: jwrap, jsize
    REAL (rprec)                 :: Bjmod_real
!_____________________________________________________________________________
!   last update: 11-28-84               by:rws
!                                                                       
!   this function subporgram returns bjmod with a value
!   between jwrap and jmax-1. In RCM, arrays in j (local time angle)
!   are dimensioned from 1 to jsize, but the grid wraps around and
!   overlaps such that array (jwrap) = array (jsize)
!   array (jwrap-1) = array (jsize-1), etc. In other words, only
!   elements from j=jwrap to j=jsize-1 are unique. This function takes
!   a non-integer j index, BJ, and makes sure that it is larger or
!   equal to jwrap but smaller than jsize-1. Then when array is interpolated
!   on two elements, j is never larger than jsize or smaller than jwrap-1.
!   For the case of jwrap = 3 and a 1-dim array, this looks like:
!
!   j-value:   1  2  3  4  5              jsize-2   jsize-1   jsize
!              x  x  x  x  x .................x        x     x
!              |  |  |                        |        |     |
!              |  |   --------->---------->---|--------|-----
!              |  --------------->------------|-->-----
!               ---------->----------->-------
!
!   Dependency:  none
!
    Bjmod_real = bj
!                                                                       
!   do_1: DO
!      IF (Bjmod_real < REAL (jsize - 1,rprec)) EXIT
!      Bjmod_real = Bjmod_real - REAL (jsize - jwrap,rprec)
!   END DO do_1
!
!   do_2: DO
!      IF (Bjmod_real >= REAL (jwrap,rprec)) EXIT
!      Bjmod_real = Bjmod_real + REAL(jsize - jwrap,rprec)
!   END DO do_2
!
bjmod_real = MODULO(bj-REAL(jwrap),REAL(jsize-jwrap-1)) + REAL(jwrap) 
    RETURN
    END FUNCTION Bjmod_real
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    FUNCTION Bjmod_int (bj, jwrap, jsize)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: bj
    INTEGER (iprec), INTENT (IN) :: jwrap, jsize
    INTEGER (iprec)              :: Bjmod_int
!_____________________________________________________________________________
!   last update: 11-28-84               by:rws
!                                                                       
!   this function subporgram returns bjmod with a value
!   between jwrap and jmax-1. In RCM, arrays in j (local time angle)
!   are dimensioned from 1 to jsize, but the grid wraps around and
!   overlaps such that array (jwrap) = array (jsize)
!   array (jwrap-1) = array (jsize-1), etc. In other words, only
!   elements from j=jwrap to j=jsize-1 are unique. This function takes
!   a non-integer j index, BJ, and makes sure that it is larger or
!   equal to jwrap but smaller than jsize-1. Then when array is interpolated
!   on two elements, j is never larger than jsize or smaller than jwrap-1.
!   For the case of jwrap = 3 and a 1-dim array, this looks like:
!
!   j-value:   1  2  3  4  5              jsize-2   jsize-1   jsize
!              x  x  x  x  x .................x        x     x
!              |  |  |                        |        |     |
!              |  |   --------->---------->---|--------|-----
!              |  --------------->------------|-->-----
!               ---------->----------->-------
!
!   Dependency:  none
!
    Bjmod_int = bj
!                                                                       
    do_1: DO
       IF (Bjmod_int > jsize - 1) THEN
          Bjmod_int = Bjmod_int - (jsize - jwrap)
       ELSE
          EXIT do_1
       END IF
    END DO do_1
!
    do_2: DO
       IF (Bjmod_int < jwrap) THEN
           Bjmod_int = Bjmod_int + (jsize - jwrap)
       ELSE
           EXIT do_2
       END IF
    END DO do_2
!
    RETURN
    END FUNCTION Bjmod_int
!
!
    END MODULE Rcm_mod_subs
