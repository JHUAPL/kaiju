! module to convert imag data to mhd data
module imag2mix_interface
    use mixdefs
    use mixgeom
    use volttypes
    use cmiutils
    use planethelper
    use gridutils
    use shellInterp

    implicit none

    contains

   subroutine CoupleIMagToMix(vApp)
      class(voltApp_T), intent(inout) :: vApp
      integer :: i,j
      real(rp) :: x1,x2,x1c,x2c
      logical :: isTasty
      real(rp) :: imP(nVars_imag2mix)
      integer, dimension(9) :: varNames = (/IM_EAVG,IM_ENFLX,IM_EFLUX,IM_GTYPE,IM_EDEN,IM_EPRE,IM_NPSP,IM_THCON,IM_PHCON/)

      ! Zero out before conditional updates
      do i=1,size(varNames)
         vApp%remixApp%ion(NORTH)%St%Vars(:,:,varNames(i)) = 0.0
         vApp%remixApp%ion(SOUTH)%St%Vars(:,:,varNames(i)) = 0.0
      enddo
      ! G%t(G%Np,G%Nt), G%p = modulo((atan2(G%y,G%x)+2*pi),(2*pi)) 
      ! G%x = sin(G%t)*cos(G%p)
      ! G%y = sin(G%t)*sin(G%p)
      associate(G=>vApp%remixApp%ion(NORTH)%G, ion=>vApp%remixApp%ion)
         !Loop over and get imag data
         !$OMP PARALLEL DO default(shared) collapse(2) &
         !$OMP private(i,j,x1,x2,x1c,x2c,imP,isTasty)
            do j=1,G%Np
               do i=1,G%Nt
                  ! NH mapping
                  x1 = G%t(j,i) ! G%t is colatitude.
                  x2 = G%p(j,i) ! need to be [0,2pi]
                  imP = 0.0_rp
                  ! Here we use the global voltron grid to first find the conjugate point in the SH.
                  call InterpShellVar_TSC_pnt(vApp%shGrid, vApp%State%tubeShell%latc, x1, x2, x1c)
                  call InterpShellVar_TSC_pnt(vApp%shGrid, vApp%State%tubeShell%lonc, x1, x2, x2c)
                  ion(NORTH)%St%Vars(j,i,IM_THCON) = PI/2.0_rp - x1c ! conjugate co-lat in radians, 0-pi ! State%thcon(i0,j0)
                  ion(NORTH)%St%Vars(j,i,IM_PHCON) = x2c ! conjugate long in radians, 0-2pi ! State%phcon(i0,j0)
                  !call vApp%imagApp%getMomentsPrecip(x1,x2,imP,isTasty)
                  call superSampleMomentsPrecip(G, i, j, 5, imP, isTasty)
                  ! gtype is 1 when isTasty is true, but can be 0 or 0.5 when isTasty is false.
                  ion(NORTH)%St%Vars(j,i,IM_GTYPE) = imP(RAI_GTYPE)      ! [0~1]
                  if (isTasty) then
                     ion(NORTH)%St%Vars(j,i,IM_EAVG ) = imP(RAI_EAVG )      ! [keV]
                     ion(NORTH)%St%Vars(j,i,IM_ENFLX) = imP(RAI_ENFLX)      ! [#/cm^2/s]
                     ion(NORTH)%St%Vars(j,i,IM_EFLUX) = imP(RAI_EFLUX)      ! [ergs/cm^2/s]
                     ion(NORTH)%St%Vars(j,i,IM_EDEN ) = imP(RAI_EDEN )      ! [#/m^3]
                     ion(NORTH)%St%Vars(j,i,IM_EPRE ) = imP(RAI_EPRE )      ! [Pa]
                     ion(NORTH)%St%Vars(j,i,IM_NPSP ) = imP(RAI_NPSP )      ! [#/m^3]
                  endif

                  ! SH mapping
                  ! Earlier we used a simple mirro mapping, i.e.:
                  ! ion(SOUTH)%St%Vars(j,i,IM_EAVG ) = ion(NORTH)%St%Vars(G%Np:1:-1,j,i,IM_EAVG )
                  x1 = PI-G%t(G%Np+1-j,i)     ! G%t is colatitude. SH colat should be 0.75pi-pi if remix llb is 0.25pi (45 deg).
                  x2 = G%p(G%Np+1-j,i)        ! need to be [0,2pi]
                  imP = 0.0_rp
                  ! Here we use the global voltron grid to first find the conjugate point in the NH.
                  call InterpShellVar_TSC_pnt(vApp%shGrid, vApp%State%tubeShell%latc, x1, x2, x1c)
                  call InterpShellVar_TSC_pnt(vApp%shGrid, vApp%State%tubeShell%lonc, x1, x2, x2c)
                  ion(SOUTH)%St%Vars(j,i,IM_THCON) = PI/2.0_rp - x1c ! conjugate co-lat in radians, 0-pi ! State%thcon(i0,j0)
                  ion(SOUTH)%St%Vars(j,i,IM_PHCON) = x2c ! conjugate long in radians, 0-2pi ! State%phcon(i0,j0)

                  ! Note x1c is the conjugate lat in the NH (rather than colat). 
                  ! x1c is returned as 0 for open field lines. 
                  ! Only interpolate if x1c is within the remix grid.
                  if (x1c>=PI/2.0_rp - G%t(1,G%Nt)) then
                     !call vApp%imagApp%getMomentsPrecip(PI/2.0_rp-x1c,x2c,imP,isTasty)
                     call superSampleMomentsPrecip(G, i, j, 5, imP, isTasty, thetaConO=PI/2.0_rp-x1c, phiConO=x2c)
                     ! gtype is 1 when isTasty is true, but can be 0 or 0.5 when isTasty is false.
                     ion(SOUTH)%St%Vars(j,i,IM_GTYPE) = imP(RAI_GTYPE)      ! [0~1]
                     if (isTasty) then
                        ion(SOUTH)%St%Vars(j,i,IM_EAVG ) = imP(RAI_EAVG )      ! [keV]
                        ion(SOUTH)%St%Vars(j,i,IM_ENFLX) = imP(RAI_ENFLX)      ! [#/cm^2/s]
                        ion(SOUTH)%St%Vars(j,i,IM_EFLUX) = imP(RAI_EFLUX)      ! [ergs/cm^2/s]
                        ion(SOUTH)%St%Vars(j,i,IM_EDEN ) = imP(RAI_EDEN )      ! [#/m^3]
                        ion(SOUTH)%St%Vars(j,i,IM_EPRE ) = imP(RAI_EPRE )      ! [Pa]
                        ion(SOUTH)%St%Vars(j,i,IM_NPSP ) = imP(RAI_NPSP )      ! [#/m^3]
                     endif
                  endif
               enddo !i
            enddo !j
      end associate


      contains

      subroutine superSampleMomentsPrecip(G, i, j, nPnts, imP_avg, isTasty, thetaConO, phiConO)
         !! Super sample a remix cell to account for disparate resolutions between remix and imag model
         !! Should ultimately be replaced with a flux-conserving interpolation method eventually
         !! Note: For conjugate points, we are assuming the same cell size as the remix point's origin, which can have some issues
         !!   In order to do better cell area, need to calculate all conjugate points first and then pass all neighbors, do fancy stuff
         type(mixGrid_T), intent(in) :: G
         !class(voltApp_T), intent(in) :: vApp  ! Some sort of intent complaint, now we just use it from parent scope
         integer, intent(in) :: i, j
         integer, intent(in) :: nPnts
         real(rp), intent(inout) :: imP_avg(nVars_imag2mix)
         logical, intent(inout) :: isTasty
         real(rp), intent(in), optional :: thetaConO, phiConO
            !! If wanting conjugate points, provide them here. Otherwise will use G%t and G%p at point i,j

         logical :: doSouth
         integer :: n, nGood
         integer :: iPnt, jPnt
         logical :: isTastyPnt
         real(rp) :: thetaIn, phiIn
            !! Theta and phi vals we are evaluating (either G%t or thetaConO, and same for phi)
         real(rp) :: thMin, thMax, phMin, phMax
         real(rp), dimension(:), allocatable :: theta_subdiv, phi_subdiv
         real(rp) :: imP_pnt(nVars_imag2mix)

         imP_avg = 0.0
         nGood = 0
         if (nPnts < 2) then
            write(*,*)"ERROR in imag2mix_interface.F90:superSampleMomentsPrecip"
            write(*,*)"nPnts must be at least 2, preferrably higher"
            stop
         endif

         if (present(thetaConO) .and. present(phiConO)) then
            thetaIn = thetaConO
            phiIn = phiConO
         else
            thetaIn = G%t(j,i)
            phiIn = G%p(j,i)
         endif

         allocate(theta_subdiv(nPnts))
         allocate(phi_subdiv  (nPnts))

         ! Calculate remix cell bounds for index i,j
         ! Lazy handling of i bounds
         if (i == 1 .or. i == G%Nt) then
            theta_subdiv = thetaIn
         else
            thMin = thetaIn - 0.5_rp*G%dt(j,i-1)
            thMax = thetaIn + 0.5_rp*G%dt(j,i  )
            theta_subdiv = genThetaSubdiv(thMin, thMax, nPnts)
         endif
         ! Gen phi vals
         if (j == 1) then
            phMin = phiIn - 0.5_rp*G%dp(G%Np,i)
         else
            phMin = phiIn - 0.5_rp*G%dp(j-1,i)
         endif
         phMax = phiIn + 0.5_rp*G%dp(j  ,i)
         do n=1,nPnts
            phi_subdiv = phMin + (n-1)*1.0_rp/(nPnts-1)*(phMax - phMin)  ! Don't forget to cast something to real so we don't do integer division
         enddo

         ! Now that we have our points we can start getting our values
         do iPnt=1,nPnts
            do jPnt=1,nPnts
               call vApp%imagApp%getMomentsPrecip(theta_subdiv(iPnt), phi_subdiv(jPnt), imP_pnt, isTastyPnt)
               if (.not. isTastyPnt) then
                  cycle
               endif
               ! Otherwise we got good data
               nGood = nGood + 1
               imP_avg = imP_avg + imP_pnt
            enddo
         enddo
         ! Finish average and make final isTasty decision
         if (nGood > 0) then
            isTasty = .true.

            ! For some variables, don't divide by nGood, divide by # of all points sampled
            ! e.g. in case where 1/nPnts^2 were good it gets applied to whole cell when really the rest are telling you there should be no precip there
            ! So you actually to wanna include some effect from bad points
            !IM_EAVG,IM_ENFLX,IM_EFLUX,IM_GTYPE,IM_EDEN,IM_EPRE,IM_NPSP,IM_THCON,IM_PHCON
            imP_avg(RAI_ENFLX) = imP_avg(RAI_ENFLX)/nPnts**2
            imP_avg(RAI_EFLUX) = imP_avg(RAI_EFLUX)/nPnts**2
            imP_avg(RAI_EDEN)  = imP_avg(RAI_EDEN )/nPnts**2
            imP_avg(RAI_EPRE)  = imP_avg(RAI_EPRE )/nPnts**2
            imP_avg(RAI_NPSP)  = imP_avg(RAI_NPSP )/nPnts**2
            imP_avg(RAI_EAVG)  = imP_avg(RAI_EFLUX) / max(TINY, imP_avg(RAI_ENFLX))
            imP_avg(RAI_GTYPE) = imP_avg(RAI_GTYPE)/nGood
            imP_avg(RAI_THCON) = imP_avg(RAI_THCON)/nGood
            imP_avg(RAI_PHCON) = imP_avg(RAI_PHCON)/nGood
            
         else
            isTasty = .false.
         endif

      end subroutine superSampleMomentsPrecip

      function genThetaSubdiv(theta0, theta1, nPnts) result(theta_subdiv)
         !! Return theta array between theta0 and theta1 where areas of cells will be equal
         !! Expecting Theta0 < Theta1
         real(rp), intent(in) :: theta0, theta1
         integer, intent(in) :: nPnts

         integer :: i
         real(rp) :: dCos, arg
         real(rp), dimension(nPnts) :: theta_subdiv

         dCos = cos(theta0) - cos(theta1)

         theta_subdiv(1) = theta0
         theta_subdiv(nPnts) = theta1
         do i=2,nPnts-1
            arg = cos(theta0) - (i-1)*1.0_rp/(nPnts-1)*dCos  ! Don't forget to cast something to real so we don't do integer division
            theta_subdiv(i) = acos(arg)
         enddo

      end function genThetaSubdiv

   end subroutine CoupleIMagToMix

end module

