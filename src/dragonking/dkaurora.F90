module dkaurora
  use kdefs
  use mixdefs
  use mixtypes
  use dkdefs
  use dktypes
!  use planethelper
!  use gcmtypes
!  use gcminterp
  use math
!  use euvhelper
  use auroralhelper
!  use kai2geo
  use raijudefs, ONLY : def_tiote
  
  implicit none

  real(rp), dimension(:,:), allocatable, private :: beta_RCM,alpha_RCM,gtype_RCM ! two-dimensional beta based on RCM fluxes.

  !Replacing some hard-coded inline values (bad) w/ module private values (slightly less bad)
  real(rp), parameter, private :: maxDrop = 20.0 !Hard-coded max potential drop [kV]
  real(rp), parameter, private :: eTINY = mixeTINY ! Floor of average energy [keV]
  real(rp), parameter, private :: Ne_floor = 0.03e6 ! minimum Ne in [/m^3] when evaluating the linearized FL relation.
  real(rp), parameter, private :: Ne_psp = 10.0e6 ! Ne threshold for the plasmasphere in [/m^3].
  real(rp), parameter, private :: GuABNF = 1.e7 ! Gussenhoven+[1983] Auroral Boundary Number Flux in [#/cm^s/s].
  real(rp), private :: RinMHD = 0.0 !Rin of MHD grid (0 if not running w/ MHD)
  real(rp), private :: MIXgamma
  real(rp), private :: beta_inp

  contains

    subroutine dragonking_init(aurora,Params,G)
      ! Initialize precipitation variables. Called in mixmain.F90/subroutine init_mix.
      type(dkParams_T) , intent(inout) :: aurora
      type(mixParams_T), intent(in)    :: Params
      type(mixGrid_T)  , intent(in)    :: G

      ! define these module-wide variables so we don't have to pass the Params object to all the aurora functions
      ! this is similar to how it was done in MIX
      aurora%alpha             = Params%alpha
      aurora%beta              = Params%beta
      aurora%R                 = Params%R
      aurora%doRamp            = Params%doRamp
      aurora%doChill           = Params%doChill
      aurora%doMR              = Params%doMR
      aurora%doAuroralSmooth   = Params%doAuroralSmooth
      aurora%aurora_model_type = Params%aurora_model_type

      if (.not. allocated(aurora%rampFactor)) allocate(aurora%rampFactor(G%Np,G%Nt))
      if (.not. allocated(aurora%deltaE)) allocate(aurora%deltaE(G%Np,G%Nt))
      if (.not. allocated(aurora%E0)) allocate(aurora%E0(G%Np,G%Nt))
      if (.not. allocated(aurora%phi0)) allocate(aurora%phi0(G%Np,G%Nt))
      if (.not. allocated(aurora%PrecipMask)) allocate(aurora%PrecipMask(G%Np,G%Nt))

      ! these arrays are global and should not be! reallocate them
      if(allocated(beta_RCM)) deallocate(beta_RCM)
      if(allocated(alpha_RCM)) deallocate(alpha_RCM)
      if(allocated(gtype_RCM)) deallocate(gtype_RCM)

      allocate(beta_RCM(G%Np,G%Nt))
      allocate(alpha_RCM(G%Np,G%Nt))
      allocate(gtype_RCM(G%Np,G%Nt))

      call SetMIXgamma(Params%gamma)
      ! MHD inner boundary, used to calc mirror ratio.
      RinMHD = Params%RinMHD
      ! Te/Tmhd
      alpha_RCM = 1.0/(def_tiote+1.0)
      ! Loss cone rate
      beta_RCM  = aurora%beta
      beta_inp  = aurora%beta
      ! RCM grid weight: 1. Totally on closed RCM; 0. Totally outside RCM.
      ! if conductance_IM_GTYPE is not called, gtype_RCM has all zero, MHD values have a weight of 1.
      gtype_RCM = 0.0

    end subroutine dragonking_init

    subroutine dragonking_total(aurora,G,St,conductance)
      ! The main subroutine in module dragonking.
      ! It derives auroral precipitation.
      type(dkParams_T), intent(inout) :: aurora
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      ! Compute EUV though because it's used in fedder
!      call conductance_euv(conductance,G,St) ! moved to mixmain

      select case ( aurora%aurora_model_type )
        case (FEDDER)
          call dragonking_fedder95(aurora,G,St,conductance)
        case (SUNNY)
          call dragonking_sunny(aurora,G,St)
        case (LINMRG)
          call dragonking_linmrg(aurora,G,St)
        case default
          stop "The auroral precipitation model type entered is not supported."
      end select

      ! Correct for multiple reflections if you're so inclined
      if (aurora%doMR) call dragonking_mr(aurora,G,St)

      ! Smooth precipitation energy and flux before calculating conductance.
      if (aurora%doAuroralSmooth) call dragonking_auroralsmooth(St,G,aurora)

    end subroutine dragonking_total

    subroutine dragonking_fedder95(aurora,G,St,conductance)
      ! Derive electron precipitation energy flux and avg energy using Fedder95 formula.
      type(dkParams_T), intent(inout) :: aurora
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      type(mixConductance_T), intent(inout) :: conductance
      
      real(rp) :: signOfJ
      real(rp) :: Rout = 6.D0, Rin = 1.2D0
      real(rp) :: rhoFactor = 3.3D-24*0.5D0
      real(rp) :: rPolarBound = 20.0D0*pi/180.D0, rEquatBound = 30.0D0*pi/180.0D0, rLowLimit = 0.02D0
      real(rp) :: D,Cs,Pe,Ne,dE,phi0,E0,aRes,dEc
      integer :: i,j

      if (St%hemisphere==NORTH) then
         signOfJ = -1  
      elseif (St%hemisphere==SOUTH) then
         signOfJ = 1
      else
         stop 'Wrong hemisphere label. Stopping...'
      endif

      ! fills in rampFactor
      if (aurora%doRamp) then
         call dragonking_ramp(aurora,G,rPolarBound,rEquatBound,rLowLimit)
      else 
         aurora%rampFactor = 1.0D0
      endif
      
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,D,Cs,dE,phi0,E0,aRes,dEc)
      do j=1,G%Nt
         do i=1,G%Np
            if (aurora%doChill) then
               ! MHD density replaced with gallagher where it's lower      
               ! and temperature changed correspondingly
               D  = max(G%D0(i,j)*Mp_cgs, St%Vars(i,j,DENSITY)) ! [g/cm^3]
               Cs = St%Vars(i,j,SOUND_SPEED)*sqrt(St%Vars(i,j,DENSITY)/D) ! [cm/s]
            else
               D  = St%Vars(i,j,DENSITY)
               Cs = St%Vars(i,j,SOUND_SPEED)
            endif
            ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
            E0 = aurora%alpha*Mp_cgs*heFrac*erg2kev*Cs**2*aurora%rampFactor(i,j)
            aurora%E0  (i,j) = E0
            ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
            phi0 = sqrt(kev2erg)/(heFrac*Mp_cgs)**1.5D0*aurora%beta*D*sqrt(E0)*aurora%rampFactor(i,j)
            aurora%phi0(i,j) = phi0
            
            ! resistence out of the ionosphere is 2*rout resistence into the
            ! ionosphere is 2*rin outward current is positive
            if( signOfJ*St%Vars(i,j,FAC) >=0. ) then
              aRes = 2.D0*Rout
            else
              aRes = 2.D0*Rin
            endif

            ! Density floor to limit characteristic energy.  See Wiltberger et al. 2009 for details.
            D = min(D,rhoFactor*conductance%euvSigmaP(i,j))
            ! Original form for ref: 
            ! (heFrac*Mp_cgs)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*aurora%R*aRes*signOfJ*(St%Vars(:,:,FAC)*1.D-6)*sqrt(E0)/D
            dE = (heFrac*Mp_cgs)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*aurora%R*aRes*signOfJ*(St%Vars(i,j,FAC)*1.D-6)*sqrt(E0)/D
            aurora%deltaE(i,j) = dE
            dEc = min(dE,maxDrop)

            St%Vars(i,j,AVG_ENG) = max( E0+dEc, eTINY)
            if(dEc>0.) then
              St%Vars(i,j,NUM_FLUX) = phi0*(8.D0-7.D0*exp(-dEc/7.D0/E0))
            else
              St%Vars(i,j,NUM_FLUX) = phi0*exp(dEc/E0)
            endif
            St%Vars(i,j,DELTAE) = aurora%deltaE(i,j) ! [kV]
         enddo ! i
      enddo ! j

    end subroutine dragonking_fedder95

    subroutine dragonking_sunny(aurora,G,St)
      ! Assign zero precipitation when only EUV conductance is used.
      type(dkParams_T), intent(inout) :: aurora
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      ! Make sure zero Eavg/Nflux does not cause issues in conductance calculation.
      St%Vars(:,:,AUR_TYPE) = AT_NoPre
      St%Vars(:,:,AVG_ENG) = 0.D0
      St%Vars(:,:,NUM_FLUX) = 0.D0
      St%Vars(:,:,DELTAE) = 0.D0

    end subroutine dragonking_sunny

    subroutine dragonking_linmrg(aurora,G,St)
      ! Derive mono-diffuse electron precipitation where mono is based on linearized FL relation,
      ! and diffuse is a combination of MHD and RCM precipitation.
      type(dkParams_T), intent(inout) :: aurora
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      integer :: i,j
      real(rp) :: wRCM,wMHD
      real(rp) :: mhd_nflx,mhd_eflx,mhd_eavg
      real(rp) :: rcm_nflx,rcm_eflx,rcm_eavg
      real(rp) :: mix_nflx,mhd_SigPH,rcm_SigPH
      logical :: isMono,isPSP

      St%Vars(:,:,AUR_TYPE) = 0

      ! derive spatially varying beta and gtype 
      ! using RCM precipitation and thermal fluxes.
      call dragonking_beta_gtype(G,St)

      ! Derive mono using the linearized FL relation.
      call dragonking_linmono(aurora,G,St)

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,isMono,isPSP,wRCM,wMHD) &
      !$OMP private(mhd_nflx,mhd_eflx,mhd_eavg) &
      !$OMP private(rcm_nflx,rcm_eflx,rcm_eavg) &
      !$OMP private(mix_nflx,mhd_SigPH,rcm_SigPH)
      do j=1,G%Nt
         do i=1,G%Np
            !Start by figuring out where we are
            isPSP = St%Vars(i,j,IM_NPSP)>=Ne_psp .or. beta_RCM(i,j)<=0.01
            isMono = St%Vars(i,j,DELTAE)>eTINY
            !Grab values
            mhd_nflx = St%Vars(i,j,NUM_FLUX)
            if(mhd_nflx>GuABNF) then
               mhd_eavg = St%Vars(i,j,AVG_ENG)
            else
               mhd_eavg = eTINY
            endif
            mhd_eflx = mhd_eavg*mhd_nflx*kev2erg

            rcm_nflx = St%Vars(i,j,IM_ENFLX)
            if(rcm_nflx>GuABNF) then
               rcm_eavg = St%Vars(i,j,IM_EFLUX)/(rcm_nflx*kev2erg)
            else
               rcm_eavg = eTINY
            endif
            rcm_eflx = rcm_eavg*rcm_nflx*kev2erg

            if(isPSP) then
               ! Set auroral type to diffuse or no precipitation. 
               ! Use RCM values for diffuse nflux and eavg in the plasmasphere.
               St%Vars(i,j,NUM_FLUX) = rcm_nflx
               St%Vars(i,j,AVG_ENG)  = rcm_eavg
               if(rcm_nflx>GuABNF) then
                  St%Vars(i,j,AUR_TYPE) = AT_RMnoE
               else
                  St%Vars(i,j,AUR_TYPE) = AT_NoPre
               endif
            elseif(isMono .and. .not.isPSP) then
               ! Set auroral type to mono if mono is above Gussenhove threshold and gives higher Robinson L2.
               ! Else set diffuse or no precipitation depending on rcm values relative to the threshold.
               mhd_SigPH = SigmaP_Robinson(mhd_eavg,mhd_eflx)**2+SigmaH_Robinson(mhd_eavg,mhd_eflx)**2
               rcm_SigPH = SigmaP_Robinson(rcm_eavg,rcm_eflx)**2+SigmaH_Robinson(rcm_eavg,rcm_eflx)**2
               St%Vars(i,j,AUR_TYPE) = AT_RMono
               if(mhd_nflx>GuABNF .and. mhd_SigPH>rcm_SigPH) then
                  St%Vars(i,j,AUR_TYPE) = AT_RMono
               else
                  St%Vars(i,j,NUM_FLUX) = rcm_nflx
                  St%Vars(i,j,AVG_ENG)  = rcm_eavg
                  if(rcm_nflx>GuABNF) then
                     St%Vars(i,j,AUR_TYPE) = AT_RMfnE
                  else
                     St%Vars(i,j,AUR_TYPE) = AT_NoPre
                  endif
               endif
            else
               ! Linearly merge MHD and RCM diffuse nflux and eflux.
               ! Note where deltaE>eTINY but beta_RCM<=0.01, gtype will be near 1.
               wRCM = gtype_RCM(i,j)
               wMHD = 1.0-wRCM
               mix_nflx = wMHD*mhd_nflx + wRCM*rcm_nflx
               St%Vars(i,j,NUM_FLUX) = mix_nflx
               if(mix_nflx>GuABNF) then
                  St%Vars(i,j,AVG_ENG) = (wMHD*mhd_eflx+wRCM*rcm_eflx)/(mix_nflx*kev2erg)
                  St%Vars(i,j,AUR_TYPE)= AT_RMnoE
               else
                  St%Vars(i,j,AVG_ENG) = eTINY
                  St%Vars(i,j,AUR_TYPE)= AT_NoPre
               endif
            endif
         enddo
      enddo

    end subroutine dragonking_linmrg

    subroutine dragonking_linmono(aurora,G,St)
      ! Linearized Fridman-Lemaire relation for mono.
      type(dkParams_T), intent(inout) :: aurora
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      
      real(rp) :: signOfJ
      integer :: i,j
      real(rp) :: D,Cs,Pe,Ne,J2eF0,eV2kT,dE,phi0,kT
      real(rp) :: Pe_mhd,Ne_mhd,kT_mhd,phi0_mhd
      real(rp) :: Pe_rcm,Ne_rcm,kT_rcm,phi0_rcm
      real(rp) :: wMHD,wRCM
      logical :: hasMono

      if (St%hemisphere==NORTH) then
         signOfJ = -1  
      elseif (St%hemisphere==SOUTH) then
         signOfJ = 1
      else
         stop 'Wrong hemisphere label. Stopping...'
      endif

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,D,Cs,Pe,Ne,J2eF0,eV2kT,dE,phi0,kT) &
      !$OMP private(Pe_mhd,Ne_mhd,kT_mhd,phi0_mhd) &
      !$OMP private(Pe_rcm,Ne_rcm,kT_rcm,phi0_rcm) &
      !$OMP private(wMHD, wRCM, hasMono)
      do j=1,G%Nt
         do i=1,G%Np
            if (aurora%doChill) then
               ! MHD density replaced with gallagher where it's lower      
               ! and temperature changed correspondingly
               D  = max(G%D0(i,j)*Mp_cgs, St%Vars(i,j,DENSITY)) ! [g/cm^3]
               Cs = St%Vars(i,j,SOUND_SPEED)*sqrt(St%Vars(i,j,DENSITY)/D) ! [cm/s]
            else
               D  = St%Vars(i,j,DENSITY)
               Cs = St%Vars(i,j,SOUND_SPEED)
            endif
            ! electron pressure from MHD side in [Pa]
            ! 0.1 is to convert [g/cm^3]*[cm/s]^2=[g/cm/s^2] to [Pa].
            Pe_mhd = 0.1/MIXgamma*alpha_RCM(i,j)*D*(Cs**2)
            ! electron number density from MHD side in [/m^3].
            Ne_mhd = D/(Mp_cgs*heFrac)*1.0D6
            ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
            kT_mhd = Pe_mhd/(Ne_mhd*kev2J)
            ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
            ! Treat MHD flux on closed field lines as trapped flux.
            ! Beta now stands for the ratio between precipitation and trapped.
            ! Note phi0 can be zero in the polar cap.
            phi0_mhd = beta_RCM(i,j)*sqrt(Pe_mhd*Ne_mhd/(2.0D-3*PI*Me_cgs))*1.0D-4

            ! RCM flux
            Pe_rcm = St%Vars(i,j,IM_EPRE)
            Ne_rcm = St%Vars(i,j,IM_EDEN)
            phi0_rcm = St%Vars(i,j,IM_ENFLX)

            ! Merged flux
            wRCM = gtype_RCM(i,j)
            wMHD = 1.0-wRCM
            phi0 = wMHD*phi0_mhd + wRCM*phi0_rcm
            Ne   = wMHD*Ne_mhd + wRCM*Ne_rcm
            Pe   = wMHD*Pe_mhd + wRCM*Pe_rcm
            if(Ne>TINY) then
               kT = Pe/(Ne*kev2J)
            else
               kT = eTINY
            endif

            hasMono = .false.

            aurora%E0  (i,j) = 2.0*kT
            aurora%phi0(i,j) = phi0
            
            ! Note JF0 may go inf where phi0/Pe_MHD/Ne_MHD is zero.
            if (phi0 > TINY) then
               J2eF0 = signOfJ*(St%Vars(i,j,FAC)*1.0D-6)/(eCharge*phi0*1.0D4)
               ! Linearize the Fridman-Lemaire relation: 
               ! eV = kB*Te*(RM-1)*ln((RM-1)/(RM-J/e/F0)) when 1<=J/e/F0<=RM.
               ! eV ~ kB*Te*(J/e/F0-1) when RM>>J/e/F0
               if (J2eF0>1.0 .and. Ne>Ne_floor) then
                  hasMono = .true.
               endif
            endif

            if (hasMono) then
               dE = kT*(J2eF0-1.0)
               aurora%deltaE(i,j) = dE
               eV2kT = min(dE,maxDrop)/kT + 1.0
               St%Vars(i,j,Z_NFLUX) = phi0*eV2kT
               St%Vars(i,j,Z_EAVG) = kT*(eV2kT+1.0/eV2kT)
            else
               ! When there is no acceleration, MHD flux represents the diffuse precipitation.
               aurora%deltaE(i,j) = 0.0
               St%Vars(i,j,Z_NFLUX) = phi0
               St%Vars(i,j,Z_EAVG) = 2.0*kT
            endif

            St%Vars(i,j,DELTAE)   = aurora%deltaE(i,j)
            St%Vars(i,j,AVG_ENG)  = max(St%Vars(i,j,Z_EAVG),eTINY)   ! [keV]
            St%Vars(i,j,NUM_FLUX) = St%Vars(i,j,Z_NFLUX)  ! [#/cm^2/s]
         enddo ! i
      enddo ! j
    end subroutine dragonking_linmono

    subroutine dragonking_mr(aurora,G,St)
      ! George Khazanov's multiple reflection(MR) corrections
      ! Modify the diffuse precipitation energy and mean energy based on equation (5) in Khazanov et al. [2019JA026589]
      ! Kc = 3.36-exp(0.597-0.37*Eavg+0.00794*Eavg^2)
      ! Eavg_c = 0.073+0.933*Eavg-0.0092*Eavg^2
      ! Nflx_c = Eflx_c/Eavg_c = Eflx*Kc/Eavg_c = Nflx*Eavg*Kc/Eavg_c
      ! where Eavg is the mean energy in [keV] before MR, Eavg_c is the modified mean energy. 
      ! Kc is the ratio between modified enflux and unmodified enflux.
      ! No longer use the conductance modification factor in Khazanov et al. [2018SW001837] (eq 7).
      type(dkParams_T), intent(in) :: aurora
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      real(rp) :: eavg,nflx,eavgMR,nflxMR
      integer :: i,j

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,eavg,nflx,eavgMR,nflxMR)
      do j=1,G%Nt
        do i=1,G%Np
          if(aurora%deltaE(i,j)<=0.0) then
            ! Only modify diffuse precipitation. May modify mono precipitation when a formula is available.
            eavg = St%Vars(i,j,AVG_ENG)
            nflx = St%Vars(i,j,NUM_FLUX)
            call AugmentMR(eavg,nflx,eavgMR,nflxMR)
            St%Vars(i,j,AVG_ENG ) = eavgMR
            St%Vars(i,j,NUM_FLUX) = nflxMR
          endif
        enddo
      enddo
    end subroutine dragonking_mr

    subroutine AugmentMR(eavg,nflx,eavgMR,nflxMR)
      ! Cleaned up routine to calculate MR-augmented eavg and nflx
      real(rp), intent(in)  :: eavg  ,nflx
      real(rp), intent(out) :: eavgMR,nflxMR
      real(rp) :: eflux,Kc

      if ( (eavg<=30.0) .and. (eavg>=0.5) ) then
         ! Kc becomes negative when E_avg > 47-48 keV.
         Kc = 3.36 - exp(0.597-0.37*eavg+0.00794*eavg**2)
         eflux = eavg*nflx
         eavgMR = 0.073+0.933*eavg-0.0092*eavg**2 ! [keV]
         nflxMR = Kc*eflux/eavgMR
      else
         !Nothing to do
         eavgMR = eavg
         nflxMR = nflx
      endif
    end subroutine AugmentMR

    subroutine dragonking_auroralsmooth(St,G,aurora)
      type(mixState_T), intent(inout) :: St
      type(mixGrid_T), intent(in) :: G      
      type(dkParams_T), intent(inout) :: aurora
      integer :: i, j
      logical :: smthDEPonly = .false.
      integer :: smthE
      real(rp) :: tmpC(G%Np,G%Nt),tmpD(G%Np,G%Nt),tmpE(G%Np+4,G%Nt+4),tmpF(G%Np+4,G%Nt+4)
      smthE = 1 ! 1. smooth EFLUX; 2. smooth EAVG

      tmpC = 0.D0
      tmpD = 0.D0
      tmpE = 0.D0
      tmpF = 0.D0

      ! Test only smoothing diffuse precipitation.
      ! Diffuse mask is St%Vars(:,:,Z_NFLUX)<0.
      ! Get AVG_ENG*mask and NUM_FLUX*mask for smoothing.
      ! Fill AVG_ENG and NUM_FLUX with smoothed where mask is true.
      ! this domain of conductance is never used. Here it's used for diffuse precipitation mask. Be careful if not in RCMONO mode.
      aurora%PrecipMask = 1.D0 
      if(smthDEPonly) then
        where(St%Vars(:,:,DELTAE)>0)
          aurora%PrecipMask = 0.D0
        end where
        ! back up mono
        tmpC = St%Vars(:,:,AVG_ENG)
        tmpD = St%Vars(:,:,NUM_FLUX)
      endif

      ! get expanded diffuse with margin grid (ghost cells)
      if(smthE==1) then
        call dragonking_margin(G,St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)*aurora%PrecipMask, tmpE)
      else
        call dragonking_margin(G,St%Vars(:,:,AVG_ENG)*aurora%PrecipMask, tmpE)
      end if
      call dragonking_margin(G,St%Vars(:,:,NUM_FLUX)*aurora%PrecipMask, tmpF)

      !$OMP PARALLEL DO default(shared) collapse(2) &
      !$OMP private(i,j)
      do j=1,G%Nt
        do i=1,G%Np
          ! SmoothOperator55(Q,isGO)
          St%Vars(i,j,NUM_FLUX)=SmoothOperator55(tmpF(i:i+4,j:j+4))
          if(smthE==1) then
            St%Vars(i,j,AVG_ENG) = max(SmoothOperator55(tmpE(i:i+4,j:j+4))/St%Vars(i,j,NUM_FLUX),1.D-8)
          else
            St%Vars(i,j,AVG_ENG) = max(SmoothOperator55(tmpE(i:i+4,j:j+4)),1.D-8)
          end if
        end do
      end do

      if(smthDEPonly) then
        ! combine smoothed diffuse and unsmoothed mono
        where(St%Vars(:,:,DELTAE)>0)
          St%Vars(:,:,AVG_ENG)  = tmpC
          St%Vars(:,:,NUM_FLUX) = tmpD
        end where
      endif
    end subroutine dragonking_auroralsmooth

    subroutine dragonking_ramp(aurora,G,rPolarBound,rEquatBound,rLowLimit)
      type(dkParams_T), intent(inout) :: aurora
      type(mixGrid_T), intent(in) :: G      
      real(rp), intent(in) :: rPolarBound,rEquatBound,rLowLimit
      
      where (G%r <= rPolarBound)
        aurora%rampFactor = 1.0D0
      elsewhere ( (G%r > rPolarBound).and.(G%r <= rEquatBound) )
        aurora%rampFactor = 1.0D0+(rLowLimit - 1.0D0)*(G%r - rPolarBound)/(rEquatBound-rPolarBound)
      elsewhere ! G%r > rEquatBound
        aurora%rampFactor = rLowLimit
      end where
    end subroutine dragonking_ramp

    subroutine dragonking_beta_gtype(G,St)
      ! Use RCM precipitation and source population to derive the loss cone rate beta.
      ! Assume beta is one in the polar cap and smooth across RCM boundary.
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      logical :: isAncB(G%Np,G%Nt),isAncG(G%Np,G%Nt)
      real(rp) :: phi0_rcm
      integer :: i,j

      ! Initialize IM_BETA with xml input value.
      ! IM_GTYPE is interpolated from RCM: 1=RCM and 0=MHD.
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j)
      do j=1,G%Nt
         do i=1,G%Np
            St%Vars(i,j,IM_BETA) = beta_inp
         enddo
      enddo
      isAncB = .false. ! beta is smoothed everywhere.
      isAncG = .false.

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,phi0_rcm)
      do j=1,G%Nt
         do i=1,G%Np
            ! Set grids outside RCM as anchors that won't be smoothed.
            if(St%Vars(i,j,IM_GTYPE)<=0.01) then
               isAncG(i,j) = .true.
            endif

            ! Derive beta as the ratio between RCM precipitation and trapped flux.
            ! sqrt([Pa]*[#/m^3]/[kg]) = sqrt([#/m^4/s^2]) = 1e-4*[#/cm^2/s]
            phi0_rcm = sqrt(St%Vars(i,j,IM_EPRE)*St%Vars(i,j,IM_EDEN)/(Me_cgs*1e-3*2*pi))*1.0e-4
            if(phi0_rcm>TINY) then
               ! This criterion includes all RCM grid. Note beta is 0 where IM_ENFLX is zero.
               St%Vars(i,j,IM_BETA) = St%Vars(i,j,IM_ENFLX)/phi0_rcm
               if(St%Vars(i,j,IM_BETA)>1.0) then
                  ! Near RCM high lat boundary, trapped flux is not well constrained.
                  ! The beta there can be easily >>1 and is thus not as reliable.
                  ! Reset grid weight to MHD if RCM beta>1.
                  St%Vars(i,j,IM_GTYPE) = 0.0 !1.0/St%Vars(i,j,IM_BETA)**2
                  St%Vars(i,j,IM_BETA)  = 0.5*(1.0+beta_inp)
                  isAncG(i,j) = .false.
               endif
            endif
         enddo
      enddo

      ! Smooth IM_BETA and save in beta_RCM
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j)
      do j=1,G%Nt
         do i=1,G%Np
            beta_RCM(i,j) = St%Vars(i,j,IM_BETA)
         enddo
      enddo
      call dragonking_smooth(G,beta_RCM,isAncB)
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j)
      do j=1,G%Nt
         do i=1,G%Np
            St%Vars(i,j,IM_BETA) = beta_RCM(i,j)
            if(beta_RCM(i,j)>1.D0) beta_RCM(i,j) = 1.D0
            if(beta_RCM(i,j)<0.D0) beta_RCM(i,j) = 0.D0
         enddo
      enddo

      ! Smooth IM_GTYPE and save in gtype_RCM
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j)
      do j=1,G%Nt
         do i=1,G%Np
            gtype_RCM(i,j) = St%Vars(i,j,IM_GTYPE)
         enddo
      enddo
      call dragonking_smooth(G,gtype_RCM,isAncG)
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j)
      do j=1,G%Nt
         do i=1,G%Np
            St%Vars(i,j,IM_GTYPE) = gtype_RCM(i,j)
            if(gtype_RCM(i,j)>1.D0) gtype_RCM(i,j) = 1.D0
            if(gtype_RCM(i,j)<0.D0) gtype_RCM(i,j) = 0.D0
         enddo
      enddo
    end subroutine dragonking_beta_gtype

    subroutine dragonking_smooth(Gr,Q,isAnchor)
      ! Do smoothing window on ReMIX grid quantity
      ! Skip anchor points
      type(mixGrid_T), intent(in) :: Gr
      real(rp), intent(inout) :: Q(Gr%Np,Gr%Nt)
      logical, intent(in) :: isAnchor(Gr%Np,Gr%Nt)
      real(rp) :: temp(Gr%Np,Gr%Nt)
      real(rp) :: thres,mad,Ttmp,madij
      integer :: i,j,it,im1,ip1,jm1,jp1,MaxIter

      thres = 0.025
      MaxIter = 15
      call FixPole(Gr,Q)

      do it=1,MaxIter
        mad = 0.D0 ! max abs difference from last iteration.
        temp = Q

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,jm1,jp1,im1,ip1,Ttmp,madij) &
        !$OMP reduction(max:mad)
        do j=1,Gr%Nt ! use open BC for lat.
          do i=1,Gr%Np ! use periodic BC for lon.
            ! Do not smooth anchor points.
            if(.not. isAnchor(i,j)) then
              jm1 = j-1
              jp1 = j+1
              if (j == 1)     jm1 = 1
              if (j == Gr%Nt) jp1 = Gr%Nt

              im1 = i-1
              ip1 = i+1
              if (i == 1)     im1 = Gr%Np
              if (i == Gr%Np) ip1 = 1

              Ttmp =(temp(im1,jm1)+temp(im1,j)+temp(im1,jp1) &
                   + temp(i  ,jm1)+temp(i  ,j)+temp(i  ,jp1) &
                   + temp(ip1,jm1)+temp(ip1,j)+temp(ip1,jp1))/9.D0
              madij = abs(Q(i,j)-Ttmp)
              Q(i,j) = Ttmp
            else
              madij = 0.D0
            endif
            mad = max(mad,madij)
          enddo
        enddo
        call FixPole(Gr,Q)
        if(mad<thres) exit
      enddo
    end subroutine dragonking_smooth

    !Enforce pole condition that all values at the same point are equal
    subroutine FixPole(Gr,Q)
      type(mixGrid_T), intent(in)    :: Gr
      real(rp)       , intent(inout) :: Q(Gr%Np,Gr%Nt)

      real(rp) :: Qavg

      Qavg = sum(Q(:,2))/Gr%Np
      Q(:,1) = Qavg
    end subroutine FixPole

    subroutine dragonking_margin(G,array,arraymar)
      ! Conductance boundary processing.
      type(mixGrid_T), intent(in) :: G      
      real(rp),dimension(G%Np,G%Nt), intent(in) :: array
      real(rp),dimension(G%Np+4,G%Nt+4), intent(out) :: arraymar
      arraymar(3:G%Np+2,3:G%Nt+2) = array
      ! reflective BC for lat
      arraymar(3:G%Np+2,1) = array(:,3)
      arraymar(3:G%Np+2,2) = array(:,2)
      arraymar(3:G%Np+2,G%Nt+3) = array(:,G%Nt-1)
      arraymar(3:G%Np+2,G%Nt+4) = array(:,G%Nt-2)
      ! periodic BC for lon
      arraymar(1:2,:) = arraymar(G%Np+1:G%Np+2,:)
      arraymar(G%Np+3:G%Np+4,:) = arraymar(3:4,:)
    end subroutine dragonking_margin

    subroutine SetMIXgamma(gamma)
      !Routine to change MIX-gamma on the fly if necessary
      real(rp), intent(in) :: gamma
      MIXgamma = gamma
    end subroutine SetMIXgamma

  end module dkaurora
