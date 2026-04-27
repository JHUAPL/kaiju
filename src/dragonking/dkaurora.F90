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

  !Replacing some hard-coded inline values (bad) w/ module private values (slightly less bad)
  real(rp), parameter, private :: maxDrop = 20.0 !Hard-coded max potential drop [kV]
  real(rp), parameter, private :: eTINY = mixeTINY ! Floor of average energy [keV]
  real(rp), parameter, private :: Ne_floor = 0.03e6 ! minimum Ne in [/m^3] when evaluating the linearized FL relation.
  real(rp), parameter, private :: Ne_psp = 10.0e6 ! Ne threshold for the plasmasphere in [/m^3].
  real(rp), parameter, private :: GuABNF = 1.e7 ! Gussenhoven+[1983] Auroral Boundary Number Flux in [#/cm^s/s].

  contains

    subroutine precipitation_total(aurora,grid,vars,conductance,params)
      ! The main subroutine in module dragonking.
      ! It derives auroral precipitation.
      type(dkAurora_T), intent(inout) :: aurora
      type(mixGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkConductance_T), intent(inout) :: conductance
      type(dkParameters_T), intent(in) :: params

      integer :: i, j
      logical :: isAncB(grid%Np, grid%Nt)
      logical :: isAncG(grid%Np, grid%Nt)
      real(rp) :: signOfJ
      integer :: Np, Nt
      real(rp) :: r(grid%Np, grid%Nt), D0(grid%Np, grid%Nt)

      Np = grid%Np
      Nt = grid%Nt
      D0 = grid%D0
      r  = grid%r

      isAncB = .false.
      isAncG = .false.

      select case (params%aurora_model_type)
        case (FEDDER)
         !$OMP PARALLEL DO collapse(2) default(shared)
         do j=1,Nt
            do i=1,Np
               call dragonking_fedder95(i,j,aurora,vars,conductance,params,r(i,j),D0(i,j))
            end do
         end do
         !$OMP END PARALLEL DO

        case (SUNNY)
          call dragonking_sunny(aurora,vars)
        
        case (LINMRG)
            if (vars%hemisphere == NORTH) then
               signOfJ = -1.0_rp
            elseif (vars%hemisphere == SOUTH) then
               signOfJ =  1.0_rp
            else
               stop 'Wrong hemisphere label. Stopping...'
            endif

            !$OMP PARALLEL DO collapse(2) default(shared)
            do j = 1,Nt
               do i = 1,Np
                  call dragonking_beta_gtype(i,j,vars,aurora,params,isAncG(i,j))
               enddo
            enddo
            !$OMP END PARALLEL DO

            call dragonking_smooth(Np,Nt,aurora%beta_RCM,isAncB)
            call dragonking_smooth(Np,Nt,aurora%gtype_RCM,isAncG)

            !$OMP PARALLEL DO collapse(2) default(shared)
            do j = 1,Nt
               do i = 1,Np
                  vars%Im_Beta(i,j)  = aurora%beta_RCM(i,j)
                  aurora%beta_RCM(i,j)  = min(max(aurora%beta_RCM(i,j),  0.0_rp), 1.0_rp)
                  
                  vars%Im_GTYPE(i,j) = aurora%gtype_RCM(i,j)
                  aurora%gtype_RCM(i,j) = min(max(aurora%gtype_RCM(i,j), 0.0_rp), 1.0_rp)

                  call dragonking_linmrg(i, j, aurora, vars, params, D0(i,j), signOfJ)
               enddo
            enddo
            !$OMP END PARALLEL DO

        case default
          stop "The auroral precipitation model type entered is not supported."
      end select

      ! Correct for multiple reflections if you're so inclined
      if (params%doMR) then
         !$OMP PARALLEL DO collapse(2) default(shared)
         do j = 1, Nt
            do i = 1, Np
               call dragonking_mr(aurora%deltaE(i,j),vars%Avg_Eng(i,j),vars%Num_Flux(i,j))
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

      ! Smooth precipitation energy and flux before calculating conductance.
      if (params%doAuroralSmooth) call dragonking_auroralsmooth(Np,Nt,vars,aurora)

    end subroutine precipitation_total

    subroutine dragonking_fedder95(i,j,aurora,vars,conductance,params,r,D0)
      ! Derive electron precipitation energy flux and avg energy using Fedder95 formula.
      integer, intent(in) :: i,j
      type(dkAurora_T), intent(inout) :: aurora
      type(dkVariables_T), intent(inout) :: vars
      type(dkConductance_T), intent(inout) :: conductance
      type(dkParameters_T), intent(in) :: params
      real(rp), intent(in) :: r, D0
      
      real(rp) :: signOfJ
      real(rp) :: Rout = 6.D0, Rin = 1.2D0
      real(rp) :: rhoFactor = 3.3D-24*0.5D0
      real(rp) :: rPolarBound = 20.0D0*pi/180.D0, rEquatBound = 30.0D0*pi/180.0D0, rLowLimit = 0.02D0
      real(rp) :: D,Cs,Pe,Ne,dE,phi0,E0,aRes,dEc

      if (vars%hemisphere==NORTH) then
         signOfJ = -1  
      elseif (vars%hemisphere==SOUTH) then
         signOfJ = 1
      else
         stop 'Wrong hemisphere label. Stopping...'
      endif

      ! fills in rampFactor
      if (params%doRamp) then
         call dragonking_ramp(aurora%rampFactor(i,j),r,rPolarBound,rEquatBound,rLowLimit)
      else 
         aurora%rampFactor(i,j) = 1.0D0
      endif

      if (params%doChill) then
         ! MHD density replaced with gallagher where it's lower      
         ! and temperature changed correspondingly
         D  = max(D0*Mp_cgs, vars%Density(i,j)) ! [g/cm^3]
         Cs = vars%Sound_Speed(i,j)*sqrt(vars%Density(i,j)/D) ! [cm/s]
      else
         D  = vars%Density(i,j)
         Cs = vars%Sound_Speed(i,j)
      endif
      ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
      E0 = params%alpha*Mp_cgs*heFrac*erg2kev*Cs**2*aurora%rampFactor(i,j)
      aurora%E0(i,j) = E0
      ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
      phi0 = sqrt(kev2erg)/(heFrac*Mp_cgs)**1.5D0*params%beta*D*sqrt(E0)*aurora%rampFactor(i,j)
      aurora%phi0(i,j) = phi0
      
      ! resistence out of the ionosphere is 2*rout resistence into the
      ! ionosphere is 2*rin outward current is positive
      if( signOfJ*vars%FAC(i,j) >=0. ) then
         aRes = 2.D0*Rout
      else
         aRes = 2.D0*Rin
      endif

      ! Density floor to limit characteristic energy.  See Wiltberger et al. 2009 for details.
      D = min(D,rhoFactor*conductance%euvSigmaP(i,j))
      ! Original form for ref: 
      ! (heFrac*Mp_cgs)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*aurora%R*aRes*signOfJ*(St%Vars(:,:,FAC)*1.D-6)*sqrt(E0)/D
      dE = (heFrac*Mp_cgs)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*params%R*aRes*signOfJ*(vars%FAC(i,j)*1.D-6)*sqrt(E0)/D
      aurora%deltaE(i,j) = dE
      dEc = min(dE,maxDrop)

      vars%Avg_Eng(i,j) = max( E0+dEc, eTINY)
      if(dEc>0.) then
         vars%Num_Flux(i,j) = phi0*(8.D0-7.D0*exp(-dEc/7.D0/E0))
      else
         vars%Num_Flux(i,j) = phi0*exp(dEc/E0)
      endif
      vars%DeltaE(i,j) = aurora%deltaE(i,j) ! [kV]

    end subroutine dragonking_fedder95

    subroutine dragonking_sunny(aurora,vars)
      ! Assign zero precipitation when only EUV conductance is used.
      type(dkAurora_T), intent(inout) :: aurora
      type(dkVariables_T), intent(inout) :: vars

      ! Make sure zero Eavg/Nflux does not cause issues in conductance calculation.
      vars%Aur_Type(:,:) = AT_NoPre
      vars%Avg_Eng(:,:) = 0.D0
      vars%Num_Flux(:,:) = 0.D0
      vars%DeltaE(:,:) = 0.D0

    end subroutine dragonking_sunny

    subroutine dragonking_linmrg(i,j,aurora,vars,params,D0,signOfJ)
      ! Derive mono-diffuse electron precipitation where mono is based on linearized FL relation,
      ! and diffuse is a combination of MHD and RCM precipitation.
      integer, intent(in) :: i,j
      type(dkAurora_T), intent(inout) :: aurora
      type(dkVariables_T), intent(inout) :: vars
      type(dkParameters_T), intent(in) :: params
      ! type(mixGrid_T), intent(in) :: G
      real(rp), intent(in) :: D0, signOfJ

      real(rp) :: E0, phi0, deltaE, avg_eng, num_flux
      integer  :: aur_type

      vars%Aur_Type(i,j) = 0

      call dragonking_linmono(i, j, vars, params, aurora, D0, signOfJ, &
                     E0, phi0, deltaE, avg_eng, num_flux)
      aurora%E0  (i,j)      = E0
      aurora%phi0(i,j)      = phi0
      vars%DeltaE(i,j)   = deltaE
      vars%Avg_Eng(i,j)  = avg_eng
      vars%Num_Flux(i,j) = num_flux

      call dragonking_lindiff(i, j, vars, aurora, num_flux, avg_eng, aur_type)
      vars%Num_Flux(i,j) = num_flux
      vars%Avg_Eng(i,j)  = avg_eng
      vars%Aur_Type(i,j) = aur_type
    end subroutine dragonking_linmrg

    subroutine dragonking_linmono(i,j,vars,params,aurora,D0,signOfJ,E0,phi0,deltaE,avg_eng,num_flux)
      ! Linearized Fridman-Lemaire relation for mono.
      integer, intent(in) :: i,j
      type(dkVariables_T), intent(inout) :: vars
      type(dkParameters_T), intent(in) :: params
      type(dkAurora_T), intent(inout) :: aurora
      real(rp), intent(in) :: D0, signOfJ
      real(rp), intent(out) :: E0, phi0, deltaE, avg_eng, num_flux


      real(rp) :: D,Cs,Pe,Ne,J2eF0,eV2kT,dE,kT
      real(rp) :: Pe_mhd,Ne_mhd,phi0_mhd
      real(rp) :: Pe_rcm,Ne_rcm,phi0_rcm
      real(rp) :: wMHD,wRCM

      if (params%doChill) then
         ! MHD density replaced with gallagher where it's lower      
         ! and temperature changed correspondingly
         D  = max(D0*Mp_cgs, vars%Density(i,j)) ! [g/cm^3]
         Cs = vars%Sound_Speed(i,j)*sqrt(vars%Density(i,j)/D) ! [cm/s]
      else
         D  = vars%Density(i,j)
         Cs = vars%Sound_Speed(i,j)
      endif

      ! electron pressure from MHD side in [Pa]
      ! 0.1 is to convert [g/cm^3]*[cm/s]^2=[g/cm/s^2] to [Pa].
      Pe_mhd = 0.1/params%gamma*aurora%alpha_RCM(i,j)*D*(Cs**2)
      ! electron number density from MHD side in [/m^3].
      Ne_mhd = D/(Mp_cgs*heFrac)*1.0D6
      
      ! ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
      ! kT_mhd = Pe_mhd/(Ne_mhd*kev2J)

      ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
      ! Treat MHD flux on closed field lines as trapped flux.
      ! Beta now stands for the ratio between precipitation and trapped.
      ! Note phi0 can be zero in the polar cap.
      phi0_mhd = aurora%beta_RCM(i,j)*sqrt(Pe_mhd*Ne_mhd/(2.0D-3*PI*Me_cgs))*1.0D-4

      ! RCM flux
      Pe_rcm = vars%Im_Epre(i,j)
      Ne_rcm = vars%Im_Eden(i,j)
      phi0_rcm = vars%Im_Enflx(i,j)

      ! Merged flux
      wRCM = aurora%gtype_RCM(i,j)
      wMHD = 1.0_rp - wRCM
      phi0 = wMHD*phi0_mhd + wRCM*phi0_rcm
      Ne   = wMHD*Ne_mhd + wRCM*Ne_rcm
      Pe   = wMHD*Pe_mhd + wRCM*Pe_rcm

      if(Ne>TINY) then
         kT = Pe/(Ne*kev2J)
      else
         kT = eTINY
      endif
      E0 = 2.0_rp*kT
      
      ! Note JF0 may go inf where phi0/Pe_MHD/Ne_MHD is zero.
      J2eF0 = signOfJ*(vars%FAC(i,j)*1.0D-6)/(eCharge*phi0*1.0D4)
      ! Linearize the Fridman-Lemaire relation: 
      ! eV = kB*Te*(RM-1)*ln((RM-1)/(RM-J/e/F0)) when 1<=J/e/F0<=RM.
      ! eV ~ kB*Te*(J/e/F0-1) when RM>>J/e/F0
      if (J2eF0>1.0_rp .and. Ne>Ne_floor) then
         dE = kT*(J2eF0-1.0_rp)
         deltaE = dE
         eV2kT = min(dE,maxDrop)/kT + 1.0_rp
         num_flux = phi0*eV2kT
         avg_eng = kT*(eV2kT+1.0_rp/eV2kT)
      else
         ! When there is no acceleration, MHD flux represents the diffuse precipitation.
         deltaE = 0.0_rp
         num_flux = phi0
         avg_eng = 2.0_rp*kT
      endif

      avg_eng  = max(avg_eng,eTINY)   ! [keV]
   
    end subroutine dragonking_linmono

   subroutine dragonking_lindiff(i, j, vars, aurora, num_flux, avg_eng, aur_type)
   ! Derive diffuse electron precipitation as a combination of MHD and RCM precipitation.

      integer, intent(in)  :: i, j
      type(dkVariables_T), intent(in)  :: vars
      type(dkAurora_T), intent(in) :: aurora
      real(rp), intent(out) :: num_flux, avg_eng
      integer, intent(out) :: aur_type

      logical  :: isMono, isPSP
      real(rp) :: mhd_nflx,mhd_eflx,mhd_eavg
      real(rp) :: rcm_nflx,rcm_eflx,rcm_eavg
      real(rp) :: mix_nflx,mhd_SigPH,rcm_SigPH
      real(rp) :: wMHD,wRCM

      !Start by figuring out where we are
      isPSP = vars%Im_Npsp(i,j)>=Ne_psp .or. aurora%beta_RCM(i,j)<=0.01_rp
      isMono = vars%DeltaE(i,j)>eTINY

      !Grab values
      mhd_nflx = vars%Num_Flux(i,j)
      if(mhd_nflx>GuABNF) then
         mhd_eavg = vars%Avg_Eng(i,j)
      else
         mhd_eavg = eTINY
      endif
      mhd_eflx = mhd_eavg*mhd_nflx*kev2erg

      rcm_nflx = vars%Im_Enflx(i,j)
      if(rcm_nflx>GuABNF) then
         rcm_eavg = vars%Im_Eflux(i,j)/(rcm_nflx*kev2erg)
      else
         rcm_eavg = eTINY
      endif
      rcm_eflx = rcm_eavg*rcm_nflx*kev2erg

      if(isPSP) then
         ! Set auroral type to diffuse or no precipitation. 
         ! Use RCM values for diffuse nflux and eavg in the plasmasphere.
         num_flux = rcm_nflx
         avg_eng = rcm_eavg
         if(rcm_nflx>GuABNF) then
            aur_type = AT_RMnoE
         else
            aur_type = AT_NoPre
         endif
      elseif(isMono) then
         ! Set auroral type to mono if mono is above Gussenhove threshold and gives higher Robinson L2.
         ! Else set diffuse or no precipitation depending on rcm values relative to the threshold.
         mhd_SigPH = SigmaP_Robinson(mhd_eavg,mhd_eflx)**2+SigmaH_Robinson(mhd_eavg,mhd_eflx)**2
         rcm_SigPH = SigmaP_Robinson(rcm_eavg,rcm_eflx)**2+SigmaH_Robinson(rcm_eavg,rcm_eflx)**2
         ! vars%Aur_Type(i,j) = AT_RMono

         if(mhd_nflx>GuABNF .and. mhd_SigPH>rcm_SigPH) then
            aur_type = AT_RMono
         else
            num_flux = rcm_nflx
            avg_eng = rcm_eavg
            if(rcm_nflx>GuABNF) then
               aur_type = AT_RMfnE
            else
               aur_type = AT_NoPre
            endif
         endif
      else
         ! Linearly merge MHD and RCM diffuse nflux and eflux.
         ! Note where deltaE>eTINY but beta_RCM<=0.01, gtype will be near 1.
         wRCM = aurora%gtype_RCM(i,j)
         wMHD = 1.0_rp-wRCM
         mix_nflx = wMHD*mhd_nflx + wRCM*rcm_nflx
         num_flux = mix_nflx
         if(mix_nflx>GuABNF) then
            avg_eng = (wMHD*mhd_eflx+wRCM*rcm_eflx)/(mix_nflx*kev2erg)
            aur_type = AT_RMnoE
         else
            avg_eng = eTINY
            aur_type = AT_NoPre
         endif
      endif
   end subroutine dragonking_lindiff

    subroutine dragonking_mr(deltaE,eavg,nflx)
      ! George Khazanov's multiple reflection(MR) corrections
      ! Modify the diffuse precipitation energy and mean energy based on equation (5) in Khazanov et al. [2019JA026589]
      ! Kc = 3.36-exp(0.597-0.37*Eavg+0.00794*Eavg^2)
      ! Eavg_c = 0.073+0.933*Eavg-0.0092*Eavg^2
      ! Nflx_c = Eflx_c/Eavg_c = Eflx*Kc/Eavg_c = Nflx*Eavg*Kc/Eavg_c
      ! where Eavg is the mean energy in [keV] before MR, Eavg_c is the modified mean energy. 
      ! Kc is the ratio between modified enflux and unmodified enflux.
      ! No longer use the conductance modification factor in Khazanov et al. [2018SW001837] (eq 7).
      real(rp), intent(in) :: deltaE
      real(rp), intent(inout) :: eavg,nflx
      real(rp) :: eavgMR,nflxMR

      if(deltaE<=0.0) then
         ! Only modify diffuse precipitation. May modify mono precipitation when a formula is available.
         call AugmentMR(eavg,nflx,eavgMR,nflxMR)
         eavg = eavgMR
         nflx = nflxMR
      endif
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

    subroutine dragonking_auroralsmooth(Np,Nt,vars,aurora)
      integer, intent(in) :: Np,Nt
      type(dkVariables_T), intent(inout) :: vars
      type(dkAurora_T), intent(inout) :: aurora
      integer :: i, j
      logical :: smthDEPonly = .false.
      integer :: smthE
      real(rp) :: tmpC(Np,Nt),tmpD(Np,Nt),tmpE(Np+4,Nt+4),tmpF(Np+4,Nt+4)
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
        where(vars%DeltaE>0)
          aurora%PrecipMask = 0.D0
        end where
        ! back up mono
        tmpC = vars%Avg_Eng
        tmpD = vars%Num_Flux
      endif

      ! get expanded diffuse with margin grid (ghost cells)
      if(smthE==1) then
        call dragonking_margin(Np,Nt,vars%Avg_Eng*vars%Num_Flux*aurora%PrecipMask, tmpE)
      else
        call dragonking_margin(Np,Nt,vars%Avg_Eng*aurora%PrecipMask, tmpE)
      end if
      call dragonking_margin(Np,Nt,vars%Num_Flux*aurora%PrecipMask, tmpF)

      !$OMP PARALLEL DO default(shared) collapse(2) &
      !$OMP private(i,j)
      do j=1,Nt
        do i=1,Np
          ! SmoothOperator55(Q,isGO)
          vars%Num_Flux(i,j)=SmoothOperator55(tmpF(i:i+4,j:j+4))
          if(smthE==1) then
            vars%Avg_Eng(i,j) = max(SmoothOperator55(tmpE(i:i+4,j:j+4))/vars%Num_Flux(i,j),1.D-8)
          else
            vars%Avg_Eng(i,j) = max(SmoothOperator55(tmpE(i:i+4,j:j+4)),1.D-8)
          end if
        end do
      end do

      if(smthDEPonly) then
        ! combine smoothed diffuse and unsmoothed mono
        where(vars%DeltaE>0)
          vars%Avg_Eng  = tmpC
          vars%Num_Flux = tmpD
        end where
      endif
    end subroutine dragonking_auroralsmooth

    subroutine dragonking_ramp(rampFactor,r,rPolarBound,rEquatBound,rLowLimit)
      real(rp), intent(inout) :: rampFactor
      ! type(mixGrid_T), intent(in) :: G      
      real(rp), intent(in) :: r, rPolarBound,rEquatBound,rLowLimit
      
      if (r <= rPolarBound) then
        rampFactor = 1.0D0
      else if ( (r > rPolarBound).and.(r <= rEquatBound) ) then
        rampFactor = 1.0D0+(rLowLimit - 1.0D0)*(r - rPolarBound)/(rEquatBound-rPolarBound)
      else ! r > rEquatBound
        rampFactor = rLowLimit
      end if
    end subroutine dragonking_ramp

    subroutine dragonking_beta_gtype(i,j,vars,aurora,params,isAncG)
      ! Use RCM precipitation and source population to derive the loss cone rate beta.
      ! Assume beta is one in the polar cap and smooth across RCM boundary.
      integer, intent(in) :: i,j
      type(dkVariables_T), intent(inout) :: vars
      type(dkAurora_T), intent(inout) :: aurora
      type(dkParameters_T), intent(in) :: params
      logical, intent(inout) :: isAncG

      real(rp) :: phi0_rcm

      ! Initialize IM_BETA with xml input value.
      ! IM_GTYPE is interpolated from RCM: 1=RCM and 0=MHD.
      vars%Im_Beta(i,j) = params%beta_inp
      if (vars%Im_Gtype(i,j) <= 0.01) isAncG = .true.

      ! Derive beta as the ratio between RCM precipitation and trapped flux.
      ! sqrt([Pa]*[#/m^3]/[kg]) = sqrt([#/m^4/s^2]) = 1e-4*[#/cm^2/s]
      phi0_rcm = sqrt(vars%Im_Epre(i,j)*vars%Im_Eden(i,j)/(Me_cgs*1e-3*2*pi))*1.0e-4
      if(phi0_rcm>TINY) then
         ! This criterion includes all RCM grid. Note beta is 0 where IM_ENFLX is zero.
         vars%Im_Beta(i,j) = vars%Im_Enflx(i,j)/phi0_rcm
         if(vars%Im_Beta(i,j)>1.0) then
            ! Near RCM high lat boundary, trapped flux is not well constrained.
            ! The beta there can be easily >>1 and is thus not as reliable.
            ! Reset grid weight to MHD if RCM beta>1.
            vars%Im_Gtype(i,j) = 0.0 !1.0/St%Vars(i,j,IM_BETA)**2
            vars%Im_Beta(i,j)  = 0.5*(1.0+params%beta_inp)
            isAncG = .false.
         endif
      endif

      aurora%beta_RCM(i,j) = vars%Im_Beta(i,j)
      aurora%gtype_RCM(i,j) = vars%Im_Gtype(i,j)

    end subroutine dragonking_beta_gtype


    subroutine dragonking_smooth(Np,Nt,Q,isAnchor)
      ! Do smoothing window on ReMIX grid quantity
      ! Skip anchor points
      integer, intent(in) :: Np,Nt
      real(rp), intent(inout) :: Q(Np,Nt)
      logical, intent(in) :: isAnchor(Np,Nt)
      real(rp) :: temp(Np,Nt)
      real(rp) :: thres,mad,Ttmp,madij
      integer :: i,j,it,im1,ip1,jm1,jp1,MaxIter

      thres = 0.025
      MaxIter = 15

      call FixPole(Np,Nt,Q)
      do it=1,MaxIter
        mad = 0.D0 ! max abs difference from last iteration.
        temp = Q

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,jm1,jp1,im1,ip1,Ttmp,madij) &
        !$OMP reduction(max:mad)
        do j=1,Nt ! use open BC for lat.
          do i=1,Np ! use periodic BC for lon.
            ! Do not smooth anchor points.
            if(.not. isAnchor(i,j)) then
              jm1 = j-1
              jp1 = j+1
              if (j == 1)     jm1 = 1
              if (j == Nt) jp1 = Nt

              im1 = i-1
              ip1 = i+1
              if (i == 1)     im1 = Np
              if (i == Np) ip1 = 1

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
        call FixPole(Np,Nt,Q)
        if(mad<thres) exit
      enddo
    end subroutine dragonking_smooth

    !Enforce pole condition that all values at the same point are equal
    subroutine FixPole(Np,Nt,Q)
      integer, intent(in) :: Np,Nt
      real(rp), intent(inout) :: Q(Np,Nt)
      real(rp) :: Qavg

      Qavg = sum(Q(:,2))/Np
      Q(:,1) = Qavg
    end subroutine FixPole

    subroutine dragonking_margin(Np,Nt,array,arraymar)
      ! Conductance boundary processing.
      integer, intent(in) :: Np,Nt     
      real(rp),dimension(Np,Nt), intent(in) :: array
      real(rp),dimension(Np+4,Nt+4), intent(out) :: arraymar
      arraymar(3:Np+2,3:Nt+2) = array
      ! reflective BC for lat
      arraymar(3:Np+2,1) = array(:,3)
      arraymar(3:Np+2,2) = array(:,2)
      arraymar(3:Np+2,Nt+3) = array(:,Nt-1)
      arraymar(3:Np+2,Nt+4) = array(:,Nt-2)
      ! periodic BC for lon
      arraymar(1:2,:) = arraymar(Np+1:Np+2,:)
      arraymar(Np+3:Np+4,:) = arraymar(3:4,:)
    end subroutine dragonking_margin

  end module dkaurora
