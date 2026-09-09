module dkaurora
  use kdefs
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

  real(rp), parameter, private :: Ne_floor = 0.03e6 ! minimum Ne in [/m^3] when evaluating the linearized FL relation.
  real(rp), parameter, private :: Ne_psp = 1.0e6 ! Ne threshold for the plasmasphere in [/m^3].
  real(rp), parameter, private :: GuABNF = 1.e7 ! Gussenhoven+[1983] Auroral Boundary Number Flux in [#/cm^s/s].
  ! real(rp), private :: RinMHD = 0.0 !Rin of MHD grid (0 if not running w/ MHD)

  private :: hemiSignJ

  contains

    subroutine hemiSignJ(grid,sgnJ)
      !! Sign convention for the field-aligned current. Precomputed for performance
      type(dkGrid_T), intent(in) :: grid
      real(rp), intent(out) :: sgnJ(grid%Nt)
      integer :: j

      do j=1,grid%Nt
         sgnJ(j) = -sign(1.0_rp, cos(grid%thc(j)))
      enddo
    end subroutine hemiSignJ

    subroutine precipitation_total(aurora,grid,vars,conductance,params)
      ! The main subroutine in module dragonking.
      ! It derives auroral precipitation.
      type(dkAurora_T), intent(inout) :: aurora
      type(dkGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkConductance_T), intent(inout) :: conductance
      type(dkModel_T), intent(in) :: params

      select case ( params%aurora_model_type )
        case (DKOpts%FEDDER)
          call dragonking_fedder95(aurora,grid,vars,conductance,params)
        case (DKOpts%SUNNY)
          call dragonking_sunny(vars)
        case (DKOpts%LINMRG)
          call dragonking_linmrg(aurora,grid,vars,params)
        case default
          stop "The auroral precipitation model type entered is not supported."
      end select

      ! Correct for multiple reflections if you're so inclined
      if (params%doMR) call dragonking_mr(aurora,grid,vars,params)

      ! Smooth precipitation energy and flux before calculating conductance.
      if (params%doAuroralSmooth) call dragonking_auroralsmooth(vars,grid,aurora)

    end subroutine precipitation_total

    subroutine dragonking_fedder95(aurora,grid,vars,conductance,params)
      ! Derive electron precipitation energy flux and avg energy using Fedder95 formula.
      type(dkAurora_T), intent(inout) :: aurora
      type(dkGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkConductance_T), intent(inout) :: conductance
      type(dkModel_T), intent(in) :: params

      real(rp) :: sgnJ(grid%Nt)
      real(rp) :: Rout = 6.D0, Rin = 1.2D0
      real(rp) :: rhoFactor = 3.3D-24*0.5D0
      real(rp) :: rPolarBound = 20.0D0*pi/180.D0, rEquatBound = 30.0D0*pi/180.0D0, rLowLimit = 0.02D0
      real(rp) :: D,Cs,Pe,Ne,dE,phi0,E0,aRes,dEc
      integer :: i,j

      !Default to no precip (diffuse conductance branch)
      vars%Aur_Type%data = DKOpts%AT_NoPre

      ! fills in rampFactor
      if (params%doRamp) then
         call dragonking_ramp(aurora,grid,rPolarBound,rEquatBound,rLowLimit)
      else
         aurora%rampFactor = 1.0D0
      endif

      call hemiSignJ(grid,sgnJ)

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,D,Cs,dE,phi0,E0,aRes,dEc)
      do i=1,grid%Np
         do j=1,grid%Nt

            if (params%doChill) then
               ! MHD density replaced with gallagher where it's lower
               ! and temperature changed correspondingly
               D  = max(grid%D0(j,i)*Mp_cgs, vars%Density%data(j,i)) ! [g/cm^3]
               Cs = vars%Sound_Speed%data(j,i)*sqrt(vars%Density%data(j,i)/D) ! [cm/s]
            else
               D  = vars%Density%data(j,i)
               Cs = vars%Sound_Speed%data(j,i)
            endif
            ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
            E0 = params%alpha*Mp_cgs*heFrac*erg2kev*Cs**2*aurora%rampFactor(j,i)
            aurora%E0  (j,i) = E0
            ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
            phi0 = sqrt(kev2erg)/(heFrac*Mp_cgs)**1.5D0*params%beta*D*sqrt(E0)*aurora%rampFactor(j,i)
            aurora%phi0(j,i) = phi0

            ! resistence out of the ionosphere is 2*rout resistence into the
            ! ionosphere is 2*rin outward current is positive
            if( sgnJ(j)*vars%FAC%data(j,i) >=0. ) then
              aRes = 2.D0*Rout
            else
              aRes = 2.D0*Rin
            endif

            ! Density floor to limit characteristic energy.  See Wiltberger et al. 2009 for details.
            D = min(D,rhoFactor*conductance%euvSigmaP(j,i))
            dE = (heFrac*Mp_cgs)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*params%R*aRes*sgnJ(j)*(vars%FAC%data(j,i)*1.D-6)*sqrt(E0)/D
            aurora%deltaE(j,i) = dE
            dEc = min(dE,maxDrop)

            vars%Avg_Eng%data(j,i) = max( E0+dEc, eTINY)
            if(dEc>0.) then
              vars%Num_Flux%data(j,i) = phi0*(8.D0-7.D0*exp(-dEc/7.D0/E0))
            else
              vars%Num_Flux%data(j,i) = phi0*exp(dEc/E0)
            endif
            vars%DeltaE%data(j,i) = aurora%deltaE(j,i) ! [kV]
         enddo ! j (theta)
      enddo ! i (phi)

    end subroutine dragonking_fedder95

    subroutine dragonking_sunny(vars)
      ! Assign zero precipitation when only EUV conductance is used.
      type(dkVariables_T), intent(inout) :: vars

      ! Make sure zero Eavg/Nflux does not cause issues in conductance calculation.
      vars%Aur_Type%data = DKOpts%AT_NoPre
      vars%Avg_Eng%data = 0.D0
      vars%Num_Flux%data = 0.D0
      vars%DeltaE%data = 0.D0

    end subroutine dragonking_sunny

    subroutine dragonking_linmrg(aurora,grid,vars,params)
      ! Derive mono-diffuse electron precipitation where mono is based on linearized FL relation,
      ! and diffuse is a combination of MHD and RCM precipitation.
      type(dkAurora_T), intent(inout) :: aurora
      type(dkGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkModel_T), intent(in) :: params

      integer :: i,j
      real(rp) :: wRCM,wMHD
      real(rp) :: mhd_nflx,mhd_eflx,mhd_eavg
      real(rp) :: rcm_nflx,rcm_eflx,rcm_eavg
      real(rp) :: mix_nflx,mhd_SigPH,rcm_SigPH
      logical :: isMono,isIMAG

      vars%Aur_Type%data = 0

      ! derive spatially varying beta and gtype
      ! using RCM precipitation and thermal fluxes.
      call dragonking_beta_gtype(aurora,grid,vars,params)

      ! Derive mono using the linearized FL relation.
      call dragonking_linmono(aurora,grid,vars,params)

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,isMono,isIMAG,wRCM,wMHD) &
      !$OMP private(mhd_nflx,mhd_eflx,mhd_eavg) &
      !$OMP private(rcm_nflx,rcm_eflx,rcm_eavg) &
      !$OMP private(mix_nflx,mhd_SigPH,rcm_SigPH)
      do i=1,grid%Np
         do j=1,grid%Nt
            !Start by figuring out where we are
!            isPSP = vars%Im_Npsp%data(j,i)>=Ne_psp .or. aurora%beta_RCM(j,i)<=0.01
            isIMAG = vars%Im_Gtype%data(j,i)>=0.5 ! defined in raiCplTypesSub.F90: "-1=Inactive, 0=Buffer, 1=Active" -> 0, 0.5, and 1.0.
            isMono = vars%DeltaE%data(j,i)>eTINY
            !Grab values
            mhd_nflx = vars%Num_Flux%data(j,i)
            if(mhd_nflx>GuABNF) then
               mhd_eavg = vars%Avg_Eng%data(j,i)
            else
               mhd_eavg = eTINY
            endif
            mhd_eflx = mhd_eavg*mhd_nflx*kev2erg

            rcm_nflx = vars%Im_Enflx%data(j,i)
            if(rcm_nflx>GuABNF) then
               rcm_eavg = vars%Im_Eflux%data(j,i)/(rcm_nflx*kev2erg)
            else
               rcm_eavg = eTINY
            endif
            rcm_eflx = rcm_eavg*rcm_nflx*kev2erg

            if(isMono) then
               vars%Num_Flux%data(j,i) = mhd_nflx
               vars%Avg_Eng%data(j,i)  = mhd_eavg
               vars%Aur_Type%data(j,i) = DKOpts%AT_RMono
               if(isIMAG) then
                  ! Set auroral type to mono if mono is above Gussenhove threshold and gives higher Robinson L2.
                  ! Else set diffuse or no precipitation depending on rcm values relative to the threshold.
                  mhd_SigPH = SigmaP_Robinson(mhd_eavg,mhd_eflx)**2+SigmaH_Robinson(mhd_eavg,mhd_eflx)**2
                  rcm_SigPH = SigmaP_Robinson(rcm_eavg,rcm_eflx)**2+SigmaH_Robinson(rcm_eavg,rcm_eflx)**2
                  if(mhd_SigPH<rcm_SigPH) then
                     vars%Num_Flux%data(j,i) = rcm_nflx
                     vars%Avg_Eng%data(j,i)  = rcm_eavg
                     vars%Aur_Type%data(j,i) = DKOpts%AT_RMfnE
                  endif
               endif
            else
               ! Linearly merge MHD and RCM diffuse nflux and eflux.
               ! Note where deltaE>eTINY but aurora%beta_RCM<=0.01, gtype will be near 1.
               wRCM = aurora%gtype_RCM(j,i)
               wMHD = 1.0-wRCM
               mix_nflx = wMHD*mhd_nflx + wRCM*rcm_nflx
               vars%Num_Flux%data(j,i) = mix_nflx
               if(mix_nflx>GuABNF) then
                  vars%Avg_Eng%data(j,i) = (wMHD*mhd_eflx+wRCM*rcm_eflx)/(mix_nflx*kev2erg)
                  vars%Aur_Type%data(j,i)= DKOpts%AT_RMnoE
               else
                  vars%Avg_Eng%data(j,i) = eTINY
                  vars%Aur_Type%data(j,i)= DKOpts%AT_NoPre
               endif
            endif
         enddo
      enddo

    end subroutine dragonking_linmrg

    subroutine dragonking_linmono(aurora,grid,vars,params)
      ! Linearized Fridman-Lemaire relation for mono.
      type(dkAurora_T), intent(inout) :: aurora
      type(dkGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkModel_T), intent(in) :: params
      
      real(rp) :: sgnJ(grid%Nt)
      integer :: i,j
      real(rp) :: D,Cs,Pe,Ne,J2eF0,eV2kT,dE,phi0,kT
      real(rp) :: Pe_mhd,Ne_mhd,kT_mhd,phi0_mhd
      real(rp) :: Pe_rcm,Ne_rcm,kT_rcm,phi0_rcm
      real(rp) :: wMHD,wRCM
      real(rp) :: alpha_loc,beta_loc,gamma_loc
      real(rp) :: dens_loc,csnd_loc,d0_loc
      logical  :: hasMono

      gamma_loc = params%gamma

      call hemiSignJ(grid,sgnJ)

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,D,Cs,Pe,Ne,J2eF0,eV2kT,dE,phi0,kT) &
      !$OMP private(Pe_mhd,Ne_mhd,kT_mhd,phi0_mhd) &
      !$OMP private(Pe_rcm,Ne_rcm,kT_rcm,phi0_rcm) &
      !$OMP private(wMHD, wRCM, hasMono) &
      !$OMP private(alpha_loc,beta_loc,dens_loc,csnd_loc,d0_loc)
      do i=1,grid%Np
         do j=1,grid%Nt
            ! Copy struct fields to local scalars to force consistent
            ! FP evaluation order across module structure changes
            alpha_loc = aurora%alpha_RCM(j,i)
            beta_loc  = aurora%beta_RCM(j,i)
            dens_loc  = vars%Density%data(j,i)
            csnd_loc  = vars%Sound_Speed%data(j,i)
            d0_loc    = grid%D0(j,i)

            if (params%doChill) then
               ! MHD density replaced with gallagher where it's lower
               ! and temperature changed correspondingly
               D  = max(d0_loc*Mp_cgs,dens_loc) ! [g/cm^3]
               Cs = csnd_loc*sqrt(dens_loc/D) ! [cm/s]
            else
               D  = dens_loc
               Cs = csnd_loc
            endif
            ! electron pressure from MHD side in [Pa]
            ! 0.1 is to convert [g/cm^3]*[cm/s]^2=[g/cm/s^2] to [Pa].
            Pe_mhd = 0.1/gamma_loc*alpha_loc*D*(Cs**2)
            ! electron number density from MHD side in [/m^3].
            Ne_mhd = D/(Mp_cgs*heFrac)*1.0D6
            ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
            kT_mhd = Pe_mhd/(Ne_mhd*kev2J)
            ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
            ! Treat MHD flux on closed field lines as trapped flux.
            ! Beta now stands for the ratio between precipitation and trapped.
            ! Note phi0 can be zero in the polar cap.
            phi0_mhd = beta_loc*sqrt(Pe_mhd*Ne_mhd/(2.0D-3*PI*Me_cgs))*1.0D-4

            ! RCM flux
            Pe_rcm = vars%Im_Epre%data(j,i)
            Ne_rcm = vars%Im_Eden%data(j,i)
            phi0_rcm = vars%Im_Enflx%data(j,i)

            ! Merged flux
            wRCM = aurora%gtype_RCM(j,i)
            wMHD = 1.0-wRCM
            phi0 = wMHD*phi0_mhd + wRCM*phi0_rcm
            Ne   = wMHD*Ne_mhd + wRCM*Ne_rcm
            Pe   = wMHD*Pe_mhd + wRCM*Pe_rcm
            if(Ne>TINY) then
               kT = Pe/(Ne*kev2J)
            else
               kT = eTINY
            endif

            aurora%E0  (j,i) = 2.0*kT
            hasMono = .false.
            aurora%phi0(j,i) = phi0
            ! phi0 = 0 in the polar cap; guard the divide
            if (phi0 > TINY) then
               J2eF0 = sgnJ(j)*(vars%FAC%data(j,i)*1.0D-6)/(eCharge*phi0*1.0D4)
               ! Linearize the Fridman-Lemaire relation:
               ! eV = kB*Te*(RM-1)*ln((RM-1)/(RM-J/e/F0)) when 1<=J/e/F0<=RM.
               ! eV ~ kB*Te*(J/e/F0-1) when RM>>J/e/F0
               if (J2eF0>1.0 .and. Ne>Ne_floor .and. vars%Im_Npsp%data(j,i)<Ne_psp) hasMono = .true.
            endif
            if (hasMono) then
               dE = kT*(J2eF0-1.0)
               aurora%deltaE(j,i) = dE
               eV2kT = min(dE,maxDrop)/kT + 1.0
               vars%Z_Nflux%data(j,i) = phi0*eV2kT
               vars%Z_Eavg%data(j,i) = kT*(eV2kT+1.0/eV2kT)
            else
               ! When there is no acceleration, MHD flux represents the diffuse precipitation.
               aurora%deltaE(j,i) = 0.0
               vars%Z_Nflux%data(j,i) = phi0
               vars%Z_Eavg%data(j,i) = 2.0*kT
            endif

            vars%DeltaE%data(j,i)   = aurora%deltaE(j,i)
            vars%Avg_Eng%data(j,i)  = max(vars%Z_Eavg%data(j,i),eTINY)   ! [keV]
            vars%Num_Flux%data(j,i) = vars%Z_Nflux%data(j,i)  ! [#/cm^2/s]

         enddo ! j (theta)
      enddo ! i (phi)
    end subroutine dragonking_linmono

    subroutine dragonking_mr(aurora,grid,vars,params)
      ! George Khazanov's multiple reflection(MR) corrections
      ! Modify the diffuse precipitation energy and mean energy based on equation (5) in Khazanov et al. [2019JA026589]
      ! Kc = 3.36-exp(0.597-0.37*Eavg+0.00794*Eavg^2)
      ! Eavg_c = 0.073+0.933*Eavg-0.0092*Eavg^2
      ! Nflx_c = Eflx_c/Eavg_c = Eflx*Kc/Eavg_c = Nflx*Eavg*Kc/Eavg_c
      ! where Eavg is the mean energy in [keV] before MR, Eavg_c is the modified mean energy. 
      ! Kc is the ratio between modified enflux and unmodified enflux.
      ! No longer use the conductance modification factor in Khazanov et al. [2018SW001837] (eq 7).
      type(dkAurora_T), intent(in) :: aurora
      type(dkGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkModel_T), intent(in) :: params
      real(rp) :: eavg,nflx,eavgMR,nflxMR
      integer :: i,j

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,eavg,nflx,eavgMR,nflxMR)
      do i=1,grid%Np
        do j=1,grid%Nt
          if(aurora%deltaE(j,i)<=0.0) then
            ! Only modify diffuse precipitation. May modify mono precipitation when a formula is available.
            eavg = vars%Avg_Eng%data(j,i)
            nflx = vars%Num_Flux%data(j,i)
            call AugmentMR(eavg,nflx,eavgMR,nflxMR)
            vars%Avg_Eng%data(j,i) = eavgMR
            vars%Num_Flux%data(j,i) = nflxMR
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

    subroutine dragonking_auroralsmooth(vars,grid,aurora)
      type(dkVariables_T), intent(inout) :: vars
      type(dkGrid_T), intent(in) :: grid
      type(dkAurora_T), intent(inout) :: aurora
      integer :: i, j
      logical :: smthDEPonly = .false.
      integer :: smthE
      real(rp) :: tmpC(grid%Nt,grid%Np),tmpD(grid%Nt,grid%Np),tmpE(grid%Nt+4,grid%Np+4),tmpF(grid%Nt+4,grid%Np+4)
      smthE = 1 ! 1. smooth EFLUX; 2. smooth EAVG

      tmpC = 0.D0
      tmpD = 0.D0
      tmpE = 0.D0
      tmpF = 0.D0

      ! Test only smoothing diffuse precipitation.
      ! Diffuse mask is Z_NFLUX<0.
      ! Get AVG_ENG*mask and NUM_FLUX*mask for smoothing.
      ! Fill AVG_ENG and NUM_FLUX with smoothed where mask is true.
      ! this domain of conductance is never used. Here it's used for diffuse precipitation mask. Be careful if not in RCMONO mode.
      aurora%PrecipMask = 1.D0
      if(smthDEPonly) then
        where(vars%DeltaE%data>0)
          aurora%PrecipMask = 0.D0
        end where
        ! back up mono
        tmpC = vars%Avg_Eng%data
        tmpD = vars%Num_Flux%data
      endif

      ! get expanded diffuse with margin grid (ghost cells)
      if(smthE==1) then
        call dragonking_margin(grid,vars%Avg_Eng%data*vars%Num_Flux%data*aurora%PrecipMask, tmpE)
      else
        call dragonking_margin(grid,vars%Avg_Eng%data*aurora%PrecipMask, tmpE)
      end if
      call dragonking_margin(grid,vars%Num_Flux%data*aurora%PrecipMask, tmpF)

      !$OMP PARALLEL DO default(shared) collapse(2) &
      !$OMP private(i,j)
      do i=1,grid%Np
        do j=1,grid%Nt
          ! SmoothOperator55(Q,isGO)
          vars%Num_Flux%data(j,i)=SmoothOperator55(tmpF(j:j+4,i:i+4))
          if(smthE==1) then
            vars%Avg_Eng%data(j,i) = max(SmoothOperator55(tmpE(j:j+4,i:i+4))/max(vars%Num_Flux%data(j,i), eTINY),1.D-8)
          else
            vars%Avg_Eng%data(j,i) = max(SmoothOperator55(tmpE(j:j+4,i:i+4)),1.D-8)
          end if
        end do
      end do

      if(smthDEPonly) then
        ! combine smoothed diffuse and unsmoothed mono
        where(vars%DeltaE%data>0)
          vars%Avg_Eng%data  = tmpC
          vars%Num_Flux%data = tmpD
        end where
      endif
    end subroutine dragonking_auroralsmooth

    subroutine dragonking_ramp(aurora,grid,rPolarBound,rEquatBound,rLowLimit)
      type(dkAurora_T), intent(inout) :: aurora
      type(dkGrid_T), intent(in) :: grid
      real(rp), intent(in) :: rPolarBound,rEquatBound,rLowLimit
      
      where (grid%r <= rPolarBound)
        aurora%rampFactor = 1.0D0
      elsewhere ( (grid%r > rPolarBound).and.(grid%r <= rEquatBound) )
        aurora%rampFactor = 1.0D0+(rLowLimit - 1.0D0)*(grid%r - rPolarBound)/(rEquatBound-rPolarBound)
      elsewhere ! grid%r > rEquatBound
        aurora%rampFactor = rLowLimit
      end where
    end subroutine dragonking_ramp

    subroutine dragonking_beta_gtype(aurora,grid,vars,params)
      ! Use RCM precipitation and source population to derive the loss cone rate beta.
      ! Assume beta is one in the polar cap and smooth across RCM boundary.
      type(dkAurora_T), intent(inout) :: aurora
      type(dkGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkModel_T), intent(in) :: params

      logical :: isAncB(grid%Nt,grid%Np),isAncG(grid%Nt,grid%Np)
      real(rp) :: phi0_rcm
      integer :: i,j

      ! Initialize IM_BETA with xml input value.
      ! IM_GTYPE is interpolated from RCM: 1=RCM and 0=MHD.
      vars%Im_Beta%data = params%beta_inp
      isAncB = .false. ! beta is smoothed everywhere.
      isAncG = .false.

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,phi0_rcm)
      do i=1,grid%Np
         do j=1,grid%Nt
            ! Set grids outside RCM as anchors that won't be smoothed.
            if(vars%Im_Gtype%data(j,i)<=0.01) then
               isAncG(j,i) = .true.
            endif

            ! Derive beta as the ratio between RCM precipitation and trapped flux.
            ! sqrt([Pa]*[#/m^3]/[kg]) = sqrt([#/m^4/s^2]) = 1e-4*[#/cm^2/s]
            phi0_rcm = sqrt(vars%Im_Epre%data(j,i)*vars%Im_Eden%data(j,i)/(Me_cgs*1e-3*2*pi))*1.0e-4
            if(phi0_rcm>TINY) then
               ! This criterion includes all RCM grid. Note beta is 0 where IM_ENFLX is zero.
               vars%Im_Beta%data(j,i) = vars%Im_Enflx%data(j,i)/phi0_rcm
               !if(vars%Im_Beta%data(j,i)>1.0) then
                  ! Near RCM high lat boundary, trapped flux is not well constrained.
                  ! The beta there can be easily >>1 and is thus not as reliable.
                  ! Reset grid weight to MHD if RCM beta>1.
               !   vars%Im_Gtype%data(j,i) = 0.0 !1.0/vars%Im_Beta%data(j,i)**2
               !   vars%Im_Beta%data(j,i)  = 0.5*(1.0+params%beta_inp)
               !   isAncG(j,i) = .false.
               !endif
            endif
         enddo
      enddo

      ! Smooth IM_BETA and save in aurora%beta_RCM
      aurora%beta_RCM = vars%Im_Beta%data
      call dragonking_smooth(grid,aurora%beta_RCM,isAncB)
      ! IM_BETA keeps the unclamped smoothed value; beta_RCM is clamped to [0,1]
      vars%Im_Beta%data = aurora%beta_RCM
      aurora%beta_RCM  = min(max(aurora%beta_RCM,0.D0),1.D0)

      ! Smooth IM_GTYPE and save in gtype_RCM
      aurora%gtype_RCM = vars%Im_Gtype%data
      call dragonking_smooth(grid,aurora%gtype_RCM,isAncG)
      ! IM_GTYPE keeps the unclamped smoothed value; gtype_RCM is clamped to [0,1]
      vars%Im_Gtype%data = aurora%gtype_RCM
      aurora%gtype_RCM  = min(max(aurora%gtype_RCM,0.D0),1.D0)
    end subroutine dragonking_beta_gtype

    subroutine dragonking_smooth(grid,Q,isAnchor)
      ! Do smoothing window on ReMIX grid quantity
      ! Skip anchor points
      type(dkGrid_T), intent(in) :: grid
      real(rp), intent(inout) :: Q(grid%Nt,grid%Np)
      logical, intent(in) :: isAnchor(grid%Nt,grid%Np)
      real(rp) :: temp(grid%Nt+4,grid%Np+4)
        !Haloed copy of Q: cross-pole in theta, periodic in phi
      real(rp) :: thres,mad,Ttmp,madij
      integer :: i,j,ih,jh,it,MaxIter

      thres = 0.025
      MaxIter = 15

      do it=1,MaxIter
        mad = 0.D0 ! max abs difference from last iteration.
        ! dragonking_margin supplies both BCs the stencil needs.
        ! On a cell-centered grid the neighbour across the pole is the same
        ! near-pole ring at antipodal longitude
        call dragonking_margin(grid,Q,temp)

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,ih,jh,Ttmp,madij) &
        !$OMP reduction(max:mad)
        do i=1,grid%Np
          ih = i+2   ! interior of the haloed array starts at index 3
          do j=1,grid%Nt
            ! Do not smooth anchor points.
            if(.not. isAnchor(j,i)) then
              jh = j+2

              Ttmp =(temp(jh-1,ih-1)+temp(jh,ih-1)+temp(jh+1,ih-1) &
                   + temp(jh-1,ih  )+temp(jh,ih  )+temp(jh+1,ih  ) &
                   + temp(jh-1,ih+1)+temp(jh,ih+1)+temp(jh+1,ih+1))/9.D0
              madij = abs(Q(j,i)-Ttmp)
              Q(j,i) = Ttmp
            else
              madij = 0.D0
            endif
            mad = max(mad,madij)
          enddo
        enddo
        if(mad<thres) exit
      enddo
    end subroutine dragonking_smooth

    subroutine dragonking_margin(grid,array,arraymar)
      ! Conductance boundary processing.
      type(dkGrid_T), intent(in) :: grid
      real(rp),dimension(grid%Nt,grid%Np), intent(in) :: array
      real(rp),dimension(grid%Nt+4,grid%Np+4), intent(out) :: arraymar

      if (mod(grid%Np,2) /= 0) then
         write(*,*) "dragonking_margin: cross-pole BC needs even Np, got ", grid%Np
         stop
      endif

      arraymar(3:grid%Nt+2,3:grid%Np+2) = array
      ! Cross-pole lat BC: neighbor is the near-pole cell at antipodal lon
      arraymar(2,3:grid%Np+2)         = cshift(array(1,:),         grid%Np/2)
      arraymar(1,3:grid%Np+2)         = cshift(array(2,:),         grid%Np/2)
      arraymar(grid%Nt+3,3:grid%Np+2) = cshift(array(grid%Nt,:),   grid%Np/2)
      arraymar(grid%Nt+4,3:grid%Np+2) = cshift(array(grid%Nt-1,:), grid%Np/2)
      ! periodic BC for lon (phi, dim2)
      arraymar(:,1:2) = arraymar(:,grid%Np+1:grid%Np+2)
      arraymar(:,grid%Np+3:grid%Np+4) = arraymar(:,3:4)
    end subroutine dragonking_margin

  end module dkaurora
