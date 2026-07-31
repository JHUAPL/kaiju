module dkconductance
  use dktypes
  use planethelper
  use math
  use euvhelper
  use auroralhelper
  use kai2geo
  
  implicit none

  contains

    subroutine conductance_total(conductance,grid,vars,params)
      ! It derives auroral precipitation and conductance caused by both solar EUV and aurora.
      type(dkConductance_T), intent(inout) :: conductance
      type(dkGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkParameters_T), intent(in) :: params

      real(rp), dimension(:,:), allocatable :: SigH0,SigP0,MagE !Old values of SigH/SigP
      real(rp) :: dT,upTau,dnTau,wAvgU,wAvgD
      integer :: i,j

      !Initialize delta conductances to zero to avoid using uninitialized memory
      conductance%deltaSigmaP = 0.0_rp
      conductance%deltaSigmaH = 0.0_rp

      !Save old Sigs
      allocate(SigH0(grid%Nt,grid%Np))
      allocate(SigP0(grid%Nt,grid%Np))
      SigP0 = vars%SigmaP%data
      SigH0 = vars%SigmaH%data

      if (conductance%const_sigma) then
         ! Use constant conductance.
         vars%SigmaP%data = params%ped0
         vars%SigmaH%data = 0.D0
      else
         ! Use vector addition of auroral and EUV conductance.
         call conductance_aurora(conductance,grid,vars,params)

         if (params%doET) then
            !Add ET enhancement to embiggen auroral conductance
            allocate(MagE(grid%Nt,grid%Np))
            MagE = 0.0
            call CalcMagE(grid,vars%Pot%data,MagE) !Get magnitude of E field [mV/m]
            !Apply ET in diffuse regions where auroral conductance is above floor
            where ( (conductance%deltaSigmaP > params%pedmin + TINY) .and. (vars%DeltaE%data < TINY) )
               conductance%deltaSigmaP = SigET_P(MagE)*conductance%deltaSigmaP
               conductance%deltaSigmaH = SigET_H(MagE)*conductance%deltaSigmaH
            endwhere
            deallocate(MagE)
         endif

         vars%SigmaP%data = sqrt( conductance%euvSigmaP**2 + conductance%deltaSigmaP**2)
         vars%SigmaH%data = sqrt( conductance%euvSigmaH**2 + conductance%deltaSigmaH**2)

      endif

      !Before applying cap, do optional exponential smoothing
      if ( (.not. conductance%const_sigma) .and. params%doEMA ) then
         !Want 95% of weight to come from last tau seconds
         !Lazily setting values here
         dT = 5.D0

         dnTau = 30.D0 ![s], fall timescale (eg recombination)
         upTau = 30.D0  ![s], increase timescale (eg state averaging)

         wAvgD = 1.D0 - exp(-3.D0*dT/max(dT,dnTau))
         wAvgU = 1.D0 - exp(-3.D0*dT/max(dT,upTau))
         !Throttle how fast conductance drops (ie lazy recombination timescale)
         !$OMP PARALLEL DO default(shared) &
         !$OMP private(i,j)
         do j=1,grid%Nt
            do i=1,grid%Np
               if(vars%SigmaP%data(j,i) < SigP0(j,i)) then
                  !Local conductance dropping
                  vars%SigmaP%data(j,i) = wAvgD*vars%SigmaP%data(j,i) + (1-wAvgD)*SigP0(j,i)
               else
                  !Local conductance increasing
                  vars%SigmaP%data(j,i) = wAvgU*vars%SigmaP%data(j,i) + (1-wAvgU)*SigP0(j,i)
               endif
               if(vars%SigmaH%data(j,i) < SigH0(j,i)) then
                  !Local conductance dropping
                  vars%SigmaH%data(j,i) = wAvgD*vars%SigmaH%data(j,i) + (1-wAvgD)*SigH0(j,i)
               else
                  !Local conductance increasing
                  vars%SigmaH%data(j,i) = wAvgU*vars%SigmaH%data(j,i) + (1-wAvgU)*SigH0(j,i)
               endif
            enddo
         enddo

      endif
      deallocate(SigP0)
      deallocate(SigH0)

      ! Apply cap
      if ((params%apply_cap).and.(.not. conductance%const_sigma)) then
         vars%SigmaP%data = max(params%pedmin,vars%SigmaP%data)
         vars%SigmaH%data = min(max(params%hallmin,vars%SigmaH%data),vars%SigmaP%data*params%sigma_ratio)
      endif

    end subroutine conductance_total

    subroutine conductance_euv(conductance,grid,vars,params)
      ! Derive solar EUV conductance: euvSigmaP, euvSigmaH.
      type(dkConductance_T), intent(inout) :: conductance
      type(dkGrid_T), intent(in)         :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkParameters_T), intent(in)   :: params

      integer :: j_euv
      ! cos(zenith) = x*cos(tilt) + z*sin(tilt); z must be SIGNED cos(theta) on the full sphere
      do j_euv = 1, grid%Nt
        conductance%coszen(j_euv,:) = grid%x(j_euv,:)*cos(vars%tilt) + &
            cos(grid%thc(j_euv))*sin(vars%tilt)
      end do
      conductance%zenith = acos(conductance%coszen)

      select case ( params%euv_model_type )
         case (DKModel%LOMPE)
            conductance%euvSigmaP = SigP_EUV_LOMPE(conductance%zenith,params%f107)
            conductance%euvSigmaH = SigH_EUV_LOMPE(conductance%zenith,params%f107)

         case default
            stop "The EUV model type entered is not supported."
      end select

      if (params%doStarlight) then ! Makes sense to turn off apply_cap if this is on
         ! add background conductance quadratically instead of a sharp cutoff
         ! first need to remove negative conductances resulting form AMIE for strong tilt
         ! otherwise, the quadratic addition results in conductances growing toward nightside.
         conductance%euvSigmaP = max(conductance%euvSigmaP,0.)
         conductance%euvSigmaH = max(conductance%euvSigmaH,0.)
         
         conductance%euvSigmaP = sqrt(conductance%euvSigmaP**2 + params%pedmin**2)
         conductance%euvSigmaH = sqrt(conductance%euvSigmaH**2 + params%hallmin**2)
      else
         ! otherwise, default to the standard way of applying the floor
         ! I'll leave it as default for now (doStarlight=.false.) but we can reconsider later
         conductance%euvSigmaP = max(conductance%euvSigmaP,params%pedmin)
         conductance%euvSigmaH = max(conductance%euvSigmaH,params%hallmin)
      end if
      
    end subroutine conductance_euv

    subroutine conductance_aurora(conductance,grid,vars,params)
      ! Use Robinson formula to get Pedersen and Hall conductance from electron precipitation.
      type(dkConductance_T), intent(inout) :: conductance
      type(dkGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkParameters_T), intent(in) :: params

      integer :: i,j
      real(rp) :: eavg,eflux,dE
!      real(rp), dimension(:,:), allocatable :: engFlux
!      allocate(engFlux(G%Np,G%Nt))
!      engFlux = kev2erg*St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)  ! Energy flux in ergs/cm^2/s

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,eflux,eavg,dE)
      do j=1,grid%Nt
         do i=1,grid%Np
            eflux = kev2erg*vars%Avg_Eng%data(j,i)*vars%Num_Flux%data(j,i)  ! Energy flux in ergs/cm^2/s
            eavg  = vars%Avg_Eng%data(j,i) ! Average energy in keV
            if (vars%Aur_Type%data(j,i) == DKModel%AT_RMono) then
               dE = vars%DeltaE%data(j,i) ! Energy width in keV
            else
               dE = 0.0_rp
            endif

            select case ( params%sigma_model_type )
               case (DKModel%KAE23)
                  ! Kaeppler et al. [2023] model, mono/diffuse specific parameters.
                  conductance%deltaSigmaP(j,i) = SigmaP_Kaep23(eavg,eflux,dE)
                  conductance%deltaSigmaH(j,i) = SigmaH_Kaep23(eavg,eflux,dE)
               case (DKModel%ROB87)
                  ! Robinson formula
                  conductance%deltaSigmaP(j,i) = SigmaP_Robinson(eavg,eflux)
                  conductance%deltaSigmaH(j,i) = SigmaH_Robinson(eavg,eflux)
               case default
                  stop "The SIGMA model type entered is not supported."
            end select
         enddo
      enddo

    end subroutine conductance_aurora
    
end module dkconductance