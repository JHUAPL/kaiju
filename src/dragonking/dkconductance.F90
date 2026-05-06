module dkconductance
  use dkdefs
  use dktypes
  use mixtypes
  use planethelper
  use gcmtypes
!   use gcminterp
  use math
  use euvhelper
  use auroralhelper
  use kai2geo
  
  implicit none

  contains

   !  subroutine conductance_init(conductance,params,grid)
   !    ! Initialize conductance variables. Called in mixmain.F90/subroutine init_mix.
   !    type(dkConductance_T), intent(inout)   :: conductance
   !    type(dkParameters_T) , intent(in)      :: params
   !    type(mixGrid_T)      , intent(in)      :: grid

   !    ! define these module-wide variables so we don't have to pass the Params object to all the conductance functions
   !    ! this is similar to how it was done in MIX
   !    conductance%euv_model_type    = params%euv_model_type
   !    conductance%et_model_type     = params%et_model_type
   !    conductance%sigma_model_type  = params%sigma_model_type
   !    conductance%F107              = params%F107
   !    conductance%pedmin            = params%pedmin
   !    conductance%hallmin           = params%hallmin
   !    conductance%sigma_ratio       = params%sigma_ratio
   !    conductance%ped0              = params%ped0
   !    conductance%const_sigma       = params%const_sigma
   !    conductance%doStarlight       = params%doStarlight
   !    conductance%apply_cap         = params%apply_cap
   !    conductance%doEMA             = params%doEMA
   !    conductance%doET              = params%doET

   !    if (.not. allocated(conductance%zenith)) allocate(conductance%zenith(grid%Np,grid%Nt))
   !    if (.not. allocated(conductance%coszen)) allocate(conductance%coszen(grid%Np,grid%Nt))
   !    if (.not. allocated(conductance%euvSigmaP)) allocate(conductance%euvSigmaP(grid%Np,grid%Nt))
   !    if (.not. allocated(conductance%euvSigmaH)) allocate(conductance%euvSigmaH(grid%Np,grid%Nt))
   !    if (.not. allocated(conductance%deltaSigmaP)) allocate(conductance%deltaSigmaP(grid%Np,grid%Nt))
   !    if (.not. allocated(conductance%deltaSigmaH)) allocate(conductance%deltaSigmaH(grid%Np,grid%Nt))

   !  end subroutine conductance_init

    subroutine conductance_total(conductance,grid,vars,params)
      ! The main subroutine in module mixconductance.
      ! It derives auroral precipitation and conductance caused by both solar EUV and aurora.
      type(dkConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkParameters_T), intent(in) :: params

      real(rp), dimension(:,:), allocatable :: SigH0,SigP0,MagE !Old values of SigH/SigP
      real(rp) :: dT,upTau,dnTau,wAvgU,wAvgD
      integer :: i,j

      !Initialize delta conductances to zero to avoid using uninitialized memory
      conductance%deltaSigmaP = 0.0_rp
      conductance%deltaSigmaH = 0.0_rp

      !Save old Sigs
      allocate(SigH0(grid%Np,grid%Nt))
      allocate(SigP0(grid%Np,grid%Nt))
      SigP0 = vars%SigmaP
      SigH0 = vars%SigmaH
      
      if (params%doGCM) then
         ! Use GCM conductance.
         vars%SigmaP = max(params%pedmin,vars%SigmaPgcm)
         vars%SigmaH = max(params%hallmin,vars%SigmaHgcm)
      else if (params%const_sigma) then
         ! Use constant conductance.
         vars%SigmaP = params%ped0
         vars%SigmaH = 0.D0
      else
         ! Use vector addition of auroral and EUV conductance.
         call conductance_aurora(conductance,grid,vars,params)

         if (params%doET) then
            !Add ET enhancement to embiggen auroral conductance
            allocate(MagE(grid%Np,grid%Nt))
            MagE = 0.0
            call CalcMagE(grid,vars%Pot,MagE) !Get magnitude of E field [mV/m]
            !Apply ET in diffuse regions where auroral conductance is above floor
            where ( (conductance%deltaSigmaP > params%pedmin + TINY) .and. (vars%DeltaE < TINY) )
               conductance%deltaSigmaP = SigET_P(MagE)*conductance%deltaSigmaP
               conductance%deltaSigmaH = SigET_H(MagE)*conductance%deltaSigmaH
            endwhere
            deallocate(MagE)
         endif

         vars%SigmaP = sqrt( conductance%euvSigmaP**2 + conductance%deltaSigmaP**2) 
         vars%SigmaH = sqrt( conductance%euvSigmaH**2 + conductance%deltaSigmaH**2)
         
      endif

      !Before applying cap, do optional exponential smoothing
      if ( (.not. params%const_sigma) .and. params%doEMA ) then
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
               if(vars%SigmaP(i,j) < SigP0(i,j)) then
                  !Local conductance dropping
                  vars%SigmaP(i,j) = wAvgD*vars%SigmaP(i,j) + (1-wAvgD)*SigP0(i,j)
               else
                  !Local conductance increasing
                  vars%SigmaP(i,j) = wAvgU*vars%SigmaP(i,j) + (1-wAvgU)*SigP0(i,j)
               endif
               if(vars%SigmaH(i,j) < SigH0(i,j)) then
                  !Local conductance dropping
                  vars%SigmaH(i,j) = wAvgD*vars%SigmaH(i,j) + (1-wAvgD)*SigH0(i,j)
               else
                  !Local conductance increasing
                  vars%SigmaH(i,j) = wAvgU*vars%SigmaH(i,j) + (1-wAvgU)*SigH0(i,j)
               endif
            enddo
         enddo

      endif
      deallocate(SigP0)
      deallocate(SigH0)

      ! Apply cap
      if ((params%apply_cap).and.(.not. params%const_sigma).and.(.not. params%doGCM)) then
         vars%SigmaP = max(params%pedmin,vars%SigmaP)
         vars%SigmaH = min(max(params%hallmin,vars%SigmaH),vars%SigmaP*params%sigma_ratio)
      endif

    end subroutine conductance_total

    subroutine conductance_euv(conductance,grid,vars,params)
      ! Derive solar EUV conductance: euvSigmaP, euvSigmaH.
      type(dkConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in)         :: grid
      type(dkVariables_T), intent(inout) :: vars
      type(dkParameters_T), intent(in)   :: params

      conductance%coszen = grid%x*cos(vars%tilt)+sqrt(1.-grid%x**2-grid%y**2)*sin(vars%tilt) ! as it should be
      conductance%zenith = acos(conductance%coszen)

      select case ( params%euv_model_type )
         case (LOMPE) 
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
      type(mixGrid_T), intent(in) :: grid
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
            eflux = kev2erg*vars%Avg_Eng(i,j)*vars%Num_Flux(i,j)  ! Energy flux in ergs/cm^2/s
            eavg  = vars%Avg_Eng(i,j) ! Average energy in keV
            if (vars%Aur_Type(i,j) == AT_RMono) then
               dE = vars%DeltaE(i,j) ! Energy width in keV
            else
               dE = 0.0_rp
            endif

            select case ( params%sigma_model_type )
               case (KAE23)
                  ! Kaeppler et al. [2023] model, mono/diffuse specific parameters.
                  conductance%deltaSigmaP(i,j) = SigmaP_Kaep23(eavg,eflux,dE)
                  conductance%deltaSigmaH(i,j) = SigmaH_Kaep23(eavg,eflux,dE)
               case (ROB87)
                  ! Robinson formula
                  conductance%deltaSigmaP(i,j) = SigmaP_Robinson(eavg,eflux)
                  conductance%deltaSigmaH(i,j) = SigmaH_Robinson(eavg,eflux)
               case default
                  stop "The SIGMA model type entered is not supported."
            end select
         enddo
      enddo

    end subroutine conductance_aurora
    
end module dkconductance