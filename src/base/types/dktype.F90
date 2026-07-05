module dktypes
  use kdefs
  use shellGrid

  implicit none

  type dkParams_T
    integer  :: aurora_model_type
    real(rp) :: alpha, beta, R
    logical  :: doRamp, doChill, doMR, doAuroralSmooth

    ! arrays on the grid
    real(rp), dimension(:,:), allocatable :: E0, phi0, deltaE, rampFactor, PrecipMask
  end type dkParams_T

  type auroraState_T
  end type auroraState_T

  type auroraGrid_T
  end type auroraGrid_T
  
  ! used to store an entire instance of MIX (e.g., one per hemisphere)
  type auroraIon_T
!     type(auroraState_T)       :: St
!     type(auroraGrid_T)        :: G       ! G - primary MIX grid used for the solver. 
!     type(dkParams_T)          :: P
  end type auroraIon_T

  ! used to store all instances of mixIon type, i.e., all hemispheres
  type auroraApp_T
     type(auroraIon_T), dimension(:), allocatable :: dragonking
  end type auroraApp_T

  ! use this as a container to store the variables read from a previous H5 file
  type auroraIO_T
  end type auroraIO_T

  type precip_t  
      integer :: nVars ! change together with the enumerator above
      !type(ShellGridVar_T),target :: FAC,SIGMAP,SIGMAH,SOUND_SPEED,DENSITY,AVG_ENG,NUM_FLUX,NEUTRAL_WIND,EFIELD,IM_EAVG,IM_EFLUX,IM_IAVG,IM_IFLUX,&
      !  Z_EAVG,Z_NFLUX,CRPOT,TPOT,IM_GTYPE,AUR_TYPE,IM_BETA,IM_EDEN,IM_EPRE,IM_ENFLX,IM_INFLX,DELTAE,IM_NPSP,IM_THCON,IM_PHCON
      !type(ShellGridVar_T), pointer :: data(:)
      type(ShellGridVar_T), allocatable :: data(:)
  end type precip_T

end module dktypes

