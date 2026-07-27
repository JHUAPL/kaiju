module dktypes
  use kdefs
  use shellGrid
  use dkdefs
  use basetypes

  implicit none

  type dkAurora_T
    ! arrays on the grid
    real(rp), dimension(:,:), allocatable :: E0, phi0, deltaE, rampFactor, PrecipMask
    real(rp), dimension(:,:), allocatable :: beta_RCM, alpha_RCM, gtype_RCM
  end type dkAurora_T

  type dkConductance_T
    logical :: const_sigma
    
    ! arrays on the grid
    real(rp), dimension(:,:), allocatable :: zenith, coszen
    real(rp), dimension(:,:), allocatable :: euvSigmaP, euvSigmaH
    real(rp), dimension(:,:), allocatable :: deltaSigmaP, deltaSigmaH
  end type dkConductance_T

  ! DK variables live natively on the DK shell grid as ShellGridVar_T (theta,phi),
  ! matching the ShellGridVar convention used by gam2shell/imag2shell staging.
  ! This removes the transpose at every coupling hop (shellDK, mapDKToMix).
  ! All %data are cell-centered (SHGR_CC) on the ghost-free DK grid: (1:Nt,1:Np).
  type dkVariables_T
    type(ShellGridVar_T) :: Pot, FAC, SigmaP, SigmaH
    type(ShellGridVar_T) :: Sound_Speed, Density, Avg_Eng, Num_Flux, Aur_Type
    type(ShellGridVar_T) :: Im_Eflux, Im_Gtype, Im_Beta, Im_Eden, Im_Epre
    type(ShellGridVar_T) :: Im_Enflx, DeltaE, Im_Npsp, Z_Eavg, Z_Nflux
    type(ShellGridVar_T) :: Im_Eavg

    real(rp) :: tilt=0.
  end type dkVariables_T
  
  ! used to store an entire instance of MIX (e.g., one per hemisphere)
  type dkParameters_T
    ! aurora type parameters
    integer  :: aurora_model_type
    real(rp) :: alpha, beta, R, RinMHD, beta_inp, gamma
    logical  :: doRamp, doChill, doMR, doAuroralSmooth

    ! conductance type parameters
    integer :: euv_model_type, et_model_type, sigma_model_type
    real(rp) :: F107, pedmin, hallmin, sigma_ratio, ped0
    logical :: const_sigma, doStarlight, apply_cap, doEMA, doET
  end type dkParameters_T

  type dkIO_T
    ! IO parameters
    real(rp) :: time, mjd, tilt
    character(len=strLen) :: RunID = ""   ! run tag for native DK restart/output filenames
    real(rp), dimension(:,:), allocatable :: x,y
    real(rp), dimension(:,:,:,:), allocatable :: vars !(Np,Nt,Nvars,Nhemispheres)
  end type dkIO_T

  type, extends(ShellGrid_T) :: dkGrid_T
  ! grid elements needed for DK calculations
    real(rp), dimension(:,:), allocatable :: x, y    ! SM coordinates: sin(th)*cos(ph), sin(th)*sin(ph)
    real(rp), dimension(:,:), allocatable :: r       ! sin(thc) — colatitude mapping for FEDDER ramp
    real(rp), dimension(:,:), allocatable :: D0      ! Gallagher background density (earthhelper.F90)
    real(rp), dimension(:,:), allocatable :: dt, dp  ! cell spacings for CalcMagE (auroralhelper.F90)
    ! Np, Nt inherited from ShellGrid_T — do not redeclare
  end type dkGrid_T


  type, extends(BaseApp_T) :: dkApp_T
    type(dkGrid_T)          :: grid
    type(dkAurora_T)        :: aurora
    type(dkConductance_T)   :: conductance
    type(dkVariables_T)     :: vars
    type(dkParameters_T)    :: params
    type(dkIO_T)            :: IO
    
    contains

    procedure :: InitModel => dkInitModel
    procedure :: InitIO => dkInitIO
    procedure :: WriteRestart => dkWriteRestart
    procedure :: ReadRestart => dkReadRestart
    procedure :: WriteConsoleOutput => dkWriteConsoleOutput
    procedure :: WriteFileOutput => dkWriteFileOutput
    procedure :: WriteSlimFileOutput => dkWriteSlimFileOutput
    procedure :: AdvanceModel => dkAdvanceModel
    procedure :: Cleanup => dkCleanup

  end type dkApp_T

  interface
    module subroutine dkInitModel(App,Xml)
      class(dkApp_T), intent(inout) :: App
      type(XML_Input_T), intent(inout) :: Xml
    end subroutine dkInitModel

    module subroutine dkInitIO(App, Xml)
      class(dkApp_T), intent(inout) :: App
      type(XML_Input_T), intent(inout) :: Xml
    end subroutine dkInitIO

    module subroutine dkWriteRestart(App, nRes)
      class(dkApp_T), intent(inout) :: App
      integer, intent(in) :: nRes
    end subroutine dkWriteRestart

    module subroutine dkReadRestart(App, resId, nRes)
      class(dkApp_T), intent(inout) :: App
      character(len=*), intent(in) :: resId
      integer, intent(in) :: nRes
    end subroutine dkReadRestart

    module subroutine dkWriteConsoleOutput(App)
      class(dkApp_T), intent(inout) :: App
    end subroutine dkWriteConsoleOutput

    module subroutine dkWriteFileOutput(App, nStep)
      class(dkApp_T), intent(inout) :: App
      integer, intent(in) :: nStep
    end subroutine dkWriteFileOutput

    module subroutine dkWriteSlimFileOutput(App, nStep)
      class(dkApp_T), intent(inout) :: App
      integer, intent(in) :: nStep
    end subroutine dkWriteSlimFileOutput

    module subroutine dkAdvanceModel(App, dt)
      class(dkApp_T), intent(inout) :: App
      real(rp), intent(in) :: dt
    end subroutine dkAdvanceModel

    module subroutine dkCleanup(App)
      class(dkApp_T), intent(inout) :: App
    end subroutine dkCleanup

  end interface

  type precip_t  
      integer :: nVars ! change together with the enumerator above
      !type(ShellGridVar_T),target :: FAC,SIGMAP,SIGMAH,SOUND_SPEED,DENSITY,AVG_ENG,NUM_FLUX,NEUTRAL_WIND,EFIELD,IM_EAVG,IM_EFLUX,IM_IAVG,IM_IFLUX,&
      !  Z_EAVG,Z_NFLUX,CRPOT,TPOT,IM_GTYPE,AUR_TYPE,IM_BETA,IM_EDEN,IM_EPRE,IM_ENFLX,IM_INFLX,DELTAE,IM_NPSP,IM_THCON,IM_PHCON
      !type(ShellGridVar_T), pointer :: data(:)
      type(ShellGridVar_T), allocatable :: data(:)
  end type precip_T

end module dktypes

