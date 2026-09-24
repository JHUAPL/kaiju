module dktypes
  use kdefs
  use shellGrid
  use basetypes

  implicit none

  ! === DragonKing model / output codes (formerly module dkdefs)
  type :: dkOpts_T
     ! precipitation/aurora_model_type
     integer :: FEDDER = 1, LINMRG = 2, SUNNY = 3
     ! conductance/euv_model_type
     integer :: LOMPE = 1
     ! conductance/et_model_type  (electron temperature)
     integer :: ET_NONE = 1, ET_ADHOC = 2, ET_CMIT = 3
     ! conductance/sigma_model_type  (Robinson variant)
     integer :: KAE23 = 1, ROB87 = 2
     ! auroral output type codes (written to h5 as the Aur_Type variable)
     ! 0: No precipitation; 1: Precipitation without mono acceleration; 
     ! 2: Precipitation with finite mono acceleration but not producing stronger conductance than diffuse; 
     ! 3: Mono-energetic precipitation.
     integer :: AT_NoPre = 0, AT_RMnoE = 1, AT_RMfnE = 2, AT_RMono = 3
  end type dkOpts_T
  type(dkOpts_T), parameter :: DKOpts = dkOpts_T()

  real(rp), parameter :: eTINY = 1.D-8 ! Floor of average energy [keV]

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
  
  ! DragonKing run configuration
  type dkModel_T
    ! aurora type parameters
    integer  :: aurora_model_type
    real(rp) :: alpha, beta, R, RinMHD, beta_inp, gamma
    logical  :: doRamp, doChill, doMR, doAuroralSmooth

    ! conductance type parameters
    integer :: euv_model_type, et_model_type, sigma_model_type
    real(rp) :: F107, pedmin, hallmin, sigma_ratio, ped0
    logical :: const_sigma, doStarlight, apply_cap, doEMA, doET
    logical :: doSWF107     !Take F107 from the SW instead of XML
  end type dkModel_T

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
    type(dkModel_T)         :: model
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

end module dktypes

