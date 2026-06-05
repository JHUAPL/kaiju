module dktypes
  use kdefs
  use mixdefs
  use shellGrid
  use mixtypes
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

  type dkVariables_T
    real (rp), dimension(:,:),allocatable :: Pot, FAC, SigmaP, SigmaH, SigmaPgcm, SigmaHgcm
    real (rp), dimension(:,:),allocatable :: Sound_Speed, Density, Avg_Eng, Num_Flux, Aur_Type
    real (rp), dimension(:,:),allocatable :: Im_Eflux, Im_Gtype, Im_Beta, Im_Eden, Im_Epre
    real (rp), dimension(:,:),allocatable :: Im_Enflx, DeltaE, Im_Npsp, Z_Eavg, Z_Nflux

    integer :: hemisphere=NORTH
    real(rp) :: tilt=0.

    ! type(ShellGridVar_T)   :: pot_shGr
  end type dkVariables_T

  
  ! used to store an entire instance of MIX (e.g., one per hemisphere)
  type dkParameters_T
    ! aurora type parameters
    integer  :: aurora_model_type
    real(rp) :: alpha, beta, R, RinMHD, beta_inp, gamma
    logical  :: doRamp, doChill, doMR, doAuroralSmooth, doGCM

    ! conductance type parameters
    integer :: euv_model_type, et_model_type, sigma_model_type
    real(rp) :: F107, pedmin, hallmin, sigma_ratio, ped0
    logical :: const_sigma, doStarlight, apply_cap, doEMA, doET
  end type dkParameters_T

  type dkIO_T
    ! IO parameters
    real(rp) :: time, mjd, tilt
    real(rp), dimension(:,:), allocatable :: x,y
    real(rp), dimension(:,:,:,:), allocatable :: vars !(Np,Nt,Nvars,Nhemispheres)
  end type dkIO_T


  type, extends(BaseApp_T) :: dkApp_T
    type(mixGrid_T), dimension(:), allocatable        :: grid
    type(dkAurora_T), dimension(:), allocatable       :: aurora
    type(dkConductance_T), dimension(:), allocatable  :: conductance
    type(dkVariables_T), dimension(:), allocatable    :: vars
    type(dkParameters_T)                              :: params
    type(dkIO_T)                                      :: IO
    
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

