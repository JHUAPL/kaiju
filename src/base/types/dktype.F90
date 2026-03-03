module dktypes
  use kdefs
  use shellGrid
  
  use basetypes

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

  type dkState_T
    type(ShellGridVar_T) :: 
  end type dkState_T


  type, extends(BaseApp_T) :: dkApp_T
    ! figure out what dk has
    type(ShellGrid_T)      :: Gr
    
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

