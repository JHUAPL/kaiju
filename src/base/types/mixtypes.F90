module mixtypes
  use mixdefs
  ! use dktypes
  use shellGrid

  use basetypes

  implicit none

  type mixParams_T
    ! constants to be read from XML input file

    ! DK reads now the conductance/precipitation section
    real(rp) :: pedmin
        !! Floor applied to GCM Pedersen conductance in mapOneGCMGridToMix
    real(rp) :: hallmin
        !! Floor applied to GCM Hall conductance in mapOneGCMGridToMix
    logical  :: doSWF107
        !! Take F107 from the solar wind file rather than the input deck

    ! solver
    integer :: maxitr
    integer :: mr
    real(rp) :: tol_abs
    real(rp) :: tol_rel
    real(rp) :: llbc_value

    ! grid
    real(rp) :: LowLatBoundary
    integer :: Np, Nt
    real(rp) :: RinMHD = 0.0 !MHD Rinner if defined
    
    ! init
    logical :: init_from_file

    ! IO
    real(rp) :: dtOut
    integer :: nRes

    ! debug
    integer :: mklmsglvl

  end type mixParams_T

  type mixState_T
    real(rp), dimension(:,:,:),allocatable :: Vars
    integer :: hemisphere=NORTH
    real(rp) :: tilt=0.
    logical :: isIMAG = .false.

    type(ShellGridVar_T)   :: pot_shGr
  end type mixState_T

  type mixGrid_T
    integer :: Np, Nt
    logical :: isPeriodic=.true.
    
    real(rp), dimension(:,:), allocatable :: x,y,t,p,r
    real(rp), dimension(:,:), allocatable :: dt,dp
    real(rp), dimension(:,:), allocatable :: ft,fp
    real(rp), dimension(:,:), allocatable :: dtdt,dpdp
    real(rp), dimension(:,:), allocatable :: cosd
    real(rp), dimension(:,:), allocatable :: D0  ! background density, usually plasmasphere, to apply to precipitation
    real(rp), dimension(:,:,:,:), allocatable :: Interpolant
    integer, dimension(:,:), allocatable :: mask
  end type mixGrid_T

  type Map_T
    real(rp), dimension(:,:,:), allocatable :: M  ! map matrix
    integer, dimension(:,:), allocatable :: I1, J1  ! indices of grid
                                                    ! cells on G1 to
                                                    ! which cells on
                                                    ! G2 belong where
                                                    ! G1(G2) are
                                                    ! grids to(from)
                                                    ! which we
                                                    ! interpolate
  end type Map_T

  type Solver_T
    real(rp), dimension(:), allocatable :: RHS,data
    integer,dimension(:), allocatable :: II,JJ,rowI
    real(rp), dimension(:,:), allocatable :: F11,F22,F12 ! (assuming F21=-F12)
    integer :: nnz  ! number of non-zeros
    ! low lat boundary array
    real(rp),dimension(:),allocatable :: LLBC 
    ! solution storage
    real(rp),dimension(:),allocatable :: solution
  end type Solver_T

!   type mixModel_T
!      type(mixParams_T)      :: P
!      type(dkParams_T)       :: aurora
!      type(mixConductance_T) :: conductance
!      real(rp)               :: rad_iono_m   ! Ionosphere radius in meters
!      real(rp)               :: rad_planet_m ! Planet     radius in meters
!   end type mixModel_T

!   type mixGrid_T
!      type(mixGgrid_T)        :: G       ! G - primary MIX grid used for the solver. 
!      type(mixGgrid_T)        :: mixGfpd ! mixGfpd - flipped grid for mapping from MHD, moved from mhd2mix type.  
!      type(ShellGrid_T)      :: shGr    ! ShellGrid representation of grid
!      type(mixGgrid_T)        :: Ggeo    ! Ggeo - G grid converted to geo updated every coupling step
!      type(mixGgrid_T)        :: Gapx    ! Gapx - G grid converted to apex updated every coupling step
!   end type mixGrid_T
  
!   ! used to store an entire instance of MIX (e.g., one per hemisphere)
!   type mixIon_T
!      type(mixState_T)       :: State
!      type(mixGrid_T)        :: Grid
!      type(mixModel_T)       :: Model
!      type(Solver_T)         :: Solver 
!   end type mixIon_T

  type mixIon_T
    type(mixState_T)       :: St
    type(mixGrid_T)        :: G       ! G - primary MIX grid used for the solver. 
    type(mixGrid_T)        :: mixGfpd ! mixGfpd - flipped grid for mapping from MHD, moved from mhd2mix type.  
    type(ShellGrid_T)      :: shGr    ! ShellGrid representation of grid
    type(mixGrid_T)        :: Ggeo    ! Ggeo - G grid converted to geo updated every coupling step
    type(mixGrid_T)        :: Gapx    ! Gapx - G grid converted to apex updated every coupling step
    type(mixParams_T)      :: P
    type(Solver_T)         :: S
    ! type(dkAurora_T)       :: aurora
    real(rp)               :: rad_iono_m   ! Ionosphere radius in meters
    real(rp)               :: rad_planet_m ! Planet     radius in meters
  end type mixIon_T

   ! use this as a container to store the variables read from a previous H5 file
  type mixIO_T
    real(rp) :: time, mjd, tilt
    real(rp), dimension(:,:), allocatable :: x,y
    real(rp), dimension(:,:,:,:), allocatable :: vars !(Np,Nt,Nvars,Nhemispheres)
  end type mixIO_T


  ! used to store all instances of mixIon type, i.e., all hemispheres
  type, extends(BaseApp_T) :: mixApp_T
    type(mixIon_T), dimension(:), allocatable :: ion
    type(mixIO_T) :: IO
    
    contains

    procedure :: InitModel => mixInitModel
    procedure :: InitIO => mixInitIO
    procedure :: WriteRestart => mixWriteRestart
    procedure :: ReadRestart => mixReadRestart
    procedure :: WriteConsoleOutput => mixWriteConsoleOutput
    procedure :: WriteFileOutput => mixWriteFileOutput
    procedure :: WriteSlimFileOutput => mixWriteSlimFileOutput
    procedure :: AdvanceModel => mixAdvanceModel
    procedure :: Cleanup => mixCleanup

  end type mixApp_T

  interface
    module subroutine mixInitModel(App,Xml)
      class(mixApp_T), intent(inout) :: App
      type(XML_Input_T), intent(inout) :: Xml
    end subroutine mixInitModel

    module subroutine mixInitIO(App, Xml)
      class(mixApp_T), intent(inout) :: App
      type(XML_Input_T), intent(inout) :: Xml
    end subroutine mixInitIO

    module subroutine mixWriteRestart(App, nRes)
      class(mixApp_T), intent(inout) :: App
      integer, intent(in) :: nRes
    end subroutine mixWriteRestart

    module subroutine mixReadRestart(App, resId, nRes)
      class(mixApp_T), intent(inout) :: App
      character(len=*), intent(in) :: resId
      integer, intent(in) :: nRes
    end subroutine mixReadRestart

    module subroutine mixWriteConsoleOutput(App)
      class(mixApp_T), intent(inout) :: App
    end subroutine mixWriteConsoleOutput

    module subroutine mixWriteFileOutput(App, nStep)
      class(mixApp_T), intent(inout) :: App
      integer, intent(in) :: nStep
    end subroutine mixWriteFileOutput

    module subroutine mixWriteSlimFileOutput(App, nStep)
      class(mixApp_T), intent(inout) :: App
      integer, intent(in) :: nStep
    end subroutine mixWriteSlimFileOutput

    module subroutine mixAdvanceModel(App, dt)
      class(mixApp_T), intent(inout) :: App
      real(rp), intent(in) :: dt
    end subroutine mixAdvanceModel

    module subroutine mixCleanup(App)
      class(mixApp_T), intent(inout) :: App
    end subroutine mixCleanup

  end interface



end module mixtypes

