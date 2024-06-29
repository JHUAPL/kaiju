module mixtypes
  use mixdefs

  implicit none

  type mixParams_T
     ! constants to be read from XML input file

     ! conductance model
     integer  :: euv_model_type
     integer  :: et_model_type
     integer  :: aurora_model_type
     real(rp) :: alphaZ
     real(rp) :: betaZ
     real(rp) :: alpha
     real(rp) :: beta
     real(rp) :: gamma
     real(rp) :: r
     real(rp) :: f107
     real(rp) :: pedmin 
     real(rp) :: hallmin
     real(rp) :: sigma_ratio
     real(rp) :: ped0
     logical  :: const_sigma
     logical  :: doRamp
     logical  :: doChill
     logical  :: doStarlight
     logical  :: doMR
     logical  :: doAuroralSmooth
     logical  :: apply_cap
     logical  :: doSWF107
     logical  :: doEMA

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

  type mixConductance_T
    integer :: euv_model_type, et_model_type, aurora_model_type
    real(rp) :: alpha, beta, R, F107,pedmin,hallmin,sigma_ratio,ped0, alphaZ, betaZ
    logical :: const_sigma, doRamp, doChill, doStarlight, apply_cap, doMR, doAuroralSmooth, doEMA

    ! arrays on the grid
    real(rp), dimension(:,:), allocatable :: zenith, coszen
    real(rp), dimension(:,:), allocatable :: euvSigmaP, euvSigmaH
    real(rp), dimension(:,:), allocatable :: deltaSigmaP, deltaSigmaH
    real(rp), dimension(:,:), allocatable :: E0, phi0, deltaE, aRes
    real(rp), dimension(:,:), allocatable :: rampFactor, AuroraMask, PrecipMask, drift
    real(rp), dimension(:,:), allocatable :: engFlux, avgEng
  end type mixConductance_T
  
  ! used to store an entire instance of MIX (e.g., one per hemisphere)
  type mixIon_T
     type(mixState_T)       :: St
     type(mixGrid_T)        :: G       ! G - primary MIX grid used for the solver. 
     type(mixGrid_T)        :: mixGfpd ! mixGfpd - flipped grid for mapping from MHD, moved from mhd2mix type.  
     type(mixGrid_T)        :: Ggeo    ! Ggeo - G grid converted to geo updated every coupling step
     type(mixGrid_T)        :: Gapx    ! Gapx - G grid converted to apex updated every coupling step
     type(mixParams_T)      :: P
     type(Solver_T)         :: S
     type(mixConductance_T) :: conductance
     real(rp)               :: rad_iono_m   ! Ionosphere radius in meters
     real(rp)               :: rad_planet_m ! Planet     radius in meters
  end type mixIon_T

  ! used to store all instances of mixIon type, i.e., all hemispheres
  type mixApp_T
     type(mixIon_T), dimension(:), allocatable :: ion
  end type mixApp_T

  ! use this as a container to store the variables read from a previous H5 file
  type mixIO_T
     real(rp) :: time, mjd, tilt
     real(rp), dimension(:,:), allocatable :: x,y
     real(rp), dimension(:,:,:,:), allocatable :: vars !(Np,Nt,Nvars,Nhemispheres)
  end type mixIO_T

end module mixtypes

