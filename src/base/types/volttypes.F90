!Voltron types

module volttypes
    use kdefs
    use cmidefs
    use kronos
    use ioclock
    use mixtypes
    use ebtypes
    use gcmtypes
    use helpertypes
    use basetypes
    use gamtypes
    use raijutypes
    use shellGrid
    use voltCplTypes

    implicit none

    enum, bind(C)
        enumerator :: IMAGRCM=1,IMAGRCMX,IMAGSST,IMAGRAIJU
    endenum

    !Projection types
    enum, bind(C) 
        !L-phi (equatorial), Lat-Lon (northern ionospheric)
        enumerator :: LPPROJ=1,LLPROJ
    endenum

    ! Data for voltron => gamera Gas0 variables
    enum, bind(C)
        enumerator :: IONEX=1,IONEY,IONEZ, & !Convection E field from remix potential
                      PROJLAT,PROJLON !Projected (to NH) lat/lon
    endenum
    integer, parameter :: NVARVOLT0 = 5

    ! Data for inner mag => gamera Gas0 variables
    enum, bind(C)
        enumerator :: IM_D_RING=NVARVOLT0+1,IM_P_RING,IM_D_COLD, IM_P_COLD, IM_TSCL
    endenum
    integer, parameter :: NVARIMAG0 = 5
    integer, parameter :: RCFLUID=1,COLDFLUID=2

    !Data for outflow => gamera Gas0 variables
    enum, bind(C)
        enumerator :: WIND_D,WIND_V,WIND_P
    endenum
    integer, parameter :: NVARWIND0 = 3

    integer, parameter :: NVARVOLTSRC = NVARVOLT0 + NVARIMAG0

    ! Data channels for inExyz electric field array
    enum, bind(C)
        enumerator :: IONO_EX=1,IONO_EY,IONO_EZ
    endenum
    integer, parameter :: NVARIONO0 = 3

    enum, bind(C)
        enumerator :: V_GRID_UNIFORM=1, V_GRID_SHAFEE
    endenum

        ! IMAG2MIX field indices
    integer, parameter :: nVars_imag2mix = 10 ! change together with the enumerator below
    enum, bind(C)
        enumerator :: RAI_GTYPE=1,RAI_THCON,RAI_PHCON,RAI_EAVG,RAI_ENFLX,RAI_EFLUX,RAI_EDEN,RAI_EPRE,RAI_NPSP,RAI_CCHF
    endenum

    ! data for imag => remix for conductance
    type imag2Mix_T
        !Assuming IMag data is coming on northern hemisphere
        logical :: isFresh = .false.
        logical :: isInit = .false.
        !Working on logically Cartesian grid w/ lat, lon 
        real(rp), dimension(:), allocatable :: gcolat,glong
        real(rp), dimension(:,:), allocatable :: eflux,eavg,enflx
        real(rp), dimension(:,:), allocatable :: iflux,iavg,inflx
        !latc/lonc are the mappings to the southern hemisphere
        real(rp), dimension(:,:), allocatable :: latc,lonc
        real(rp), dimension(:,:), allocatable :: fac
        real(rp), dimension(:,:), allocatable :: eden,epre,npsp ! add electron density, pressure, plasmasphere density channels to REMIX.
        integer , dimension(:,:), allocatable :: gtype ! RCM grid info: active, buffer, or outside
        
        logical, dimension(:,:), allocatable :: inIMag
    end type imag2Mix_T

    ! data for gamera -> remix conversion
    type mhd2Mix_T
        real(rp), dimension(:,:,:,:,:), allocatable :: mixInput
        real(rp), dimension(:,:,:,:), allocatable :: gJ,gBAvg
        type(Map_T), allocatable, dimension(:,:) :: Jmaps
        real(rp) :: dtAvg,wAvg
    end type mhd2mix_T

    ! data for chimp -> gamera conversion
    type chmp2Mhd_T
        !Ni,Nj,Nk,2 array
        !Holds mapping from cell-centered xyz => x1,x1 (projection coordinates)
        !Projection coordinates can be R,phi (cylindrical) or lat,lon (ionospheric)
        real(rp), dimension(:,:,:,:), allocatable :: xyzSquish
        logical , dimension(:,:,:)  , allocatable :: isGood   !Good projection or not
        logical , dimension(:,:,:)  , allocatable :: isEdible !Good values for MHD to eat (ingest)
        
        integer :: iMax !Possibly changing i-boundary of squish mapping
        real(rp) :: epsSquish,epsds0 !Epsilon parameter for tracing squish/default respectively

    end type chmp2Mhd_T

    ! data for gamera -> chimp conversion
    type mhd2Chmp_T
        logical :: isLonely = .true.
        real(rp) :: Rin
        real(rp) :: lowlatBC

    end type mhd2Chmp_T

    type innerMagBase_T
        contains

        procedure baseInit_Imag
        procedure baseAdvance
        procedure baseEval
        procedure baseIO
        procedure baseRestart

        ! functions to be over-written by specific inner magnetosphere implementations
        procedure :: doInit => baseInit_Imag
        procedure :: doAdvance => baseAdvance
        procedure :: doEval => baseEval
        procedure :: doIO => baseIO
        procedure :: doConIO => baseConIO
        procedure :: doRestart => baseRestart

    end type innerMagBase_T

    type, extends(BaseOptions_T) :: imagOptions_T
        character(len=strLen) :: swF
            !! Solar wind filename
        real(rp) :: mhdRin
            !! GAMERA near-Earth active cell boundary
        real(rp) :: mhdRinG
            !! GAMERA near-Earth ghost cell boundary
        real(rp) :: mjd0
            !! MJD at sim t=0
        logical :: doColdStart
            !! Whether IMAG model should cold start itself at volt%t = 0
        type(ShellGrid_T) :: voltGrid

        contains
    end type imagOptions_T

    type, abstract, extends(BaseApp_T) :: imagCoupler_T

        class(imagOptions_T), allocatable :: opt

        contains

        ! voltApp to IMAG
        procedure(toIMAG_T)  , deferred :: toIMAG
        ! Routines to help get stuff out of IMAG
        procedure(getMomentsIMAG_T), deferred :: getMoments
        procedure(getMomentsPrecipIMAG_T), deferred :: getMomentsPrecip
        
        !procedure :: InitModel           => baseInitImag
        !procedure :: InitIO              => baseInitIOImag
        !procedure :: WriteRestart        => baseWriteRestartImag
        !procedure :: ReadRestart         => baseReadRestartImag
        !procedure :: WriteConsoleOutput  => baseWriteConsoleOutputImag
        !procedure :: WriteFileOutput     => baseWriteFileOutputImag
        !procedure :: WriteSlimFileOutput => baseWriteFileOutputImag
        !procedure :: AdvanceModel        => baseAdvanceModelImag
        !procedure :: Cleanup             => baseCleanupImag

    end type imagCoupler_T

    type, extends(imagCoupler_T) :: raijuCoupler_T

        class(raijuApp_T), allocatable :: raiApp

        ! --- Options --- !
        real(rp) :: startup_blendTscl = 3600.0
            !! [s] Time scale over which we ramp up to full IM_TSCL for MHD ingestion
        logical :: doColdstartCX = .true.
            !! Whether or not we apply charge exchange to initial coldstart ion profile
        real(rp) :: tsclSm_dL = 1.0_rp
            !! [Rp] theta/L-direction smoothing scale for mhd ingestion timescale
        real(rp) :: tsclSm_dMLT = 1.0_rp
            !! [hr] phi/MLT-direction smoothing scale for mhd ingestion timescale

        ! --- Grid --- !
        type(ShellGrid_T) :: shGr
            !! Copy of raijuModel's shellGrid
        
        ! --- State --- !
        real(rp) :: tLastUpdate
            !! Time of last update, according to voltron

        ! Stuff going into raiju
        type(magLine_T), dimension(:,:), allocatable :: magLines
        type(IMAGTube_T), dimension(:,:), allocatable :: ijTubes

        type(ShellGridVar_T), dimension(:), allocatable :: Pavg
            !! (Ngi, Ngj, Ns) [nPa] Average pressure along flux tube
        type(ShellGridVar_T), dimension(:), allocatable :: Davg
            !! (Ngi, Ngj, Ns) [#/cc] Average density along flux tube
        type(ShellGridVar_T), dimension(:), allocatable :: Pstd
            !! (Ngi, Ngj, Ns) Normalized standard deviation of the species pressure along the field line
        type(ShellGridVar_T), dimension(:), allocatable :: Dstd
            !! (Ngi, Ngj, Ns) Normalized standard deviation of the species density along the field line
        type(ShellGridVar_T) :: tiote
            !! (Ngi, Ngj) Ratio of ion temperature to electron temperature
        type(ShellGridVar_T), dimension(NDIM) :: Bmin
            !! (Ngi+1, Ngj+1, NDIM) [nT] Bmin vector
        type(ShellGridVar_T), dimension(NDIM) :: xyzMin
            !! (Ngi+1, Ngj+1, 3) [Rp] bMin xyz coordinates
        type(ShellGridVar_T), dimension(NDIM) :: xyzMincc
            !! (Ngi, Ngj, 3) [Rp] cell-centered bMin xyz coordinates
        type(ShellGridVar_T) :: topo
            !! (Ngi+1, Ngj+1) Topology (0=open, 1=closed)
        type(ShellGridVar_T) :: thcon
            !! (Ngi+1, Ngj+1) Co-latitude  of conjugate points
        type(ShellGridVar_T) :: phcon
            !! (Ngi+1, Ngj+1) Longitude of conjugate points
        type(ShellGridVar_T) :: bvol
            !! (Ngi+1, Ngj+1) [Rp/nT] Flux-tube volume
        type(ShellGridVar_T) :: bvol_cc
            !! (Ngi, Ngj) [Rp/nT] Flux-tube volume averaged from corners
        type(ShellGridVar_T) :: vaFrac
            !! (Ngi+1, Ngj+1) Fraction of total velocity coming from Alfven speed
            !! Used to limit active region to tubes that can reasonably be treated as averaged and slowly-evolving
        type(ShellGridVar_T) :: Tb
            !! (Ngi, Ngj) Bounce timesale
        type(ShellGridVar_T) :: avgBeta
            !! (Ngi, Ngj) Cross-sectional flux-tube area - averaged beta (woof)

        type(ShellGridVar_T) :: pot_total
            !! Total electrostatic potential from (ionosphere + corot) [kV]
        type(ShellGridVar_T) :: pot_corot
            !! Just corotation potential, just for output purposes [kV]

        ! Post-advance stuff leaving raiju
        type(ShellGridVar_T) :: tscl_mhdIngest
            !! [s] Timescale at which we think MHD should ingest our info

        contains

        procedure :: toIMAG => volt2RAIJU
        procedure :: getMoments => getMomentsRAIJU
        procedure :: getMomentsPrecip => getMomentsPrecipRAIJU
        
        procedure :: InitModel           => raiCplInitModel
        procedure :: InitIO              => raiCplInitIO
        procedure :: WriteRestart        => raiCplWriteRestart
        procedure :: ReadRestart         => raiCplReadRestart
        procedure :: WriteConsoleOutput  => raiCplWriteConsoleOutput
        procedure :: WriteFileOutput     => raiCplWriteFileOutput
        procedure :: WriteSlimFileOutput => raiCplWriteSlimFileOutput
        procedure :: AdvanceModel        => raiCplAdvanceModel
        procedure :: Cleanup             => raiCplCleanup

    end type raijuCoupler_T


    integer, parameter :: mix2mhd_varn = 1  ! for now just the potential is sent back

    type, extends(gamApp_T) :: gamCoupler_T

        ! data for coupling remix to gamera
        real(rp), dimension(:,:,:,:,:), allocatable :: mixOutput
        real(rp), dimension(:,:,:), allocatable :: gPsi
        type(Map_T), allocatable, dimension(:,:) :: PsiMaps
        real(rp) :: rm2g
        real(rp) :: Rion

        ! data for coupler
        character(len=strLen) :: vh5File

        contains

        ! only over-riding specific functions
        !procedure :: InitModel => gamCplInitModel
        procedure :: InitIO => gamCplInitIO
        procedure :: WriteRestart => gamCplWriteRestart
        procedure :: ReadRestart => gamCplReadRestart
        procedure :: WriteConsoleOutput => gamCplWriteConsoleOutput
        procedure :: WriteFileOutput => gamCplWriteFileOutput
        procedure :: WriteSlimFileOutput => gamCplWriteSlimFileOutput
        !procedure :: AdvanceModel => gamCplAdvanceModel

        ! add new coupling function which can be over-ridden by children
        procedure :: InitMhdCoupler => gamInitMhdCoupler
        procedure :: StartUpdateMhdData => gamStartUpdateMhdData
        procedure :: PartialUpdateMhdData => gamPartialUpdateMhdData
        procedure :: FinishUpdateMhdData => gamFinishUpdateMhdData

    end type gamCoupler_T

    type, extends(gamCoupler_T) :: SHgamCoupler_T

        contains

        ! only over-riding specific functions
        procedure :: InitModel => SHgamCplInitModel
        procedure :: InitIO => SHgamCplInitIO
        procedure :: WriteRestart => SHgamCplWriteRestart
        procedure :: ReadRestart => SHgamCplReadRestart
        procedure :: WriteConsoleOutput => SHgamCplWriteConsoleOutput
        procedure :: WriteFileOutput => SHgamCplWriteFileOutput
        procedure :: WriteSlimFileOutput => SHgamCplWriteSlimFileOutput
        procedure :: AdvanceModel => SHgamCplAdvanceModel

        ! add new coupling function which can be over-ridden by children
        procedure :: InitMhdCoupler => SHgamInitMhdCoupler
        procedure :: StartUpdateMhdData => SHgamStartUpdateMhdData
        procedure :: PartialUpdateMhdData => SHgamPartialUpdateMhdData
        procedure :: FinishUpdateMhdData => SHgamFinishUpdateMhdData

    end type SHgamCoupler_T

    type, extends(BaseOptions_T) :: VoltOptions_T
        procedure(StateIC_T), pointer, nopass :: gamUserInitFunc

        contains
    end type VoltOptions_T

    type :: voltState_t

        type(Tube_T), dimension(:,:), allocatable :: ijTubes
        type(TubeShell_T)    :: tubeShell
        type(ShellGridVar_T) :: potential_total
            !! (Corners) [kV] Ionospheric ExB and corotation potential
        type(ShellGridVar_T) :: potential_corot
            !! (Corners) [kV] Corotation potential
        type(ShellGridVar_T) :: pot_drop
            !! (Corners) [kV] potential drop along field line

        type(ShellGridVar_T) :: bIonoMag
            !! [nT] Magnitude of the magnetic field at grid's ionosphere location
        type(ShellGridVar_T) :: bIonoRad
            !! [nT] Radial component of the magnetic field at grid's ionosphere location

    end type voltState_t
    type :: voltApp_T

        !Planet information
        type(planet_T) :: planet

        ! Voltron ShellGrid information
        type(ShellGrid_T) :: shGrid
        integer :: gridType
            !! Populated with enum (e.g. V_GRID_UNIFORM)

        ! Voltron state information
        type(TimeSeries_T) :: tilt,symh
        real(rp) :: time, MJD,tFin
        real(rp) :: BSDst=0.0 !Most recent bsdst calculated
        integer :: ts
        type(voltState_t) :: State

        !Voltron output/restart info
        type (IOClock_T) :: IO
        logical :: isLoud = .true. !Console output
        logical :: writeFiles = .true. !File output

        !Apps
        type(mixApp_T) :: remixApp
        type(mhd2Mix_T) :: mhd2mix

        type(ebTrcApp_T)  :: ebTrcApp
        type(mhd2Chmp_T)  :: mhd2chmp
        type(chmp2Mhd_T)  :: chmp2mhd
        type(imag2Mix_T)  :: imag2mix

        type(gcm_T) :: gcm

        !class(innerMagBase_T), allocatable :: imagApp
        class(imagCoupler_T), allocatable :: imagApp

        !Deep coupling information
        real(rp) :: DeepT ! Time of next deep coupling
        real(rp) :: DeepDT ! Time between deep couplings
        !real(rp) :: TargetDeepDT ! Desired deep step from Voltron
        logical  :: doDeep = .false. !Whether to do deep coupling
        real(rp) :: rTrc  !Radius to do tracing (ebSquish) inside of
        integer  :: nTrc = MaxFL !Max number of tracer steps to take (ebSquish)
        integer  :: iDeep  = 0 !Index of max i shell containing deep coupling radius
        integer  :: imType = 0 !Type of inner magnetosphere model (0 = None)
        integer  :: prType = 0 !Type of projection for coupling   (0 = None)
        logical  :: doQkSquish = .false. !Whether or not to do fast squishing
        integer  :: qkSquishStride = 2 ! Stride to use when fast squishing
        logical  :: doGCM = .false.

        !Dynamic coupling info
        logical :: doDynCplDT = .false. !Whether to do dynamic coupling cadence

        !Have special flag to indicate this is Earth, which is special
        logical :: isEarth = .false.

        ! Model flags
        logical :: doGeoCorot = .false.
            !! Whether or not to use geopack/kai2geo corotation potential, or simple axis-aligned potential

        !Local gamera object to couple to
        class(gamCoupler_T), allocatable :: gApp
        logical :: doSerialMHD = .true.

        !voltron specific options
        type(VoltOptions_T) :: vOptions

    end type voltApp_T


    abstract interface
        subroutine toIMAG_T(App, vApp)
            import imagCoupler_T
            import voltApp_T
            class(imagCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: vApp
        end subroutine toIMAG_T

        subroutine getMomentsIMAG_T(App,th,ph,t,imW,isEdible)
            import imagCoupler_T
            import rp, IM_D_RING,IM_TSCL
            class(imagCoupler_T), intent(inout) :: App
            real(rp), intent(in) :: th,ph,t
            real(rp), intent(out) :: imW(IM_D_RING:IM_TSCL)
            logical, intent(out) :: isEdible
        end subroutine getMomentsIMAG_T


        subroutine getMomentsPrecipIMAG_T(App,th,ph,imP,isEdible)
            import imagCoupler_T
            import rp, nVars_imag2mix
            class(imagCoupler_T), intent(inout) :: App
            real(rp), intent(in) :: th,ph
            real(rp), intent(out) :: imP(nVars_imag2mix)
            logical, intent(out) :: isEdible
        end subroutine getMomentsPrecipIMAG_T
    end interface


    ! interface for coupling functions
    interface
        module subroutine gamCplInitIO(App, Xml)
            class(gamCoupler_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine

        module subroutine gamCplWriteRestart(App, nRes)
            class(gamCoupler_T), intent(inout) :: App
            integer, intent(in) :: nRes
        end subroutine

        module subroutine gamCplReadRestart(App, resId, nRes)
            class(gamCoupler_T), intent(inout) :: App
            character(len=*), intent(in) :: resId
            integer, intent(in) :: nRes
        end subroutine

        module subroutine gamCplWriteConsoleOutput(App)
            class(gamCoupler_T), intent(inout) :: App
        end subroutine

        module subroutine gamCplWriteFileOutput(App, nStep)
            class(gamCoupler_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine

        module subroutine gamCplWriteSlimFileOutput(App, nStep)
            class(gamCoupler_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine

        module subroutine gamInitMhdCoupler(App, voltApp)
            class(gamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
        end subroutine

        module subroutine gamStartUpdateMhdData(App, voltApp)
            class(gamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
        end subroutine

        module subroutine gamPartialUpdateMhdData(App, voltApp, vDT)
            class(gamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
            real(rp), intent(in) :: vDT
        end subroutine

        module subroutine gamFinishUpdateMhdData(App, voltApp)
            class(gamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
        end subroutine

    end interface

    ! functions for squish helper specific gamera coupler
    interface
        module subroutine SHgamCplInitModel(App, Xml)
            class(SHgamCoupler_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine

        module subroutine SHgamCplInitIO(App, Xml)
            class(SHgamCoupler_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine

        module subroutine SHgamCplWriteRestart(App, nRes)
            class(SHgamCoupler_T), intent(inout) :: App
            integer, intent(in) :: nRes
        end subroutine

        module subroutine SHgamCplReadRestart(App, resId, nRes)
            class(SHgamCoupler_T), intent(inout) :: App
            character(len=*), intent(in) :: resId
            integer, intent(in) :: nRes
        end subroutine

        module subroutine SHgamCplWriteConsoleOutput(App)
            class(SHgamCoupler_T), intent(inout) :: App
        end subroutine

        module subroutine SHgamCplWriteFileOutput(App, nStep)
            class(SHgamCoupler_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine

        module subroutine SHgamCplWriteSlimFileOutput(App, nStep)
            class(SHgamCoupler_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine

        module subroutine SHgamCplAdvanceModel(App, dt)
            class(SHgamCoupler_T), intent(inout) :: App
            real(rp), intent(in) :: dt
        end subroutine

        module subroutine SHgamCplCleanup(App)
            class(SHgamCoupler_T), intent(inout) :: App
        end subroutine

        module subroutine SHgamInitMhdCoupler(App, voltApp)
            class(SHgamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
        end subroutine

        module subroutine SHgamStartUpdateMhdData(App, voltApp)
            class(SHgamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
        end subroutine

        module subroutine SHgamPartialUpdateMhdData(App, voltApp, vDT)
            class(SHgamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
            real(rp), intent(in) :: vDT
        end subroutine

        module subroutine SHgamFinishUpdateMhdData(App, voltApp)
            class(SHgamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
        end subroutine

    end interface

    ! raijuCoupler_T
    interface
        module subroutine volt2RAIJU(App, vApp)
            class(raijuCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: vApp
        end subroutine

        module subroutine getMomentsRAIJU(App,th,ph,t,imW,isEdible)
            class(raijuCoupler_T), intent(inout) :: App
            real(rp), intent(in) :: th,ph,t
            real(rp), intent(out) :: imW(IM_D_RING:IM_TSCL)
            logical, intent(out) :: isEdible
        end subroutine getMomentsRAIJU

        module subroutine getMomentsPrecipRAIJU(App,th,ph,imP,isEdible)
            class(raijuCoupler_T), intent(inout) :: App
            real(rp), intent(in) :: th,ph
            real(rp), intent(out) :: imP(nVars_imag2mix)
            logical, intent(out) :: isEdible
        end subroutine getMomentsPrecipRAIJU

        module subroutine raiCplInitModel(App, xml)
            class(raijuCoupler_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine 

        module subroutine raiCplInitIO(App, Xml)
            class(raijuCoupler_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine

        module subroutine raiCplWriteRestart(App, nRes)
            class(raijuCoupler_T), intent(inout) :: App
            integer, intent(in) :: nRes
        end subroutine

        module subroutine raiCplReadRestart(App, resId, nRes)
            class(raijuCoupler_T), intent(inout) :: App
            character(len=*), intent(in) :: resId
            integer, intent(in) :: nRes
        end subroutine

        module subroutine raiCplWriteConsoleOutput(App)
            class(raijuCoupler_T), intent(inout) :: App
        end subroutine

        module subroutine raiCplWriteFileOutput(App, nStep)
            class(raijuCoupler_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine

        module subroutine raiCplWriteSlimFileOutput(App, nStep)
            class(raijuCoupler_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine

        module subroutine raiCplAdvanceModel(App, dt)
            class(raijuCoupler_T), intent(inout) :: App
            real(rp), intent(in) :: dt
        end subroutine

        module subroutine raiCplCleanup(App)
            class(raijuCoupler_T), intent(inout) :: App
        end subroutine
        

    end interface

    contains

!    subroutine baseEvalIMAG(App,x1,x2,t,imW,isEdible)
!        class(imagCoupler_T), intent(inout) :: App
!        real(rp), intent(in) :: x1,x2,t
!        real(rp), intent(out) :: imW(NVARIMAG)
!        logical, intent(out) :: isEdible
!
!        imW = 0.0_rp
!        isEdible = .false.
!    end subroutine baseEvalIMAG
!
!    ! null default subroutines for inner mag base type
!    subroutine baseInitImag(App,xml)
!        class(imagCoupler_T), intent(inout) :: App
!        type(XML_Input_T), intent(inout) :: xml
!    end subroutine
!
!    subroutine baseInitIOImag(App,xml)
!        class(imagCoupler_T), intent(inout) :: App
!        type(XML_Input_T), intent(inout) :: xml
!    end subroutine
!
!    subroutine baseWriteRestartImag(App,nRes)
!        class(imagCoupler_T), intent(inout) :: App
!        integer, intent(in) :: nRes
!    end subroutine
!
!    subroutine baseReadRestartImag(App,resId,nRes)
!        class(imagCoupler_T), intent(inout) :: App
!        character(len=*), intent(in) :: resId
!        integer, intent(in) :: nRes
!    end subroutine
!
!    subroutine baseWriteConsoleOutputImag(App)
!        class(imagCoupler_T), intent(inout) :: App
!    end subroutine
!
!    subroutine baseWriteFileOutputImag(App, nStep)
!        class(imagCoupler_T), intent(inout) :: App
!        integer, intent(in) :: nStep
!    end subroutine
!
!    subroutine baseAdvanceModelImag(App,dt)
!        class(imagCoupler_T), intent(inout) :: App
!        real(rp), intent(in) :: dt
!    end subroutine
!
!    subroutine baseCleanupImag(App)
!        class(imagCoupler_T), intent(inout) :: App
!    end subroutine


    ! Old innermagbase stuff
    subroutine baseInit_Imag(imag,iXML,isRestart,vApp)
        class(innerMagBase_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        type(voltApp_T), intent(inout) :: vApp
    end subroutine

    subroutine baseAdvance(imag,vApp,tAdv)
        class(innerMagBase_T), intent(inout) :: imag
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv
    end subroutine

    subroutine baseEval(imag,x1,x2,t,imW,isEdible)
        class(innerMagBase_T), intent(inout) :: imag
        real(rp), intent(in) :: x1,x2,t
        real(rp), intent(out) :: imW(NVARIMAG0)
        logical, intent(out) :: isEdible

        imW = 0.0_rp
        isEdible = .false.
    end subroutine

    subroutine baseIO(imag,nOut,MJD,time)
        class(innerMagBase_T), intent(inout) :: imag
        integer, intent(in) :: nOut
        real(rp), intent(in) :: MJD,time
    end subroutine

    subroutine baseConIO(imag,MJD,time)
        class(innerMagBase_T), intent(inout) :: imag
        real(rp), intent(in) :: MJD,time
    end subroutine

    subroutine baseRestart(imag,nRes,MJD,time)
        class(innerMagBase_T), intent(inout) :: imag
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time
    end subroutine
    
end module volttypes
