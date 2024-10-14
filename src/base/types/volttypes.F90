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

    implicit none

    enum, bind(C)
        enumerator :: IMAGRCM=1,IMAGRCMX,IMAGSST
    endenum

    !Projection types
    enum, bind(C) 
        !L-phi (equatorial), Lat-Lon (northern ionospheric)
        enumerator :: LPPROJ=1,LLPROJ
    endenum

    !Data for inner mag => gamera variables
    enum, bind(C)
        enumerator :: IMDEN=1,IMX1,IMX2,IMTSCL,IMPR
    endenum
    integer, parameter :: NVARIMAG = 5

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

        procedure baseInitImag
        procedure baseAdvance
        procedure baseEval
        procedure baseIO
        procedure baseRestart

        ! functions to be over-written by specific inner magnetosphere implementations
        procedure :: doInit => baseInitImag
        procedure :: doAdvance => baseAdvance
        procedure :: doEval => baseEval
        procedure :: doIO => baseIO
        procedure :: doConIO => baseConIO
        procedure :: doRestart => baseRestart

    end type innerMagBase_T


    integer, parameter :: mix2mhd_varn = 1  ! for now just the potential is sent back

    type, extends(gamApp_T) :: gamCoupler_T

        ! data for coupling remix to gamera
        real(rp), dimension(:,:,:,:,:), allocatable :: mixOutput
        real(rp), dimension(:,:,:), allocatable :: gPsi
        type(Map_T), allocatable, dimension(:,:) :: PsiMaps
        real(rp) :: rm2g
        real(rp) :: Rion

        ! data for coupling imag to gamera
        real(rp), dimension(:,:,:,:), allocatable :: SrcNC !Node-centered source terms

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
        procedure :: FinishUpdateMhdData => SHgamFinishUpdateMhdData

    end type SHgamCoupler_T

    type, extends(BaseOptions_T) :: VoltOptions_T
        procedure(StateIC_T), pointer, nopass :: gamUserInitFunc

        contains
    end type VoltOptions_T

    type :: voltApp_T

        !Planet information
        type(planet_T) :: planet

        !Voltron state information
        type(TimeSeries_T) :: tilt,symh
        real(rp) :: time, MJD,tFin
        real(rp) :: BSDst=0.0 !Most recent bsdst calculated
        integer :: ts

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

        class(innerMagBase_T), allocatable :: imagApp

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

        !Local gamera object to couple to
        class(gamCoupler_T), allocatable :: gApp
        logical :: doSerialMHD = .true.

        !voltron specific options
        type(VoltOptions_T) :: vOptions

    end type voltApp_T

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

        module subroutine SHgamFinishUpdateMhdData(App, voltApp)
            class(SHgamCoupler_T), intent(inout) :: App
            class(voltApp_T), intent(inout) :: voltApp
        end subroutine

    end interface

    contains

    ! null default subroutines for inner mag base type
    subroutine baseInitImag(imag,iXML,isRestart,vApp)
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
        real(rp), intent(out) :: imW(NVARIMAG)
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

