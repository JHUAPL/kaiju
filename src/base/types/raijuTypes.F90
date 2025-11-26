module raijutypes

    use helpertypes
    use shellgrid
    use xml_input
    use ioclock
    use kronos

    use basetypes

    use raijudefs

    implicit none

    !  var*Var0 = CODE units --> In/Out units
    type kmUnits_T
        real(rp) :: V0   = 1.0  ! [m/s]
        real(rp) :: N0   = 1.0  ! [kg/m3]
        real(rp) :: T0   = 1.0  ! [s]
        real(rp) :: Eta0 = 1.0  ! [#/(Bs*Rx^2)] --> [#/Wb]
        real(rp) :: Pot0 = 1.0  ! [Volts]
    end type kmUnits_T

!------
! Species
!------
    type raijuSpecies_T
        !! Container for all info related to a specific species
        
        ! These are all specified by raijuconfig.h5
        character(len=strLen) :: name
        integer :: N
        integer :: flav = -1
            !! Species flavor
        integer :: numNuc_p
            !! Number of protons in nucleus
        integer :: numNuc_n
            !! Number of neutrons in nucleus
        integer :: q
            !! Net charge of species
        real(rp) :: fudge
            !! Strong-scattering limit
        real(rp), dimension(:), allocatable :: alami
            !! [eV*(Rx/nT)^(2/3)] Lambda channel cell interfaces/edges

        ! These are calculated after read-in
        integer :: kStart, kEnd
            !! Start and end indices for this species in Nk-size arrays
        real(rp) :: amu
            !! Species mass in amu (even the electrons)
        integer :: spcType
            !! Enum of species type

        ! Determined based on coupling details
        logical :: isMappedTo = .false.
        
    end type raijuSpecies_T


!------
! Losses
!------

    type :: baseRaijuLoss_T
        !! Base loss process type for specific loss implementations to extend
        logical :: isPrecip = .false.

        contains
        
        ! Default do-nothing placeholders
        procedure baseLossInit, baseLossUpdate, baseLossCalcTau, baseLossValidSpc

        procedure :: doInit   => baseLossInit
            !! Initialize any arrays and params that persist throughout run
        procedure :: doUpdate => baseLossUpdate
            !! Update any variables based on current raijuState
        procedure :: isValidSpc => baseLossValidSpc
            !! Tell someone if this loss process applies to given species
        procedure :: calcTau  => baseLossCalcTau
            !! Report instantaneous loss rate for a given energy and lat/lon
        procedure :: doOutput  => baseLossDoOutput
            !! Output any relevant state info to file
    end type baseRaijuLoss_T

    type :: raijuLPHolder_T
        !! Container for collection of loss processes
        class(baseRaijuLoss_T), allocatable :: p
    end type raijuLPHolder_T

!------
! Precipitation models
!------
!    type eLossWM_T
!        !! Parameters used in electron wave model from Dedong Wang and Shanshan Bao
!        
!        ! -- Model -- !
!        real(rp) :: NpsphHigh = def_NpsphHigh
!        real(rp) :: NpsphLow = def_NpsphLow
!        real(rp) :: ChorusLMax = def_ChorusLmax
!        real(rp) :: PsheetLMin = def_PsheetLmin
!        real(rp) :: ChorusEMin = def_ChorusEMin
!
!        logical :: doOutput = .false.
!            !! Whether or not we will be asked to provide detailed info in the output file
!        
!        
!        type(TimeSeries_T) :: KpTS
!            !! Kp data from wind file
!
!        ! -- Grid -- !
!        ! Chorus info
!        integer :: Nkp, Nmlt, Nl, Ne
!            !! Number of bins for Kp, MLT, L shell, and Energy
!        real(rp), allocatable, dimension(:) :: Kp1D
!            !! 1D array of Kp dimension for Tau4D
!        real(rp), allocatable, dimension(:) :: MLT1D
!            !! 1D array of MLT dimension for Tau4D
!        real(rp), allocatable, dimension(:) :: L1D
!            !! 1D array of L shell dimension for Tau4D [Re]
!        real(rp), allocatable, dimension(:) :: Energy1D
!            !! 1D array of energy dimension for Tau4D [MeV]
!        real(rp), allocatable, dimension(:,:,:,:) :: Tau4D
!            !! Tau(Kp, MLT, L, E) table electron scattering table [seconds]
!
!        ! -- State -- !
!        real(rp), allocatable, dimension(:,:) :: wPS
!        real(rp), allocatable, dimension(:,:) :: wHISS
!        real(rp), allocatable, dimension(:,:) :: wCHORUS
!
!    end type eLossWM_T

!------
! General helpers
!------
    type raijuColdstarter_T
        !! Hold various settings and state info to help manage raiju's cold starting

        !--- Model ---!
        real(rp) :: tEnd
            !! [s] Sim time we stop doing cold start stuff
        real(rp) :: evalCadence = 300
            !! [s] How often to check target Dst and adjust raiju's etas
        logical :: doCX = .true.
            !! Whether or not we apply charge exchange
        logical :: doUpdate = .false.
            !! Whether or not we do updates after first init


        !--- State ---!
        logical :: doCS_next_preAdv = .false.
            !! Signal to run coldstart next time raiju preAdvances
        real(rp) :: modelDst_next_preAdv = 0.0_rp
            !! Target Dst [nT] when we run coldstart next

        logical :: doneFirstCS = .false.
            !! Have we executed once already?
        real(rp) :: lastEval = -1*HUGE
            !! [s] Last eval time
        real(rp) :: lastTarget = 0
            !! [nT] Last target Dst


    end type raijuColdstarter_T

!------
! Coupling helpers
!------
    type mhd2raiSpcMap_T
        integer :: idx_mhd
            !! Index of mhd fluid (in Gas array, probably)
        integer :: flav
            !! RAIJU flavor to map to
        logical :: doExcessToPsph = .false.
            !! Whether or not we map low-energy part of distribution to zero-energy plasmasphere channel
    end type mhd2raiSpcMap_T

!------
! Main Model, Grid, State
!------
    type raijuModel_T

        ! Misc. bookkeeping stuff
        character(len=strLen) :: RunID = ""    
        character(len=strLen) :: configFName
            !! Filename of the .h5 config that holds lambda grid, wavemodel values, etc.
        character(len=strLen) :: raijuH5
            !! Filename of the h5 file we output to
        character(len=strLen) :: tsF
            !! Filename for timeseries information like mjd, solar wind (bcwind.h5)

        ! Restart info
        logical :: isRestart
            !! Whether we initialize state from a restart file
        integer :: nResIn = -1
            !! Restart number to use for state init (-1 for symlink)
        character(len=strLen) :: resId = ""
            !! Restart filename id (resID.raiju.Res.xxxxx.h5) to read from
        character(len=strLen) :: ResF = ""
            !! Restart filename to read from
        
        integer :: nSpc, nSpcMHD, nG
            !! Number of species in raiju
            !! Number of species/fluids in MHD
            !! # ghosts
        real(rp) :: t0, tFin, dt
            !! Start and end time, delta coupling time (may be replaced/ignored later by voltron setting a dynamic coupling time)
        
        ! Solver params
        integer :: maxItersPerSec
        real(rp) :: CFL
            !! CFL condition used in deciding time step
        integer :: maxOrder
            !! Max reconstruction order we will use when calculating etas at interfaces
        real(rp) :: PDMB
            !! 0 < PDMB < 1 used in PDM limiter. 1 = very diffusive
        logical :: doUseVelLRs

        ! https://media0.giphy.com/media/XfDPdSRhYFUhIU7EPw/giphy.gif
        logical :: isSA
            !! Is RAIJU running in statnd-alone mode. If so, we shouldn't expect any coupling information to exist
        logical :: fixedTimestep
            !! Fixed or dynamic timestep
        logical :: isMPI
            !! Are we in MPI mode (does not exist yet so better be no)
        logical :: doLosses
            !! Whether or not to calculate eta losses during advance
        logical :: isLoud
            !! For debug
        logical :: writeGhosts
            !! For debug
        logical :: doSmoothGrads
            !! Whether or not we smooth variables (bvol and electric potential) before taking gradients
        logical :: doClockConsoleOut
            !! If we are driving, output clock info
        logical :: doOutput_debug
            !! Dump lots of otherwise unnecessary stuff
        logical :: doOutput_3DLoss
            !! Dump extra 3D loss variables
        logical :: doOwnCorot
            !! If true, we calculate corortation potential ourselves and include it in our velocity calculation
            !! If false, we assume its already included in whatever ionospheric potential we're using
        logical :: doGeoCorot
            !! If true, calc corotation potential from Geopack
            !! If false, calc corotation potential assuming dipole and rotational axes are aligned
        logical :: doExcessToPsph
            !! Allow mapping of excess H+ to plasmasphere channel

        ! Plasmasphere settings
        logical :: doPlasmasphere
            !! Use for now to determine if we should be doing plasmasphere stuff
            !! Likely, in the future, we will determine automatically by the presence of a flavor 0
        real(rp) :: psphInitKp
        logical :: doPsphEvol
            !! Whether or not to actually evolve the plasmasphere
        real(rp) :: psphEvolRad
            !! [Rp] Radius below which plasmasphere is not evolved
        ! TODO: Extra params for refilling rate, determining initial profile, etc.

        ! Some constants
        real(rp) :: tiote  ! Ion temp over electron temp. In the future, should be fancier

        ! Active domain settings
        integer :: n_bndLim
            !! cell dt/dp slope
        real(rp) :: maxTail_buffer
            !! Maximum tailward extent of the buffer region
        real(rp) :: maxSun_buffer
            !! Maximum sunward extent of the buffer region
        real(rp) :: maxTail_active
            !! Maximum tailward extent of the active region
        real(rp) :: maxSun_active
            !! Maximum sunward extent of the active region
        real(rp) :: activeDomRad
            !! [Rp] Cells are forced to be active below this radius

        ! Active shell settings
        logical :: doActiveShell
            !! Use activeShell logic to try to boost dt
        real(rp) :: worthyFrac  
            !! Fracton that a channel must contribute to pressure or density for its i shell to be evolved

        ! Lambda controls
        real(rp) :: kappa
            !! Kappa value, used in case of Kappa mapping from moments to eta channels
        logical :: doDynamicLambdaRanges
            !! Dynamic lambda ranges so we only evolve certain energies at certain L shells

        ! Detailed information
        type(planet_T) :: planet  
            !! Planet info like radius, mag. moment, etc.

        ! Losses
        !logical :: doSS
            !! Do strong scattering
        logical :: doCC
            !! Do coulomb collisions
        logical :: doCX
            !! Do charge exchange
        !logical :: doFLC
            !! Do field-line curvature
        logical :: doEWM
            !! Do electron loss model
        character(len=strLen) :: ewmType
            !! Type of electron loss model to use
            !! C05 = Chen05 blended with Hiss
            !! BW = Bao-Wang empirical wave model with Hiss, Chorus, and Plasma sheet strong scattering
            !! SS = Schulz strong scattering only, we hate electrons

        ! Coupling info
        integer :: nFluidIn = 0
        type(mhd2raiSpcMap_T), dimension(:), allocatable :: fluidInMaps
        ! Coupling-related knobs
        !real(rp) :: vaFracThresh, bminThresh, normAngThresh, PstdThresh
        real(rp) :: nBounce
        ! Domain determination knobs
        real(rp) :: lim_vaFrac_soft, lim_vaFrac_hard
            !! soft/hard limit for flux tube energy partition
        real(rp) :: lim_bmin_soft, lim_bmin_hard
            !! [nT] soft/hard limit for minimum allowable magnetic field strength at bmin point

        character(len=strLen) :: icStr
        procedure(raijuStateIC_T     ), pointer, nopass :: initState => NULL()
        procedure(raijuUpdateV_T     ), pointer, nopass :: updateV   => NULL()
        !> TODO: Retire this and just have a single function with certain options like maxwellian or kappa
        procedure(raijuDP2EtaMap_T   ), pointer, nopass :: dp2etaMap => NULL()

        ! Hax?
        logical :: doHack_rcmEtaPush = .false.  ! If on, we explicitly do eta conservation and ignore difference of magnetic flux between cells

    end type raijuModel_T


    type raijuGrid_T
        integer :: gType  ! Enum of grid type

        type(ShellGrid_T) :: shGrid
        real(rp), dimension(:), allocatable :: thRp
            !! (Ngi+1) [radians] Corner-centered theta value at planet's surface (1 Rp)
        real(rp), dimension(:), allocatable :: thcRp
            !! (Ngi) [radians] Cell-centered theta value at planet's surface (1 Rp)
        real(rp), dimension(:), allocatable :: delTh
            !! (Ngi+1) [radians] Delta theta between cell centers. For cell i, delTh(i) = lower theta del, delTh(i+1) = higher theta del
        real(rp), dimension(:), allocatable :: delPh
            !! (Ngj+1) [radians] Delta phi between cell centers. For cell j, delPh(j) = lower phi del, delPh(j+1) = higher phi del
        real(rp), dimension(:,:), allocatable :: areaCC
            !! (Ngi, Ngj) [Rp^2] Area of each cell
        real(rp), dimension(:,:,:), allocatable :: areaFace
            !! (Ngi+1, Ngj+1, 2) [Rp^2] Estimated area at each interface. (i,j,1) = theta dir, (i,j,2) = phi dir
        real(rp), dimension(:,:,:), allocatable :: lenFace
            !! (Ngi+1, Ngj+1, 2) [Rp] arc length of each face
        real(rp), dimension(:,:), allocatable :: Bmag
            !! (Ngi, Ngj) [nT] Magnitude of B field at ionosphere (cell-centered)
        real(rp), dimension(:,:), allocatable :: cosdip
            !! (Ngi, Ngj) Cosine of the dip angle at ionosphere (cell-centered)
        real(rp), dimension(:,:,:), allocatable :: BrFace
            !! (Ngi+1, Ngj+1, 2) [nT] Radial/normal component of B field [nT] at ionosphere (cell faces)
        real(rp), dimension(:,:), allocatable :: Brcc
            !! (Ngi, Ngj) [nT] same as BrFace but at cell centers

        integer :: nB ! Number of buffer cells between open region and active domain

        ! Flags
        logical :: ignoreConfigMismatch
            !! In the case that the config file has more species than Model%nSpc,
            !! raiju will complain and die if this is false, but will carry on if true

        ! MPI things
        integer :: NumRk  ! Number of ranks in energy space
        integer :: Rk=0   ! Number of this rank

        ! (Nj) I of last active cell for each j slice
        integer, dimension(:), allocatable :: iBnd

        ! Species / lambda stuff
        integer :: Nk  ! Total number of channels for all species
        integer :: nSpc  ! Model has the main copy of this, but helpful to keep here too
        integer :: nFluidIn  ! Model has the main copy of this, but helpful to keep here too
        type(raijuSpecies_T), dimension(:), allocatable :: spc
            !! Collection of raijuSpecies that contain all relevant species info, including alami
        real(rp), dimension(:), allocatable :: alamc
            !! Cell-centered lamba channel values
        integer, dimension(:), allocatable :: k2spc
            !! Nk length mapping of k value to corresponding species index

    end type raijuGrid_T


    type raijuState_T
        logical :: isFirstCpl = .true.

        real(rp) :: t, dt
            !! Current time and last coupling dt made
        real(rp), dimension(:), allocatable :: dtk
            !! Time step for every lambda channel
        integer, dimension(:), allocatable :: nStepk
            !! Number of steps each channel has been evolved
        real(rp) :: mjd
            !! Current mjd
        integer :: ts, tss
            !! Current coupling timestep and sub-stepping timestep
        type(IOClock_T) :: IO
            !! Timers for IO operations

        ! I feel like philosophically this should be in Grid but that feels weird so its here
        type(TimeSeries_T) :: KpTS
            !! Kp timeseries from solar wind file

        ! -- Solver values -- !
        real(rp), dimension(:,:,:), allocatable :: eta
            !! (Ngi, Ngj, Nk) [#/cc * Rp/T] etas
        real(rp), dimension(:,:,:), allocatable :: eta_half
            !! (Ngi, Ngj, Nk) [#/cc * Rp/T] etas  0.5*dt forward from t
        real(rp), dimension(:,:,:), allocatable :: eta_last
            !! (Ngi, Ngj, Nk) [#/cc * Rp/T] etas from previous time step, used for halt-time calculation
        real(rp), dimension(:,:,:), allocatable :: eta_avg
            !! (Ngi, Ngj, Nk) [#/cc * Rp/T] etas averaged over time step State%dt

        logical, dimension(:,:), allocatable :: activeShells
            !! (Ngi, Nk) i shells that should be evolved for a given lambda
        real(rp), dimension(:,:,:), allocatable :: pEff
            !! (Ngi+1, Ngj+1, Nk) [V] Effective potential (ExB + corot + gradient-curvature)
            !! Not actually used in calculations, but helpful for output
        real(rp), dimension(:,:,:), allocatable :: gradPotE, gradPotCorot, gradVM
            !! (Ngi+1, Ngj+1,2) [V/m] Th/phi gradient of the ionospheric potential, corotation potential, and flux tube volume raised to -2/3
            !! units of gradVM are [FTV^(2/3)/m]. Multiplying by lambda gives units of [V/m]
        real(rp), dimension(:,:,:), allocatable :: gradPotE_cc, gradPotCorot_cc, gradVM_cc
            !! (Ngi, Ngj,2) [V/m] Th/phi gradient of the ionospheric potential, corotation potential, and flux tube volume raised to -2/3
            !! Cell-centered, calculated using Green-Gauss method
        real(rp), dimension(:,:,:,:), allocatable :: iVel
            !! (Ngi+1, Ngj+1, Nk, 2) [m/s] Edge-centered normal velocities
        real(rp), dimension(:,:,:,:), allocatable :: iVelL
            !! (Ngi+1, Ngj+1, Nk, 2) [m/s] Left edge-centered normal velocity reconstructed from cell-centered velocities
        real(rp), dimension(:,:,:,:), allocatable :: iVelR
            !! (Ngi+1, Ngj+1, Nk, 2) [m/s] Right edge-centered normal velocity reconstructed from cell-centered velocities
        real(rp), dimension(:,:,:,:), allocatable :: cVel
            !! (Ngi, Ngj, Nk, 2) [m/s] Cell-centered velocities


        ! -- Variables coming from MHD flux tube tracing -- !
        real(rp), dimension(:,:,:), allocatable :: Pavg
            !! (Ngi, Ngj, Ns) [nPa] Average pressure along flux tube
        real(rp), dimension(:,:,:), allocatable :: Davg
            !! (Ngi, Ngj, Ns) [#/cc] Average density along flux tube
        real(rp), dimension(:,:,:), allocatable :: Pstd
            !! (Ngi, Ngj, Ns) Normalized standard deviation of the species pressure along the field line
        real(rp), dimension(:,:,:), allocatable :: Dstd
            !! (Ngi, Ngj, Ns) Normalized standard deviation of the species density along the field line
        real(rp), dimension(:,:), allocatable :: tiote
            !! (Ngi, Ngj) Ratio of ion temperature to electron temperature
        real(rp), dimension(:,:,:), allocatable :: Bmin
            !! (Ngi+1, Ngj+1, NDIM) [nT] Bmin vector
        real(rp), dimension(:,:,:), allocatable :: xyzMin
            !! (Ngi+1, Ngj+1, 3) [Rp] bMin xyz coordinates
        real(rp), dimension(:,:,:), allocatable :: xyzMincc
            !! (Ngi, Ngj, 3) [Rp] cell-centered bMin xyz coordinates
        integer , dimension(:,:), allocatable :: topo   
            !! (Ngi+1, Ngj+1) Topology (0=open, 1=closed)
        real(rp), dimension(:,:), allocatable :: thcon
            !! (Ngi+1, Ngj+1) Co-latitude  of conjugate points
        real(rp), dimension(:,:), allocatable :: phcon
            !! (Ngi+1, Ngj+1) Longitude of conjugate points
        real(rp), dimension(:,:), allocatable :: bvol
            !! (Ngi+1, Ngj+1) [Rp/nT] Flux-tube volume
        real(rp), dimension(:,:), allocatable :: bvol_cc
            !! (Ngi, Ngj) [Rp/nT] Flux-tube volume averaged from corners
        real(rp), dimension(:,:), allocatable :: vaFrac
            !! (Ngi+1, Ngj+1) Fraction of total velocity coming from Alfven speed
            !! Used to limit active region to tubes that can reasonably be treated as averaged and slowly-evolving

        real(rp), dimension(:,:), allocatable :: domWeights
            !! (Ngi,Ngj) weights for incoming moments based on domain constraints

        ! -- Coming from ionosphere solve -- !
        real(rp), dimension(:,:), allocatable :: espot
            !! (Ngi+1, Ngj+1) [kV] electro-static potential
        real(rp), dimension(:,:), allocatable :: pot_corot
            !! (Ngi+1, Ngj+1) [kV] corotation potential
        
        integer, dimension(:), allocatable :: bndLoc
            !! (Ngi) i value of boundary between inactive and buffer domain
        ! (Ngi, Ngj) cell-centered values
        integer , dimension(:,:), allocatable :: active
            !! (Ngi, Ngj) (-1=inactive, 0=buffer, 1=active)
        integer , dimension(:,:), allocatable :: active_last
            !! (Ngi, Ngj) Active domain from the previous coupling time step, used for half-step eta calculation
        integer , dimension(:,:), allocatable :: OCBDist
            !! (Ngi, Ngj) Cell distance from open-closed boundary

        ! Loss-related things
        type(raijuLPHolder_T), dimension(:), allocatable :: lps
            !! Array of loss processes
        integer :: lp_cc_idx = -1
            !! If CC active, this is the index of CC object within lps array

        ! Other helpers
        type(raijuColdstarter_T) :: coldStarter
        
        ! (Ngi, Ngj, Nk) Varibles coming from RAIJU
        real(rp), dimension(:,:,:), allocatable :: lossRates
            !! (Ngi, Ngj, Nk) [1/s] Loss rates for each grid and lambda point. Generally stays the same over coupling time so we store them all here
        real(rp), dimension(:,:,:), allocatable :: lossRatesPrecip
            !! (Ngi, Ngj, Nk) [1/s] Loss rates that result in precipitation. Should be <= lossRates
        real(rp), dimension(:,:,:), allocatable :: precipType_ele
            !! (Ngi, Ngj, Nk) Prepication type used for electrons
        !real(rp), dimension(:,:,:), allocatable :: precipNFlux
        type(ShellGridVar_T), dimension(:), allocatable :: precipNFlux
            !! (Ngi, Ngj, Nk) [#/cm^2/s] Precipitation number fluxes
        !real(rp), dimension(:,:,:), allocatable :: precipEFlux
        type(ShellGridVar_T), dimension(:), allocatable :: precipEFlux
            !! (Ngi, Ngj, Nk) [erg/cm^2/s] Precipitation energy fluxes
        real(rp), dimension(:,:,:), allocatable :: dEta_dt
            !! (Ngi, Ngj, Nk) [eta units/s] Average 
        real(rp), dimension(:,:,:), allocatable :: CCHeatFlux
            !! (Ngi, Ngj, Nk) [erg/cm^2/s] Heat flux from RC ions to plasmasphere electrons due to Coulumb collisions
        
        ! (Ngi, Ngj, Nspc+1) (First Nspc index is bulk) Moments
        ! Last dimension will be D/P of different populations (not necessarily same as species)
        ! Example: Total, hot protons, cold protons, electrons, other species
        !real(rp), dimension(:,:,:), allocatable :: Den
        !real(rp), dimension(:,:,:), allocatable :: vAvg
        !real(rp), dimension(:,:,:), allocatable :: Press
        type(ShellGridVar_T), dimension(:), allocatable :: Den
            !! (Ngi, Ngj, Nspc+1) Density, evaluated using instantaneous etas  [#/cc]
        type(ShellGridVar_T), dimension(:), allocatable :: Press
            !! (Ngi, Ngj, Nspc+1) Pressure, evaluated using instantaneous etas [nPa]
        type(ShellGridVar_T), dimension(:), allocatable :: Den_avg
            !! (Ngi, Ngj, Nspc+1) Density, evaluated using etas averaged over time State%dt  [#/cc]
        type(ShellGridVar_T), dimension(:), allocatable :: Press_avg
            !! (Ngi, Ngj, Nspc+1) Pressure, evaluated using etas averaged over time State%dt [nPa]
        type(ShellGridVar_T) :: Tb
            !! (Ngi, Ngj) [s] Bounce timescale (Alfven crossing time)
        


        !> Only used when debugging
        real(rp), dimension(:,:,:,:), allocatable :: etaFaceReconL
        real(rp), dimension(:,:,:,:), allocatable :: etaFaceReconR   
        real(rp), dimension(:,:,:,:), allocatable :: etaFacePDML
        real(rp), dimension(:,:,:,:), allocatable :: etaFacePDMR
        real(rp), dimension(:,:,:,:), allocatable :: etaFlux   

    end type raijuState_T


!------
! Higher-level types, using above types
!------
    type raiAppOverride_T
        !! Individual values that, if set, are meant to override otherwise default raiju behavior
        logical, private :: isValSet = .false.
        real(rp), private :: val

        contains

        procedure :: isSet => raiAppOverride_isSet
        procedure :: set   => raiAppOverride_setVal
        procedure :: get   => raiAppOverride_getVal
    end type raiAppOverride_T

    type, extends(BaseOptions_T) :: raijuOptions_T
        ! If set, will use the set value instead of default or xml setting
        type(raiAppOverride_T) :: thetaL
            !! [deg] lower theta (poleward) boundary location
        type(raiAppOverride_T) :: thetaU
            !! [deg] upper theta (equatorward) boundary location
    end type raijuOptions_T

    type, extends(BaseApp_T) :: raijuApp_T
        type(raijuModel_T) :: Model
        type(raijuGrid_T ) :: Grid
        type(raijuState_T) :: State
        
        type(raijuOptions_T) :: opt

        contains

        procedure :: InitModel           => raiInitModel
        procedure :: InitIO              => raiInitIO
        procedure :: WriteRestart        => raiWriteRestart
        procedure :: ReadRestart         => raiReadRestart
        procedure :: WriteConsoleOutput  => raiWriteConsoleOutput
        procedure :: WriteFileOutput     => raiWriteFileOutput
        procedure :: WriteSlimFileOutput => raiWriteSlimFileOutput
        procedure :: AdvanceModel        => raiAdvanceModel
        procedure :: Cleanup             => raiCleanup
    end type raijuApp_T

!------
! Interfaces
!------

    !raijuapp function placeholders, bodies are in src/raiju/raijutypessub.F90 to prevent circular dependency
    interface
        module subroutine raiInitModel(App, Xml)
            class(raijuApp_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine raiInitModel

        module subroutine raiInitIO(App, Xml)
            class(raijuApp_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine raiInitIO

        module subroutine raiWriteRestart(App, nRes)
            class(raijuApp_T), intent(inout) :: App
            integer, intent(in) :: nRes
        end subroutine raiWriteRestart

        module subroutine raiReadRestart(App, resId, nRes)
            class(raijuApp_T), intent(inout) :: App
            character(len=*), intent(in) :: resId
            integer, intent(in) :: nRes
        end subroutine raiReadRestart

        module subroutine raiWriteConsoleOutput(App)
            class(raijuApp_T), intent(inout) :: App
        end subroutine raiWriteConsoleOutput

        module subroutine raiWriteFileOutput(App, nStep)
            class(raijuApp_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine raiWriteFileOutput

        module subroutine raiWriteSlimFileOutput(App, nStep)
            class(raijuApp_T), intent(inout) :: App
            integer, intent(in) :: nStep
        end subroutine raiWriteSlimFileOutput

        module subroutine raiAdvanceModel(App, dt)
            class(raijuApp_T), intent(inout) :: App
            real(rp), intent(in) :: dt
        end subroutine raiAdvanceModel

        module subroutine raiCleanup(App)
            class(raijuApp_T), intent(inout) :: App
        end subroutine raiCleanup

    end interface

    abstract interface
        subroutine raijuStateIC_T(Model,Grid,State,inpXML)
            Import :: raijuModel_T, raijuGrid_T, raijuState_T, strLen, XML_Input_T
            type(raijuModel_T) , intent(in) :: Model
            type(raijuGrid_T)  , intent(in) :: Grid
            type(raijuState_T) , intent(inout) :: State
            type(XML_Input_T), intent(in) :: inpXML
        end subroutine raijuStateIC_T

        subroutine raijuUpdateV_T(Model,Grid,State)
            Import :: raijuModel_T, raijuGrid_T, raijuState_T
            type(raijuModel_T) , intent(in) :: Model
            type(raijuGrid_T)  , intent(in) :: Grid
            type(raijuState_T) , intent(inout) :: State
        end subroutine raijuUpdateV_T

        function raijuDP2EtaMap_T(Model,D,kT,vm,amin,amax) result (etaK)
            Import :: rp, raijuModel_T
            type(raijuModel_T), intent(in) :: Model
            real(rp), intent(in) :: D,kT,vm,amin,amax
            real(rp) :: etaK
        end function raijuDP2EtaMap_T

        function raijuELossRate_T(Model,Grid,State,i,j,k) result (lossRate2)
            Import :: rp, raijuModel_T, raijuGrid_T, raijuState_T
            type(raijuModel_T) , intent(inout) :: Model
            type(raijuGrid_T)  , intent(in) :: Grid
            type(raijuState_T) , intent(in) :: State
            integer, intent(in) :: i,j,k
            real(rp), dimension(2) :: lossRate2
        end function raijuELossRate_T

    end interface

    contains

    ! New loss stuff
    subroutine baseLossInit(this, Model, Grid, xmlInp)
        !Import :: baseRaijuLoss_T, raijuModel_T, raijuGrid_T, XML_Input_T
        class(baseRaijuLoss_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(XML_Input_T) , intent(in) :: xmlInp
    end subroutine baseLossInit

    subroutine baseLossUpdate(this, Model, Grid, State)
        !Import :: baseRaijuLoss_T, raijuModel_T, raijuGrid_T, raijuState_T
        class(baseRaijuLoss_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
    end subroutine baseLossUpdate

    function baseLossValidSpc(this, spc) result(isValid)
        !Import :: baseRaijuLoss_T, raijuModel_T, raijuGrid_T, raijuState_T
        class(baseRaijuLoss_T), intent(in) :: this
        type(raijuSpecies_T), intent(in) :: spc
        logical :: isValid
        isValid = .false.
    end function baseLossValidSpc

    function baseLossCalcTau(this, Model, Grid, State, i, j, k) result(tau)
        !Import :: baseRaijuLoss_T, raijuModel_T, raijuGrid_T, raijuState_T, rp
        class(baseRaijuLoss_T), intent(in) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau
        tau = HUGE
    end function baseLossCalcTau

    subroutine baseLossDoOutput(this, Model, Grid, State, gStr, doGhostsO)
        !Import :: baseRaijuLoss_T, raijuModel_T, raijuGrid_T, raijuState_T, rp
        class(baseRaijuLoss_T), intent(in) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: gStr
        logical, intent(in), optional :: doGhostsO
    end subroutine baseLossDoOutput


    ! Options
    function raiAppOverride_isSet(this) result(isSet)
        class(raiAppOverride_T), intent(in) :: this
        logical :: isSet
        isSet = this%isValSet
    end function raiAppOverride_isSet

    function raiAppOverride_getVal(this) result(val)
        class(raiAppOverride_T), intent(in) :: this
        real(rp) :: val
        val = this%val
    end function raiAppOverride_getVal

    subroutine raiAppOverride_setVal(this, val)
        class(raiAppOverride_T), intent(inout) :: this
        real(rp), intent(in) :: val
        this%val = val
        this%isValSet = .true.
    end subroutine raiAppOverride_setVal

end module raijutypes








