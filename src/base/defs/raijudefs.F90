! kaimag definitions/constants

module raijudefs
    use kdefs

    implicit none

    !------
    ! Enumerators
    !------

    ! Index things
    enum, bind(C)
        enumerator :: RAI_TH=1, RAI_PH  ! Theta and phi directions
    endenum

    ! Grid settings
    enum, bind(C)
        enumerator :: RAI_G_UNISPH, RAI_G_WARPSPH, RAI_G_SHGRID
    endenum

    ! Species Flavors
    enum, bind(c)
        enumerator :: F_PSPH=0,F_HOTE,F_HOTP
            !! These flavors have reserved numbers
    endenum

    ! Species type
    enum, bind(C)
        enumerator :: RAIJUNSPC=0,RAIJUELE,RAIJUHPLUS,RAIJUOPLUS
                    ! Null species, electron, h+, o+
    endenum

    ! Topology
    enum, bind(C)
        enumerator :: RAIJUOPEN=0, RAIJUCLOSED
            !! Whether the field line corresponding to grid point is open or closed
    endenum

    ! Active/buffer/inactive cells
    enum, bind(C)
        enumerator :: RAIJUINACTIVE=-1, RAIJUBUFFER, RAIJUACTIVE
            !! Helps determine how RAIJU is going to treat the grid point
    endenum

    ! Electron loss models
    enum, bind(C)
        enumerator :: RaiELOSS_WM=1
    endenum

    ! Electron wave model loss types
    enum, bind(C)
        enumerator :: RaiEWM_HISS=1,RaiEWM_CHORUS,RaiEWM_PSHEET
    endenum

    !------
    ! Defaults
    !------

    ! Names
    character(len=strLen), parameter :: RAI_SG_NAME = "RAIJU"

    ! Units
    real(rp), parameter :: sclEta = 1.0e9  ! [1/nT -> 1/T on DkT2eta conversion]
    real(rp), parameter :: sclIntens = 1.e-4*sqrt(ev2J/(8.0*dalton))/PI ! code eta to intensity [1/(s*sr*keV*cm^2)]

    ! Solver param defaults and limits
    integer, parameter :: def_maxItersPerSec = 1E3
    real(rp), parameter :: def_pdmb = 0.75
    real(rp), parameter :: def_cfl  = 0.3
    real(rp), parameter :: cflMax = 0.3
    logical, parameter :: def_doUseVelLRs = .true.
    logical, parameter :: def_doSmoothGrads = .true.

    ! Domain limits
    ! Buffer not allowed beyond min of maxTail and maxSun
    real(rp), parameter :: def_maxTail_buffer = 15.0  ! [Rp]
    real(rp), parameter :: def_maxSun_buffer  = 10.0  ! [Rp]
    ! Active not allowed beyond min of maxTail and maxSun
    real(rp), parameter :: def_maxTail_active = 10.0  ! [Rp]
    real(rp), parameter :: def_maxSun_active  = 10.0  ! [Rp]
    ! Active is forced below activeDomRad
    real(rp), parameter :: def_activeDomRad = 3.0  ! [Rp]

    ! Settings
    integer, parameter :: raiRecLen = 8
        !! Reconstruction stencil length
    integer, parameter :: nSpacesDef = 4
        !! Number of i spaces between last good value and active i for species
    real(rp), parameter :: fracWorthyDef = 0.001
        !! Fraction that a lambda channel must contribute to total pressure or density in order to be worthy of being evolved
    real(rp), parameter :: pressFracThreshDef = 0.01
        !! If fraction of total pressure below lowest lambda bin is greater than this, we complain
    real(rp), parameter :: def_tiote = 4.0

    ! Coupling defaults
    real(rp), parameter :: def_vaFracThresh = 0.10  
        !! min allowable fracton of Alfven to total tube velocity
    real(rp), parameter :: def_PstdThresh = 0.10
        !! Max allowable standard deviation of pressure along a field line
    real(rp), parameter :: def_normAngle = 25
        !! [deg] Max allowable angle between any two normals surrounding a cell corner
    real(rp), parameter :: def_bminThresh   = 1.0  
        !! [nT] default allowable bmin strencgh
    real(rp), parameter :: def_nBounce   = 3.0
        !! Number of Alfven bounces (Tb) required to be considered a "good" flux tube for coupling


    real(rp), parameter :: def_lim_vaFrac_soft = 0.6_rp
    real(rp), parameter :: def_lim_vaFrac_hard = 0.4_rp
    real(rp), parameter :: def_lim_bmin_soft  = 5.0_rp
        !! [nT]
    real(rp), parameter :: def_lim_bmin_hard  = 2.0_rp
        !! [nT]

    ! Plasmasphere stuff
    real(rp), parameter :: def_psphEvolRad = 2.25
        !! [Rp] radius under which plasmasphere is not evolved

end module raijudefs