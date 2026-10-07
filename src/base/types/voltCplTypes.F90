module voltCplTypes
    !! Lower-level types useful for voltron coupling
    !! Putting here to remove the basic types from the voltron-dependent routines that use then
    !! So other modules can rely on the types without caring about how they're used by others

    use kdefs
    use shellGrid
    use ebtypes
    use basetypes
    use raijudefs, ONLY : def_tiote

    implicit none

    integer, parameter :: MAXTUBEFLUIDS = 5
    integer, parameter :: TUBE_CLOSED=2,TUBE_OPEN=1,TUBE_UNDEF=-1,TUBE_DISASTER=-666

!------
! Flux tubes
!------

    !New tube type for global Voltron grid
    type Tube_T
        real(rp) :: xyz0(NDIM)
        !! XYZ of seed point [Rx]
        real(rp) :: lat0,lon0,invlat
        !! Lat/lon of tube footpoint and invariant latitude [rad]
        real(rp) :: latc,lonc
        !! Lat/lon of conjugate point [rad] if topo=TUBE_CLOSED
        real(rp) :: J11i
        !! Field-aligned current calculated at inner MHD boundary and mapped (preserving J/B) to ionosphere, [micro A/m^2]

        integer :: topo
        !! topo = TUBE_CLOSED,OPEN,UNDEF
        real(rp) :: bmin, X_bmin(NDIM)
        !! Minimum B along field line [nT] and coordinates [Rx]
        real(rp) :: X_maxRad(NDIM)
        !! xyz coords of max point along field line
        real(rp) :: Jxyz_iono(NDIM)
        !! NDIM current just above ionosphere boundary
        !! Note: magnitude is from sampled point, same location as Bxyz_iono
        real(rp) :: Bxyz_iono(NDIM)
        !! NDIM magnetic field just above ionosphere boundary
        real(rp) :: bVol
        !! Flux-tube volume [Rx/nT]
        real(rp) :: Lb
        !! Length of field line [Rx]
        real(rp) :: Tb
        !! Alfven bounce time [s] if CLOSED
        real(rp) :: pot,crpot,potc
        !! Electrostatic potential (TOTAL and corotation only), and conjugate potential [kV]
        real(rp) :: wMAG
        !! Energy partition, Mag E / (Kin + Mag + Thermal). 0 <= wMAG <= 1
        !! wMAG = 1 => magnetically dominated, wMAG << 1 => not strongly magnetized
        real(rp) :: rCurv
        !!Curvature radius [Rx] @ bmin if CLOSED
        real(rp) :: avgBeta
        !! Flux-tube volume-weighted average beta
        real(rp), dimension(0:MAXTUBEFLUIDS) :: avgP,avgN,stdP,stdN
        !! Average and standard deviation of pressure [nPa] and number density [#/cc] on field line
        real(rp), dimension(0:MAXTUBEFLUIDS) :: bcN,bcP
        !! Density (#/cc) and pressure (nPa) near inner MHD boundary on this field line
        real(rp) :: losscone,lossconec
        !! Size of losscone [RAD] at footpoint and conjugate (if CLOSED)
        real(rp) :: TioTe0
        !! Empirical Ti/Te estimate if it exists
        integer :: nTrc
        !! Number of points traced for this tube
    end type Tube_T

    type TubeShell_T
        type(ShellGridVar_T), dimension(NDIM) :: xyz0
            !! XYZ of seed point [Rx]
        type(ShellGridVar_T) :: lat0, lon0, invlat
            !! Lat/lon of tube footpoint and invariant latitude [rad]
        type(ShellGridVar_T) :: latc, lonc
            !! Conjugate lat/lon
        type(ShellGridVar_T) :: J11i
            !! Ionospheric FAC [muA/m2]
        type(ShellGridVar_T) :: topo
            !! topo = TUBE_CLOSED,OPEN,UNDEF
            !! Note: integer data cast to real precision
        type(ShellGridVar_T) :: bmin
            !! Magnitude B at bmin location
        type(ShellGridVar_T), dimension(NDIM) :: X_bmin
            !! (NDIM) xyz-SM location of bmin surface
        type(ShellGridVar_T), dimension(NDIM) :: X_maxRad
            !! (NDIM) xyz-SM location of max radius along field line
        type(ShellGridVar_T), dimension(NDIM) :: Jxyz_iono
            !! NDIM current just above ionosphere boundary
            !! Note: magnitude is from sampled point, same location as Bxyz_iono
        type(ShellGridVar_T), dimension(NDIM) :: Bxyz_iono
            !! NDIM magnetic field just above ionosphere boundary
        type(ShellGridVar_T) :: bVol
            !! Flux tube volume [Rx/nT]
        type(ShellGridVar_T) :: Lb
            !! Length of field line [Rx]
        type(ShellGridVar_T) :: Tb
            !! Alfven bounce time [s] if CLOSED
        type(ShellGridVar_T) :: wMAG
            !! Energy partition, Mag E / (Kin + Mag + Thermal). 0 <= wMAG <= 1
            !! wMAG = 1 => magnetically dominated, wMAG << 1 => not strongly magnetized
        type(ShellGridVar_T) :: rCurv
            !! Curvature radius [Rx] @ bmin if CLOSED
        type(ShellGridVar_T) :: avgBeta
            !! Flux-tube volume-weighted average beta
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: avgP
            !! (N_MHDspc) average plasma pressure [nPa]
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: avgN
            !! (N_MHDspc) average plasma density [#/cc]
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: stdP
            !! (N_MHDspc) Standard deviation of pressure along field line
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: stdN
            !! (N_MHDspc) Standard deviation of density along field line
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: bcN
            !! (N_MHDspc) Near-BC density (#/cc)
        type(ShellGridVar_T), dimension(0:MAXTUBEFLUIDS) :: bcP
            !! (N_MHDspc) Near-BC pressure (nPa)
        type(ShellGridVar_T) :: losscone
            !! Size of losscone [RAD] at footpoint (if CLOSED)
        type(ShellGridVar_T) :: lossconec
            !! Size of losscone [RAD] at conjugate footpoint (if CLOSED)
        type(ShellGridVar_T) :: TioTe0
            !! Empirical Ti/Te estimate if it exists
        type(ShellGridVar_T) :: nTrc
            !! Number of points traced for each tube
    end type TubeShell_T

    type IMAGTube_T
        real(rp) :: Vol,bmin,beta_average
            !! Flux tube volume, minimum b along FL
            !! Average plasma beta
        real(rp), dimension(0:MAXTUBEFLUIDS) :: Pave, Nave
            !! average pressure, average density
        real(rp), dimension(0:MAXTUBEFLUIDS) :: Pstd, Nstd
            !! standard deviation of pressure and densiy along the field line
        real(rp) :: X_bmin(NDIM)
        integer(ip) :: topo
        real(rp) :: latc,lonc !Conjugate lat/lon
        real(rp) :: Lb, Tb !Arc length/bounce time
        real(rp) :: losscone,rCurv,wIMAG,TioTe0=4.0_rp
        real(rp) :: Veb  ! [km/s]
    end type IMAGTube_T


    ! =========================================================
    ! ionoCoupler_T — abstract ionosphere coupler interface
    ! Mirrors imagCoupler_T / raijuCoupler_T pattern.
    ! Any ionospheric model (REMIX, SMITE, ...) extends this.
    ! =========================================================

    type, extends(BaseOptions_T) :: ionoOptions_T
        type(ShellGrid_T) :: voltGrid
            !! Copy of Voltron ShellGrid passed during init

        contains
    end type ionoOptions_T

end module voltCplTypes