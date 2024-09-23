!Main kaiju definitions/globals
module kdefs
    use, intrinsic :: iso_fortran_env
    implicit none

    !Define variable precisions
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128 !Quad precision
    !integer, parameter :: ip = INT64
    integer, parameter :: ip = INT32
    
    !Compute precision, IO precision
    integer, parameter :: rp = dp
    integer, parameter :: iop = sp !Precision for IO

!Globals
    real(rp), parameter :: PI  = 4.0D0 *atan(1.0D0)
    real(qp), parameter :: QPI = 4.0_qp*atan(1.0_qp)
    real(rp), parameter :: TINY=1.0D-12, HUGE=1.0D+12
    real(rp), parameter :: rad2deg = 180.0/PI
    real(rp), parameter :: deg2rad = PI/180.0
    
!Defaults
    integer, parameter :: strLen = 1000 !Default size for strings
    integer, parameter :: ALIGN = 64 !Compiler-directive alignment
    integer, parameter :: NDIM=3 !Number of dimensions, size of vector
    integer, parameter :: rseed0=31337 !Default random seed

!Physical constants
    !CGS Constants
    real(rp), parameter :: vc_cgs = 2.99792458D+10 ![cm/s], speed of light
    real(rp), parameter :: qe_cgs = 4.80320425D-10 ![CGS], electron charge
    real(rp), parameter :: Re_cgs = 6.3781D8     ![cm]    Earth's radius
    real(rp), parameter :: Me_cgs = 9.1093837015D-28 ![g] Electron mass
    real(rp), parameter :: Mp_cgs = 1.67262192369D-24 ![g] Proton mass
    real(rp), parameter :: G_cgs  = 6.6726D-8     ![cm^3/g/s^2], Gravitational constant (per NRL plasma formulary'21)

    !MKS Constants
    real(rp), parameter :: vc_mks = vc_cgs*(1.0D-2) ![m/s], Speed of light
    real(rp), parameter :: G_mks  = 6.6726D-11     ![m^3/kg/s^2], Gravitational constant (per NRL plasma formulary'21)

    !Helper conversions
    real(rp), parameter :: G2nT = 1.0D+5 !Gauss->nT
    real(rp), parameter :: G2T = 1.0D-4 !Gauss->T
    real(rp), parameter :: kev2J = 1.602176634D-16 !keV->J
    real(rp), parameter ::  ev2J = kev2J*(1.0D-3)  ! eV->J
    real(rp), parameter ::   J2K = 7.242971666663D+22 !J -> Kelvin
    real(rp), parameter :: kev2erg = kev2J*1.0D+7
    real(rp), parameter :: erg2kev = 1.0/kev2erg
    real(rp), parameter :: Re_km = Re_cgs*(1.0D-2)*(1.0D-3) !km

    !Misc
    real(rp), parameter :: Mu0     = 4*PI*1.0D-7 ! Tm/A
    real(rp), parameter :: Kbltz   = 1.380649D-16      ![cm^2 g /s^2/K=erg/K] Boltzmann constant
    real(rp), parameter :: mec2    = (Me_cgs*vc_cgs**2.0)*(1.0D-3)/kev2erg ! [MeV] electron rest mass
    real(rp), parameter :: heFrac  = 1.16D0       ! Accounts for 4% helium
    real(rp), parameter :: eCharge = 1.60217663D-19  ! Charge of electron
    !NOTE: dalton isn't precisely Mp b/c carbon binding energy business
    real(rp), parameter :: dalton  = 1.6605390666D-27  ! Mass unit [kg]

!Planetary constants
    !Earth
    !real(rp), parameter :: EarthM0g = 0.31 !Gauss, old LFM value
    real(rp), parameter :: EarthM0g = 0.2961737_rp !Gauss, Olsen++ 2000

    real(rp), parameter :: REarth = Re_cgs*1.0D-2 !m

    real(rp), parameter :: RionE  = 6.5_rp    ! Earth Ionosphere radius in 1000 km
    ! Earth corotation potential
    !real(rp), parameter :: EarthPsi0 = 92.4 ! !!OUTDATED, based on LFM magntic field strength
    !real(rp), parameter :: EarthPsi0 = 87.62  ! Calculated using B = 0.2961737 Gauss
    real(rp), parameter :: EarthPsi0 = REarth**2 * (2.0*PI/86400.0) * (EarthM0g*G2T) * 1.0D-3 ! [kV]
        ! Everyone should ideally get Psi0 from planet object, but things specific to Earth should still be able to rely on EarthPsi0
    
    !Saturn
    real(rp), parameter :: SaturnM0g = 0.21_rp !Gauss
    real(rp), parameter :: RSaturnXE = 9.5_rp  !Rx = X*Re
    !Jupiter
    real(rp), parameter :: JupiterM0g = 4.8_rp  !Gauss
    real(rp), parameter :: RJupiterXE = 11.0_rp !Rx = X*Re
    !Mercury
    real(rp), parameter :: MercuryM0g = 0.00345_rp  !Gauss
    real(rp), parameter :: RMercuryXE = 0.31397_rp  !Rx = X*Re
    !Neptune
    real(rp), parameter :: NeptuneM0g = 0.142_rp  !Gauss
    real(rp), parameter :: RNeptuneXE = 3.860_rp  !Rx = X*Re

!Helio constants
    real(rp), parameter :: Rsolar = 6.956D5    ! [km] Solar radius
    real(rp), parameter :: Msolar = 1.98847D30 ! [kg] Solar mass
    real(rp), parameter :: Tsolar_synodic = 27.28_rp ![days] synodic Tsolar

!Numbered accessors
    !Directions
    enum, bind(C) 
        enumerator :: IDIR=1, JDIR, KDIR
    endenum 
    enum, bind(C) 
        enumerator :: XDIR=1, YDIR, ZDIR
    endenum 

    !Conserved variables
    enum, bind(C)
        enumerator :: DEN=1,MOMX,MOMY,MOMZ,ENERGY
    endenum
    !Primitive variables
    enum, bind(C)
        enumerator :: VELX=MOMX,VELY,VELZ,PRESSURE
    endenum

    ! Directions, used by ShellGrid, Remix, calcdb, etc
    enum, bind(C)
        enumerator :: NORTH=1,SOUTH,EAST,WEST
    end enum

#ifdef USECOLORTEXT
!Color options for funsies
character, parameter :: ANSIESCAPE = char(27) !Escape character
integer, parameter :: ANSILEN = 5
character(ANSILEN), parameter :: &
    ANSIRED    = ANSIESCAPE // '[31m', &
    ANSIGREEN  = ANSIESCAPE // '[32m', &
    ANSIYELLOW = ANSIESCAPE // '[33m', &
    ANSIBLUE   = ANSIESCAPE // '[34m', &
    ANSIPURPLE = ANSIESCAPE // '[35m', &
    ANSICYAN   = ANSIESCAPE // '[36m', &
    ANSIWHITE  = ANSIESCAPE // '[37m', &
    ANSIRESET  = ANSIESCAPE // '[0m'
#else
!Fake values to avoid text
integer, parameter :: ANSILEN = 0
character(ANSILEN), parameter :: &
    ANSIRED    = "", &
    ANSIGREEN  = "", &
    ANSIYELLOW = "", &
    ANSIBLUE   = "", &
    ANSIPURPLE = "", &
    ANSICYAN   = "", &
    ANSIWHITE  = "", &
    ANSIRESET  = ""
#endif

    contains

    !Dump contents of /proc/self/status to output
    subroutine printProcessInfo()
        logical :: fileEnd
        character(len=strLen) :: line
        integer :: io

        OPEN(37, file='/proc/self/status', action='read')
        fileEnd = .false.
        do while(.not. fileEnd)
            read(37,'(A)',iostat=io) line
            if (io /= 0) then
                fileEnd = .true.
            else
                write (*,*) trim(line)
            endif
        enddo
        CLOSE(37)

    end subroutine printProcessInfo


    !Generate name of output file based on tiling
    function genName(caseName,Ri,Rj,Rk,i,j,k,useOldStyle) result(fName)
        character(len=*), intent(in) :: caseName
        integer, intent(in) :: Ri,Rj,Rk,i,j,k
        logical, optional, intent(in) :: useOldStyle
        character(len=strLen) :: fName

        character(len=strLen) :: fId,fHd

        ! need to split this into two if statements to ensure ordering
        if(present(useOldStyle)) then
            if(useOldStyle) then
                fName = genName_old(caseName,Ri,Rj,Rk,i,j,k)
                return
            endif
        endif

        fId = genRunId(caseName,Ri,Rj,Rk,i,j,k)

        write(fHd ,'(a,a)') trim(fId), '.gam.h5'
        fName = trim(fHd)
        !write(*,*) 'ijk / file = ',i,j,k,trim(fName)
    end function genName

    ! Generate runID based on tiling (output file without extension)
    function genRunId(caseName, Ri, Rj, Rk, i, j, k) result(fName)
        character(len=*), intent(in) :: caseName
        integer, intent(in) :: Ri,Rj,Rk,i,j,k
        character(len=strLen) :: fName

        character(len=strLen) :: fHd,fRn,fijk

        if (Ri > 1 .or. Rj > 1 .or. Rk > 1) then
            write(fHd ,'(a,a)') trim(caseName), '_'
            write(fRn ,'(I0.4,a,I0.4,a,I0.4,a)') Ri,'_',Rj,'_',Rk,'_'
            write(fijk,'(I0.4,a,I0.4,a,I0.4)') i-1,'_',j-1,'_',k-1

            fName = trim(fHd) // trim(fRn) // trim(fijk)
        else
            fName = trim(caseName)
        endif
    end function genRunId

    ! Generate old omega format name based on tiling
    function genName_old(caseName,Ri,Rj,Rk,i,j,k) result(fName)
        character(len=*), intent(in) :: caseName
        integer, intent(in) :: Ri,Rj,Rk,i,j,k
        character(len=strLen) :: fName

        character(len=strLen) :: fHd,fRn,fijk,fTl
        integer :: n

        if (Ri > 1 .or. Rj > 1 .or. Rk > 1) then
            n = (j-1) + (i-1)*Rj + (k-1)*Ri*Rj
            write(fHd ,'(a,a)') trim(caseName), '_'
            write(fRn ,'(I0.4,a,I0.4,a,I0.4,a)') Ri,'_',Rj,'_',Rk,'_'
            write(fijk,'(I0.4,a,I0.4,a,I0.4,a)') i-1,'_',j-1,'_',k-1,'_'
            write(fTl ,'(I0.12,a)') n,'.h5'

            fName = trim(fHd) // trim(fRn) // trim(fijk) // trim(fTl)
        else
            fName = caseName
        endif
        !write(*,*) 'ijk / file = ',i,j,k,trim(fName)
    end function genName_old

end module kdefs
