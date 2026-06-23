!CHIMP definitions/constants

module chmpdefs
    use kdefs
    use math
    use xml_input
    use ioclock
    implicit none

    !Algorithmic/Run options
    type chmpModel_T
        character(len=strLen) :: RunID, uID
        real(rp) :: t,T0,tFin,dt !Various times
        !Mass/charge in units of me and |e|
        real(rp) :: m0,q0
        integer :: imeth !Integration method (IFO,IGC,IDYN)
        logical :: isDynamic=.false. !Whether or not using dynamic integrator
        integer :: nTh=1 !Number of threads per node/group (currently unused)
        integer :: ts,tsOut,nOut
        logical :: doTimer=.false.
        logical :: doPureB0=.false.,doNumB0=.false. !Background field options (Pure/Numerical)
        logical :: doEBOut=.false.,doTPOut=.false.,doFLOut=.false. !slice/particle/field line output
        logical :: doSlim=.false.,doFat=.false. !Do slim/fat output
        logical :: doLL=.false. !Output lat-lon projection
        logical :: do2D=.false. !Force 2D tp integration
        logical :: doMHD=.false. !Do full MHD variables instead of just E/B
        logical :: doJ=.false. !Deal w/ currents in addition to just E/B
        logical :: doReadSrc=.false. ! Whether or not to read src data from file
        real(rp) :: tOut,dtOut
        integer :: Nblk=1,Ng=1
        real(rp) :: epsht !Small number for timestep calculation
        real(rp) :: epsgc !Max value for adiabaticity parameter
        real(rp) :: epsds !Small number for tracing
        real(rp) :: MinK = 1.0e-2 !Min K [keV] to kill GC particles
        logical :: doEBInit=.false. !Initialize particles to have specified energy in ExB frame (instead of lab)
        logical :: doEBFix=.true. !Enforce E.B=0

        integer  :: MaxIter !Maximum iterations
        real(rp) :: TolGC !Tolerance for GC iteration

        logical :: doOldNaming = .false. !Whether to use old-style naming
        logical :: doEQProj=.false. !Force projection to EQ before output
        logical :: doTrc=.false. !Do field line topology tracing
        logical :: doStream=.false. !Inject TPs over time
        logical :: doWPI=.false. !Include wave particle interactions
        logical :: doEQScat=.false. !Do equatorial PA scattering
        real(rp) :: reqScat
        logical :: doPP=.false. !Include plasmapaus location calculations

        logical :: isRestart=.false.
        logical :: isMAGE = .false.
        logical :: doDip = .false.

        integer :: nSpc = 0 !Number of plasma species if doMHD=.true.
        !Output info
        type (IOClock_T) :: IO

        !Background field pointers
        procedure(B0_T)   , pointer, nopass :: B0
        procedure(JacB0_T), pointer, nopass :: JacB0
    end type chmpModel_T


    !Field types
    enum, bind(C)
        enumerator :: BFLD=1,EFLD,DBFLD
    endenum

    !Grid types
    enum, bind(C) 
        enumerator :: LFMGRID,EGGGRID,CARTGRID,SPHGRID
    endenum

    !Integrator types
    enum, bind(C)
        enumerator :: IFO,IGC,IDYN
    endenum
!Function pointer types for background fields
    !NOTE: Jb0(i,j) = B_{i,j}, ie derivative of B_i in direction j
    abstract interface
        function B0_T(r) result(Bxyz)
            Import :: rp,NDIM
            real(rp), intent(in) :: r(NDIM)
            real(rp) :: Bxyz(NDIM)

        end function B0_T
    end interface

    abstract interface
        function JacB0_T(r) result(Jb0)
            Import :: rp,NDIM
            real(rp), intent(in) :: r(NDIM)
            real(rp) :: Jb0(NDIM,NDIM)

        end function JacB0_T
    end interface

    contains
    
    !Generate dimension sampling from XML (w/ default bounds)
    subroutine getSample(Model,inpXML,xID,Np,Xs,xMinO,xMaxO,doLogO)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: inpXML
        character(len=*), intent(in) :: xID
        integer, intent(in) :: Np
        real(rp), allocatable, dimension(:), intent(inout) :: Xs
        real(rp), intent(in), optional :: xMinO,xMaxO
        logical, intent(in), optional :: doLogO

        character(len=strLen) :: minStr,maxStr,rStr
        real(rp) :: xMinD,xMaxD,xMin,xMax
        logical :: doLogD,doLog
        integer :: n

        !Set default bounds/random method
        xMinD = 1.0
        xMaxD = 10.0
        doLogD = .false.
        if (present(xMinO))  xMinD  =  xMinO
        if (present(xMaxO))  xMaxD  =  xMaxO
        if (present(doLogO)) doLogD = doLogO

        !Now get from XML and overwrite if necessary
        minStr = trim(xID) // "/min"
        maxStr = trim(xID) // "/max"
        rStr   = trim(xID) // "/doLog"
        call inpXML%Set_Val(xMin ,trim(minStr), xMinD)
        call inpXML%Set_Val(xMax ,trim(maxStr), xMaxD)
        call inpXML%Set_Val(doLog,trim(rStr)  ,doLogD)

        if (allocated(Xs)) deallocate(Xs)
        allocate(Xs(Np))

        if (abs(xMin-xMax)<TINY) then
            Xs(:) = xMin
            return
        endif
        
        do n=1,Np
            if (doLog) then
                Xs(n) = genRandLog(xMin,xMax)
            else
                Xs(n) = genRand(xMin,xMax)
            endif
        enddo

    end subroutine getSample

    !Generate dimension spacing (w/ default bounds)
    subroutine getDim(Model,inpXML,xID,Nc,Xi,Xc,Xd,xMinO,xMaxO,doLogO)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: inpXML
        character(len=*), intent(in) :: xID
        integer, intent(inout) :: Nc
        real(rp), allocatable, dimension(:), intent(inout) :: Xi,Xc,Xd
        real(rp), intent(in), optional :: xMinO,xMaxO
        logical, intent(in), optional :: doLogO

        character(len=strLen) :: minStr,maxStr,rStr,nStr
        real(rp) :: xMinD,xMaxD,xMin,xMax,dx,dxn
        logical :: doLogD,doLog
        integer :: n,NcD

        !Set default bounds/random method
        xMinD = 1.0
        xMaxD = 10.0
        doLogD = .false.
        NcD = Nc
        if (present(xMinO))  xMinD  =  xMinO
        if (present(xMaxO))  xMaxD  =  xMaxO
        if (present(doLogO)) doLogD = doLogO

        !Now get from XML and overwrite if necessary
        minStr = trim(xID) // "/min"
        maxStr = trim(xID) // "/max"
        rStr   = trim(xID) // "/doLog"
        nStr   = trim(xID) // "/N"
        call inpXML%Set_Val(xMin ,trim(minStr), xMinD)
        call inpXML%Set_Val(xMax ,trim(maxStr), xMaxD)
        call inpXML%Set_Val(doLog,trim(rStr)  ,doLogD)

        !Get number of cells
        call inpXML%Set_Val(Nc,trim(nStr),NcD)

        if (allocated(Xi)) deallocate(Xi)
        if (allocated(Xc)) deallocate(Xc)
        if (allocated(Xd)) deallocate(Xd)

        allocate(Xi(Nc+1)) !Cell interfaces
        allocate(Xc(Nc)) !Cell centers
        allocate(Xd(Nc)) !Cell sizes

        !Calculate interfaces
        if (doLog) then
            dx = (log10(xMax)-log10(xMin))/Nc
        else
            dx = (xMax-xMin)/Nc
        endif

        do n=1,Nc+1
            dxn = (n-1)*dx
            if (doLog) then
                Xi(n) = xMin*10**(dxn)
            else
                Xi(n) = xMin + dxn
            endif
        enddo

        !Calculate centers
        do n=1,Nc
            Xc(n) = 0.5*(Xi(n)+Xi(n+1))
            Xd(n) = Xi(n+1)-Xi(n)
        enddo

    end subroutine getDim

    !Does the obvious
    subroutine SetTrcEpsilon(Model,epsds)
        type(chmpModel_T), intent(inout) :: Model
        real(rp), intent(in) :: epsds
        
        Model%epsds = epsds

    end subroutine SetTrcEpsilon

end module chmpdefs
