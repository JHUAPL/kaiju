!Various time-series utilities to read in variables from a file and return its value at any given time

module kronos
    use ioH5,     only : IOVAR_T, ClearIO, AddInVar, ReadVars, ioExist
    use strings,  only : toUpper
    use kdefs,    only : strLen, dp, rp
 
    implicit none

    ! Global Parameters
    integer, parameter :: MAXWINDVARS = 20

    !Type for generic solar wind BC (from file or subroutine)
    !Either use discrete tF,Qw(NVAR) series or subroutine
    type TimeSeries_T
        !File ID string
        character(len=strLen) :: wID   !File name
        character(len=strLen) :: varID !Variable name
        logical :: isDiscrete=.false.

        !Holder for discrete data
        real(rp), allocatable :: tF(:), Q(:)
        real(rp) :: tMin,tMax
        integer :: NumT

        contains

        procedure :: initTS
        procedure :: getValue
        procedure :: getMax
        
    end type TimeSeries_T


    contains

    !Read variable data from file
    subroutine initTS(varTS,varStr,doLoudO)
        class(TimeSeries_T), intent(inout) :: varTS
        character(len=*),intent(in) :: varStr
        logical, intent(in), optional :: doLoudO

        integer :: N
        logical :: fExist,doLoud,vEx1,vEx2
        type(IOVAR_T), dimension(MAXWINDVARS) :: IOVars

        varTS%varID = varStr

        if (present(doLoudO)) then
            doLoud = doLoudO
        else
            doLoud = .true.
        endif

        if (trim(toUpper(varTS%wID)) .eq. "NONE") then
            write(*,*) "Error: No input file specified to be read in"
            stop
        endif

        !Check file
        inquire(file=trim(varTS%wID),exist=fExist)
        if (.not. fExist) then
            write(*,*) "Error reading ", trim(varTS%wID), " exiting ..."
            stop
        endif

        !Check variables
        vEx1 = ioExist(varTS%wID,"T")
        vEx2 = ioExist(varTS%wID,varTS%varID)
              
        if (.not. (vEx1 .and. vEx2) ) then
            write(*,*) 'Error reading time series variables from ', trim(varTS%wID), ' : T / ', trim(varTS%varID)
            stop
        endif

        if (doLoud) then
            write(*,*) "---------------"
            write(*,*) "Creating time series for ", trim(varTS%varID)
            write(*,*) "Reading data from ", trim(varTS%wID)
        endif

        !Set discrete wind function                                                                                                                                                            
        varTS%isDiscrete = .true.

        !Setup input chain
        call ClearIO(IOVars)
        call AddInVar(IOVars,"T")
        call AddInVar(IOVars,trim(varTS%varID))

        !Read data, don't use IO precision
        call ReadVars(IOVars,.false.,varTS%wID)
        N = IOVars(1)%N
        

        !Allocate holders
        varTS%NumT = N
        if(allocated(varTS%tF)) deallocate(varTS%tF)
        if(allocated(varTS%Q)) deallocate(varTS%Q)
        allocate(varTS%tF(N))
        allocate(varTS%Q(N))

        !Not normalized, must be put into code units after value is returned
        varTS%tF            = IOVars(1)%data
        varTS%Q             = IOVars(2)%data

        varTS%tMin = minval(varTS%tF)
        varTS%tMax = maxval(varTS%tF)

        if (doLoud) then
            write(*,*) "Finished reading data for ", trim(varTS%varID)
            write(*,*) "---------------"
        endif

    end subroutine initTS

    ! Interpolate between discrete points given by the file
    ! to provide a value at any given time
    ! Variable is returned without being normalized to code units
    ! Must be converted after value is returned
    subroutine getValue(varTS,t,Var)
        class(TimeSeries_T), intent(in) :: varTS
        real(rp), intent(in) :: t
        real(rp), intent(out) :: Var

        integer :: i0,i1
        real(rp) :: w0,w1,dT
        Var = 0.0

        if (t >= varTS%tMax) then
            i0 = varTS%NumT
            i1 = i0
            w0 = 1.0
            w1 = 0.0

        else if (t <= varTS%tMin) then
            i0 = 1
            i1 = i0
            w0 = 1.0
            w1 = 0.0
        else
            i0 = maxloc(varTS%tF,dim=1,mask=varTS%tF .le. t)
            i1 = i0+1
            dT = varTS%tF(i1)-varTS%tF(i0)
            w0 = (varTS%tF(i1)-t)/dT
            w1 = (t-varTS%tF(i0))/dT
        endif

        Var = w0*varTS%Q(i0) + w1*varTS%Q(i1)
        
    end subroutine getValue

    !Function to return value (just wraps getValue)
    function evalTimeSeriesAt(varTS,t) result(tVal)
        type(TimeSeries_T), intent(in) :: varTS
        real(rp), intent(in) :: t
        real(rp) :: tVal

        call varTS%getValue(t,tVal)
    end function evalTimeSeriesAt

    !Return maximum value of time series
    function getMax(varTS) result(tsMax)
        class(TimeSeries_T), intent(in) :: varTS
        real(rp) :: tsMax

        tsMax = maxval(varTS%Q)

    end function getMax

end module kronos
