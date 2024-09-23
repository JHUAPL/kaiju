!Generic output clock for various models
module ioclock
    use kdefs
    use xml_input
    use clocks
    implicit none

    type IOClock_T
        !Output/restart (cadence/next output)
        real(rp) :: dtOut,tOut ! output timing info
        real(rp) :: dtRes,tRes ! restart timing info
        real(rp) :: dtCon,tCon ! console output timing info
        integer :: nOut=0,nRes=0 !Output numbering for output/restarts

        logical :: doResOut,doConOut,doDataOut,doTimerOut !Logical flags to do various outputs
        logical :: doFat,doSlim !Not yet used generically
        
        contains

            procedure :: init => IOClockInit
            procedure :: doConsole => doConsoleIOClock
            procedure :: doOutput  => doOutputIOClock
            procedure :: doRestart => doRestartIOClock
            procedure :: doTimer => doTimerIOClock
            procedure :: nextIOTime => doNextIOTime

    end type IOClock_T

    contains

    !Initialize IOClock from an already created XML reader
    subroutine IOClockInit(this,iXML,time,ts,isResO)
        class(IOClock_T), intent(inout) :: this
        type(XML_Input_T), intent(in)   :: iXML
        real(rp), intent(in) :: time
        integer, intent(in) :: ts
        logical, intent(in), optional :: isResO
        
        logical :: slimTimers
        logical :: isRes

        if(present(isResO)) then
            isRes = isResO
        else
            isRes = .false.
        endif

        if(iXML%Exists("output/tsOut")) then
            write (*,*) "Please update the XML to use output/dtCon for console output"
            write (*,*) "Instead of the no longer used output/tsOut"
        endif

        call iXML%Set_Val(this%dtCon,'output/dtCon' ,5.0)
        call iXML%Set_Val(this%dtOut,'output/dtOut' ,10.0)
        call iXML%Set_Val(this%dtRes,'restart/dtRes',100.0)

        !Get booleans
        call iXML%Set_Val(this%doTimerOut,   'output/doTimer'   ,.false.)
        call iXML%Set_Val(this%doFat,        'output/doFat'     ,.false.)
        call iXML%Set_Val(this%doSlim,       'output/doSlim'    ,.false.)
        call iXML%Set_Val(slimTimers,        'output/slimTimers',.true.)
        if(this%doTimerOut) slimTimers = .false.
        call setSlimTimers(slimTimers)

        if (this%dtCon<0) then
            this%doConOut = .false.
            this%tCon = HUGE
        else
            this%doConOut = .true.
            this%tCon = time
        endif

        !NOTE: Setting it so that first output is @ time
        if (this%dtOut<0) then
            this%doDataOut = .false.
            this%tOut = HUGE
        else
            this%doDataOut = .true.
            if(.not. isRes) this%tOut = time
        endif

        if (this%dtRes<0) then
            this%doResOut = .false.
            this%tRes = HUGE
        else
            this%doResOut = .true.
            if(.not. isRes) this%tRes = time
        endif

    end subroutine IOClockInit
    
    !Check for various output times
    function doConsoleIOClock(this,t)
        class(IOClock_T), intent(in) :: this
        real(rp), intent(in) :: t
        logical :: doConsoleIOClock

        if (t>=this%tCon .and. this%doConOut ) then
            doConsoleIOClock = .true.
        else
            doConsoleIOClock = .false.
        endif

    end function doConsoleIOClock

    function doRestartIOClock(this,t)
        class(IOClock_T), intent(in) :: this
        real(rp), intent(in) :: t
        logical :: doRestartIOClock

        if ( (t>=this%tRes) .and. this%doResOut ) then
            doRestartIOClock = .true.
        else
            doRestartIOClock = .false.
        endif
    end function doRestartIOClock

    function doOutputIOClock(this,t)
        class(IOClock_T), intent(in) :: this
        real(rp), intent(in) :: t
        logical :: doOutputIOClock

        if ( (t>=this%tOut) .and. this%doDataOut ) then
            doOutputIOClock = .true.
        else
            doOutputIOClock = .false.
        endif
    end function doOutputIOClock

    function doTimerIOClock(this,t)
        class(IOClock_T), intent(in) :: this
        real(rp), intent(in) :: t
        logical :: doTimerIOClock

        if (t>=this%tCon) then
            doTimerIOClock = .true.
        else
            doTimerIOClock = .false.
        endif

    end function doTimerIOClock

    ! Returns the time sim time when sime kind of IO should occur
    function doNextIOTime(this)
        class(IOClock_T), intent(in) :: this
        real(rp) :: doNextIOTime

        doNextIOTime = HUGE
        ! console output time is a prediction based on current DT and timesteps until output
        if(this%doConOut) doNextIOTime = min(doNextIOTime, this%tCon)
        if(this%doResOut) doNextIOTime = min(doNextIOTime, this%tRes)
        if(this%doDataOut) doNextIOTime = min(doNextIOTime, this%tOut)
        ! don't check for clock cleaning timer, it's not important enough (?)

    end function doNextIOTime

    !Copy ioA=>ioB using tScl scaling
    subroutine IOSync(ioA,ioB,tSclO)
        class(IOClock_T), intent(in   ) :: ioA
        class(IOClock_T), intent(inout) :: ioB
        real(rp), intent(in), optional :: tSclO

        real(rp) :: tScl

        if (present(tSclO)) then
            tScl = tSclO
        else
            tScl = 1.0
        endif

        ioB%tOut = ioA%tOut*tScl
        ioB%tRes = ioA%tRes*tScl
        ioB%tCon = ioA%tCon*tScl

        ioB%dtOut = ioA%dtOut*tScl
        ioB%dtRes = ioA%dtRes*tScl
        ioB%dtCon = ioA%dtCon*tScl

    end subroutine IOSync

end module ioclock
