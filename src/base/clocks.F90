!Timer construct for Kaiju routines

module clocks
    use kdefs
    use strings
    implicit none

    !Global clock parameters
    integer, parameter :: maxClocks = 80
    integer :: clockRate=0,clockMax=0

    !Clock output (min/max depth)
    integer :: clkMinD = 1, clkMaxD = 5

    real(rp) :: clkMinT = 0.0 !Minimum time (in %) to display in clock output

    !Timer construct
    type Timer_T
        character(len=strLen) :: tName = 'aTimer' !Timer ID
        logical :: isOn = .false.
        integer :: iTic=0,iToc=0 !Integer ticks
        real(rp) :: tElap=0.0 !Elapsed time
        integer :: clockRate=0, clockMax=0, nCalls=0
        integer :: level=-1, parent=-1

        !Routines
        contains
        procedure, pass :: tic => startTimer
        procedure, pass :: toc => stopTimer
        procedure, pass :: reset => resetTimer

    end type Timer_T

    !Clock construct
    type Clock_T
        character(len=strLen) :: cID = 'NULL' !Clock ID
        logical :: isOn = .false.
        !Depth/parent, ie array entry of parent of this clock
        integer :: level=-1,parent=-1
        integer :: iTic=0,iToc=0 !Integer ticks
        real(rp) :: tElap=0.0 !Elapsed time
    end type Clock_T

    !Global clock array
    type(Clock_T), dimension(maxClocks) :: kClocks
    integer, private :: nclk=0 !Current number of clocks
    logical, private :: slimTimers = .false.

    !interface for reading clock time
    interface readClock
        module procedure readClock_str, readClock_int
    end interface

    contains

    !Initialize gobal clock array
    subroutine initClocks()
        integer :: n,itoc

        !Create "Omega" as overall parent clock
        kClocks(1)%cID = "Omega"
        kClocks(1)%level = 0
        kClocks(1)%parent = 1 !Omega is its own parent

        nclk = 1
        !Get clock info
        call system_clock(itoc,clockRate,clockMax)
    end subroutine initClocks

    subroutine setSlimTimers(isSlim)
        logical, intent(in) :: isSlim

        slimTimers = isSlim

    end subroutine

    !Start clock given by cID, create it if necessary
    subroutine Tic(cID,priority)
        character(len=*), intent(in) :: cID
        logical, optional, intent(in) :: priority

        integer :: n,iblk
        logical :: lPriority
        iblk = 0

        if(present(priority)) then
            lPriority = priority
        else
            lPriority = .false.
        endif

        ! only process priority times when doing slim timing
        if(slimTimers .and. .not. lPriority) return

        !Find timer
        do n=1,nclk
            if (toUpper(kClocks(n)%cID) == toUpper(cID)) then
                !Found it, save ID
                iblk = n
            endif
        enddo

        if (iblk == 0) then
            !Not found, create new timer
            call newClock(cID)
            iblk = nclk
        endif

        !Start timer in slot iblk
        kClocks(iblk)%isOn = .true.
        call system_clock(count=kClocks(iblk)%iTic)
    end subroutine Tic

    !Add new clock
    subroutine newClock(cID)
        character(len=*), intent(in) :: cID
        integer :: Npa(1),Np,Nl

        if (nclk == maxClocks-1) then
            write(*,*) 'Ran out of clocks, increse maxClocks ...'
            stop
        endif

        !Find parent of this clock, ie on clock w/ largest level
        kClocks(nclk+1)%cID = trim(cID)
        Npa = maxloc(kClocks(1:nclk)%level,mask=kClocks(1:nclk)%isOn)
        Np = max(Npa(1),1)
        Nl = kClocks(Np)%level+1
        kClocks(nclk+1)%parent = Np
        kClocks(nclk+1)%level = Nl

        !write(*,*) "Adding clock ",cID," w/ parent/level = ", Np,Nl

        nclk = nclk+1
    end subroutine newClock

    !Stop clock, save time.  Error if clock doesn't exist
    subroutine Toc(cID, priority)
        character(len=*), intent(in) :: cID
        logical, optional, intent(in) :: priority

        integer :: n,iblk
        real :: wclk
        logical :: lPriority
        iblk = 0

        if(present(priority)) then
            lPriority = priority
        else
            lPriority = .false.
        endif

        ! only process priority times when doing slim timing
        if(slimTimers .and. .not. lPriority) return

        !Find timer
        do n=1,nclk
            if (toUpper(kClocks(n)%cID) == toUpper(cID)) then
                !Found it, save ID
                iblk = n
            endif
        enddo

        if (iblk == 0) then
            write(*,*) '<Clock Error, Stopping Non-existent Clock!>'
            iblk = 1 !Avoid crash
        endif

        kClocks(iblk)%isOn = .false.
        call system_clock(count=kClocks(iblk)%iToc)
        wclk = real(kClocks(iblk)%iToc-kClocks(iblk)%iTic)/real(clockRate)

        kClocks(iblk)%tElap = kClocks(iblk)%tElap + wclk
    end subroutine Toc

    !Reset clocks
    subroutine cleanClocks()
        integer :: n

        do n=1,nclk
            kClocks(n)%tElap = 0
            ! if the clock is active, reset the tic to right now
            if(kClocks(n)%isOn) call Tic(kClocks(n)%cID, .true.)
        enddo

    end subroutine cleanClocks

    !Get specific elapsed time from clock w/ ID cID
    function readClock_str(cID) result(wclk)
        character(len=*), intent(in) :: cID

        integer :: n,iblk
        real :: wclk

        iblk = 0
        !Find timer
        do n=1,nclk
            if (toUpper(kClocks(n)%cID) == toUpper(cID)) then
                !Found it, save ID
                iblk = n
            endif
        enddo
        
        wclk = readClock_int(iblk)

    end function readClock_str

    function readClock_int(iblk) result(wclk)
        integer, intent(in) :: iblk

        integer :: tmpToc
        real :: wclk

        if (iblk == 0) then
            wclk = 0.0
        else
            if (kClocks(iblk)%isOn) then
                call system_clock(count=tmpToc)
                wclk = kClocks(iblk)%tElap + real(tmpToc-kClocks(iblk)%iTic)/real(clockRate)
            else
                wclk = kClocks(iblk)%tElap
            endif
        endif
    end function readClock_int

    !Output clock data
    subroutine printClocks()
        integer :: n,l
        real(rp) :: tScl
        character(len=strLen) :: tStr

        ! do n=1,nclk
        !     write(*,*) 'ID/T = ', trim(kClocks(n)%cID),kClocks(n)%tElap
        ! enddo
        !Start w/ Omega
        tScl = getClockTime(1)
        write(tStr,'(f8.3)') tScl

        write(*,'(a)') 'Timing Data (' // trim(adjustl(tStr)) // 's)'
        do n=1,maxClocks
            !Loop until find clock at min depth
            l = kClocks(n)%level
            if (l == clkMinD) then
                call printChildren(n)
            endif
        enddo
        write(*,*) ''

        contains
        !Print pID and children up to max depth
            recursive subroutine printChildren(pID)
                integer, intent(in) :: pID
        
                integer :: n,nP,nC,nInd,lev
                real(rp) :: tF,tS
                character(len=StrLen) :: sPer,sSec
        
                lev = kClocks(pID)%level
                nInd = (lev - clkMinD)+1
                
                tS = getClockTime(pID)
                tF = 100*tS/tScl
        
                !if ( (tF > clkMinT) .and. (lev <= clkMaxD) ) then
                    write(sPer,'(f6.2)') tF
                    write(sSec,'(f8.3)') tS
                    write(*,'(a)') repeat("   ",nInd) // trim(kClocks(pID)%cID) // ' : ' // &
                        trim(adjustl(sPer)) // '% / ' // trim(adjustl(sSec)) // 's'
                    do n=1,maxClocks
                        nP = kClocks(n)%parent
                        if (np == pID) then
                            call printChildren(n)
                        endif
                    enddo
                !endif
        
                !write(*,*) trim(gTime(pID)%tName), ' : ', gTime(pID)%tElap
            end subroutine printChildren

            function getClockTime(clockIndex)
                integer, intent(in) :: clockIndex
                real(rp) :: getClockTime,wclk
                integer :: curT

                getClockTime = kClocks(clockIndex)%tElap
                if(kClocks(clockIndex)%isOn) then
                    ! clock is active, include elapsed time
                    call system_clock(count=curT)
                    wclk = real(curT-kClocks(clockIndex)%iTic)/real(clockRate)
                    getClockTime = getClockTime + wclk
                endif

            end function getClockTime
    end subroutine printClocks

    !Basic clock routines
    subroutine resetClocks(nClocks)
        type(Timer_T), intent(inout), dimension(maxClocks) :: nClocks
        integer :: n

        do n=1,maxClocks
            call nClocks(n)%reset()
        enddo
    end subroutine resetClocks
    
    !Starts timer
    !Using system_clock to correctly handle threads
    !(as opposed to cpu_time)
    subroutine startTimer(self)
        implicit none
        class(Timer_T) :: self

        self%isOn = .true.
        call system_clock(self%iTic,self%clockRate,self%clockMax)

    end subroutine startTimer

    !End timer, print if optional flag is T
    subroutine stopTimer(self,optPrint)
        implicit none
        class(Timer_T) :: self
        logical, intent(in), optional :: optPrint

        logical :: doPrint
        doPrint = .false.
        if (present(optPrint)) then
            doPrint = optPrint
        endif
        !TODO Add check for calling stop w/o start
        self%isOn = .false.
        self%nCalls = self%nCalls+1

        call system_clock(self%iToc,self%clockRate,self%clockMax)
        self%tElap = self%tElap + real(self%iToc-self%iTic)/real(self%clockRate)
        
        if (doPrint) then
            write(*,*) 'Timer (', trim(self%tName), ') = ', self%tElap
        endif

    end subroutine stopTimer

    !Reset timer
    subroutine resetTimer(self)
        implicit none
        class(Timer_T) :: self

        self%nCalls = 0
        self%isOn = .false.
        self%tElap = 0.0

    end subroutine resetTimer

end module clocks
