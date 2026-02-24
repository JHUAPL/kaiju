! code for helping perform dynamic load balancing

module loadBalance
    use kdefs
    use helpertypes
    
    implicit none

    contains

    subroutine createLoadBalancer(lb, nL, nO, plusOne)
        type(loadBal_T), intent(inout) :: lb
        integer, intent(in) :: nL
        integer, intent(in) :: nO
        logical, optional, intent(in) :: plusOne ! whether instantLoads should include
        ! one extra space for root timing info that won't be used

        logical :: po = .false.

        if(present(plusOne)) po = plusOne

        lb%nL = nL
        lb%nO = nO

        if(allocated(lb%instantTimes)) deallocate(lb%instantTimes)
        if(po) then
            allocate(lb%instantTimes(lb%nL+1))
        else
            allocate(lb%instantTimes(lb%nL))
        endif

        if(allocated(lb%smoothLoads)) deallocate(lb%smoothLoads)
        allocate(lb%smoothLoads(lb%nL))

        if(allocated(lb%balStartInd)) deallocate(lb%balStartInd)
        allocate(lb%balStartInd(lb%nL))

        ! starting defaults
        lb%instantTimes(:) = 1000.0_rp
        lb%smoothLoads(:) = 1000.0_rp

        call updateLoads(lb)

    end subroutine createLoadBalancer

    subroutine updateLoads(lb)
        type(loadBal_T), intent(inout) :: lb

        real(rp) :: targetTime, defaultLoad, loadSum, thisPercent, totalPercent
        integer :: l

        ! this assumes that lb%instantLoads has new data since the last call

        ! amount of time each worker should take
        if(size(lb%instantTimes) == lb%nL) then
            targetTime = sum(lb%instantTimes)/lb%nL
        else
            ! for the case where helperTimes includes an extra space for root timing
            ! which isn't used
            targetTime = sum(lb%instantTimes(2:))/lb%nL
        endif

        ! weighted average of previous balanced loads and calculated load adjustment
        ! if a worker's instantTime is below targetTime, it's multiplier is > 1 and its work increases
        ! if a worker's instantTime is above targetTime, it's multiplier is < 1 and its work decreases
        ! this should be approximately zero sum
        if(size(lb%instantTimes) == lb%nL) then
            lb%smoothLoads = lb%smoothLoads*(lb%hAlpha + (1.0_rp-lb%hAlpha)*(targetTime/lb%instantTimes))
        else
            lb%smoothLoads = lb%smoothLoads*(lb%hAlpha + (1.0_rp-lb%hAlpha)*(targetTime/lb%instantTimes(2:)))
        endif

        ! use the mean of loads as the default
        defaultLoad = sum(lb%smoothLoads)/lb%nL

        ! calculate the relative size of each block
        loadSum = 0.0_rp
        do l=1,lb%nL
            if(lb%smoothLoads(l) <= 0) then
                ! we can't have a worker with no work
                loadSum = loadSum + defaultLoad
            else
                ! this input is valid
                loadSum = loadSum + lb%smoothLoads(l)
            endif
        enddo

        ! use the calculated percentage of each block to assign starting indices
        totalPercent = 0.0_rp
        do l=1,lb%nL
            lb%balStartInd(l) = totalPercent * lb%nO

            if(lb%smoothLoads(l) <= 0) then
                ! using default
                thisPercent = defaultLoad / loadSum
            else
                thisPercent = lb%smoothLoads(l) / loadSum
            endif
            ! thisPercent not permitted to be less than 0%
            thisPercent = MAX(0.0_rp, thisPercent)
            totalPercent = totalPercent + thisPercent
            ! totalPercent not permitted to increase past 100%
            totalPercent = MIN(1.0_rp, totalPercent)
        enddo

        ! ensure that every worker gets at least one index to work on, regardless of load
        do l=2,lb%nL
            ! loop to ensure every worker starts at least one index above the previous one
            !      but not outside the allowed index range minus one (so the last worker has work, too)
            lb%balStartInd(l) = min(max(lb%balStartInd(l),lb%balStartInd(l-1)+1),lb%nO-1)
        enddo
        do l=lb%nL-1,1,-1
            ! loop to ensure every worker starts at least one index below the next one
            !      but not outside the allowed index range
            lb%balStartInd(l) = max(min(lb%balStartInd(l),lb%balStartInd(l+1)-1),0)
        enddo

        ! test to ensure that all calculated start indices are greater than the one before them
        do l=1,lb%nL-1
            if(lb%balStartInd(l) >= lb%balStartInd(l+1)) then
                write (*,*) "Load Balancing Failed, and assigned equal or invalid starting indices. Exitting."
                stop
            endif
        enddo


    end subroutine updateloads

end module loadBalance

