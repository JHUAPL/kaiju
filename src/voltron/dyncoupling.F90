!Routines to dynamically alter coupling cadence
module dyncoupling
    
    use volttypes
    use kronos
    use earthhelper, ONLY : SetKp0
    implicit none
    enum, bind(C)
        enumerator :: NULLALERT=0,WHITEALERT,GREYALERT,YELLOWALERT,REDALERT
    end enum

    logical, private :: doInit = .true.
    type(TimeSeries_T), private :: KpTS
    integer, private :: alertStatus = NULLALERT

    contains

    !Do checks for conditions to change coupling cadence
    subroutine UpdateCouplingCadence(vApp)
        class(voltApp_T), intent(inout) :: vApp

        integer :: n,KpI,newAlert
        real(rp) :: t0,t,KpMax,dynDT
        KpMax = -1.0
        if ((.not. vApp%doDynCplDT) .or. (vApp%time<=0)) return

        if (doInit) then
            !Read Kp time series
            KpTS%wID = trim(vApp%tilt%wID)
            call KpTS%initTS("Kp",doLoudO=.false.)
            doInit = .false.
        endif

        t0 = vApp%time
        !Loop over 15min increments -1/+1 hour, find max Kp
        do n=-4,+4
            t = t0 + 15.0*60.0*n
            KpMax = max(KpMax,KpTS%evalAt(t))
        enddo

        KpI = nint(KpMax) !Cast to integer
        !Set coupling cadences based on max Kp on one hour window
        select case (KpI)
            case (1,2)
                dynDT = 15.0
                newAlert = WHITEALERT
            case (3,4)
                dynDT = 15.0
                newAlert = GREYALERT
            case (5,6)
                dynDT = 10.0
                newAlert = YELLOWALERT
            case (7:)
                dynDT = 10.0
                newAlert = REDALERT
            case default
                dynDT = 15.0

        end select

        !Whether to reset current Kp default
        !call SetKp0(KpI)
        if (newAlert /= alertStatus) then
            !Alert status has changed
            if (newAlert < alertStatus) then
                write(*,*) 'Downgrading cadence to ', dynDT
            else
                write(*,*) 'Upgrading cadence to ', dynDT
            endif
            
            call resetDeepCoupling(vApp,dynDT)
            alertStatus = newAlert
        endif

    end subroutine UpdateCouplingCadence

    !Helper routines to reset coupling cadences
    subroutine resetDeepCoupling(vApp, newDeep)
        class(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: newDeep

        ! temporary unsupported, needs to be modified to support new coupling
        !vApp%TargetDeepDT = newDeep
    end subroutine resetDeepCoupling


end module dyncoupling
