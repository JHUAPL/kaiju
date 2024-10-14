module output
    use gdefs
    use gamtypes
    use clocks
    use gioH5
    use gridutils
    use files
    use kaiomp
    
    implicit none

    character(len=strLen), private :: zcsClk = "Gamera" !Clock ID to use for ZCS calculation
    character(len=strLen), private :: zcsTot = "Omega" !Clock ID to use for ZCS calculation

    real(rp), private :: voltWait = 0.0
    integer, private :: lastTs = -1

    !ConOut_T
    !Console output function pointer
    abstract interface
        subroutine ConsoleOut_T(Model,Grid,State)
            Import :: Model_T, Grid_T, State_T
            type(Model_T), intent(inout) :: Model
            type(Grid_T), intent(in) :: Grid
            type(State_T), intent(in) :: State
        end subroutine ConsoleOut_T
    end interface

    !Set console output function here
    procedure(ConsoleOut_T), pointer :: consoleOutput => consoleOutput_STD

    !tStr_T
    !Time string function
    abstract interface
        subroutine tStr_T(T,tStr)
            Import :: rp,strLen
            real(rp), intent(in) :: T
            character(len=strLen), intent(out) :: tStr
        end subroutine tStr_T
    end interface

    procedure(tStr_T), pointer :: timeString => tStr_STD

contains

    subroutine consoleOutput_STD(Model, Grid, State)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State

        real(rp) :: ZCs_gam, ZCs_tot, wTime_gam,wTime_tot
        integer :: nTh
        character(len=strLen) :: tStr

        ! grab the first reasonable dt after the sim has been running for a bit as dt0
        if (Model%ts > 0 .and. Model%dt0 < TINY) then
            Model%dt0 = Model%dt
        endif

        wTime_gam = readClock(zcsClk)
        wTime_tot = readClock(1)

        !Calculate zone-cycles per second
        Model%kzcsMHD = 0.0
        Model%kzcsTOT = 0.0

        if (Model%ts > 0 .and. lastTs > 0) then
            ZCs_gam = (Model%ts-lastTs)*Grid%Nip*Grid%Njp*Grid%Nkp/wTime_gam
            ZCs_tot = (Model%ts-lastTs)*Grid%Nip*Grid%Njp*Grid%Nkp/wTime_tot

            voltWait = 0.8*voltWait + 0.2*(readClock('VoltSync'))/(readClock(1)+TINY) ! Weighted average to self-correct
        else
            ZCs_gam = 0.0
            ZCs_tot = 0.0
            voltWait = 0
        endif
        lastTs = Model%ts
        
        Model%kzcsMHD = ZCs_gam*1.0e-3
        Model%kzcsTOT = ZCs_tot*1.0e-3

        if (Model%isLoud) then
            nTh = NumOMP()

            write(*,*) ANSICYAN
            write(*,*) 'GAMERA'
            call timeString(Model%t,tStr)
            write(*,'(a,a)')        '      Time = ', trim(tStr)
            write (*, '(a,I8)')     '        ts = ', Model%ts
            call timeString(Model%dt,tStr)
            write (*, '(a,a)')      '        dt = ', trim(tStr)

            if (Model%dt0 > TINY) then
                write (*, '(a,f8.3,a)')      '    dt/dt0 = ', 100*Model%dt/Model%dt0, '%'     
            endif
            if (.not. isnan(voltWait)) then
                write (*, '(a,1f7.1,a)' )    '    Spent ', voltWait*100.0, '% of time waiting for Voltron'
            endif
            if (ZCs_gam>TINY) then
                write (*, '(a,f9.2,a,f9.2,a,I0,a)') '      kZCs = ', Model%kzcsMHD, ' / ', Model%kzcsTOT, ' [MHD/TOT] (',nTh,' threads)'
            endif
            write(*,'(a)',advance="no") ANSIRESET!, ''
        endif

        !Setup for next output
        Model%IO%tCon = Model%IO%tCon + Model%IO%dtCon
        
    end subroutine consoleOutput_STD

    subroutine tStr_STD(T,tStr)
        real(rp), intent(in) :: T
        character(len=strLen), intent(out) :: tStr

        if (T>1.0e-2) then
            write(tStr,'(f0.3,a)' ) T, ' [code]'
        else
            write(tStr,'(es9.2,a)') T, ' [code]'
        endif
    end subroutine tStr_STD

    subroutine fOutput(Model, Grid, State)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State

        character(len=strLen) :: gStr,tStr

        write (gStr, '(A,I0)') "Step#", Model%IO%nOut

        if (Model%isLoud) then
            call timeString(Model%t,tStr)
            write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing HDF5 DATA @ t = ', trim(tStr), ' >', ANSIRESET
        endif

        call writeSlc(Model, Grid, State, gStr)

        !Setup for next output
        Model%IO%tOut = Model%IO%tOut + Model%IO%dtOut
        Model%IO%nOut = Model%IO%nOut + 1
    end subroutine fOutput

    subroutine resOutput(Model,Grid,oState,State)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: oState,State

        character(len=strLen) :: ResF, tStr,lnResF !Name of restart file
        logical :: fExist

        write (ResF, '(A,A,I0.5,A)') trim(Model%RunID), ".gam.Res.", Model%IO%nRes, ".h5"

        call CheckAndKill(ResF)

        if (Model%isLoud) then
            call timeString(Model%t,tStr)
            write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing HDF5 RESTART @ t = ', trim(tStr), ' >', ANSIRESET
        endif

        call writeH5Res(Model,Grid,oState,State,ResF)

        !Setup for next restart
        Model%IO%tRes = Model%IO%tRes + Model%IO%dtRes
        Model%IO%nRes = Model%IO%nRes + 1

        write (lnResF, '(A,A,A,A)') trim(Model%RunID), ".gam.Res.", "XXXXX", ".h5"

        call MapSymLink(ResF,lnResF)
    end subroutine resOutput

end module output
