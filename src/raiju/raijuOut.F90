module raijuOut
    use raijuIO
    use timeHelpers
    use dates

    implicit none

    contains

    subroutine raijuOutput(Model, Grid, State)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        character(len=strLen) :: gStr,tStr

        write (gStr, '(A,I0)') "Step#", State%IO%nOut

        if (Model%isLoud) then
            call timeStrFmt(State%t, tStr)
            write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing RAIJU HDF5 DATA @ t = ', trim(tStr), ' >', ANSIRESET
        endif
        call WriteRAIJU(Model, Grid, State, gStr, Model%writeGhosts)

        !Setup for next output
        State%IO%nOut = State%IO%nOut + 1
        State%IO%tOut = State%t + State%IO%dtOut
    end subroutine raijuOutput


    subroutine raijuResOutput(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        character(len=strLen) :: ResF, tStr,lnResF !Name of restart file
        logical :: fExist

        write (ResF, '(A,A,I0.5,A)') trim(Model%RunID), ".raiju.Res.", State%IO%nRes, ".h5"

        call CheckAndKill(ResF)

        if (Model%isLoud) then
            call timeStrFmt(State%t, tStr)
            write (*, '(a,a,a,a,a)') ANSIGREEN, '<Writing RAIJU HDF5 RESTART @ t = ', trim(tStr), ' >', ANSIRESET
        endif

        call WriteRaijuRes(Model, Grid, State, ResF)

        ! Prep for next restart
        State%IO%tRes = State%IO%tRes + State%IO%dtRes
        State%IO%nRes = State%IO%nRes + 1

        write (lnResF, '(A,A,A,A)') trim(Model%RunID), ".raiju.Res.", "XXXXX", ".h5"

        call MapSymLink(ResF,lnResF)

    end subroutine raijuResOutput


    subroutine raijuResInput(Model, Grid, State)
        !! Mirror of raijuResOutput to make sure state is restored to same exact state as above
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        character(len=strLen) :: tStr

        if (Model%isLoud) then
            call timeStrFmt(State%t, tStr)
            write (*, '(a,a,I0,a,a)') ANSIGREEN, '<Reading RAIJU HDF5 RESTART @ nRes = ', Model%nResIn, ' >', ANSIRESET
            !write (*, '(a,a,a)') ANSIGREEN, '<Reading RAIJU HDF5 RESTART >', ANSIRESET
        endif

        call ReadRaijuResState(Model, Grid, State, Model%ResF)

        ! Prep for next restart
        State%IO%tRes = State%IO%tRes + State%IO%dtRes
        State%IO%nRes = State%IO%nRes + 1

    end subroutine raijuResInput


    subroutine raijuConsoleOut(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        character(len=strLen) :: utStr, tStr, tStr2
        integer :: minDtLoc, maxDtLoc
        real(rp) :: minDt, maxDt
        integer :: s, sIdx
        integer, dimension(2) :: mpLoc
        real(rp) :: maxPress, maxDen, maxP_Xmin, maxP_Ymin, maxP_L, maxP_MLT

        call mjd2utstr(State%mjd,utStr)
        minDtLoc = minloc(State%dtk,dim=1)
        maxDtLoc = maxloc(State%dtk,dim=1)


        write(*,*) ANSIPURPLE
        write(*,*) 'RAIJU'
        write(*,'(a,a)')        '      UT   = ', trim(utStr)
        call timeStrFmt(State%t, tStr)
        write(*,'(a,a)')        '      Time = ', trim(tStr)
        call timeStrFmt(State%dt, tStr)
        write(*,'(a,a)')        '     dtCpl = ', trim(tStr)
        call timeStrFmt(State%dtk(maxDtLoc), tStr )
        call timeStrFmt(State%dtk(minDtLoc), tStr2)
        write(*,'(a)'  )        '     Max/Min dt @ k:'
        write(*,'(a,a,a,I0)') '        Max', trim(tStr), ' @ ', maxDtLoc
        write(*,'(a,a,a,I0)') '        Min', trim(tStr2) , ' @ ', minDtLoc
        write(*,'(a)'  )        '     Flav : max Press/Den @ L,MLT; DPS-dst'
        do s=1, Model%nSpc
            sIdx = spcIdx(Grid, Grid%spc(s)%flav)
            if (Grid%spc(s)%flav == F_PSPH) then
                mPLoc = maxloc(State%Den(sIdx)%data)  ! First element = 1, need to account for index offset after
            else
                mPLoc = maxloc(State%Press(sIdx)%data)  ! First element = 1, need to account for index offset after
            endif
            mpLoc(1) = mpLoc(1) + Grid%shGrid%isg - 1
            mpLoc(2) = mpLoc(2) + Grid%shGrid%jsg - 1

            maxPress = State%Press(sIdx)%data(mPLoc(1), mpLoc(2))
            maxDen   = State%Den  (sIdx)%data(mPLoc(1), mpLoc(2))!/Grid%spc(sIdx)%amu
            maxP_Xmin = State%xyzMincc(mPLoc(1), mpLoc(2),XDIR)
            maxP_Ymin = State%xyzMincc(mPLoc(1), mpLoc(2),YDIR)

            maxP_L = sqrt(maxP_Xmin**2 + maxP_Ymin**2)
            maxP_MLT = atan2(maxP_Ymin, maxP_Xmin)/PI*12D0 + 12D0
            if (maxP_MLT > 24) maxP_MLT = maxP_MLT - 24D0
            write(*,'(a,I0,a,f7.2,a,f7.2,a,f5.2,a,f5.2,a,f7.2)') '        ', &
                Grid%spc(s)%flav, ': P =', maxPress,', D =',maxDen,' @ ',maxP_L,' Rp,',maxP_MLT, &
                " MLT; DPS:",spcEta2DPS(Model, Grid, State%bvol_cc, State%eta_avg, Grid%spc(sIdx), State%active .eq. RAIJUACTIVE)

        enddo
        write(*,'(a)',advance="no") ANSIRESET

        State%IO%tCon = State%IO%tCon + State%IO%dtCon

    end subroutine raijuConsoleOut

!------
! Helpers
!------

    subroutine genResInFname(Model, ResF, runIdO)
        !!! Using Model members, determine the restart name to read from
        type(raijuModel_T), intent(in) :: Model
        character(len=strLen), intent(out) :: ResF
        character(len=*), optional, intent(in) :: runIdO

        character(len=strLen) :: runId
        character(len=strLen) :: nStr

        if (present(runIdO)) then
            runId = trim(runIdO)
        else
            runId = Model%RunID
        endif

        if (Model%nResIn == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') Model%nResIn
        endif

        write (ResF, '(A,A,A,A)') trim(runId), ".raiju.Res.", trim(nStr), ".h5"
    end subroutine genResInFname

end module raijuOut