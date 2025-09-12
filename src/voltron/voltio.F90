!Routines to write restart/outputs from Voltron
module voltio
    use gamapp
    use volttypes
    use cmiutils
    use mixio
    use clocks
    use wind
    use dyncoupling
    use dstutils
    use planethelper
    use shellGridIO
    use voltappHelper

    use tubehelper, only : writeTubeShellRestart, readTubeShellRestart
    
    implicit none

    ! these will be assigned values during voltron console output so that the output
    ! will be 90% of the correct value after about 1 model minute with the fast weight, and
    ! 90% of the correct value after about 5 model minutes with the slow weight.
    ! these are just default values
    real(rp), private :: fastWeight = 0.8
    real(rp), private :: slowWeight = 0.95

    integer , parameter, private :: MAXVOLTIOVAR = 50
    real(rp), parameter, private :: dtWallMax = 1.0 !How long between timer resets[hr]
    logical , private :: isConInit = .false.
    real(rp), private ::  oMJD = 0.0
    integer , private :: oTime = 0.0
    real(rp), private :: gamWait = 0.0
    real(rp), private :: mixWait = 0.0
    real(rp), private :: imagWait = 0.0
    real(rp), private :: chimpWait = 0.0
    real(rp), private :: tubesWait = 0.0
    real(rp), private :: ioWait = 0.0
    real(rp), private :: simRate = 0.0
    character(len=strLen), private :: vh5File

    contains

    subroutine consoleOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Using console output from Gamera
        call gApp%WriteConsoleOutput()

        !Using console output from Voltron
        call consoleOutputVOnly(vApp,gApp,gApp%Model%MJD0)

    end subroutine consoleOutputV

    subroutine consoleOutputVOnly(vApp,gApp,MJD0)
        class(voltApp_T), intent(inout) :: vApp
        class(gamApp_T) , intent(in) :: gApp
        real(rp), intent(in) :: MJD0

        real(rp) :: cpcp(2) = 0.0

        real(rp) :: dpT,dtWall,cMJD,dMJD

        integer :: nTh,curCount,countMax
        real(rp) :: clockRate
        character(len=strLen) :: utStr
        real(rp) :: BSDst0,AvgBSDst,DPSDst,symh,BSSMRs(4)

        dpT = vApp%tilt%evalAt(vApp%time)*180.0/PI

        !Figure out some perfromance info
        cMJD = T2MJD(vApp%time,MJD0) !Current MJD

        if (isConInit) then
            !Console output has been initialized
            dMJD = cMJD - oMJD !Elapsed MJD since first recorded value
            oMJD = cMJD ! clock every output separately
            call system_clock(curCount,clockRate,countMax)
            dtWall = (curCount - oTime)/clockRate
            oTime = curCount
            if(dtWall < 0) dtWall = dtWall + countMax / clockRate
            simRate = slowWeight*simRate + (1.0-slowWeight)*dMJD*24.0*60.0*60.0/dtWall !Model seconds per wall second
            ! Use weight average to self-correct model timing
            ! chimp timing is a little confusing, but it combines local time spent squishing (no helpers) with
            !   actual helper delay, which is tricky to estimate
            gamWait   = fastWeight*gamWait   + (1.0-fastWeight)*readClock('GameraSync')/(readClock(1)+TINY)
            chimpWait = fastWeight*chimpWait + (1.0-fastWeight)*(readClock('Squish')+readClock('VoltHelpers'))/(readClock(1)+TINY)
            imagWait  = fastWeight*imagWait  + (1.0-fastWeight)*readClock('InnerMag')/(readClock(1)+TINY)
            mixWait   = fastWeight*mixWait   + (1.0-fastWeight)*readClock('ReMIX')/(readClock(1)+TINY)
            tubesWait = fastWeight*tubesWait + (1.0-fastWeight)*readClock('VoltTubes')/(readClock(1)+TINY)
            ioWait    = fastWeight*ioWait    + (1.0-fastWeight)*readClock('IO')/(readClock(1)+TINY)
        else
            ! calculate weights for simple exponential smoothing
            ! halve the time constant to approximate 90% of correct value in formula
            ! fast is 90% correct after 1 model minute
            fastWeight = exp(-1.0*vApp%IO%dtCon/30.0_rp)
            ! slow is 90% correct after 5 model minutes
            slowWeight = exp(-1.0*vApp%IO%dtCon/150.0_rp)
            simRate = 0.0
            oMJD = cMJD
            call system_clock(count=oTime)
            isConInit = .true.
            dtWall = 0.0
            gamWait = 0.0
            mixWait = 0.0
            imagWait = 0.0
            chimpWait = 0.0
            tubesWait = 0.0
            ioWait = 0.0
        endif

        !Get MJD info
        call mjd2utstr(cMJD,utStr)

        !Get Dst estimate: DPS, center of earth, MLT avg of equatorial stations
        call EstDST(gApp%Model,gApp%Grid,gApp%State,BSDst0=BSDst0,AvgBSDst=AvgBSDst,BSSMRs=BSSMRs,DPSDst=DPSDst)

        vApp%BSDst = AvgBSDst
        
        !Get symh from input time series
        symh = vApp%symh%evalAt(vApp%time)

        if (vApp%isLoud) then
            write(*,*) ANSIBLUE
            write(*,*) 'VOLTRON'
            write (*,'(a,a)')                    '      UT   = ', trim(utStr)
            write (*, '(a,1f8.3,a)')             '      tilt = ' , dpT, ' [deg]'
            write (*, '(a, f8.3,a)')             '    Sym-H  = ' , symh  , ' [nT]'
            write (*, '(a, f8.3,a)')             '    BSDst  ~ ' , AvgBSDst , ' [nT]'
            write (*, '(a,4f8.2,a)')             '           dSMRs  ~ ' , BSSMRs-AvgBSDst, ' [nT, 12/18/00/06]'
            !write (*, '(a,4f8.2,a)')             '   BSSMRs  ~ ' , BSSMRs, ' [nT, 12/18/00/06]'

            if (vApp%doDeep .and. (vApp%time>0.0)) then
                write (*, '(a, f8.3,a)')             '   DPSDst  ~ ' , DPSDst, ' [nT]'
                !write (*,'(a,1f8.3,I6,a)')           '      xTrc = ', vApp%rTrc,vApp%nTrc, ' [r/n]'
            endif
            write (*, '(a,1f7.1,a)' ) '   Spent ', gamWait*100.0,   '% of time waiting for Gamera'
            write (*, '(a,1f7.1,a)' ) '         ', chimpWait*100.0, '% of time processing Chimp(Helpers)'
            write (*, '(a,1f7.1,a)' ) '         ', imagWait*100.0,  '% of time processing IMAG'
            write (*, '(a,1f7.1,a)' ) '         ', mixWait*100.0,   '% of time processing Remix'
            write (*, '(a,1f7.1,a)' ) '         ', tubesWait*100.0, '% of time processing Tubes(Helpers)'
            write (*, '(a,1f7.1,a)' ) '         ', ioWait*100.0,    '% of time writing files'
            if (simRate>TINY) then
                nTh = NumOMP()
                write (*, '(a,1f8.3,a,I0,a)')             '    Running @ ', simRate*100.0, '% of real-time (',nTh,' threads)'  
            endif
            
            write(*,'(a)',advance="no") ANSIRESET!, ''
        endif

        !Write inner mag console IO if needed
        if (vApp%doDeep) then
            !call vApp%imagApp%doConIO(vApp%MJD,vApp%time)
            call vApp%imagApp%WriteConsoleOutput()
        endif

        !Setup for next output
        vApp%IO%tCon = vApp%IO%tCon + vApp%IO%dtCon
        
        if (vApp%doDynCplDT) then
            call UpdateCouplingCadence(vApp)
        endif

    end subroutine consoleOutputVOnly

    !Given vector, get clock/cone angle and magnitude
    function ClockConeMag(V) result(aVec)
        real(rp), dimension(NDIM), intent(in) :: V
        real(rp), dimension(NDIM) :: aVec

        real(rp) :: MagV
        MagV = norm2(V)
        aVec(2) = atan2(V(YDIR),V(ZDIR) )*180.0/PI !Clock angle
        if (aVec(2) < 0) then
            aVec(2) = aVec(2) + 360.0
        endif
        
        if (MagV>TINY) then
            aVec(3) = acos (V(XDIR)/MagV)*180.0/PI
        else
            aVec(3) = 0.0
        endif
        aVec(1) = MagV
    end function ClockConeMag

    subroutine resOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Write Gamera restart
        call gApp%WriteRestart(vApp%IO%nRes)

        !Write Voltron restart data
        call resOutputVOnly(vApp,gApp)

    end subroutine resOutputV

    subroutine resOutputVOnly(vApp, gApp)
        class(voltApp_T), intent(inout) :: vApp
        class(gamApp_T) , intent(inout) :: gApp

        if (vApp%writeFiles) then
            call writeMIXRestart(vApp%remixApp%ion,vApp%IO%nRes,mjd=vApp%MJD,time=vApp%time)
            !Write inner mag restart
            if (vApp%doDeep) then
                !call vApp%imagApp%doRestart(vApp%IO%nRes,vApp%MJD,vApp%time)
                call vApp%imagApp%WriteRestart(vApp%IO%nRes)
            endif
            call writeVoltRestart(vApp,gApp)
        endif

        vApp%IO%tRes = vApp%IO%tRes + vApp%IO%dtRes
        vApp%IO%nRes = vApp%IO%nRes + 1

    end subroutine resOutputVOnly

    subroutine writeVoltRestart(vApp,gApp)
        class(voltApp_T), intent(in) :: vApp
        class(gamApp_T) , intent(in) :: gApp

        character(len=strLen) :: ResF,lnResF
        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars

        write (ResF, '(A,A,I0.5,A)') trim(gApp%Model%RunID), ".volt.Res.", vApp%IO%nRes, ".h5"
        call CheckAndKill(ResF)

        call ClearIO(IOVars)

        !Main attributes
        call AddOutVar(IOVars,"nOut",vApp%IO%nOut)
        call AddOutVar(IOVars,"nRes",vApp%IO%nRes)
        call AddOutVar(IOVars,"ts"  ,vApp%ts)
        call AddOutVar(IOVars,"MJD" ,vApp%MJD)
        call AddOutVar(IOVars,"time",vApp%time)

        !Coupling info
        call AddOutVar(IOVars,"CoupleT", vApp%DeepT)
        call AddOutVar(IOVars,"gBAvg",   vApp%mhd2Mix%gBAvg)
        
        ! State object
        call AddOutSGV(IOVars, "potential_total", vApp%State%potential_total, doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "potential_corot", vApp%State%potential_corot, doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "bIonoMag", vApp%State%bIonoMag, doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "bIonoRad", vApp%State%bIonoRad, doWriteMaskO=.true.)

        call WriteVars(IOVars,.false.,ResF)

        ! Save shellGrid
        call writeShellGrid(vApp%shGrid, ResF)

        ! And all TubeShell vars
        call writeTubeShellRestart(vApp%State%tubeShell, ResF)

        !Create link to latest restart
        write (lnResF, '(A,A,A,A)') trim(gApp%Model%RunID), ".volt.Res.", "XXXXX", ".h5"
        call MapSymLink(ResF,lnResF)

    end subroutine writeVoltRestart

    subroutine readVoltronRestart(vApp,resID, nRes)
        class(voltApp_T), intent(inout) :: vApp
        character(len=*), intent(in) :: resID
        integer, intent(in) :: nRes

        character(len=strLen) :: ResF,nStr
        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars
        logical :: fExist
        integer :: n0

        !Get number string
        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') nRes
        endif

        write (ResF, '(A,A,A,A)') trim(resID), ".volt.Res.", trim(nStr), ".h5"
        write(*,*) 'Reading Voltron restart from ', trim(ResF)
        inquire(file=ResF,exist=fExist)
        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open input voltron restart file, exiting'
            stop
        endif

        ! Re-init our shellGrid from file
        call GenShellGridFromFile(vApp%shGrid, "VOLTRON", ResF)
        ! Init our state now that we have grid info back
        call initVoltState(vApp)

        !! TODO: Load new voltState vars from restart file

        call ClearIO(IOVars)
        call AddInVar(IOVars,"gBAvg")

        call AddInVar(IOVars,"nOut"    ,vTypeO=IOINT)
        call AddInVar(IOVars,"nRes"    ,vTypeO=IOINT)
        call AddInVar(IOVars,"ts"      ,vTypeO=IOINT)
        call AddInVar(IOVars,"MJD"     ,vTypeO=IOREAL)
        call AddInVar(IOVars,"time"    ,vTypeO=IOREAL)
        call AddInVar(IOVars,"CoupleT" ,vTypeO=IOREAL)


        !Get data
        call ReadVars(IOVars,.false.,ResF)

        !Check to see if CoupleT is present
        n0 = FindIO(IOVars,"CoupleT")
        if (.not. IOVars(n0)%isDone) then
            write(*,*) "CoupleT not found in Voltron restart."
            write(*,*) "This restart must have been written with an older"
            write(*,*) " version of the code that used ShallowT and DeepT."
            write(*,*) "This code is not compatible with this restart file."
            write(*,*) "Please regenerate the restart file."
            stop
        endif

        vApp%IO%nOut  = GetIOInt(IOVars,"nOut")
        vApp%IO%nRes  = GetIOInt(IOVars,"nRes") + 1
        vApp%ts       = GetIOInt(IOVars,"ts")
        vApp%MJD      = GetIOReal(IOVars,"MJD")
        vApp%time     = GetIOReal(IOVars,"time")
        vApp%DeepT    = GetIOReal(IOVars,"CoupleT")

        !Check to see if gB0 is present
        n0 = FindIO(IOVars,"gBAvg")
        if (IOVars(n0)%isDone) then
            if (.not. allocated(vApp%mhd2Mix%gBAvg)) then
                !Allocate this if necessary
                allocate( vApp%mhd2Mix%gBAvg(IOVars(n0)%dims(1),IOVars(n0)%dims(2),IOVars(n0)%dims(3),IOVars(n0)%dims(4)) )
            endif
            call IOArray4DFill(IOVars,"gBAvg",vApp%mhd2Mix%gBAvg)
        else
            write(*,*) "gBAvg not found in Voltron restart, assuming dipole ..."
        endif

        call ReadInSGV(vApp%State%potential_total, ResF, "potential_total", doIOpO=.false.)
        call ReadInSGV(vApp%State%potential_corot, ResF, "potential_corot", doIOpO=.false.)
        call ReadInSGV(vApp%State%bIonoMag, ResF, "bIonoMag", doIOpO=.false.)
        call ReadInSGV(vApp%State%bIonoRad, ResF, "bIonoRad", doIOpO=.false.)

        ! And all of tubeShell
        call readTubeShellRestart(vApp%State%tubeShell, ResF)

    end subroutine readVoltronRestart


    subroutine fOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Write gamera data
        call gApp%WriteFileOutput(vApp%IO%nOut)

        !Write voltron data
        call fOutputVOnly(vApp,gApp)

    end subroutine fOutputV

    subroutine fOutputVOnly(vApp,gApp)
        class(voltApp_T), intent(inout) :: vApp
        class(gamApp_T) , intent(inout) :: gApp
        
        if(vApp%writeFiles) then
            !Write ReMIX data
            call writeMix(vApp%remixApp%ion,vApp%IO%nOut,mjd=vApp%MJD,time=vApp%time)

            !Write inner mag IO if needed
            if (vApp%doDeep) then
                !call vApp%imagApp%doIO(vApp%IO%nOut,vApp%MJD,vApp%time)
                call vApp%imagApp%WriteFileOutput(vApp%IO%nOut)
            endif

            call WriteVolt(vApp,gApp,vApp%IO%nOut)
        endif

        vApp%IO%tOut = vApp%IO%tOut + vApp%IO%dtOut
        vApp%IO%nOut = vApp%IO%nOut + 1

    end subroutine fOutputVOnly

    !Output voltron data
    subroutine WriteVolt(vApp,gApp,nOut,doGhostsO)
        class(voltApp_T), intent(inout) :: vApp
        type(gamApp_T)  , intent(inout) :: gApp
        integer, intent(in) :: nOut
        logical, intent(in), optional :: doGhostsO

        character(len=strLen) :: gStr
        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars
        real(rp) :: symh

        integer :: is,ie,js,je,nClkSteps
        real(rp) :: Csijk,Con(NVAR)
        real(rp) :: BSDst0,AvgBSDst,DPSDst,BSSMRs(4)
        integer, dimension(4) :: outSGVBnds_corner
        logical :: doGhosts

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        if (doGhosts) then
            is = vApp%shGrid%isg
            ie = vApp%shGrid%ieg
            js = vApp%shGrid%jsg
            je = vApp%shGrid%jeg
        else
            is = vApp%shGrid%is
            ie = vApp%shGrid%ie
            js = vApp%shGrid%js
            je = vApp%shGrid%je
        endif

        !Get symh from input time series
        symh = vApp%symh%evalAt(vApp%time)

        outSGVBnds_corner = (/is,ie+1,js,je+1/)

        !Get Dst estimate: DPS, center of earth, MLT avg of equatorial stations
        call EstDST(gApp%Model,gApp%Grid,gApp%State,BSDst0,AvgBSDst,BSSMRs,DPSDst)
        vApp%BSDst = AvgBSDst

        write(gStr,'(A,I0)') "Step#", nOut

        !Reset IO chain
        call ClearIO(IOVars)

        !---------------------
        !Do attributes

        call AddOutVar(IOVars,"BSDst" ,AvgBSDst)
        call AddOutVar(IOVars,"BSDst0",BSDst0)
        call AddOutVar(IOVars,"DPSDst",DPSDst)
        call AddOutVar(IOVars,"BSSMR12",BSSMRs(1))
        call AddOutVar(IOVars,"BSSMR18",BSSMRs(2))
        call AddOutVar(IOVars,"BSSMR00",BSSMRs(3))
        call AddOutVar(IOVars,"BSSMR06",BSSMRs(4))

        call AddOutVar(IOVars,"SymH",symh)

        call AddOutVar(IOVars,"time" ,vApp%time)
        call AddOutVar(IOVars,"MJD"  ,vApp%MJD)
        call AddOutVar(IOVars,"timestep",vApp%ts)

        !Performance metrics
        nClkSteps = readNCalls('DeepUpdate')
        call AddOutVar(IOVars,"_perf_stepTime",readClock(1)/nClkSteps)
        call AddOutVar(IOVars,"_perf_deepUpdateTime",readClock(1)/nClkSteps)
        call AddOutVar(IOVars,"_perf_gamTime", readClock('GameraSync')/nClkSteps)
        call AddOutVar(IOVars,"_perf_squishTime", (readClock('Squish')+readClock('VoltHelpers'))/nClkSteps)
        call AddOutVar(IOVars,"_perf_imagTime", readClock('InnerMag')/nClkSteps)
        call AddOutVar(IOVars,"_perf_mixTime", readClock('ReMIX')/nClkSteps)
        call AddOutVar(IOVars,"_perf_tubesTime", readClock('VoltTubes')/nClkSteps)
        call AddOutVar(IOVars,"_perf_ioTime", readClock('IO')/nClkSteps)

        ! voltState stuff
        call AddOutSGV(IOVars, "Potential_total", vApp%State%potential_total, &
                        uStr="kV", dStr="Ionospheric electrostatic potential (ExB + corotation)", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Potential_corot", vApp%State%potential_corot, &
                        uStr="kV", dStr="Corotation potential", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Potential_drop", vApp%State%pot_drop, &
                        uStr='"kV', dStr="Potential drop between given location and its conjugate footpoint", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.true.)
        associate(tubeShell=>vApp%State%tubeShell)
        call AddOutSGV(IOVars, "xMin", tubeShell%X_bmin(XDIR), &
                        uStr="Rx", dStr="SM-X location of flux tube bMin", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "yMin", tubeShell%X_bmin(YDIR), &
                        uStr="Rx", dStr="SM-Y location of flux tube bMin", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "zMin", tubeShell%X_bmin(ZDIR), &
                        uStr="Rx", dStr="SM-Z location of flux tube bMin", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "bMin", tubeShell%bmin, &
                        uStr="nT", dStr="Field strength at magnetic equator", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "topo", tubeShell%topo, &
                        dStr="Magnetic field topology (0=closed,1=open,2=undefined)", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "bVol", tubeShell%bVol, &
                        uStr="Rx/nT", dStr="Flux tube volume (if closed)", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "lat0", tubeShell%lat0, &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "lon0", tubeShell%lon0, &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "x0", tubeShell%xyz0(XDIR), &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "y0", tubeShell%xyz0(YDIR), &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "z0", tubeShell%xyz0(ZDIR), &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "latc", tubeShell%latc, &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "lonc", tubeShell%lonc, &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "rCuravture", tubeShell%rCurv, &
                        uStr="Rx", dStr="Radius of curvature @ bmin if this is a closed field line", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "tube_nTrc", tubeShell%nTrc, &
                        dStr="Number of points along field line", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "tube_wIMAG", tubeShell%wMAG, &
                        dStr="'Good' weighting for IMAG model based on Alfven speed vs. magnetofast and ExB flow", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "tube_avgP", tubeShell%avgP, &
                        ustr="nPa",dStr="Flux tube averaged pressure", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "tube_avgD", tubeShell%avgN, &
                        uStr="#/cc",dStr="Flux tube averaged density", &
                        outBndsO=outSGVBnds_corner, doWriteMaskO=.false.)
        end associate

        call WriteVars(IOVars,.true.,vh5File,gStr)


    end subroutine WriteVolt


    !Initialize Voltron-unique IO
    subroutine InitVoltIO(vApp,gApp,doGhostsO)
        class(voltApp_T), intent(inout) :: vApp
        type(gamApp_T)  , intent(inout) :: gApp
        logical, intent(in), optional :: doGhostsO

        character(len=strLen) :: RunID
        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars
        logical :: fExist, isRestart
        integer :: i,j,is,ie,js,je
        logical :: doGhosts

        real(rp), dimension(:,:), allocatable :: colat2D, lon2D
        real(rp), dimension(:,:), allocatable :: X2D, Y2D, Z2D
        real(rp), dimension(:,:), allocatable :: dLat, areaCC, bMag, bRad
        real(rp), dimension(:), allocatable :: cosThc

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        if (doGhosts) then
            is = vApp%shGrid%isg
            ie = vApp%shGrid%ieg
            js = vApp%shGrid%jsg
            je = vApp%shGrid%jeg
        else
            is = vApp%shGrid%is
            ie = vApp%shGrid%ie
            js = vApp%shGrid%js
            je = vApp%shGrid%je
        endif

        isRestart = gApp%Model%isRestart
        RunID = trim(gApp%Model%RunID)

        !Create filename
        vh5File = trim(RunID) // ".volt.h5" !Voltron output
        fExist = CheckFile(vh5File)
        write(*,*) 'Voltron outputting to ',trim(vh5File)

        if ( (.not. isRestart) .or. (isRestart .and. (.not. fExist)) ) then
            !Not a restart or it is a restart and no file
            call CheckAndKill(vh5File) !For non-restart but file exists

            !Reset IO chain
            call ClearIO(IOVars)

            ! Add shell grid as spatial array
            allocate(colat2D(is:ie+1,js:je+1))  ! +1 because we're doing corners
            allocate(lon2D(is:ie+1,js:je+1))
            allocate(X2D(is:ie+1,js:je+1))
            allocate(Y2D(is:ie+1,js:je+1))
            allocate(Z2D(is:ie+1,js:je+1))
            do i=js,je+1
                colat2D(:,i) = vApp%shGrid%th(is:ie+1)
            enddo
            do i=is,ie+1
                lon2D(i,:) = vApp%shGrid%ph(js:je+1)
            enddo
            !call AddOutVar(IOVars,"X",colat2D,uStr="radians")
            !call AddOutVar(IOVars,"Y",  lon2D,uStr="radians")

            X2D = vApp%shGrid%radius*sin(colat2D)*cos(lon2D)
            Y2D = vApp%shGrid%radius*sin(colat2D)*sin(lon2D)
            Z2D = vApp%shGrid%radius*cos(colat2D)
            call AddOutVar(IOVars,"X",X2D,uStr="Rp")
            call AddOutVar(IOVars,"Y",Y2D,uStr="Rp")
            call AddOutVar(IOVars,"Z",Z2D,uStr="Rp")


            ! Some derived cell-centered stuff
            associate(shGr=>vApp%shGrid)

            allocate(dLat(is:ie,js:je))
            allocate(areaCC(is:ie,js:je))
            allocate(cosThc(is:ie))
            allocate(bMag(is:ie,js:je))
            allocate(bRad(is:ie,js:je))

            do i=is,ie
                dLat(i,:) = vApp%shGrid%th(i+1) - vApp%shGrid%th(i)
            enddo
            cosThc = cos(shGr%thc)
            do i=is,ie
                do j=js,je
                    ! r^2 * sin(th) * dTh * dPh
                    areaCC(i,j) = (shGr%radius)**2 &
                                    * sin(shGr%thc(i)) &
                                    * (shGr%th(i+1) - shGr%th(i)) &
                                    * (shGr%ph(j+1) - shGr%ph(j))
                    bMag(:,j) = vApp%planet%magMoment*G2nT &
                                /(shGr%radius)**3.0 &
                                * sqrt(1.0+3.0*cosThc**2.0)  ! [nT]
                    bRad(:,j) = -1*vApp%planet%magMoment*G2nT &
                                /(shGr%radius)**3.0 &
                                * 2*cosThc  ! [nT]
                enddo
            enddo
            call AddOutVar(IOVars,'dLat',dLat,uStr="rad")
            call AddOutVar(IOVars,'areaCC',areaCC,uStr="Rp^2")
            call AddOutVar(IOVars,'bMag',bMag,uStr="nT",dStr="Total magnetic field strength at cell center")
            call AddOutVar(IOVars,'bRad',bRad,uStr="nT",dStr="Radial component of magnetic field strength at cell center")
            
            end associate

            call AddOutVar(IOVars,"UnitsID","VOLTRON")
            call WriteVars(IOVars,.true.,vh5File)

            call writeShellGrid(vApp%shGrid, vh5File, "/ShellGrid")
        endif

    end subroutine InitVoltIO

end module voltio

