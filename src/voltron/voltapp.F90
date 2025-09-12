! Collection of data and objects for the voltron middle man

module voltapp
    use mixtypes
    use mixdefs
    use mixmain
    use ebtypes
    use chmpdefs
    use starter
    use mhd2mix_interface
    use mhd2chmp_interface
    use chmp2mhd_interface
    use mix2mhd_interface
    use imag2mhd_interface
    use imag2mix_interface
    use ebsquish
    use innermagsphere
    use dates
    use kronos
    use voltio
    use msphutils, only : RadIonosphere
    use gcminterp
    use gcmtypes
    use planethelper
    use shellGrid
    use shellGridGen
    use voltappHelper
    use voltCplHelper
    use apex, only : init_apex
    
    implicit none

    contains

    !Initialize Voltron (after Gamera has already been initialized)
    subroutine initVoltron(vApp,optFilename)
        class(voltApp_T), intent(inout) :: vApp
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: inpXML, kaijuRoot, resID
        type(XML_Input_T) :: xmlInp
        type(TimeSeries_T) :: tsMJD
        real(rp) :: tSpin,tIO
        logical :: doSpin,isK,doRestart
        integer :: nRes

        if(.not. allocated(vApp%gApp)) then
            ! non-mpi voltron uses non-mpi local coupled gamera
            ! but don't over-ride if someone else allocated first
            allocate(gamCoupler_T :: vApp%gApp)
        endif

        associate(gApp=>vApp%gApp)

        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
            call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")
        else
            !Find input deck
            call getIDeckStr(inpXML)
        endif

        if(vApp%isLoud) write(*,*) 'Voltron Reading input deck from ', trim(inpXML)

        !Create XML reader
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Voltron',.true.)

        !Start by shutting up extra ranks
        if (.not. vApp%isLoud) call xmlInp%BeQuiet()

        ! try to verify that the XML file has "Kaiju" as a root element
        kaijuRoot = ""
        call xmlInp%Get_Key_Val("/Gamera/sim/H5Grid",kaijuRoot, .false.)
        if(len(trim(kaijuRoot)) /= 0) then
            write(*,*) "The input XML appears to be of an old style."
            write(*,*) "As of June 12th, 2021 it needs a root element of <Kaiju>."
            write(*,*) "Please modify your XML config file by adding this line at the top:"
            write(*,*) "<Kaiju>"
            write(*,*) "and this line at the bottom:"
            write(*,*) "</Kaiju>"
            write(*,*) "OR (preferred) convert your configuration to an INI file and use"
            write(*,*) " the XMLGenerator.py script to create conforming XML files."
            write(*,*) "Please refer to the python script or"
            write(*,*) " the [Generating XML Files] wiki page for additional info."
            stop
        endif

        call xmlInp%SetVerbose(.true.)

        !Setup OMP
        call SetOMP(xmlInp)

        !initialize coupled Gamera
        call xmlInp%SetRootStr('Kaiju/Gamera')
        gApp%gOptions%userInitFunc => vApp%vOptions%gamUserInitFunc
        call gApp%InitModel(xmlInp)
        call gApp%InitIO(xmlInp)

        ! adjust XMl reader root
        call xmlInp%SetRootStr('Kaiju/Voltron')

        !Initialize planet information
        call getPlanetParams(vApp%planet, xmlInp)
        if (vApp%isLoud) then
            call printPlanetParams(vApp%planet)
        endif
        call xmlInp%Set_Val(vApp%doGeoCorot,"/Kaiju/Gamera/prob/doGeoCorot",vApp%doGeoCorot)

        !Initialize state information
        !Check for Earth to decide what things need to happen
        if (trim(gApp%Model%gamOut%uID) == "EARTH") then
            vApp%isEarth = .true.
            if (vApp%isLoud) write(*,*) "Going into geospace mode ..."
        else
            vApp%isEarth = .false.
            if (vApp%isLoud) write(*,*) "Not using geospace mode ..."
        endif

        !Set file to read from and pass desired variable name to initTS
        call xmlInp%Set_Val(vApp%tilt%wID,"/Kaiju/Gamera/wind/tsfile","NONE")
        call vApp%tilt%initTS("tilt",doLoudO=.false.)
        vApp%symh%wID = vApp%tilt%wID
        call vApp%symh%initTS("symh",doLoudO=.false.)
        if (vApp%isEarth) then
            !Initialize TM03 model in case we wanna use it
            call InitTM03(vApp%tilt%wID,0.0_rp)
        endif

        !Time options
        call xmlInp%Set_Val(vApp%tFin,'time/tFin',1.0_rp)
        
        call vApp%IO%init(xmlInp,vApp%time,vApp%ts)
        
        !Deep coupling
        if (xmlInp%Exists("coupling/dt") .or. xmlInp%Exists("coupling/dtDeep")) then
                write(*,*) 'Please remove all instances of voltron/coupling/dt and voltron/coupling/dtDeep'
                write(*,*) '   from the input XML. They have been replaced with a single unified input named'
                write(*,*) '   voltron/coupling/dtCouple, which controls all coupling and is set to 5 seconds by default'
                stop
        endif

        call xmlInp%Set_Val(vApp%DeepDT, "coupling/dtCouple", 5.0_rp)
        call xmlInp%Set_Val(vApp%rTrc,   "coupling/rTrc"  , 40.0)

        !Termination can have issues if tFin is too close to a coupling time
        if(MODULO(vApp%tFin,vApp%DeepDT) < 0.1_rp .or. (vApp%DeepDT-MODULO(vApp%tFin,vApp%DeepDT)) < 0.1_rp) then
            write (*,*) "Ending a simulation too close to a coupling interval can cause synchronization issues"
            write (*,*) "Increasing the ending time by a fraction of a second to create a buffer"
            vApp%tFin = vApp%tFin + 0.25_rp
        endif

        !Correct Gamera's end time
        gApp%Model%tFin = vApp%tFin/gApp%Model%Units%gT0

        !Coupling is unified, so adding a separate XML option to control "deep" parts
        call xmlInp%Set_Val(vApp%doDeep, "coupling/doDeep", .true.)

        call xmlInp%Set_Val(vApp%doGCM, "coupling/doGCM",.false.)
        if (vApp%isEarth) then
            call xmlInp%Set_Val(vApp%mhd2mix%dtAvg ,"coupling/dtAvgB0",900.0)
        else
            call xmlInp%Set_Val(vApp%mhd2mix%dtAvg ,"coupling/dtAvgB0",-1.0)
        endif
        !Figure out weighting for exponential moving average (EMA)
        !Want weighting such that ~95% of the weight comes from the last dtAvg seconds
        if (vApp%mhd2mix%dtAvg > 0) then
            vApp%mhd2mix%wAvg = 1.0 - exp(-3*vApp%DeepDT/max(vApp%mhd2mix%dtAvg,vApp%DeepDT))
        else
            vApp%mhd2mix%wAvg = 0.0 !Ignore any corrections after initial dipole value
        endif

        call xmlInp%Set_Val(doRestart,"/Kaiju/gamera/restart/doRes",.false.)
        if(doRestart) then
            call xmlInp%Set_Val(resID,"/Kaiju/gamera/restart/resID","msphere")
            call xmlInp%Set_Val(nRes,"/Kaiju/gamera/restart/nRes" ,-1)
            call readVoltronRestart(vApp, resID, nRes)
            vApp%IO%tOut = floor(vApp%time/vApp%IO%dtOut)*vApp%IO%dtOut + vApp%IO%dtOut
            vApp%IO%tRes = floor(vApp%time/vApp%IO%dtRes)*vApp%IO%dtRes + vApp%IO%dtRes

            vApp%IO%tCon = vApp%time
        else
            ! non-restart initialization
            !Check for spinup info
            call xmlInp%Set_Val(doSpin,"spinup/doSpin",.true.)
            ! Deep enabled, not restart, not spinup is an error. Restart or spinup is required
            if (vApp%doDeep .and. (.not. doSpin) ) then
                write(*,*) 'Spinup is required with deep coupling. Please enable the spinup/doSpin option. At least 1 minute of spinup is recommended.'
                stop
            endif
            if (doSpin) then
                call xmlInp%Set_Val(tSpin,"spinup/tSpin",tSpinDef)
                !Rewind time to negative tSpin (seconds)
                vApp%time =-tSpin
                call xmlInp%Set_Val(tIO,"spinup/tIO",0.0) !Time of first restart and output
                vApp%IO%tRes = tIO
                vApp%IO%tOut = tIO
            endif
            !Use MJD from time series
            tsMJD%wID = vApp%tilt%wID
            call tsMJD%initTS("MJD",doLoudO=.false.)
            vApp%MJD = T2MJD(vApp%time,tsMJD%evalAt(0.0_rp))
            vApp%gApp%Model%MJD0 = tsMJD%evalAt(0.0_rp)
            !Set first deep coupling (defaulting to coupling immediately)
            call xmlInp%Set_Val(vApp%DeepT, "coupling/tCouple", vApp%time)
            vApp%IO%tCon = vApp%time

            ! correct Gamera start time
            gApp%Model%t = vApp%time / gApp%Model%Units%gT0
            gApp%State%time = gApp%Model%t

            call genVoltShellGrid(vApp, xmlInp)
            call initVoltState(vApp)

        endif

        !Recalculate timestep after correcting Gamera's start and end time
        gApp%Model%dt = CalcDT(gApp%Model,gApp%Grid,gApp%State)
        if (gApp%Model%dt0<TINY) gApp%Model%dt0 = gApp%Model%dt

        if (vApp%doDeep) then
            !Whether to do fast eb-squishing
            call xmlInp%Set_Val(vApp%doQkSquish,"coupling/doQkSquish",.false.)
            call xmlInp%Set_Val(vApp%qkSquishStride,"coupling/qkSquishStride", 2)
            if(vApp%doQkSquish .and. popcnt(vApp%qkSquishStride) /= 1) then
                write(*,*) 'Quick Squish Stride must be a power of 2'
                stop
            endif
            call xmlInp%Set_Val(vApp%chmp2mhd%epsSquish,"ebsquish/epsSquish",0.05)

            !Verify that Gamera has location to hold source info
            if (.not. gApp%Model%doSource) then
                write(*,*) 'Must have GAMERA/source/doSource="T" when running inner magnetosphere model'
                stop
            endif

            !Verify CHIMP data has been set
            if (      (.not. xmlInp%Exists("chimp/units/uid"))     &
               & .or. (.not. xmlInp%Exists("chimp/fields/grType")) &
               & .or. (.not. xmlInp%Exists("chimp/domain/dtype")) ) then   
                write(*,*) 'Necessary CHIMP XML paramters not found, sort that out ...'
                stop
            endif
            
            !Initialize deep coupling type/inner magnetosphere model
            call InitInnerMag(vApp,gApp,xmlInp)

            if(gApp%Model%isRestart) then
                call vApp%imagApp%ReadRestart(gApp%Model%RunID, nRes)
            endif

        endif

        !Check for dynamic coupling cadence
        call xmlInp%Set_Val(vApp%doDynCplDT ,"coupling/doDynDT" , .false.)
        
        if(present(optFilename)) then
            ! read from the prescribed file
            call initializeFromGamera(vApp, gApp, xmlInp, optFilename)
        else
            call initializeFromGamera(vApp, gApp, xmlInp)
        endif

        ! now that remix is initialized, check if precipitation model is OK with deep choice
        if(.not. vApp%doDeep) then
            ! if we aren't using "deep" parts such as RCM, we need to use a
            !    precipitation model that doesn't rely on them
            if(vApp%remixApp%ion(NORTH)%P%aurora_model_type /= FEDDER .and. vApp%remixApp%ion(NORTH)%P%aurora_model_type /= SUNNY) then
                write(*,*) 'The "FEDDER" or SUNNY precipitation model MUST be used when deep coupling is disabled.'
                write(*,*) 'Please either enable the "voltron/coupling/doDeep" option, or'
                write(*,*) ' set "remix/precipitation/aurora_model_type" to "FEDDER" or "SUNNY"'
                stop
            endif
        endif

        if (gApp%Grid%Nkp>=512) then
        !Hex or above, check for sabotage
            !For now disabling hex res for people too lazy to grep this error message
            call xmlInp%Set_Val(isK,"sabotage/isKareem" , .false.)
            if (.not. isK) then
                write(*,*) 'Womp womp womp ...'
                stop
            endif
        endif

        if (vApp%time>=vApp%DeepT) then
            call Tic("DeepCoupling", .true.)
            call DeepUpdate(vApp,gApp)
            call Toc("DeepCoupling", .true.)
        endif

        !Bring overview info
        if (vApp%isLoud) call printConfigStamp()

        !Finally do first output stuff
        !console output
        call consoleOutputV(vApp,gApp)
        !file output
        if (.not. gApp%Model%isRestart) then
            ! write initialization as first output
            ! but save and restore the initial output time so that we don't skip it
            tIO = vApp%IO%tOut
            call fOutputV(vApp, gApp)
            vApp%IO%tOut = tIO
        endif

        ! Do some final checks
        if (abs(gApp%Model%t - gApp%State%time) > TINY) then
            write(*,*)"GAMERA Model%t and State%time don't match, dying"
            write(*,*)"Model: ",gApp%Model%t,", State: ",gApp%State%time
            stop
        endif

        end associate

    end subroutine initVoltron

    !Step Voltron one coupling interval
    subroutine stepVoltron(vApp, dt)
        class(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: dt

        real(rp) :: stepEndTime, remainingStep

        stepEndTime = vApp%time + dt

        do while(stepEndTime .ge. vApp%DeepT)
            ! Finish Gamera processing up to DeepT
            call vApp%gApp%FinishUpdateMhdData(vApp)

            vApp%time = vApp%DeepT
            vApp%MJD = T2MJD(vApp%time,vApp%gApp%Model%MJD0)

            ! update the next predicted coupling interval
            vApp%DeepT = vApp%DeepT + vApp%DeepDT

            call Tic("DeepUpdate",.true.)
            call DeepUpdate(vApp, vApp%gApp)
            call Toc("DeepUpdate",.true.)
            vApp%ts = vApp%ts + 1

            call vApp%gApp%StartUpdateMhdData(vApp)

        enddo

        remainingStep = stepEndTime - vApp%time

        if(remainingStep > TINY) then
            ! partially step Gamera with dts less than coupling interval
            call vApp%gApp%PartialUpdateMhdData(vApp, remainingStep)
        endif

        ! step end time is greater than, or equal to, the current DeepT
        ! advance to that partial deep step time
        vApp%time = max(stepEndTime,vApp%gApp%Model%t*vApp%gApp%Model%Units%gT0)
        vApp%MJD = T2MJD(vApp%time,vApp%gApp%Model%MJD0)

    end subroutine stepVoltron

    !Initialize Voltron app based on Gamera data
    subroutine initializeFromGamera(vApp, gApp, xmlInp, optFilename)
        type(voltApp_T), intent(inout) :: vApp
        class(gamApp_T), intent(inout) :: gApp
        type(XML_Input_T), intent(inout) :: xmlInp
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: RunID, resID
        type(TimeSeries_T) :: f107

        logical :: isRestart
        real(rp) :: maxF107,Rin
        integer :: n, nRes

        isRestart = gApp%Model%isRestart
        RunID = trim(gApp%Model%RunID)
        
        if(vApp%writeFiles) call InitVoltIO(vApp,gApp)
        
    !Remix from Gamera
        !Set mix default grid before initializing
        Rin = norm2(gApp%Grid%xyz(1,1,1,:)) !Inner radius
        call SetMixGrid0(Rin,gApp%Grid%Nkp)
        
        if (gApp%Grid%Nkp>=512) then
            !Hex or above
            call DisableSymLinks()
        endif

        ! Setting apex reference altitude
        call init_apex(vApp%MJD,vApp%gcm%altmax)

        if(present(optFilename)) then
            ! read from the prescribed file
            call init_mix(vApp%remixApp%ion,[NORTH, SOUTH],optFilename=optFilename,RunID=RunID,isRestart=isRestart,nRes=vApp%IO%nRes,optIO=vApp%writeFiles)
        else
            call init_mix(vApp%remixApp%ion,[NORTH, SOUTH],RunID=RunID,isRestart=isRestart,nRes=vApp%IO%nRes,optIO=vApp%writeFiles)
        endif
        
        vApp%remixApp%ion%rad_iono_m  = vApp%planet%ri_m
        vApp%remixApp%ion%rad_planet_m = vApp%planet%rp_m
        !Ensure remix and voltron restart numbers match
        if (isRestart .and. vApp%IO%nRes /= vApp%remixApp%ion(1)%P%nRes) then
            write(*,*) "Voltron and Remix disagree on restart number, you should sort that out."
            write(*,*) "Error code: A house divided cannot stand"
            write(*,*) "   Voltron nRes = ", vApp%IO%nRes
            write(*,*) "   Remix   nRes = ", vApp%remixApp%ion(1)%P%nRes
            stop
        endif

        if (vApp%remixApp%ion(NORTH)%P%doSWF107 .neqv. vApp%remixApp%ion(SOUTH)%P%doSWF107) then
            write(*,*) 'Something is wrong. doSWf107 is set differently for the two hemispheres.'
            write(*,*) 'Exiting ...'
            stop
        endif

        ! read f107 from the SW file and overwrite what's been read from .xml above (in init_mix)
        ! note, only checking for NORTH, because both hemispheres read the same xml file        
        if (vApp%remixApp%ion(NORTH)%P%doSWF107) then
            !Set F10.7 from time series (using max)
            f107%wID = vApp%tilt%wID
            call f107%initTS("f10.7",doLoudO=.false.)
            maxF107 = f107%getMax()

            call updateF107(vApp%remixApp%ion,maxF107)
    
            if (vApp%isLoud) write(*,*) 'Using F10.7 = ', maxF107        
        endif                

        if (vApp%isLoud) write(*,*) 'Using MJD0  = ', gApp%Model%MJD0

        ! initialize remix to gamera structures
        call init_mix2MhdCoupler(vApp%gApp, vApp%remixApp)
        ! initialize additional coupled gamera data
        call vApp%gApp%InitMhdCoupler(vApp)
        if(isRestart) then
            call xmlInp%Set_Val(resID,"/Kaiju/gamera/restart/resID","msphere")
            call xmlInp%Set_Val(nRes,"/Kaiju/gamera/restart/nRes" ,-1)
            call vApp%gApp%ReadRestart(resID, nRes)
        endif

        call init_mhd2Mix(vApp%mhd2mix, gApp, vApp%remixApp)
        !vApp%mix2mhd%mixOutput = 0.0
        
    !CHIMP (TRC) from Gamera
        if (vApp%doDeep) then
            !Verify that there's some place to put deep coupling info
            if (.not. gApp%Model%doSource) then
                write(*,*) 'For deep coupling, GAMERA/source/doSouce="T" must be set!'
                write(*,*) 'Exiting ...'
                stop
            endif

            ! initialize chimp
            if (present(optFilename)) then
                call init_volt2Chmp(vApp,gApp,optFilename=optFilename)
            else
                call init_volt2Chmp(vApp,gApp)
            endif

            !Ensure chimp and voltron restart numbers match
            ! Actually chimp doesn't write restart files right now
            !if (isRestart .and. vApp%IO%nRes /= ebTrcApp%ebModel%IO%nRes) then
            !    write(*,*) "Voltron and Chimp disagree on restart number, you should sort that out."
            !    write(*,*) "Error code: A house divided cannot stand"
            !    write(*,*) "   Voltron nRes = ", vApp%IO%nRes
            !    write(*,*) "   Chimp   nRes = ", ebTrcApp%ebModel%IO%nRes
            !    stop
            !endif

            call init_mhd2Chmp(vApp%mhd2chmp, gApp, vApp%ebTrcApp)
            call init_chmp2Mhd(vApp%chmp2mhd, vApp%ebTrcApp, gApp)
            call init_raiju_mix(vApp%imagApp,vApp%remixApp)

            vApp%iDeep = gApp%Grid%ie-1
            
        endif !doDeep

    end subroutine initializeFromGamera

!----------
    subroutine runRemix(vApp)
        class(voltApp_T), intent(inout) :: vApp
        real(rp) :: curTilt

        ! convert gamera inputs to remix
        call MJDRecalc(vApp%MJD)
        if (vApp%doDeep) then
!            call mapIMagToRemix(vApp%imag2mix,vApp%remixApp) ! rcm style
!            call mapRaijuToRemix(vApp)
            call CoupleIMagToMix(vApp)
        endif
        call mapGameraToRemix(vApp%mhd2mix, vApp%remixApp)

        ! determining the current dipole tilt
        call vApp%tilt%getValue(vApp%time,curTilt)

        ! solve for remix output
        if (vApp%time<=0) then
            call run_mix(vApp%remixApp%ion,curTilt,doModelOpt=.false.)
        else if (vApp%doGCM) then
            call run_mix(vApp%remixApp%ion,curTilt,gcm=vApp%gcm)
        else
            call run_mix(vApp%remixApp%ion,curTilt,doModelOpt=.true.)
        endif

    end subroutine runRemix

!----------
!Deep coupling stuff (time coming from vApp%time, so in seconds)
    subroutine DeepUpdate(vApp, gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Remix code moved from old shallow coupling
        ! convert gamera data to mixInput
        call Tic("G2R")
        call convertGameraToRemix(vApp%mhd2mix, gApp, vApp%remixApp)
        call Toc("G2R")

        ! run remix
        call Tic("ReMIX", .true.)
        call runRemix(vApp)
        call Toc("ReMIX", .true.)

        ! Update potential on voltron's grid
        call updateVoltPotential(vApp)


        call Tic("R2G")
        call CouplePotentialToMhd(vApp)
        call Toc("R2G")

        ! only do imag after spinup with deep enabled
        if(vApp%doDeep .and. vApp%time >= 0) then
            ! DoImag needs PreDeep
            ! DoImag can be in || with squish
            ! PostDeep needs everything
            call PreDeep(vApp, gApp)
              call DoImag(vApp)
              call SquishStart(vApp)
                call Squish(vApp) ! do all squish blocks here
              call SquishEnd(vApp)
            call PostDeep(vApp, gApp)
        elseif(vApp%doDeep) then
            !Doing deep but t<0, call spinup Gas0 (empirical TM plasma sheet)
            !Will kick back w/o doing anything if there's no appetizer

            gApp%Grid%Gas0 = 0
            !Load TM03 into Gas0 for ingestion during spinup
            !Note: Using vApp%time instead of gamera time units
            call LoadSpinupGas0(gApp%Model,gApp%Grid,vApp%time)            
            
        endif

    end subroutine DeepUpdate

    subroutine PreDeep(vApp, gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Update i-shell to trace within in case rTrc has changed
        vApp%iDeep = gApp%Grid%ie-1
        
        !Pull in updated fields to CHIMP
        call Tic("G2C")
        call convertGameraToChimp(vApp%mhd2chmp,gApp,vApp%ebTrcApp)
        call Toc("G2C")
        call Tic("VoltTubes",.true.)
        call genVoltTubes(vApp)
        call Toc("VoltTubes",.true.)

        ! Diagnostics that need fresh voltron information
        call calcPotDrop(vApp)


    end subroutine

    subroutine DoImag(vApp)
        class(voltApp_T), intent(inout) :: vApp

        !Advance inner magnetosphere model to tAdv
        call Tic("InnerMag", .true.)
        ! Do field line tracing using chmpApp and packinto  vApp%imagApp%ijTubes here
        call vApp%imagApp%toIMAG(vApp)
        !call vApp%imagApp%doAdvance(vApp,vApp%DeepT)
        call vApp%imagApp%AdvanceModel(vApp%DeepDT)
        call Toc("InnerMag", .true.)

    end subroutine

    subroutine PostDeep(vApp, gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Now use imag model and squished coordinates to fill Gamera source terms
        call Tic("IM2G")
        call CoupleSourceToMhd(vApp)
        call Toc("IM2G")

    end subroutine

    subroutine CheckQuickSquishError(vApp, gApp, Nbase, Nx2, Nx4, x2Err, x4Err)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp
        integer, intent(out) :: Nbase, Nx2, Nx4
        real(rp), intent(out) :: x2Err, x4Err

        integer :: i,j,k,baseQkSqStr
        logical :: baseDoQkSq
        real(rp), dimension(:,:,:,:), allocatable :: baseXyzSquish
        real(rp), dimension(2) :: posErr
        real(rp) :: dErr

        associate(ebGr=>vApp%ebTrcApp%ebState%ebGr)

        allocate(baseXyzSquish,  MOLD=vApp%chmp2mhd%xyzSquish)
        baseQkSqStr = vApp%qkSquishStride
        baseDoQkSq = vApp%doQkSquish

        ! squish with no quick squish stride
        vApp%qkSquishStride = 1
        vApp%doQkSquish = .false.
        call DeepUpdate(vApp, gApp)
        Nbase = count(NORM2(vApp%chmp2mhd%xyzSquish(ebGr%is:vApp%iDeep+1,ebGr%js:ebGr%je+1,ebGr%ks:ebGr%ke+1,:),dim=4) > TINY)
        baseXyzSquish = vApp%chmp2mhd%xyzSquish

        ! squish with 2x quick squish stride
        vApp%qkSquishStride = 2
        vApp%doQkSquish = .true.
        call DeepUpdate(vApp, gApp)
        Nx2 = count(NORM2(vApp%chmp2mhd%xyzSquish(ebGr%is:vApp%iDeep+1,ebGr%js:ebGr%je+1,ebGr%ks:ebGr%ke+1,:),dim=4) > TINY)
        x2Err = 0
        do i=ebGr%is,vApp%iDeep+1
            do j=ebGr%js,ebGr%je+1
                do k=ebGr%ks,ebGr%ke+1
                    if(NORM2(baseXyzSquish(i,j,k,:)) > TINY .and. NORM2(vApp%chmp2mhd%xyzSquish(i,j,k,:)) > TINY) then
                        dErr = HaverDist(baseXyzSquish(i,j,k,:),vApp%chmp2mhd%xyzSquish(i,j,k,:))
                        x2Err = x2Err + dErr
                    endif
                enddo
            enddo
        enddo

        ! squish with 4x quick squish stride
        vApp%qkSquishStride = 4
        call DeepUpdate(vApp, gApp)
        Nx4 = count(NORM2(vApp%chmp2mhd%xyzSquish(ebGr%is:vApp%iDeep+1,ebGr%js:ebGr%je+1,ebGr%ks:ebGr%ke+1,:),dim=4) > TINY)
        x4Err = 0
        do i=ebGr%is,vApp%iDeep+1
            do j=ebGr%js,ebGr%je+1
                do k=ebGr%ks,ebGr%ke+1
                    if(NORM2(baseXyzSquish(i,j,k,:)) > TINY .and. NORM2(vApp%chmp2mhd%xyzSquish(i,j,k,:)) > TINY) then
                        dErr = HaverDist(baseXyzSquish(i,j,k,:),vApp%chmp2mhd%xyzSquish(i,j,k,:))
                        x4Err = x4Err + dErr
                    endif
                enddo
            enddo
        enddo

        !Rescale to err/pt
        x2Err = x2Err/Nx2
        x4Err = x4Err/Nx4
        
        vApp%qkSquishStride = baseQkSqStr
        vApp%doQkSquish = baseDoQkSq
        deallocate(baseXyzSquish)

        end associate

        contains

        function HaverDist(latlon1,latlon2) result(D)
            real(rp), dimension(2) :: latlon1,latlon2
            real(rp) :: D
            real(rp) :: lat1,lon1,lat2,lon2,dLat,dLon,hArg

            lat1 = latlon1(1)
            lon1 = latlon1(2)
            lat2 = latlon2(1)
            lon2 = latlon2(2)

            dLat = 0.5*(lat2-lat1)
            dLon = 0.5*(lon2-lon1)

            hArg = sin(dLat)**2.0 + cos(lat1)*cos(lat2)*sin(dLon)**2.0
            D = 2*asin(sqrt(hArg))
        end function HaverDist

    end subroutine

    !Initialize CHIMP data structure
    subroutine init_volt2Chmp(vApp,gApp,optFilename)
        class(voltApp_T), intent(inout) :: vApp
        class(gamApp_T), intent(in) :: gApp
        character(len=*), intent(in), optional     :: optFilename

        integer :: b
        character(len=strLen) :: xmlStr
        type(XML_Input_T) :: inpXML
        real(rp) :: xyz0(NDIM)

    !Create input XML object
        if (present(optFilename)) then
            xmlStr = trim(optFilename)
        else
            call getIDeckStr(xmlStr)
        endif
        inpXML = New_XML_Input(trim(xmlStr),"Kaiju/Chimp",.true.)

    !Initialize model
        associate(Model=>vApp%ebTrcApp%ebModel,ebState=>vApp%ebTrcApp%ebState,ebGr=>vApp%ebTrcApp%ebState%ebGr,Gr=>gApp%Grid)
        Model%isMAGE = .true. !Let chimp know it's part of mage
        call setChimpUnitsVoltron(Model,vApp%planet,inpXML)
        Model%T0   = 0.0
        Model%tFin = 0.0
        Model%dt   = 0.0
        Model%t    = 0.0
        call setInterpolation(Model,inpXML)
        Model%doMHD = .true.
        call inpXML%Set_Val(Model%epsds,'tracer/epsds',1.0e-2)    
        call setBackground(Model,inpXML)
        call inpXML%Set_Val(Model%doDip,'tracer/doDip',.false.)
        call setStreamline(Model,inpXML)
        
        !Initialize ebState
        if (gApp%Model%doMultiF) then
            write(*,*) "Initializing MF-Chimp ..."
            !Set proper number of species for chimp
            Model%nSpc = gApp%Model%nSpc        
        endif
        !CHIMP grid is initialized from Gamera's active corners
        call ebInit_fromMHDGrid(Model,ebState,inpXML,Gr%xyz(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,1:NDIM))
        !Replace CHIMP 8-point average centers w/ more accurate Gamera quadrature centers        
        ebGr%xyzcc(ebGr%is:ebGr%ie,ebGr%js:ebGr%je,ebGr%ks:ebGr%ke,:) = Gr%xyzcc(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:)

        call InitLoc(Model,ebState%ebGr,inpXML)

        !Initialize squish indices
        allocate(vApp%ebTrcApp%ebSquish%blockStartIndices(vApp%ebTrcApp%ebSquish%numSquishBlocks))
        do b=1,vApp%ebTrcApp%ebSquish%numSquishBlocks
            vApp%ebTrcApp%ebSquish%blockStartIndices(b) = &
                GetAdjustedSquishStart(vApp,((b-1)*(ebGr%ke+1))/vApp%ebTrcApp%ebSquish%numSquishBlocks)
        enddo

        !Do simple test to make sure locator is reasonable
        xyz0 = Gr%xyz(Gr%is+1,Gr%js,Gr%ks,:)
        if (.not. inDomain(xyz0,Model,ebState%ebGr) ) then
            write(*,*) 'Configuration error: CHIMP Domain incorrect'
            stop
        endif

        end associate

    end subroutine init_volt2Chmp


    subroutine genVoltShellGrid(vApp, xmlInp)
        class(voltApp_T) , intent(inout) :: vApp
        type(XML_Input_T), intent(in) :: xmlInp

        character(len=strLen) :: gType
        integer :: Nt, Np
            !! Number of active cells in theta and phi
        integer :: Ng
            !! Number of ghosts in every direction
        integer, dimension(4) :: nGhosts
            !! Number of ghosts in N, S, E, W directions
        real(rp), dimension(:), allocatable :: phi
            !! Active cell array in phi direction
        real(rp), dimension(:), allocatable :: theta_hemi
            !! Active cell array in theta direction for 1 hemisphere
        real(rp), dimension(:), allocatable :: theta_global
            !! Active cell array in theta direction for both hemispheres
        real(rp) :: sh_radius
            !! ShellGrid radius in Re
        ! Stuff for warp grid
        integer :: nPow
            !! Power applied to non-linear scaling term
        real(rp) :: hWgt
            !! Weight between linear and non-linear term. 1=linear, 0=non-linear
        real(rp) :: xLow, xHigh
            !! Bounds for generating x values between 0 and 1
            !! FIXME: replace with theta_center and x_scale, and calculate our x_low and x_high from that
        integer :: i

        ! Note: Nt is for a single hemisphere, we will manually double it in a minute
        ! TODO: This means we will always have even number of total cells, and a cell interfce right on the equator
        !  Can upgrade to allow for odd number later
        call xmlInp%Set_Val(Nt, "grid/Nt", 180 )  ! 1 deg res default for uniform grid
        call xmlInp%Set_Val(Np, "grid/Np", 360)  ! 1 deg res default
        ! Ghost cells
        call xmlInp%Set_Val(Ng, "grid/Ng", 4)  ! # of ghosts in every direction
        nGhosts = 0
        nGhosts(EAST) = Ng
        nGhosts(WEST) = Ng

        ! Allocate arrays
        allocate(theta_global(Nt*2+1))
        allocate(theta_hemi(Nt+1))
        allocate(phi(Np+1))

        call xmlInp%Set_Val(gType, "grid/gType","UNISPH")
        select case(gType)
            case("UNISPH")
                vApp%gridType = V_GRID_UNIFORM
                call genThetaPhi_uniform(Nt, Np, 0.0_rp, 90.0_rp, theta_hemi, phi)
            case ("WARPSPH")
                vApp%gridType = V_GRID_SHAFEE
                call xmlInp%Set_Val(nPow , "grid/nPow" , 5   )  ! Controls sharpness/bluntness. Must be odd. 5 or 7 are good
                call xmlInp%Set_Val(hWgt , "grid/h"    , 0.2 )  ! Controls width of high-res region
                call xmlInp%Set_Val(xLow , "grid/xLow" , 0.2 )  ! Controls zoom / shift towards pole
                call xmlInp%Set_Val(xHigh, "grid/xHigh", 0.90)  ! Controls zoom / shift towards equator
                call genThetaPhi_Shafee2008(Nt, Np, 0.0_rp, 90.0_rp, nPow,hWgt,xLow,xHigh,theta_hemi,phi)
        end select

        theta_global(1:Nt+1) = theta_hemi
        theta_global(Nt+1:2*Nt+1) = PI - theta_hemi(Nt+1:1:-1)  ! Reverse theta_hemi to mirror into Southern hemisphere
        sh_radius = vApp%planet%ri_m/vApp%planet%rp_m

        ! Now we can make our global ShellGrid
        call GenShellGrid(vApp%shGrid, theta_global, phi, "VOLTRON", nGhosts, sh_radius)

    end subroutine genVoltShellGrid

end module voltapp

