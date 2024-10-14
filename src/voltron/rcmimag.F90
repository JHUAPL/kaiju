!Routines to handle RCM inner magnetosphere model
!NOTES: 
!-Figure out flux-tube volume units
!-Work on upating legacy Fortran
!-Work on OMP bindings
!-Streamline console noise

module rcmimag
    use volttypes
    use files
    use earthhelper
    use imagtubes
    use imaghelper
    use planethelper
    use rcm_mhd_interfaces
    use rcm_mix_interface
    use clocks
    use kronos
    use rcm_mhd_mod, ONLY : rcm_mhd
    use rcm_mhd_io
    use gdefs, only : dFloor,pFloor
    use rcmdefs, only : DenPP0
    use rcmeval
    
    implicit none

    integer, parameter, private :: MHDPad = 0 !Number of padding cells between RCM domain and MHD ingestion
    logical , private :: doTrickyTubes = .true.  !Whether to poison bad flux tubes
    real(rp), private :: imagScl = 1.5 !Safety factor for RCM=>ebsquish
    
    !Whether to call smooth tubes routine at all, see imagtubes for specific options
    logical , private :: doSmoothTubes = .false. 
    
    !Whether to send MHD buffer information to remix
    logical , private :: doBigIMag2Ion = .false. 

    !Whether to use MHD Alfven bounce or RCM hot population bounce
    logical, private :: doHotBounce = .true.
    real(rp), dimension(:,:), allocatable, private :: mixPot

    type, extends(innerMagBase_T) :: rcmIMAG_T

        ! rcm coupling variable
        type(rcm_mhd_T) :: rcmCpl

        ! Holder for field line data
        type(fLine_T), dimension(:,:), allocatable :: rcmFLs

        logical :: doFakeTube=.false. !Only for testing

        contains

        ! over-ride the base functions with RCM versions
        procedure :: doInit => initRCM
        procedure :: doAdvance => advanceRCM
        procedure :: doEval => EvalRCM
        procedure :: doIO => doRCMIO
        procedure :: doConIO => doRCMConIO
        procedure :: doRestart => doRCMRestart

    end type

    contains

    !Initialize RCM inner magnetosphere model
    subroutine initRCM(imag,iXML,isRestart,vApp)
        class(rcmIMAG_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        type(voltApp_T), intent(inout) :: vApp

        character(len=strLen) :: RunID
        real(rp) :: t0

        associate(RCMApp => imag%rcmCpl, & !type rcm_mhd_T
                  imag2mix => vApp%imag2mix, &
                  t0 => vApp%time, &
                  dtCpl => vApp%DeepDT)

        !Set radii in RCMApp
        RCMApp%planet_radius = vApp%planet%rp_m
        RCMApp%iono_radius = vApp%planet%ri_m
        Rp_m = vApp%planet%rp_m  ! For local use
        RIonRCM = vApp%planet%ri_m/vApp%planet%rp_m

        planetM0g = vApp%planet%magMoment

        if (vApp%isLoud) then
            write(*,*) '---------------'
            write(*,*) 'RCM planet params'
            write(*,*) 'Rp        [m]  = ', Rp_m
            write(*,*) 'RIon      [Rp] = ', RIonRCM
            write(*,*) 'RIon      [m]  = ', RIonRCM*Rp_m
            write(*,*) 'MagMoment [G]  = ', planetM0g
            write(*,*) '---------------'
        endif
                
        call iXML%Set_Val(RunID,"/Kaiju/gamera/sim/runid","sim")
        RCMApp%rcm_runid = trim(RunID)

        call iXML%Set_Val(doWolfLim ,"/Kaiju/gamera/source/doWolfLim" ,doWolfLim )
        if (doWolfLim) then
            call iXML%Set_Val(doWolfNLim ,"/Kaiju/gamera/source/doWolfNLim" ,doWolfNLim )
        else
            doWolfNLim = .false.
        endif

        call iXML%Set_Val(doBounceDT,"/Kaiju/gamera/source/doBounceDT",doBounceDT)
        call iXML%Set_Val(doHotBounce,"/Kaiju/gamera/source/doHotBounce",doHotBounce)
        call iXML%Set_Val(nBounce   ,"/Kaiju/gamera/source/nBounce"   ,nBounce   )
        call iXML%Set_Val(maxBetaLim,"/Kaiju/gamera/source/betamax"   ,maxBetaLim)
        call iXML%Set_Val(doBigIMag2Ion ,"imag2ion/doBigIMag2Ion",doBigIMag2Ion)

        call iXML%Set_Val(imagScl ,"imag/safeScl",imagScl)
        call iXML%Set_Val(bMin_C  ,"imag/bMin_C" ,bMin_C )
        call iXML%Set_Val(wImag_C ,"imag/wImag_C",wImag_C)

        call iXML%Set_Val(imag%doFakeTube, "imag/doFakeTube",.false.)

        if (isRestart) then

            !Get t0 and nRes necessary for RCM restart
            call RCMRestartInfo(RCMApp,iXML,t0)

            write(*,*) 'Restarting RCM @ t = ', t0
            vApp%time = t0 !Set vApp's time to correct value from restart
            call rcm_mhd(t0,dtCpl,RCMApp,RCMRESTART,iXML=iXML)
            !Check if we need to do coldstart, assuming coldstart happens at T=0
            if (t0 <= 0) then
                !Still haven't got to T=0 even w/ restart so still need to cold start
                doColdStart = .true.
                call InitRCMICs(imag,vApp,iXML)
            else
                doColdstart = .false. ! set to false if it is a restart
            endif
            call ReadMHD2IMagRestart(imag%rcmCpl,imag%rcmCpl%rcm_nRes-1) !Subtract 1 for the one to read
        else
            t0 = vApp%time
            write(*,*) 'Initializing RCM ...'
            call InitRCMICs(imag,vApp,iXML)
            call rcm_mhd(t0,dtCpl,RCMApp,RCMINIT,iXML=iXML)
            doColdStart = .true.
        endif

        call init_rcm_mix(RCMApp,imag2mix)

        !Allocate any memory needed
        allocate(imag%rcmFLs(RCMApp%nLat_ion,RCMApp%nLon_ion))

        !Start up IO
        if(vApp%writeFiles) call initRCMIO(RCMApp,isRestart)

        end associate
            
    end subroutine initRCM

    !Setup ICs to pass to RCM if asked to
    subroutine InitRCMICs(imag,vApp,iXML)
        class(rcmIMAG_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        type(voltApp_T), intent(in) :: vApp

        real(rp) :: t0
        
        call iXML%Set_Val(RCMICs%doIC,"imag/doInit",.false.)
        t0 = TINY
        if (RCMICs%doIC) then
            !Want initial dst0
            RCMICs%dst0 = GetSWVal("symh",vApp%tilt%wID,t0)
            call iXML%Set_Val(RCMICs%ktRC,"imag/ktRC",30.0)
            RCMICs%vSW = abs(GetSWVal("Vx",vApp%tilt%wID,t0))
            RCMICs%dSW = GetSWVal("D",vApp%tilt%wID,t0)

            !Set PS values (see Borovsky paper)
            RCMICs%dPS  = 0.292*(RCMICs%dSW**0.49)
            RCMICs%kTPS = -3.65 + 0.0190*RCMICs%vSW*1.0e-3 !m/s=>km/s
            RCMICs%kTPS = max(RCMICs%kTPS,TINY)

            !Tune RC pressure profile, using just dst(T=0) (will try to do better later)
            call SetQTRC(RCMICs%dst0)
        else
            !Zero out any additional ring current
            call SetQTRC(0.0_rp)
        endif

        !Also initialize TM03
        call InitTM03(vApp%tilt%wID,t0)

        contains

        function GetSWVal(vID,fID,t0) result(qSW)
            character(len=*), intent(in) :: vID,fID
            real(rp), intent(in) :: t0
            real(rp) :: qSW

            type(TimeSeries_T) :: tsQ
            tsQ%wID = trim(fID)
            call tsQ%initTS(trim(vID),doLoudO=.false.)
            qSW = tsQ%evalAt(t0)
        end function GetSWVal

    end subroutine InitRCMICs

    !Advance RCM from Voltron data
    subroutine AdvanceRCM(imag,vApp,tAdv)
        class(rcmIMAG_T), intent(inout) :: imag
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        integer :: i,j,n,nStp,maxNum
        real(rp) :: colat,lat,lon
        real(rp) :: dtAdv
        type(RCMTube_T) :: ijTube

        real(rp) :: maxRad
        logical :: isLL,doHackIC

        if (vApp%isEarth) then
            call UpdateTM03(vApp%time) !Update plasma sheet model for MP finding and such
            call MJDRecalc(vApp%MJD)
        else
            write(*,*) "You need to do something about RCM for not Earth!"
            stop
        endif

        associate(RCMApp => imag%rcmCpl)

        RCMApp%llBC  = vApp%mhd2chmp%lowlatBC
        RCMApp%dtCpl = vApp%DeepDT
        RCMApp%pFloor = pFloor
        
        call Tic("MAP_RCMMIX")
    !Get potential from mix
        call map_rcm_mix(vApp,mixPot)
        call Toc("MAP_RCMMIX")

        call Tic("RCM_TUBES")
        if (imag%doFakeTube) write(*,*) "Using fake flux tubes for testing!"
            
    !Load RCM tubes
       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP schedule(dynamic) &
       !$OMP private(i,j,colat,lat,lon,isLL,ijTube)
        do j=1,RCMApp%nLon_ion
            do i=1,RCMApp%nLat_ion
                call CleanStream(imag%rcmFLs(i,j)) !Wipe old field line info

                colat = RCMApp%gcolat(i)
                lat = PI/2 - colat
                lon = RCMApp%glong(j)
                
                !Decide if we're below low-lat BC or not
                isLL = (lat <= RCMApp%llBC)
                if (isLL) then
                    !Use mocked up values
                    call DipoleTube(vApp,lat,lon,ijTube,imag%rcmFLs(i,j))
                else
                    if (imag%doFakeTube) then
                        call FakeTube   (vApp,lat,lon,ijTube,imag%rcmFLs(i,j))
                    else
                        !Trace through MHD
                        call MHDTube   (vApp,lat,lon,ijTube,imag%rcmFLs(i,j),vApp%nTrc)
                    endif
                endif
                
                !Stuff data into RCM
                RCMApp%Vol(i,j)          = ijTube%Vol
                RCMApp%bmin(i,j)         = ijTube%bmin
                RCMApp%iopen(i,j)        = ijTube%iopen
                RCMApp%beta_average(i,j) = ijTube%beta_average
                RCMApp%Pave(i,j)         = ijTube%Pave
                RCMApp%Nave(i,j)         = ijTube%Nave
                RCMApp%X_bmin(i,j,:)     = ijTube%X_bmin

                RCMApp%latc(i,j)         = ijTube%latc
                RCMApp%lonc(i,j)         = ijTube%lonc
                RCMApp%losscone(i,j)     = ijTube%losscone
                RCMApp%Lb(i,j)           = ijTube%Lb
                RCMApp%radcurv(i,j)      = ijTube%rCurv
                RCMApp%Tb(i,j)           = ijTube%Tb
                RCMApp%wIMAG(i,j)        = ijTube%wIMAG
                RCMApp%nTrc(i,j)         = imag%rcmFLs(i,j)%Nm+imag%rcmFLs(i,j)%Np
                !mix variables are stored in this order (longitude,colatitude), hence the index flip
                RCMApp%pot(i,j)          = mixPot(j,i)

                !Set composition
                RCMApp%oxyfrac(i,j)      = 0.0

            enddo
        enddo

        doHackIC = (vApp%time <= vApp%DeepDT) .and. RCMICs%doIC !Whether to hack MHD/RCM coupling for ICs

        call Toc("RCM_TUBES")
    !Do any tube hacking needed before sending tubes to RCM
        if (doHackIC) then
        !Tune values to send to RCM for its cold start
            !Setup quiet time ring current to hit target using both current BSDst and target dst
            !Replacing first RC estimate w/ Dst at end of blow-in period
            call SetQTRC(RCMICs%dst0-vApp%BSDst)
            call HackTubes(RCMApp,vApp)
        endif

        if (doTrickyTubes) then
            !Coverup some bad tubes
            call TrickyTubes(RCMApp)
        endif

        if (doSmoothTubes) then
            !Smooth out FTV/potential on tubes b/c RCM will take gradient
            call SmoothTubes(RCMApp,vApp)

        endif

    !Advance from vApp%time to tAdv
        call Tic("AdvRCM")
        dtAdv = tAdv-vApp%time !RCM-DT
        if (doColdstart) then
            write(*,*) 'Cold-starting RCM @ t = ', vApp%time
            call rcm_mhd(vApp%time,dtAdv,RCMApp,RCMCOLDSTART)
            doColdstart = .false.
        else
            call rcm_mhd(vApp%time,dtAdv,RCMApp,RCMADVANCE)
        end if

        call Toc("AdvRCM")

    !Set ingestion region
        if (doHackIC) then
            !For ICs ingest from entire closed field region just this once
            RCMApp%toMHD = .not. (RCMApp%iopen == RCMTOPOPEN)
        else
            call SetIngestion(RCMApp)
            !Try to tailor region to do projections over
            ! !Find maximum extent of RCM domain (RCMTOPCLOSED but not RCMTOPNULL)
            !maxRad = maxval(norm2(RCMApp%X_bmin,dim=3),mask=(RCMApp%iopen == RCMTOPCLOSED))
            maxRad = maxval(norm2(RCMApp%X_bmin,dim=3),mask=.not. (RCMApp%iopen == RCMTOPOPEN))
            maxNum = maxval(      RCMApp%nTrc              ,mask=.not. (RCMApp%iopen == RCMTOPOPEN))
            maxRad = maxRad/Rp_m
            vApp%rTrc = imagScl*maxRad
            vApp%nTrc = min( nint(imagScl*maxNum),MaxFL )
        endif

    !Pull data from RCM state for conductance calculations
        !NOTE: this is not the closed field region, this is actually the RCM domain
        !How much RCM info to use
        if (doBigIMag2Ion) then
            !Pass buffer region to remix
            vApp%imag2mix%inIMag = .not. (RCMApp%iopen == RCMTOPOPEN)
        else    
            vApp%imag2mix%inIMag = (RCMApp%iopen == RCMTOPCLOSED)
        endif

        ! Pass RCM grid type info for remix conductance merge.
        ! Gtype is integer and will be converted to numbers for interpolation in rcm_mix_interface.
        where(RCMApp%iopen == RCMTOPCLOSED)
            vApp%imag2mix%gtype = IMactive
        elsewhere((.not. RCMApp%iopen == RCMTOPCLOSED) .and. (.not. RCMApp%iopen == RCMTOPOPEN))
            vApp%imag2mix%gtype = IMbuffer
        elsewhere
            vApp%imag2mix%gtype = IMoutside
        endwhere
        
        vApp%imag2mix%latc = RCMApp%latc
        vApp%imag2mix%lonc = RCMApp%lonc

    ! electrons precipitation
        vApp%imag2mix%enflx = RCMApp%nflx   (:,:,RCMELECTRON)
        vApp%imag2mix%eflux = RCMApp%flux   (:,:,RCMELECTRON)
        vApp%imag2mix%eavg  = RCMApp%eng_avg(:,:,RCMELECTRON)
    ! ion precipitation
        vApp%imag2mix%inflx = RCMApp%nflx   (:,:,RCMPROTON)
        vApp%imag2mix%iflux = RCMApp%flux   (:,:,RCMPROTON)
        vApp%imag2mix%iavg  = RCMApp%eng_avg(:,:,RCMPROTON)

    ! Pass RCM hot electron density and pressure to REMIX.
        vApp%imag2mix%eden  = RCMApp%Nrcm
        vApp%imag2mix%epre  = RCMApp%Percm
        vApp%imag2mix%npsp  = RCMApp%Npsph

        vApp%imag2mix%isFresh = .true.

        end associate

    end subroutine AdvanceRCM

    !Set region of RCM grid that's "good" for MHD ingestion
    subroutine SetIngestion(RCMApp)
        type(rcm_mhd_T), intent(inout) :: RCMApp

        integer , dimension(:), allocatable :: jBnd
        integer :: i,j
        logical :: inMHD,isClosed
        real(rp) :: Drc,bEq,Lb,Prc

        RCMApp%toMHD(:,:) = .false.
        !Testing lazy quick boundary
        allocate(jBnd (  RCMApp%nLon_ion  ))
        !Now find nominal current boundary
        jBnd(:) = RCMApp%nLat_ion-1

       !$OMP PARALLEL DO default(shared) &
       !$OMP private(i,j,inMHD,isClosed,Drc,bEq,Lb,Prc)
        do j=1,RCMApp%nLon_ion
            do i = RCMApp%nLat_ion,1,-1
                inMHD = RCMApp%toMHD(i,j)
                isClosed = (RCMApp%iopen(i,j) == RCMTOPCLOSED)
                if ( .not. isClosed ) then
                    jBnd(j) = min(i+1+MHDPad,RCMApp%nLat_ion)
                    exit
                endif

            enddo !i loop
            RCMApp%toMHD(:,j) = .false.
            RCMApp%toMHD(jBnd(j):,j) = .true.

            !Replace bounce timescale w/ one using RCM hot population and equatorial B
            if (doHotBounce) then
                !Calculate ingestion timescale in this longitude
                do i = jBnd(j),RCMApp%nLat_ion
                    Drc = rcmNScl*RCMApp%Nrcm (i,j) !#/cc
                    Prc = rcmPScl*RCMApp%Prcm (i,j) !nPa
                    Drc = max(Drc,TINY)
                    bEq = rcmBScl*RCMApp%Bmin (i,j) !Mag field [nT]
                    Lb  = (RCMApp%planet_radius)*(1.0e-3)*RCMApp%Lb(i,j) !Lengthscale [km]
                    !RCMApp%Tb(i,j) = AlfBounce(Drc,bEq,Lb)
                    RCMApp%Tb(i,j) = FastBounce(Drc,Prc,bEq,Lb)
                enddo
            endif

        enddo

        contains
        !Calculate Alfven bounce timescale
        !D = #/cc, B = nT, L = km
        function AlfBounce(Dcc,BnT,Lkm) result(dTb)
            real(rp), intent(in) :: Dcc,BnT,Lkm
            real(rp) :: dTb

            real(rp) :: Va
            if ( (Dcc<TINY) .or. (Lkm<TINY) ) then
                dTb = 0.0
                return
            endif
            Va = 22.0*BnT/sqrt(Dcc) !km/s, from NRL plasma formulary
            dTb = Lkm/Va
        end function AlfBounce

        !Calculate "fast wave" bounce timescale
        !D = #/cc, P = nPa, B = nT, L = km
        function FastBounce(Dcc,Pnpa,Bnt,Lkm) result(dTb)
            real(rp), intent(in) :: Dcc,Pnpa,BnT,Lkm
            real(rp) :: dTb

            real(rp) :: Va,Cs,Tev
            if ( (Dcc<TINY) .or. (Lkm<TINY) ) then
                dTb = 0.0
                return
            endif
            Va = 22.0*BnT/sqrt(Dcc) !km/s, from NRL plasma formulary
            Tev = (1.0e+3)*DP2kT(Dcc,Pnpa) !Temp in eV
            !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
            Cs = 9.79*sqrt( (5.0/3)*Tev )
            dTb = Lkm/sqrt(Va**2.0 + Cs**2.0)
        end function FastBounce

    end subroutine SetIngestion

!-------------
!Eval Routines


    !Evaluate eq map at a given point
    !Returns density (#/cc) and pressure (nPa)
    !x1,x2 = lat,lon
    subroutine EvalRCM(imag,x1,x2,t,imW,isEdible)
        class(rcmIMAG_T), intent(inout) :: imag
        real(rp), intent(in) :: x1,x2,t
        real(rp), intent(out) :: imW(NVARIMAG)
        logical, intent(out) :: isEdible

        !Set defaults
        imW(:) = 0.0
        imW(IMDEN ) = 0.0
        imW(IMPR  ) = 0.0
        imW(IMTSCL) = 0.0
        isEdible = .false.

        !Do quick short circuit test
        if (x1 < TINY) return !lat < TINY

        !Call wrapper for RCM interpolatioon
        call InterpRCM(imag%rcmCpl,x1,x2,t,imW,isEdible)

    end subroutine EvalRCM

!IO wrappers
    subroutine doRCMIO(imag,nOut,MJD,time)
        class(rcmIMAG_T), intent(inout) :: imag
        integer, intent(in) :: nOut
        real(rp), intent(in) :: MJD,time

        !RCM-MHD output
        call WriteRCM   (imag%rcmCpl,nOut,MJD,time)
        !RCM output
        imag%rcmCpl%rcm_nOut = nOut
        call rcm_mhd(time,TINY,imag%rcmCpl,RCMWRITEOUTPUT)
        call WriteRCMFLs(imag%rcmFLs,nOut,MJD,time,imag%rcmCpl%nLat_ion,imag%rcmCpl%nLon_ion)
        
    end subroutine doRCMIO

    subroutine doRCMRestart(imag,nRes,MJD,time)
        class(rcmIMAG_T), intent(inout) :: imag
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time

        imag%rcmCpl%rcm_nRes = nRes
        call rcm_mhd(time,TINY,imag%rcmCpl,RCMWRITERESTART)
        call WriteMHD2IMagRestart(imag%rcmCpl,nRes,MJD,time)
    end subroutine doRCMRestart
    

    !Some quick console diagnostics for RCM info
    subroutine doRCMConIO(imag,MJD,time)
        class(rcmIMAG_T), intent(inout) :: imag
        real(rp), intent(in) :: MJD, time

        integer :: i0,j0,maxIJ(2),maxNum

        real(rp) :: maxPRCM,maxD,maxDP,maxPMHD,maxDMHD,maxL,maxMLT,maxBeta
        real(rp) :: limP,limD,wTrust,wTMin,maxT,maxWT,maxLam,maxLen

        associate(RCMApp => imag%rcmCpl)
    !Start by getting some data
        !Pressure peak info
        maxIJ = maxloc(RCMApp%Prcm,mask=RCMApp%toMHD)
        i0 = maxIJ(1); j0 = maxIJ(2)

        maxPRCM  = RCMApp%Prcm (i0,j0)*rcmPScl
        maxPMHD  = RCMApp%Pave (i0,j0)*rcmPScl
        maxBeta  = RCMApp%beta_average(i0,j0)
        maxD     = RCMApp%Nrcm (i0,j0)*rcmNScl
        maxDMHD  = RCMApp%Nave (i0,j0)*rcmNScl
        maxDP = RCMApp%Npsph(i0,j0)*rcmNScl

        maxL = norm2(RCMApp%X_bmin(i0,j0,XDIR:YDIR))/Rp_m
        maxMLT = atan2(RCMApp%X_bmin(i0,j0,YDIR),RCMApp%X_bmin(i0,j0,XDIR))*180.0/PI
        if (maxMLT<0) maxMLT = maxMLT+360.0
        maxWT = 100*RCMApp%wImag(i0,j0)

        maxLam = (1.0e-3)*RCMApp%MaxAlam/( (RCMApp%vol(i0,j0)*1.0e-9)**(2.0/3.0) ) !Max RCM energy channel [keV]

        !Get pressure weighted confidence
        wTrust = sum(RCMApp%Prcm*RCMApp%wIMAG,mask=RCMApp%toMHD)/sum(RCMApp%Prcm,mask=RCMApp%toMHD)
        wTrust = 100.0*wTrust
        !Get min confidence in MHD domain
        wTMin = 100.0*minval(RCMApp%wIMAG,mask=RCMApp%toMHD)
    !Get some info about size of closed field domain
        maxNum = maxval(RCMApp%nTrc,mask=.not. (RCMApp%iopen == RCMTOPOPEN))
        maxLen = maxval(RCMApp%Lb  ,mask=.not. (RCMApp%iopen == RCMTOPOPEN))

    !Do some output
        if ((maxPRCM<TINY) .or. (time<0)) return

        write(*,*) ANSIYELLOW
        write(*,*) 'RCM'
        !write (*, '(a, f8.2,a,f6.2,a,f6.2,a)')      '  Trust    = ' , wTrust, '% (P-AVG) / ', wTMin, '% (MIN) / ', maxWT, '% (@ MAX)'
        if (doWolfLim) then
            call WolfLimit(maxD,maxPRCM,maxDP,maxDMHD,maxPMHD,maxBeta,limD,limP)
            write (*, '(a, f8.3,a,f8.3,a)')      '  Max RC-P = ' , maxPRCM, ' (RCM) / ', limP, ' (LIM) [nPa]'
            !Get temperature from RCM raw (unlimited)
            maxT = DP2kT(maxD,maxPRCM)
            !maxT = DP2kT(limD,limP)
        else
            write (*, '(a,1f8.3,a)')             '  Max RC-P = ' , maxPRCM, ' [nPa]'
            maxT = DP2kT(maxD,maxPRCM)
        endif

        write (*, '(a,2f8.3,a)')             '   @ L/MLT = ' , maxL, maxMLT, ' [deg]'
        if (doWolfLim) then
            !Add limited density
            write (*, '(a, f8.3,a,f8.3,a,f8.3,a)')      '      w/ D = ' , maxD, ' (RC) / ', maxDP, ' (PSPH) / ', limD, ' (LIM) [#/cc]' 
        else
            write (*, '(a, f8.3,a,f8.3,a)')      '      w/ D = ' , maxD, ' (RC) / ', maxDP, ' (PSPH) [#/cc]'
        endif
        write (*, '(a,1f8.3,a)')             '      w/ T = ' , maxT, ' [keV]'

        write (*, '(a,1f8.3,a)')               '  Max RC-D = ' , maxval(RCMApp%Nrcm,mask=RCMApp%toMHD)*rcmNScl,' [#/cc]'
        write (*,'(a,1f8.3,I6,a)')             '  Max Tube = ', maxLen, maxNum, ' [Re,pts]'
        write(*,'(a,I4,a,I4)')                 '  Channels: ', RCMApp%NkT, ' / ', RCMSIZEK

        write(*,'(a)',advance="no") ANSIRESET!, ''

        end associate

    end subroutine doRCMConIO
end module rcmimag
