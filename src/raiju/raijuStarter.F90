! Initialize kaimag
module raijustarter
    
    ! Base
    use shellgrid
    use xml_input
    use planethelper
    use arrayutil

    ! Raiju
    use raijudefs
    use raijutypes
    use raijugrids
    use raijuetautils
    use raijuRecon, only : maxOrderSupported
    use raijuout
    use raijuICHelpers
    use raijulosses, only : initRaiLosses
    use raijuColdStartHelper, only : initRaijuColdStarter

    ! Cmake points to this
    use raijuuseric

    implicit none

    contains


!------
! Main Initialization Routines
!------

    subroutine raijuInit(app, inpXML, shGridO)
        class(raijuApp_T), intent(inout) :: app
        type(XML_Input_T), intent(in) :: inpXML
        type(ShellGrid_T), intent(in), optional :: shGridO 
            !! If provided, and xml arg is set, we make a child grid for ourselves off of this one

        type(XML_Input_T) :: iXML, iXML_volt
        character(len=strLen) :: tmpStr

        ! Make sure root is Kaiju/raiju
        call inpXML%GetFileStr(tmpStr)
        ! Create new XML reader w/ RAIJU as root
        iXML = New_XML_Input(trim(tmpStr),'Kaiju/RAIJU',.true.)
        
        ! Init model, grid, state
        call raijuInitModel(app%Model, iXML)

        if(present(shGridO)) then
            call raijuInitGrid(app%Model, app%Grid, app%opt, iXML, shGridO)
        else
            call raijuInitGrid(app%Model, app%Grid, app%opt, iXML)
        endif

        call raijuInitState(app%Model,app%Grid,app%State,iXML)

        ! Do losses after everything else has been set, just in case they need something from it
        call initRaiLosses(app%Model, app%Grid, app%State, iXML)

        ! Initialize IOCLOCK
        if (app%Model%isSA) then
            ! If we are driving, we look to raiju section for IO timing
            call app%State%IO%init(iXML,app%State%t,app%State%ts)
        else
            ! If we are not driving, we look to voltron section for IO timing
            iXML_volt = New_XML_Input(trim(tmpStr),'Kaiju/VOLTRON',.true.)
            call app%State%IO%init(iXML_volt,app%State%t,app%State%ts)
        endif


        ! TODO: Add some final checks before returning
        ! e.g., if we are mapping fluids, are we mapping something to HOTP?
        !  (totally just thought of it randomly and not because I spent 4 hours realizing this was my problem)

    end subroutine raijuInit


    ! Sets up Model; Grid and State must be set up separately
    subroutine raijuInitModel(Model, iXML)
        type(raijuModel_T), intent(inout) :: Model
        type(XML_Input_T) , intent(in)    :: iXML
         
        character(len=strLen) :: tmpStr
        character(len=strLen) :: tmpResId
        real(rp) :: cfl0

        write(*,*) "raijuInitModel is starting"

        ! Assume we are running with others, but allow for solo run
        call iXML%Set_Val(Model%isSA, "driver/isSA", .false.)

        ! Assuming that if being controlled by e.g. Voltron, someone else will set RunID accordingly
        ! If we are in SA mode, need to set it ourselves
        if (Model%isSA) then
            call iXML%Set_Val(Model%RunID, "prob/RunID","raijuSA")  ! raiju stand-alone
        else
            call iXML%Set_Val(Model%RunID, "/Kaiju/Gamera/sim/runid","msphere") 
        endif

        ! Timing info, if provided
        call iXML%Set_Val(Model%t0  ,'time/T0',0.0)
        call iXML%Set_Val(Model%tFin,'time/tFin',60.0)
        call iXML%Set_Val(Model%dt  ,'time/dt',1.0)

        ! Config file
        call iXML%Set_Val(Model%configFName, "config/fname","raijuconfig.h5")
        call CheckFileOrDie(Model%configFName,"RAIJU unable to open config file")

        ! Solar wind filename
        call iXML%Set_Val(Model%tsF, '/Kaiju/Gamera/wind/tsFile', 'bcwind.h5')

        call iXML%Set_Val(Model%isMPI, "mpi/isMPI",.false.)
        if (Model%isMPI) then
            write(*,*) "(RAIJU) MPI not implemented yet, dying."
            stop
        endif

        !--- Restarts ---!
        if (Model%isSA) then
            tmpStr = "restart/doRes"
        else
            tmpStr = "/Kaiju/Gamera/restart/doRes"
        endif
        call iXML%Set_Val(Model%isRestart, trim(tmpStr),.false.)
        if (Model%isRestart) then
            ! Restart number
            if (Model%isSA) then
                tmpStr = "restart/nRes"
            else
                tmpStr = "/Kaiju/Gamera/restart/nRes"
            endif
            call iXML%Set_Val(Model%nResIn, trim(tmpStr), Model%nResIn)
            ! Restart id (allows it to be different from our current runId)
            if (Model%isSA) then
                tmpStr = "restart/resId"
            else
                tmpStr = "/Kaiju/Gamera/restart/resId"
            endif
            call iXML%Set_Val(tmpResId, trim(tmpStr), Model%RunID)
            call genResInFname(Model, Model%ResF, runIdO=tmpResId)  ! Determine filename to read from
        endif
        
        !--- Plasmasphere ---!
        call iXML%Set_Val(Model%doPlasmasphere, "plasmasphere/doPsphere",.true.)
        call iXML%Set_Val(Model%doPsphEvol, 'plasmasphere/doEvol',.true.)
        call iXML%Set_Val(Model%psphEvolRad, 'plasmasphere/evolRad', def_psphEvolRad)
        ! Determine number of species. First set default, then read from xml to overwrite if present
        if (Model%doPlasmasphere) then
            Model%nSpc = 3
            call iXML%Set_Val(Model%psphInitKp, "plasmasphere/initKp",3.0)
        else
            Model%nSpc = 2
        endif
        call iXML%Set_Val(Model%doExcessToPsph, "prob/doExcessMap",.true.)
            !! Allow mapping of excess H+ to plasmasphere channel
        if (Model%doExcesstoPsph .and. .not. Model%doPlasmasphere) then
            write(*,*)"Error: Can't map excess to psph because psph doesn't exist!"
        endif

        call iXML%Set_Val(Model%nSpc, "prob/nSpc",Model%nSpc)

        !--- Domain ---!
        call iXML%Set_Val(Model%n_bndLim      , "domain/bndLim"     , 3)
        call iXML%Set_Val(Model%maxTail_buffer, "domain/tail_buffer", def_maxTail_buffer)
        call iXML%Set_Val(Model%maxSun_buffer , "domain/sun_buffer" , def_maxSun_buffer)
        call iXML%Set_Val(Model%maxTail_active, "domain/tail_active", def_maxTail_active)
        call iXML%Set_Val(Model%maxSun_active , "domain/sun_active" , def_maxSun_active)
        call iXML%Set_Val(Model%activeDomRad  , "domain/activeRad"  , def_activeDomRad)
        ! Store all distances as positive values, we'll add signs as needed later
        Model%maxTail_buffer = abs(Model%maxTail_buffer)
        Model%maxSun_buffer  = abs(Model%maxSun_buffer)
        Model%maxTail_active = abs(Model%maxTail_active)
        Model%maxSun_active  = abs(Model%maxSun_active)
        Model%activeDomRad   = abs(Model%activeDomRad)

        !---Solver ---!
        call iXML%Set_Val(Model%doUseVelLRs,'sim/useVelLRs',def_doUseVelLRs)
        call iXML%Set_Val(Model%maxItersPerSec,'sim/maxIter',def_maxItersPerSec)
        call iXML%Set_Val(Model%maxOrder,'sim/maxOrder',7)
        if (Model%maxOrder > maxOrderSupported) then
            write(*,*) "(RAIJU) Too much order, allowable orders <= ",maxOrderSupported
            stop
        endif
        call iXML%Set_Val(Model%PDMB,'sim/pdmb',def_pdmb)
        cfl0 = min(0.5/(Model%PDMB+0.5), cflMax) !Set CFL based on PDM

        !Set CFL from XML
        call iXML%Set_Val(Model%CFL ,'sim/cfl'  ,cfl0)
        if (Model%CFL > cfl0) then
            write(*,*) '-------------------------------------'
            write(*,*) '(RAIJU) WARNING, CFL is above critical value!'
            write(*,*) 'CFL/Critical/PDMB = ', Model%CFL,cfl0,Model%PDMB
            write(*,*) '-------------------------------------'
        else
            write(*,*) 'CFL # = ', Model%CFL
        endif

        ! Certain physical constants that shouldn't be constants
        call iXML%Set_Val(Model%tiote, "prob/tiote",4.0)

        ! Active shell settings
        call iXML%Set_Val(Model%doActiveShell, "activeShell/doAS",.false.)
        call iXML%Set_Val(Model%worthyFrac, "activeShell/worthyFrac",fracWorthyDef)

        ! Corotation
        if (Model%isSA) then
            call iXML%Set_Val(Model%doOwnCorot, "prob/doCorot",.true.)
        else
            call iXML%Set_Val(Model%doOwnCorot, "prob/doCorot",.false.)
        endif
        if (Model%doOwnCorot) then
            call iXML%Set_Val(Model%doGeoCorot, "prob/doGeoCorot",.false.)
        endif

        ! Lambda channel settings
        call iXML%Set_Val(Model%doDynamicLambdaRanges, "lambdas/dynamicRanges",.false.)

        ! Determine which DP2eta mapping should be used (e.g. Maxwellian, Kappa, user-defined)
        call iXML%Set_Val(tmpStr, "moments/DP2EtaMap","MAXWELL")
        select case(tmpStr)
            case("MAXWELL")
                Model%dp2etaMap => Maxwell2Eta
            case("KAPPA")
                Model%dp2etaMap => Kappa2Eta
                call iXML%Set_Val(Model%kappa, "moments/kappa",6.0)
            ! TODO: Add user option which will point to whatever is in the chosen IC file
            case DEFAULT
                write(*,*) "(RAIJU) Received unavailable DP2Eta mapping: ",tmpStr
                write(*,*) " Dying."
                stop
        end select


        !--- Losses ---!
        call iXML%Set_Val(Model%doLosses, "losses/doLosses",.true.)
        !call iXML%Set_Val(Model%doSS    , "losses/doSS" ,.true. )
        call iXML%Set_Val(Model%doCC    , "losses/doCC" ,.true. )
        call iXML%Set_Val(Model%doCX    , "losses/doCX" ,.true.)
        !call iXML%Set_Val(Model%doFLC   , "losses/doFLC",.false.)
        call iXML%Set_Val(Model%doEWM   , "losses/doEWM",.true.)
        call iXML%Set_Val(Model%ewmType , "losses/ewmType","SS")

        !--- Output ---!
        call iXML%Set_Val(Model%isLoud           , "output/loudConsole",.false.)
        call iXML%Set_Val(Model%writeGhosts      , "output/writeGhosts",.false.)
        call iXML%Set_Val(Model%doOutput_potGrads, "output/doFat"      ,.false.)
        call iXML%Set_Val(Model%doOutput_3DLoss  , "output/doLossExtras" ,.false.)  ! Several (Ni,Nj,Nk) arrays
        call iXML%Set_Val(Model%doOutput_debug   , "output/doDebug"    ,.false.)

        ! Model Hax?
        call iXML%Set_Val(Model%doHack_rcmEtaPush, "hax/rcmEtaPush",.false.)
        

        ! Determine coupling information
        call setMHD2RaiInfo(Model, iXML)

        ! Set planet params
        call getPlanetParams(Model%planet, iXML, doLoudO=.true.)

    end subroutine raijuInitModel


    subroutine setMHD2RaiInfo(Model, iXML)
        type(raijuModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: iXML

        integer :: nFluids
        integer :: f
        character(len=strLen) :: xmlFluidTag


        ! Domain determination
        call iXML%Set_Val(Model%lim_vaFrac_soft, "domain/vaFrac_soft", def_lim_vaFrac_soft)
        call iXML%Set_Val(Model%lim_vaFrac_hard, "domain/vaFrac_hard", def_lim_vaFrac_hard)
        call iXML%Set_Val(Model%lim_bmin_soft  , "domain/bmin_soft"  , def_lim_bmin_soft  )
        call iXML%Set_Val(Model%lim_bmin_hard  , "domain/bmin_hard"  , def_lim_bmin_hard  )
        !call iXML%Set_Val(Model%vaFracThresh , "cpl/vaFracThresh" , def_vaFracThresh)
        !call iXML%Set_Val(Model%bminThresh   , "cpl/bminThresh"   , def_bminThresh)
        !call iXML%Set_Val(Model%PstdThresh   , "cpl/Pstd"         , def_PstdThresh)
        !call iXML%Set_Val(Model%normAngThresh, "cpl/normAngThresh", def_normAngle)
        ! Note: Not used now within raiju, just passing through. But should probably include as a constraint
        call iXML%Set_Val(Model%nBounce, "cpl/nBounce", def_nBounce)
        ! Convert degrees to dot product value
        !Model%normAngThresh = cos(Model%normAngThresh*PI/180.0)

        ! Fluid mapping
        call iXML%Set_Val(nFluids, "cpl/nFluidsIn",0)
        Model%nFluidIn = nFluids
        
        if (nFluids > 0) then
            allocate(Model%fluidInMaps(nFluids))
            do f=1,nFluids
                write(xmlFluidTag,'(A,I0,A)')"fluidIn",f,"/imhd"
                call iXML%Set_Val(Model%fluidInMaps(f)%idx_mhd, xmlFluidTag, 0)
                write(xmlFluidTag,'(A,I0,A)')"fluidIn",f,"/flav"
                call iXML%Set_Val(Model%fluidInMaps(f)%flav, xmlFluidTag, -1)
                write(xmlFluidTag,'(A,I0,A)')"fluidIn",f,"/excessToPsph"
                call iXML%Set_Val(Model%fluidInMaps(f)%doExcessToPsph, xmlFluidTag, .false.)
            enddo
        endif

    end subroutine setMHD2RaiInfo


    subroutine raijuInitGrid(Model, Grid, opts, iXML, shGridO)
        type(raijuModel_T) , intent(inout) :: Model
        type(raijuGrid_T)  , intent(inout) :: Grid
        type(raijuOptions_T), intent(in) :: opts
        type(XML_Input_T), intent(in)   :: iXML
        type(ShellGrid_T), optional, intent(in)    :: shGridO

        character(len=strLen) :: tmpStr
        type(ShellGrid_T) :: shGrid
        logical :: doHack_constB

        ! Set grid params
        call iXML%Set_Val(Grid%nB, "grid/Nbnd", 4      )  ! Number of cells between open boundary and active domain
        call iXML%Set_Val(tmpStr, "grid/gType","UNISPH")

        ! Fill out Grid object depending on chosen method
        select case(tmpStr)
            case("UNISPH")
                Grid%gType = RAI_G_UNISPH
                ! Generate our own grid from scratch
                call raijuGenUniSphGrid(Model, Grid, iXML)
            case("WARPSPH")
                Grid%gType = RAI_G_WARPSPH
                ! Generate our own grid from scratch
                !call raijuGenWarpSphGrid_Fok2021(Model, Grid, iXML)
                call raijuGenWarpSphGrid_Shafee2008(Model, Grid, iXML)
            case("SHGRID")
                Grid%gType = RAI_G_SHGRID
                ! Then we should be receiving a predefined ShellGrid that Voltron has set up
                if(present(shGridO)) then
                    shGrid = shGridO
                    call raijuGenGridFromShGrid(Grid%shGrid, shGrid, iXML, opts)
                else
                    write(*,*) "RAIJU expecting a ShellGrid_T but didn't receive one. Dying."
                endif
            case DEFAULT
                write(*,*) "RAIJU Received invalid grid definition: ",Grid%gType
                write(*,*) " Dying."
                stop
        end select
        ! Make sure we aren't going beyond the equator with our ghost cells, that wouldn't make sense
        if (any(Grid%shGrid%th(Grid%shGrid%ie+1:) > PI/2.0_rp)) then
            write(*,*)"ERROR: RAIJU ghosts going beyond 0 deg lat. Probably you are pushing it in too close to planet's surface"
            stop
        endif

        ! Finalize the spatial part of the grid
        call finalizeLLGrid(Grid, Model%planet)

        ! Hax?
        call iXML%Set_Val(doHack_constB, "hax/doConstB",.false.)
        if (doHack_constB) then
            Grid%Bmag = Model%planet%magMoment*G2nT
            Grid%Brcc = Model%planet%magMoment*G2nT
            Grid%BrFace = Model%planet%magMoment*G2nT
        endif

        ! Now handle lambda grid
        Grid%nSpc = Model%nSpc  ! Make a copy of nSpc
        Grid%nFluidIn = Model%nFluidIn  ! Make a copy of nFluidIn
        call iXML%Set_Val(Grid%ignoreConfigMismatch, "config/ignoreMismatch",.false.)
        call initLambdaGrid(Model, Grid, Model%configFName)

    end subroutine raijuInitGrid


    subroutine raijuInitState(Model, Grid, State, iXML)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T) , intent(in)    :: Grid
        type(raijuState_T), intent(inout) :: State
        type(XML_Input_T), intent(in)   :: iXML

        integer :: s

        ! Allocate arrays

        associate(sh=>Grid%shGrid)

            ! dt for every lambda channel
            allocate( State%dtk (Grid%Nk) )
            call fillArray(State%dtk, 0.0_rp)
            ! nSteps for each channel
            allocate( State%nStepk(Grid%Nk) )
            call fillArray(State%nStepk, 0)
            ! Where we keep all our stuff
            allocate( State%eta      (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%eta, 0.0_rp)
            ! Where we keep all our stuff but a half-step ahead of now
            allocate( State%eta_half (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%eta_half, 0.0_rp)
            ! Where we kept all our stuff one step ago
            allocate( State%eta_last (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%eta_last, 0.0_rp)
            ! Where all the stuff sorta was over the last State%dt
            allocate( State%eta_avg  (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%eta_avg, 0.0_rp)
            ! I shells shat should be evolved for each k
            allocate( State%activeShells (sh%isg:sh%ieg, Grid%Nk) )
            State%activeShells = .false.
            ! Effective potential (used for output only)
            allocate( State%pEff(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk) )
            call fillArray(State%pEff, 0.0_rp)
            ! Gradient of ionspheric potential
            allocate( State%gradPotE     (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 2) )
            call fillArray(State%gradPotE, 0.0_rp)
            ! Gradient of corotation potential
            allocate( State%gradPotCorot (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 2) )
            call fillArray(State%gradPotCorot, 0.0_rp)
            ! Gradient of (flux tube volume ^ -2/3)
            allocate( State%gradVM      (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 2) )
            call fillArray(State%gradVM, 0.0_rp)
            ! Interface and cell velocities
            allocate( State%gradPotE_cc    (sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
            call fillArray(State%gradPotE_cc, 0.0_rp)
            allocate( State%gradPotCorot_cc(sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
            call fillArray(State%gradPotCorot_cc, 0.0_rp)
            allocate( State%gradVM_cc      (sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
            call fillArray(State%gradVM_cc, 0.0_rp)
            allocate( State%iVel (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
            call fillArray(State%iVel, 0.0_rp)
            allocate( State%iVelL(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
            call fillArray(State%iVelL, 0.0_rp)
            allocate( State%iVelR(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
            call fillArray(State%iVelR, 0.0_rp)
            allocate( State%cVel (sh%isg:sh%ieg  , sh%jsg:sh%jeg  , Grid%Nk, 2) )
            call fillArray(State%cVel, 0.0_rp)
            
            ! Coupling input moments
            allocate( State%Pavg(sh%isg:sh%ieg  , sh%jsg:sh%jeg, 0:Grid%nFluidIn) )
            call fillArray(State%Pavg, 0.0_rp)
            allocate( State%Davg(sh%isg:sh%ieg  , sh%jsg:sh%jeg, 0:Grid%nFluidIn) )
            call fillArray(State%Davg, 0.0_rp)
            allocate( State%Pstd(sh%isg:sh%ieg  , sh%jsg:sh%jeg, 0:Grid%nFluidIn) )
            call fillArray(State%Pstd, 0.0_rp)
            allocate( State%Dstd(sh%isg:sh%ieg  , sh%jsg:sh%jeg, 0:Grid%nFluidIn) )
            call fillArray(State%Dstd, 0.0_rp)
            allocate( State%domWeights(sh%isg:sh%ieg  , sh%jsg:sh%jeg) )
            call fillArray(State%domWeights, 0.0_rp)
            allocate( State%tiote(sh%isg:sh%ieg  , sh%jsg:sh%jeg) )
            call fillArray(State%tiote, 0.0_rp)
            call initShellVar(Grid%shGrid, SHGR_CC, State%Tb)
            State%Tb%data = 0.0
            State%Tb%mask = .false.
            ! Bmin surface
            allocate( State%Bmin    (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 3 ) )
            call fillArray(State%Bmin, 0.0_rp)
            allocate( State%xyzMin  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 3 ) )
            call fillArray(State%xyzMin, 0.0_rp)
            allocate( State%xyzMincc(sh%isg:sh%ieg  , sh%jsg:sh%jeg  , 3 ) )
            call fillArray(State%xyzMincc, 0.0_rp)
            allocate( State%thcon   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1    ) )
            call fillArray(State%thcon, 0.0_rp)
            allocate( State%phcon   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1    ) )
            call fillArray(State%phcon, 0.0_rp)
            ! 2D corner quantities
            allocate( State%topo     (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            call fillArray(State%topo, 0)
            allocate( State%espot    (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            call fillArray(State%espot, 0.0_rp)
            allocate( State%pot_corot(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            call fillArray(State%pot_corot, 0.0_rp)
            allocate( State%bvol     (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            call fillArray(State%bvol, 0.0_rp)
            allocate( State%bvol_cc  (sh%isg:sh%ieg  , sh%jsg:sh%jeg  ) )
            call fillArray(State%bvol_cc, 0.0_rp)
            allocate( State%vaFrac   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            call fillArray(State%vaFrac, 0.0_rp)
            ! 1D cell-centered quantities
            allocate( State%bndLoc(sh%jsg:sh%jeg) )
            call fillArray(State%bndLoc, 0)
            ! 2D cell-centered quantities
            allocate( State%active      (sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            call fillArray(State%active, 0)
            allocate( State%active_last (sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            call fillArray(State%active_last, 0)
            allocate( State%OCBDist(sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            call fillArray(State%OCBDist, 0)
            
            allocate( State%lossRates      (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%lossRates, 0.0_rp)
            allocate( State%precipType_ele (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%precipType_ele, 0.0_rp)
            allocate( State%lossRatesPrecip(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%lossRatesPrecip, 0.0_rp)
            !allocate( State%precipNFlux    (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            !allocate( State%precipEFlux    (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%dEta_dt        (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%dEta_dt, 0.0_rp)
            allocate( State%CCHeatFlux     (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            call fillArray(State%CCHeatFlux, 0.0_rp)
            ! Coupling output data
            allocate(State%Den  (0:Model%nSpc))
            allocate(State%Press(0:Model%nSpc))
            allocate(State%Den_avg  (0:Model%nSpc))
            allocate(State%Press_avg(0:Model%nSpc))
            do s=0,Model%nSpc
                call initShellVar(Grid%shGrid, SHGR_CC, State%Den  (s))
                call initShellVar(Grid%shGrid, SHGR_CC, State%Press(s))
                State%Den  (s)%data = 0.0
                State%Press(s)%data = 0.0
                State%Den  (s)%mask = .false.
                State%Press(s)%mask = .false.

                call initShellVar(Grid%shGrid, SHGR_CC, State%Den_avg  (s))
                call initShellVar(Grid%shGrid, SHGR_CC, State%Press_avg(s))
                State%Den_avg  (s)%data = 0.0
                State%Press_avg(s)%data = 0.0
                State%Den_avg  (s)%mask = .false.
                State%Press_avg(s)%mask = .false.
            enddo
            allocate(State%precipNFlux(Grid%Nk))
            allocate(State%precipEFlux(Grid%Nk))
            do s=1,Grid%Nk
                call initShellVar(Grid%shGrid, SHGR_CC, State%precipNFlux(s))
                call initShellVar(Grid%shGrid, SHGR_CC, State%precipEFlux(s))
                State%precipNFlux(s)%data = 0.0
                State%precipEFlux(s)%data = 0.0
                State%precipNFlux(s)%mask = .false.
                State%precipEFlux(s)%mask = .false.
            enddo

            ! Only bother allocating persistent versions of debug stuff if we need them
            if (Model%doOutput_debug) then
                allocate( State%etaFaceReconL(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                call fillArray(State%etaFaceReconL, 0.0_rp)
                allocate( State%etaFaceReconR(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                call fillArray(State%etaFaceReconR, 0.0_rp)
                allocate( State%etaFacePDML  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                call fillArray(State%etaFacePDML, 0.0_rp)
                allocate( State%etaFacePDMR  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                call fillArray(State%etaFacePDMR, 0.0_rp)
                allocate( State%etaFlux      (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                call fillArray(State%etaFlux, 0.0_rp)
            endif

            State%KpTS%wID = Model%tsF
            call State%KpTS%initTS("Kp")
        
        end associate
        
        ! For now, just set t to tStart and ts to 0
        State%t = Model%t0
        State%ts = 0
        State%nStepk = 0

        ! Set all activeShells to true. There should be doActiveShell checks anywhere it is used, but just in case, make sure default has all shells on
        State%activeShells = .true.
        ! Similarly, set vaFrac to safe value in case stand-alone never writes to it
        State%vaFrac = 1.0

        ! Init State sub-modules
        if (Model%isSA) then
            ! If we are standalone, this is the place to get coldStarter settings
            ! Otherwise, raiCpl should call it so we can use other models' info
            call initRaijuColdStarter(Model, iXML, State%coldStarter)
        endif

    ! Init state
        call iXML%Set_Val(Model%icStr, "prob/IC","DIP")
        select case(trim(Model%icStr))
            case("DIP")
                !! Simple Maxwellian in a dipole field
                !Model%initState => initRaijuIC_DIP
                call initRaijuIC_DIP(Model, Grid, State, iXML)
            case("TTA")
                !! Test Theta Advection
                !Model%initState => initRaijuIC_testThetaAdvection
                call initRaijuIC_testThetaAdvection(Model, Grid, State, iXML)
            case("USER")
                ! Call the IC in the module raijuuseric
                ! This module is set in cmake via the RAIJUIC variable
                !Model%initState => userInitStateFunc
                call raijuInitState_useric(Model, Grid, State, iXML)       
            case DEFAULT
                write(*,*)"Invalid IC name to RAIJU, see raijuStarter.F90:raijuInitState. Bye."
                stop
        end select

        !call Model%initState(Model, Grid, State, iXML)

    end subroutine raijuInitState


end module raijustarter
