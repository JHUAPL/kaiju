module raijuIO
    use ioh5
    use files
    use planethelper
    use kai2geo
    use shellGridIO
    
    use raijutypes
    use raijuetautils
    use shellGrid
    use shellInterp
    use shellGridIO
    use raijuDomain, only : calcMapJacNorm

    implicit none

    integer, parameter, private :: MAXIOVAR = 200
    logical, private :: doRoot = .true. !Whether root variables need to be written to main output file
    logical, private :: doFat = .false. !Whether to output lots of extra datalogical, private :: doRoot = .true. !Whether root variables need to be written

    contains

    subroutine raijuInitIO(Model, Grid, doGhostsO)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        logical, optional, intent(in) :: doGhostsO

        logical :: doGhosts
        logical :: fExist
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

        doRoot = .false. ! Don't call again

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        Model%raijuH5 = trim(Model%RunID) // ".raiju.h5"

        fExist = CheckFile(Model%raijuH5)
        write(*,*) "RAIJU outputting to ",trim(Model%raijuH5)

        if(.not. Model%isRestart) then
            ! Remove all old files, start fresh
            call CheckAndKill(Model%raijuH5)
        endif

        if (Model%isRestart .and. fExist) then
            ! No need to do anything because file already exist
            !!TODO: When more output files are active, maybe put in a check for each of them
            return
        endif

        ! If still here, proceed to init
        call raijuWriteRootH5Data(Model, Grid, fnameO=Model%raijuH5, doGhostsO=doGhosts)
        
    end subroutine raijuInitIO


    subroutine raijuWriteRootH5Data(Model, Grid, fnameO, doGhostsO)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        character(len=strLen), optional, intent(in) :: fnameO
        logical, optional, intent(in) :: doGhostsO

        integer :: is, ie, js, je
        logical :: doGhosts
        integer :: i
        logical :: fExist
        real(rp), dimension(:,:), allocatable :: lat2D, lon2D
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        character(len=strLen) :: fStr, gStr

        if(present(fnameO)) then
            fStr = trim(fnameO)
        else
            fStr = trim(Model%raijuH5)
        endif

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        if (doGhosts) then
            is = Grid%shGrid%isg
            ie = Grid%shGrid%ieg
            js = Grid%shGrid%jsg
            je = Grid%shGrid%jeg
        else
            is = Grid%shGrid%is
            ie = Grid%shGrid%ie
            js = Grid%shGrid%js
            je = Grid%shGrid%je
        endif

        fExist = CheckFile(fStr)

        associate(sh=>Grid%shGrid, spc=>Grid%spc)
        ! Add spatial grid as 2D array
        allocate(lat2D(is:ie+1,js:je+1))  ! +1 because we're doing corners
        allocate(lon2D(is:ie+1,js:je+1))

        do i=js,je+1
            lat2D(:,i) = sh%th(is:ie+1)
        enddo

        do i=is,ie+1
            lon2D(i,:) = sh%ph(js:je+1)
        enddo

        ! Some model setting info
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"doDebugOutput" ,Model%doOutput_debug   )  ! Attr
        call AddOutVar(IOVars,"doWriteGhosts" ,Model%writeGhosts      )  ! Attr
        call AddOutVar(IOVars,"doGeoCorot"    ,Model%doGeoCorot       )  ! Attr
        call AddOutVar(IOVars,"doOwnCorot"    ,Model%doOwnCorot       )  ! Attr
        call AddOutVar(IOVars,"doExcesstoPsph",Model%doExcesstoPsph   )  ! Attr
        call AddOutVar(IOVars,"doPlasmasphere",Model%doPlasmasphere   )  ! Attr
        call AddOutVar(IOVars,"doActiveShell",Model%doActiveShell     )  ! Attr
        call AddOutVar(IOVars,"nSpcMHD",Model%nSpcMHD)  ! Attr
        call AddOutVar(IOVars,"nSpcRAIJU",Model%nSpc)  ! Attr

        call AddOutVar(IOVars,"X",lat2D,uStr="radians")
        call AddOutVar(IOVars,"Y",lon2D,uStr="radians")
        call WriteVars(IOVars,.true.,fStr)
        
        ! Grid data
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"th_surf"  ,Grid%thRp  (is:ie+1),uStr="radians",dStr="(corner) Thetas mapped to 1 Rp")
        call AddOutVar(IOVars,"thcc_surf",Grid%thcRp (is:ie  ),uStr="radians",dStr="(cell) Thetas mapped to 1 Rp"  )
        call AddOutVar(IOVars,"Bmag"     ,Grid%Bmag  (is:ie,js:je),uStr="nT")
        call AddOutVar(IOVars,"BrCC"     ,Grid%BrCC  (is:ie,js:je),uStr="nT")
        call AddOutVar(IOVars,"areaCC"   ,Grid%areaCC(is:ie,js:je),uStr="Rp^2")
        call AddOutVar(IOVars,"BrFace"   ,Grid%BrFace(is:ie+1,js:je+1,:),uStr="nT")
        call AddOutVar(IOVars,"alamc"    ,Grid%alamc,uStr="eV * (Rp/nT)^(2/3)")
        call AddOutVar(IOVars,"areaFace",Grid%areaFace(is:ie+1,js:je+1,:),uStr="Rp^2")
        call AddOutVar(IOVars,"lenFace" ,Grid%lenFace (is:ie+1,js:je+1,:),uStr="Rp")
        call WriteVars(IOVars,.true.,fStr,"Grid")

        call writeShellGrid(Grid%shGrid, fStr, "Grid/ShellGrid")

        ! Output detailed lambda grid info
        call writeSpeciesInfo(Model, Grid, fStr)

        ! Output planet info
        call writePlanetParams(Model%planet, .true., fStr)
        end associate

    end subroutine raijuWriteRootH5Data


    subroutine WriteRaiju(Model, Grid, State, gStr, doGhostsO)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: gStr
        logical, optional, intent(in) :: doGhostsO

        character(len=strLen) :: tmpStr
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        integer :: i,j,k,s, nClkSteps
        integer :: is, ie, js, je, ks, ke
        integer, dimension(4) :: outBnds2D, outBnds2D_corner
        logical :: doGhosts
        real(rp) :: axT, axP
        real(rp), dimension(:,:), allocatable :: outTmp2D
        real(rp), dimension(:,:,:), allocatable :: outTmp3D
        real(rp), dimension(:,:,:), allocatable :: outPrecipN, outPrecipE, outDeltaPressure, outCCHeatFlux
        real(rp), dimension(Model%nSpc) :: dpsdst_spc

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        if (doGhosts) then
            is = Grid%shGrid%isg
            ie = Grid%shGrid%ieg
            js = Grid%shGrid%jsg
            je = Grid%shGrid%jeg
        else
            is = Grid%shGrid%is
            ie = Grid%shGrid%ie
            js = Grid%shGrid%js
            je = Grid%shGrid%je
        endif
        outBnds2D = (/is,ie,js,je/)
        outBnds2D_corner = (/is,ie+1,js,je+1/)

        ! First, make sure root variables are there
        if (doRoot) then
            call raijuInitIO(Model, Grid, doGhosts)
            doRoot = .false.
        endif
        !Reset IO chain
        call ClearIO(IOVars)

        ! Calc dpsdst to add to attributes
        do s=1,Model%nSpc
            dpsdst_spc(s) = spcEta2DPS(Model, Grid, State%bvol_cc, State%eta, Grid%spc(s), State%active .eq. RAIJUACTIVE)
            write(tmpStr, '(A,I0)')"DPSDst_",Grid%spc(s)%flav
            call AddOutVar(IOVars, trim(tmpStr), dpsdst_spc(s))  ! Save as attribute
        enddo
        call AddOutVar(IOVars,"DPSDst_total",sum(dpsdst_spc)) !Save total over all species
        
        ! Add attributes
        call AddOutVar(IOVars,"time",State%t)
        call AddOutVar(IOVars,"MJD",State%mjd)
        call AddOutVar(IOVars,"nOut",State%IO%nOut)
        call AddOutVar(IOVars,"nRes",State%IO%nRes)
        call AddOutVar(IOVars,"isFirstCpl", merge(1,0,State%isFirstCpl))
        call WriteRaijuColdStartInfo(State, IOVars)

        ! Save axis of rotation as attributes
        if (Model%doGeoCorot) then
            call GeoAxisTP(axT, axP)
        else
            axT = 0.0
            axP = 0.0
        endif
        call AddOutVar(IOVars,"rotAxisT",axT)  ! Radians
        call AddOutVar(IOVars,"rotAxisP",axP)  ! Radians

        ! Add State variables
        call AddOutVar(IOVars, "dtk", State%dtk, uStr="s")
        call AddOutVar(IOVars, "nStepk", State%nStepk*1.0_rp, uStr="#", dStr="Number of steps each channel has taken")
        call AddOutVar(IOVars, "bminZ",State%Bmin(is:ie+1,js:je+1,ZDIR),uStr="nT")

        ! Calc intensity
        allocate(outTmp3D(is:ie,js:je,Grid%Nk))
        outTmp3D = 0.0
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,s)
        do j=js,je
            do i=is,ie
                if (State%active(i,j) .eq. RAIJUINACTIVE) then
                    cycle
                endif
                do s=1,Grid%nSpc
                    if (Grid%spc(s)%flav==F_PSPH) then
                        cycle  ! Skip plasmasphere since it has zero energy
                    endif
                    
                    outTmp3D(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd) = &
                        eta2intensity(Grid%spc(s),   &
                        State%bvol_cc(i,j),   &
                        State%eta_avg(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd))
                enddo
            enddo
        enddo
        call AddOutVar(IOVars,"intensity",outTmp3D(is:ie,js:je, :),uStr="1/(s*sr*keV*cm^2)")
        deallocate(outTmp3D)

        ! Coupling things
        call AddOutVar(IOVars,"dtCpl"   ,State%dt)  ! Attr
        call AddOutVar(IOVars,"xmin"    ,State%xyzMin (is:ie+1,js:je+1,XDIR)  ,uStr="Rx",dStr="(corners) X location of bmin surface")
        call AddOutVar(IOVars,"ymin"    ,State%xyzMin (is:ie+1,js:je+1,YDIR)  ,uStr="Rx",dStr="(corners) Y location of bmin surface")
        call AddOutVar(IOVars,"zmin"    ,State%xyzMin (is:ie+1,js:je+1,ZDIR)  ,uStr="Rx",dStr="(corners) Z location of bmin surface")
        call AddOutVar(IOVars,"radMax"   ,State%radMax    (is:ie+1,js:je+1)   ,uStr="Rx",dStr="(corners) max radius of field line")
        call AddOutVar(IOVars,"radMax_cc",State%radMax_cc (is:ie  ,js:je  )   ,uStr="Rx",dStr="max radius of field line")
        call AddOutVar(IOVars,"rCurv"   ,State%rCurv  (is:ie,js:je)           ,uStr="Rx",dStr="Min curvature radius")
        call AddOutVar(IOVars,"topo"    ,State%topo   (is:ie+1,js:je+1)*1.0_rp,uStr="0=Open, 1=Closed",dStr="(corners) Magnetic topology")
        call AddOutVar(IOVars,"active"  ,State%active (is:ie,js:je)*1.0_rp    ,uStr="-1=Inactive, 0=Buffer, 1=Active")
        call AddOutVar(IOVars,"OCBDist" ,State%OCBDist(is:ie,js:je)*1.0_rp    ,dStr="Cell distance from an open closed boundary")
        call AddOutVar(IOVars,"bVol"    ,State%bvol   (is:ie+1,js:je+1)       ,uStr="Rx/nT",dStr="(corners) Flux Tube Volume")
        call AddOutVar(IOVars,"vaFrac"  ,State%vaFrac (is:ie+1,js:je+1)       ,uStr="fraction",dStr="Fraction of Alfven speed over magnetofast + ExB speed")
        call AddOutVar(IOVars,"Pavg_in" ,State%Pavg   (is:ie,js:je, :)        ,uStr="nPa" ,dStr="Pressures from imagtubes")
        call AddOutVar(IOVars,"Davg_in" ,State%Davg   (is:ie,js:je, :)        ,uStr="#/cc",dStr="Densities from imagtubes")
        call AddOutVar(IOVars,"Pstd_in" ,State%Pstd   (is:ie,js:je, :)        ,uStr="normalized" ,dStr="Std. dev. of species pressure from imagtubes")
        call AddOutVar(IOVars,"Dstd_in" ,State%Dstd   (is:ie,js:je, :)        ,uStr="normalized" ,dStr="Std. dev. of species density from imagtubes")
        call AddOutVar(IOVars,"domWeights",State%domWeights(is:ie,js:je)     ,                   dStr="Weights used for smoothing moments for questionable cells along boundary")
        call AddOutSGV(IOVars,"Tbounce",State%Tb, outBndsO=outBnds2D_corner, uStr="[s]", dStr="Bounce timescale along field line (Alfven crossing time)", doWriteMaskO=.false.)
        call AddOutVar(IOVars,"bndLoc"    ,State%bndLoc    (js:je)*1.0_rp    ,dStr="Inactive cell closest to planet")
        call AddOutVar(IOVars,"activeLoc" ,State%activeLoc (js:je)*1.0_rp    ,dStr="Active cell furthest from planet")
        if (Model%doOwnCorot) then
            ! If we are addng our own corotation potential, then its not included in espot
            call AddOutVar(IOVars,"pot_iono",State%espot(is:ie+1,js:je+1)                                  ,uStr="kV",dStr="(corners) Ionospheric potential")
        else
            ! Otherwise, its in espot, so substract it out
            call AddOutVar(IOVars,"pot_iono",State%espot(is:ie+1,js:je+1)-State%pot_corot(is:ie+1,js:je+1) ,uStr="kV",dStr="(corners) Ionospheric potential")
        endif
        call AddOutVar(IOVars,"pot_corot",State%pot_corot(is:ie+1,js:je+1),uStr="kV",dStr="(corners) Corotation potential")
        call AddOutVar(IOVars,"pot_eb",State%pot_eb(is:ie+1,js:je+1),uStr="V",dStr="(corners) Smoothed electric potential used for RAIJU drift")
        call AddOutVar(IOVars,"pot_ftv",State%pot_ftv(is:ie+1,js:je+1),dStr="(corners) Smoothed flux-tube-volume potential piece used for RAIJU drift")

        ! Idk about you but I did not expect true to equal -1
        allocate(outTmp2D(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%Nk))
        where (State%activeShells)
            outTmp2D = 1.0
        elsewhere
            outTmp2D = 0.0
        end where
        call AddOutVar(IOVars,"activeShells",outTmp2D,uStr="[Ni, Nk]")
        deallocate(outTmp2D)

        ! Moments
        call AddOutSGV(IOVars, "Pressure", State%Press_avg, outBndsO=outBnds2D, uStr="nPa" , dStr="Step-averaged Pressure from RAIJU flavors", doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Density" , State%Den_avg  , outBndsO=outBnds2D, uStr="#/cc", dStr="Step-averaged Density from RAIJU flavors" , doWriteMaskO=.false.)
        

        ! Calculate flux tube entropy using bulk pressure
        allocate(outTmp2D(is:ie,js:je))
        outTmp2D = -1.0
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j)
        do j=js,je
            do i=is,ie
                if (State%active(i,j) .eq. RAIJUINACTIVE) then
                    cycle
                endif
                outTmp2D(i,j) = State%Press(0)%data(i,j)*State%bvol_cc(i,j)**(5./3.)
            enddo
        enddo
        call AddOutVar(IOVars,"FTEntropy",outTmp2D,uStr="nPa*(Rp/nT)^(5/3)")
        deallocate(outTmp2D)



        !--- Losses ---!
        if (Model%doLosses) then
            ! (Ni, Nj, Nk)
            call AddOutVar(IOVars, "lossRate", State%lossRates(is:ie,js:je,:), uStr="1/s")
            !call AddOutVar(IOVars, "lossRatePrecip", State%lossRatesPrecip(is:ie,js:je,:), uStr="1/s")

            ! Calculate accumulated precipitation for each species
            allocate(outPrecipN   (is:ie,js:je,Grid%nSpc))
            allocate(outPrecipE   (is:ie,js:je,Grid%nSpc))
            allocate(outCCHeatFlux(is:ie,js:je,Grid%nSpc))
            outPrecipN    = 0.0_rp
            outPrecipE    = 0.0_rp
            outCCHeatFlux = 0.0_rp
            do s=1,Grid%nSpc
                ks = Grid%spc(s)%kStart
                ke = Grid%spc(s)%kEnd
                do k=ks,ke !1,Grid%Nk
                    outPrecipN(:,:,s) = outPrecipN(:,:,s) + State%precipNFlux(k)%data(is:ie,js:je)
                    outPrecipE(:,:,s) = outPrecipE(:,:,s) + State%precipEFlux(k)%data(is:ie,js:je)
                enddo
                outCCHeatFlux(:,:,s) = sum(State%CCHeatFlux (is:ie,js:je,kS:kE), dim=3)

            enddo
            call AddOutVar(IOVars, "precipNFlux", outPrecipN   , uStr="#/cm^2/s")
            call AddOutVar(IOVars, "precipEFlux", outPrecipE   , uStr="erg/cm^2/s")
            call AddOutVar(IOVars, "CCHeatFlux" , outCCHeatFlux, uStr="eV/cm^2/s")
            deallocate(outPrecipN)
            deallocate(outPrecipE)
            deallocate(outCCHeatFlux)

            allocate(outDeltaPressure(is:ie,js:je,Grid%nSpc))
            outDeltaPressure = 0.0_rp
            do j=js,je
                do i=is,ie
                    do s=1,Grid%nSpc
                        outDeltaPressure(i,j,s) = outDeltaPressure(i,j,s) + SpcEta2Press(grid%spc(s), State%dEta_dt(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd), State%bVol_cc(i,j))
                    enddo
                enddo
            enddo
            call AddOutVar(IOVars, "dPressure_dt" , outDeltaPressure, uStr="nPa/s")
            deallocate(outDeltaPressure)

            if (Model%doOutput_3DLoss .or. Model%doOutput_debug) then
                call AddOutVar(IOVars, "dEta_dt" , State%dEta_dt(is:ie,js:je,:), uStr="eta_units/s")
                call AddOutVar(IOVars, "CCHeatFlux_Nk" , State%CCHeatFlux (is:ie,js:je,:), uStr="eV/cm^2/s")
                call AddOutSGV(IOVars, "precipNFlux_Nk", State%precipNFlux, outBndsO=outBnds2D, uStr="#/cm^2/s" , dStr="precipNFlux from RAIJU flavors", doWriteMaskO=.false.)
                call AddOutSGV(IOVars, "precipEFlux_Nk", State%precipEFlux, outBndsO=outBnds2D, uStr="erg/cm^2/s" , dStr="precipEFlux from RAIJU flavors", doWriteMaskO=.false.)
            endif
        endif

        if (Model%doOutput_velocities .or. Model%doOutput_debug) then
            call AddOutVar(IOVars, "iVel", State%iVel(is:ie+1,js:je+1,:,:), uStr="m/s")
            call AddOutVar(IOVars, "cVel", State%cVel(is:ie  ,js:je  ,:,:), uStr="m/s")
        endif
        
        
        ! (Ni, Nj, 2)
        if (Model%doOutput_debug) then
            ! New members, being tucked away into debug only
            call AddOutVar(IOVars,"bminX",State%Bmin(is:ie+1,js:je+1,XDIR),uStr="nT")
            call AddOutVar(IOVars,"bminY",State%Bmin(is:ie+1,js:je+1,YDIR),uStr="nT")
            call AddOutVar(IOVars,"bVol_cc",State%bvol_cc(is:ie,js:je)           ,uStr="Rx/nT",dStr="(corners) Flux Tube Volume")
            call AddOutVar(IOVars,"colatc" ,State%thcon  (is:ie+1,js:je+1)       ,uStr="radians",dStr="(corners) Congugate latitude")
            call AddOutVar(IOVars,"lonc"   ,State%phcon  (is:ie+1,js:je+1)       ,uStr="radians",dStr="(corners) Congugate longitude")
            call AddOutVar(IOVars,"tiote"  ,State%tiote  (is:ie,js:je)           ,uStr="normalized" ,dStr="Ratio of ion temperature over electron temperature")
            call AddOutVar(IOVars,"eta"    ,State%eta    (is:ie,js:je, :),uStr="#/cm^3 * Rx/T") !! TODO: Maybe swap with just intensity instead
            call AddOutVar(IOVars,"eta_avg",State%eta_avg(is:ie,js:je, :),uStr="#/cm^3 * Rx/T") !! TODO: Maybe swap with just intensity instead
            call AddOutSGV(IOVars,"Pressure_inst", State%Press, outBndsO=outBnds2D, uStr="nPa" , dStr="Instantaneous pressure from RAIJU flavors", doWriteMaskO=.false.)
            call AddOutSGV(IOVars,"Density_inst" , State%Den  , outBndsO=outBnds2D, uStr="#/cc", dStr="Instantaneous density from RAIJU flavors" , doWriteMaskO=.false.)

            call AddOutVar(IOVars, "eta_half"     , State%eta_half     (is:ie  ,js:je  ,:)  , uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFaceReconL", State%etaFaceReconL(is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFaceReconR", State%etaFaceReconR(is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFacePDML"  , State%etaFacePDML  (is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFacePDMR"  , State%etaFacePDMR  (is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFlux"      , State%etaFlux      (is:ie+1,js:je+1,:,:), uStr="(#/cm^3 * Rx/T) * m^2 / s")
            
            ! Mapping diagnostic
            allocate(outTmp2D(Grid%shGrid%isg:Grid%shGrid%ieg  ,Grid%shGrid%jsg:Grid%shGrid%jeg))
            call calcMapJacNorm(Grid, State%xyzMin, outTmp2D)
            call AddOutVar(IOVars, "mapJacNorm", outTmp2D(is:ie,js:je), dStr="L_(2,1) norm of lat/lon => xyzMin Jacobian")
            deallocate(outTmp2D)
        endif

        ! Species/flavor info
        call WriteRaijuSpeciesStepInfo(Model, Grid, State, IOVars)

        !Performance Metrics
        nClkSteps = readNCalls('DeepUpdate')
        call AddOutVar(IOVars, "_perf_stepTime", readClock(1)/nClkSteps)
        call AddOutVar(IOVars, "_perf_preAdvance", readClock("Pre-Advance")/nClkSteps)
        call AddOutVar(IOVars, "_perf_advanceState", readClock("AdvanceState")/nClkSteps)
        call AddOutVar(IOVars, "_perf_moments", readClock("Moments Eval")/nClkSteps)

        call WriteVars(IOVars,.true.,Model%raijuH5, gStr)

        ! Any extra groups to add
        if (Model%doLosses .and. Model%doOutput_3DLoss) then
            do i=1,size(State%lps)
                call State%lps(i)%p%doOutput(Model, Grid, State, gStr, doGhostso=doGhosts)
            enddo
        endif

    end subroutine WriteRaiju


    subroutine WriteRaijuSpeciesStepInfo(Model, Grid, State, IOVars)
        !! Calculates and writes some higher-level species/flavor info to help with postprocessing and analysis
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        type(IOVAR_T), dimension(MAXIOVAR), intent(inout) :: IOVars

        integer :: s, sIdx
        integer, dimension(2) :: mpLoc
        real(rp), dimension(Grid%nSpc) :: maxPress, maxDen

        ! Max Density and Pressure
        do s=1, Model%nSpc
            sIdx = spcIdx(Grid, Grid%spc(s)%flav)
            if (Grid%spc(s)%flav == F_PSPH) then
                mPLoc = maxloc(State%Den(sIdx)%data, mask=State%active==RAIJUACTIVE)  ! First element = 1, need to account for index offset after
            else
                mPLoc = maxloc(State%Press(sIdx)%data, mask=State%active==RAIJUACTIVE)  ! First element = 1, need to account for index offset after
            endif
            mpLoc(1) = mpLoc(1) + Grid%shGrid%isg - 1
            mpLoc(2) = mpLoc(2) + Grid%shGrid%jsg - 1

            maxPress(s) = State%Press(sIdx)%data(mPLoc(1), mpLoc(2))
            maxDen(s)   = State%Den  (sIdx)%data(mPLoc(1), mpLoc(2))!/Grid%spc(sIdx)%amu
        enddo

        call AddOutVar(IOVars, "spc_maxPress", maxPress, uStr="[nPa]", dStr="Max pressure for each species/flavor")
        call AddOutVar(IOVars, "spc_maxDen"  , maxDen  , uStr="[#/cc]", dStr="Max density for each species/flavor")

    end subroutine WriteRaijuSpeciesStepInfo


    subroutine WriteRaijuColdStartInfo(State, IOVars)
        type(raijuState_T), intent(in) :: State
        type(IOVAR_T), dimension(MAXIOVAR), intent(inout) :: IOVars

        call AddOutVar(IOVars, "cs_doneFirstCS", merge(1,0,State%coldStarter%doneFirstCS))
        call AddOutVar(IOVars, "cs_doCS_next_preAdv", merge(1,0,State%coldStarter%doCS_next_preAdv))
        call AddOutVar(IOVars, "cs_modelDst_next_preAdv", State%coldStarter%modelDst_next_preAdv)
        call AddOutVar(IOVars, "cs_lastEval", State%coldStarter%lastEval)
        call AddOutVar(IOVars, "cs_lastTarget", State%coldStarter%lastTarget)

    end subroutine WriteRaijuColdStartInfo


    subroutine ReadRaijuColdStartInfo(State, IOVars)
        type(raijuState_T), intent(inout) :: State
        type(IOVAR_T), dimension(MAXIOVAR), intent(inout) :: IOVars
        integer :: tmpInt

        State%coldStarter%lastEval = GetIOReal(IOVars, "cs_lastEval")
        State%coldStarter%lastTarget = GetIOReal(IOVars, "cs_lastTarget")
        State%coldStarter%modelDst_next_preAdv = GetIOReal(IOVars, "cs_modelDst_next_preAdv")

        tmpInt = GetIOInt(IOVars, "cs_doneFirstCS")
        State%coldStarter%doneFirstCS = tmpInt .eq. 1
        tmpInt = GetIOInt(IOVars, "cs_doCS_next_preAdv")
        State%coldStarter%doCS_next_preAdv = tmpInt .eq. 1

    end subroutine ReadRaijuColdStartInfo


    subroutine WriteRaijuPotentialInfo(Model, State, IOVars, is, ie, js, je)
        type(raijuModel_T), intent(in) :: Model
        type(raijuState_T), intent(in) :: State
        type(IOVAR_T), dimension(MAXIOVAR), intent(inout) :: IOVars
        integer, intent(in) :: is, ie, js, je
        ! If doOwnCorot=True, then espot will just be ionosphere E-field
        ! If doOwnCorot=False, we expect espot to already have iono + corot in it, so to get the iono contribution we need to take out corot
        if (Model%doOwnCorot) then
            call AddOutVar(IOVars,"pot_iono",State%espot(is:ie+1,js:je+1),uStr="kV",dStr="(corners) Ionospheric potential")
        else
            call AddOutVar(IOVars,"pot_iono",State%espot(is:ie+1,js:je+1)-State%pot_corot(is:ie+1,js:je+1),uStr="kV",dStr="(corners) Ionospheric potential")
        endif
        call AddOutVar(IOVars,"pot_corot",State%pot_corot(is:ie+1,js:je+1),uStr="kV",dStr="(corners) Corotation potential")

    end subroutine WriteRaijuPotentialInfo


    subroutine ReadRaijuPotentialInfo(Model, State, IOVars)
        type(raijuModel_T), intent(in) :: Model
        type(raijuState_T), intent(inout) :: State
        type(IOVAR_T), dimension(MAXIOVAR), intent(inout) :: IOVars

        real(rp), dimension(:,:), allocatable :: pot_iono

        allocate(pot_iono(lbound(State%espot,1):ubound(State%espot,1), &
                          lbound(State%espot,2):ubound(State%espot,2)))

        call IOArray2DFill(IOVars, "pot_iono", pot_iono)
        call IOArray2DFill(IOVars, "pot_corot", State%pot_corot(:,:))

        if (Model%doOwnCorot) then
            State%espot = pot_iono
        else
            State%espot = pot_iono + State%pot_corot
        endif

        deallocate(pot_iono)

    end subroutine ReadRaijuPotentialInfo


    subroutine WriteRaijuRes(Model, Grid, State, opt, ResF)
        !! Writes RAIJU restart info to provided path ResF
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        type(raijuOptions_T), intent(in) :: opt
        character(len=strLen), intent(in) :: ResF

        ! If a restart already exists, get rid of old one
        call CheckAndKill(ResF)

        call raijuWriteRootH5Data(Model, Grid, ResF, doGhostsO=.true.)
        ! Save ShellGrid to root of file (Not needed now, done in raijuWriteRootH5Data and saved to Grid/ShellGrid)
        !call writeShellGrid(Grid%shGrid, ResF)
        ! And species info (also done in raijuWriteRootH5Data now)
        !call writeSpeciesInfo(Model, Grid, ResF)
        ! App options
        call WriteRaijuResOpts(opt, ResF)
        ! All necessary State info
        call WriteRaijuResState(Model, Grid, State, ResF)

    end subroutine WriteRaijuRes


    subroutine WriteRaijuResOpts(opt, ResF)
        type(raijuOptions_T), intent(in) :: opt
        character(len=strLen), intent(in) :: ResF

        type(IOVAR_T), dimension(20) :: IOVars

        call AddOutVar(IOVars, "thetaL", opt%thetaL%get())
        call AddOutVar(IOVars, "thetaU", opt%thetaU%get())

        call WriteVars(IOVars,.false.,ResF,"Options")

    end subroutine

    subroutine WriteRaijuResState(Model, Grid, State, ResF)
        !! Writes RAIJU State restart info to provided path ResF
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: ResF

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        integer :: i,j,k,s
        integer :: is, ie, js, je
        character(len=strLen) :: tmpStr
        real(rp), dimension(:,:), allocatable :: outActiveShell
        real(rp), dimension(:,:), allocatable :: outFTEntropy
        real(rp), dimension(Model%nSpc) :: dpsdst_spc

        ! Always do ghosts
        is = Grid%shGrid%isg
        ie = Grid%shGrid%ieg
        js = Grid%shGrid%jsg
        je = Grid%shGrid%jeg

        !Reset IO chain
        call ClearIO(IOVars)

        call AddOutVar(IOVars,"time",State%t)
        call AddOutVar(IOVars,"ts"  ,State%ts)
        call AddOutVar(IOVars,"MJD" ,State%mjd)
        if (State%isFirstCpl) then
            call AddOutVar(IOVars,"isFirstCpl", 1)
        else
            call AddOutVar(IOVars,"isFirstCpl", 0)
        endif
        ! IOClock data
        call AddOutVar(IOVars,"nOut" ,State%IO%nOut)
        call AddOutVar(IOVars,"nRes" ,State%IO%nRes)
        call AddOutVar(IOVars,"tRes" ,State%IO%tRes)
        call WriteRaijuColdStartInfo(State, IOVars)

        do s=1,Model%nSpc
            dpsdst_spc(s) = spcEta2DPS(Model, Grid, State%bvol_cc, State%eta, Grid%spc(s), State%active .eq. RAIJUACTIVE)
            write(tmpStr, '(A,I0)')"DPSDst_",Grid%spc(s)%flav
            call AddOutVar(IOVars, trim(tmpStr), dpsdst_spc(s))
        enddo
        call AddOutVar(IOVars,"DPSDst_total",sum(dpsdst_spc))

        ! Add State variables

        ! Coupling things
        call AddOutVar(IOVars,"dtCpl"  ,State%dt)  ! Attr
        call AddOutVar(IOVars,"xmin"   ,State%xyzMin (:,:,XDIR)  ,uStr="Rx",dStr="(corners) X location of bmin surface")
        call AddOutVar(IOVars,"ymin"   ,State%xyzMin (:,:,YDIR)  ,uStr="Rx",dStr="(corners) Y location of bmin surface")
        call AddOutVar(IOVars,"zmin"   ,State%xyzMin (:,:,ZDIR)  ,uStr="Rx",dStr="(corners) Z location of bmin surface")
        call AddOutVar(IOVars,"radMax"   ,State%radMax    (is:ie+1,js:je+1)   ,uStr="Rx",dStr="(corners) max radius of field line")
        call AddOutVar(IOVars,"radMax_cc",State%radMax_cc (is:ie  ,js:je  )   ,uStr="Rx",dStr="max radius of field line")
        call AddOutVar(IOVars,"bminX"  ,State%Bmin   (:,:,XDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"bminY"  ,State%Bmin   (:,:,YDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"bminZ"  ,State%Bmin   (:,:,ZDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"topo"   ,State%topo   (:,:)*1.0_rp,uStr="0=Open, 1=Closed",dStr="(corners) Magnetic topology")
        call AddOutVar(IOVars,"colatc" ,State%thcon  (:,:)       ,uStr="radians",dStr="(corners) Congugate latitude")
        call AddOutVar(IOVars,"lonc"   ,State%phcon  (:,:)       ,uStr="radians",dStr="(corners) Congugate longitude")
        call WriteRaijuPotentialInfo(Model, State, IOVars, is, ie, js, je)
        call AddOutVar(IOVars,"bVol"   ,State%bvol   (:,:)       ,uStr="Rx/nT",dStr="(corners) Flux Tube Volume")
        call AddOutVar(IOVars,"bVol_cc",State%bvol_cc(:,:)       ,uStr="Rx/nT",dStr="(centers) Flux Tube Volume")
        call AddOutVar(IOVars,"rCurv"  ,State%rCurv   (:,:)      ,uStr="Rx",dStr="Min curvature radius")
        call AddOutVar(IOVars,"vaFrac" ,State%vaFrac (:,:)       ,uStr="fraction",dStr="Fraction of Alfven speed over magnetofast + ExB speed")

        ! If only 1 element in 3rd position, AddOutVar will write as 2D array
        ! We know it should be 3D, so force it
        call AddOutVar(IOVars,"Pavg_in",State%Pavg   (:,:,:)     ,uStr="nPa" ,dStr="Pressures from imagtubes", doSqzO=.false.)
        call AddOutVar(IOVars,"Davg_in",State%Davg   (:,:,:)     ,uStr="#/cc",dStr="Densities from imagtubes", doSqzO=.false.)
        call AddOutVar(IOVars,"Pstd_in",State%Pstd   (:,:,:)     ,uStr="normalized",dStr="Std. dev. of species pressure from imagtubes", doSqzO=.false.)
        call AddOutVar(IOVars,"Dstd_in",State%Dstd   (:,:,:)     ,uStr="normalized",dStr="Std. dev. of species density from imagtubes", doSqzO=.false.)
        call AddOutVar(IOVars,"domWeights",State%domWeights(:,:)  ,dStr="Weights used for smoothing moments for questionable cells along boundary")
        call AddOutVar(IOVars,"bndLoc"    ,State%bndLoc(:)*1.0_rp ,dStr="Inactive cell closest to planet")
        call AddOutVar(IOVars,"activeLoc" ,State%activeLoc(:)*1.0_rp,dStr="Active cell furthest from planet")

        ! Core variables
        call AddOutVar(IOVars,"eta"        ,State%eta         (:,:,:)     ,uStr="#/cm^3 * Rx/T")
        call AddOutVar(IOVars,"eta_last"   ,State%eta_last    (:,:,:)     ,uStr="#/cm^3 * Rx/T")
        call AddOutVar(IOVars,"eta_avg"    ,State%eta_avg     (:,:,:)     ,uStr="#/cm^3 * Rx/T")
        call AddOutVar(IOVars,"active"     ,State%active      (:,:)*1.0_rp,uStr="-1=Inactive, 0=Buffer, 1=Active")
        call AddOutVar(IOVars,"active_last",State%active_last (:,:)*1.0_rp,uStr="-1=Inactive, 0=Buffer, 1=Active")
        call AddOutVar(IOVars,"OCBDist"    ,State%OCBDist     (:,:)*1.0_rp,dStr="Cell distance from an open closed boundary")
        ! Note: no need to convert to more intuitive 0,1 here cause we're just gonna un-do it on readResState
        allocate(outActiveShell(is:ie, Grid%Nk))
        where (State%activeShells)
            outActiveShell = 1.0
        elsewhere
            outActiveShell = 0.0
        end where
        call AddOutVar(IOVars,"activeShells",outActiveShell,uStr="[Ni, Nk]")
        ! Moments
        call AddOutSGV(IOVars, "Pressure", State%Press_avg, doWriteMaskO=.true., uStr="nPa")
        call AddOutSGV(IOVars, "Density" , State%Den_avg  , doWriteMaskO=.true., uStr="#/cc")
        call AddOutSGV(IOVars, "Pressure_inst", State%Press, doWriteMaskO=.true., uStr="nPa")
        call AddOutSGV(IOVars, "Density_inst" , State%Den  , doWriteMaskO=.true., uStr="#/cc")
        call AddOutSGV(IOVars, "Tbounce", State%Tb, doWriteMaskO=.true., uStr="[s]")
        ! Precip
        call AddOutSGV(IOVars, "precipNFlux_Nk", State%precipNFlux, doWriteMaskO=.true., uStr="#/cm^2/s")
        call AddOutSGV(IOVars, "precipEFlux_Nk", State%precipEFlux, doWriteMaskO=.true., uStr="erg/cm^2/s")
        call AddOutVar(IOVars, "lossRate"      , State%lossRates      (:,:,:), uStr="1/s")
        call AddOutVar(IOVars, "lossRatePrecip", State%lossRatesPrecip(:,:,:), uStr="1/s")
        ! More solver stuff
        call AddOutVar(IOVars, "dtk"   , State%dtk(:), uStr="s")
        call AddOutVar(IOVars, "nStepk", State%nStepk(:)*1.0_rp, uStr="#", dStr="Number of steps each channel has taken")
        call AddOutVar(IOVars, "iVel"   , State%iVel  (:,:,:,:)     , uStr="m/s")
        call AddOutVar(IOVars, "cVel"   , State%cVel  (:,:,:,:), uStr="m/s")

        allocate(outFTEntropy(is:ie,js:je))
        outFTEntropy = -1.0_rp
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j)
        do j=js,je
            do i=is,ie
                if (State%active(i,j) .eq. RAIJUINACTIVE) then
                    cycle
                endif
                outFTEntropy(i,j) = State%Press(0)%data(i,j)*State%bvol_cc(i,j)**(5./3.)
            enddo
        enddo
        call AddOutVar(IOVars,"FTEntropy",outFTEntropy,uStr="nPa*(Rp/nT)^(5/3)")
        deallocate(outFTEntropy)

        call WriteRaijuSpeciesStepInfo(Model, Grid, State, IOVars)

        call WriteVars(IOVars,.false.,ResF,"State")

    end subroutine WriteRaijuResState

    
    subroutine ReadRaijuResState(Model, Grid, State, inH5)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        character(len=*), intent(in) :: inH5

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        integer :: Ntg, Npg
            !! Number of theta and phi cells, including ghosts
        integer :: is,ie,js,je, s
        integer :: tmpInt
        real(rp), dimension(:)  , allocatable :: tmpReal1D
        real(rp), dimension(:,:), allocatable :: tmpReal2D

        is = Grid%shGrid%isg
        ie = Grid%shGrid%ieg
        js = Grid%shGrid%jsg
        je = Grid%shGrid%jeg

        !Reset IO chain
        call ClearIO(IOVars)

        call AddInVar(IOVars,"time")
        call AddInVar(IOVars,"ts"  )
        call AddInVar(IOVars,"MJD" )
        call AddInVar(IOVars,"nOut")
        call AddInVar(IOVars,"nRes")
        call AddInVar(IOVars,"tRes")
        call AddInVar(IOVars,"isFirstCpl")
        call AddInVar(IOVars,"cs_doCS_next_preAdv")
        call AddInVar(IOVars,"cs_modelDst_next_preAdv")
        call AddInVar(IOVars,"cs_doneFirstCS")
        call AddInVar(IOVars,"cs_lastEval")
        call AddInVar(IOVars,"cs_lastTarget")

        call AddInVar(IOVars,"dtCpl"  )
        call AddInVar(IOVars,"xmin"   )
        call AddInVar(IOVars,"ymin"   )
        call AddInVar(IOVars,"zmin"   )
        call AddInVar(IOVars,"radMax" )
        call AddInVar(IOVars,"radMax_cc")
        call AddInVar(IOVars,"bminX"  )
        call AddInVar(IOVars,"bminY"  )
        call AddInVar(IOVars,"bminZ"  )
        call AddInVar(IOVars,"topo"   )
        call AddInVar(IOVars,"colatc" )
        call AddInVar(IOVars,"lonc"   )
        call AddInVar(IOVars,"pot_iono")
        call AddInVar(IOVars,"pot_corot")
        call AddInVar(IOVars,"bVol"   )
        call AddInVar(IOVars,"bVol_cc")
        call AddInVar(IOVars,"rCurv"  )
        call AddInVar(IOVars,"vaFrac" )
        call AddInVar(IOVars,"Pavg_in")
        call AddInVar(IOVars,"Davg_in")
        call AddInVar(IOVars,"Pstd_in")
        call AddInVar(IOVars,"Dstd_in")
        call AddInVar(IOVars,"domWeights")
        call AddInVar(IOVars,"bndLoc")
        call AddInVar(IOVars,"activeLoc")

        call AddInVar(IOVars,"eta"        )
        call AddInVar(IOVars,"eta_last"   )
        call AddInVar(IOVars,"eta_avg"    )
        call AddInVar(IOVars,"active"     )
        call AddInVar(IOVars,"active_last")
        call AddInVar(IOVars,"OCBDist"    )
        call AddInVar(IOVars,"activeShells")

        call AddInVar(IOVars,"lossRate")
        call AddInVar(IOVars,"lossRatePrecip")

        call AddInVar(IOVars, "dtk"    )
        call AddInVar(IOVars, "nStepk" )
        call AddInVar(IOVars, "iVel"   )
        call AddInVar(IOVars, "cVel")

        call ReadVars(IOVars,.false.,inH5,"State")

        State%t       = GetIOReal(IOVars, "time" )
        State%ts      = GetIOInt (IOVars, "ts"   )
        State%mjd     = GetIOReal(IOVars, "MJD"  )
        State%dt      = GetIOReal(IOVars, "dtCpl")
        State%IO%nOut = GetIOInt (IOVars, "nOut" )
        State%IO%nRes = GetIOInt (IOVars, "nRes" )
        State%IO%tRes = GetIOReal(IOVars, "tRes" )

        call ReadRaijuColdStartInfo(State, IOVars)
        
        ! Handle boolean attributes
        tmpInt = GetIOInt(IOVars, "isFirstCpl")
        State%isFirstCpl = tmpInt .eq. 1
        call IOArray2DFill(IOVars, "xmin" , State%xyzMin(:,:,XDIR))
        call IOArray2DFill(IOVars, "ymin" , State%xyzMin(:,:,YDIR))
        call IOArray2DFill(IOVars, "zmin" , State%xyzMin(:,:,ZDIR))
        call IOArray2DFill(IOVars, "radMax"   , State%radMax(:,:))
        call IOArray2DFill(IOVars, "radMax_cc", State%radMax_cc(:,:))
        call IOArray2DFill(IOVars, "bminX", State%Bmin  (:,:,XDIR))
        call IOArray2DFill(IOVars, "bminY", State%Bmin  (:,:,YDIR))
        call IOArray2DFill(IOVars, "bminZ", State%Bmin  (:,:,ZDIR))
        call IOArray2DFill(IOVars, "colatc", State%thcon  (:,:))
        call IOArray2DFill(IOVars, "lonc"  , State%phcon  (:,:))
        call ReadRaijuPotentialInfo(Model, State, IOVars)
        call IOArray2DFill(IOVars, "bVol"  , State%bvol   (:,:))
        call IOArray2DFill(IOVars, "bVol_cc", State%bvol_cc(:,:))
        call IOArray2DFill(IOVars, "rCurv" , State%rCurv   (:,:))
        call IOArray2DFill(IOVars, "vaFrac", State%vaFrac (:,:))
        call IOArray3DFill(IOVars, "Pavg_in", State%Pavg(:,:,:))
        call IOArray3DFill(IOVars, "Davg_in", State%Davg(:,:,:))
        call IOArray3DFill(IOVars, "Pstd_in", State%Pstd(:,:,:))
        call IOArray3DFill(IOVars, "Dstd_in", State%Dstd(:,:,:))
        call IOArray2DFill(IOVars, "domWeights", State%domWeights(:,:))
        

        call IOArray3DFill(IOVars, "eta"     , State%eta     (:,:,:))
        call IOArray3DFill(IOVars, "eta_last", State%eta_last(:,:,:))
        call IOArray3DFill(IOVars, "eta_avg" , State%eta_avg (:,:,:))

        call IOArray3DFill(IOVars, "lossRate", State%lossRates(:,:,:))
        call IOArray3DFill(IOVars, "lossRatePrecip", State%lossRatesPrecip(:,:,:))
        
        call IOArray1DFill(IOVars, "dtk", State%dtk(:))

        call IOArray4DFill(IOVars, "iVel", State%iVel(:,:,:,:))
        call IOArray4DFill(IOVars, "cVel", State%cVel(:,:,:,:))

        ! Handle real -> int arrays
        associate(sh=>Grid%shGrid)
        Ntg = sh%Nt+sh%Ngn+sh%Ngs
        Npg = sh%Np+sh%Nge+sh%Ngw
        ! Cell corner integer variables
        allocate(tmpReal2D(Ntg+1,Npg+1))
        call IOArray2DFill(IOVars, "topo",tmpReal2D)
        State%topo  (:,:)      = INT(tmpReal2D)
        ! Cell center integer variables
        deallocate(tmpReal2D)
        allocate(tmpReal2D(Ntg, Npg))
        call IOArray2DFill(IOVars, "active",tmpReal2D)
        State%active(:,:)      = INT(tmpReal2D)
        call IOArray2DFill(IOVars, "active_last",tmpReal2D)
        State%active_last(:,:)      = INT(tmpReal2D)
        call IOArray2DFill(IOVars, "OCBDist",tmpReal2D)
        State%OCBDist(:,:)      = INT(tmpReal2D)
        ! Weird [Ni, Nk] activeShells
        deallocate(tmpReal2D)
        allocate(tmpReal2D(Ntg, Grid%Nk))
        call IOArray2DFill(IOVars, "activeShells",tmpReal2D)
        State%activeShells = merge(.true., .false., tmpReal2D .eq. 1.0)

        ! 1D Nk
        allocate(tmpReal1D(lbound(State%bndLoc,1):ubound(State%bndLoc,1)))
        call IOArray1DFill(IOVars, "bndLoc", tmpReal1D)
        State%bndLoc = INT(tmpReal1D)
        call IOArray1DFill(IOVars, "activeLoc", tmpReal1D)
        State%activeLoc = INT(tmpReal1D)
        deallocate(tmpReal1D)

        allocate(tmpReal1D(Grid%Nk))
        call IOArray1DFill(IOVars, "nStepk", tmpReal1D)
        State%nStepk = INT(tmpReal1D)

        ! ShellGridVars
        call ReadInSGV(State%Press_avg, inH5, "Pressure", "State", doIOpO=.false.)
        call ReadInSGV(State%Den_avg  , inH5, "Density" , "State", doIOpO=.false.)
        call ReadInSGV(State%Press, inH5, "Pressure_inst", "State", doIOpO=.false.)
        call ReadInSGV(State%Den  , inH5, "Density_inst" , "State", doIOpO=.false.)
        call ReadInSGV(State%Tb   , inH5, "Tbounce", "State", doIOpO=.false.)
        call ReadInSGV(State%precipNFlux, inH5, "precipNFlux_Nk", "State", doIOpO=.false.)
        call ReadInSGV(State%precipEFlux, inH5, "precipEFlux_Nk", "State", doIOpO=.false.)

        end associate

    end subroutine ReadRaijuResState


    subroutine WriteRaijuBlackBox(raiApp)
        class(raijuApp_T), intent(in) :: raiApp

        character(len=strLen) :: outF
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        integer :: is, ie, js, je, ks, ke
        integer :: i,j,k,s, nClkSteps
        integer, dimension(4) :: outBnds2D, outBnds2D_corner
        real(rp) :: axT, axP
        real(rp), dimension(:,:), allocatable :: outTmp2D

        write(outF,"(A,A)") "CRASH",trim(raiApp%Model%raijuH5)

        write(*,*)"Writing raiju black box to ",trim(outF)

        call raijuWriteRootH5Data(raiApp%Model, raiApp%Grid, fnameO=outF, doGhostsO=.true.)

        call ClearIO(IOVars)

        associate(Model=>raiApp%Model, Grid=>raiApp%Grid, State=>raiApp%State)

        is = Grid%shGrid%isg
        ie = Grid%shGrid%ieg
        js = Grid%shGrid%jsg
        je = Grid%shGrid%jeg
        outBnds2D = (/is,ie,js,je/)
        outBnds2D_corner = (/is,ie+1,js,je+1/)

        call AddOutVar(IOVars,"time",State%t)
        call AddOutVar(IOVars,"dtCpl"  ,State%dt)
        call AddOutVar(IOVars,"MJD",State%mjd)
        call AddOutVar(IOVars,"nOut",State%IO%nOut)
        call AddOutVar(IOVars,"nRes",State%IO%nRes)
        call AddOutVar(IOVars,"isFirstCpl", merge(1,0,State%isFirstCpl))
        call WriteRaijuColdStartInfo(State, IOVars)
        ! Save axis of rotation as attributes
        if (Model%doGeoCorot) then
            call GeoAxisTP(axT, axP)
        else
            axT = 0.0
            axP = 0.0
        endif
        call AddOutVar(IOVars,"rotAxisT",axT)  ! Radians
        call AddOutVar(IOVars,"rotAxisP",axP)  ! Radians

        call AddOutVar(IOVars, "dtk", State%dtk, uStr="s")
        call AddOutVar(IOVars, "nStepk", State%nStepk*1.0_rp, uStr="#", dStr="Number of steps each channel has taken")

        call AddOutVar(IOVars,"eta"     , State%eta     (is:ie,js:je,:), uStr="#/cm^3 * Rx/T")
        call AddOutVar(IOVars,"eta_half", State%eta_half(is:ie,js:je,:), uStr="#/cm^3 * Rx/T")
        call AddOutVar(IOVars,"eta_last", State%eta_last(is:ie,js:je,:), uStr="#/cm^3 * Rx/T")
        call AddOutVar(IOVars,"eta_avg" , State%eta_avg (is:ie,js:je,:), uStr="#/cm^3 * Rx/T")

        allocate(outTmp2D(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%Nk))
        where (State%activeShells)
            outTmp2D = 1.0
        elsewhere
            outTmp2D = 0.0
        end where
        call AddOutVar(IOVars,"activeShells",outTmp2D,uStr="[Ni, Nk]")
        deallocate(outTmp2D)
        call AddOutVar(IOVars, "iVel"           , State%iVel           (is:ie+1,js:je+1,:,:), uStr="m/s")
        call AddOutVar(IOVars, "cVel"           , State%cVel           (is:ie  ,js:je  ,:,:), uStr="m/s")

        ! Coupling things
        call AddOutVar(IOVars,"Pavg_in"   ,State%Pavg      (is:ie  ,js:je  , :)    ,uStr="nPa" ,dStr="Pressures from imagtubes")
        call AddOutVar(IOVars,"Davg_in"   ,State%Davg      (is:ie  ,js:je  , :)    ,uStr="#/cc",dStr="Densities from imagtubes")
        call AddOutVar(IOVars,"Pstd_in"   ,State%Pstd      (is:ie  ,js:je  , :)    ,uStr="normalized" ,dStr="Std. dev. of species pressure from imagtubes")
        call AddOutVar(IOVars,"Dstd_in"   ,State%Dstd      (is:ie  ,js:je  , :)    ,uStr="normalized" ,dStr="Std. dev. of species density from imagtubes")
        call AddOutVar(IOVars,"tiote"     ,State%tiote     (is:ie  ,js:je  )       ,uStr="normalized" ,dStr="Ratio of ion temperature over electron temperature")
        call AddOutVar(IOVars,"bminX"     ,State%Bmin      (is:ie+1,js:je+1,XDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"bminY"     ,State%Bmin      (is:ie+1,js:je+1,YDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"bminZ"     ,State%Bmin      (is:ie+1,js:je+1,ZDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"xmin"      ,State%xyzMin    (is:ie+1,js:je+1,XDIR)  ,uStr="Rx",dStr="(corners) X location of bmin surface")
        call AddOutVar(IOVars,"ymin"      ,State%xyzMin    (is:ie+1,js:je+1,YDIR)  ,uStr="Rx",dStr="(corners) Y location of bmin surface")
        call AddOutVar(IOVars,"zmin"      ,State%xyzMin    (is:ie+1,js:je+1,ZDIR)  ,uStr="Rx",dStr="(corners) Z location of bmin surface")
        call AddOutVar(IOVars,"xmin_cc"   ,State%xyzMincc  (is:ie  ,js:je  ,XDIR)  ,uStr="Rx",dStr="(centers) X location of bmin surface")
        call AddOutVar(IOVars,"ymin_cc"   ,State%xyzMincc  (is:ie  ,js:je  ,YDIR)  ,uStr="Rx",dStr="(centers) Y location of bmin surface")
        call AddOutVar(IOVars,"zmin_cc"   ,State%xyzMincc  (is:ie  ,js:je  ,ZDIR)  ,uStr="Rx",dStr="(centers) Z location of bmin surface")
        call AddOutVar(IOVars,"radMax"    ,State%radMax    (is:ie+1,js:je+1)       ,uStr="Rx",dStr="(corners) Maximum radius of flux tube")
        call AddOutVar(IOVars,"radMax_cc" ,State%radMax_cc (is:ie  ,js:je  )       ,uStr="Rx",dStr="(centers) Maximum radius of flux tube")
        call AddOutVar(IOVars,"topo"      ,State%topo      (is:ie+1,js:je+1)*1.0_rp,uStr="0=Open, 1=Closed",dStr="(corners) Magnetic topology")
        call AddOutVar(IOVars,"colatc"    ,State%thcon     (is:ie+1,js:je+1)       ,uStr="radians",dStr="(corners) Congugate latitude")
        call AddOutVar(IOVars,"lonc"      ,State%phcon     (is:ie+1,js:je+1)       ,uStr="radians",dStr="(corners) Congugate longitude")
        call AddOutVar(IOVars,"bVol"      ,State%bvol      (is:ie+1,js:je+1)       ,uStr="Rx/nT",dStr="(corners) Flux Tube Volume")
        call AddOutVar(IOVars,"bVol_cc"   ,State%bvol_cc   (is:ie  ,js:je  )       ,uStr="Rx/nT",dStr="(centers) Flux Tube Volume")
        call AddOutVar(IOVars,"rCurv"     ,State%rCurv     (is:ie+1,js:je+1)       ,uStr="Rx",dStr="(corners) Radius of curvature")
        call AddOutVar(IOVars,"vaFrac"    ,State%vaFrac    (is:ie+1,js:je+1)       ,uStr="fraction",dStr="Fraction of Alfven speed over magnetofast + ExB speed")
        call AddOutVar(IOVars,"domWeights",State%domWeights(is:ie  ,js:je  )       ,                   dStr="Weights used for smoothing moments for questionable cells along boundary")

        call WriteRaijuPotentialInfo(Model, State, IOVars, is, ie, js, je)
        call AddOutVar(IOVars,"pot_eb",State%pot_eb(is:ie+1,js:je+1),uStr="V",dStr="(corners) Smoothed electric potential used for RAIJU drift")
        call AddOutVar(IOVars,"pot_ftv",State%pot_ftv(is:ie+1,js:je+1),dStr="(corners) Smoothed flux-tube-volume potential piece used for RAIJU drift")

        call AddOutVar(IOVars,"bndLoc"     ,State%bndLoc     (      js:je)*1.0_rp    ,dStr="Inactive cell closest to planet")
        call AddOutVar(IOVars,"activeLoc"  ,State%activeLoc  (      js:je)*1.0_rp    ,dStr="Active cell furthest from planet")
        call AddOutVar(IOVars,"active"     ,State%active     (is:ie,js:je)*1.0_rp    ,uStr="-1=Inactive, 0=Buffer, 1=Active")
        call AddOutVar(IOVars,"active_last",State%active_last(is:ie,js:je)*1.0_rp    ,uStr="-1=Inactive, 0=Buffer, 1=Active")
        call AddOutVar(IOVars,"OCBDist"    ,State%OCBDist    (is:ie,js:je)*1.0_rp    ,dStr="Cell distance from an open closed boundary")

        call AddOutVar(IOVars, "lossRate"      , State%lossRates(is:ie,js:je,:), uStr="1/s")
        call AddOutVar(IOVars, "lossRatePrecip", State%lossRatesPrecip(is:ie,js:je,:), uStr="1/s")
        call AddOutSGV(IOVars, "precipNFlux_Nk", State%precipNFlux, outBndsO=outBnds2D, uStr="#/cm^2/s" , dStr="precipNFlux from RAIJU flavors", doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "precipEFlux_Nk", State%precipEFlux, outBndsO=outBnds2D, uStr="erg/cm^2/s" , dStr="precipEFlux from RAIJU flavors", doWriteMaskO=.false.)
        call AddOutVar(IOVars, "dEta_dt"       , State%dEta_dt(is:ie,js:je,:), uStr="eta_units/s")
        call AddOutVar(IOVars, "CCHeatFlux_Nk" , State%CCHeatFlux (is:ie,js:je,:), uStr="eV/cm^2/s")

        call AddOutSGV(IOVars, "Pressure_inst", State%Press    , outBndsO=outBnds2D, uStr="nPa" , dStr="Instantaneous pressure from RAIJU flavors", doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Density_inst" , State%Den      , outBndsO=outBnds2D, uStr="#/cc", dStr="Instantaneous density from RAIJU flavors" , doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Pressure"     , State%Press_avg, outBndsO=outBnds2D, uStr="nPa" , dStr="Step-averaged Pressure from RAIJU flavors", doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Density"      , State%Den_avg  , outBndsO=outBnds2D, uStr="#/cc", dStr="Step-averaged Density from RAIJU flavors" , doWriteMaskO=.false.)
        call AddOutSGV(IOVars,"Tbounce"       , State%Tb       , outBndsO=outBnds2D_corner, uStr="[s]", dStr="Bounce timescale along field line (Alfven crossing time)", doWriteMaskO=.false.)

        if (allocated(State%etaFacePDML)) then
            call AddOutVar(IOVars, "etaFaceReconL", State%etaFaceReconL(is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFaceReconR", State%etaFaceReconR(is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFacePDML"  , State%etaFacePDML  (is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFacePDMR"  , State%etaFacePDMR  (is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFlux"      , State%etaFlux      (is:ie+1,js:je+1,:,:), uStr="(#/cm^3 * Rx/T) * m^2 / s")
        endif

        nClkSteps = readNCalls('DeepUpdate')
        call WriteRaijuSpeciesStepInfo(Model, Grid, State, IOVars)
        call AddOutVar(IOVars, "_perf_stepTime", readClock(1)/nClkSteps)
        call AddOutVar(IOVars, "_perf_preAdvance", readClock("Pre-Advance")/nClkSteps)
        call AddOutVar(IOVars, "_perf_advanceState", readClock("AdvanceState")/nClkSteps)
        call AddOutVar(IOVars, "_perf_moments", readClock("Moments Eval")/nClkSteps)

        call WriteVars(IOVars,.true.,outF,gstrO="Step#0")

        end associate

    end subroutine WriteRaijuBlackBox
end module raijuIO
