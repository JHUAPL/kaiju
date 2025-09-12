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
    ! Functions/subs we borrow for debug/diagnostic output
    use raijuPreAdvancer, only : calcEffectivePotential, calcGradVM_cc, calcVelocity_cc
    use raijuDomain, only : calcMapJacNorm

    implicit none

    integer, parameter, private :: MAXIOVAR = 70
    logical, private :: doRoot = .true. !Whether root variables need to be written
    logical, private :: doFat = .false. !Whether to output lots of extra datalogical, private :: doRoot = .true. !Whether root variables need to be written

    contains

    subroutine raijuInitIO(Model, Grid, doGhostsO)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        logical, optional, intent(in) :: doGhostsO

        integer :: is, ie, js, je
        logical :: doGhosts
        integer :: i
        logical :: fExist
        real(rp), dimension(:,:), allocatable :: lat2D, lon2D
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        character(len=strLen) :: gStr

        doRoot = .false. ! Don't call again

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

        associate(sh => Grid%shGrid, spc=>Grid%spc)
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

        ! Add spatial grid as 2D array
        allocate(lat2D(is:ie+1,js:je+1))  ! +1 because we're doing corners
        allocate(lon2D(is:ie+1,js:je+1))

        do i=js,je+1
            lat2D(:,i) = sh%th(is:ie+1)
        enddo

        do i=is,ie+1
            lon2D(i,:) = sh%ph(js:je+1)
        enddo

        ! Ready for output
        ! Some model setting info
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"doFatOutput"   ,Model%doOutput_potGrads)  ! Attr
        call AddOutVar(IOVars,"doDebugOutput" ,Model%doOutput_debug   )  ! Attr
        call AddOutVar(IOVars,"doWriteGhosts" ,Model%writeGhosts      )  ! Attr
        call AddOutVar(IOVars,"doGeoCorot"    ,Model%doGeoCorot       )  ! Attr
        call AddOutVar(IOVars,"doExcesstoPsph",Model%doExcesstoPsph   )  ! Attr
        call AddOutVar(IOVars,"doPlasmasphere",Model%doPlasmasphere   )  ! Attr
        call AddOutVar(IOVars,"doActiveShell",Model%doActiveShell     )  ! Attr
        call AddOutVar(IOVars,"nSpcIn",Model%nSpcMHD)  ! Attr
        call AddOutVar(IOVars,"nSpcRAIJU",Model%nSpc)  ! Attr
        call WriteVars(IOVars,.true.,Model%raijuH5)

        ! Root data
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",lat2D,uStr="radians")
        call AddOutVar(IOVars,"Y",lon2D,uStr="radians")
        call WriteVars(IOVars,.true.,Model%raijuH5)
        ! Grid data
        call AddOutVar(IOVars,"th_surf"  ,Grid%thRp  (is:ie+1),uStr="radians",dStr="(corner) Thetas mapped to 1 Rp")
        call AddOutVar(IOVars,"thcc_surf",Grid%thcRp (is:ie  ),uStr="radians",dStr="(cell) Thetas mapped to 1 Rp"  )
        call AddOutVar(IOVars,"Bmag"     ,Grid%Bmag  (is:ie,js:je),uStr="nT")
        call AddOutVar(IOVars,"BrCC"     ,Grid%BrCC  (is:ie,js:je),uStr="nT")
        call AddOutVar(IOVars,"areaCC"   ,Grid%areaCC(is:ie,js:je),uStr="Rp^2")
        call AddOutVar(IOVars,"BrFace"   ,Grid%BrFace(is:ie+1,js:je+1,:),uStr="nT")
        call AddOutVar(IOVars,"alamc"    ,Grid%alamc,uStr="eV * (Rp/nT)^(2/3)")
        call AddOutVar(IOVars,"areaFace",Grid%areaFace(is:ie+1,js:je+1,:),uStr="Rp^2")
        call AddOutVar(IOVars,"lenFace" ,Grid%lenFace (is:ie+1,js:je+1,:),uStr="Rp")
        call WriteVars(IOVars,.true.,Model%raijuH5,"Grid")

        call writeShellGrid(Grid%shGrid, Model%raijuH5, "Grid/ShellGrid")

        ! Output detailed lambda grid info
        do i=1,Grid%nSpc
            call ClearIO(IOVars)
            ! Attrs
            call AddOutVar(IOVars,"flav"     ,spc(i)%flav    )
            call AddOutVar(IOVars,"N"        ,spc(i)%N       )
            call AddOutVar(IOVars,"fudge"    ,spc(i)%fudge   )
            call AddOutVar(IOVars,"numNuc_p" ,spc(i)%numNuc_p)
            call AddOutVar(IOVars,"numNuc_n" ,spc(i)%numNuc_n)
            call AddOutVar(IOVars,"spcType"  ,spc(i)%spcType )
            call AddOutVar(IOVars,"q"        ,spc(i)%q       )
            call AddOutVar(IOVars,"amu"      ,spc(i)%amu     )
            call AddOutVar(IOVars,"kStart"   ,spc(i)%kStart-1)  ! Change to assume zero-based
            call AddOutVar(IOVars,"kEnd"     ,spc(i)%kEnd  -1)  ! Change to assume zero-based
            ! Datasets
            call AddOutVar(IOVars,"alami",spc(i)%alami,uStr="eV * (Rx/nT)^(2/3)")
            write(gStr,'(I0)') spc(i)%flav  ! Idk if this is the easiest way to format ints as strings
            call WriteVars(IOVars,.true.,Model%raijuH5,"Species",gStr)
        enddo

        ! Output planet info
        call writePlanetParams(Model%planet, .true., Model%raijuH5)
        end associate

    end subroutine raijuInitIO


    subroutine WriteRaiju(Model, Grid, State, gStr, doGhostsO)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: gStr
        logical, optional, intent(in) :: doGhostsO

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        integer :: i,j,k,s, nClkSteps
        integer :: is, ie, js, je, ks, ke
        integer, dimension(4) :: outBnds2D
        logical :: doGhosts
        real(rp) :: axT, axP
        real(rp) :: bVol_cc
        real(rp), dimension(:,:), allocatable :: outTmp2D
        logical , dimension(:,:), allocatable :: tmpGood2D
        real(rp), dimension(:,:,:), allocatable :: outTmp3D
        real(rp), dimension(:,:,:,:), allocatable :: outTmp4D
        logical , dimension(:,:,:), allocatable :: tmpGood3D
        real(rp), dimension(:,:,:), allocatable :: outPrecipN, outPrecipE, outPrecipAvgE, outCCHeatFlux

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

        ! First, make sure root variables are there
        if (doRoot) then
            call raijuInitIO(Model, Grid, doGhosts)
            doRoot = .false.
        endif
        !Reset IO chain
        call ClearIO(IOVars)

        ! Add attributes
        call AddOutVar(IOVars,"time",State%t)
        call AddOutVar(IOVars,"MJD",State%mjd)

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
        call AddOutVar(IOVars,"bminX",State%Bmin(is:ie+1,js:je+1,XDIR),uStr="nT")
        call AddOutVar(IOVars,"bminY",State%Bmin(is:ie+1,js:je+1,YDIR),uStr="nT")
        call AddOutVar(IOVars,"bminZ",State%Bmin(is:ie+1,js:je+1,ZDIR),uStr="nT")

        call AddOutVar(IOVars,"eta"    ,State%eta    (is:ie,js:je, :),uStr="#/cm^3 * Rx/T") !! TODO: Maybe swap with just intensity instead
        call AddOutVar(IOVars,"eta_avg",State%eta_avg(is:ie,js:je, :),uStr="#/cm^3 * Rx/T") !! TODO: Maybe swap with just intensity instead

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
                        State%eta (i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd))
                enddo
            enddo
        enddo
        call AddOutVar(IOVars,"intensity",outTmp3D(is:ie,js:je, :),uStr="1/(s*sr*keV*cm^2)")
        deallocate(outTmp3D)

        ! Coupling things
        call AddOutVar(IOVars,"dtCpl"  ,State%dt)  ! Attr
        call AddOutVar(IOVars,"xmin"   ,State%xyzMin (is:ie+1,js:je+1,XDIR)  ,uStr="Rx",dStr="(corners) X location of bmin surface")
        call AddOutVar(IOVars,"ymin"   ,State%xyzMin (is:ie+1,js:je+1,YDIR)  ,uStr="Rx",dStr="(corners) Y location of bmin surface")
        call AddOutVar(IOVars,"zmin"   ,State%xyzMin (is:ie+1,js:je+1,ZDIR)  ,uStr="Rx",dStr="(corners) Z location of bmin surface")
        call AddOutVar(IOVars,"topo"   ,State%topo   (is:ie+1,js:je+1)*1.0_rp,uStr="0=Open, 1=Closed",dStr="(corners) Magnetic topology")
        call AddOutVar(IOVars,"colatc" ,State%thcon  (is:ie+1,js:je+1)       ,uStr="radians",dStr="(corners) Congugate latitude")
        call AddOutVar(IOVars,"lonc"   ,State%phcon  (is:ie+1,js:je+1)       ,uStr="radians",dStr="(corners) Congugate longitude")
        call AddOutVar(IOVars,"active" ,State%active (is:ie,js:je)*1.0_rp    ,uStr="-1=Inactive, 0=Buffer, 1=Active")
        call AddOutVar(IOVars,"OCBDist",State%OCBDist(is:ie,js:je)*1.0_rp    ,dStr="Cell distance from an open closed boundary")
        call AddOutVar(IOVars,"bVol"   ,State%bvol   (is:ie+1,js:je+1)       ,uStr="Rx/nT",dStr="(corners) Flux Tube Volume")
        call AddOutVar(IOVars,"bVol_cc",State%bvol_cc(is:ie,js:je)           ,uStr="Rx/nT",dStr="(corners) Flux Tube Volume")
        call AddOutVar(IOVars,"vaFrac" ,State%vaFrac (is:ie+1,js:je+1)       ,uStr="fraction",dStr="Fraction of Alfven speed over magnetofast + ExB speed")
        call AddOutVar(IOVars,"Pavg_in",State%Pavg   (is:ie,js:je, :)        ,uStr="nPa" ,dStr="Pressures from imagtubes")
        call AddOutVar(IOVars,"Davg_in",State%Davg   (is:ie,js:je, :)        ,uStr="#/cc",dStr="Densities from imagtubes")
        call AddOutVar(IOVars,"Pstd_in",State%Pstd   (is:ie,js:je, :)        ,uStr="normalized" ,dStr="Std. dev. of species pressure from imagtubes")
        call AddOutVar(IOVars,"Dstd_in",State%Dstd   (is:ie,js:je, :)        ,uStr="normalized" ,dStr="Std. dev. of species density from imagtubes")
        call AddOutVar(IOVars,"tiote"  ,State%tiote  (is:ie,js:je)           ,uStr="normalized" ,dStr="Ratio of ion temperature over electron temperature")
        call AddOutVar(IOVars,"domWeights",State%domWeights(is:ie,js:je)     ,                   dStr="Weights used for smoothing moments for questionable cells along boundary")
        call AddOutSGV(IOVars,"Tbounce",State%Tb, outBndsO=outBnds2D, uStr="[s]", dStr="Bounce timescale along field line (Alfven crossing time)", doWriteMaskO=.false.)

        if (Model%doOwnCorot) then
            ! If we are addng our own corotation potential, then its not included in espot
            call AddOutVar(IOVars,"pot_iono",State%espot(is:ie+1,js:je+1)                                  ,uStr="kV",dStr="(corners) Ionospheric potential")
        else
            ! Otherwise, its in espot, so substract it out
            call AddOutVar(IOVars,"pot_iono",State%espot(is:ie+1,js:je+1)-State%pot_corot(is:ie+1,js:je+1) ,uStr="kV",dStr="(corners) Ionospheric potential")
        endif
        call AddOutVar(IOVars,"pot_corot",State%pot_corot(is:ie+1,js:je+1),uStr="kV",dStr="(corners) Corotation potential")

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
        call AddOutSGV(IOVars, "Pressure", State%Press, outBndsO=outBnds2D, uStr="nPa" , dStr="Pressure from RAIJU flavors", doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Density" , State%Den  , outBndsO=outBnds2D, uStr="#/cc", dStr="Density from RAIJU flavors" , doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Pressure_avg", State%Press_avg, outBndsO=outBnds2D, uStr="nPa" , dStr="Step-averaged Pressure from RAIJU flavors", doWriteMaskO=.false.)
        call AddOutSGV(IOVars, "Density_avg" , State%Den_avg  , outBndsO=outBnds2D, uStr="#/cc", dStr="Step-averaged Density from RAIJU flavors" , doWriteMaskO=.false.)
        

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
            ! Calculate accumulated precipitation for each species
            allocate(outPrecipN   (is:ie,js:je,Grid%nSpc))
            allocate(outPrecipE   (is:ie,js:je,Grid%nSpc))
            allocate(outPrecipAvgE(is:ie,js:je,Grid%nSpc))
            allocate(outCCHeatFlux(is:ie,js:je,Grid%nSpc))
            outPrecipN    = 0.0_rp
            outPrecipE    = 0.0_rp
            outPrecipAvgE = 0.0_rp
            outCCHeatFlux = 0.0_rp
            do s=1,Grid%nSpc
                ks = Grid%spc(s)%kStart
                ke = Grid%spc(s)%kEnd
                do k=ks,ke !1,Grid%Nk
                    outPrecipN(:,:,s) = outPrecipN(:,:,s) + State%precipNFlux(k)%data(is:ie,js:je)
                    outPrecipE(:,:,s) = outPrecipE(:,:,s) + State%precipEFlux(k)%data(is:ie,js:je)
                enddo
                outCCHeatFlux(:,:,s) = sum(State%CCHeatFlux (is:ie,js:je,kS:kE), dim=3)

                where (outPrecipN(:,:,s) > TINY)
                    outPrecipAvgE(:,:,s) = outPrecipE(:,:,s)/outPrecipN(:,:,s) * erg2kev  ! Avg E [keV]
                elsewhere
                    outPrecipAvgE(:,:,s) = 0.0
                end where
            enddo
            call AddOutVar(IOVars, "precipNFlux", outPrecipN   , uStr="#/cm^2/s")
            call AddOutVar(IOVars, "precipEFlux", outPrecipE   , uStr="erg/cm^2/s")
            call AddOutVar(IOVars, "precipAvgE" , outPrecipAvgE, uStr="keV")
            call AddOutVar(IOVars, "CCHeatFlux" , outCCHeatFlux, uStr="eV/cm^2/s")
            deallocate(outPrecipN)
            deallocate(outPrecipE)
            deallocate(outPrecipAvgE)
            deallocate(outCCHeatFlux)
            ! (Ni, Nj, Nk)
            call AddOutVar(IOVars, "lossRate", State%lossRates(is:ie,js:je,:), uStr="1/s")
            call AddOutVar(IOVars, "lossRatePrecip", State%lossRatesPrecip(is:ie,js:je,:), uStr="1/s")

            if (Model%doOutput_3DLoss) then
                call AddOutVar(IOVars, "dEta_dt" , State%dEta_dt(is:ie,js:je,:), uStr="eta_units/s")
                call AddOutVar(IOVars, "CCHeatFlux_Nk" , State%CCHeatFlux (is:ie,js:je,:), uStr="eV/cm^2/s")
                call AddOutSGV(IOVars, "precipNFlux_Nk", State%precipNFlux, outBndsO=outBnds2D, uStr="#/cm^2/s" , dStr="precipNFlux from RAIJU flavors", doWriteMaskO=.false.)
                call AddOutSGV(IOVars, "precipEFlux_Nk", State%precipEFlux, outBndsO=outBnds2D, uStr="erg/cm^2/s" , dStr="precipEFlux from RAIJU flavors", doWriteMaskO=.false.)
            endif
        endif
        
        if (Model%doOutput_potGrads) then
            ! (Ni, Nj, 2)
            call AddOutVar(IOVars, "gradPotE"       , State%gradPotE    (is:ie+1,js:je+1,:), uStr="V/m")
            call AddOutVar(IOVars, "gradPotCorot"   , State%gradPotCorot(is:ie+1,js:je+1,:), uStr="V/m")
            call AddOutVar(IOVars, "gradVM"         , State%gradVM      (is:ie+1,js:je+1,:), uStr="V/m/lambda")
            call AddOutVar(IOVars, "gradPotE_cc"    , State%gradPotE_cc    (is:ie,js:je,:) , uStr="V/m")
            call AddOutVar(IOVars, "gradPotCorot_cc", State%gradPotCorot_cc(is:ie,js:je,:) , uStr="V/m")
            call AddOutVar(IOVars, "gradVM_cc"      , State%gradVM_cc      (is:ie,js:je,:) , uStr="V/m/lambda")

            ! Calc pEffective based on current state
            ! Make full ghost size since that's what the subroutine expects
            allocate(outTmp3D(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1,Grid%Nk))
            call calcEffectivePotential(Model, Grid, State, outTmp3D)
            call AddOutVar(IOVars, "pEffective", outTmp3D(is:ie+1,js:je+1,:)*1e-3, uStr="kV")
            deallocate(outTmp3D)
        endif

        if (Model%doOutput_debug) then
            call AddOutVar(IOVars, "eta_half"     , State%eta_half     (is:ie  ,js:je  ,:)  , uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "iVel"         , State%iVel         (is:ie+1,js:je+1,:,:), uStr="m/s")
            call AddOutVar(IOVars, "iVelL"        , State%iVelL        (is:ie+1,js:je+1,:,:), uStr="m/s")
            call AddOutVar(IOVars, "iVelR"        , State%iVelR        (is:ie+1,js:je+1,:,:), uStr="m/s")
            call AddOutVar(IOVars, "cVel"         , State%cVel         (is:ie  ,js:je  ,:,:), uStr="m/s")
            call AddOutVar(IOVars, "etaFaceReconL", State%etaFaceReconL(is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFaceReconR", State%etaFaceReconR(is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFacePDML"  , State%etaFacePDML  (is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFacePDMR"  , State%etaFacePDMR  (is:ie+1,js:je+1,:,:), uStr="#/cm^3 * Rx/T")
            call AddOutVar(IOVars, "etaFlux"      , State%etaFlux      (is:ie+1,js:je+1,:,:), uStr="(#/cm^3 * Rx/T) * m^2 / s")
            
            ! Unsmoothed gradVMcc
            allocate(tmpGood2D(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1   ))
            allocate(outTmp3D (Grid%shGrid%isg:Grid%shGrid%ieg  ,Grid%shGrid%jsg:Grid%shGrid%jeg  ,2))
            allocate(outTmp4D (Grid%shGrid%isg:Grid%shGrid%ieg  ,Grid%shGrid%jsg:Grid%shGrid%jeg,Grid%Nk ,2))
            tmpGood2D = .false.
            where (State%topo .eq. RAIJUCLOSED)
                tmpGood2D = .true.                
            end where
            call calcGradVM_cc(Model%planet%rp_m, Model%planet%ri_m, Model%planet%magMoment, &
                                Grid, tmpGood2D, State%bvol, outTmp3D, doSmoothO=.false., doLimO=.false.)
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(k)
            do k=1,Grid%Nk
                call calcVelocity_cc(Model, Grid, State, k, outTmp4D(:,:,k,:), outTmp3D)
            enddo
            call AddOutVar(IOVars, "gradVM_cc_nosm", outTmp3D(is:ie,js:je,:)  , uStr="V/m/lambda")
            call AddOutVar(IOVars, "cVel_nosm"     , outTmp4D(is:ie,js:je,:,:), uStr="m/s")
            call calcGradVM_cc(Model%planet%rp_m, Model%planet%ri_m, Model%planet%magMoment, &
                                Grid, tmpGood2D, State%bvol, outTmp3D, doSmoothO=.true., doLimO=.false.)
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(k)
            do k=1,Grid%Nk
                call calcVelocity_cc(Model, Grid, State, k, outTmp4D(:,:,k,:), outTmp3D)
            enddo
            call AddOutVar(IOVars, "gradVM_cc_sm_nolim", outTmp3D(is:ie,js:je,:)  , uStr="V/m/lambda")
            call AddOutVar(IOVars, "cVel_sm_nolim"     , outTmp4D(is:ie,js:je,:,:), uStr="m/s")
            deallocate(tmpGood2D)
            deallocate(outTmp3D)
            deallocate(outTmp4D)

            ! Mapping diagnostic
            allocate(outTmp2D(Grid%shGrid%isg:Grid%shGrid%ieg  ,Grid%shGrid%jsg:Grid%shGrid%jeg))
            call calcMapJacNorm(Grid, State%xyzMin, outTmp2D)
            call AddOutVar(IOVars, "mapJacNorm", outTmp2D(is:ie,js:je), dStr="L_(2,1) norm of lat/lon => xyzMin Jacobian")
            deallocate(outTmp2D)
        endif

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


    subroutine WriteRaijuRes(Model, Grid, State, ResF)
        !! Writes RAIJU restart info to provided path ResF
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: ResF

        ! If a restart already exists, get rid of old one
        call CheckAndKill(ResF)

        ! Save ShellGrid to root of file
        call writeShellGrid(Grid%shGrid, ResF)
        ! And species info
        call writeSpeciesInfo(Model, Grid, ResF)
        ! All necessary State info
        call WriteRaijuResState(Model, Grid, State, ResF)

    end subroutine WriteRaijuRes


    subroutine WriteRaijuResState(Model, Grid, State, ResF)
        !! Writes RAIJU State restart info to provided path ResF
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: ResF

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        integer :: i,j,k,s
        integer :: is, ie, js, je, ks, ke
        real(rp), dimension(:,:), allocatable :: outActiveShell
        real(rp), dimension(:,:,:), allocatable :: tmpOut3D

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
        ! ColdStarter state
        if (State%coldStarter%doneFirstCS) then
            call AddOutVar(IOVars,"cs_doneFirstCS", 1)
        else
            call AddOutVar(IOVars,"cs_doneFirstCS", 0)
        endif
        call AddOutVar(IOVars, "cs_lastEval", State%coldStarter%lastEval)
        call AddOutVar(IOVars, "cs_lastTarget", State%coldStarter%lastTarget)

        ! Add State variables

        ! Coupling things
        call AddOutVar(IOVars,"dtCpl"  ,State%dt)  ! Attr
        call AddOutVar(IOVars,"xmin"   ,State%xyzMin (:,:,XDIR)  ,uStr="Rx",dStr="(corners) X location of bmin surface")
        call AddOutVar(IOVars,"ymin"   ,State%xyzMin (:,:,YDIR)  ,uStr="Rx",dStr="(corners) Y location of bmin surface")
        call AddOutVar(IOVars,"zmin"   ,State%xyzMin (:,:,ZDIR)  ,uStr="Rx",dStr="(corners) Z location of bmin surface")
        call AddOutVar(IOVars,"bminX"  ,State%Bmin   (:,:,XDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"bminY"  ,State%Bmin   (:,:,YDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"bminZ"  ,State%Bmin   (:,:,ZDIR)  ,uStr="nT")
        call AddOutVar(IOVars,"topo"   ,State%topo   (:,:)*1.0_rp,uStr="0=Open, 1=Closed",dStr="(corners) Magnetic topology")
        call AddOutVar(IOVars,"colatc" ,State%thcon  (:,:)       ,uStr="radians",dStr="(corners) Congugate latitude")
        call AddOutVar(IOVars,"lonc"   ,State%phcon  (:,:)       ,uStr="radians",dStr="(corners) Congugate longitude")
        call AddOutVar(IOVars,"espot"  ,State%espot  (:,:)       ,uStr="kV",dStr="(corners) Electrostatic potential")
        call AddOutVar(IOVars,"bVol"   ,State%bvol   (:,:)       ,uStr="Rx/nT",dStr="(corners) Flux Tube Volume")
        call AddOutVar(IOVars,"bVolcc" ,State%bvol_cc(:,:)       ,uStr="Rx/nT",dStr="(centers) Flux Tube Volume")
        call AddOutVar(IOVars,"vaFrac" ,State%vaFrac (:,:)       ,uStr="fraction",dStr="Fraction of Alfven speed over magnetofast + ExB speed")

        ! If only 1 element in 3rd position, AddOutVar will write as 2D array
        ! We know it should be 3D, so force it
        call AddOutVar(IOVars,"Pavg_in",State%Pavg   (:,:,:)     ,uStr="nPa" ,dStr="Pressures from imagtubes", doSqzO=.false.)
        call AddOutVar(IOVars,"Davg_in",State%Davg   (:,:,:)     ,uStr="#/cc",dStr="Densities from imagtubes", doSqzO=.false.)

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
        call AddOutSGV(IOVars, "Pressure", State%Press, doWriteMaskO=.true., uStr="nPa")
        call AddOutSGV(IOVars, "Density" , State%Den  , doWriteMaskO=.true., uStr="#/cc")
        call AddOutSGV(IOVars, "Pressure_avg", State%Press_avg, doWriteMaskO=.true., uStr="nPa")
        call AddOutSGV(IOVars, "Density_avg" , State%Den_avg  , doWriteMaskO=.true., uStr="#/cc")
        ! Precip
        call AddOutSGV(IOVars, "precipNFlux_Nk", State%precipNFlux, doWriteMaskO=.true., uStr="#/cm^2/s")
        call AddOutSGV(IOVars, "precipEFlux_Nk", State%precipEFlux, doWriteMaskO=.true., uStr="erg/cm^2/s")
        call AddOutVar(IOVars,"precipLossRates_Nk", State%lossRates(:,:,:), uStr="1/s")
        ! (Probably not needed but we will save anyways)
        call AddOutVar(IOVars, "gradPotE"    , State%gradPotE     (:,:,:), uStr="V/m")
        call AddOutVar(IOVars, "gradPotCorot", State%gradPotCorot (:,:,:), uStr="V/m")
        call AddOutVar(IOVars, "gradVM"      , State%gradVM       (:,:,:), uStr="V/m/lambda")
        ! More solver stuff
        call AddOutVar(IOVars, "dtk"   , State%dtk(:), uStr="s")
        call AddOutVar(IOVars, "nStepk", State%nStepk(:)*1.0_rp, uStr="#", dStr="Number of steps each channel has taken")
        call AddOutVar(IOVars, "iVel"   , State%iVel  (:,:,:,:)     , uStr="m/s")
        call AddOutVar(IOVars, "cVel"   , State%cVel  (:,:,:,:), uStr="m/s")

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
        real(rp), dimension(:,:,:), allocatable :: tmpReal3D

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
        call AddInVar(IOVars,"cs_doneFirstCS")
        call AddInVar(IOVars,"cs_lastEval")
        call AddInVar(IOVars,"cs_lastTarget")

        call AddInVar(IOVars,"dtCpl"  )
        call AddInVar(IOVars,"xmin"   )
        call AddInVar(IOVars,"ymin"   )
        call AddInVar(IOVars,"zmin"   )
        call AddInVar(IOVars,"bminX"  )
        call AddInVar(IOVars,"bminY"  )
        call AddInVar(IOVars,"bminZ"  )
        call AddInVar(IOVars,"topo"   )
        call AddInVar(IOVars,"colatc" )
        call AddInVar(IOVars,"lonc"   )
        call AddInVar(IOVars,"espot"  )
        call AddInVar(IOVars,"bVol"   )
        call AddInVar(IOVars,"bVolcc" )
        call AddInVar(IOVars,"vaFrac" )
        call AddInVar(IOVars,"Pavg_in")
        call AddInVar(IOVars,"Davg_in")

        call AddInVar(IOVars,"eta"        )
        call AddInVar(IOVars,"eta_last"   )
        call AddInVar(IOVars,"eta_avg"    )
        call AddInVar(IOVars,"active"     )
        call AddInVar(IOVars,"active_last")
        call AddInVar(IOVars,"OCBDist"    )
        call AddInVar(IOVars,"activeShells")

        call AddInVar(IOVars,"precipLossRates_Nk")

        call AddInVar(IOVars, "gradPotE"    )
        call AddInVar(IOVars, "gradPotCorot")
        call AddInVar(IOVars, "gradVM"      )

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

        ! Coldstarter
        State%coldStarter%lastEval   = GetIOReal(IOVars, "cs_lastEval")
        State%coldStarter%lastTarget = GetIOReal(IOVars, "cs_lastTarget")
        
        ! Handle boolean attributes
        tmpInt = GetIOInt(IOVars, "isFirstCpl")
        State%isFirstCpl = tmpInt .eq. 1
        tmpInt = GetIOInt(IOVars, "cs_doneFirstCS")
        State%coldStarter%doneFirstCS = tmpInt .eq. 1

        call IOArray2DFill(IOVars, "xmin" , State%xyzMin(:,:,XDIR))
        call IOArray2DFill(IOVars, "ymin" , State%xyzMin(:,:,YDIR))
        call IOArray2DFill(IOVars, "zmin" , State%xyzMin(:,:,ZDIR))
        call IOArray2DFill(IOVars, "bminX", State%Bmin  (:,:,XDIR))
        call IOArray2DFill(IOVars, "bminY", State%Bmin  (:,:,YDIR))
        call IOArray2DFill(IOVars, "bminZ", State%Bmin  (:,:,ZDIR))
        call IOArray2DFill(IOVars, "colatc", State%thcon  (:,:))
        call IOArray2DFill(IOVars, "lonc"  , State%phcon  (:,:))
        call IOArray2DFill(IOVars, "espot" , State%espot  (:,:))
        call IOArray2DFill(IOVars, "colatc", State%thcon  (:,:))
        call IOArray2DFill(IOVars, "lonc"  , State%phcon  (:,:))
        call IOArray2DFill(IOVars, "espot" , State%espot  (:,:))
        call IOArray2DFill(IOVars, "bVol"  , State%bvol   (:,:))
        call IOArray2DFill(IOVars, "bVolcc", State%bvol_cc(:,:))
        call IOArray2DFill(IOVars, "vaFrac", State%vaFrac (:,:))
        call IOArray3DFill(IOVars, "Pavg_in", State%Pavg(:,:,:))
        call IOArray3DFill(IOVars, "Davg_in", State%Davg(:,:,:))
        

        call IOArray3DFill(IOVars, "eta"     , State%eta     (:,:,:))
        call IOArray3DFill(IOVars, "eta_last", State%eta_last(:,:,:))
        call IOArray3DFill(IOVars, "eta_avg" , State%eta_avg (:,:,:))

        call IOArray3DFill(IOVars, "precipLossRates_Nk", State%lossRates(:,:,:))

        call IOArray3DFill(IOVars, "gradPotE"    , State%gradPotE    (:,:,:))
        call IOArray3DFill(IOVars, "gradPotCorot", State%gradPotCorot(:,:,:))
        call IOArray3DFill(IOVars, "gradVM"      , State%gradVM      (:,:,:))
        
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
        allocate(tmpReal1D(Grid%Nk))
        call IOArray1DFill(IOVars, "nStepk", tmpReal1D)
        State%nStepk = INT(tmpReal1D)

        ! ShellGridVars
        call ReadInSGV(State%Press, inH5, "Pressure", "State", doIOpO=.false.)
        call ReadInSGV(State%Den  , inH5, "Density" , "State", doIOpO=.false.)
        call ReadInSGV(State%Press_avg, inH5, "Pressure_avg", "State", doIOpO=.false.)
        call ReadInSGV(State%Den_avg  , inH5, "Density_avg" , "State", doIOpO=.false.)
        call ReadInSGV(State%precipNFlux, inH5, "precipNFlux_Nk", "State", doIOpO=.false.)
        call ReadInSGV(State%precipEFlux, inH5, "precipEFlux_Nk", "State", doIOpO=.false.)

        end associate

    end subroutine ReadRaijuResState

end module raijuIO
