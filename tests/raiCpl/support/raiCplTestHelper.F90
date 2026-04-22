module raiCplTestHelper

    use kdefs, only : rp, PI, strLen, TINY, NDIM, XDIR, YDIR, ZDIR, rad2deg

    use shellGrid
    use xml_input
    use planethelper, only : DipoleShift
    use volttypes, only : raijuCoupler_T, imagOptions_T, IMAG_nMHD_R_NEAREARTH, voltApp_T
    use voltCplTypes, only : TubeShell_T, MAXTUBEFLUIDS, TUBE_CLOSED, TUBE_OPEN
    use raijutypes, only : raijuApp_T
    use raijustarter, only : raijuInit
    use raijuCplHelper, only : raijuCpl_init, setRaijuILowLat_MHD
    use tubehelper, only : init_TubeShell

    implicit none

    integer, parameter :: raiCplTestNG = 4
        !! Assume 4 radial ghost cells in the LFM grid files used by raiCpl tests
    integer, parameter :: raiCplTestNGb = raiCplTestNG + 1
        !! Need 5 corner radii to bound 4 ghost cells
    integer, parameter :: nRaiCplTestProfiles = 6
    character(len=9), parameter :: raiCplTestProfileNames(nRaiCplTestProfiles) = [character(len=9) :: &
        "lfmD_1.25", "lfmD_2p0", "lfmO_1p25", "lfmO_2p0", "lfmQ_1p25", "lfmQ_2p0" ]

contains

    subroutine genGlobalVoltParentSG(parentSG, NtO, NpO, NgPhiO, radO, nameO)
        !! Build a Voltron-style global parent shell grid:
        !! full theta span, no theta ghosts, phi ghosts only.
        type(ShellGrid_T), intent(out) :: parentSG
        integer, intent(in), optional :: NtO, NpO, NgPhiO
        real(rp), intent(in), optional :: radO
        character(len=*), intent(in), optional :: nameO

        integer :: Nt, Np, NgPhi, i
        integer, dimension(4) :: nGhosts
        real(rp) :: rad
        real(rp), allocatable :: theta(:), phi(:)
        character(len=strLen) :: gridName

        Nt = 135
        if (present(NtO)) Nt = NtO

        Np = 28
        if (present(NpO)) Np = NpO

        NgPhi = 4
        if (present(NgPhiO)) NgPhi = NgPhiO

        rad = 1.1_rp
        if (present(radO)) rad = radO

        gridName = "voltParentSG"
        if (present(nameO)) gridName = nameO

        allocate(theta(Nt+1), phi(Np+1))

        do i=1,Nt+1
            theta(i) = PI*(i-1)/Nt
        enddo

        do i=1,Np+1
            phi(i) = 2.0_rp*PI*(i-1)/Np
        enddo

        nGhosts = (/0, 0, NgPhi, NgPhi/)
        call GenShellGrid(parentSG, theta, phi, trim(gridName), nGhosts=nGhosts, radO=rad)
    end subroutine genGlobalVoltParentSG


    subroutine getRaiCplTestMhdProfile(profileName, mhdRAll, mhdRinGcc, mhdRNearEarth, mhdRin)
        !! Return one LFM-derived radial profile along the x-axis, laid out as:
        !! 5 ghost-bounding corner radii followed by IMAG_nMHD_R_NEAREARTH active corner radii.
        !! The active list starts at mhdRin, so the ghost/active interface corner is duplicated
        !! to keep the profile easy for humans to read.
        character(len=*), intent(in) :: profileName
        real(rp), dimension(raiCplTestNGb + IMAG_nMHD_R_NEAREARTH), intent(out) :: mhdRAll
        real(rp), dimension(raiCplTestNG), intent(out) :: mhdRinGcc
        real(rp), dimension(IMAG_nMHD_R_NEAREARTH), intent(out) :: mhdRNearEarth
        real(rp), intent(out) :: mhdRin
        integer :: i

        select case (trim(profileName))
        case ("lfmD_1.25")
            mhdRAll = [ &
                1.1500000099999999_rp, 1.1750000074999998_rp, 1.2000000050000000_rp, 1.2250000024999999_rp, 1.2500000000000004_rp, &
                1.2500000000000004_rp, 1.5255736251434022_rp, 1.8302668923949814_rp, 2.1657695719101710_rp, &
                2.5336715532140666_rp, 2.9356098099740900_rp, 3.3732553968872160_rp, 3.8483185151678770_rp, &
                4.3625437543068810_rp, 4.9177082923334140_rp ]
        case ("lfmD_2p0")
            mhdRAll = [ &
                1.1500000099999999_rp, 1.3625000075000000_rp, 1.5750000049999997_rp, 1.7875000024999996_rp, 2.0000000000000018_rp, &
                2.0000000000000018_rp, 2.4012156664761326_rp, 2.8160939941865150_rp, 3.2418104973610213_rp, &
                3.6754498544192820_rp, 4.1140536736815250_rp, 4.5546815639024745_rp, 4.9943837997034430_rp, &
                5.4302169696985740_rp, 5.8592369004197264_rp ]
        case ("lfmO_1p25")
            mhdRAll = [ &
                1.1500000099999999_rp, 1.1750000074999998_rp, 1.2000000050000000_rp, 1.2250000024999999_rp, 1.2500000000000004_rp, &
                1.2500000000000004_rp, 1.3163451108595015_rp, 1.3842989289111058_rp, 1.4540348666045437_rp, &
                1.5255736251434022_rp, 1.5989929496902840_rp, 1.6742154783069854_rp, 1.7512747320248350_rp, &
                1.8302668923949814_rp, 1.9112255175308180_rp ]
        case ("lfmO_2p0")
            mhdRAll = [ &
                1.6020335503721244_rp, 1.7015251627790930_rp, 1.8010167751860620_rp, 1.9005083875930306_rp, 2.0000000000000018_rp, &
                2.0000000000000018_rp, 2.0991431188506957_rp, 2.1989475811182880_rp, 2.2996286965727517_rp, &
                2.4012156664761326_rp, 2.5037519486306685_rp, 2.6071129165499474_rp, 2.7112370379813204_rp, &
                2.8160939941865150_rp, 2.9216406148843164_rp ]
        case ("lfmQ_1p25")
            mhdRAll = [ &
                1.1500000099999999_rp, 1.1750000074999998_rp, 1.2000000050000000_rp, 1.2250000024999999_rp, 1.2500000000000004_rp, &
                1.2500000000000004_rp, 1.3842989289111058_rp, 1.5255736251434022_rp, 1.6742154783069854_rp, &
                1.8302668923949814_rp, 1.9940958226591960_rp, 2.1657695719101710_rp, 2.3455788965803777_rp, &
                2.5336715532140666_rp, 2.7302708471347470_rp ]
        case ("lfmQ_2p0")
            mhdRAll = [ &
                1.1976734237411777_rp, 1.3982550678058832_rp, 1.5988367118705886_rp, 1.7994183559352940_rp, 2.0000000000000018_rp, &
                2.0000000000000018_rp, 2.1989475811182880_rp, 2.4012156664761326_rp, 2.6071129165499474_rp, &
                2.8160939941865150_rp, 3.0277973211550386_rp, 3.2418104973610213_rp, 3.4578309010436160_rp, &
                3.6754498544192820_rp, 3.8943167941395957_rp ]
        case default
            write(*,*) "ERROR in getRaiCplTestMhdProfile: unsupported profileName=", trim(profileName)
            stop
        end select

        do i=1,raiCplTestNG
            mhdRinGcc(i) = 0.5_rp*(mhdRAll(i) + mhdRAll(i+1))
        enddo
        mhdRNearEarth = mhdRAll(raiCplTestNGb:raiCplTestNGb+IMAG_nMHD_R_NEAREARTH-1)
        mhdRin = 0.5_rp*(mhdRAll(raiCplTestNGb) + mhdRAll(raiCplTestNGb+1))
    end subroutine getRaiCplTestMhdProfile


    subroutine initRaiCplTestOpt(raiCpl, parentSG, mhdRin, mhdRNearEarthO, mhdRinGccO, mhdDrO, mjd0O, swFO, doColdStartO)
        !! Populate the same imag option fields that InitInnerMag sets at runtime.
        !! The early raiCpl routine tests mainly need voltGrid, mhdRin, mhdR_nearEarth_arr, and mjd0.
        class(raijuCoupler_T), intent(inout) :: raiCpl
        type(ShellGrid_T), intent(in) :: parentSG
        real(rp), intent(in) :: mhdRin
        real(rp), dimension(IMAG_nMHD_R_NEAREARTH), intent(in), optional :: mhdRNearEarthO
        real(rp), dimension(:), intent(in), optional :: mhdRinGccO
        real(rp), intent(in), optional :: mhdDrO, mjd0O
        character(len=*), intent(in), optional :: swFO
        logical, intent(in), optional :: doColdStartO

        integer :: i
        real(rp) :: mhdDr, mjd0
        logical :: doColdStart
        character(len=strLen) :: swF

        if (allocated(raiCpl%opt)) deallocate(raiCpl%opt)
        allocate(imagOptions_T :: raiCpl%opt)

        mhdDr = 0.1_rp
        if (present(mhdDrO)) mhdDr = mhdDrO

        mjd0 = 60000.0_rp
        if (present(mjd0O)) mjd0 = mjd0O

        swF = "NONE"
        if (present(swFO)) swF = swFO

        doColdStart = .false.
        if (present(doColdStartO)) doColdStart = doColdStartO

        raiCpl%opt%swF = swF
        raiCpl%opt%mhdRin = mhdRin
        raiCpl%opt%mjd0 = mjd0
        raiCpl%opt%doColdStart = doColdStart
        raiCpl%opt%voltGrid = parentSG

        if (present(mhdRNearEarthO)) then
            raiCpl%opt%mhdR_nearEarth_arr = mhdRNearEarthO
        else
            do i=1,IMAG_nMHD_R_NEAREARTH
                raiCpl%opt%mhdR_nearEarth_arr(i) = mhdRin + (i-1)*mhdDr
            enddo
        endif

        if (allocated(raiCpl%opt%mhdRinGcc_arr)) deallocate(raiCpl%opt%mhdRinGcc_arr)
        if (present(mhdRinGccO)) then
            allocate(raiCpl%opt%mhdRinGcc_arr(size(mhdRinGccO)))
            raiCpl%opt%mhdRinGcc_arr = mhdRinGccO
        else
            allocate(raiCpl%opt%mhdRinGcc_arr(raiCplTestNG))
            do i=1,raiCplTestNG
                raiCpl%opt%mhdRinGcc_arr(i) = mhdRin - (raiCplTestNG-i+1)*mhdDr
            enddo
        endif
    end subroutine initRaiCplTestOpt


    subroutine initReadyRaiCpl(raiCpl, xmlInp, parentSG, profileName)
        !! Build a fully initialized raiCpl fixture for routine-level tests.
        !! This follows the same high-level startup path as production, but uses
        !! a deterministic shell grid and one of the canned near-Earth MHD profiles.
        class(raijuCoupler_T), intent(inout) :: raiCpl
        type(XML_Input_T), intent(inout) :: xmlInp
        type(ShellGrid_T), intent(in) :: parentSG
        character(len=*), intent(in) :: profileName

        real(rp), dimension(raiCplTestNGb + IMAG_nMHD_R_NEAREARTH) :: mhdRAll
        real(rp), dimension(raiCplTestNG) :: mhdRinGcc
        real(rp), dimension(IMAG_nMHD_R_NEAREARTH) :: mhdRNearEarth
        real(rp) :: mhdRin

        if (allocated(raiCpl%raiApp)) deallocate(raiCpl%raiApp)
        allocate(raijuApp_T :: raiCpl%raiApp)

        call getRaiCplTestMhdProfile(profileName, mhdRAll, mhdRinGcc, mhdRNearEarth, mhdRin)
        call initRaiCplTestOpt(raiCpl, parentSG, mhdRin, mhdRNearEarthO=mhdRNearEarth, mhdRinGccO=mhdRinGcc, &
                               mjd0O=60321.25_rp, swFO="bcwind_2step.h5")
        call applyRuntimeThetaUOverride(raiCpl)
        call raijuInit(raiCpl%raiApp, xmlInp, parentSG, doVerbXMLO=.false.)
        call raijuCpl_init(raiCpl, xmlInp)
        call setRaijuILowLat_MHD(raiCpl, doVerboseO=.false.)
    end subroutine initReadyRaiCpl


    subroutine initVoltAppTestState(vApp, parentSG)
        !! Prepare just enough Voltron-side state for pack/coupling tests:
        !! a parent shell grid, an initialized TubeShell, and corner potentials.
        type(voltApp_T), intent(inout) :: vApp
        type(ShellGrid_T), intent(in) :: parentSG

        vApp%shGrid = parentSG
        call init_TubeShell(vApp%shGrid, vApp%State%tubeShell)
        call initShellVar(vApp%shGrid, SHGR_CORNER, vApp%State%potential_total)
        call initShellVar(vApp%shGrid, SHGR_CORNER, vApp%State%potential_corot)
        vApp%time = 0.0_rp
    end subroutine initVoltAppTestState


    subroutine fillSyntheticCornerVar(var, offset)
        !! Fill a shell variable with a simple affine pattern in i/j so tests can
        !! recognize interpolation/copy behavior without ambiguous repeated values.
        type(ShellGridVar_T), intent(inout) :: var
        real(rp), intent(in) :: offset

        integer :: i, j

        do j=var%jsv,var%jev
            do i=var%isv,var%iev
                var%data(i,j) = offset + 10.0_rp*real(i, rp) + real(j, rp)
            enddo
        enddo
    end subroutine fillSyntheticCornerVar


    subroutine fillSyntheticTubeShell(parentSG, tubeShell, firstClosedIO)
        !! Fill a TubeShell with deterministic synthetic data.
        !! Most fields use simple smooth patterns, while the optional topology
        !! front creates a phi-dependent poleward-open / equatorward-closed mask.
        type(ShellGrid_T), intent(in) :: parentSG
        type(TubeShell_T), intent(inout) :: tubeShell
        integer, intent(in), optional :: firstClosedIO

        integer :: i, j, n, s
        integer :: firstClosedI
        real(rp) :: phase, frontShift

        call init_TubeShell(parentSG, tubeShell)

        do n=1,NDIM
            call fillSyntheticCornerVar(tubeShell%xyz0(n), 10.0_rp*real(n, rp))
            call fillSyntheticCornerVar(tubeShell%X_bmin(n), 100.0_rp*real(n, rp))
            call fillSyntheticCornerVar(tubeShell%X_maxRad(n), 200.0_rp*real(n, rp))
            call fillSyntheticCornerVar(tubeShell%Jxyz_iono(n), 300.0_rp*real(n, rp))
            call fillSyntheticCornerVar(tubeShell%Bxyz_iono(n), 400.0_rp*real(n, rp))
        enddo

        call fillSyntheticCornerVar(tubeShell%J11i, 500.0_rp)
        call fillSyntheticCornerVar(tubeShell%bmin, 600.0_rp)
        call fillSyntheticCornerVar(tubeShell%bVol, 700.0_rp)
        call fillSyntheticCornerVar(tubeShell%Lb, 800.0_rp)
        call fillSyntheticCornerVar(tubeShell%Tb, 900.0_rp)
        call fillSyntheticCornerVar(tubeShell%rCurv, 1000.0_rp)
        call fillSyntheticCornerVar(tubeShell%avgBeta, 1100.0_rp)
        call fillSyntheticCornerVar(tubeShell%losscone, 1200.0_rp)
        call fillSyntheticCornerVar(tubeShell%lossconec, 1300.0_rp)
        call fillSyntheticCornerVar(tubeShell%nTrc, 1400.0_rp)

        do s=0,MAXTUBEFLUIDS
            call fillSyntheticCornerVar(tubeShell%avgP(s), 2000.0_rp + 100.0_rp*real(s, rp))
            call fillSyntheticCornerVar(tubeShell%avgN(s), 3000.0_rp + 100.0_rp*real(s, rp))
            call fillSyntheticCornerVar(tubeShell%stdP(s), 4000.0_rp + 100.0_rp*real(s, rp))
            call fillSyntheticCornerVar(tubeShell%stdN(s), 5000.0_rp + 100.0_rp*real(s, rp))
            call fillSyntheticCornerVar(tubeShell%bcP(s), 6000.0_rp + 100.0_rp*real(s, rp))
            call fillSyntheticCornerVar(tubeShell%bcN(s), 7000.0_rp + 100.0_rp*real(s, rp))
        enddo

        do j=tubeShell%lat0%jsv,tubeShell%lat0%jev
            do i=tubeShell%lat0%isv,tubeShell%lat0%iev
                tubeShell%lat0%data(i,j) = 0.25_rp - 0.002_rp*real(i - tubeShell%lat0%isv, rp)
                tubeShell%lon0%data(i,j) = 0.10_rp + 0.01_rp*real(j - tubeShell%lon0%jsv, rp)
                tubeShell%invlat%data(i,j) = 0.40_rp + 0.001_rp*real(i - tubeShell%invlat%isv, rp)
                tubeShell%latc%data(i,j) = -0.30_rp + 0.0015_rp*real(i - tubeShell%latc%isv, rp) - 0.0005_rp*real(j - tubeShell%latc%jsv, rp)
                tubeShell%lonc%data(i,j) = 0.30_rp + 0.002_rp*real(i - tubeShell%lonc%isv, rp) + 0.01_rp*real(j - tubeShell%lonc%jsv, rp)
                tubeShell%wMAG%data(i,j) = 0.20_rp + 0.001_rp*real(i - tubeShell%wMAG%isv, rp) + 0.002_rp*real(j - tubeShell%wMAG%jsv, rp)
                tubeShell%TioTe0%data(i,j) = 2.0_rp + 0.01_rp*real(i - tubeShell%TioTe0%isv, rp) + 0.005_rp*real(j - tubeShell%TioTe0%jsv, rp)
                tubeShell%topo%data(i,j) = real(TUBE_CLOSED, rp)
            enddo
        enddo

        if (present(firstClosedIO)) then
            do j=tubeShell%topo%jsv,tubeShell%topo%jev
                ! Make a phi-periodic open/closed front with a few harmonics so the
                ! synthetic topology looks more like the real wavering OCB shape.
                phase = modulo(parentSG%ph(j), 2.0_rp*PI)
                frontShift = 2.5_rp*sin(phase - 0.35_rp) + 1.0_rp*sin(3.0_rp*phase + 0.60_rp)
                firstClosedI = nint(real(firstClosedIO, rp) + frontShift)
                firstClosedI = max(tubeShell%topo%isv, min(firstClosedI, tubeShell%topo%iev+1))

                ! Open field is on the poleward/smaller-theta side of the front.
                if (firstClosedI > tubeShell%topo%isv) then
                    tubeShell%topo%data(tubeShell%topo%isv:firstClosedI-1, j) = real(TUBE_OPEN, rp)
                endif
                if (firstClosedI <= tubeShell%topo%iev) then
                    tubeShell%topo%data(firstClosedI:tubeShell%topo%iev, j) = real(TUBE_CLOSED, rp)
                endif
            enddo
        endif

        call applyClosedTubeMasks(tubeShell)
    end subroutine fillSyntheticTubeShell


    subroutine applyRuntimeThetaUOverride(raiCpl)
        !! Reproduce the runtime thetaU logic used in raiCplInitModel so tests can
        !! initialize the child grid from mhdRin without depending on full startup.
        type(raijuCoupler_T), intent(inout) :: raiCpl

        real(rp), dimension(NDIM) :: rxyz
        real(rp) :: thetaU

        rxyz(XDIR) = max(raiCpl%opt%mhdRin, 1.0_rp)
        rxyz(YDIR) = 0.0_rp
        rxyz(ZDIR) = 0.0_rp
        rxyz = DipoleShift(rxyz, raiCpl%opt%voltGrid%radius)
        thetaU = atan2(rxyz(XDIR), rxyz(ZDIR))
        thetaU = min(thetaU, raiCpl%opt%voltGrid%th(raiCpl%opt%voltGrid%Nt/2))
        call raiCpl%raiApp%opt%thetaU%set(thetaU*rad2deg)
    end subroutine applyRuntimeThetaUOverride


    subroutine applyClosedTubeMasks(tubeShell)
        !! Mirror the current production convention that most TubeShell quantities are
        !! only valid for IMAG on closed field lines and therefore share the topo mask.
        type(TubeShell_T), intent(inout) :: tubeShell

        logical, dimension(tubeShell%topo%isv:tubeShell%topo%iev, tubeShell%topo%jsv:tubeShell%topo%jev) :: topoMask
        integer :: n, s

        where (abs(tubeShell%topo%data - real(TUBE_CLOSED, rp)) < TINY)
            topoMask = .true.
        elsewhere
            topoMask = .false.
        end where

        do n=1,NDIM
            tubeShell%X_bmin(n)%mask = topoMask
            tubeShell%X_maxRad(n)%mask = topoMask
        enddo

        tubeShell%latc%mask = topoMask
        tubeShell%lonc%mask = topoMask
        tubeShell%bmin%mask = topoMask
        tubeShell%bVol%mask = topoMask
        tubeShell%Lb%mask = topoMask
        tubeShell%Tb%mask = topoMask
        tubeShell%wMAG%mask = topoMask
        tubeShell%rCurv%mask = topoMask
        tubeShell%avgBeta%mask = topoMask
        do s=0,MAXTUBEFLUIDS
            tubeShell%avgP(s)%mask = topoMask
            tubeShell%avgN(s)%mask = topoMask
            tubeShell%stdP(s)%mask = topoMask
            tubeShell%stdN(s)%mask = topoMask
        enddo
        tubeShell%TioTe0%mask = topoMask
    end subroutine applyClosedTubeMasks

end module raiCplTestHelper
