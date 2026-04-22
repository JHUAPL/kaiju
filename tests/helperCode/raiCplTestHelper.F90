module raiCplTestHelper

    use kdefs, only : rp, PI, strLen

    use shellGrid
    use volttypes, only : raijuCoupler_T, imagOptions_T, IMAG_nMHD_R_NEAREARTH

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

end module raiCplTestHelper
