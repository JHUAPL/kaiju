!remixCoupler_T procedures

submodule (volttypes) remixCplTypesSub

    use mixtypes
    use mixdefs
    use mhd2mix_interface
    use mix2mhd_interface
    use imag2mix_interface
    use mix2shell_interface,  only: mapMixToShell
    use voltappHelper
    use kai2geo
    use mixmain,   only: updateF107
    use shellutils,   only: wrapJ_SGV
    use planethelper, only: DP2Cs
    use kdefs,        only: Mp_cgs, TINY, PI
    use gdefs,        only: BLK

    implicit none

    contains

    module subroutine remixCplInitModel(App, xml)
        class(remixCoupler_T), intent(inout) :: App
        type(XML_Input_T),     intent(inout) :: xml

        call App%mixApp%InitModel(xml)

    end subroutine remixCplInitModel


    module subroutine remixCpl_InitCouplers(App, vApp, xml)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp
        type(XML_Input_T),     intent(inout) :: xml

        type(TimeSeries_T)    :: f107
        character(len=strLen) :: RunID, resID
        logical               :: isRestart
        real(rp)              :: maxF107
        integer               :: nRes

        isRestart = vApp%gApp%Model%isRestart
        RunID     = trim(vApp%gApp%Model%RunID)

        if (isRestart) then
            call App%ReadRestart(resId=RunID, nRes=vApp%IO%nRes)
        endif

        if (vApp%writeFiles) then
            call App%InitIO(xml)
        endif

        App%mixApp%ion%rad_iono_m   = vApp%planet%ri_m
        App%mixApp%ion%rad_planet_m = vApp%planet%rp_m

        ! Ensure remix and voltron restart numbers match
        if (isRestart .and. vApp%IO%nRes /= App%mixApp%ion(1)%P%nRes) then
            write(*,*) "Voltron and Remix disagree on restart number, you should sort that out."
            write(*,*) "Error code: A house divided cannot stand"
            write(*,*) "   Voltron nRes = ", vApp%IO%nRes
            write(*,*) "   Remix   nRes = ", App%mixApp%ion(1)%P%nRes
            stop
        endif

        if (App%mixApp%ion(NORTH)%P%doSWF107 .neqv. App%mixApp%ion(SOUTH)%P%doSWF107) then
            write(*,*) 'Something is wrong. doSWf107 is set differently for the two hemispheres.'
            write(*,*) 'Exiting ...'
            stop
        endif

        ! read f107 from the SW file and overwrite what's been read from .xml above (in init_mix)
        ! note, only checking for NORTH, because both hemispheres read the same xml file
        if (App%mixApp%ion(NORTH)%P%doSWF107) then
            f107%wID = vApp%tilt%wID
            call f107%initTS("f10.7", doLoudO=.false.)
            maxF107 = f107%getMax()
            call updateF107(App%mixApp%ion, maxF107)
            vApp%dkApp%params%f107 = maxF107
            if (vApp%isLoud) write(*,*) 'Using F10.7 = ', maxF107
        endif

        if (vApp%isLoud) write(*,*) 'Using MJD0  = ', vApp%gApp%Model%MJD0

        call init_mix2MhdCoupler(vApp%gApp, App%mixApp)
        call vApp%gApp%InitMhdCoupler(vApp)
        if (isRestart) then
            call xml%Set_Val(resID, "/Kaiju/gamera/restart/resID", "msphere")
            call xml%Set_Val(nRes,  "/Kaiju/gamera/restart/nRes", -1)
            call vApp%gApp%ReadRestart(resID, nRes)
        endif

    end subroutine remixCpl_InitCouplers


    module subroutine remixCpl_toIono(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        real(rp) :: curTilt

        call MJDRecalc(vApp%MJD)

        call App%packInputs(vApp)

        call vApp%tilt%getValue(vApp%time, curTilt)
        App%mixApp%IO%tilt = curTilt
        App%mixApp%IO%time = vApp%time
        App%mixApp%IO%mjd  = vApp%MJD

        call mapDKToMix(vApp%shGrid, vApp%dkApp%vars, App%mixApp)
        if (vApp%doDeep .and. (vApp%time > 0)) call mapNonDKToMix(vApp, App%mixApp)

        ! GCM conductance overwrites DK when enabled
        if (vApp%doGCM) call mapGCMToMix(vApp%gcm, App%mixApp)

    end subroutine remixCpl_toIono


    module subroutine remixCpl_fromIono(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        call mapMixToShell(App%mixApp, vApp)

        ! Corner potential for DK doET
        call mixVarToVoltron(App%mixApp, POT, vApp%shGrid, vApp%pot_iono)

    end subroutine remixCpl_fromIono


    module subroutine remixCplAdvanceModel(App, dt)
        class(remixCoupler_T), intent(inout) :: App
        real(rp),              intent(in)    :: dt

        call App%mixApp%AdvanceModel(dt)

    end subroutine remixCplAdvanceModel


    module subroutine remixCplInitIO(App, xml)
        class(remixCoupler_T), intent(inout) :: App
        type(XML_Input_T),     intent(inout) :: xml

        call App%mixApp%InitIO(xml)

    end subroutine remixCplInitIO

    module subroutine remixCplWriteRestart(App, nRes)
        class(remixCoupler_T), intent(inout) :: App
        integer,               intent(in)    :: nRes

        call App%mixApp%WriteRestart(nRes)

    end subroutine remixCplWriteRestart

    module subroutine remixCplReadRestart(App, resId, nRes)
        class(remixCoupler_T), intent(inout) :: App
        character(len=*),      intent(in)    :: resId
        integer,               intent(in)    :: nRes

        call App%mixApp%ReadRestart(resId, nRes)

    end subroutine remixCplReadRestart

    module subroutine remixCplWriteConsoleOutput(App)
        class(remixCoupler_T), intent(inout) :: App

        call App%mixApp%WriteConsoleOutput()

    end subroutine remixCplWriteConsoleOutput

    module subroutine remixCplWriteFileOutput(App, nStep)
        class(remixCoupler_T), intent(inout) :: App
        integer,               intent(in)    :: nStep

        call App%mixApp%WriteFileOutput(nStep)

    end subroutine remixCplWriteFileOutput

    module subroutine remixCplWriteSlimFileOutput(App, nStep)
        class(remixCoupler_T), intent(inout) :: App
        integer,               intent(in)    :: nStep

        call App%mixApp%WriteSlimFileOutput(nStep)

    end subroutine remixCplWriteSlimFileOutput

    module subroutine remixCplCleanup(App)
        class(remixCoupler_T), intent(inout) :: App

        call App%mixApp%Cleanup()

    end subroutine remixCplCleanup


    !packInputs: tubeShell (J11i/bcN/bcP) -> REMIX hemisphere grids
    module subroutine remix_packInputs(App, vApp)
        class(remixCoupler_T), intent(inout) :: App
        class(voltApp_T),      intent(inout) :: vApp

        type(ShellGridVar_T)                 :: facV, denV, csV
        logical, dimension(:,:), allocatable :: mask
        real(rp) :: floorD(2), floorCs(2)
        logical  :: doDeepFeed
        integer  :: i, j, h, Nth, is_, ie_, js_, je_

        doDeepFeed = (vApp%doDeep .and. vApp%time >= 0)

        call initShellVar(vApp%shGrid, SHGR_CORNER, facV)
        call initShellVar(vApp%shGrid, SHGR_CORNER, denV)
        call initShellVar(vApp%shGrid, SHGR_CORNER, csV)

        facV%data = vApp%State%tubeShell%J11i%data                  ! uA/m^2
        denV%data = vApp%State%tubeShell%bcN(BLK)%data * Mp_cgs      ! #/cc -> g/cm^3
        do j = lbound(csV%data,2), ubound(csV%data,2)
            do i = lbound(csV%data,1), ubound(csV%data,1)
                csV%data(i,j) = DP2Cs(vApp%State%tubeShell%bcN(BLK)%data(i,j), &
                                      vApp%State%tubeShell%bcP(BLK)%data(i,j)) * 1.0e5_rp  ! km/s -> cm/s
            enddo
        enddo
        ! bcN/bcP = 0 marks unsampled tubes; mask them, FAC 0 is physical
        facV%mask = .true.
        denV%mask = (denV%data > 0._rp)
        csV%mask  = (vApp%State%tubeShell%bcN(BLK)%data > 0._rp) .and. &
                    (vApp%State%tubeShell%bcP(BLK)%data > 0._rp)

        call mapShellVarToMix(vApp%shGrid, facV, App%mixApp%ion, FAC)
        call mapShellVarToMix(vApp%shGrid, denV, App%mixApp%ion, DENSITY)
        call mapShellVarToMix(vApp%shGrid, csV,  App%mixApp%ion, SOUND_SPEED)

        ! Per-hemisphere D/Cs floor over sampled corners
        is_ = vApp%shGrid%is; ie_ = vApp%shGrid%ie
        js_ = vApp%shGrid%js; je_ = vApp%shGrid%je
        floorD = 0.0_rp; floorCs = 0.0_rp
        if (doDeepFeed) then
            do h = NORTH, SOUTH
                floorD(h)  = 1.0e30_rp
                floorCs(h) = 1.0e30_rp
                do j = js_, je_+1
                    do i = is_, ie_+1
                        if (sqrt(App%mhd_Rin)*sin(vApp%shGrid%th(i)) > 1._rp) cycle
                        if ((h == NORTH) .neqv. (cos(vApp%shGrid%th(i)) > 0._rp)) cycle
                        if (denV%data(i,j) <= 0._rp) cycle
                        floorD(h)  = min(floorD(h),  denV%data(i,j))
                        floorCs(h) = min(floorCs(h), csV%data(i,j))
                    enddo
                enddo
                if (floorD(h)  >= 1.0e30_rp) floorD(h)  = 0.0_rp
                if (floorCs(h) >= 1.0e30_rp) floorCs(h) = 0.0_rp
            enddo
        endif

        ! Equatorward of the MHD inner boundary: FAC->0, D/Cs->per-hemi floor.
        do h = NORTH, SOUTH
            Nth = App%mixApp%ion(h)%G%Nt
            allocate(mask(App%mixApp%ion(h)%G%Np, App%mixApp%ion(h)%G%Nt))
            mask = (sqrt(App%mhd_Rin) * sin(App%mixApp%ion(h)%G%t) > 1._rp - TINY)
            where (mask) App%mixApp%ion(h)%St%Vars(:,:,FAC)         = 0._rp
            where (mask) App%mixApp%ion(h)%St%Vars(:,:,DENSITY)     = floorD(h)
            where (mask) App%mixApp%ion(h)%St%Vars(:,:,SOUND_SPEED) = floorCs(h)
            App%mixApp%ion(h)%St%Vars(:, Nth-1:Nth, FAC) = 0._rp
            deallocate(mask)
        enddo

    end subroutine remix_packInputs


    !DK vars (cell-centered, active) -> both REMIX hemispheres
    subroutine mapDKToMix(voltGrid, dkVars, mixApp)
        type(ShellGrid_T),   intent(in)    :: voltGrid
        type(dkVariables_T), intent(in)    :: dkVars
        class(mixApp_T),     intent(inout) :: mixApp
        type(ShellGridVar_T) :: sgv
        integer :: is_, ie_, js_, je_

        is_ = voltGrid%is; ie_ = voltGrid%ie
        js_ = voltGrid%js; je_ = voltGrid%je

        call mapOneDK(dkVars%SigmaP%data,   SIGMAP)
        call mapOneDK(dkVars%SigmaH%data,   SIGMAH)
        call mapOneDK(dkVars%Avg_Eng%data,  AVG_ENG)
        call mapOneDK(dkVars%Num_Flux%data, NUM_FLUX)
        call mapOneDK(dkVars%DeltaE%data,   DELTAE)
        call mapOneDK(dkVars%Z_Eavg%data,   Z_EAVG)
        call mapOneDK(dkVars%Z_Nflux%data,  Z_NFLUX)
        call mapOneDK(dkVars%Im_Beta%data,  IM_BETA)
        call mapOneDK(dkVars%Im_Eden%data,  IM_EDEN)
        call mapOneDK(dkVars%Im_Epre%data,  IM_EPRE)
        call mapOneDK(dkVars%Im_Enflx%data, IM_ENFLX)
        call mapOneDK(dkVars%Im_Eflux%data, IM_EFLUX)
        call mapOneDK(dkVars%Im_Npsp%data,  IM_NPSP)
        call mapOneDK(dkVars%Im_Eavg%data,  IM_EAVG)
        call mapOneDK(dkVars%Im_Gtype%data, IM_GTYPE)
        call mapOneDKcat(dkVars%Aur_Type%data, AUR_TYPE)

    contains

        subroutine mapOneDK(dkField, mixID)
            real(rp), dimension(:,:), intent(in) :: dkField
            integer,                  intent(in) :: mixID
            integer :: i, j
            call initShellVar(voltGrid, SHGR_CC, sgv)
            do i = is_, ie_
                do j = js_, je_
                    sgv%data(i,j) = dkField(i - is_ + 1, j - js_ + 1)
                enddo
            enddo
            call wrapJ_SGV(voltGrid, sgv)
            sgv%mask = .true.
            call mapShellVarToMix(voltGrid, sgv, mixApp%ion, mixID)
        end subroutine mapOneDK

        ! Nearest-neighbour for categorical fields
        subroutine mapOneDKcat(dkField, mixID)
            real(rp), dimension(:,:), intent(in) :: dkField
            integer,                  intent(in) :: mixID
            integer :: i, j
            call initShellVar(voltGrid, SHGR_CC, sgv)
            do i = is_, ie_
                do j = js_, je_
                    sgv%data(i,j) = dkField(i - is_ + 1, j - js_ + 1)
                enddo
            enddo
            call wrapJ_SGV(voltGrid, sgv)
            sgv%mask = .true.
            call mapShellVarToMixNN(voltGrid, sgv, mixApp%ion, mixID)
        end subroutine mapOneDKcat

    end subroutine mapDKToMix

    subroutine mapNonDKToMix(vApp, mixApp)
        class(voltApp_T), intent(inout) :: vApp
        class(mixApp_T),  intent(inout) :: mixApp

        type(ShellGridVar_T) :: sinLon, cosLon
        integer  :: i, j, jM
        real(rp) :: x1, x2, latc_v, lonc_v, s, c

        ! No proton source yet
        mixApp%ion(NORTH)%St%Vars(:,:,IM_IFLUX) = 0.0_rp
        mixApp%ion(NORTH)%St%Vars(:,:,IM_IAVG ) = 0.0_rp
        mixApp%ion(NORTH)%St%Vars(:,:,IM_INFLX) = 0.0_rp
        mixApp%ion(SOUTH)%St%Vars(:,:,IM_IFLUX) = 0.0_rp
        mixApp%ion(SOUTH)%St%Vars(:,:,IM_IAVG ) = 0.0_rp
        mixApp%ion(SOUTH)%St%Vars(:,:,IM_INFLX) = 0.0_rp

        ! Conjugate lat/lon from tubeShell; lon via sin/cos to avoid the wrap
        sinLon = vApp%State%tubeShell%lonc; sinLon%data = sin(vApp%State%tubeShell%lonc%data)
        cosLon = vApp%State%tubeShell%lonc; cosLon%data = cos(vApp%State%tubeShell%lonc%data)

        associate(GN => mixApp%ion(NORTH)%G, GS => mixApp%ion(SOUTH)%G, &
                  ion => mixApp%ion, sh => vApp%shGrid, ts => vApp%State%tubeShell)
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,jM,x1,x2,latc_v,lonc_v,s,c)
            do j = 1, GN%Np
                do i = 1, GN%Nt
                    x1 = GN%t(j,i); x2 = GN%p(j,i)
                    call InterpShellVar_TSC_pnt(sh, ts%latc, x1, x2, latc_v)
                    call InterpShellVar_TSC_pnt(sh, sinLon,  x1, x2, s)
                    call InterpShellVar_TSC_pnt(sh, cosLon,  x1, x2, c)
                    lonc_v = atan2(s, c); if (lonc_v < 0.0_rp) lonc_v = lonc_v + 2.0_rp*PI
                    ion(NORTH)%St%Vars(j,i,IM_THCON) = PI/2.0_rp - latc_v
                    ion(NORTH)%St%Vars(j,i,IM_PHCON) = lonc_v

                    ! SH: reflect theta about the equator
                    jM = GS%Np + 1 - j
                    x1 = PI - GS%t(jM,i); x2 = GS%p(jM,i)
                    call InterpShellVar_TSC_pnt(sh, ts%latc, x1, x2, latc_v)
                    call InterpShellVar_TSC_pnt(sh, sinLon,  x1, x2, s)
                    call InterpShellVar_TSC_pnt(sh, cosLon,  x1, x2, c)
                    lonc_v = atan2(s, c); if (lonc_v < 0.0_rp) lonc_v = lonc_v + 2.0_rp*PI
                    ion(SOUTH)%St%Vars(j,i,IM_THCON) = PI/2.0_rp - latc_v
                    ion(SOUTH)%St%Vars(j,i,IM_PHCON) = lonc_v
                end do
            end do
        end associate

    end subroutine mapNonDKToMix

    subroutine mapGCMToMix(gcm, mixApp)
        class(gcmCoupler_T),    intent(inout) :: gcm
        class(mixApp_T),intent(inout) :: mixApp

        if (gcm%GEO%gcm2mix_nvar > 0) call mapOneGCMGridToMix(gcm%GEO, mixApp)
        if (gcm%APEX%gcm2mix_nvar > 0) call mapOneGCMGridToMix(gcm%APEX, mixApp)
    end subroutine mapGCMToMix

    subroutine mapOneGCMGridToMix(gcmGrid, mixApp)
        type(gcm_grid_T), intent(inout) :: gcmGrid
        class(mixApp_T),  intent(inout) :: mixApp

        integer :: v, h, mixVar

        if (gcmGrid%gcm2mix_nvar == 0) return

        do v = 1, gcmGrid%gcm2mix_nvar

            select case (trim(gcmGrid%inlist(v)))
            case ("SIGMAP")
                mixVar = SIGMAP
            case ("SIGMAH")
                mixVar = SIGMAH
            case default
                write(*,*) "MIXCPL WARNING: Unknown GCM inlist variable: ", trim(gcmGrid%inlist(v))
                cycle
            end select

            call mapShellVarToMix(gcmGrid%grid, gcmGrid%gcmInput(v), mixApp%ion, mixVar)
            do h = NORTH, SOUTH
                select case (mixVar)
                case (SIGMAP)
                    mixApp%ion(h)%St%Vars(:,:,SIGMAP) = max(mixApp%ion(h)%P%pedmin, &
                                                             mixApp%ion(h)%St%Vars(:,:,SIGMAP))
                case (SIGMAH)
                    mixApp%ion(h)%St%Vars(:,:,SIGMAH) = max(mixApp%ion(h)%P%hallmin, &
                                                             mixApp%ion(h)%St%Vars(:,:,SIGMAH))
                end select
            enddo
            write(*,*) "Mapped GCM variable ", trim(gcmGrid%inlist(v)), " to REMIX variable ", mixVar
        end do

        write(*,*) "Val Check NORTH: ", maxval(mixApp%ion(NORTH)%St%Vars(:,:,SIGMAP)), minval(mixApp%ion(NORTH)%St%Vars(:,:,SIGMAP))
        write(*,*) "Val Check SOUTH: ", maxval(mixApp%ion(SOUTH)%St%Vars(:,:,SIGMAP)), minval(mixApp%ion(SOUTH)%St%Vars(:,:,SIGMAP))
    end subroutine mapOneGCMGridToMix
end submodule remixCplTypesSub
