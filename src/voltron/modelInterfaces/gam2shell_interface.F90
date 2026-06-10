!Gamera innermost J-shell -> voltron ShellGrid packer
!(single global grid; hemisphere logic lives in shell2mix)

module gam2shell_interface
    use kdefs
    use cmidefs,     only: JpSh, JpSt
    use mixdefs,     only: MHDJ, MHDD, MHDC
    use volttypes
    use gamtypes
    use gamutils
    use cmiutils
    use mixgeom,     only: init_grid_fromTP
    use mixinterp,   only: mix_map_grids
    use shellGrid,   only: ShellGrid_T, ShellGridVar_T, initShellVar, SHGR_CORNER
    use shellUtils,  only: wrapJ_SGV
    use shellMixMap, only: mix_to_coords_set_map
    use msphutils,   only: RadIonosphere

    implicit none
    private
    public :: init_gam2shell, convertGameraToShell, mapGameraToShell

    ! Number of mhd2mix variables; kept local for module independence
    integer, parameter :: mhd2mix_varn = 3

    ! Populated in init_gam2shell; used by BIon2Mag.
    real(rp), private :: Rion

    contains

    subroutine init_gam2shell(g2s, gApp, shGrid)
        type(gam2shell_T), intent(inout) :: g2s
        type(gamApp_T),    intent(in)    :: gApp
        type(ShellGrid_T), intent(in)    :: shGrid

        integer :: l, i, j, k, iG
        integer :: is_, ie_, js_, je_
        real(rp) :: xc, yc, zc, xv, yv, zv, thV, phV, s_mhd
        real(rp) :: mhd_Rin
        real(rp), dimension(:,:,:,:), allocatable :: mhdJGrid
        real(rp), dimension(:,:), allocatable :: tFpd, pFpd
        real(rp), dimension(:,:), allocatable :: tFpd_v, pFpd_v
        type(mixGrid_T) :: mhdG

        Rion = RadIonosphere()

        if (.not. allocated(g2s%gJ)) then
            allocate(g2s%gJ(JpSh, gApp%Grid%js:gApp%Grid%je, gApp%Grid%ks:gApp%Grid%ke, NDIM))
        endif
        if (.not. allocated(g2s%gBAvg)) then
            allocate(g2s%gBAvg(JpSh, gApp%Grid%js:gApp%Grid%je, gApp%Grid%ks:gApp%Grid%ke, NDIM))
        endif
        if (.not. allocated(g2s%mhdBuf)) then
            allocate(g2s%mhdBuf(JpSh, gApp%Grid%js:gApp%Grid%je, &
                                gApp%Grid%ks:gApp%Grid%ke, mhd2mix_varn))
        endif
        if (.not. allocated(g2s%JMaps)) then
            allocate(g2s%JMaps(JpSh))
        endif

        allocate(mhdJGrid(JpSh, gApp%Grid%js:gApp%Grid%je, &
                          gApp%Grid%ks:gApp%Grid%ke, 3))

        ! Cartesian cell centers for every J-shell cell (no hemisphere split).
        do k = gApp%Grid%ks, gApp%Grid%ke
            do j = gApp%Grid%js, gApp%Grid%je
                do i = 1, JpSh
                    iG = JpSt + i - 1
                    call cellCenter(gApp%Grid, iG, j, k, xc, yc, zc)
                    g2s%gBAvg(i, j, k, :) = gApp%Grid%B0(iG, j, k, :)
                    mhdJGrid(i, j, k, :) = [xc, yc, zc] / Rion
                enddo
            enddo
        enddo

        !Mean Gamera inner-boundary radius (Rion units) for the dipole trace.
        mhd_Rin = sum(sqrt(mhdJGrid(JpSh,:,:,XDIR)**2 + &
                           mhdJGrid(JpSh,:,:,YDIR)**2 + &
                           mhdJGrid(JpSh,:,:,ZDIR)**2)) / &
                  size(mhdJGrid(JpSh,:,:,XDIR))

        !Trace voltron corners along dipole to MHD inner boundary in flipped frame:
        !sin(theta_MHD) = sqrt(R_in)*sin(theta_iono), phi preserved.
        allocate(tFpd_v(shGrid%Np+1, shGrid%Nt+1))
        allocate(pFpd_v(shGrid%Np+1, shGrid%Nt+1))
        is_ = shGrid%is; ie_ = shGrid%ie
        js_ = shGrid%js; je_ = shGrid%je
        do i = 1, shGrid%Nt+1
            thV = shGrid%th(is_ + i - 1)
            do j = 1, shGrid%Np+1
                phV = shGrid%ph(js_ + j - 1)
                s_mhd = sqrt(mhd_Rin) * sin(thV)
                if (s_mhd >= 1._rp) then
                    ! Ionosphere corner inside the dipole-reach cone: stick
                    ! it at the flipped equator. shell2mix masks these out.
                    xv = cos(phV)
                    yv = sin(phV)
                    zv = 0._rp
                else
                    xv = s_mhd * cos(phV)
                    yv = s_mhd * sin(phV)
                    zv = sign(sqrt(1._rp - s_mhd*s_mhd), cos(thV))
                endif
                tFpd_v(j, i) = acos(xv)
                pFpd_v(j, i) = modulo(atan2(zv, yv), 2._rp * PI)
            enddo
        enddo

        !isPeriodic=.true.: Gamera's native phi IS the flipped pFpd; without it,
        !voltron's equator corners fall in the half-cell seam gap and zero out.
        do l = 1, JpSh
            call mhdgrid_to_flipped(mhdJGrid(l,:,:,:), tFpd, pFpd, mhd_Rin)
            call init_grid_fromTP(mhdG, tFpd, pFpd, &
                                  isSolverGrid=.false., isPeriodic=.true.)
            call mix_to_coords_set_map(mhdG, pFpd_v, tFpd_v, g2s%JMaps(l))
            deallocate(tFpd, pFpd)
        enddo

        deallocate(tFpd_v, pFpd_v)

        !shell2mix consumes this for the per-REMIX-grid boundary mask.
        g2s%mhd_Rin = mhd_Rin

        ! Per-variable, per-hemisphere floor values. Populated every step by
        ! convertGameraToShell; shell2mix consumes to floor masked cells.
        if (.not. allocated(g2s%floorVal)) allocate(g2s%floorVal(mhd2mix_varn, 2))

        ! Output ShellGridVars live on voltron shGrid corners (matches Map dest layout)
        call initShellVar(shGrid, SHGR_CORNER, g2s%fac)
        call initShellVar(shGrid, SHGR_CORNER, g2s%density)
        call initShellVar(shGrid, SHGR_CORNER, g2s%soundSpd)

        deallocate(mhdJGrid)
    end subroutine init_gam2shell


    !Port of convertGameraToRemix unit conversions; no hemisphere branch.
    subroutine convertGameraToShell(g2s, gApp)
        type(gam2shell_T), intent(inout) :: g2s
        class(gamApp_T),   intent(in)    :: gApp

        real(rp) :: Bi2m
        integer  :: i, j, k, iG, v
        integer  :: ksN, keN, ksS, keS
        real(rp) :: xc, yc, zc
        real(rp) :: Con(NVAR)
        real(rp) :: Cs, gB0, gv0, gx0, J11, B0(NDIM)

        gB0 = gApp%Model%Units%gB0
        gv0 = gApp%Model%Units%gv0
        gx0 = gApp%Model%Units%gx0

        call GetShellJ(gApp%Model, gApp%Grid, gApp%State%Bxyz, g2s%gJ)

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,iG,j,k,Bi2m,xc,yc,zc,Con,J11,B0,Cs)
        do k = gApp%Grid%ks, gApp%Grid%ke
            do j = gApp%Grid%js, gApp%Grid%je
                do i = 1, JpSh
                    iG = JpSt + i - 1

                    B0 = gApp%Grid%B0(iG, j, k, :)
                    g2s%gBAvg(i, j, k, :) = B0

                    B0  = g2s%gBAvg(i, j, k, :)
                    J11 = dot_product(g2s%gJ(i, j, k, :), normVec(B0))

                    call cellCenter(gApp%Grid, iG, j, k, xc, yc, zc)
                    Bi2m = BIon2Mag(xc, yc, zc)

                    g2s%mhdBuf(i, j, k, MHDJ) = &
                        J11 * Bi2m * (gB0/gx0 * 1.e4_rp / 4._rp / PI)
                    g2s%mhdBuf(i, j, k, MHDD) = &
                        gApp%State%Gas(iG, j, k, DEN, BLK) * &
                        (gB0 * 1.e-7_rp / gv0)**2 / 4._rp / PI
                    Con = gApp%State%Gas(iG, j, k, :, BLK)
                    call CellPress2Cs(gApp%Model, Con, Cs)
                    g2s%mhdBuf(i, j, k, MHDC) = Cs * gv0 * 1.e2_rp
                enddo
            enddo
        enddo

        ! Per-hemi floors from the last-shell slab (only l=JpSh propagates
        ! through mapGameraToShell). Matches mapGameraToRemix's per-hemi
        ! `minval(mhd2mix%mixInput(l,:,:,v,h))` for l=JpSh.
        ksN = gApp%Grid%ks;           keN = gApp%Grid%ke / 2
        ksS = gApp%Grid%ke / 2 + 1;   keS = gApp%Grid%ke
        do v = 1, mhd2mix_varn
            g2s%floorVal(v, NORTH) = minval(g2s%mhdBuf(JpSh, :, ksN:keN, v))
            g2s%floorVal(v, SOUTH) = minval(g2s%mhdBuf(JpSh, :, ksS:keS, v))
        enddo
    end subroutine convertGameraToShell


    !Apply per-shell map to each MHD variable; only last shell propagates (matches mapGameraToRemix).
    subroutine mapGameraToShell(g2s, shGrid)
        type(gam2shell_T), intent(inout) :: g2s
        type(ShellGrid_T), intent(in)    :: shGrid

        real(rp), dimension(:,:), allocatable :: F
        integer  :: l, v
        integer  :: is_, ie_, js_, je_
        logical, save :: firstCall = .true.

        is_ = shGrid%is; ie_ = shGrid%ie
        js_ = shGrid%js; je_ = shGrid%je

        ! No boundary-floor application here -- shell2mix applies the floor on
        ! REMIX grids using g2s%mhd_Rin and g2s%floorVal.
        do v = 1, mhd2mix_varn
            do l = 1, JpSh
                ! mhdBuf(l,:,:,v) is (j=theta, k=phi); map-apply needs (phi, theta)
                call mix_map_grids(g2s%JMaps(l), &
                                   transpose(g2s%mhdBuf(l, :, :, v)), F)
            enddo

            ! Zero ghosts, write active corners in (theta, phi) layout, wrap phi.
            select case (v)
            case (MHDJ)
                g2s%fac%data = 0._rp
                g2s%fac%data(is_:ie_+1, js_:je_+1) = transpose(F)
                call wrapJ_SGV(shGrid, g2s%fac)
                g2s%fac%mask = .true.
            case (MHDD)
                g2s%density%data = 0._rp
                g2s%density%data(is_:ie_+1, js_:je_+1) = transpose(F)
                call wrapJ_SGV(shGrid, g2s%density)
                g2s%density%mask = .true.
            case (MHDC)
                g2s%soundSpd%data = 0._rp
                g2s%soundSpd%data(is_:ie_+1, js_:je_+1) = transpose(F)
                call wrapJ_SGV(shGrid, g2s%soundSpd)
                g2s%soundSpd%mask = .true.
            end select
        enddo

        if (allocated(F)) deallocate(F)

        if (firstCall) then
            write(*,'(A)') '[gam2shell DIAG] First-call stats (voltron shGrid):'
            write(*,'(A,1P,2E14.5,A,2E14.5)') &
                '  mhdBuf(JpSh,:,:,MHDJ) min/max=', &
                minval(g2s%mhdBuf(JpSh,:,:,MHDJ)), maxval(g2s%mhdBuf(JpSh,:,:,MHDJ)), &
                '   on shell fac min/max=', &
                minval(g2s%fac%data), maxval(g2s%fac%data)
            write(*,'(A,1P,2E14.5,A,2E14.5)') &
                '  mhdBuf(JpSh,:,:,MHDD) min/max=', &
                minval(g2s%mhdBuf(JpSh,:,:,MHDD)), maxval(g2s%mhdBuf(JpSh,:,:,MHDD)), &
                '   on shell den min/max=', &
                minval(g2s%density%data), maxval(g2s%density%data)
            write(*,'(A,1P,2E14.5,A,2E14.5)') &
                '  mhdBuf(JpSh,:,:,MHDC) min/max=', &
                minval(g2s%mhdBuf(JpSh,:,:,MHDC)), maxval(g2s%mhdBuf(JpSh,:,:,MHDC)), &
                '   on shell cs  min/max=', &
                minval(g2s%soundSpd%data), maxval(g2s%soundSpd%data)
            write(*,'(A,1P,E14.5)') '  mhd_Rin = ', g2s%mhd_Rin
            write(*,'(A,1P,6E12.4)') '  floorVal(MHDJ,MHDD,MHDC) N,S = ', &
                g2s%floorVal(MHDJ,NORTH), g2s%floorVal(MHDJ,SOUTH), &
                g2s%floorVal(MHDD,NORTH), g2s%floorVal(MHDD,SOUTH), &
                g2s%floorVal(MHDC,NORTH), g2s%floorVal(MHDC,SOUTH)
            write(*,'(A,4I6)') '  shGrid is/ie js/je = ', is_, ie_, js_, je_
            write(*,'(A,2I6)') '  shGrid Nt Np = ', shGrid%Nt, shGrid%Np
            write(*,'(A,1P,2E14.5)') '  shGrid th(is) th(ie+1) = ', &
                shGrid%th(is_), shGrid%th(ie_+1)
            firstCall = .false.
        endif
    end subroutine mapGameraToShell


    !Cartesian (Rion units) -> flipped-frame (tFpd, pFpd) with x-axis polar;
    !inner J-shell becomes a disk, no hole. Output in mix layout (phi, theta).
    subroutine mhdgrid_to_flipped(mhdg, tFpd, pFpd, Rinner)
        real(rp), dimension(:,:,:), intent(in) :: mhdg   ! (nj, nk, 3)
        real(rp), dimension(:,:), allocatable, intent(out) :: tFpd, pFpd
        real(rp), intent(out) :: Rinner

        real(rp), dimension(:,:), allocatable :: r
        integer :: nj, nk

        nj = size(mhdg, 1); nk = size(mhdg, 2)
        allocate(r(nj, nk))
        allocate(tFpd(nk, nj), pFpd(nk, nj))

        r = sqrt(mhdg(:,:,XDIR)**2 + mhdg(:,:,YDIR)**2 + mhdg(:,:,ZDIR)**2)

        tFpd = transpose(acos(mhdg(:,:,XDIR) / r))
        pFpd = transpose(modulo(atan2(mhdg(:,:,ZDIR), mhdg(:,:,YDIR)), 2._rp * PI))

        Rinner = sum(r) / size(r)

        deallocate(r)
    end subroutine mhdgrid_to_flipped


    !Ion->Mag dipole field scaling, xyz in Rion units.
    !Verbatim copy from mhd2mix_interface.F90 (avoids cross-module Rion dep).
    function BIon2Mag(x, y, z) result(Bi2m)
        real(rp), intent(in) :: x, y, z
        real(rp) :: Bi2m

        real(rp) :: rMag, rRatio, zor2, mZ, pZ

        rMag   = norm2([x, y, z])
        rRatio = rMag / Rion
        zor2   = (z / rMag)**2
        mZ     = (1._rp - zor2) / rRatio
        pZ     = 1._rp + 3._rp * zor2

        Bi2m = (rRatio**3) * sqrt(1._rp + 3._rp * (1._rp - mZ)) / sqrt(pZ)
    end function BIon2Mag

end module gam2shell_interface
