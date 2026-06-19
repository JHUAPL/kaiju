! Routines bridging mixGrid_T and ShellGrid_T, mirroring mix_set_map /
! mix_interpolant / mix_search. Output Map_T matches mix_set_map's, so
! mix_map_grids applies it unchanged. Map storage stays in mix (phi,theta)
! layout; transposes at ShellGridVar boundaries are the caller's job.

module shellMixMap
    use kdefs
    use mixtypes,  only: mixGrid_T, Map_T
    use shellGrid, only: ShellGrid_T
    use math,      only: matinv4
    use mixinterp, only: mix_search

    implicit none
    private
    public :: shell_to_mix_set_map, mix_to_coords_set_map

    contains

    ! Bisection-ish locate for a ShellGrid active corner range.
    ! Mirrors mix_search but uses 1D shGr%th / shGr%ph.
    subroutine shell_search(shGr, p, t, j, i)
        type(ShellGrid_T), intent(in) :: shGr
        real(rp), intent(in) :: p, t
        integer, intent(out) :: j, i
        integer :: is_, ie_, js_, je_

        is_ = shGr%is; ie_ = shGr%ie
        js_ = shGr%js; je_ = shGr%je

        i = minloc(abs(shGr%th(is_:ie_+1) - t), 1) + is_ - 1
        j = minloc(abs(shGr%ph(js_:je_+1) - p), 1) + js_ - 1

        ! Use TINY tolerance, not strict > 0: shGr corners carry round-trip
        ! noise (th(1) ~4.5e-17, not exact 0), so a strict test at the pole
        ! decrements i to 0 and fills row 0 with dead cells (pole artifact).
        if ((shGr%th(i) - t) > TINY) i = i - 1
        if ((shGr%ph(j) - p) > TINY) j = j - 1

        ! Snap upper-corner hit back into the last active cell. Without this,
        ! t = shGr%th(ie_+1) exactly (e.g. SOUTH pole via virtG = PI - 0)
        ! leaves i at ie_+1, which the caller rejects as out-of-range.
        if (i > ie_) i = ie_

        ! periodic wrap on phi
        if (shGr%isPeriodic .and. (p - shGr%ph(je_+1)) > 0 .and. p < 2*PI) j = je_+1
    end subroutine shell_search


    ! Precompute the 4x4 inverse Vandermonde per active phi-theta cell of
    ! shGr. Mirror of mix_interpolant. Output layout matches mix:
    !   Interpolant(Np+1, Nt+1, 4, 4)   — shGr-cell index (j, i)
    ! Poles are left as zero (same convention as mix_interpolant).
    subroutine shell_interpolant(shGr, Interpolant)
        type(ShellGrid_T), intent(in) :: shGr
        real(rp), dimension(:,:,:,:), allocatable, intent(out) :: Interpolant
        real(rp), dimension(4,4) :: A, Ainv
        real(rp) :: dp, thI, thIp1, phJ, phJp1
        integer :: i, j, jp1
        integer :: is_, ie_, js_, je_

        is_ = shGr%is; ie_ = shGr%ie
        js_ = shGr%js; je_ = shGr%je

        allocate(Interpolant(shGr%Np+1, shGr%Nt+1, 4, 4))
        Interpolant = 0.0_rp

        ! Cells: 1..Nt in theta, 1..Np in phi (with wrap at last phi if periodic)
        do i = 1, shGr%Nt
            thI   = shGr%th(is_ + i - 1)
            thIp1 = shGr%th(is_ + i)
            do j = 1, shGr%Np
                if (j == shGr%Np .and. .not. shGr%isPeriodic) cycle

                if (j /= shGr%Np) then
                    jp1 = j + 1
                    dp  = 0.0_rp
                else
                    jp1 = 1
                    dp  = 2*PI
                endif
                phJ   = shGr%ph(js_ + j - 1)
                phJp1 = shGr%ph(js_ + jp1 - 1) + dp

                A = reshape( &
                    (/ 1.0_rp, phJ,   thI,   phJ   * thI,   &
                       1.0_rp, phJp1, thI,   phJp1 * thI,   &
                       1.0_rp, phJ,   thIp1, phJ   * thIp1, &
                       1.0_rp, phJp1, thIp1, phJp1 * thIp1 /), shape(A))

                ! Skip poles (matches mix_interpolant)
                if (.not. ( &
                     (i == 1         .and. thI        <= TINY) .or. &
                     (i == shGr%Nt   .and. (PI - thIp1) <= TINY) &
                     )) then
                    call matinv4(A, Ainv)
                    Interpolant(j, i, :, :) = Ainv
                endif
            enddo
        enddo
    end subroutine shell_interpolant


    ! Builds a mix-source bilinear map, but destination coords are supplied as
    ! explicit 2D arrays (destP, destT) in mix layout (Np_dst, Nt_dst)
    ! instead of being read from a ShellGrid_T. Used by gam2shell to
    ! drive the interpolation in Gamera's flipped (x-polar) frame: the
    ! Gamera source G1 is built in flipped coords, and destP/destT are
    ! the voltron corners converted to the same flipped frame.
    subroutine mix_to_coords_set_map(G1, destP, destT, Map)
        type(mixGrid_T),          intent(in)  :: G1
        real(rp), dimension(:,:), intent(in)  :: destP, destT
        type(Map_T),              intent(out) :: Map
        real(rp), dimension(4) :: x
        real(rp) :: p2, t2
        integer :: i1, j1, i2, j2
        integer :: Np_dst, Nt_dst

        Np_dst = size(destP, 1)
        Nt_dst = size(destP, 2)

        if (.not. allocated(Map%M))  allocate(Map%M (Np_dst, Nt_dst, 4))
        if (.not. allocated(Map%I1)) allocate(Map%I1(Np_dst, Nt_dst))
        if (.not. allocated(Map%J1)) allocate(Map%J1(Np_dst, Nt_dst))

        do i2 = 1, Nt_dst
            do j2 = 1, Np_dst
                p2 = destP(j2, i2)
                t2 = destT(j2, i2)

                call mix_search(G1, p2, t2, j1, i1)

                ! Periodic wrap at the low-p edge. mix_search returns j1=0
                ! when p2 is below the first cell center; for a periodic
                ! grid the cell straddling p2 in [0, G%p(1)) is cell Np,
                ! whose interpolant was built with jp1=1, dp=2*pi. Re-anchor
                ! j1=Np and shift the query by 2*pi so it lands inside the
                ! interpolant's [G%p(Np), G%p(1)+2*pi] support.
                if (G1%isPeriodic .and. j1 == 0) then
                    j1 = G1%Np
                    p2 = p2 + 2.0_rp * PI
                endif

                Map%I1(j2, i2) = i1
                Map%J1(j2, i2) = j1
                x = (/ 1.0_rp, p2, t2, p2*t2 /)

                if (i1 == 0 .or. i1 == G1%Nt) then
                    Map%M(j2, i2, :) = 0.0_rp
                    cycle
                endif
                if (.not. G1%isPeriodic .and. (j1 == 0 .or. j1 == G1%Np)) then
                    Map%M(j2, i2, :) = 0.0_rp
                    cycle
                endif

                if ((i1 == 1) .and. (G1%t(1, 1) < TINY)) then
                    Map%M(j2, i2, :) = 1.0_rp/3.0_rp * (/ 0.5_rp, 0.5_rp, 1.0_rp, 1.0_rp /)
                else if ((i1 == G1%Nt - 1) .and. ((PI - G1%t(1, G1%Nt-1)) < TINY)) then
                    Map%M(j2, i2, :) = 1.0_rp/3.0_rp * (/ 1.0_rp, 1.0_rp, 0.5_rp, 0.5_rp /)
                else
                    Map%M(j2, i2, :) = matmul(G1%Interpolant(j1, i1, :, :), x)
                endif
            enddo
        enddo

        Map%M = max(Map%M, 0.0_rp)
        Map%M = min(Map%M, 1.0_rp)
    end subroutine mix_to_coords_set_map


    ! Build bilinear map from shGr source onto mix destination G2 corner
    ! points. Shell-side interpolant computed here and discarded after use.
    subroutine shell_to_mix_set_map(shGr, G2, Map)
        type(ShellGrid_T), intent(in)  :: shGr
        type(mixGrid_T),   intent(in)  :: G2
        type(Map_T),       intent(out) :: Map
        real(rp), dimension(:,:,:,:), allocatable :: Interp
        real(rp), dimension(4) :: x
        real(rp) :: p2, t2
        integer :: i1, j1, i2, j2
        integer :: is_, ie_, js_, je_

        is_ = shGr%is; ie_ = shGr%ie
        js_ = shGr%js; je_ = shGr%je

        call shell_interpolant(shGr, Interp)

        if (.not. allocated(Map%M))  allocate(Map%M (G2%Np, G2%Nt, 4))
        if (.not. allocated(Map%I1)) allocate(Map%I1(G2%Np, G2%Nt))
        if (.not. allocated(Map%J1)) allocate(Map%J1(G2%Np, G2%Nt))

        do i2 = 1, G2%Nt
            do j2 = 1, G2%Np
                t2 = G2%t(j2, i2)
                p2 = G2%p(j2, i2)

                call shell_search(shGr, p2, t2, j1, i1)
                ! shell_search returns active-frame corners, but mix_map_grids
                ! reads the ghost-inclusive transpose(SGV%data) whose index 1 is
                ! the first west ghost. Rebase the stored indices to that frame
                ! (subtract the ghost start); otherwise every cell reads Ngw deg
                ! west, a rigid azimuthal rotation. Local i1/j1 stay active-frame
                ! for the interpolant and boundary tests below.
                Map%I1(j2, i2) = i1 - shGr%isg + 1
                Map%J1(j2, i2) = j1 - shGr%jsg + 1
                x = (/ 1.0_rp, p2, t2, p2*t2 /)

                ! outside source theta range (using shGr active corners)
                if (i1 < is_ .or. i1 > ie_) then
                    Map%M(j2, i2, :) = 0.0_rp
                    cycle
                endif

                ! non-periodic phi edge on shell source
                if (.not. shGr%isPeriodic .and. (j1 < js_ .or. j1 > je_)) then
                    Map%M(j2, i2, :) = 0.0_rp
                    cycle
                endif

                ! poles
                if ((i1 == is_) .and. (shGr%th(is_) < TINY)) then
                    Map%M(j2, i2, :) = 1.0_rp/3.0_rp * (/ 0.5_rp, 0.5_rp, 1.0_rp, 1.0_rp /)
                else if ((i1 == ie_) .and. ((PI - shGr%th(ie_+1)) < TINY)) then
                    Map%M(j2, i2, :) = 1.0_rp/3.0_rp * (/ 1.0_rp, 1.0_rp, 0.5_rp, 0.5_rp /)
                else
                    ! Interpolant is indexed from 1 in both axes; translate to cell index
                    Map%M(j2, i2, :) = matmul(Interp(j1 - js_ + 1, i1 - is_ + 1, :, :), x)
                endif
            enddo
        enddo

        Map%M = max(Map%M, 0.0_rp)
        Map%M = min(Map%M, 1.0_rp)

        deallocate(Interp)
    end subroutine shell_to_mix_set_map

end module shellMixMap
