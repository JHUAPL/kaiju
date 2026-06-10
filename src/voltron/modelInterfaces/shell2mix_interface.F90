!Voltron ShellGrid -> REMIX hemisphere-grid packer.
!NORTH map is direct; SOUTH uses a virtual grid mirroring SOUTH cells into
!voltron's frame (PI-t, 2PI-p) so mix_map_grids output lands in St%Vars
!without reshuffling. Applies inner-boundary mask: FAC->0, D/Cs->floorVal.

module shell2mix_interface
    use kdefs
    use mixdefs,     only: MHDJ, MHDD, MHDC, FAC, DENSITY, SOUND_SPEED
    use mixtypes,    only: mixGrid_T, Map_T
    use volttypes,   only: gam2shell_T, shell2mix_T, remixCoupler_T, mixApp_T
    use mixinterp,   only: mix_map_grids
    use shellGrid,   only: ShellGrid_T
    use shellMixMap, only: shell_to_mix_set_map

    implicit none
    private
    public :: init_shell2mix, mapShellToMix

    integer, parameter :: mhd2mix_varn = 3

    contains

    subroutine init_shell2mix(s2m, shGr, mixApp)
        type(shell2mix_T), intent(inout) :: s2m
        type(ShellGrid_T), intent(in)    :: shGr
        type(mixApp_T),    intent(in)    :: mixApp

        type(mixGrid_T) :: virtGS

        !NORTH grid already in voltron frame; SOUTH needs a mirrored virtual grid.
        call shell_to_mix_set_map(shGr, mixApp%ion(NORTH)%G, s2m%maps(NORTH))

        call buildSouthVirtualGrid(mixApp%ion(SOUTH)%G, virtGS)
        call shell_to_mix_set_map(shGr, virtGS, s2m%maps(SOUTH))
        call releaseMixGrid(virtGS)
    end subroutine init_shell2mix


    !Apply maps + inner-boundary floor; call after App%fac/density/soundSpd staging.
    subroutine mapShellToMix(s2m, App, g2s, mixApp)
        type(shell2mix_T),     intent(in)    :: s2m
        class(remixCoupler_T), intent(in)    :: App
        type(gam2shell_T),     intent(in)    :: g2s
        type(mixApp_T),        intent(inout) :: mixApp

        real(rp), dimension(:,:), allocatable :: F, src
        logical,  dimension(:,:), allocatable :: mask
        integer :: h, Nth, v, varID, mhdID

        do h = NORTH, SOUTH
            Nth = mixApp%ion(h)%G%Nt

            !Per-hemi mask: cells equatorward of MHD inner boundary (matches flip_grid).
            allocate(mask(mixApp%ion(h)%G%Np, mixApp%ion(h)%G%Nt))
            mask = (sqrt(g2s%mhd_Rin) * sin(mixApp%ion(h)%G%t) > 1._rp - TINY)

            do v = 1, mhd2mix_varn
                select case (v)
                case (1); mhdID = MHDJ; varID = FAC
                case (2); mhdID = MHDD; varID = DENSITY
                case (3); mhdID = MHDC; varID = SOUND_SPEED
                end select

                !ShellGridVar is (theta, phi); map_grids wants (phi, theta).
                if (allocated(src)) deallocate(src)
                select case (mhdID)
                case (MHDJ)
                    allocate(src(size(App%fac%data, 2), size(App%fac%data, 1)))
                    src = transpose(App%fac%data)
                case (MHDD)
                    allocate(src(size(App%density%data, 2), size(App%density%data, 1)))
                    src = transpose(App%density%data)
                case (MHDC)
                    allocate(src(size(App%soundSpd%data, 2), size(App%soundSpd%data, 1)))
                    src = transpose(App%soundSpd%data)
                end select

                call mix_map_grids(s2m%maps(h), src, F)

                !Boundary floor: FAC->0, D/Cs->floorVal (matches mapGameraToRemix).
                select case (mhdID)
                case (MHDJ)
                    where (mask) F = 0._rp
                case (MHDD, MHDC)
                    where (mask) F = g2s%floorVal(mhdID, h)
                end select

                mixApp%ion(h)%St%Vars(:, :, varID) = F
                if (mhdID == MHDJ) then
                    mixApp%ion(h)%St%Vars(:, Nth - 1:Nth, FAC) = 0._rp
                endif
            enddo

            deallocate(mask)
        enddo

        if (allocated(src)) deallocate(src)
        if (allocated(F))   deallocate(F)
    end subroutine mapShellToMix


    !SOUTH virtual grid: same storage dims as ion(SOUTH)%G, (t,p) reflected into voltron frame.
    subroutine buildSouthVirtualGrid(Gs, virtG)
        type(mixGrid_T), intent(in)    :: Gs
        type(mixGrid_T), intent(inout) :: virtG

        virtG%Np = Gs%Np
        virtG%Nt = Gs%Nt
        virtG%isPeriodic = Gs%isPeriodic

        if (.not. allocated(virtG%t)) allocate(virtG%t(Gs%Np, Gs%Nt))
        if (.not. allocated(virtG%p)) allocate(virtG%p(Gs%Np, Gs%Nt))

        virtG%t = PI - Gs%t
        virtG%p = modulo(2._rp * PI - Gs%p, 2._rp * PI)
    end subroutine buildSouthVirtualGrid


    subroutine releaseMixGrid(G)
        type(mixGrid_T), intent(inout) :: G
        if (allocated(G%t)) deallocate(G%t)
        if (allocated(G%p)) deallocate(G%p)
    end subroutine releaseMixGrid

end module shell2mix_interface
