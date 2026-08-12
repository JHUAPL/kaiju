module voltappHelper

    use kdefs
    use voltTypes
    use kai2geo
    use mixinterfaceutils
    use tubehelper
    use arrayutil
    
    implicit none

    contains

    subroutine initVoltState(vApp)
        type(voltApp_T), intent(inout) :: vApp

        associate(sh=>vApp%shGrid, State=>vApp%State)

            ! Allocations
            allocate(State%ijTubes(sh%is:sh%ie+1, sh%js:sh%je+1))
            call init_TubeShell(sh, State%tubeShell)
            call initShellVar(sh, SHGR_CORNER, State%potential_total)
            call initShellVar(sh, SHGR_CORNER, State%potential_corot)
            call initShellVar(sh, SHGR_CORNER, State%potential_base)
            call initShellVar(sh, SHGR_CORNER, State%pot_drop)

            call initShellVar(sh, SHGR_CC, State%bIonoMag)
            call initShellVar(sh, SHGR_CC, State%bIonoRad)

            ! DK precipitation outputs, mirrored from dkApp each step (shellDK)
            call initShellVar(sh, SHGR_CC, State%Num_Flux)
            call initShellVar(sh, SHGR_CC, State%Avg_Eng)
            call initShellVar(sh, SHGR_CC, State%Aur_Type)

            ! Default values
            call initMagVars(sh, vApp%planet, State)
        end associate 


        contains

        subroutine initMagVars(sh, planet, State)
            !! Default routine to init magnetic field-related variables on voltron grid. 
            !! TODO: Will make this more formal / someone else's job once we have different mag field models like IGRF
            type(ShellGrid_T), intent(in) :: sh
            type(planet_T   ), intent(in) :: planet
            type(voltState_t), intent(inout) :: State

            integer :: j
            real(rp), dimension(sh%isg:sh%ieg) :: cosThc
                !! cos(th) at cell centers
            real(rp), dimension(sh%isg:sh%ieg) :: BMagTh
                !! [nT] mag field along theta direction

            cosThc = cos(sh%thc)
            do j=sh%jsg,sh%jeg
                State%bIonoMag%data(:,j) = planet%magMoment*G2nT &
                                         /sh%radius**3.0 &
                                         * sqrt(1.0+3.0*cosThc**2.0)  ! [nT]
                State%bIonoRad%data(:,j) = planet%magMoment*G2nT &
                                         /sh%radius**3.0 &
                                         * 2*cosThc  ! [nT]
            enddo


        end subroutine initMagVars
        
    end subroutine initVoltState

    subroutine calcCorotPotential(planet, sh, pCorot,doGeoCorotO)
        ! Calculate corotation potential to store on voltron ShellGrid
        type(planet_T), intent(in) :: planet
        type(ShellGrid_T), intent(in) :: sh
        type(ShellGridVar_T), intent(inout) :: pCorot
            !! [kV] corotation potential we calculate
        logical, intent(in), optional :: doGeoCorotO
        
        integer :: i,j

        pCorot%data = 0.0
        pCorot%mask = .true.

        if (pCorot%loc .ne. SHGR_CORNER) then
            write(*,*)"Error in voltappHelper.F90's calcCorotPotential:"
            write(*,*)"Only handling corner-located corot var for now"
            stop
        endif

        if (present(doGeoCorotO)) then
            if (doGeoCorotO) then
                do j=sh%jsg,sh%jeg+1
                    do i=sh%isg,sh%ieg+1
                        ! geopack corotation returns in kV
                        call geocorotation_from_SMTP(sh%th(i), sh%ph(j), pCorot%data(i,j))
                    enddo
                enddo

                return

            endif
        endif

        ! Default: corotation for an aligned dipole + rotation axis
        do j=sh%jsg,sh%jeg+1
            pCorot%data(:,j) = -planet%psiCorot*(planet%rp_m/planet%ri_m)*sin(sh%th)**2  ! [kV]
        enddo

    end subroutine calcCorotPotential


    subroutine calcPotDrop(vApp)
        !! Calculate as a diagnostic quantity the potential drop along field lines by comparing potentail at either end of field line
        class(voltApp_T), intent(inout) :: vApp

        integer :: i,j
        real(rp) :: latc, colatc, lonc, potc
            !! Conjugate latitude, colatitude, longitude, and potential
        logical, dimension(vApp%shGrid%isg:vApp%shGrid%ieg+1,vApp%shGrid%jsg:vApp%shGrid%jeg+1) :: isGood

        associate(sh=>vApp%shGrid, State=>vApp%State)

        where (State%tubeShell%topo%data == TUBE_CLOSED)
            isGood = .true.
        elsewhere
            isGood = .false.
        end where

        State%pot_drop%data = 0.0_rp
        State%pot_drop%mask = .false.  ! We will turn on valid points in our loop

        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i, j, colatc, lonc, potc)
        do j=sh%jsg,sh%jeg+1
            do i=sh%isg,sh%ieg+1
                if (.not. isGood(i,j)) then
                    cycle
                endif

                ! Get conjugate points in the right coordinates
                colatc = PI/2.0_rp - State%tubeShell%latc%data(i,j)
                lonc = State%tubeShell%lonc%data(i,j)
                call InterpShellVar_TSC_pnt(sh, State%potential_total, colatc, lonc, potc, srcMaskO=isGood)
                State%pot_drop%data(i,j) = State%potential_total%data(i,j) - potc
                State%pot_drop%mask(i,j) = .true.

            enddo
        enddo

        end associate

    end subroutine

end module voltappHelper