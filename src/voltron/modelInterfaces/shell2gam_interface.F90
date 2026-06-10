! Map the ionospheric potential from the voltron ShellGrid to Gamera's inner boundary.
! TSC-interpolate the electrostatic part only (potential_total - potential_corot, which
! is smooth) and add corotation back analytically per corner, so TSC never softens its
! sharp sin^2(theta) gradient. The gPsi -> E-field write matches convertRemixToGamera.
module shell2gam_interface
    use kdefs
    use volttypes
    use cmiutils                   ! Ion2MHD  (+ PsiSh/PsiSt, Model_T/Grid_T transitively)
    use uservoltic                 ! IonInnerBC_T, FindBC, INI
    use shellGrid,    only: initShellVar, ShellGridVar_T, SHGR_CORNER
    use shellInterp,  only: InterpShellVar_TSC_pnt
    use planethelper, only: InvLatitude, sInvLatitude
    use math,         only: katan2
    use kai2geo,      only: geocorotation_from_SMTP

    implicit none
    private
    public :: convertVolttoGamera

    contains

    subroutine convertVolttoGamera(vApp)
        class(voltApp_T), intent(inout) :: vApp

        type(ShellGridVar_T), save :: dPot      ! cached across calls (shell grid is fixed)
        logical,              save :: dPotInit = .false.
        real(rp) :: xyz(NDIM), x2, x1c, x1N, x1S, gI, gN, gS, invlat
        integer  :: i, j, k, iG, nbc

        associate(Gr=>vApp%gApp%Grid, sh=>vApp%shGrid)

        ! Electrostatic-only potential on the shell (corotation removed); corner-centered.
        ! dPot is cached (initShellVar once); its data is fully overwritten each call.
        if (.not. dPotInit) then
            call initShellVar(sh, SHGR_CORNER, dPot)
            dPotInit = .true.
        endif
        dPot%data = vApp%State%potential_total%data - vApp%State%potential_corot%data
        dPot%mask = .true.

        vApp%gApp%gPsi = 0.0_rp

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,iG,j,k,xyz,x2,x1c,x1N,x1S,gI,gN,gS,invlat)
        do k = Gr%ks, Gr%ke+1
            do j = Gr%js, Gr%je+1
                do i = 1, PsiSh+1
                    iG  = PsiSt + i - 1
                    xyz = Gr%xyz(iG, j, k, :)     ! cell corner (Rx; scale-free for lat/atan2)
                    x2  = katan2(xyz(YDIR), xyz(XDIR))
                    if (abs(xyz(ZDIR)) < TINY) then
                        ! On the magnetic equator: average both conjugate footpoints.
                        invlat = InvLatitude(xyz)
                        x1N = PI/2._rp - invlat
                        x1S = PI/2._rp + invlat
                        call InterpShellVar_TSC_pnt(sh, dPot, x1N, x2, gN)
                        call InterpShellVar_TSC_pnt(sh, dPot, x1S, x2, gS)
                        gN = gN + pointCorotPotential(x1N, x2)
                        gS = gS + pointCorotPotential(x1S, x2)
                        vApp%gApp%gPsi(i, j, k) = 0.5_rp*(gN + gS)
                    else
                        x1c = PI/2._rp - sInvLatitude(xyz)
                        call InterpShellVar_TSC_pnt(sh, dPot, x1c, x2, gI)
                        vApp%gApp%gPsi(i, j, k) = gI + pointCorotPotential(x1c, x2)
                    endif
                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO

        ! Enforce the z=0 phi seam (k=ke+1 corner is geometrically k=ks).
        vApp%gApp%gPsi(:, :, Gr%ke+1) = vApp%gApp%gPsi(:, :, Gr%ks)

        ! Potential -> inner-boundary E-field (identical to convertRemixToGamera).
        nbc = FindBC(vApp%gApp%Model, Gr, INI)
        SELECT type(iiBC=>Gr%externalBCs(nbc)%p)
            TYPE IS (IonInnerBC_T)
                call Ion2MHD(vApp%gApp%Model, Gr, vApp%gApp%gPsi, &
                             iiBC%inEijk, iiBC%inExyz, vApp%gApp%rm2g)
            CLASS DEFAULT
                write(*,*) 'convertVolttoGamera: could not find Ion Inner BC'
                stop
        END SELECT

        end associate

    contains

        ! Analytic corotation potential at colatitude th, longitude ph [rad]; same
        ! formula as calcCorotPotential, so it re-adds exactly what was subtracted into dPot.
        function pointCorotPotential(th, ph) result(pc)
            real(rp), intent(in) :: th, ph
            real(rp) :: pc
            if (vApp%doGeoCorot) then
                call geocorotation_from_SMTP(th, ph, pc)
            else
                pc = -vApp%planet%psiCorot*(vApp%planet%rp_m/vApp%planet%ri_m)*sin(th)**2
            endif
        end function pointCorotPotential

    end subroutine convertVolttoGamera

end module shell2gam_interface
