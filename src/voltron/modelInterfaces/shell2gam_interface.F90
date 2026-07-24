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
    use planethelper, only: InvLatitude
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
        real(rp),             save :: thAvgFull ! colat where averaging reaches full weight
        real(rp),             save :: rampW     ! transition width (2 local grid cells)
        real(rp) :: xyz(NDIM), x2, x1c, x1N, x1S, gI, gN, gS, invlat, wAvg
        real(rp) :: thFoot, s
        integer  :: i, j, k, iG, nbc, iC

        associate(Gr=>vApp%gApp%Grid, sh=>vApp%shGrid)

        ! Electrostatic-only potential on the shell (corotation removed); corner-centered.
        ! dPot is cached (initShellVar once); its data is fully overwritten each call.
        if (.not. dPotInit) then
            call initShellVar(sh, SHGR_CORNER, dPot)
            invlat = InvLatitude(Gr%xyz(PsiSt+PsiSh, Gr%js, Gr%ks, :))
            thFoot = PI/2._rp - invlat    ! outermost shell
            iC = lbound(sh%th, 1)
            do while (iC < ubound(sh%th, 1) - 1 .and. sh%th(iC+1) < thFoot)
                iC = iC + 1
            enddo
            rampW     = 2.0_rp*(sh%th(iC+1) - sh%th(iC))
            thAvgFull = thFoot - rampW    ! TSC-stencil margin
            dPotInit = .true.
        endif
        dPot%data = vApp%State%potential_total%data - vApp%State%potential_corot%data
        dPot%mask = .true.

        vApp%gApp%gPsi = 0.0_rp

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,iG,j,k,xyz,x2,x1c,x1N,x1S,gI,gN,gS,invlat,wAvg,s)
        do k = Gr%ks, Gr%ke+1
            do j = Gr%js, Gr%je+1
                do i = 1, PsiSh+1
                    iG  = PsiSt + i - 1
                    xyz = Gr%xyz(iG, j, k, :)     ! cell corner (Rx; scale-free for lat/atan2)
                    x2  = katan2(xyz(YDIR), xyz(XDIR))
                    invlat = InvLatitude(xyz)
                    x1N = PI/2._rp - invlat       
                    x1S = PI/2._rp + invlat       
                    if (xyz(ZDIR) < 0.0_rp) then
                        x1c = x1S                 
                    else
                        x1c = x1N
                    endif
                    ! 0 = pure sign(z) sampling, 1 = full average.
                    s    = min(1.0_rp, max(0.0_rp, (x1N - (thAvgFull - rampW))/rampW))
                    wAvg = s*s*(3.0_rp - 2.0_rp*s)
                    if (wAvg <= 0.0_rp) then
                        call InterpShellVar_TSC_pnt(sh, dPot, x1c, x2, gI)
                        vApp%gApp%gPsi(i, j, k) = gI + pointCorotPotential(x1c, x2)
                    else
                        call InterpShellVar_TSC_pnt(sh, dPot, x1N, x2, gN)
                        call InterpShellVar_TSC_pnt(sh, dPot, x1S, x2, gS)
                        gN = gN + pointCorotPotential(x1N, x2)
                        gS = gS + pointCorotPotential(x1S, x2)
                        if (xyz(ZDIR) < 0.0_rp) then
                            gI = gS             
                        else
                            gI = gN
                        endif
                        vApp%gApp%gPsi(i, j, k) = wAvg*0.5_rp*(gN + gS) &
                                                + (1.0_rp - wAvg)*gI
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
