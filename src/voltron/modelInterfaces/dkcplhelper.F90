! Voltron -> DragonKing coupling: adapts shell state (tubeShell feed + imag2shell
! precip) into dkApp%vars and advances DK. The DK -> REMIX hop is remixCplTypesSub's mapDKToMix.
module dkcplhelper

    use volttypes
    use kdefs, only: Mp_cgs          ! #/cc -> mass density for the tubeShell feed
    use gdefs, only: BLK             ! bulk fluid index for tubeShell bcN/bcP
    use planethelper, only: DP2Cs    ! (D,P) -> sound speed
    use shellGrid, only: ShellGrid_T, SHGR_CC, SHGR_CORNER, initShellVar

    implicit none

    contains

    ! DragonKing conductance on the full-sphere shell from the tubeShell feed (J11i/bcN/bcP),
    ! corner->cell averaged. Density/Sound_Speed skip bcN=0 no-sample corners; FAC 0 is physical.
    ! Runs before toIono (mapDKToMix reads dkApp%vars).
    subroutine shellDK(vApp)
        class(voltApp_T), intent(inout) :: vApp
        real(rp) :: curTilt
        real(rp) :: dn(4), cs(4), bn(4)
        real(rp) :: dmin, cmin
        integer  :: i, j, nGood

        ! DK IO/tilt update (DKcompute forces const_sigma when IO%time<=0, so spinup is safe).
        call vApp%tilt%getValue(vApp%time, curTilt)
        vApp%dkApp%IO%time   = vApp%time
        vApp%dkApp%IO%mjd    = vApp%MJD
        vApp%dkApp%IO%tilt   = curTilt
        vApp%dkApp%vars%tilt = vApp%dkApp%IO%tilt

        associate(sg  => vApp%shGrid, &
                  j11 => vApp%State%tubeShell%J11i%data, &
                  bcN => vApp%State%tubeShell%bcN(BLK)%data, &
                  bcP => vApp%State%tubeShell%bcP(BLK)%data)
            do i = sg%is, sg%ie          ! theta (active)
                do j = sg%js, sg%je      ! phi (active)
                    ! FAC: plain 4-corner average [uA/m^2]
                    vApp%dkApp%vars%FAC(j-sg%js+1, i-sg%is+1) = &
                        0.25_rp*( j11(i,j)+j11(i+1,j)+j11(i,j+1)+j11(i+1,j+1) )

                    ! Density/Sound_Speed: mask-aware 4-corner average (bcN>0 closed corners only).
                    bn = [ bcN(i,j), bcN(i+1,j), bcN(i,j+1), bcN(i+1,j+1) ]
                    dn = bn * Mp_cgs                                            ! g/cm^3
                    cs = [ DP2Cs(bcN(i,j),    bcP(i,j)),    DP2Cs(bcN(i+1,j),  bcP(i+1,j)),   &
                           DP2Cs(bcN(i,j+1),  bcP(i,j+1)),  DP2Cs(bcN(i+1,j+1),bcP(i+1,j+1)) ] * 1.0e5_rp
                    nGood = count(bn > 0.0_rp)
                    if (nGood > 0) then
                        vApp%dkApp%vars%Density(j-sg%js+1, i-sg%is+1)     = sum(dn, mask=bn>0.0_rp) / nGood
                        vApp%dkApp%vars%Sound_Speed(j-sg%js+1, i-sg%is+1) = sum(cs, mask=bn>0.0_rp) / nGood
                    else
                        vApp%dkApp%vars%Density(j-sg%js+1, i-sg%is+1)     = 0.0_rp
                        vApp%dkApp%vars%Sound_Speed(j-sg%js+1, i-sg%is+1) = 0.0_rp
                    endif
                enddo
            enddo
        end associate

        !Backfill no-sample cells so FEDDER never sees zero Density/Sound_Speed (dE = 0/0 -> NaN)
        if (any(vApp%dkApp%vars%Density > 0.0_rp)) then
            dmin = minval(vApp%dkApp%vars%Density, mask = vApp%dkApp%vars%Density > 0.0_rp)
        else
            dmin = 0.01_rp*Mp_cgs   ! feed entirely empty: 0.01 /cc
        endif
        if (any(vApp%dkApp%vars%Sound_Speed > 0.0_rp)) then
            cmin = minval(vApp%dkApp%vars%Sound_Speed, mask = vApp%dkApp%vars%Sound_Speed > 0.0_rp)
        else
            cmin = 1.0e6_rp          ! feed entirely empty: 10 km/s [cm/s]
        endif
        where (vApp%dkApp%vars%Density     <= 0.0_rp) vApp%dkApp%vars%Density     = dmin
        where (vApp%dkApp%vars%Sound_Speed <= 0.0_rp) vApp%dkApp%vars%Sound_Speed = cmin

        ! IM precip: imag2shell is SHGR_CC (theta,phi); dkVars is (phi,theta) -> transpose
        associate(sg => vApp%shGrid)
            vApp%dkApp%vars%Im_Eflux = transpose(vApp%imag2shell%Im_Eflux%data(sg%is:sg%ie, sg%js:sg%je))
            vApp%dkApp%vars%Im_Enflx = transpose(vApp%imag2shell%Im_Enflx%data(sg%is:sg%ie, sg%js:sg%je))
            vApp%dkApp%vars%Im_Gtype = transpose(vApp%imag2shell%Im_Gtype%data(sg%is:sg%ie, sg%js:sg%je))
            vApp%dkApp%vars%Im_Eden  = transpose(vApp%imag2shell%Im_Eden%data (sg%is:sg%ie, sg%js:sg%je))
            vApp%dkApp%vars%Im_Epre  = transpose(vApp%imag2shell%Im_Epre%data (sg%is:sg%ie, sg%js:sg%je))
            vApp%dkApp%vars%Im_Npsp  = transpose(vApp%imag2shell%Im_Npsp%data (sg%is:sg%ie, sg%js:sg%je))
        end associate

        call vApp%dkApp%AdvanceModel(0.0_rp)
    end subroutine shellDK

    !Average corner potential (pot_iono) to DK cell centers for doET
    subroutine shellPotToDK(vApp)
        class(voltApp_T), intent(inout) :: vApp
        integer :: i, j
        associate(sg => vApp%shGrid, pot => vApp%pot_iono%data)
            do j = sg%js, sg%je
                do i = sg%is, sg%ie
                    vApp%dkApp%vars%Pot(j-sg%js+1, i-sg%is+1) = &
                        0.25_rp * ( pot(i  ,j  ) + pot(i+1,j  ) + &
                                    pot(i  ,j+1) + pot(i+1,j+1) )
                enddo
            enddo
        end associate
    end subroutine shellPotToDK

    ! Allocate the IMAG->shell precip staging ShellGridVars (cell-centered, full-sphere)
    subroutine init_imag2shell(i2s, shGrid)
        type(imag2shell_T), intent(inout) :: i2s
        type(ShellGrid_T),  intent(in)    :: shGrid

        call initShellVar(shGrid, SHGR_CC, i2s%Im_Eflux)
        call initShellVar(shGrid, SHGR_CC, i2s%Im_Enflx)
        call initShellVar(shGrid, SHGR_CC, i2s%Im_Gtype)
        call initShellVar(shGrid, SHGR_CC, i2s%Im_Eden)
        call initShellVar(shGrid, SHGR_CC, i2s%Im_Epre)
        call initShellVar(shGrid, SHGR_CC, i2s%Im_Npsp)
    end subroutine init_imag2shell

end module dkcplhelper
