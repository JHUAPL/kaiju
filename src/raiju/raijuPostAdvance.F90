module raijuPostAdvancer

    use kdefs
    use math
    use planethelper
    use raijutypes
    use raijuSpeciesHelper
    use raijuetautils

    use, intrinsic :: ieee_arithmetic

    implicit none

    contains


    subroutine raijuPostAdvanceChecks(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State

        integer :: s
        ! NaN checks
        if (any(ieee_is_nan(State%eta_avg       ))) call handleRaiPostCheckFailed("NaNs in RAIJU's eta_avg. Writing black box and exiting")
        do s=1,Model%nSpc
            if (any(ieee_is_nan(State%Press_avg  (s)%data))) call handleRaiPostCheckFailed("NaNs in RAIJU's Press_avg. Writing black box and exiting")
            if (any(ieee_is_nan(State%Den_avg    (s)%data))) call handleRaiPostCheckFailed("NaNs in RAIJU's Den_avg. Writing black box and exiting")
            if (any(ieee_is_nan(State%precipEFlux(s)%data))) call handleRaiPostCheckFailed("NaNs in RAIJU's precipEFlux. Writing black box and exiting")
            if (any(ieee_is_nan(State%precipNFlux(s)%data))) call handleRaiPostCheckFailed("NaNs in RAIJU's precipNFlux. Writing black box and exiting")
        enddo

        contains

        subroutine handleRaiPostCheckFailed(msg)
            character(len=*), intent(in) :: msg

            write(*,*)msg
            call State%myRaiApp%WriteBlackBox()
            stop
        end subroutine handleRaiPostCheckFailed
    
    end subroutine raijuPostAdvanceChecks

    subroutine initTepidBnd(Model, Grid, eta_tepid, LFullO, LEmptyO)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        real(rp), dimension(:,:,:), allocatable, intent(out) :: eta_tepid
        real(rp), intent(in), optional :: LFullO, LEmptyO

        integer :: i,j,iSpc
        real(rp) :: P = 0.02  ! [nPa]
        real(rp) :: kT = 10.0  ! [keV]
        real(rp), dimension(2) :: L_ramp
            !! Bounds over which to ramp from 1 to zero
        real(rp) :: N, L, Pscl, vm
        write(*,*)"InitTepidBnd"

        associate(sh=>Grid%shGrid, spc_ele=>Grid%spc(spcIdx(Grid, F_HOTE)), spc_p=>Grid%spc(spcIdx(Grid, F_HOTP)))

        if(present(LFullO)) then
            L_ramp(1) = LFullO
        else
            L_ramp(1) = DipColat2L(sh%thc(sh%ie)) + 1.0_rp
        endif
        if(present(LEmptyO)) then
            L_ramp(2) = LEmptyO
        else
            L_ramp(2) = L_ramp(1) + 1.0_rp
        endif

        if (allocated(eta_tepid)) deallocate(eta_tepid)

        allocate(eta_tepid(sh%isg:sh%ieg,sh%jsg:sh%jeg,1:Grid%Nk))
        eta_tepid = 0.0_rp

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,iSpc,L,vm,Pscl,N)
        do i=sh%isg,sh%ieg
            L = DipColat2L(sh%thc(i))
            vm = DipFTV_colat(sh%thc(i))**(-2./3.)
            Pscl = RampDown(L, L_ramp(1), L_ramp(2))
            N = PkT2Den(Pscl, kT)
            
            do j=sh%jsg,sh%jeg
                do iSpc=1,Grid%nSpc
                    associate(spc=>Grid%spc(iSpc))
                
                    if (spc%flav == F_PSPH) cycle  ! Don't touch plasmasphere
                    call DkT2SpcEta(Model, spc, eta_tepid(i,j,spc%kStart:spc%kEnd),&
                                    N, kT, vm)
                
                    end associate
                enddo
            enddo
        enddo

        end associate

    end subroutine initTepidBnd


    subroutine applyTepidBnd_k(Model, sh, eta_tepid, eta)
        !! Apply tepid boundary for a single energy channel
        type(raijuModel_T), intent(in) :: Model
        type(ShellGrid_T) , intent(in) :: sh
        real(rp), dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(in) :: eta_tepid
            !! Eta values for tepid floor
        real(rp), dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(inout) :: eta
            !! Eta values we are potentially adjusting

        integer :: i,j

        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg
                if (eta(i,j) < eta_tepid(i,j)) then
                    eta(i,j) = eta_tepid(i,j)
                endif
            enddo
        enddo
    end subroutine applyTepidBnd_k

end module raijuPostAdvancer