module dkmain
    use kdefs
    use gcmtypes
    use dkdefs
    use dktypes
    use dkconductance
    use dkaurora
    use auroralhelper
    use math
    use kronos
    
    implicit none

    contains

!---------- Compute precipitation and conductance/conductivity
    subroutine DKcompute(dkApp)
        type(dkApp_T),intent(inout) :: dkApp
        logical :: doModel
        
        doModel = .true. ! always default to xml input deck unless doModelOpt is present and on

        if (dkApp%params%doGCM) then
            doModel = .true.
        else
            if (dkApp%IO%time<=0) doModel = .false.
        end if

        if (doModel) then
            dkApp%conductance%const_sigma = dkApp%params%const_sigma
        else
            dkApp%conductance%const_sigma = .true.
        end if

        call Tic("DK-COND")
        call conductance_euv(dkApp%conductance,dkApp%grid,dkApp%vars,dkApp%params)
        call precipitation_total(dkApp%aurora,dkApp%grid,dkApp%vars,dkApp%conductance,dkApp%params) !to adapt?
        call conductance_total(dkApp%conductance,dkApp%grid,dkApp%vars,dkApp%params)
        call Toc("DK-COND")
    end subroutine DKcompute

    ! subroutine DKdogcm(dkApp, gcm)
    !     type(dkApp_T), intent(inout) :: dkApp
    !     type(gcm_T), intent(in) :: gcm
    !     integer :: h

    !     do h = 1, size(dkApp%vars)
    !         call apply_gcm2dk(gcm, dkApp%vars(h), h)   ! writes into vars%SigmaP/SigmaH
    !     end do
    !     dkApp%params%doGCM = .true.

    ! end subroutine DKdogcm

    ! subroutine apply_gcm2dk(gcm,vars,h)
    !   type(gcm_T),intent(in) :: gcm
    !   type(dkVariables_T),intent(inout) :: vars
    !   integer,intent(in) :: h

    !   vars%SigmaPgcm = gcm%APEX%mixInput(h,1)%var !sigmaP
    !   vars%SigmaHgcm = gcm%APEX%mixInput(h,2)%var !sigmaH
    ! end subroutine apply_gcm2dk

end module dkmain