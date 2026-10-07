module dkmain
    use kdefs
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

        ! DK is GCM-agnostic: it always computes its own precipitation and conductance.
        ! When GCM coupling is active, GCM conductance overwrites DK's directly in REMIX
        ! (mapGCMToMix after mapDKToMix), and DK precip (Avg_Eng/Num_Flux) is exported to
        ! GCM -- so DK needs nothing from GCM and behaves identically here.
        if (dkApp%IO%time<=0) doModel = .false.

        if (doModel) then
            dkApp%conductance%const_sigma = dkApp%model%const_sigma
        else
            dkApp%conductance%const_sigma = .true.
        end if

        call Tic("DK-COND")
        call conductance_euv(dkApp%conductance,dkApp%grid,dkApp%vars,dkApp%model)
        call precipitation_total(dkApp%aurora,dkApp%grid,dkApp%vars,dkApp%conductance,dkApp%model) !to adapt?
        call conductance_total(dkApp%conductance,dkApp%grid,dkApp%vars,dkApp%model)
        call Toc("DK-COND")
    end subroutine DKcompute

end module dkmain
