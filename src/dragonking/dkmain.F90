module dkmain
    use kdefs
    use mixtypes
    use mixdefs
    ! use mixmain
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

!----------
    subroutine getMixToDK(mixApp,dkApp)
        class(mixApp_T), intent(in) :: mixApp
        class(dkApp_T), intent(inout) :: dkApp

        integer :: h,i,j

        dkApp%params%RinMHD = mixApp%ion(1)%P%RinMHD
        dkApp%params%gamma = mixApp%ion(1)%P%gamma

        do h=1,size(mixApp%ion)
            do j=1,mixApp%ion(h)%G%Nt
                do i=1,mixApp%ion(h)%G%Np
                    ! dkApp%vars(h)%Neutral_Wind(i,j) = mixApp%ion(h)%St%Vars(i,j,NEUTRAL_WIND)
                    ! dkApp%vars(h)%Efield(i,j) = mixApp%ion(h)%St%Vars(i,j,EFIELD)
                    ! dkApp%vars(h)%Im_Eavg(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_EAVG)
                    ! dkApp%vars(h)%Im_Iavg(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_IAVG)
                    ! dkApp%vars(h)%Im_Iflux(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_IFLUX)
                    ! dkApp%vars(h)%Z_Eavg(i,j) = mixApp%ion(h)%St%Vars(i,j,Z_EAVG)
                    ! dkApp%vars(h)%Z_Nflux(i,j) = mixApp%ion(h)%St%Vars(i,j,Z_NFLUX)
                    ! dkApp%vars(h)%CrPot(i,j) = mixApp%ion(h)%St%Vars(i,j,CR_POT)
                    ! dkApp%vars(h)%TPot(i,j) = mixApp%ion(h)%St%Vars(i,j,TPOT)
                    ! dkApp%vars(h)%Im_Inflx(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_INFLX)
                    ! dkApp%vars(h)%Im_Thcon(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_THCON)
                    ! dkApp%vars(h)%Im_Phcon(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_PHCON)

                    dkApp%vars(h)%Pot(i,j) = mixApp%ion(h)%St%Vars(i,j,POT)
                    dkApp%vars(h)%FAC(i,j) = mixApp%ion(h)%St%Vars(i,j,FAC)
                    dkApp%vars(h)%SigmaP(i,j) = mixApp%ion(h)%St%Vars(i,j,SIGMAP)
                    dkApp%vars(h)%SigmaH(i,j) = mixApp%ion(h)%St%Vars(i,j,SIGMAH)
                    dkApp%vars(h)%Sound_Speed(i,j) = mixApp%ion(h)%St%Vars(i,j,SOUND_SPEED)
                    dkApp%vars(h)%Density(i,j) = mixApp%ion(h)%St%Vars(i,j,DENSITY)
                    dkApp%vars(h)%Avg_Eng(i,j) = mixApp%ion(h)%St%Vars(i,j,AVG_ENG)
                    dkApp%vars(h)%Num_Flux(i,j) = mixApp%ion(h)%St%Vars(i,j,NUM_FLUX)
                    dkApp%vars(h)%Im_Eflux(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_EFLUX)
                    dkApp%vars(h)%Im_Gtype(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_GTYPE)
                    dkApp%vars(h)%Aur_Type(i,j) = mixApp%ion(h)%St%Vars(i,j,AUR_TYPE)
                    dkApp%vars(h)%Im_Beta(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_BETA)
                    dkApp%vars(h)%Im_Eden(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_EDEN)
                    dkApp%vars(h)%Im_Epre(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_EPRE)
                    dkApp%vars(h)%Im_Enflx(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_ENFLX)
                    dkApp%vars(h)%DeltaE(i,j) = mixApp%ion(h)%St%Vars(i,j,DELTAE)
                    dkApp%vars(h)%Im_Npsp(i,j) = mixApp%ion(h)%St%Vars(i,j,IM_NPSP)
                end do
            end do    
        end do

    end subroutine getMixToDK

!----------
    subroutine getDKToMix(mixApp,dkApp)
        class(mixApp_T), intent(inout) :: mixApp
        class(dkApp_T), intent(in) :: dkApp

        ! type(mixIon_T), dimension(:) :: ion
        integer :: h,i,j

        do h=1,size(mixApp%ion)
            do j=1,mixApp%ion(h)%G%Nt
                do i=1,mixApp%ion(h)%G%Np    
                    ! mixApp%ion(h)%St%Vars(i,j,NEUTRAL_WIND) = dkApp%vars(h)%Neutral_Wind(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,EFIELD) = dkApp%vars(h)%Efield(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,IM_EAVG) = dkApp%vars(h)%Im_Eavg(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,IM_IAVG) = dkApp%vars(h)%Im_Iavg(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,IM_IFLUX) = dkApp%vars(h)%Im_Iflux(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,Z_EAVG) = dkApp%vars(h)%Z_Eavg(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,Z_NFLUX) = dkApp%vars(h)%Z_Nflux(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,CR_POT) = dkApp%vars(h)%CrPot(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,TPOT) = dkApp%vars(h)%TPot(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,IM_INFLX) = dkApp%vars(h)%St%Vars(i,j,IM_Inflx)
                    ! mixApp%ion(h)%St%Vars(i,j,IM_THCON) = dkApp%vars(h)%Im_Thcon(i,j)
                    ! mixApp%ion(h)%St%Vars(i,j,IM_PHCON) = dkApp%vars(h)%Im_Phcon(i,j)

                    mixApp%ion(h)%St%Vars(i,j,POT) = dkApp%vars(h)%Pot(i,j)
                    mixApp%ion(h)%St%Vars(i,j,FAC) = dkApp%vars(h)%FAC(i,j)
                    mixApp%ion(h)%St%Vars(i,j,SIGMAP) = dkApp%vars(h)%SigmaP(i,j)
                    mixApp%ion(h)%St%Vars(i,j,SIGMAH) = dkApp%vars(h)%SigmaH(i,j)
                    mixApp%ion(h)%St%Vars(i,j,SOUND_SPEED) = dkApp%vars(h)%Sound_Speed(i,j)
                    mixApp%ion(h)%St%Vars(i,j,DENSITY) = dkApp%vars(h)%Density(i,j)
                    mixApp%ion(h)%St%Vars(i,j,AVG_ENG) = dkApp%vars(h)%Avg_Eng(i,j)
                    mixApp%ion(h)%St%Vars(i,j,NUM_FLUX) = dkApp%vars(h)%Num_Flux(i,j)
                    mixApp%ion(h)%St%Vars(i,j,IM_EFLUX) = dkApp%vars(h)%Im_Eflux(i,j)
                    mixApp%ion(h)%St%Vars(i,j,IM_GTYPE) = dkApp%vars(h)%Im_Gtype(i,j)
                    mixApp%ion(h)%St%Vars(i,j,AUR_TYPE) = dkApp%vars(h)%Aur_Type(i,j)
                    mixApp%ion(h)%St%Vars(i,j,IM_BETA) = dkApp%vars(h)%Im_Beta(i,j)
                    mixApp%ion(h)%St%Vars(i,j,IM_EDEN) = dkApp%vars(h)%Im_Eden(i,j)
                    mixApp%ion(h)%St%Vars(i,j,IM_EPRE) = dkApp%vars(h)%Im_Epre(i,j)
                    mixApp%ion(h)%St%Vars(i,j,IM_ENFLX) = dkApp%vars(h)%Im_Enflx(i,j)
                    mixApp%ion(h)%St%Vars(i,j,DELTAE) = dkApp%vars(h)%DeltaE(i,j)
                    mixApp%ion(h)%St%Vars(i,j,IM_NPSP) = dkApp%vars(h)%Im_Npsp(i,j)
                end do
            end do
        end do

    end subroutine getDKToMix


!---------- Compute precipitation and conductance/conductivity
    subroutine DKcompute(dkApp,gcm)
        type(dkApp_T),intent(inout) :: dkApp
        ! type(mixIon_T),dimension(:),intent(inout) :: I
        ! type(mixIO_T), intent(in) :: mixIOobj  
        type(gcm_T),optional,intent(inout) :: gcm

        ! logical :: doModel=.true.   ! always default to xml input deck unless doModelOpt is present and on
        logical :: doModel
        integer :: h,NumH

        ! here we should check if Remix or GPS
        NumH = size(dkApp%vars)
        doModel = .true. ! always default to xml input deck unless doModelOpt is present and on

        if (dkApp%params%doGCM) then
            doModel = .true.
        else
            if (dkApp%IO%time<=0) doModel = .false.
        end if

        do h=1,NumH
            if (dkApp%vars(h)%hemisphere.eq.NORTH) then
                dkApp%vars(h)%tilt = dkApp%IO%tilt
            else
                dkApp%vars(h)%tilt = -dkApp%IO%tilt
            end if

            if (doModel) then
                dkApp%conductance(h)%const_sigma = dkApp%params%const_sigma         
            else
                dkApp%conductance(h)%const_sigma = .true.            
            end if

            ! if (.not. doModel) dkApp%params%const_sigma = .true. 

            call Tic("DK-COND")
            call conductance_euv(dkApp%conductance(h),dkApp%grid(h),dkApp%vars(h),dkApp%params)
            call precipitation_total(dkApp%aurora(h),dkApp%grid(h),dkApp%vars(h),dkApp%conductance(h),dkApp%params) !to adapt?
            call conductance_total(dkApp%conductance(h),dkApp%grid(h),dkApp%vars(h),dkApp%params)
            call Toc("DK-COND")
        end do

    end subroutine DKcompute

    subroutine DKdogcm(dkApp, gcm)
        type(dkApp_T), intent(inout) :: dkApp
        type(gcm_T), intent(in) :: gcm
        integer :: h

        do h = 1, size(dkApp%vars)
            call apply_gcm2dk(gcm, dkApp%vars(h), h)   ! writes into vars%SigmaP/SigmaH
        end do
        dkApp%params%doGCM = .true.

    end subroutine DKdogcm

    subroutine apply_gcm2dk(gcm,vars,h)
      type(gcm_T),intent(in) :: gcm
      type(dkVariables_T),intent(inout) :: vars
      integer,intent(in) :: h

      vars%SigmaPgcm = gcm%APEX%mixInput(h,1)%var !sigmaP
      vars%SigmaHgcm = gcm%APEX%mixInput(h,2)%var !sigmaH
    end subroutine apply_gcm2dk

end module dkmain