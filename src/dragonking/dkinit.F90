module dkinit
    use kdefs
    use dkdefs
    use dktypes
    use dkconductance
    use dkaurora
    use auroralhelper
    use math
    use kronos
    use xml_input
    use raijudefs, ONLY : def_tiote
    use earthhelper, only: GallagherRP
    use shellGrid, only: initShellVar, SHGR_CC
    
    implicit none

    contains

   !  subroutine initDKmesh(dkApp, I)
   !    type(dkApp_T), intent(inout) :: dkApp
   !    type(mixIon_T),dimension(:),intent(in) :: I  
   !    integer :: h

   !    if (.not. allocated(dkApp%grid)) allocate(dkApp%grid(size(I)))
      
   !    do h=1,size(I)
   !      dkApp%grid(h) = I(h)%G
   !    end do

   !    dkApp%params%R = I(1)%P%R
   !    dkApp%params%RinMHD = I(1)%P%RinMHD

   !  end subroutine initDKmesh

    subroutine initDK(dkApp, inpXML)
      type(dkApp_T), intent(inout) :: dkApp
      type(XML_Input_T), intent(in) :: inpXML
      
      call initDKparams(dkApp%params, inpXML)
      call initDKvars(dkApp%aurora, dkApp%conductance, dkApp%vars, dkApp%grid)

      ! Cold-start for the conductance EMA: begin from the const_sigma reference
      dkApp%vars%SigmaP%data = dkApp%params%ped0
      dkApp%vars%SigmaH%data = 0.0_rp

      dkApp%params%R = dkApp%grid%radius
      dkApp%aurora%alpha_RCM = 1.0/(def_tiote+1.0)
      dkApp%aurora%beta_RCM = dkApp%params%beta
      dkApp%aurora%gtype_RCM = 0.0
         
    end subroutine initDK

   subroutine initDKgrid(dkApp, shGrid)
      type(dkApp_T), intent(inout) :: dkApp
      type(ShellGrid_T), intent(in) :: shGrid
      integer :: Np, Nt, j, is, ie, js, je
      real(rp), dimension(:), allocatable :: tmp

      ! DK runs on the active grid only: ghost cells are dropped
      dkApp%grid%ShellGrid_T = shGrid

      is = shGrid%is; ie = shGrid%ie
      js = shGrid%js; je = shGrid%je
      Nt = ie - is + 1   ! active theta cells
      Np = je - js + 1   ! active phi cells

      ! DK grid ghost-free
      tmp = dkApp%grid%phc(js:je);   call move_alloc(tmp, dkApp%grid%phc)
      tmp = dkApp%grid%ph (js:je+1); call move_alloc(tmp, dkApp%grid%ph)
      tmp = dkApp%grid%thc(is:ie);   call move_alloc(tmp, dkApp%grid%thc)
      tmp = dkApp%grid%th (is:ie+1); call move_alloc(tmp, dkApp%grid%th)

      dkApp%grid%Nt = Nt;   dkApp%grid%Np  = Np
      dkApp%grid%is = 1;    dkApp%grid%ie  = Nt
      dkApp%grid%js = 1;    dkApp%grid%je  = Np
      dkApp%grid%isg = 1;   dkApp%grid%ieg = Nt
      dkApp%grid%jsg = 1;   dkApp%grid%jeg = Np
      dkApp%grid%Ngn = 0;   dkApp%grid%Ngs = 0
      dkApp%grid%Ngw = 0;   dkApp%grid%Nge = 0

      ! DK arrays are shell-native (theta,phi) = (Nt,Np), matching ShellGridVar%data.
      allocate(dkApp%grid%x (Nt,Np), dkApp%grid%y (Nt,Np), dkApp%grid%r (Nt,Np), &
             dkApp%grid%D0(Nt,Np), dkApp%grid%dt(Nt,Np), dkApp%grid%dp(Nt,Np))

      ! SM coordinates: x = sin(th)*cos(ph), y = sin(th)*sin(ph)
      dkApp%grid%x  = spread(sin(dkApp%grid%thc), 2, Np) * spread(cos(dkApp%grid%phc), 1, Nt)
      dkApp%grid%y  = spread(sin(dkApp%grid%thc), 2, Np) * spread(sin(dkApp%grid%phc), 1, Nt)

      ! FEDDER ramp coordinate: sin(th), symmetric about equator, phi-independent
      dkApp%grid%r  = spread(sin(dkApp%grid%thc), 2, Np)

      ! Gallagher background density: function of radius and colatitude, phi-independent
      do j=1,Nt
         dkApp%grid%D0(j,:) = GallagherRP(dkApp%grid%radius, dkApp%grid%thc(j))
      end do

      ! Cell spacings for CalcMagE: differences of corner arrays
      dkApp%grid%dt = spread(dkApp%grid%th(2:Nt+1) - dkApp%grid%th(1:Nt), 2, Np)
      dkApp%grid%dp = spread(dkApp%grid%ph(2:Np+1) - dkApp%grid%ph(1:Np), 1, Nt)

    end subroutine initDKgrid


    subroutine initDKparams(params, inpXML)
      type(dkParameters_T), intent(inout) :: params
      type(XML_Input_T), intent(in) :: inpXML
      character(len=strLen) :: tmpStr
      logical :: doSuicide
      type(XML_Input_T) :: xmlInp
      
      doSuicide = .false.

      ! Make sure root is Kaiju/Remix
      call inpXML%GetFileStr(tmpStr)
      ! Create new XML reader w/ REMIX as root
      xmlInp = New_XML_Input(trim(tmpStr),'Kaiju/REMIX',.true.)

      ! =========== CONDUCTANCE MODEL PARAMTERS =================== !
      ! EUV_MODEL_TYPE
      !NOTE: Changing default EUV to LOMPE, it's 2021 ffs
      call xmlInp%Set_Val(tmpStr,"conductance/euv_model_type","LOMPE")
      select case (tmpSTR)
         case ("LOMPE")
            params%euv_model_type = LOMPE
         case default 
            stop "The EUV model type entered is not supported."              
      end select

      ! ET_MODEL_TYPE
      call xmlInp%Set_Val(tmpStr,"conductance/et_model_type","NONE")
      select case (tmpSTR)
         case ("NONE")
            params%et_model_type = NONE
         case ("ADHOC")
            params%et_model_type = ADHOC
         case ("CMIT")
            params%et_model_type = CMIT
         case default 
            stop "The ET model type entered is not supported."              
      end select

      ! SIGMA_MODEL_TYPE
      call xmlInp%Set_Val(tmpStr,"conductance/sigma_model_type","KAE23")
      select case (tmpSTR)
         case ("KAE23")
            params%sigma_model_type = KAE23
         case ("ROB87")
            params%sigma_model_type = ROB87
         case default
            stop "The SIGMA model type entered is not supported."
      end select

      ! AURORA_MODEL_TYPE
      call xmlInp%Set_Val(tmpStr,"precipitation/aurora_model_type","FEDDER")
      select case (tmpSTR)
         case ("FEDDER")
            params%aurora_model_type = FEDDER
         case ("SUNNY")
            params%aurora_model_type = SUNNY
         case ("LINMRG")
            params%aurora_model_type = LINMRG
         case default 
            stop "The auroral model type entered is not supported (Available options: FEDDER, SUNNY, LINMRG)."
      end select

      ! =========== CONDUCTANCE MODEL PARAMETERS =================== !
      call xmlInp%Set_Val(params%alpha,"precipitation/alpha",1.0332467_rp)
      call xmlInp%Set_Val(params%beta,"precipitation/beta",0.4362323_rp)

      call xmlInp%Set_Val(params%doAuroralSmooth,"precipitation/doAuroralSmooth",.false.)        
      call xmlInp%Set_Val(params%F107,"conductance/F107",120._rp)
      call xmlInp%Set_Val(params%pedmin,"conductance/pedmin",2.0_rp)
      call xmlInp%Set_Val(params%hallmin,"conductance/hallmin",1.0_rp)
      call xmlInp%Set_Val(params%sigma_ratio,"conductance/sigma_ratio",6.0_rp)
      call xmlInp%Set_Val(params%ped0,"conductance/ped0",10.0_rp)
      call xmlInp%Set_Val(params%const_sigma,"conductance/const_sigma",.false.)
      call xmlInp%Set_Val(params%doRamp,"conductance/doRamp",.true.)
      call xmlInp%Set_Val(params%doChill,"conductance/doChill",.false.)
      call xmlInp%Set_Val(params%doStarlight,"conductance/doStarlight",.false.)        
      call xmlInp%Set_Val(params%doMR,"conductance/doMR",.false.)        
      call xmlInp%Set_Val(params%apply_cap,"conductance/apply_cap",.true.)
      call xmlInp%Set_Val(params%gamma,"precipitation/gamma",5.0/3.0)
      call xmlInp%Set_Val(params%doEMA,"conductance/doEMA",.true.)
      call xmlInp%Set_Val(params%doET,"conductance/doET",.false.)

      params%beta_inp = params%beta

      ! if this is on, params%F107 read from xml file above will be overwritten
      ! and F107 will be set from the SW file (i.e., from OMNI data)
      ! call xmlInp%Set_Val(params%doSWF107,"conductance/doSWF107",.true.)      !not sure this applies to DK
      ! =========== CONDUCTANCE MODEL PARAMETERS =================== !


      ! Check for old-style precipitation quantities
      write(*,*) ANSIRED
      if (xmlInp%Exists("conductance/alpha")) then
         write(*,*) "conductance/alpha is gone. You're entering a world of pain. Use precipitation/alpha"
         doSuicide=.true.
      else if (xmlInp%Exists("conductance/beta")) then
         write(*,*) "conductance/beta is gone. You're entering a world of pain. Use precipitation/beta"
         doSuicide=.true.
      else if (xmlInp%Exists("conductance/R")) then
         write(*,*) "conductance/R is gone. You're entering a world of pain. Use precipitation/R"
         doSuicide=.true.
      else if (xmlInp%Exists("conductance/aurora_model_type")) then
         write(*,*) "conductance/aurora_model_type is gone. You're entering a world of pain. Use precipitation/aurora_model_type"
         doSuicide=.true.
      end if
      write(*,*) ANSIRESET, ''
      if (doSuicide) stop

    end subroutine initDKparams

    subroutine initDKvars(aurora,conductance,vars,grid)
      type(dkAurora_T) , intent(inout) :: aurora
      type(dkConductance_T) , intent(inout) :: conductance
      type(dkVariables_T) , intent(inout) :: vars
      type(dkGrid_T) , intent(in) :: grid
      integer :: Np, Nt

      Np = grid%Np
      Nt = grid%Nt

      ! Aurora/conductance scratch: plain (theta,phi)=(Nt,Np) arrays (internal, never coupled out)
      if (.not. allocated(aurora%E0))         allocate(aurora%E0(Nt,Np))
      if (.not. allocated(aurora%phi0))       allocate(aurora%phi0(Nt,Np))
      if (.not. allocated(aurora%deltaE))     allocate(aurora%deltaE(Nt,Np))
      if (.not. allocated(aurora%rampFactor)) allocate(aurora%rampFactor(Nt,Np))
      if (.not. allocated(aurora%PrecipMask)) allocate(aurora%PrecipMask(Nt,Np))

      if (.not. allocated(aurora%beta_RCM))  allocate(aurora%beta_RCM(Nt,Np))
      if (.not. allocated(aurora%alpha_RCM)) allocate(aurora%alpha_RCM(Nt,Np))
      if (.not. allocated(aurora%gtype_RCM)) allocate(aurora%gtype_RCM(Nt,Np))

      if (.not. allocated(conductance%zenith))      allocate(conductance%zenith(Nt,Np))
      if (.not. allocated(conductance%coszen))      allocate(conductance%coszen(Nt,Np))
      if (.not. allocated(conductance%euvSigmaP))   allocate(conductance%euvSigmaP(Nt,Np))
      if (.not. allocated(conductance%euvSigmaH))   allocate(conductance%euvSigmaH(Nt,Np))
      if (.not. allocated(conductance%deltaSigmaP)) allocate(conductance%deltaSigmaP(Nt,Np))
      if (.not. allocated(conductance%deltaSigmaH)) allocate(conductance%deltaSigmaH(Nt,Np))

      ! DK variables as ShellGridVar_T (cell-centered) on the ghost-free DK grid.
      ! initShellVar allocates %data(1:Nt,1:Np) and zero-fills; %mask defaults .true.
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Pot)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%FAC)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%SigmaP)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%SigmaH)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Sound_Speed)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Density)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Avg_Eng)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Num_Flux)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Im_Eflux)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Z_Eavg)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Z_Nflux)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Im_Gtype)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Aur_Type)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Im_Beta)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Im_Eden)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Im_Epre)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Im_Enflx)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%DeltaE)
      call initShellVar(grid%ShellGrid_T, SHGR_CC, vars%Im_Npsp)

      vars%Aur_Type%data = AT_NoPre

      vars%SigmaP%data = 0.0_rp
      vars%SigmaH%data = 0.0_rp

    end subroutine initDKvars



end module dkinit