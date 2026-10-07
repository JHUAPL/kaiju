module mixparams
  use mixdefs
  use mixtypes
  use xml_input

  implicit none

  !Grid default parameter type
  type MixGrid0_T
    integer  :: Np=256,Nt=128
    real(rp) :: LowLatBC = 45.0
    real(rp) :: R0 = 0.083567956_rp
    real(rp) :: Rin = 0.0
  end type MixGrid0_T

  type(MixGrid0_T), private :: MixGrid0

  contains

    !Set default mix grid parameters based on inner MHD boundary and number of k cells
    subroutine SetMixGrid0(Rin,Nk)
      real(rp), intent(in) :: Rin
      integer , intent(in) :: Nk

      real(rp) :: dDeg

      MixGrid0%LowLatBC = asin(sqrt(1.0/Rin))*180.0/PI
      !Set default grid resolution (in degrees)
      if (Nk<=64 ) then
        !DBL
        dDeg = 1.0
        MixGrid0%R0 = 0.08
      else if (Nk<=128) then
        !QUAD
        dDeg = 1.0
        MixGrid0%R0 = 0.04
      else if (Nk<=256) then
        !OCT
        dDeg = 0.5
        MixGrid0%R0 = 0.02
      else
        !HEX or above
        dDeg = 0.25
        MixGrid0%R0 = 0.02 ! Maybe?
      endif

      MixGrid0%Np = nint(360.0/dDeg)
      MixGrid0%Nt = nint(MixGrid0%LowLatBC/dDeg)
      MixGrid0%Rin = Rin
    end subroutine SetMixGrid0

    subroutine initMIXParams(Params, inpXML)
      type(mixParams_T), intent(out) :: Params
      type(XML_Input_T), intent(in) :: inpXML
      ! character(len=*), optional, intent(in) :: optFilename
      character(len=strLen) :: tmpStr
      integer :: Narg
      logical :: fExist,doSuicide
      type(XML_Input_T) :: xmlInp
      
      doSuicide = .false.

      ! if(present(optFilename)) then
      !    ! read from the prescribed file
      !    inpXML = optFilename
      ! else
      !    !Find input deck
      !    Narg = command_argument_count()
      !    if (Narg .eq. 0) then
      !       write(*,*) 'No input deck specified, defaulting to mixparams.xml'
      !       inpXML = "mixparams.xml"
      !    else
      !       call get_command_argument(1,inpXML)
      !    endif
      ! endif

      ! Make sure root is Kaiju/Remix
      call inpXML%GetFileStr(tmpStr)
      ! Create new XML reader w/ REMIX as root
      xmlInp = New_XML_Input(trim(tmpStr),'Kaiju/REMIX',.true.)

      ! NOTE: precipitation (auroral) and EUV/ET/SIGMA conductance model
      !       selection now live in DragonKing only; the mix-side copies were
      !       unused and have been removed.

      ! =========== CONDUCTANCE MODEL PARAMETERS =================== !
      call xmlInp%Set_Val(Params%alpha,"precipitation/alpha",1.0332467_rp)
      call xmlInp%Set_Val(Params%beta,"precipitation/beta",0.4362323_rp)
      call xmlInp%Set_Val(Params%R,"precipitation/R",MixGrid0%R0)
      call xmlInp%Set_Val(Params%doAuroralSmooth,"precipitation/doAuroralSmooth",.false.)        
      call xmlInp%Set_Val(Params%F107,"conductance/F107",120._rp)
      call xmlInp%Set_Val(Params%pedmin,"conductance/pedmin",2.0_rp)
      call xmlInp%Set_Val(Params%hallmin,"conductance/hallmin",1.0_rp)
      call xmlInp%Set_Val(Params%sigma_ratio,"conductance/sigma_ratio",6.0_rp)
      call xmlInp%Set_Val(Params%ped0,"conductance/ped0",10.0_rp)
      call xmlInp%Set_Val(Params%const_sigma,"conductance/const_sigma",.false.)
      call xmlInp%Set_Val(Params%doRamp,"conductance/doRamp",.true.)
      call xmlInp%Set_Val(Params%doChill,"conductance/doChill",.false.)
      call xmlInp%Set_Val(Params%doStarlight,"conductance/doStarlight",.false.)        
      call xmlInp%Set_Val(Params%doMR,"conductance/doMR",.false.)        
      call xmlInp%Set_Val(Params%apply_cap,"conductance/apply_cap",.true.)
      call xmlInp%Set_Val(Params%gamma,"precipitation/gamma",5.0/3.0)
      call xmlInp%Set_Val(Params%doEMA,"conductance/doEMA",.true.)
      call xmlInp%Set_Val(Params%doET,"conductance/doET",.false.)

      ! if this is on, Params%F107 read from xml file above will be overwritten
      ! and F107 will be set from the SW file (i.e., from OMNI data)
      call xmlInp%Set_Val(Params%doSWF107,"conductance/doSWF107",.true.)
      ! =========== CONDUCTANCE MODEL PARAMETERS =================== !

      ! =========== SOLVER PARAMETERS =================== !
      call xmlInp%Set_Val(Params%maxitr,"solver/maxitr",400)
      call xmlInp%Set_Val(Params%mr,"solver/mr",30)
      call xmlInp%Set_Val(Params%tol_abs,"solver/tol_abs",0.000001_rp)
      call xmlInp%Set_Val(Params%tol_rel,"solver/tol_rel",0.000001_rp)
      call xmlInp%Set_Val(Params%llbc_value,"solver/llbc_value",0.0_rp)
      ! =========== SOLVER PARAMETERS =================== !

      ! =========== GRID PARAMETERS =================== !
      call xmlInp%Set_Val(Params%Np,"grid/Np",MixGrid0%Np)
      call xmlInp%Set_Val(Params%Nt,"grid/Nt",MixGrid0%Nt)
      call xmlInp%Set_Val(Params%LowLatBoundary,"grid/LowLatBoundary",MixGrid0%LowLatBC)
      Params%RinMHD = MixGrid0%Rin !May be zero if not connected to MHD
      ! =========== GRID PARAMETERS =================== !

      ! =========== INIT PARAMETERS =================== !
      call xmlInp%Set_Val(Params%init_from_file,"grid/init_from_file",.false.)
      ! =========== INIT PARAMETERS =================== !

      ! =========== IO PARAMETERS =================== !
      call xmlInp%Set_Val(Params%dtOut,"output/dtOut",1.0_rp)
      call xmlInp%Set_Val(Params%nRes,"/Kaiju/gamera/restart/nRes",-1)
      ! =========== IO PARAMETERS =================== !

      ! =========== DEBUG PARAMETERS ================ !
      call xmlInp%Set_Val(Params%mklmsglvl,"debug/mklmsglvl", 0)
      ! =========== DEBUG PARAMETERS   ============== !

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
      end if
      write(*,*) ANSIRESET, ''
      if (doSuicide) stop

    end subroutine initMIXParams
end module mixparams
