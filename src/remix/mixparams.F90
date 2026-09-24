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

      ! NOTE: precipitation/* and conductance/* are read by DragonKing from its own
      !       <DRAGONKING> section (initDKparams).

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
