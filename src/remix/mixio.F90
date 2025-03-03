  ! meant to replace mixio.F90 by using Kareems pipeline in base/ioH5.F90

module mixio
  use ioH5
  use mixdefs
  use mixtypes
  use files
  implicit none
  
  integer, parameter :: MAXMIXIOVAR = 100
  type(IOVAR_T), dimension(MAXMIXIOVAR), private :: IOVars
  character(len=strLen), private :: h5File,h5RunID
  character(len=strLen), dimension(nVars), private :: mixVarNames
  character(len=strLen), dimension(nVars), private :: mixUnitNames   

contains

  ! since Fortran doesn't allow arrays of strings of variable length
  ! use this function to define these arrays
  ! note, this needs to be called before use (e.g., in writeMIX and readMIX)
  subroutine initMIXNames()
    ! NOTE: these have to be in the same order as the
    ! variable enumerator in mixdefs
    mixVarNames(POT)           = "Potential"
    mixUnitNames(POT)          = "kV"
    mixVarNames(FAC)           = "Field-aligned current"
    mixUnitNames(FAC)          = "muA/m**2"
    mixVarNames(SIGMAP)        = "Pedersen conductance"
    mixUnitNames(SIGMAP)       = "S"
    mixVarNames(SIGMAH)        = "Hall conductance"
    mixUnitNames(SIGMAH)       = "S"
    mixVarNames(SOUND_SPEED)   = "Sound speed"
    mixUnitNames(SOUND_SPEED)  = "cm/s"
    mixVarNames(DENSITY)       = "Density"    
    mixUnitNames(DENSITY)      = "g/cm^3"
    mixVarNames(AVG_ENG)       = "Average energy"
    mixUnitNames(AVG_ENG)      = "keV"
    mixVarNames(NUM_FLUX)      = "Number flux"
    mixUnitNames(NUM_FLUX)     = "1/cm^2 s"
    mixVarNames(NEUTRAL_WIND)  = "Neutral wind"
    mixUnitNames(NEUTRAL_WIND) = "cm/s"
    mixVarNames(EFIELD)        = "Electric field"
    mixUnitNames(EFIELD)       = "mV/m" 
    mixVarNames(IM_EFLUX)      = "IM Energy flux"
    mixUnitNames(IM_EFLUX)     = "ergs/cm^2 s"
    mixVarNames(IM_EAVG)       = "IM average energy"
    mixUnitNames(IM_EAVG)      = "keV" ! add *1e-3 in rcm_mix_interface.F90
    mixVarNames(IM_IFLUX)      = "IM Energy flux proton"
    mixUnitNames(IM_IFLUX)     = "ergs/cm^2 s"
    mixVarNames(IM_IAVG)       = "IM average energy proton"
    mixUnitNames(IM_IAVG)      = "keV" ! add *1e-3 in rcm_mix_interface.F90
    mixVarNames(Z_NFLUX)       = "MHD number flux"
    mixUnitNames(Z_NFLUX)      = "1/cm^2 s"
    mixVarNames(Z_EAVG)        = "MHD average energy"
    mixUnitNames(Z_EAVG)       = "keV"
    mixVarNames(CRPOT)         = "Corotation Potential"
    mixUnitNames(CRPOT)        = "kV"
    mixVarNames(TPOT)          = "Total Potential"
    mixUnitNames(TPOT)         = "kV"
    mixVarNames(IM_GTYPE)      = "RCM grid type"
    mixUnitNames(IM_GTYPE)     = "0-1"
    mixVarNames(AUR_TYPE)      = "Auroral model type"
    mixUnitNames(AUR_TYPE)     = "0-3"
    mixVarNames(IM_BETA)       = "RCM beta"
    mixUnitNames(IM_BETA)      = "0-1"
    mixVarNames(IM_EDEN)       = "RCM electron density"
    mixUnitNames(IM_EDEN)      = "#/m^3"
    mixVarNames(IM_EPRE)       = "RCM electron pressure"
    mixUnitNames(IM_EPRE)      = "Pa"
    mixVarNames(IM_ENFLX)      = "IM Number flux"
    mixUnitNames(IM_ENFLX)     = "1/cm^2 s"
    mixVarNames(IM_INFLX)      = "IM Number flux proton"
    mixUnitNames(IM_INFLX)     = "1/cm^2 s"
    mixVarNames(DELTAE)        = "Mono potential drop"
    mixUnitNames(DELTAE)       = "kV"
    mixVarNames(IM_NPSP)       = "RCM plasmasphere density"
    mixUnitNames(IM_NPSP)      = "#/m^3"
  end subroutine initMIXNames

  subroutine initMIXIO(I,RunID,isRestart,nRes)
    type(mixIon_T),dimension(:),intent(inout) :: I
    character(len=*),optional, intent(in) :: RunID   ! these two things are needed when we're coupled with Gamera
    logical,optional, intent(in) :: isRestart
    integer,optional, intent(inout) :: nRes
    
    logical :: fExist
    real(rp), dimension(:,:),allocatable :: xc,yc
    integer :: n0

    !Setup remix output file
    h5File = trim(RunID) // ".mix.h5"
    if (present(RunID))  h5RunID = trim(RunID)

    inquire(file=h5File,exist=fExist)

    call initMIXNames() !Always do this ... maybe?

    if (isRestart .and. present(nRes)) then
      !Read the damn restart if you gotta
      call readMIXrestart(trim(RunID),nRes,I)
    endif

    if (isRestart .and. fExist) then
      !File is already here and have read restart, let's bounce
      return
    endif

    !If we're still here, then not restart or don't have output file
    !Either way, let's do it
    call CheckAndKill(h5File)
    !Reset IO chain
    call ClearIO(IOVars)

    call genOutGrid(I(NORTH)%G%x,I(NORTH)%G%y,xc,yc)
    
    ! save grid only for north
    call AddOutVar(IOVars,"X",xc,uStr="Ri")
    call AddOutVar(IOVars,"Y",yc,uStr="Ri")
   ! Save units of Ri and Rp as root attributes
    call AddOutVar(IOVars,"Ri_m",I(1)%rad_iono_m  )
    call AddOutVar(IOVars,"Rp_m",I(1)%rad_planet_m)


    call AddOutVar(IOVars,"UnitsID","ReMIX")
             
    !Write out the chain (to root)
    call WriteVars(IOVars,.false.,h5File)

  end subroutine initMIXIO

  subroutine writeMIX(I,Step,mjd,time)
    type(mixIon_T),dimension(:),intent(in) :: I    
    integer, intent(in) :: Step
    real(rp), optional, intent(in) :: time, mjd
    character(len=strLen) :: vStr

    integer :: v,h,n0
    character(len=strLen) :: gStr,uStr,hStr
    logical :: doDump = .true.
    real(rp) :: cpcp = 0.0

    !Reset IO chain
    call ClearIO(IOVars)

    ! I don't know how to handle the IO when the I object has more than one hemisphere
    ! because they're labeled by appending the hemisphere stirng to the data set name
    ! so kill it for now and worry about it when a case like this actually appears
    if (size(I)>2) then
       write(*,*) "writeMIX: Wrong hemisphere identifier. Stopping..."
       stop
    end if

    do h=1,size(I)
       ! hemisphere should be set up properly by now but still check just in case
       if (I(h)%St%hemisphere.eq.NORTH) then
          hStr = "NORTH"          
       else if (I(h)%St%hemisphere.eq.SOUTH) then
          hStr = "SOUTH"          
       else
          write(*,*) "writeMIX: Wrong hemisphere identifier. Stopping..."
          stop
       end if
       
       do v=1,nVars
          select case (v)
          case (POT)
             doDump = .true.
          case (FAC)
             doDump = .true.             
          case (SIGMAP)
             doDump = .true.             
          case (SIGMAH)
             doDump = .true.             
          case (SOUND_SPEED)
             doDump = .true.             
          case (DENSITY)
             doDump = .true.             
          case (AVG_ENG)
             doDump = .true.             
          case (NUM_FLUX)
             doDump = .true.    
          case (NEUTRAL_WIND) 
             doDump = .false.         
          case (EFIELD) 
             doDump = .false.         
          case (IM_EFLUX)
             doDump = .true.
          case (IM_EAVG)
             doDump = .true.
          case (IM_IFLUX)
             doDump = .true.
          case (IM_IAVG)
             doDump = .true.
          case (Z_NFLUX)
             doDump = .true.
          case (Z_EAVG)
             doDump = .true.         
          case (CRPOT) 
             doDump = .false.       
          case (TPOT) 
             doDump = .false.
          case (IM_GTYPE)
             doDump = .true.
          case (AUR_TYPE)
             doDump = .true.
          case (IM_BETA)
             doDump = .true.
          case (IM_EDEN)
             doDump = .true.
          case (IM_EPRE)
             doDump = .true.
          case (IM_ENFLX)
             doDump = .true.
          case (IM_INFLX)
             doDump = .true.
          case (DELTAE)
             doDump = .true.
          case (IM_NPSP)
             doDump = .true.
          case DEFAULT
             doDump = .false. ! only dump the variables explicitely set to be dumped above
          end select

          ! NOTE: assuming initMIXNames got called before
          vStr = trim(mixVarNames(v)) // " "//trim(hStr)
          uStr = trim(mixUnitNames(v))
          
          if (doDump) then
             call AddOutVar(IOVars,vStr,I(h)%St%Vars(:,2:,v))
             ! inelegantly specifying the units       
             n0 = FindIO(IOVars,vStr)
             IOVars(n0)%unitStr = uStr
          endif
       enddo
      ! add geographic latgeo and longeo
      vstr = trim("Latitude (Geographic)") // " "//trim(hStr)
      ustr = trim("radians")
      call AddOutVar(IOVars,vstr,I(h)%Ggeo%t,uStr)
      vstr = trim("Longitude (Geographic)") // " "//trim(hStr)
      call AddOutVar(IOVars,vstr,I(h)%Ggeo%p,uStr)
      ! add apex lat and lon
      vstr = trim("Latitude (Apex)") // " "//trim(hStr)
      call AddOutVar(IOVars,vstr,I(h)%Ggeo%t,uStr)
      vstr = trim("Longitude (Apex)") // " "//trim(hStr)
      call AddOutVar(IOVars,vstr,I(h)%Ggeo%p,uStr)
      ! add geographic X and Y
      vstr = trim("X (Geographic)") // " "//trim(hStr)
      ustr = trim("Ri")
      call AddOutVar(IOVars,vstr,I(h)%Ggeo%x,uStr)
      vstr = trim("Y (Geographic)") // " "//trim(hStr)
      call AddOutVar(IOVars,vstr,I(h)%Ggeo%y,uStr)
      ! add Apex X and Y
      vstr = trim("X (Apex)") // " "//trim(hStr)
      call AddOutVar(IOVars,vstr,I(h)%Ggeo%x,uStr)
      vstr = trim("Y (Apex)") // " "//trim(hStr)
      call AddOutVar(IOVars,vstr,I(h)%Ggeo%y,uStr)
    enddo

    ! now add time
    if (present(time)) call AddOutVar(IOVars,"time",time)
    if (present(mjd))  call AddOutVar(IOVars,"MJD",mjd)
    ! also add tilt
    call AddOutVar(IOVars,"tilt",I(NORTH)%St%tilt)

    ! add cpcp
    call AddOutVar(IOVars,"nCPCP",maxval(I(NORTH)%St%Vars(:,:,POT))-minval(I(NORTH)%St%Vars(:,:,POT)))
    call AddOutVar(IOVars,"sCPCP",maxval(I(SOUTH)%St%Vars(:,:,POT))-minval(I(SOUTH)%St%Vars(:,:,POT)))    
    
    !Write out the chain (to root)
    write(gStr,'(A,I0)') "Step#", Step
    call WriteVars(IOVars,.false.,h5File,gStr)
  end subroutine writeMIX

  subroutine readMIX(inH5,Step,mixIOobj)
    character(len=*), intent(in) :: inH5
    integer, intent(in) :: Step
    type(mixIO_T), intent(inout) :: mixIOobj
    
    integer,parameter :: hmsphrs(2) = [NORTH,SOUTH]
    real(rp),dimension(:,:),allocatable :: xc,yc
    integer :: h, v, n0
    integer :: dims(2)    
    
    character(len=strLen) :: gStr,hStr,uStr,vStr

    ! filling in var and unit names
    call initMIXNames()    

    call CheckFileOrDie(inH5,"Input H5 file does not exist.")

    !Reset IO chain
    call ClearIO(IOVars)
    
    ! read grid corners from root
    call AddInVar(IOVars,"X")
    call AddInVar(IOVars,"Y")
    call ReadVars(IOVars,.true.,inH5)

    dims = IOVars(1)%dims(1:2)

    if (.not.allocated(xc)) allocate(xc(dims(1),dims(2)))
    if (.not.allocated(yc)) allocate(yc(dims(1),dims(2)))

    call IOArray2DFill(IOVars,"X",xc)
    call IOArray2DFill(IOVars,"Y",yc)

    ! convert to original mix grid
    ! and fill in mixIOobj%x,y
    call genInGrid(xc,yc,mixIOobj%x,mixIOobj%y)

    ! now read from step

    !Reset IO chain
    call ClearIO(IOVars)
    
    call AddInVar(IOVars,"time")
    call AddInVar(IOVars,"MJD")
    call AddInVar(IOVars,"tilt")    

    if ( (.not.(size(hmsphrs).eq.2)) ) then
       write(*,*) 'Code is only implemented to do two hemispheres (north,south)'
       stop
    end if 

    do h=1,size(hmsphrs)
       if (h.eq.NORTH) then
          hStr = "NORTH"          
       else if (h.eq.SOUTH) then
          hStr = "SOUTH"          
       else
          write(*,*) "readMIX: Wrong hemisphere identifier. Stopping..."
          stop
       end if

       do v=1,nVars
          vStr = trim(mixVarNames(v)) // " "// trim(hStr)
          call AddInVar(IOVars,vStr)
       end do
    end do 
      

    write(gStr,'(A,I0)') "Step#", Step
    call ReadVars(IOVars,.true.,inH5,trim(gSTr)) ! note, this checks if step exists

    ! finally fill in the mixIO object for passing to calling program
    ! (mixIOobj%x,y already filled in above by genInGrid

    ! time & mjd
    mixIOobj%time = GetIOReal(IOVars,"time")
    mixIOobj%mjd  = GetIOReal(IOVars,"MJD")

    ! allow for no tilt in the restart file
    ! for backward compatibility
    if (IOVars(3)%isDone) then
       mixIOobj%tilt  = IOVars(3)%data(1)
    else
       mixIOobj%tilt  = 0
    end if

    ! allocate as necessary
    ! NOTE: we're assuming readMIX is not called in a loop
    ! and is just used for 1-step calculation
    ! thus allocating here
    ! also, remember dims below was taken from the x,y arrays
    ! those have one extra point in phi but THE SAME size in theta as our target array here
    ! (although the actual stored arrays had -1 point in the theta direction as well)
    ! since we cut out the pole but extrapolated the low lat boundary
    if (.not.allocated(mixIOobj%Vars)) allocate(mixIOobj%Vars(dims(1)-1,dims(2),nVars,size(hmsphrs)))
    mixIOobj%Vars = 0 ! and initialize to zero

    do h=1,size(hmsphrs)
       if (h.eq.NORTH) then
          hStr = "NORTH"          
       else if (h.eq.SOUTH) then
          hStr = "SOUTH"          
       end if

       do v=1,nVars
          vStr = trim(mixVarNames(v)) // " "// trim(hStr)
          n0 = FindIO(IOVars,vStr)

          ! check whether the variable exists in the file
          ! since we didn't necessarily dump all of them in writeMIX
          if (IOVars(n0)%isDone) then
             mixIOobj%Vars(:,2:dims(2),v,h) = reshape(IOVars(n0)%data,dims-1)
             ! fix pole
             mixIOobj%Vars(:,1,v,h) = sum(mixIOobj%Vars(:,2,v,h))/size(mixIOobj%Vars(:,2,v,h))
          end if
       end do
    end do

  end subroutine readMIX

  subroutine readMIXrestart(inH5,nRes,I)
    character(len=*), intent(in) :: inH5
    integer, intent(in) :: nRes
    type(mixIon_T), dimension(:), intent(inout) :: I
    
    integer,parameter :: hmsphrs(2) = [NORTH,SOUTH]
    integer :: h, v, n0
    integer :: dims(2)    
    
    character(len=strLen) :: gStr,hStr,uStr,vStr,nStr,h5Str

    ! filling in var and unit names
    call initMIXNames() 
    
    !Get number string
    if (nRes == -1) then
        nStr = "XXXXX"
    else
        write (nStr,'(I0.5)') nRes
    endif 
    
    h5Str = trim(inH5)//'.mix.Res.'// trim(nStr)//'.h5'
    write(*,*) "Restarting from: ", trim(h5Str)

    call CheckFileOrDie(h5Str,"Restart file not found ...")

    !Reset IO chain
    call ClearIO(IOVars)
    
    ! read grid corners from root
    call AddInVar(IOVars,"nRes")
    call ReadVars(IOVars,.false.,h5Str)

    ! record the nRes number
    if (ioExist(h5Str,"nRes")) then
        call ClearIO(IOVars)
        call AddInVar(IOVars,"nRes")
        call ReadVars(IOVars,.false.,h5Str)
        I(NORTH)%P%nRes = GetIOInt(IOVars,"nRes") + 1
        I(SOUTH)%P%nRes = GetIOInt(IOVars,"nRes") + 1
    else
        write (*,*) 'Remix now requires that nRes be in its restart files'
        write (*,*) '  in order to ensure consistent restarts.'
        write (*,*) 'Please either regenerate the remix restart files, or'
        write (*,*) '  manually add the nRes value.'
        write (*,*) 'File: ', trim(h5Str)
        stop
    endif

    ! now read from step

    !Reset IO chain
    call ClearIO(IOVars)
    

    if ( (.not.(size(hmsphrs).eq.2)) ) then
       write(*,*) 'Code is only implemented to do two hemispheres (north,south)'
       stop
    end if 

    do h=1,size(hmsphrs)
       if (h.eq.NORTH) then
          hStr = "NORTH"          
       else if (h.eq.SOUTH) then
          hStr = "SOUTH"          
       else
          write(*,*) "readMIX: Wrong hemisphere identifier. Stopping..."
          stop
       end if

       do v=1,nVars
          vStr = trim(mixVarNames(v)) // " "// trim(hStr)
          call AddInVar(IOVars,vStr)
       end do
    end do 
      
    !write(gStr,'(A,I0)') "Step#", Step
    call ReadVars(IOVars,.false.,h5Str) ! note, this checks if step exists


    !if (.not.allocated(mixIOobj%Vars)) allocate(mixIOobj%Vars(dims(1)-1,dims(2),nVars,size(hmsphrs)))
    !mixIOobj%Vars = 0 ! and initialize to zero
    I(NORTH)%St%Vars = 0
    I(SOUTH)%St%Vars = 0
    do h=1,size(hmsphrs)
       if (h.eq.NORTH) then
          hStr = "NORTH"          
       else if (h.eq.SOUTH) then
          hStr = "SOUTH"          
       end if

       do v=1,nVars
          
          vStr = trim(mixVarNames(v)) // " "// trim(hStr)
          n0 = FindIO(IOVars,vStr)
          dims = IOVars(n0)%dims(1:2)
          !write(*,*) trim(vStr),maxval(IOVars(n0)%data),minval(IOVars(n0)%data)
          ! check whether the variable exists in the file
          ! since we didn't necessarily dump all of them in writeMIX
          if (IOVars(n0)%isDone) then
             call IOArray2DFill(IOVars,vStr,I(h)%St%Vars(:,1:dims(2),v))
             write(*,*) trim(vStr),maxval(I(h)%St%Vars(:,:,v)),minval(I(h)%St%Vars(:,:,v))
          else
             write(*,*) "Variable not found in restart, skipping: ", trim(vStr)
          end if
       end do
    end do
    write(*,*) "Done readMIXrestart"

  end subroutine readMIXrestart

  subroutine writeMIXRestart(I,nRes,mjd,time)
    type(mixIon_T),dimension(:),intent(in) :: I    
    integer, intent(in) :: nRes
    real(rp), optional, intent(in) :: time, mjd
    character(len=strLen) :: vStr

    integer :: v,h,n0
    character(len=strLen) :: gStr,uStr,hStr,h5Str,nStr
    logical :: doDump = .true.
    real(rp) :: cpcp = 0.0
    character(len=strLen) :: ResF, tStr,lnResF !Name of restart file
    logical :: fExist

    if (nRes == -1) then
        write(*,*) 'nRes is -1 in writeMIXRestart, now die'
        stop
    else
        write (nStr,'(I0.5)') nRes
    endif

    write (ResF, '(A,A,A,A)') trim(h5RunID), '.mix.Res.', trim(nStr), '.h5'

    call CheckAndKill(ResF)

    !Reset IO chain
    call ClearIO(IOVars)

    ! I don't know how to handle the IO when the I object has more than one hemisphere
    ! because they're labeled by appending the hemisphere stirng to the data set name
    ! so kill it for now and worry about it when a case like this actually appears
    if (size(I)>2) then
       write(*,*) "writeMIX: Wrong hemisphere identifier. Stopping..."
       stop
    end if

    h5Str = ResF

    do h=1,size(I)
       ! hemisphere should be set up properly by now but still check just in case
       if (I(h)%St%hemisphere.eq.NORTH) then
          hStr = "NORTH"          
       else if (I(h)%St%hemisphere.eq.SOUTH) then
          hStr = "SOUTH"          
       else
          write(*,*) "writeMIX: Wrong hemisphere identifier. Stopping..."
          stop
       end if
       
       do v=1,nVars
          select case (v)
          case (POT)
             doDump = .true.
          case (FAC)
             doDump = .true.             
          case (SIGMAP)
             doDump = .true.             
          case (SIGMAH)
             doDump = .true.             
          case (SOUND_SPEED)
             doDump = .true.             
          case (DENSITY)
             doDump = .true.             
          case (AVG_ENG)
             doDump = .true.             
          case (NUM_FLUX)
             doDump = .true.    
          case (IM_EFLUX)
             doDump = .true.
          case (IM_EAVG)
             doDump = .true.
          case (IM_IFLUX)
             doDump = .true.
          case (IM_IAVG)
             doDump = .true.
          case (Z_NFLUX)
             doDump = .true.
          case (Z_EAVG)
             doDump = .true.   
          case (IM_GTYPE)
             doDump = .true.
          case (AUR_TYPE)
             doDump = .true.
          case (IM_BETA)
             doDump = .true.
          case (IM_EDEN)
             doDump = .true.
          case (IM_EPRE)
             doDump = .true.
          case (IM_ENFLX)
             doDump = .true.
          case (IM_INFLX)
             doDump = .true.
          case (DELTAE)
             doDump = .true.
          case (IM_NPSP)
             doDump = .true.
          case (NEUTRAL_WIND) 
             doDump = .false.
          case (EFIELD)
             ! we never compute it
             doDump = .false.
          case DEFAULT
             doDump = .false. ! only dump the variables explicitely set to be dumped above
          end select

          ! NOTE: assuming initMIXNames got called before
          vStr = trim(mixVarNames(v)) // " "//trim(hStr)
          uStr = trim(mixUnitNames(v))
          
          if (doDump) then
             call AddOutVar(IOVars,vStr,I(h)%St%Vars(:,:,v))
             ! inelegantly specifying the units       
             n0 = FindIO(IOVars,vStr)
             IOVars(n0)%unitStr = uStr
          endif
       enddo
    enddo

    ! now add time
    if (present(time)) call AddOutVar(IOVars,"time",time)
    if (present(mjd))  call AddOutVar(IOVars,"MJD",mjd)
    ! also add tilt
    call AddOutVar(IOVars,"tilt",I(NORTH)%St%tilt)

    ! add cpcp
    call AddOutVar(IOVars,"nCPCP",maxval(I(NORTH)%St%Vars(:,:,POT))-minval(I(NORTH)%St%Vars(:,:,POT)))
    call AddOutVar(IOVars,"sCPCP",maxval(I(SOUTH)%St%Vars(:,:,POT))-minval(I(SOUTH)%St%Vars(:,:,POT)))    
    
    ! add nres
    call AddOutVar(IOVars,"nRes",nRes)

    !Write out the chain (to root)
    call WriteVars(IOVars,.false.,h5Str)

    write (lnResF, '(A,A,A,A)') trim(h5RunID), ".mix.Res.", "XXXXX", ".h5"
    call MapSymLink(ResF,lnResF)

  end subroutine writeMIXRestart

  subroutine genOutGrid(x,y,xc,yc)  
    real(rp), dimension(:,:),intent(in) :: x,y
    real(rp), dimension(:,:),allocatable,intent(out) :: xc,yc ! with corners 1/2-cell shifted from original
    real(rp), dimension(:,:),allocatable :: xtmp,ytmp ! temporary arrays
    integer, dimension(2) :: dims
    integer :: Np, Nt

    dims = shape(x)
    Np = dims(1); Nt = dims(2)

    if (.not.allocated(xc)) allocate(xc(Np+1,Nt))
    if (.not.allocated(yc)) allocate(yc(Np+1,Nt))
    if (.not.allocated(xtmp)) allocate(xtmp(Np+2,Nt))
    if (.not.allocated(ytmp)) allocate(ytmp(Np+2,Nt))

    ! first fix up periodic
    ! overlap two lines and roll back one
    xtmp(2:Np+1,:) = x
    xtmp(1,:)    = x(Np,:)
    xtmp(Np+2,:) = x(1,:)
    ytmp(2:Np+1,:) = y
    ytmp(1,:)    = y(Np,:)
    ytmp(Np+2,:) = y(1,:)

    ! move to cell centers
    xc(:,1:Nt-1) = 0.25*(xtmp(2:Np+2,2:Nt)+xtmp(1:Np+1,2:Nt)+xtmp(2:Np+2,1:Nt-1)+xtmp(1:Np+1,1:Nt-1))
    yc(:,1:Nt-1) = 0.25*(ytmp(2:Np+2,2:Nt)+ytmp(1:Np+1,2:Nt)+ytmp(2:Np+2,1:Nt-1)+ytmp(1:Np+1,1:Nt-1))

    ! finally extrapolate below equatorward boundary
    xc(:,Nt) = 2*xc(:,Nt-1)-xc(:,Nt-2)
    yc(:,Nt) = 2*yc(:,Nt-1)-yc(:,Nt-2)      
  end subroutine genOutGrid

  ! mirror of genOutGrid
  ! take staggered grid from what's stored in the H5 file
  ! and turn it into the original remix grid
  subroutine genInGrid(xc,yc,x,y)
    real(rp), dimension(:,:),intent(in) :: xc,yc ! with corners 1/2-cell shifted from original    
    real(rp), dimension(:,:),allocatable,intent(out) :: x,y

    real(rp), dimension(:,:),allocatable :: xtmp,ytmp ! temporary arrays
    integer, dimension(2) :: dims
    integer :: Np, Nt

    dims = shape(xc)
    Np = dims(1)-1; Nt = dims(2)

    if (.not.allocated(x)) allocate(x(Np,Nt))
    if (.not.allocated(y)) allocate(y(Np,Nt))

    x(:,2:Nt) = 0.25*(xc(2:Np+1,2:Nt)+xc(1:Np,2:Nt)+xc(2:Np+1,1:Nt-1)+xc(1:Np,1:Nt-1))
    y(:,2:Nt) = 0.25*(yc(2:Np+1,2:Nt)+yc(1:Np,2:Nt)+yc(2:Np+1,1:Nt-1)+yc(1:Np,1:Nt-1))

    ! fix pole
    x(:,1) = 0
    y(:,1) = 0
  end subroutine genInGrid
  
end module mixio
