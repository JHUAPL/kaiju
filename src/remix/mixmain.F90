module mixmain
  use xml_input
  use clocks
  use mixtypes
  use mixparams
  use mixgeom
  use mixconductance
  use mixstate
  use mixsolver
  use mixio

  implicit none

  contains

    subroutine init_mix(I,hmsphrs,optFilename,RunID,isRestart,mixIOobj,nRes,optIO)
      type(mixIon_T),dimension(:),allocatable,intent(inout) :: I ! I for ionosphere (is an array of 1 or 2 elements for north and south) or it can be artibrarily many, e.g., for different solves done in loop
      integer, dimension(:), intent(in) :: hmsphrs       
      character(len=*), optional, intent(in) :: optFilename
      character(len=*),optional, intent(in) :: RunID   ! these two things are needed when we're coupled with Gamera
      logical,optional, intent(in) :: isRestart
      type(mixIO_T),optional, intent(in) :: mixIOobj  
      integer,optional,intent(in) :: nRes
      logical,optional,intent(in) :: optIO
      integer :: mixnRes    
      logical :: doRestart
      integer :: h
      logical :: doIO ! do IO initialziaiton

      if (.not.allocated(I)) allocate(I(size(hmsphrs)))
      
      I%rad_iono_m = RIonE*1.e+6 ! Default to RIonE. To change, overwrite directly after the init_mix call
      I%rad_planet_m = REarth
      if (present(isRestart)) doRestart = isRestart
      if (present(optIO)) then
          doIO = optIO
      else
          doIO = .true. ! IO enabled default
      endif

      do h=1,size(I)
         I(h)%St%hemisphere = hmsphrs(h)
         if(present(optFilename)) then
            call initMIXParams(I(h)%P, optFilename)
         else
            call initMIXParams(I(h)%P)
         endif

         if (present(mixIOobj)) then
            call init_grid(I(h),mixIOobj)
            call init_state(I(h)%G,I(h)%St,mixIOobj%vars(:,:,:,h)) ! passing only the hemisphere variables            
         else
            call init_grid(I(h))
            call init_state(I(h)%G,I(h)%St)            
         end if
         
         call conductance_init(I(h)%conductance,I(h)%P,I(h)%G)

         ! check that hemisphere makes sense.
         if ((I(h)%St%hemisphere.ne.NORTH).and.(I(h)%St%hemisphere.ne.SOUTH)) then
            write(*,*) 'Hemisphere is set to an unallowable value: ',I(h)%St%hemisphere
            write(*,*) 'Stopping...'
            stop
         end if

         ! initialize solver for each ionosphere instance
         call init_solver(I(h)%P,I(h)%G,I(h)%S)

         ! flip the sign of cosd for south in this init_mix function
         ! -- which is the only place where it makes sense.  Could've
         ! done it in the mixgeom.F90 functions but they don't share
         ! the info of the entire mixIon_T object.
         if (I(h)%St%hemisphere.eq.SOUTH) I(h)%G%cosd = -I(h)%G%cosd
      end do

      call initMIXNames() !Just for luck?
      mixnRes = I(1)%P%nRes
      !NOTE: mixnres comes from I(1)%P%nRes and is the one that should be used not subroutine argument nRes!
      !Why is nRes even passed to this function then, you might ask.
      !Good question hypothetical question person, good question.

      ! initialize the mix dump file (AND RESTART)
      if (present(isRestart) .and. doRestart) then
        !call xmlInp%Set_Val(mixnRes ,"/Kaiju/gamera/restart/nRes" ,-1)
        
        !Just restart here and be done with it
        call readMIXrestart(trim(RunID),mixnRes,I)
      endif

      if (doIO) then
        call initMIXIO(I,RunID,isRestart)
      endif
      
    end subroutine init_mix

    subroutine updateF107(I,F107)
        implicit none
        type(mixIon_T),dimension(:),intent(inout) :: I 
        real(rp), intent(in) :: F107
        integer :: h
    
        do h=1,size(I)
            I(h)%conductance%f107 = F107
        enddo
    end subroutine updateF107


    subroutine get_potential(I)
      type(mixIon_T),intent(inout) :: I

      I%St%Vars(:,:,POT) = reshape(I%S%solution,[I%G%Np,I%G%Nt])*(I%rad_iono_m*1.e-6)**2*1.D3 ! in kV
    end subroutine get_potential

    subroutine run_mix(I,tilt,doModelOpt,gcm,mjd)
      type(mixIon_T),dimension(:),intent(inout) :: I 
      type(gcm_T),optional,intent(inout) :: gcm
      real(rp),intent(in) :: tilt
      logical, optional, intent(in) :: doModelOpt  ! allow to change on the fly whether we use conductance model
      real(rp), optional, intent(in) :: mjd ! used to calculate the remix GEO grid as a function of time

      logical :: doModel=.true.   ! always default to xml input deck unless doModelOpt is present and on
      logical :: isRestart = .false.
      integer :: h,NumH

      if (present(doModelOpt)) doModel = doModelOpt
      if (present(gcm)) isRestart = gcm%isRestart
      if (present(mjd)) call MJDRecalc(mjd)

      NumH = size(I)

      !Removing OMP loop due to weird runtime errors
      do h=1,NumH
        if (I(h)%St%hemisphere.eq.NORTH) then
          I(h)%St%tilt = tilt
        else
          I(h)%St%tilt = -tilt
        end if

        if (doModel) then
          I(h)%conductance%const_sigma = I(h)%P%const_sigma         
        else
          I(h)%conductance%const_sigma = .true.            
        end if

        if (present(mjd) .and. .not. present(gcm)) then
          !write(*,*) "GRID TRANSFORM!"
          if (I(h)%St%hemisphere .eq. NORTH) then
            call transform_grid(I(h)%G,I(h)%Ggeo,iSMtoGEO,h,ym1=1)
            call transform_grid(I(h)%G,I(h)%Gapx,iSMtoAPEX,h,ym1=1)
          else
            call transform_grid(I(h)%G,I(h)%Ggeo,iSMtoGEO,h,ym1=-1)
            call transform_grid(I(h)%G,I(h)%Gapx,iSMtoAPEX,h,ym1=-1)
          endif
          !write(*,*) "GRID TRANSFORM!!"
        end if

        if (present(gcm) .and. isRestart) then
          !write(*,*) "conductance: restart"
          !we read the conductance from file, so we're going to skip
          gcm%isRestart = .false.
          !write(*,*) "Get rePOT: ", maxval(I(h)%St%Vars(:,:,POT)),minval(I(h)%St%Vars(:,:,POT))
        else
          call Tic("MIX-COND")
          if (present(gcm)) then
            !write(*,*) 'doGCM!'
            call conductance_total(I(h)%conductance,I(h)%G,I(h)%St,gcm,h)
          else
            !write(*,*) "conductance: total"
            call conductance_total(I(h)%conductance,I(h)%G,I(h)%St)
          end if
          call Toc("MIX-COND")

          call Tic("MIX-SOLVE")
          call run_solver(I(h)%P,I(h)%G,I(h)%St,I(h)%S)
          call Toc("MIX-SOLVE")

          call Tic("MIX-POT")
          call get_potential(I(h))
          call Toc("MIX-POT")
        end if

      end do
    end subroutine run_mix

end module mixmain
