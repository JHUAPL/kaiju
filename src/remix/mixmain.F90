module mixmain
  use xml_input
  use clocks
  use mixtypes
  use mixparams
  use mixgeom
  use mixstate
  use mixsolver
  use mixio

  implicit none

  contains

    subroutine init_mix(I,hmsphrs,inpXML,mixIOobj)
    ! subroutine init_mix(I,hmsphrs,optFilename,RunID,isRestart,mixIOobj,nRes,optIO)
    ! subroutine init_mix(mixApp,hmsphrs,inpXML)
      type(mixIon_T),dimension(:),allocatable,intent(inout) :: I
      integer, dimension(:), intent(in) :: hmsphrs       
      type(XML_Input_T), intent(in) :: inpXML
      type(mixIO_T), intent(in) :: mixIOobj  

      integer :: h

      if (.not.allocated(I)) allocate(I(size(hmsphrs)))

      I%rad_iono_m = RIonE*1.e+6 ! Default to RIonE. To change, overwrite directly after the init_mix call
      I%rad_planet_m = REarth

      do h=1,size(I)
         I(h)%St%hemisphere = hmsphrs(h)
         call initMIXParams(I(h)%P, inpXML)

         if (I(h)%P%init_from_file) then
            call init_grid(I(h),mixIOobj)
            call init_state(I(h)%G,I(h)%St,mixIOobj%vars(:,:,:,h), shGrO=I(h)%shGr) ! passing only the hemisphere variables            
         else
            call init_grid(I(h))
            call init_state(I(h)%G,I(h)%St, shGrO=I(h)%shGr)            
         end if

        !  call conductance_init(I(h)%conductance,I(h)%P,I(h)%G)
        !  call dragonking_init(I(h)%aurora,I(h)%P,I(h)%G)

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

      ! initialize the mix dump file (AND RESTART)
      ! if (doRestart) then        
      !   ! call readMIXrestart(trim(RunID),mixnRes,I)
      !   call mixApp%ReadRestart(mixApp,trim(RunID),I(1)%P%nRes)
      ! endif

      ! if (doIO) then
      !   ! call initMIXIO(I,RunID,isRestart)
      !   call mixApp%InitIO(mixApp,inpXML)
      ! endif
      
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

    subroutine run_mix(I)
      type(mixIon_T),dimension(:),intent(inout) :: I 
      integer :: h,NumH
      logical :: isGoodSolve

      NumH = size(I)
        
      do h=1,NumH
          call Tic("MIX-SOLVE")
          call run_solver(I(h)%P,I(h)%G,I(h)%St,I(h)%S,isGoodSolveO=isGoodSolve)
          call Toc("MIX-SOLVE")

          if (.not. isGoodSolve) then
            Write(*,*)"Didn't get a good solve, writing state and dying"
            call writeMIXCRASH(I)
            stop
          endif

          call Tic("MIX-POT")
          call get_potential(I(h))
          call Toc("MIX-POT")
      end do

    end subroutine run_mix

end module mixmain
