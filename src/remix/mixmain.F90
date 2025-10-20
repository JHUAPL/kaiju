module mixmain
  use xml_input
  use clocks
  use mixtypes
  use mixparams
  use mixgeom
  use dkaurora
  use mixconductance
  use mixstate
  use mixsolver
  use mixio

  implicit none

  contains

    subroutine init_mix(mixApp,hmsphrs,inpXML)
    ! subroutine init_mix(I,hmsphrs,optFilename,RunID,isRestart,mixIOobj,nRes,optIO)
      class(mixApp_T), intent(inout) :: mixApp
      integer, dimension(:), intent(in) :: hmsphrs       
      type(XML_Input_T), intent(in) :: inpXML

      type(mixIon_T),dimension(:),allocatable :: I ! I for ionosphere (is an array of 1 or 2 elements for north and south) or it can be artibrarily many, e.g., for different solves done in loop
      type(mixIO_T) :: mixIOobj
      ! character(len=*) :: RunID   ! these two things are needed when we're coupled with Gamera
      ! logical :: doRestart

      I = mixApp%ion
      mixIOobj = mixApp%IO

      ! call inpXML%Set_Val(RunID, "/Kaiju/Gamera/sim/runid","msphere")
      ! call inpXML%Set_Val(doRestart, "/Kaiju/Gamera/restart/doRes","false")

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

         call conductance_init(I(h)%conductance,I(h)%P,I(h)%G)
         call dragonking_init(I(h)%aurora,I(h)%P,I(h)%G)

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

    subroutine run_mix(I,tilt,doModelOpt,gcm)
      type(mixIon_T),dimension(:),intent(inout) :: I 
      type(gcm_T),optional,intent(inout) :: gcm
      real(rp),intent(in) :: tilt
      logical, optional, intent(in) :: doModelOpt  ! allow to change on the fly whether we use conductance model

      logical :: doModel=.true.   ! always default to xml input deck unless doModelOpt is present and on
      logical :: isRestart = .false.
      integer :: h,NumH

      if (present(doModelOpt)) doModel = doModelOpt
      if (present(gcm)) isRestart = gcm%isRestart

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

        if (.not. present(gcm)) then
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

        !if (present(gcm) .and. isRestart) then
        !  !write(*,*) "conductance: restart"
        !  !we read the conductance from file, so we're going to skip
        !  gcm%isRestart = .false.
        !  !write(*,*) "Get rePOT: ", maxval(I(h)%St%Vars(:,:,POT)),minval(I(h)%St%Vars(:,:,POT))
        !else
          call Tic("MIX-COND")
          ! Compute auroral precipitation flux regardless of coupling mode.
          ! Note conductance_euv is used inside dragonking_total.
          call conductance_euv(I(h)%conductance,I(h)%G,I(h)%St)
          call dragonking_total(I(h)%aurora,I(h)%G,I(h)%St,I(h)%conductance)

          if (present(gcm)) then
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
        !end if

      end do
    end subroutine run_mix

end module mixmain
