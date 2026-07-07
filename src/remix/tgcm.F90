module gcminterp
  use gcmtypes
  use math
  use ioh5
  use mixio
  use kai2geo
  use clocks
  use shellInterp

  implicit none

  contains

      subroutine init_gcm(gcm,isRestart)
      type(gcm_T),intent(inout) :: gcm
      logical, intent(in) :: isRestart

      write(*,*) "start init_gcm"

      call init_gcm_grid(gcm)

      write(*,*) "done init_gcm_grid"

      if (.not. allocated(gcm%GEO%outlist)) allocate(gcm%GEO%outlist(gcm%GEO%mix2gcm_nvar))
      if (.not. allocated(gcm%APEX%outlist)) allocate(gcm%APEX%outlist(gcm%APEX%mix2gcm_nvar))
      
      if (gcm%isWACCMX) then
        gcm%APEX%outlist(1) = "POT"
        gcm%APEX%outlist(2) = "AVG_ENG"
        gcm%APEX%outlist(3) = "NUM_FLUX"
      else ! is TIEGCM
        gcm%APEX%outlist(1) = "POT"
        gcm%APEX%outlist(2) = "NUM_FLUX"
        gcm%GEO%outlist(1) = "AVG_ENG"
        gcm%GEO%outlist(2) = "NUM_FLUX"
      endif

      if ( .not. allocated(gcm%GEO%inlist)) allocate(gcm%GEO%inlist(gcm%GEO%gcm2mix_nvar))
      if ( .not. allocated(gcm%APEX%inlist)) allocate(gcm%APEX%inlist(gcm%APEX%gcm2mix_nvar))
      gcm%APEX%inlist(1) = "SIGMAP"
      gcm%APEX%inlist(2) = "SIGMAH"

      gcm%isRestart = isRestart

      write(*,*) "done init_gcm_mpi"
      
    end subroutine init_gcm

    subroutine init_gcm_grid(gcm)
      type(gcm_T), intent(inout) :: gcm

      if (gcm%GEO%mix2gcm_nvar > 0) call create_gcm_grid(gcm%GEO, "GEO")
      if (gcm%APEX%mix2gcm_nvar > 0) call create_gcm_grid(gcm%APEX, "APEX")
      
    end subroutine init_gcm_grid

    subroutine init_gcm_vars(gcm)
      type(gcm_T), intent(inout) :: gcm

      if (gcm%isWACCMX) then
        write(*,*) "MIXCPL Coupling to WACCMX: ",gcm%isWACCMX
        gcm%GEO%gcm2mix_nvar = 0 ! Neutral Winds(eventually)
        gcm%GEO%mix2gcm_nvar = 0 ! Nothing passed on GEO

        gcm%APEX%gcm2mix_nvar = 2 ! SigmaP, sigmaH
        gcm%APEX%mix2gcm_nvar = 3 ! pot, eng, flux
      else ! is TIEGCM
        write(*,*) "MIXCPL Coupling to TIEGCM: ", (.not. gcm%isWACCMX)
        gcm%GEO%gcm2mix_nvar = 0 ! Neutral Winds (eventually)
        gcm%GEO%mix2gcm_nvar = 2 ! eng, flux

        gcm%APEX%gcm2mix_nvar = 2 ! SigmaP, sigmaH
        gcm%APEX%mix2gcm_nvar = 2 ! pot, flux (for auroral boundary)

      endif

    end subroutine init_gcm_vars


    subroutine create_gcm_grid(gcm,coord)
      type(gcm_grid_T), intent(inout) :: gcm
      character(len=*), intent(in) :: coord

      integer :: nGhosts(4)
      integer :: Np, Nt, Nh, i,j,h,v

      ! nlon includes periodic ghost

      if (.not.allocated(gcm%gclat)) allocate(gcm%gclat(gcm%nlat))
      if (.not.allocated(gcm%glon)) allocate(gcm%glon(gcm%nlon))

      ! Ascending Lat starts in southern hemisphere
      ! Ascending Clat starts in northern hemisphere
      gcm%gclat = 90. - gcm%inlat
      gcm%gclat = gcm%gclat(gcm%nlat:1:-1) ! reverse order. Ascending.

      ! Check first and last point of clat to make sure it's 0 and 180
      if (.not. (gcm%gclat(1) < TINY .and. gcm%gclat(gcm%nlat) > 180.-TINY)) write(*,*) 'CLAT RANGE ERROR: ',gcm%gclat(1),gcm%gclat(gcm%nlat)
      
      ! find the point where longitude is 0, so we know how much we need to shift by
      gcm%lonshift = minloc(abs(gcm%inlon(:gcm%nlon-1)), 1) - 1

      gcm%glon(:gcm%nlon-1) = modulo(cshift(gcm%inlon(1:gcm%nlon-1),gcm%lonshift)+360.,360.)

      if (abs(gcm%glon(1)) > 0.01) write(*,*) "WARNING: GCM lonshift is not near 0.0 degrees: ",gcm%glon(1)

      gcm%glon(1) = 0.0_rp ! Clamp to 0.0 just in case
      gcm%glon(gcm%nlon) = 360.0 ! Clamp final value exactly to 360.0

      write(*,*) "GCM Grid GLON: ", gcm%glon
      write(*,*) "GCM Grid GCLAT: ", gcm%gclat

      gcm%glon = gcm%glon * deg2rad
      gcm%gclat = gcm%gclat * deg2rad

      nGhosts = (/0,0,0,0/)
      call GenShellGrid(gcm%grid, gcm%gclat, gcm%glon, "GCM", nGhosts=nGhosts, coord=coord)

      if (.not. allocated(gcm%gcmInput)) allocate(gcm%gcmInput(gcm%gcm2mix_nvar))
      do i=1,gcm%gcm2mix_nvar
        call initShellVar(gcm%grid, SHGR_CORNER, gcm%gcmInput(i))
      end do

      if (.not. allocated(gcm%gcmOutput)) allocate(gcm%gcmOutput(gcm%mix2gcm_nvar))
      do i=1,gcm%mix2gcm_nvar
        call initShellVar(gcm%grid, SHGR_CORNER, gcm%gcmOutput(i))
      end do

      write(*,*) "Created GCM grid with nlat,nlon,nhlat: ",gcm%nlat,gcm%nlon
      write(*,*) "GCM Grid Check: ", shape(gcm%gcmOutput(1)%data), gcm%grid%Nt, gcm%grid%Np
      write(*,*) "GCM Lower bounds: ", lbound(gcm%gcmOutput(1)%data)
      write(*,*) "GCM Upper bounds: ", ubound(gcm%gcmOutput(1)%data)
    end subroutine create_gcm_grid

    subroutine transform_gcm_export(gcm)
        type(gcm_grid_T),intent(inout) :: gcm

        integer :: v

        if (gcm%mix2gcm_nvar .eq. 0) return

        if (.not. allocated(gcm%outvar2d)) allocate(gcm%outvar2d(gcm%nlat,gcm%nlon,gcm%mix2gcm_nvar))

        gcm%outvar2d = 0._rp
        do v=1,gcm%mix2gcm_nvar
          ! Invert latitude to map colat -> lat
          ! undo the lonshift by applying it backwards
          gcm%outvar2d(gcm%nlat:1:-1,:gcm%nlon-1,v) = cshift(gcm%gcmOutput(v)%data(:,:gcm%nlon-1),-1*gcm%lonshift, dim=2)
          gcm%outvar2d(:,gcm%nlon,v) = gcm%outvar2d(:,1,v)
        end do

    end subroutine transform_gcm_export

    subroutine process_gcm_import(gcm)
        type(gcm_grid_T),intent(inout) :: gcm

        integer :: v

        ! gcmInput lives in colat space on the Shell Grid
        ! invar2d is the direct lat space from GCM exchange
        do v=1,gcm%gcm2mix_nvar
          gcm%gcmInput(v)%data = 0._rp
        
          gcm%gcmInput(v)%data(gcm%nlat:1:-1,:gcm%nlon-1) = cshift(gcm%invar2d(:,:gcm%nlon-1,v),gcm%lonshift, dim=2)
          gcm%gcmInput(v)%data(:,gcm%nlon) = gcm%gcmInput(v)%data(:,1)
        end do

    end subroutine process_gcm_import

!    subroutine calculate_auroral_boundary(ion,bc)
! Fix this to calculate auroral boundary on remix's grid if necessary
!      type(mixIon_T),dimension(:),intent(in) :: ion
!      real(rp),allocatable,dimension(:,:),intent(out) :: bc
!
!      integer :: v,i,j
!      real(rp) :: nfluxllb = 1.0e7
!
!      if (.not. allocated(bc)) allocate(bc(gcm%nlon,2))
!
!        do v=1,gcm%mix2gcm_nvar
!            if (gcm%outlist(v) .eq. NUM_FLUX) then
!
!            bc = 0.D0
!
!            do j=1,gcm%nlon   ! NORTH
!                do i=gcm%nlat/2+1,gcm%nlat
!                if(gcm%outvar2d(i,j,v)>=nfluxllb) exit ! Find the lowest lat where numflux is above 1e6/cm^2/s
!                enddo
!                i = min(i,gcm%nlat)
!                bc(i,1) = max(90.-gcm%inlat(j,i),15.) ! aurllbj is co-lat.
!            enddo
!
!            do j=1,gcm%nlon   ! SOUTH
!                do i=gcm%nlat/2,1,-1
!                if(gcm%outvar2d(i,j,v)>=nfluxllb) exit ! Find the lowest lat where numflux is above 1e6/cm^2/s
!                enddo
!                i = max(i,1)
!                bc(i,2) = max(90.+gcm%inlat(j,i),15.) ! aurllbj is co-lat from south pole! Backwards.
!            enddo
!
!        endif
!      enddo
!
!    end subroutine

end module gcminterp

