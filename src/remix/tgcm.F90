module gcminterp
  use gcmtypes
  use math
  use ioh5
  use mixio
  use mixshellmap, only: mapGlobalShellVarToMix
  use kai2geo
  use clocks
  use shellInterp

  implicit none

  contains

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
      Np = gcm%nlon
      Nt = gcm%nhlat

      if (.not.allocated(gcm%gclat)) allocate(gcm%gclat(gcm%nlat))
      if (.not.allocated(gcm%glon)) allocate(gcm%glon(gcm%nlon))

      ! Ascending Lat starts in southern hemisphere
      ! Ascending Clat starts in northern hemisphere
      gcm%gclat = 90. - gcm%inlat
      gcm%gclat = gcm%gclat(gcm%nlat:1:-1) ! reverse order. Ascending.

      ! Check first and last point of clat to make sure it's 0 and 180
      if (.not. (gcm%gclat(1) < TINY .or. gcm%gclat(gcm%nlat) > 180.-TINY)) write(*,*) 'CLAT RANGE ERROR: ',gcm%gclat(1),gcm%gclat(gcm%nlat)
      
      ! find the point where longitude is 0, so we know how much we need to shift by
      gcm%lonshift = findloc(gcm%inlon(:gcm%nlon-1),0.,1)-1

      gcm%glon(:gcm%nlon-1) = modulo(cshift(gcm%inlon(1:gcm%nlon-1),gcm%lonshift)+360.,360.)
      gcm%glon(gcm%nlon) = 360.0

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

      write(*,*) "Created GCM grid with nlat,nlon,nhlat: ",gcm%nlat,gcm%nlon,gcm%nhlat
      write(*,*) "GCM Grid Check: ", shape(gcm%gcmInput(1)%data), gcm%grid%Nt, gcm%grid%Np
    end subroutine create_gcm_grid

    subroutine process_gcmimport(gcm,ion)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      
      integer :: g

      call Tic("Processing",.true.)
      ! Split global GCM array into two hemispheres
      ! then shift the array so that longitude starts at 0
      do g=1,2
        select case (g)
        case (1)
          if (gcm%GEO%mix2gcm_nvar > 0) then
            call mapGCM2MIX(gcm%GEO,ion,g)
          end if
        case (2)
          if (gcm%APEX%mix2gcm_nvar > 0) then
            call mapGCM2MIX(gcm%APEX,ion,g)
          end if
        end select
      end do
      call Toc("Processing",.true.)

    end subroutine process_gcmimport

    subroutine transform_gcm_export(gcm,g)
        type(gcm_grid_T),intent(inout) :: gcm
        integer,intent(in) :: g

        integer :: h,ymod2,v

        if (gcm%mix2gcm_nvar .eq. 0) return

        if (.not. allocated(gcm%outvar2d)) allocate(gcm%outvar2d(gcm%nlat,gcm%nlon,gcm%mix2gcm_nvar))

        gcm%outvar2d = 0._rp
        ! Before we start, we collapse to 1 globe instead of 2 hemispheres
        do v=1,gcm%mix2gcm_nvar
          ! Invert latitude to map colat -> lat
          ! undo the lonshift by applying it backwards
          gcm%outvar2d(gcm%nlat:1:-1,:gcm%nlon-1,v) = transpose(cshift(gcm%gcmOutput(v)%data(:gcm%nlon-1,:),-1*gcm%lonshift))
          gcm%outvar2d(:,gcm%nlon,v) = gcm%outvar2d(:,1,v)
        end do

    end subroutine transform_gcm_export

    subroutine mapGCM2MIX(gcm,ion,g)
      type(gcm_grid_T),intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      integer,intent(in) :: g

      integer :: v,i,mixVar

      ! Skip if not importing on this grid type
      if ( gcm%gcm2mix_nvar .eq. 0 ) return

      call Tic("Shifting")
      do v=1,gcm%gcm2mix_nvar
        gcm%gcmInput(v)%data = 0._rp
        do i=1,gcm%nlat
          gcm%gcmInput(v)%data(:gcm%nlon-1,i) = cshift(gcm%invar2d(gcm%nlat-i+1,:gcm%nlon-1,v),gcm%lonshift)
          gcm%gcmInput(v)%data(gcm%nlon,i) = gcm%gcmInput(v)%data(1,i)
        enddo
      end do
      call Toc("Shifting")

      call Tic("Mapping")
      do v=1,gcm%gcm2mix_nvar
        select case (trim(gcm%inlist(v)))
        case ("SIGMAP")
          mixVar = SIGMAP
        case ("SIGMAH")
          mixVar = SIGMAH
        case default
          write(*,*) "MIXCPL WARNING: Unknown GCM inlist variable: ", trim(gcm%inlist(v))
          cycle
        end select

        call mapGlobalShellVarToMix(gcm%grid, gcm%gcmInput(v), ion, mixVar)
      end do

      call Toc("Mapping")
    end subroutine mapGCM2MIX

    subroutine apply_gcm2mix(gcm,St,h)
      type(gcm_T),intent(in) :: gcm
      type(mixState_T),intent(inout) :: St
      integer :: v,h,g,mixVar
      
      do g=1,2
        select case (g)
        case (1)
          if (gcm%GEO%gcm2mix_nvar > 0) then
            do v=1,gcm%GEO%gcm2mix_nvar
              select case (trim(gcm%GEO%inlist(v)))
              case ("SIGMAP")
                mixVar = SIGMAP
              case ("SIGMAH")
                mixVar = SIGMAH
              case default
                cycle
              end select
              St%Vars(:,:,mixVar) = gcm%GEO%mixInput(h,v)%var(:,:)
            end do
          end if
        case (2)
          if (gcm%APEX%gcm2mix_nvar > 0) then
            do v=1,gcm%APEX%gcm2mix_nvar
              select case (trim(gcm%APEX%inlist(v)))
              case ("SIGMAP")
                mixVar = SIGMAP
              case ("SIGMAH")
                mixVar = SIGMAH
              case default
                cycle
              end select
              St%Vars(:,:,mixVar) = gcm%APEX%mixInput(h,v)%var(:,:)
            end do
          end if
        end select
      end do
      
    end subroutine apply_gcm2mix

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

