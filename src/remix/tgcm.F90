module gcminterp
  use mixdefs
  use mixtypes
  use gcmtypes
  use mixgeom
  use math
  use ioh5
  use mixio
  use kai2geo
  use clocks

  implicit none

  contains

    subroutine init_gcm_grid(gcm,ion)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(in) :: ion

      gcm%GEO%Coord = 'GEO'
      gcm%APEX%Coord = 'APEX'

      gcm%GEO%gcm2mix_nvar = 0 ! Neutral Winds (eventually)
      gcm%GEO%mix2gcm_nvar = 2 ! eng, flux

      gcm%APEX%gcm2mix_nvar = 2 ! SigmaP, sigmaH
      gcm%APEX%mix2gcm_nvar = 2 ! pot, flux (for auroral boundary)

      call create_gcm_grid(gcm%GEO,ion)
      call create_gcm_grid(gcm%APEX,ion)
      
    end subroutine init_gcm_grid

    subroutine create_gcm_grid(grid,ion)
      type(gcm_grid_T), intent(inout) :: grid
      type(mixIon_T),dimension(:),intent(in) :: ion

      integer :: Np, Nt, Nh, i,j,h,v
      Np = grid%nlon
      Nt = grid%nhlat
      Nh = GCMhemispheres

      ! Let's try something
      ! Build a list of index that says these are N and these are S
      ! Then use that to map to N/S in gcm
      grid%t2N = pack([(i,i=1,grid%nlat)], grid%inlat > 10)
      grid%t2S = pack([(i,i=1,grid%nlat)], grid%inlat < -10)
      grid%lonshift = findloc(grid%inlon(:grid%nlon-1),0.,1)-1

      ! Going to assume the hemispheres are symetrically sized for now
      if ( size(grid%t2N) /= size(grid%t2S) ) write(*,*) 'WE GOT ASYMMETRY'
      grid%nhlat = size(grid%t2N)

      if (.not.allocated(grid%gclat)) allocate(grid%gclat(grid%nlon,grid%nhlat,Nh))
      if (.not.allocated(grid%glon)) allocate(grid%glon(grid%nlon,grid%nhlat,Nh))

      ! Convert things from latitude to colatitude (and funky colat for south)
      do i=1,grid%nlon
        grid%gclat(i,:,GCMNORTH) = 90.-grid%inlat(grid%t2N)
        grid%gclat(i,:,GCMNORTH) = grid%gclat(i,grid%nhlat:1:-1,GCMNORTH) ! reverse order. Ascending.
        grid%gclat(i,:,GCMSOUTH) = 90.+grid%inlat(grid%t2S) !For mapping to work, we're going to use remix's funky southern colat
      enddo

      ! This complicated looking thing converts [-180,180] longitude into [0,360] and also shifts it so that longitude starts at 0
      do j=1,grid%nhlat
        grid%glon(:grid%nlon-1,j,GCMNORTH) = modulo(cshift(grid%inlon(:grid%nlon-1),grid%lonshift)+360.,360.)!modulo((atan2(G2%y,G2%x)+2*pi),(2*pi))
        grid%glon(:grid%nlon-1,j,GCMSOUTH) = modulo(cshift(grid%inlon(:grid%nlon-1),grid%lonshift)+360.,360.)
        ! cover periodic longitude point
        grid%glon(grid%nlon,j,GCMNORTH) = 360. ! no matter what, last point is periodic point
        grid%glon(grid%nlon,j,GCMSOUTH) = 360. ! hard coding the value of last point for now
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Now do remix mapping
      if (.not.allocated(grid%MyGrid)) allocate(grid%MyGrid(Nh))
      if (.not.allocated(grid%SM)) allocate(grid%SM(Nh))
      if (.not.allocated(grid%r2tMaps)) allocate(grid%r2tMaps(Nh))
      if (.not.allocated(grid%t2rMaps)) allocate(grid%t2rMaps(Nh))

      call init_grid_fromTP(grid%MyGrid(GCMNORTH),grid%gclat(:,:,GCMNORTH)*deg2rad,grid%glon(:,:,GCMNORTH)*deg2rad,isSolverGrid=.true.)
      call init_grid_fromTP(grid%MyGrid(GCMSOUTH),grid%gclat(:,:,GCMSOUTH)*deg2rad,grid%glon(:,:,GCMSOUTH)*deg2rad,isSolverGrid=.true.)

      ! Create input and output data storage arrays
      if (.not.allocated(grid%mixInput)) allocate(grid%mixInput(Nh,grid%gcm2mix_nvar))
      if (.not.allocated(grid%gcmInput)) allocate(grid%gcmInput(Nh,grid%gcm2mix_nvar))
      if (.not.allocated(grid%gcmOutput)) allocate(grid%gcmOutput(Nh,grid%mix2gcm_nvar))

      do h=1,Nh
        do v=1,grid%gcm2mix_nvar
          !write(*,*) "MIXCPL1: ",h,v,grid%nlon,grid%nhlat
          if (.not.allocated(grid%mixInput(h,v)%var)) allocate(grid%mixInput(h,v)%var(ion(h)%G%Np,ion(h)%G%Nt))
          if (.not.allocated(grid%gcmInput(h,v)%var)) allocate(grid%gcmInput(h,v)%var(grid%nlon,grid%nhlat))
        end do
        do v=1,grid%mix2gcm_nvar
          !write(*,*) "MIXCPL2: ",h,v,grid%nlon,grid%nhlat
          if (.not. allocated(grid%gcmOutput(h,v)%var)) allocate(grid%gcmOutput(h,v)%var(grid%nlon,grid%nhlat))
        end do
      enddo

    end subroutine create_gcm_grid

    subroutine process_gcmimport(gcm,ion)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      
      integer :: g

      call Tic("Processing")
      ! Split global GCM array into two hemispheres
      ! then shift the array so that longitude starts at 0
      do g=1,2
        select case (g)
        case (1)
          call mapGCM2MIX(gcm%GEO,ion,g)
        case (2)
          call mapGCM2MIX(gcm%APEX,ion,g)
        end select
      end do
      !Map the data to MIX grid
      !call mapGCM2MIX(gcm,ion)
      call Toc("Processing")

    end subroutine process_gcmimport

    subroutine transform_gcm_export(gcm,ion,g)
        type(gcm_grid_T),intent(inout) :: gcm
        type(mixIon_T),dimension(:),intent(inout) :: ion
        integer,intent(in) :: g

        integer :: h,ymod2,v

        if (gcm%mix2gcm_nvar .eq. 0) return

        ! The weird ymod here is to undo the funky southern hemisphere colat issue.
        do h=1,GCMhemispheres
        if (h == GCMSOUTH) then
            ymod2 = -1
        else
            ymod2 = 1
        endif
        !ymod2 = 1
        !write(*,*) "GEO -> SM START:",h,ymod2
        select case(g)
        case (1)
            call transform_grid(gcm%MyGrid(h),gcm%SM(h),iGEOtoSM,h,ym2=ymod2)
        case (2)
            call transform_grid(gcm%MyGrid(h),gcm%SM(h),iAPEXtoSM,h,ym2=ymod2)
        end select
        !write(*,*) "GEO -> SM END: ",h
        end do

        ! Map from mix grid to gcm grid
        call mapMIX2GCM(ion,gcm)

        if (.not. allocated(gcm%outvar2d)) allocate(gcm%outvar2d(gcm%nlat,gcm%nlon,gcm%mix2gcm_nvar))

        gcm%outvar2d = 0._rp
        ! Before we start, we collapse to 1 globe instead of 2 hemispheres
        do v=1,gcm%mix2gcm_nvar
          gcm%outvar2d(gcm%t2N(gcm%nhlat:1:-1),:gcm%nlon-1,v) = transpose(cshift(gcm%gcmOutput(GCMNORTH,v)%var(:gcm%nlon-1,:),-1*gcm%lonshift))
          gcm%outvar2d(gcm%t2S,:gcm%nlon-1,v) = transpose(cshift(gcm%gcmOutput(GCMSOUTH,v)%var(:gcm%nlon-1,:),-1*gcm%lonshift))
          gcm%outvar2d(gcm%t2N,gcm%nlon,v) = gcm%outvar2d(gcm%t2N,1,v)
          gcm%outvar2d(gcm%t2S,gcm%nlon,v) = gcm%outvar2d(gcm%t2S,1,v)
        end do

    end subroutine transform_gcm_export

    subroutine mapGCM2MIX(gcm,ion,g)
      type(gcm_grid_T),intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      integer,intent(in) :: g

      type(Map_T) :: Map
      
      real(rp), dimension(:,:), allocatable :: F
      integer :: h,v,i,ymod1

      ! Skip if not importing on this grid type
      if ( gcm%gcm2mix_nvar .eq. 0 ) return

      call Tic("Shifting")
      do v=1,gcm%gcm2mix_nvar
        gcm%gcmInput(GCMNORTH,v)%var = 0._rp
        do i=1,gcm%nhlat
          gcm%gcmInput(GCMNORTH,v)%var(:gcm%nlon-1,i) = cshift(gcm%invar2d(gcm%t2N(gcm%nhlat-i+1),:gcm%nlon-1,v),gcm%lonshift)
          gcm%gcmInput(GCMSOUTH,v)%var(:gcm%nlon-1,i) = cshift(gcm%invar2d(gcm%t2S(i),:gcm%nlon-1,v),gcm%lonshift)
          gcm%gcmInput(GCMNORTH,v)%var(gcm%nlon,i) = gcm%gcmInput(GCMNORTH,v)%var(1,i)
          gcm%gcmInput(GCMSOUTH,v)%var(gcm%nlon,i) = gcm%gcmInput(GCMSOUTH,v)%var(1,i)
        enddo
      end do
      call Toc("Shifting")

      call Tic("Mapping")
      ! The weird ymod here is to undo the funky southern hemisphere colat issue.
      do h=1,size(ion)
        if (h == GCMSOUTH) then
          ymod1 = -1
        else
          ymod1 = 1
        endif
        !ymod1 = 1
        !write(*,*) "SM -> GEO START:",h,ymod1
        select case(g)
        case (1)
          call transform_grid(ion(h)%G,ion(h)%Ggeo,iSMtoGEO,h,ym1=ymod1)
        case (2)
          call transform_grid(ion(h)%G,ion(h)%Gapx,iSMtoAPEX,h,ym1=ymod1)
        end select
        !write(*,*) "SM -> GEO END: ",h
      end do

      do h=1,GCMhemispheres
        select case(g)
        case (1)
        call mix_set_map(gcm%MyGrid(h),ion(h)%Ggeo,gcm%t2rMaps(h))
        case (2)
        call mix_set_map(gcm%MyGrid(h),ion(h)%Gapx,gcm%t2rMaps(h))
        end select
        do v=1,gcm%gcm2mix_nvar
          call mix_map_grids(gcm%t2rMaps(h),gcm%gcmInput(h,v)%var(:,:),F)
          gcm%mixInput(h,v)%var(:,:) = F
          !select case (v)
          !case (1)
          !  gcm%mixInput(h,v)%var(:,:) = F
          !case (2)
          !  gcm%mixInput(h,v)%var(:,:) = F
          !end select
        end do
      end do
      call Toc("Mapping")
    end subroutine mapGCM2MIX
    
    subroutine mapMIX2GCM(ion,gcm)
      type(mixIon_T),dimension(:),intent(inout) :: ion
      type(gcm_grid_T), intent(inout) :: gcm
      type(Map_T) :: Map
      
      real(rp), dimension(:,:), allocatable :: F,F1,F2
      integer :: h,v

      do h=1,GCMhemispheres
        call mix_set_map(ion(h)%G,gcm%SM(h),gcm%r2tMaps(h))
        do v=1,gcm%mix2gcm_nvar
          if (gcm%outlist(v) .eq. AVG_ENG) then
            ! If we want average energy, do the interpolation and mapping on energy flux instead
            call mix_map_grids(gcm%r2tMaps(h),ion(h)%St%Vars(:,:,AVG_ENG)*ion(h)%St%Vars(:,:,NUM_FLUX),F1)
            call mix_map_grids(gcm%r2tMaps(h),ion(h)%St%Vars(:,:,NUM_FLUX),F2)

            ! Catch any dvide by 0
            where (F2 > 0.0_rp)
                F = F1 / F2
            elsewhere
                F = 0.0_rp  ! or handle it as needed
            end where
          else
            call mix_map_grids(gcm%r2tMaps(h),ion(h)%St%Vars(:,:,gcm%outlist(v)),F)
          endif

          gcm%gcmOutput(h,v)%var(:,:) = F
        end do
      end do

    end subroutine mapMIX2GCM
    
    subroutine apply_gcm2mix(gcm,St,h)
      type(gcm_T),intent(in) :: gcm
      type(mixState_T),intent(inout) :: St
      integer :: v,h,g
      
      do g=1,2
        select case (g)
        case (1)
          do v=1,gcm%GEO%gcm2mix_nvar
            St%Vars(:,:,gcm%GEO%inlist(v)) = gcm%GEO%mixInput(h,v)%var(:,:)
          end do
        case (2)
          do v=1,gcm%APEX%gcm2mix_nvar
            St%Vars(:,:,gcm%APEX%inlist(v)) = gcm%APEX%mixInput(h,v)%var(:,:)
          end do
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

