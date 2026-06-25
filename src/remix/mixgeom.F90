module mixgeom
  use mixdefs
  use mixtypes
  use mixinterp
  use earthhelper
  use shellGrid

  implicit none
  
  contains
      real(rp) function cosDipAngle(t)
          ! NOTE: sign needs to be flipped for the south. It's done in init_mix.
      real(rp),intent(in) :: t
      cosDipAngle = -2._rp*cos(t)/sqrt(1._rp+3._rp*cos(t)**2)
    end function cosDipAngle

    subroutine generate_uniformTP(Np,Nt,LowLatBoundary,HighLatBoundary,t,p)
      integer, intent(in) :: Np, Nt
      real(rp), intent(in) :: LowLatBoundary,HighLatBoundary ! in degrees
      integer :: i,j
      real(rp), dimension(:,:), allocatable, intent(out) :: t,p

      if (.not.allocated(t)) allocate(t(Np,Nt))
      if (.not.allocated(p)) allocate(p(Np,Nt))

      ! Note, we assume the number of vertices on the grid with
      ! periodic boundary cut out is Np
      ! That means that the number of cells in 2*pi is also Np
      do j=1,Np
         p(j,:) = 2*pi/Np*(j-1)
      end do

      do i=1,Nt
         t(:,i) = (LowLatBoundary-HighLatBoundary)/(Nt-1)*(i-1) + HighLatBoundary
      end do
    end subroutine generate_uniformTP


    subroutine init_grid_fromTP(G,t,p,isSolverGrid,isPeriodic)
      type(mixGrid_T),intent(out) :: G
      real(rp), dimension(:,:), intent(in) :: t,p
      logical, intent(in) :: isSolverGrid
      logical, optional, intent(in) :: isPeriodic
      integer :: i,j
      
      G%Np = size(p,1)
      G%Nt = size(p,2)

      if (present(isPeriodic)) G%isPeriodic = isPeriodic   ! otherwise, always periodic as is set in mix types
      
      if (.not.allocated(G%t)) allocate(G%t(G%Np,G%Nt))
      if (.not.allocated(G%p)) allocate(G%p(G%Np,G%Nt)) 

      ! interpolant
      ! True size (Np,Nt-1) (note, we include the cell between Np and 1 for periodic; 
      ! otherwise, it's unused and interpolant set to 0. below)
      if (.not.allocated(G%Interpolant)) allocate(G%Interpolant(G%Np,G%Nt,4,4)) 

      G%Interpolant = 0.0

      G%t = t 
      G%p = p

      ! define x, y
      if (.not.allocated(G%x)) allocate(G%x(G%Np,G%Nt))
      if (.not.allocated(G%y)) allocate(G%y(G%Np,G%Nt))

      G%x = sin(G%t)*cos(G%p)
      G%y = sin(G%t)*sin(G%p)

      ! calculate interpolant
      call mix_interpolant(G)

      if (isSolverGrid) call set_grid(G)
    end subroutine init_grid_fromTP

    subroutine init_grid(I,mixIOobj)
      type(mixIon_T),intent(inout) :: I
      type(mixIO_T),optional,intent(in) :: mixIOobj

      real(rp) :: highLatBoundary = 0.  ! always default the remix grid to start at the pole: that's the only thing that it knows how to do

      if (present(mixIOobj)) then
         call init_grid_fromXY(I%G,mixIOobj%x,mixIOobj%y,isSolverGrid=.true.)
      else
         ! for now always do uniform if grid not passed via mixIOobj
         call init_uniform(I%G,I%P%Np,I%P%Nt,I%P%LowLatBoundary*pi/180._rp,highLatBoundary*pi/180._rp,isSolverGrid=.true.,shGrO=I%shGr)
      endif
      call setD0(I%G,I%St%hemisphere)  ! pay attention: if the functional form is not axisymmetric, should make sure north and south are treated correctly.
    end subroutine init_grid
    
    subroutine init_uniform(G,Np,Nt,LowLatBoundary,HighLatBoundary,isSolverGrid, shGrO)
      type(mixGrid_T),intent(inout) :: G
      integer, intent(in) :: Np, Nt
      real(rp), intent(in) :: LowLatBoundary,HighLatBoundary ! in degrees
      logical, intent(in) :: isSolverGrid
      type(ShellGrid_T), intent(inout), optional :: shGrO
      real(rp), dimension(:,:), allocatable :: t,p
      real(rp), dimension(:), allocatable :: p_periodic

      call generate_uniformTP(Np,Nt,LowLatBoundary,HighLatBoundary,t,p)
      call init_grid_fromTP(G,t,p,isSolverGrid)
      if (present(shGrO)) then
         allocate(p_periodic(Np+1))
         p_periodic(1:Np) = p(:,1)
         p_periodic(Np+1) = 2*PI
         call GenShellGrid(shGrO, t(1,:), p_periodic, "REMIX", nGhosts=(/0,0,0,0/))
      endif
    end subroutine init_uniform

    subroutine init_grid_fromXY(G,x,y,isSolverGrid,isPeriodic)
      type(mixGrid_T),intent(inout) :: G
      real(rp), dimension(:,:), intent(in) :: x,y
      logical, intent(in) :: isSolverGrid
      logical, optional, intent(in) :: isPeriodic
      integer, dimension(2) :: dims

      ! set grid size
      dims = shape(x); G%Nt = dims(2); G%Np = dims(1)

      if (present(isPeriodic)) G%isPeriodic = isPeriodic   ! otherwise, always periodic as is set in mix types

      if (.not.allocated(G%x)) allocate(G%x(G%Np,G%Nt))
      if (.not.allocated(G%y)) allocate(G%y(G%Np,G%Nt))

      G%x = x  
      G%y = y

      if (.not.allocated(G%t)) allocate(G%t(G%Np,G%Nt))
      if (.not.allocated(G%p)) allocate(G%p(G%Np,G%Nt))

      ! interpolant
      ! True size (Np,Nt-1) (note, we include the cell between Np and 1 for periodic; 
      ! otherwise, it's unused and interpolant set to 0. below)      
      if (.not.allocated(G%Interpolant)) allocate(G%Interpolant(G%Np,G%Nt,4,4)) 

      G%Interpolant = 0.0

      ! define spherical angular coordinates
      G%t = asin(sqrt(G%x**2+G%y**2))
      G%p = modulo((atan2(G%y,G%x)+2*pi),(2*pi)) 
      ! note, this mangles phi at theta=0; Fix it, although we don't
      ! need it, since pole coordinates are never used
      G%p(:,1)=G%p(:,2)

      ! calculate interpolant
      call mix_interpolant(G)

      if (isSolverGrid) call set_grid(G)
    end subroutine init_grid_fromXY

    
    subroutine setD0(G,hemisphere)
      ! pay attention: if the functional form is not axisymmetric, should make sure north and south are treated correctly.
      type(mixGrid_T),intent(inout) :: G
      integer,intent(in) :: hemisphere
      integer :: i,j
      real(rp) :: r,L,phi
      
      ! allocate space for background density
      if (.not.allocated(G%D0)) allocate(G%D0(G%Np,G%Nt))
      
      do i=1,G%Nt
         do j=1,G%Np
            ! compute invariant latitude
            r = sqrt(G%x(j,i)**2+G%y(j,i)**2)   ! =cos(lambda)
            L = (1./r**2)*(RionE/REarth)*1.0e+6 ! correcting for the Ri/Re difference

            ! computing phi and flipping in the southern hemisphere
            if (hemisphere.eq.SOUTH) then
               phi = 2*pi-G%p(j,i)
            else
               phi = G%p(j,i)
            endif

            G%D0(j,i) = GallagherRP(L,phi)  ! use Gallagher from earthhelper
         enddo
      enddo
      
    end subroutine setD0
    

    ! NOTE, periodic boundary already cut out in mix2h5.py
    ! VGM 10142019: deprecating and renaming this function
    ! reserving the name for initializing grid from x,y arrays in remix convention (above)
    subroutine init_grid_fromXY_oldMIX(G,x,y,isSolverGrid)
      type(mixGrid_T),intent(inout) :: G
      real(rp), dimension(:,:), intent(in) :: x,y
      logical, intent(in) :: isSolverGrid
      integer, dimension(2) :: dims

      ! set grid size
      dims = shape(x); G%Nt = dims(2); G%Np = dims(1)

      if (.not.allocated(G%x)) allocate(G%x(G%Np,G%Nt))
      if (.not.allocated(G%y)) allocate(G%y(G%Np,G%Nt))

      G%x = x  ! note the transposes to conform to the new definition (Np,Nt)
      G%y = y

      if (.not.allocated(G%t)) allocate(G%t(G%Np,G%Nt))
      if (.not.allocated(G%p)) allocate(G%p(G%Np,G%Nt))

      ! interpolant
      if (.not.allocated(G%Interpolant)) allocate(G%Interpolant(G%Np,G%Nt,4,4)) ! True size (Np,Nt-1) (note, we include the cell between Np and 1)

      G%Interpolant = 0.0

      ! define spherical angular coordinates
      G%t = asin(sqrt(G%x**2+G%y**2))
      G%p = modulo((atan2(G%y,G%x)+2*pi),(2*pi)) 
      ! note, this mangles phi at theta=0; Fix it, although we don't
      ! need it, since pole coordinates are never used
      G%p(:,1)=G%p(:,2)

      ! calculate interpolant
      call mix_interpolant(G)

      if (isSolverGrid) call set_grid(G)
    end subroutine init_grid_fromXY_oldMIX

    !Simple MIX allocation routine to try to avoid memchecker flags
    subroutine AllocMIX(Q,Np,Nt)
        real(rp), dimension(:,:), allocatable, intent(inout) :: Q
        integer, intent(in) :: Np,Nt

        if (.not. allocated(Q)) allocate(Q(Np,Nt))
        Q(:,:) = 0.0
    end subroutine AllocMIX
    
    subroutine set_grid(G)
      type(mixGrid_T),intent(inout) :: G
      integer :: i,j

      call AllocMIX(G%r,G%Np,G%Nt)

      ! note allocating everything with the same size
      ! careful with unused points
      call AllocMIX(G%dp,G%Np,G%Nt) ! True size (Np-1,Nt)
      call AllocMIX(G%dt,G%Np,G%Nt) ! True size (Np,Nt-1)


      ! things needed for the solver and dependent only on the grid, not conductances
      ! same caveats about the sizes
      call AllocMIX(G%ft  ,G%Np,G%Nt)
      call AllocMIX(G%fp  ,G%Np,G%Nt)
      call AllocMIX(G%dtdt,G%Np,G%Nt)
      call AllocMIX(G%dpdp,G%Np,G%Nt)

      ! cyllindrical radius
      G%r = sqrt(G%x**2+G%y**2)

      ! note, keep explicit size on the LHS to avoid compiler-dependent
      ! problems down the road
      G%dt(:,1:G%Nt-1) = G%t(:,2:G%Nt)-G%t(:,1:G%Nt-1)  
      G%dp(1:G%Np-1,:) = G%p(2:G%Np,:)-G%p(1:G%Np-1,:)
      G%dp(G%Np,:) = modulo(G%p(1,:)-G%p(G%Np,:),2*pi)      ! fix up periodic

      ! note, unlike dp and dt above that are edge-centered, the things
      ! below are vertex centered; we just don't define them on the ends
      ! (e.g., ft(1,:) where we don't need them)
      G%ft(:,2:G%Nt) = 1.0D0/(G%dt(:,2:G%Nt)+G%dt(:,1:G%Nt-1))/sin(G%t(:,2:G%Nt))
      G%fp(2:G%Np,:) = 1.0D0/(G%dp(2:G%Np,:)+G%dp(1:G%Np-1,:)) ! note, dp(Np) defined above
      G%fp(1,:) = 1.0D0/(G%dp(1,:)+G%dp(G%Np,:))  ! fix up periodic

      ! this should work but since we don't use the last element (G%Nt) let's use the next line instead to avoid possibly dividing by zero (G%dt(:,2:G%Nt))
!      G%dtdt(:,2:G%Nt) = G%dt(:,2:G%Nt)/G%dt(:,1:G%Nt-1)-G%dt(:,1:G%Nt-1)/G%dt(:,2:G%Nt)
      G%dtdt(:,2:G%Nt-1) = G%dt(:,2:G%Nt-1)/G%dt(:,1:G%Nt-2)-G%dt(:,1:G%Nt-2)/G%dt(:,2:G%Nt-1)
      G%dpdp(2:G%Np,:) = G%dp(2:G%Np,:)/G%dp(1:G%Np-1,:)-G%dp(1:G%Np-1,:)/G%dp(2:G%Np,:)
      G%dpdp(1,:) = G%dp(1,:)/G%dp(G%Np,:)-G%dp(G%Np,:)/G%dp(1,:) ! fix up periodic

      ! dip angle: compute and store for access later
      call AllocMIX(G%cosd,G%Np,G%Nt)

      do i=1,G%Nt
         do j=1,G%Np
            G%cosd(j,i) = cosDipAngle(G%t(j,i))
         enddo
      enddo
    end subroutine set_grid

    subroutine flip_grid(G,GFpd,Rinner)
      type(mixGrid_T), intent(in) :: G
      type(mixGrid_T), intent(out) :: GFpd
      real(rp), dimension(:,:), allocatable :: xout,yout
      real(rp), intent(in) :: Rinner

      GFpd%Np = G%Np
      GFpd%Nt = G%Nt

      if (.not.allocated(xout)) allocate(xout(G%Np,G%Nt))    
      if (.not.allocated(yout)) allocate(yout(G%Np,G%Nt))    
      if (.not.allocated(GFpd%mask)) allocate(GFpd%mask(G%Np,G%Nt))
      if (.not.allocated(GFpd%t)) allocate(GFpd%t(G%Np,G%Nt))
      if (.not.allocated(GFpd%p)) allocate(GFpd%p(G%Np,G%Nt))
      GFpd%mask = 1  ! set interpolation mask to 1 by default everywhere


      ! map out to the magnetosphere
      where (sqrt(Rinner)*sin(G%t) > 1._rp-TINY )
         GFpd%mask=-1
         xout = 0._rp
         yout = 0._rp
      elsewhere
         xout = sqrt(Rinner)*sin(G%t)*cos(G%p)
         yout = sqrt(Rinner)*sin(G%t)*sin(G%p)
      end where

      ! note, the grid G that was created from gamera coordinates assumed R=1
      ! this doesn't matter as the interpolation happens in the (t,p) coordinate space
      ! and we've mapped xout,yout along dipole already above

      GFpd%t = acos(xout)
      GFpd%p = modulo(atan2(sqrt(1.-xout**2-yout**2),yout),2*pi)
    end subroutine flip_grid

end module mixgeom
