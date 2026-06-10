! Miscellaneous routines used by the gamera<->remix coupling modules

module mixinterfaceutils
    use mixdefs
    use volttypes
    use shellInterp

    implicit none

    contains

  ! This function gets the MHD grid from MHD directly
  subroutine mix_mhd_grid(mhdg,t,p,tFpd,pFpd,Rinner)
    real(rp), dimension(:,:,:), intent(in) :: mhdg ! MHD grid cell coords
    real(rp), dimension(:,:), allocatable, intent(out) :: t,p,tFpd,pFpd
    real(rp), dimension(:,:), allocatable :: r
    real(rp), intent(out) :: Rinner
    real(rp), dimension(:,:), allocatable :: dipRatio
    integer :: nj, nk2

    ! note, mhdg here already has the i-dimension stripped
    ! so, dimension 1 = j and dimension 2 = k below
    nj = size(mhdg,1); nk2 = size(mhdg,2) ! nk2+1 for Psi shells but keep the notation for brevity

    ! Also define t, p variables in the gamera coordinate space (x-axis is the spherical axis)
    ! "Fpd" suffix stands for "flipped"
    if (.not.allocated(t)) allocate(t(nj,nk2))
    if (.not.allocated(p)) allocate(p(nj,nk2))
    if (.not.allocated(r)) allocate(r(nj,nk2))
    if (.not.allocated(dipRatio)) allocate(dipRatio(nj,nk2))

    if (.not.allocated(tFpd)) allocate(tFpd(nk2,nj))
    if (.not.allocated(pFpd)) allocate(pFpd(nk2,nj))

    ! NOTE: this is already in units of Rion
    r = sqrt(mhdg(:,:,XDIR)**2+mhdg(:,:,YDIR)**2+mhdg(:,:,ZDIR)**2)

    ! capture for the rare (but possible!) case where the psiShell
    ! we're mapping to is below the ionosphere
    ! this may happen, e.g., for the bottom-most ghost cell inside the Gamera inner boundary.
    ! In this case, we simply set theta=pi/2 so mix_map_grids (in mixinterp.F90) will simply
    ! extrapolate equatorward from the MIX low lat boundary
    dipRatio = sqrt(mhdg(:,:,XDIR)**2+mhdg(:,:,YDIR)**2)/r**1.5
    where( dipRatio < 1. )
       t = asin(dipRatio)
    elsewhere
       t = PI/2.
    end where

    p = modulo(atan2(mhdg(:,:,YDIR),mhdg(:,:,XDIR)),2*pi)
    Rinner = sum(r)/size(r)

    ! NOTE, transposition is necessary to make sure the phi coordinate goes first
    ! because all the codes, particulalry, interpolation (mixinterp) assume that
    tFpd = transpose(acos(mhdg(:,:,XDIR)/r))
    pFpd = transpose(modulo(atan2(mhdg(:,:,ZDIR),mhdg(:,:,YDIR)),2*pi))

    ! Fixes for south
    if (mhdg(nj/2,nk2/2,ZDIR).lt.0) then ! pick the pole and see if
       ! z<0. this works even if nk2
       ! is nk2+1 (for PsiShells)
       ! For south, phi obtained this way goes from pi to 2pi.  so we
       ! subtract pi because we always assume one hemisphere when
       ! mapping
       pFpd = pFpd-pi ! this is phi in the lfm/gamera space

       ! FIXME: make sure this actually works (comment saved for historical reasons)
       ! we did and it does
       p = modulo(atan2(-mhdg(:,:,YDIR),mhdg(:,:,XDIR)),2*pi)
    end if
  end subroutine mix_mhd_grid


   subroutine mixVarToVoltron(mixApp, mixvarID, voltGrid, varOut)
      !! Map a given remix variable onto the voltron grid
      !! Convert certain remix quantities to Voltron State variables
      class(mixApp_T  )   , intent(inout) :: mixApp
      integer             , intent(in)    :: mixvarID
      type(ShellGrid_T)   , intent(in)    :: voltGrid
      type(ShellGridVar_T), intent(inout) :: varOut
         !! (w.r.t. voltGrid) Potential we write to 

      !real(rp), dimension(mixApp%ion(NORTH)%shGr%Nt,mixApp%ion(NORTH)%shGr%Np) :: tmpPot
      integer :: iLat
      ! South grid stuff
      real(rp), dimension(:), allocatable :: thS, phS
      type(ShellGrid_T), save :: mixS          ! built once; the REMIX south grid is fixed
      logical,           save :: mixS_init = .false.
      type(ShellGridVar_T) :: mixSGV, tmpN, tmpS

      ! Init a ShellGridVar w.r.t. remix's ShellGrid to map its variables onto
      call initShellVar(mixApp%ion(NORTH)%shGr, SHGR_CORNER, mixSGV)

      ! Init temp arrays to map to per hemisphere
      call initShellVar(voltGrid, SHGR_CORNER, tmpN)
      call initShellVar(voltGrid, SHGR_CORNER, tmpS)

      varOut%data = 0.0
      varOut%mask = .false.
      associate(rmHemi=>mixApp%ion(NORTH), shGr=>mixApp%ion(NORTH)%shGr) !, Nt=>mixApp%ion(NORTH)%shGr%Nt, Np=>mixApp%ion(NORTH)%shGr%Np)       
         mixSGV%data(:,shGr%js:shGr%je) = transpose(rmHemi%St%Vars(:,:,mixVarID))
         mixSGV%data(:,shGr%je+1) = mixSGV%data(:,shGr%js)
         mixSGV%mask = .true.
         call InterpShellVar_TSC_SG(rmHemi%shGr, mixSGV, voltGrid, tmpN)
         
         ! Hacky version for now. Not needed if mix's sg is child of voltron's
         iLat = voltGrid%is
         do while (voltGrid%th(iLat+1) < rmHemi%shGr%maxTheta)
            iLat = iLat + 1
         enddo
         varOut%mask(voltGrid%is:iLat,:) = .true.

      end associate


      associate(rmHemi=>mixApp%ion(SOUTH), shGr=>mixApp%ion(SOUTH)%shGr)       
         ! Set up a proper ShellGrid for Southern hemisphere in SM coordinates
         ! Build the SM-frame south grid once (the REMIX south grid is fixed).
         if (.not. mixS_init) then
            allocate(thS(shGr%Nt+1))
            allocate(phS(shGr%Np+1))
            thS = PI - shGr%th(shGr%ie+1:shGr%is:-1)  ! Flip direction so we go from eq to pole in memory
            phS = shGr%ph(shGr%js:shGr%je+1)  ! Positions are same values, but handedness is flipped. Handle on data mashing
            call GenShellGrid(mixS, thS,phS,"REMIX_SOUTH",nGhosts=(/0,0,0,0/),radO=shGr%radius)
            deallocate(thS, phS)
            mixS_init = .true.
         endif
         call initShellVar(mixS, SHGR_CORNER, mixSGV)  ! blow up old mixSGV and use new southern grid (even thouth they should be the same)
         mixSGV%mask = .true.
         ! Now map Southern hemi mix var onto a SM shellgridvar
         mixSGV%data(:,mixS%js:mixS%je) = transpose(rmHemi%St%Vars(1:shGr%Np,:,mixvarID))
         mixSGV%data(:,mixS%js:mixS%je) = mixSGV%data(mixS%ie+1:mixS%is:-1,mixS%je:mixS%js:-1)
         mixSGV%data(:,mixS%je+1) = mixSGV%data(:,mixS%js)
         ! Ready to interp onto the voltron grid
         call InterpShellVar_TSC_SG(mixS, mixSGV, voltGrid, tmpS)
         
         ! Hacky version for now. Not needed if mix's sg is child of voltron's
         iLat = voltGrid%ie
         do while (voltGrid%th(iLat-1) > shGr%minTheta)
            iLat = iLat - 1
            if (ilat == 1) exit ! ** FIX THIS
         enddo
         varOut%mask(iLat:voltGrid%ie+1,:) = .true.

      end associate

      ! Add both hemis together to get final global potential
      varOut%data = tmpN%data + tmpS%data

   end subroutine mixVarToVoltron

end module mixinterfaceutils

