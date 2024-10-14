! Miscellaneous routines used by the gamera<->remix coupling modules

module mixinterfaceutils
    use mixdefs

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

end module mixinterfaceutils

