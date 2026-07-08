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

   ! Author: Anthony M. Sciola (origin/backport/namelessking-fixes; clone of imag2mix's
   ! superSampleMomentsPrecip, NOT used by imag2mix -- packInputs is the only caller).
   ! git blame shows the splice commit 4a2169405 (cookieenick), not the author: use
   ! `git log -S "subroutine voltVarToMix"` or blame the backport branch. Only local
   ! change: nGood (not nPnts**2) in superSampleVoltToMix.
   subroutine mapShellVarToMix(srcGrid, srcVar, ion, mixVarID)
        type(ShellGrid_T), intent(in) :: srcGrid
        type(ShellGridVar_T), intent(in) :: srcVar
        type(mixIon_T), dimension(:), intent(inout) :: ion
        integer, intent(in) :: mixVarID

        integer :: i, j
        real(rp) :: varOut
        logical :: isTasty
        integer :: nPnts_sample

        nPnts_sample = 3

        ion(NORTH)%St%Vars(:,:,mixVarID) = 0.0_rp
        ion(SOUTH)%St%Vars(:,:,mixVarID) = 0.0_rp

        associate(GN => ion(NORTH)%G, GS => ion(SOUTH)%G)
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,varOut,isTasty)
            do j = 1, GN%Np
                do i = 1, GN%Nt
                    varOut = 0.0_rp
                    call superSampleShellToMix(srcGrid, srcVar, GN, i, j, nPnts_sample, varOut, isTasty)
                    if (isTasty) ion(NORTH)%St%Vars(j,i,mixVarID) = varOut

                    varOut = 0.0_rp
                    call superSampleShellToMix(srcGrid, srcVar, GS, i, j, nPnts_sample, varOut, isTasty, doSouthO=.true.)
                    if (isTasty) ion(SOUTH)%St%Vars(j,i,mixVarID) = varOut
                end do
            end do
        end associate

    end subroutine mapShellVarToMix

    ! Nearest-neighbour analogue of mapShellVarToMix, for categorical / integer-code fields
    subroutine mapShellVarToMixNN(srcGrid, srcVar, ion, mixVarID)
        type(ShellGrid_T),    intent(in)    :: srcGrid
        type(ShellGridVar_T), intent(in)    :: srcVar
        type(mixIon_T), dimension(:), intent(inout) :: ion
        integer,              intent(in)    :: mixVarID

        integer  :: i, j, jM
        real(rp) :: varOut, thIn, phIn
        logical  :: isTasty

        ion(NORTH)%St%Vars(:,:,mixVarID) = 0.0_rp
        ion(SOUTH)%St%Vars(:,:,mixVarID) = 0.0_rp

        associate(GN => ion(NORTH)%G, GS => ion(SOUTH)%G)
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,jM,varOut,isTasty,thIn,phIn)
            do j = 1, GN%Np
                do i = 1, GN%Nt
                    ! North (REMIX cell center), sampled directly on the shell
                    thIn = GN%t(j,i); phIn = GN%p(j,i)
                    call InterpShellVar_NN_pnt(srcGrid, srcVar, thIn, phIn, varOut, goodInterpO=isTasty)
                    if (isTasty) ion(NORTH)%St%Vars(j,i,mixVarID) = varOut

                    ! South, reflect theta about the equator
                    jM = GS%Np + 1 - j
                    thIn = PI - GS%t(jM,i); phIn = GS%p(jM,i)
                    call InterpShellVar_NN_pnt(srcGrid, srcVar, thIn, phIn, varOut, goodInterpO=isTasty)
                    if (isTasty) ion(SOUTH)%St%Vars(j,i,mixVarID) = varOut
                end do
            end do
        end associate

    end subroutine mapShellVarToMixNN

    subroutine superSampleShellToMix(srcGrid, srcVar, G, i, j, nPnts, varOut, isTasty, doSouthO)
        type(ShellGrid_T), intent(in) :: srcGrid
        type(ShellGridVar_T), intent(in) :: srcVar
        type(mixGrid_T), intent(in) :: G
        integer, intent(in) :: i, j
        integer, intent(in) :: nPnts
        real(rp), intent(inout) :: varOut
        logical, intent(inout) :: isTasty
        logical, optional, intent(in) :: doSouthO

        integer :: iPnt, jPnt, nGood, jM, imap
        logical :: doSouth, isTastyPnt
        real(rp) :: thetaIn, phiIn
        real(rp), dimension(2) :: thetaBnds, phiBnds
        real(rp), dimension(nPnts) :: theta_subdiv, phi_subdiv
        real(rp) :: varPnt

        varOut = 0.0_rp
        nGood = 0
        if (nPnts < 2) then
            write(*,*) "ERROR in remixCplTypesSub.F90:superSampleShellToMix"
            write(*,*) "nPnts must be at least 2"
            stop
        end if

        if (present(doSouthO)) then
            doSouth = doSouthO
        else
            doSouth = .false.
        end if

        if (.not. doSouth) then
            thetaIn = G%t(j,i)
            phiIn = G%p(j,i)
            if (i == 1) then
                thetaBnds(1) = thetaIn
                thetaBnds(2) = thetaIn + 0.5_rp*G%dt(j,i)
            else if (i == G%Nt) then
                thetaBnds(1) = thetaIn - 0.5_rp*G%dt(j,i-1)
                thetaBnds(2) = thetaIn + 0.5_rp*G%dt(j,i-1)
            else
                thetaBnds(1) = thetaIn - 0.5_rp*G%dt(j,i-1)
                thetaBnds(2) = thetaIn + 0.5_rp*G%dt(j,i)
            end if

            if (j == 1 .or. j == G%Np) then
                phiBnds(1) = phiIn - 0.5_rp*G%dp(G%Np-1,i)
                phiBnds(2) = phiIn + 0.5_rp*G%dp(1,i)
            else
                phiBnds(1) = phiIn - 0.5_rp*G%dp(j-1,i)
                phiBnds(2) = phiIn + 0.5_rp*G%dp(j,i)
            end if
        else
            jM = G%Np + 1 - j
            thetaIn = PI - G%t(jM,i)
            phiIn = G%p(jM,i)
            if (i == 1) then
                thetaBnds(1) = thetaIn - 0.5_rp*G%dt(jM,i)
                thetaBnds(2) = thetaIn
            else if (i == G%Nt) then
                thetaBnds(1) = thetaIn - 0.5_rp*G%dt(jM,i-1)
                thetaBnds(2) = thetaIn + 0.5_rp*G%dt(jM,i-1)
            else
                thetaBnds(1) = thetaIn - 0.5_rp*G%dt(jM,i)
                thetaBnds(2) = thetaIn + 0.5_rp*G%dt(jM,i-1)
            end if

            if (j == 1 .or. j == G%Np) then
                phiBnds(1) = phiIn - 0.5_rp*G%dp(1,i)
                phiBnds(2) = phiIn + 0.5_rp*G%dp(G%Np-1,i)
            else
                phiBnds(1) = phiIn - 0.5_rp*G%dp(jM-1,i)
                phiBnds(2) = phiIn + 0.5_rp*G%dp(jM,i)
            end if
        end if

        call genMixGridSubdivs(thetaBnds, phiBnds, nPnts, theta_subdiv, phi_subdiv)

        ! The destination is Remix's SM grid
        call get_coord_mapping(srcGrid%coord, "SM", imap)

        do iPnt = 1, nPnts
            do jPnt = 1, nPnts
                call InterpShellVar_TSC_pnt(srcGrid, srcVar, theta_subdiv(iPnt), phi_subdiv(jPnt), varPnt, goodInterpO=isTastyPnt, imap_in=imap)
                 if (.not. isTastyPnt) then
               cycle
            endif
            ! Otherwise we got good data
            nGood = nGood + 1
            varOut = varOut + varPnt
         enddo
      enddo
      ! Finish average and make final isTasty decision
      if (nGood > 0) then
         isTasty = .true.

         ! If dividing by nPnts**2, we are including zero-vals from bad points
         ! If dividing by nGood, we might be letting just a few points represent a whole cell
         ! Since we are getting stuff from voltron grid, all points should be well-defined over the remix domain, so there shouldn't be a difference
         varOut = varOut / nGood
      else
         isTasty = .false.
      endif
    end subroutine superSampleShellToMix

    subroutine genMixGridSubdivs(thetaBnds, phiBnds, nPnts, theta_subdiv, phi_subdiv)
        real(rp), dimension(2), intent(in) :: thetaBnds
        real(rp), dimension(2), intent(in) :: phiBnds
        integer, intent(in) :: nPnts
        real(rp), dimension(nPnts), intent(out) :: theta_subdiv
        real(rp), dimension(nPnts), intent(out) :: phi_subdiv

        integer :: n
        real(rp) :: dCos, arg

        dCos = cos(thetaBnds(1)) - cos(thetaBnds(2))

        theta_subdiv(1) = min(max(thetaBnds(1), 0.0_rp), PI)
        theta_subdiv(nPnts) = min(max(thetaBnds(2), 0.0_rp), PI)

        do n = 2, nPnts-1
            arg = cos(thetaBnds(1)) - dCos*real(n-1,rp)/real(nPnts-1,rp)
            arg = min(max(arg, -1.0_rp), 1.0_rp)
            theta_subdiv(n) = min(max(acos(arg), 0.0_rp), PI)
        end do

        phi_subdiv = [(modulo(phiBnds(1)+(n-1)*(phiBnds(2)-phiBnds(1))/real(nPnts-1,rp), 2*PI), n=1,nPnts)]
    end subroutine genMixGridSubdivs

end module mixinterfaceutils
