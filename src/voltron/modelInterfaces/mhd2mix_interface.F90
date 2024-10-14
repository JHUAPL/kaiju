!Various routines to handle coupling between Gamera magnetosphere and remix ionosphere
!General coupling consists of XXX

module mhd2mix_interface
    use kdefs
    use gamtypes
    use gamutils
    use volttypes
    use cmiutils
    use msphutils, only : RadIonosphere
    use gamapp
    use mixgeom
    use mixinterfaceutils
    
    implicit none

    ! how many variables are we sending (should be consistent with the enumerator in mixdefs.F90)
    integer, parameter :: mhd2mix_varn = 3
    real(rp), private :: Rion

    contains

    subroutine init_mhd2Mix(mhd2mix, gameraApp, remixApp)
        type(mhd2Mix_T), intent(inout) :: mhd2mix
        type(gamApp_T), intent(in) :: gameraApp
        type(mixApp_T), intent(inout) :: remixApp

        integer :: h,l,i,j,k,iG
        real(rp) :: xc,yc,zc, mhd_Rin
        real(rp), allocatable, dimension(:,:,:,:,:) :: mhdJGrid
        real(rp), allocatable, dimension(:,:) :: mhdt, mhdp, mhdtFpd, mhdpFpd
        type(mixGrid_T) :: mhdGfpd
        type(Map_T) :: Map

        Rion = RadIonosphere()
        
        ! allocate remix arrays (adding allocation checks as some of these are now coming from restart file)
        if (.not. allocated(mhd2Mix%gJ)) then
            allocate(mhd2Mix%gJ   (1:JpSh, gameraApp%Grid%js:gameraApp%Grid%je, gameraApp%Grid%ks:GameraApp%Grid%ke, 1:NDIM))
        endif
        
        if (.not. allocated(mhd2Mix%gBAvg)) then
            allocate(mhd2Mix%gBAvg(1:JpSh, gameraApp%Grid%js:gameraApp%Grid%je, gameraApp%Grid%ks:GameraApp%Grid%ke, 1:NDIM))
        endif
        
        if (.not. allocated(mhdJGrid)) then
            allocate(mhdJGrid(1:JpSh, gameraApp%Grid%js:gameraApp%Grid%je, gameraApp%Grid%ks:gameraApp%Grid%ke/2, 1:3, 1:2))
        endif
        
        if (.not. allocated(mhd2Mix%mixInput)) then
            allocate(mhd2Mix%mixInput(1:JpSh, gameraApp%Grid%js:gameraApp%Grid%je, gameraApp%Grid%ks:gameraApp%Grid%ke/2, 1:mhd2mix_varn, 1:2))
        endif
        
        if (.not. allocated(mhd2Mix%JMaps)) then
            allocate(mhd2Mix%JMaps(JpSh,size(remixApp%ion)))
        endif
        
        ! get those grid coordinates (cell centers for Jp)
        do k=gameraApp%Grid%ks,gameraApp%Grid%ke
            do j=gameraApp%Grid%js,gameraApp%Grid%je
                do i=1,JpSh
                    iG = JpSt+i-1
                    call cellCenter(gameraApp%Grid,iG,j,k,xc,yc,zc)
                    mhd2Mix%gBAvg(i,j,k,:) = gameraApp%Grid%B0(iG,j,k,:)
                    ! note conversion to Rion units which are expected on the remix
                    ! side
                    if (k<=gameraApp%Grid%ke/2) then
                        mhdJGrid(i,j,k,:,NORTH) = [xc,yc,zc]/Rion
                    else
                        mhdJGrid(i,j,k-gameraApp%Grid%ke/2,:,SOUTH) = [xc,yc,zc]/Rion
                    endif
                enddo
            enddo
        enddo

        do h=1,size(remixApp%ion)
           ! set up interpolation map(s) for mhd2mix
           do l=1,JpSh
               call mix_mhd_grid(mhdJGrid(l,:,:,:,h),mhdt,mhdp,mhdtFpd,mhdpFpd,mhd_Rin)
               call init_grid_fromTP(mhdGfpd,mhdtFpd,mhdpFpd,isSolverGrid=.false.,isPeriodic=.false.)
               call flip_grid(remixApp%ion(h)%G,remixApp%ion(h)%mixGfpd,mhd_Rin) ! storing flipped
               ! grid for MIX only to
               ! use in mhd2mix
               ! below for zeroing
               ! out equatorward of
               ! MHD boundary, if
               ! necessary
               call mix_set_map(mhdGfpd,remixApp%ion(h)%mixGfpd,Map)
               mhd2mix%JMaps(l,h) = Map
           end do
        end do
        
        deallocate(mhdJGrid)

    end subroutine init_mhd2Mix

    subroutine convertGameraToRemix(mhd2mix, gameraApp, remixApp)
        type(mhd2Mix_T), intent(inout) :: mhd2mix
        type(mixApp_T), intent(inout) :: remixApp
        class(gamApp_T), intent(in) :: gameraApp

        ! convert incoming gamera data to the "remixInputs" variable
        real(rp) ::  Bi2m
        integer  :: i,j,k,iG
        real(rp) :: xc,yc,zc

        real(rp) :: Con(NVAR)
        real(rp) :: Cs,gB0,gv0,gx0,J11,B0(NDIM)

        gB0 = gameraApp%Model%Units%gB0
        gv0 = gameraApp%Model%Units%gv0
        gx0 = gameraApp%Model%Units%gx0

        !Only working on Bxyz from perturbation
        !B0 in inner region (where we care for remix)
        !Assumd to be current-free
        call GetShellJ(gameraApp%Model, gameraApp%Grid, gameraApp%State%Bxyz, mhd2mix%gJ)

        !Now loop over only shells and populate mhdvars to be sent to mix
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,iG,j,k,Bi2m,xc,yc,zc,Con,J11,B0,Cs)
        do k=gameraApp%Grid%ks,gameraApp%Grid%ke
            do j=gameraApp%Grid%js,gameraApp%Grid%je
                do i=1,JpSh
                    iG = JpSt+i-1

                    !Get new bAvg
                    B0 = gameraApp%State%Bxyz(iG,j,k,:) + gameraApp%Grid%B0(iG,j,k,:)
                    mhd2mix%gBAvg(i,j,k,:) = mhd2mix%wAvg*B0 + (1-mhd2mix%wAvg)*mhd2mix%gBAvg(i,j,k,:)

                    !Get FAC
                    B0 = mhd2mix%gBAvg(i,j,k,:)
                    !B0 = gameraApp%Grid%B0(iG,j,k,:)
                    J11 = dot_product(mhd2mix%gJ(i,j,k,:),normVec(B0))

                    ! get cell centers
                    call cellCenter(gameraApp%Grid,iG,j,k,xc,yc,zc)
                    !Get ion2mag scaling factor for field
                    Bi2m = BIon2Mag(xc,yc,zc)

                    if (k<=gameraApp%Grid%ke/2) then
                        !!! NOTE: assuming gB0 in nT and gx0 in m and gv0 in m/s

                        ! note conversion to microA/m^2
                        mhd2mix%mixInput(i,j,k,MHDJ,NORTH) = J11*Bi2m*(gB0/gx0*1.e4/4/PI)
                        ! note conversion to g/cm^3
                        mhd2mix%mixInput(i,j,k,MHDD,NORTH) = gameraApp%State%Gas(iG,j,k,DEN,BLK)*(gB0*1.e-7/gv0)**2/4/pi
                        ! get sound speed first
                        Con = gameraApp%State%Gas(iG,j,k,:,BLK)
                        call CellPress2Cs(gameraApp%Model,Con,Cs)
                        mhd2mix%mixInput(i,j,k,MHDC,NORTH) = Cs*gv0*1.e2
                    else
                        mhd2mix%mixInput(i,j,k-gameraApp%Grid%ke/2,MHDJ,SOUTH) = J11*Bi2m*(gB0/gx0*1.e4/4/PI)
                        mhd2mix%mixInput(i,j,k-gameraApp%Grid%ke/2,MHDD,SOUTH) = gameraApp%State%Gas(iG,j,k,DEN,BLK)*(gB0*1.e-7/gv0)**2/4/pi
                        ! get sound speed first
                        Con = gameraApp%State%Gas(iG,j,k,:,BLK)
                        call CellPress2Cs(gameraApp%Model,Con,Cs)
                        mhd2mix%mixInput(i,j,k-gameraApp%Grid%ke/2,MHDC,SOUTH) = Cs*gv0*1.e2
                    endif
                enddo
            enddo
        enddo

    end subroutine convertGameraToRemix

  ! assume what's coming here is mhdvars(i,j,k,var,hemisphere)
  ! thus the transposes below
  subroutine mapGameraToRemix(mhd2mix, remixApp)
    type(mhd2Mix_T), intent(inout) :: mhd2mix
    type(mixApp_T), intent(inout) :: remixApp

    real(rp), dimension(:,:), allocatable :: F
    integer :: l,h ! hemisphere
    integer :: v ! mhd var
    integer :: Nth
    integer :: i,j
    real(rp) :: mixImin

    if (size(mhd2mix%mixInput,5).ne.size(remixApp%ion)) then
       write(*,*) 'The number of hemispheres in mhdvars is different from the size of the MIX ionosphere object. I am stopping.'
       stop
    end if

    do h=1,size(remixApp%ion)
       Nth = remixApp%ion(h)%G%Nt !Theta cells
       do v=1,size(mhd2mix%mixInput,4)
          do l=1,JpSh ! here we loop over Jshells but always use the last one (F)
             ! note the transpose to conform to the MIX layout (phi,theta)
             call mix_map_grids(mhd2mix%JMaps(l,h),transpose(mhd2mix%mixInput(l,:,:,v,h)),F)

             ! note, cleaning MHD vars equatorward of the MHD boundary
             ! if the MIX boundary is equatorward of MHD boundary
             select case (v)
             case (MHDJ)
                ! zero out the current
!                where (remixApp%ion(h)%mixGfpd%mask.eq.-1) F=0._rp
                do i=1,remixApp%ion(h)%mixGfpd%Np
                   do j=1,remixApp%ion(h)%mixGfpd%Nt
                      if(remixApp%ion(h)%mixGfpd%mask(i,j) .eq. -1) F(i,j) = 0._rp
                   enddo
                enddo
             case (MHDD, MHDC)
                ! set density and sound speed to min values
!                where (remixApp%ion(h)%mixGfpd%mask.eq.-1) F=minval(mhd2mix%mixInput(l,:,:,v,h))
                mixImin = minval(mhd2mix%mixInput(l,:,:,v,h))
                do i=1,remixApp%ion(h)%mixGfpd%Np
                   do j=1,remixApp%ion(h)%mixGfpd%Nt
                      if(remixApp%ion(h)%mixGfpd%mask(i,j) .eq. -1) F(i,j) = mixImin
                   enddo
                enddo
             end select
          end do

          select case (v)
          case (MHDJ)
             remixApp%ion(h)%St%Vars(:,:,FAC) = F
             !This shouldn't matter b/c it's already zero, but if it's not oh man that's bad
             !Probably if it's not this should return an error and the user should rethink what they're doing
             !But i've come to realize that the user probably doesn't want the warning, they just want a mommy who can make problems magically go away
             remixApp%ion(h)%St%Vars(:,Nth-1:Nth,FAC) = 0.0
          case (MHDD)
             remixApp%ion(h)%St%Vars(:,:,DENSITY) = F
          case (MHDC)
             remixApp%ion(h)%St%Vars(:,:,SOUND_SPEED) = F
          end select
       end do
    end do

  end subroutine mapGameraToRemix


    !Calculate Ion->Mag scaling factor, xyz are in units of Re (like Gamera)
    ! note for Bi2m to be correct, xc,yc,zc,r need to be in
    ! units of Rion where we solve for the potential
    ! Rion is defined in msphutils in code units

    function BIon2Mag(x,y,z) result(Bi2m)
        real(rp), intent(in) :: x,y,z
        real(rp) :: Bi2m

        real(rp) :: rMag,Rp,zor2,mZ,pZ

        rMag = norm2([x,y,z])
        Rp = rMag/Rion
        zor2 = (z/rMag)**2.0
        mZ = (1.0-zor2)/Rp
        pZ = 1.0+3.0*zor2

        Bi2m = (Rp**3.0)*sqrt( 1.0 + 3.0*(1.0-mZ) )/sqrt(pZ)

        ! bion2bmag in LFM-para/src/interfaces/MHDBoundaryInterface.C
        !Bi2m  = r**3*sqrt(1+3*(1-(1-(zc/r)**2)/r))/sqrt(1+3*(zc/r)**2) 


    end function BIon2Mag

    function need_remix_BC(Model,Grid)
        type(Model_T) ,intent(in) :: Model
        type(Grid_T)  ,intent(in) :: Grid
        logical :: need_remix_BC
        
        need_remix_BC = Grid%ijkShift(IDIR) == 0

    end function need_remix_BC

end module mhd2mix_interface
