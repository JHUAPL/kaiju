! module to convert remix data to mhd data
module mix2mhd_interface
    use volttypes
    use mixgeom, only : init_grid_fromTP
    use msphutils, only : RadIonosphere
    use mixinterfaceutils
    use mixinterp
    use cmiutils
    use uservoltic ! required to have IonInnerBC_T defined


    implicit none

    contains

    ! have this short interface function in case we want to alter how we fill potential data
    subroutine CouplePotentialToMhd(vApp)
        class(voltApp_T), intent(inout) :: vApp

        call mapRemixToGamera(vApp%gApp, vApp%remixApp)
        call convertRemixToGamera(vApp%gApp, vApp%remixApp)

    end subroutine

    subroutine init_mix2MhdCoupler(gameraApp, remixApp)
        class(gamCoupler_T), intent(inout) :: gameraApp
        class(mixApp_T), intent(inout) :: remixApp

        integer :: i,j,k,iG,h,l
        real(rp) :: mhd_Rin
        real(rp), allocatable, dimension(:,:,:,:,:) :: mhdPsiGrid
        real(rp), allocatable, dimension(:,:) :: mhdt, mhdp, mhdtFpd, mhdpFpd
        type(mixGrid_T) :: mhdG
        type(Map_T) :: Map
        real(rp) :: gB0,gv0,gx0
        logical :: isRestart

        gB0 = gameraApp%Model%Units%gB0
        gv0 = gameraApp%Model%Units%gv0
        gx0 = gameraApp%Model%Units%gx0

        gameraApp%rm2g = gB0*gV0*gx0*1.0e-12 !Scaling factor for remix potential [kV]
        gameraApp%Rion = RadIonosphere()

        ! allocate remix arrays
        allocate(gameraApp%gPsi(1:PsiSh+1,gameraApp%Grid%js:gameraApp%Grid%je+1,gameraApp%Grid%ks:gameraApp%Grid%ke+1))
        allocate(mhdPsiGrid(1:PsiSh+1, gameraApp%Grid%js:gameraApp%Grid%je+1, gameraApp%Grid%ks:gameraApp%Grid%ke/2+1, 1:3, 1:2))
        allocate(gameraApp%mixOutput(1:PsiSh+1, gameraApp%Grid%js:gameraApp%Grid%je+1, gameraApp%Grid%ks:gameraApp%Grid%ke/2+1, 1:mix2mhd_varn, 1:2))
        allocate(gameraApp%PsiMaps(PsiSh,size(remixApp%ion)))
        gameraApp%gPsi = 0.0

        ! get those grid coordinates (corner centers for Psi)
        do k=gameraApp%Grid%ks,gameraApp%Grid%ke+1
            do j=gameraApp%Grid%js,gameraApp%Grid%je+1
                ! note, PsiShells give shell numbers based on cell centers per our
                ! convenion
                ! thus, no -1 below
                do i=1,PsiSh+1
                    iG = PsiSt+i-1

                    ! note conversion to Rion units which are expected on the remix
                    ! side
                    if (k<=gameraApp%Grid%ke/2+1) then
                        mhdPsiGrid(i,j,k,:,NORTH) = gameraApp%Grid%xyz(iG,j,k,:)/gameraApp%Rion
                    endif
                    if (k>=gameraApp%Grid%ke/2+1) then
                        mhdPsiGrid(i,j,k-gameraApp%Grid%ke/2,:,SOUTH) = gameraApp%Grid%xyz(iG,j,k,:)/gameraApp%Rion
                    endif
                enddo
            enddo
        enddo

        do h=1,size(remixApp%ion)
            ! set up interpolation map(s) for mix2mhd
            do l=1,PsiSh
                call mix_mhd_grid(mhdPsiGrid(l,:,:,:,h),mhdt,mhdp,mhdtFpd,mhdpFpd,mhd_Rin)
                call init_grid_fromTP(mhdG,mhdt,mhdp,isSolverGrid=.false.)
                call mix_set_map(remixApp%ion(h)%G,mhdG,Map)
                gameraApp%PsiMaps(l,h) = Map
            enddo
        enddo

        !Initialize the mixOutput mapping
        gameraApp%mixOutput = 0.0
        isRestart = gameraApp%Model%isRestart
        if (isRestart) then
            !We have data from mix restart file
            write(*,*) "mapRemixToGamera"
            call mapRemixToGamera(gameraApp, remixApp)
        end if

        deallocate(mhdPsiGrid, mhdt, mhdp, mhdtFpd, mhdpFpd)

    end subroutine

    subroutine mapRemixToGamera(gameraApp, remixApp)
        class(gamCoupler_T), intent(inout) :: gameraApp
        class(mixApp_T), intent(inout) :: remixApp
        integer :: l,h ! hemisphere
        integer :: v ! mhd var

        real(rp), dimension(:,:), allocatable, save :: gPsi_tmp  ! gamera potential

        ! map to potential shells
        do h=1,size(remixApp%ion)
            do v=1,size(gameraApp%mixOutput,4)  ! mirroring mix2mhd here, but for now only one variable
                do l=1,PsiSh
                    call mix_map_grids(gameraApp%PsiMaps(l,h),remixApp%ion(h)%St%Vars(:,:,POT),gPsi_tmp)
                    ! this is going back to gamera
                    select case (v)
                    case (MHDPSI)
                        gameraApp%mixOutput(l,:,:,v,h) = gPsi_tmp   ! north
                    end select
                enddo
            enddo
        enddo
    end subroutine mapRemixToGamera

    subroutine convertRemixToGamera(gameraApp, remixApp, doCorotO)
        class(gamCoupler_T), intent(inout) :: gameraApp
        class(mixApp_T), intent(inout) :: remixApp
        logical, intent(in), optional :: doCorotO

        ! convert the "remixOutputs" variable to inEijk and inExyz, which are in
        ! Gamera coordinates
        integer :: i,nbc
        logical :: doCorot

         ! populate potential on gamera grid
         gameraApp%gPsi = 0.0
         do i=1,PsiSh+1
            gameraApp%gPsi(i,:,gameraApp%Grid%ks:gameraApp%Grid%ke/2+1)   = gameraApp%mixOutput(i,:,:,MHDPSI,NORTH)
            gameraApp%gPsi(i,:,gameraApp%Grid%ke/2+1:gameraApp%Grid%ke+1) = gameraApp%mixOutput(i,:,:,MHDPSI,SOUTH)
         enddo

        if (present(doCorotO)) then
          doCorot = doCorotO
        else
          doCorot = .true.
        endif

        ! add corotation
        if (doCorot) call CorotationPot(gameraApp%Model, gameraApp%Grid, gameraApp%gPsi)

        ! find the remix BC to write data into
        nbc = FindBC(gameraApp%Model,gameraApp%Grid,INI)
        SELECT type(iiBC=>gameraApp%Grid%externalBCs(nbc)%p)
            TYPE IS (IonInnerBC_T)
                call Ion2MHD(gameraApp%Model,gameraApp%Grid,gameraApp%gPsi,iiBC%inEijk,iiBC%inExyz,gameraApp%rm2g)
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in remix IC'
                stop
        END SELECT

    end subroutine convertRemixToGamera

end module

