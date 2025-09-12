! module to convert imag data to mhd data
module imag2mhd_interface
    use volttypes
    use cmiutils
    use planethelper
    use gridutils
    use shellInterp

    implicit none

    contains

    ! have this short interface function in case we want to alter how we fill source data, or use something other than chimp and imagApp
    subroutine CoupleSourceToMhd(vApp)
        class(voltApp_T), intent(inout) :: vApp

        real(rp), dimension(:,:,:)  , allocatable :: gPsi
        real(rp), dimension(:,:,:,:), allocatable :: Eijk
        logical , dimension(:,:,:)  , allocatable :: isGoodCC

        integer :: n,i,j,k,Nk
        real(rp) :: x1,x2,x1c,di,dj,dk,Rion,t,eScl
        real(rp), dimension(NDIM) :: xyz
        real(rp) :: Qs(8) !Squish all the corners of a single cell
        real(rp) :: imW(IM_D_RING:IM_TSCL) !All imag variables
        logical :: isTasty,isG


        !Add projection to ion instead of invlat
        Rion = vApp%shGrid%radius !Voltron grid radius in Rx
        eScl = vApp%gApp%rm2g !Scaling for remix kV

        associate(Gr=>vApp%gApp%Grid,gModel=>vApp%gApp%Model)
            Gr%Gas0 = 0.0

        !Clean up axis projections
        !---
            !Force agreement at k-seam
            vApp%chmp2mhd%xyzSquish(:,:,Gr%ke+1,1:2) = vApp%chmp2mhd%xyzSquish(:,:,Gr%ks,1:2)
            vApp%chmp2mhd%isGood   (:,:,Gr%ke+1)     = vApp%chmp2mhd%isGood   (:,:,Gr%ks)

            !Clean up axis
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,k)
            do i=Gr%is,Gr%ie+1
                !If any nodes around the ring are bad then poison all of them
                !+X pole
                if (.not. all(vApp%chmp2mhd%isGood(i,Gr%js,Gr%ks:Gr%ke+1))) then
                    vApp%chmp2mhd%isGood(i,Gr%js,:) = .false.
                else
                    !All are good, define all to be average
                    vApp%chmp2mhd%xyzSquish(i,Gr%js,:,1) = ArithMean(vApp%chmp2mhd%xyzSquish(i,Gr%js,Gr%ks:Gr%ke,1))
                    vApp%chmp2mhd%xyzSquish(i,Gr%js,:,2) =  CircMean(vApp%chmp2mhd%xyzSquish(i,Gr%js,Gr%ks:Gr%ke,2))
                endif

                !-X pole
                if (.not. all(vApp%chmp2mhd%isGood(i,Gr%je+1,Gr%ks:Gr%ke+1))) then
                    vApp%chmp2mhd%isGood(i,Gr%je+1,:) = .false.
                else
                    !All are good, define all to be average
                    vApp%chmp2mhd%xyzSquish(i,Gr%je+1,:,1) = ArithMean(vApp%chmp2mhd%xyzSquish(i,Gr%je+1,Gr%ks:Gr%ke,1))
                    vApp%chmp2mhd%xyzSquish(i,Gr%je+1,:,2) =  CircMean(vApp%chmp2mhd%xyzSquish(i,Gr%je+1,Gr%ks:Gr%ke,2))
                endif
            enddo

        !Get cell-centered projections
        !---
            allocate(isGoodCC(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg))
            isGoodCC = .false.
            

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,xyz,x1,x2,Qs)
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%isg,Gr%ie
                        
                        if (i<Gr%is) then
                            !Below inner boundary, do dipole projection
                            isGoodCC(i,j,k) = .true.
                            xyz = Gr%xyzcc(i,j,k,:) !Gamera grid center
                            call Proj2Rad(xyz,Rion,x1,x2)
                            Gr%Gas0(i,j,k,PROJLAT) = x1
                            Gr%Gas0(i,j,k,PROJLON) = x2

                        else
                            !Get value from xyzsquish

                            if ( all(vApp%chmp2mhd%isGood(i:i+1,j:j+1,k:k+1)) ) then
                                !All values are good, so just do this thing
                                call SquishCorners(vApp%chmp2mhd%xyzSquish(i:i+1,j:j+1,k:k+1,1),Qs)
                                Gr%Gas0(i,j,k,PROJLAT) = ArithMean(Qs)
                                call SquishCorners(vApp%chmp2mhd%xyzSquish(i:i+1,j:j+1,k:k+1,2),Qs)
                                Gr%Gas0(i,j,k,PROJLON) = CircMean(Qs)
                                isGoodCC(i,j,k) = .true.
                            else
                                Gr%Gas0(i,j,k,PROJLAT) = 0.0
                                Gr%Gas0(i,j,k,PROJLON) = 0.0
                                isGoodCC(i,j,k) = .false.
                            endif !corner projection => center
                        endif !inner-i vs. not

                    enddo
                enddo !j
            enddo !k

        !Get convective E field for ingestion
        !---

            call FillGhostsCC(gModel,Gr,Gr%Gas0(:,:,:,PROJLAT))
            call FillGhostsCC(gModel,Gr,Gr%Gas0(:,:,:,PROJLON))

            allocate(gPsi(Gr%isg:Gr%ieg+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1)) !All nodes
            gPsi = 0.0

            !Loop over "real" nodes and get total potential from voltron grid
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,x1,x2,x1c,xyz,isG)
            do k=Gr%ks,Gr%ke+1
                do j=Gr%js,Gr%je+1
                    do i=Gr%isg,Gr%ie+1

                        if (i < Gr%is) then
                            !Use dipole projection
                            xyz = Gr%xyz(i,j,k,:) !Gamera grid corner
                            call Proj2Rad(xyz,Rion,x1,x2)
                            isG = .true.
                        else
                            x1  = vApp%chmp2mhd%xyzSquish(i,j,k,1)
                            x2  = vApp%chmp2mhd%xyzSquish(i,j,k,2)
                            isG = vApp%chmp2mhd%isGood(i,j,k)
                        endif
                        x1c = PI/2 - x1
                        if (isG) then
                            !Get potential in kV
                            call InterpShellVar_TSC_pnt(vApp%shGrid,vApp%State%potential_total,&
                                                        x1c,x2,gPsi(i,j,k) )
                        endif
                    enddo
                enddo !j
            enddo !k

            !Create edge (ijk) E and then convert to cell-centered Exyz
            allocate(Eijk(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM))
            Eijk = 0.0

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,di,dj,dk)
            do k=Gr%ks,Gr%ke+1
                do j=Gr%js,Gr%je+1
                    do i=Gr%isg,Gr%ie+1
                        di = Gr%edge(i,j,k,IDIR)
                        Eijk(i,j,k,IDIR) = -( gPsi(i+1,j,k) - gPsi(i,j,k) )/di

                        dj = Gr%edge(i,j,k,JDIR)
                        Eijk(i,j,k,JDIR) = -( gPsi(i,j+1,k) - gPsi(i,j,k) )/dj

                        dk = Gr%edge(i,j,k,KDIR)
                        Eijk(i,j,k,KDIR) = -( gPsi(i,j,k+1) - gPsi(i,j,k) )/dk

                        Eijk(i,j,k,:) = Eijk(i,j,k,:)/eScl
                        
                    enddo
                enddo
            enddo !k

            call FixEFieldLFM(gModel,Gr,Eijk) !Force agreement at k-seam

            !Convert edge E fields to xyz and store
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k)
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%isg,Gr%ieg
                        if (isGoodCC(i,j,k)) then
                            Gr%Gas0(i,j,k,IONEX:IONEZ) = CellExyz(gModel,Gr,Eijk,i,j,k)
                        endif
                    enddo
                enddo !j
            enddo !k

            do n=IONEX,IONEZ
                call FillGhostsCC(gModel,Gr,Gr%Gas0(:,:,:,n))
            enddo


            !t = (vApp%gApp%Model%t)*(vApp%gApp%Model%Units%gT0)

            !Loop over and get imag data
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,x1c,x2,imW,isTasty)
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%isg,Gr%ie
                        Gr%Gas0(i,j,k,IM_D_RING:IM_P_RING) = 0.0
                        Gr%Gas0(i,j,k,IM_D_COLD:IM_P_COLD) = 0.0
                        Gr%Gas0(i,j,k,IM_TSCL) = 0.0

                        if (isGoodCC(i,j,k)) then
                            !Get imag value here
                            imW = 0.0
                            x1c = PI/2 - Gr%Gas0(i,j,k,PROJLAT)
                            x2  = Gr%Gas0(i,j,k,PROJLON)

                            !call vApp%imagApp%getMoments(x1c,x2,t,imW,isTasty)
                            call vApp%imagApp%getMoments(x1c,x2,vApp%time,imW,isTasty)
                            if (isTasty) then
                                !Density/pressure coming back in #/cc and nPa
                                !Ingestion timescale coming back in seconds
                                Gr%Gas0(i,j,k,IM_D_RING) = imW(IM_D_RING)
                                Gr%Gas0(i,j,k,IM_D_COLD) = imW(IM_D_COLD)
                                Gr%Gas0(i,j,k,IM_P_RING) = imW(IM_P_RING)/vApp%gApp%Model%Units%gP0
                                Gr%Gas0(i,j,k,IM_P_COLD) = imW(IM_P_COLD)/vApp%gApp%Model%Units%gP0
                                Gr%Gas0(i,j,k,IM_TSCL  ) = imW(IM_TSCL  )/vApp%gApp%Model%Units%gT0
                            endif
                            
                        endif
                    enddo !i
                enddo
            enddo !k

            !Do averaging for first cell next to singularity
            Nk = Gr%ke-Gr%ks+1

            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i)
            do i=Gr%is,Gr%ie
                !+X pole
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IM_D_RING),Nk)
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IM_P_RING),Nk)
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IM_D_COLD),Nk)
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IM_P_COLD),Nk)
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IM_TSCL  ),Nk)

                !-X pole
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IM_D_RING),Nk)
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IM_P_RING),Nk)
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IM_D_COLD),Nk)
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IM_P_COLD),Nk)
                call ReAverage(gModel,Gr,Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IM_TSCL  ),Nk)
                
            enddo

            do n=IM_D_RING,IM_TSCL
                call FillGhostsCC(gModel,Gr,Gr%Gas0(:,:,:,n))
            enddo
        end associate

        !---
        contains

            !Take Nk ring of cell-centered data and average
            subroutine ReAverage(Model,Gr,GasK,Nk)
                type(Model_T), intent(in)      :: Model
                type(Grid_T) , intent(in)      :: Gr
                integer      , intent(in)      :: Nk
                real(rp)     , intent(inout)   :: GasK(Nk)

                real(rp) :: Q
                Q = AvgOverGood(GasK,Nk)
                GasK(:) = Q
            end subroutine 

            function AvgOverGood(Q,Nk) result(Qavg)
                real(rp), intent(in), dimension(Nk) :: Q
                integer , intent(in) :: Nk

                real(rp) :: Qavg
                integer :: Nkg

                if ( all(Q>TINY) ) then
                    Nkg = count(Q>TINY)
                    Qavg = sum(Q,mask=(Q>TINY))/Nkg
                else
                    Qavg = 0.0
                endif

            end function AvgOverGood

            !Take a full-sized cell-centered grid w/ active values defined and fill ghosts
            subroutine FillGhostsCC(Model,Gr,Q)
                type(Model_T), intent(in) :: Model
                type(Grid_T) , intent(in) :: Gr
                real(rp), intent(inout) :: Q(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg)

                integer :: i,j,k,ip,jp,kp
                logical :: isActive

                !!$OMP PARALLEL DO default(shared) collapse(2) &
                !!$OMP private(i,j,k,isActive,ip,jp,kp)
                ! this causes a race condition copying values between ghost cells
                ! probably a false positive since some of the cells are just copying
                ! values onto themselves, but easier to remove for now
                do k=Gr%ksg,Gr%keg
                    do j=Gr%jsg,Gr%jeg
                        do i=Gr%isg,Gr%ieg
                            isActive = (j >= Gr%js) .and. (j <= Gr%je) .and. &
                                       (k >= Gr%ks) .and. (k <= Gr%ks)
                            if(.not. isActive) then
                                !If still here map this ghost to active and set value based on active
                                call lfmIJKcc(Model,Gr,i,j,k,ip,jp,kp)
                                Q(i,j,k) = Q(ip,jp,kp)
                            endif
                        enddo
                    enddo !j
                enddo !k

            end subroutine FillGhostsCC

            !Project xyz along dipole to R0 and return lat (x1) and lon (x2)
            subroutine Proj2Rad(xyz,R0,x1,x2)
                real(rp), intent(in ) :: xyz(NDIM), R0
                real(rp), intent(out) :: x1,x2

                real(rp), dimension(NDIM) :: xyz0

                xyz0 = DipoleShift(xyz,R0)
                x1 = asin(xyz0(ZDIR)/R0) !Lat
                x2 = katan2(xyz0(YDIR),xyz0(XDIR)) !katan => [0,2pi] instead of [-pi,pi]

            end subroutine Proj2Rad 

    end subroutine

end module

