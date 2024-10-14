!Various routines to map to/calculate on a Grid
module gridutils
    use gamtypes
    use math
    use gamutils
    use quadrature
    use metric
    use planethelper
    use ringutils
    use multifluid
    
    implicit none


    !Routine to perform flux->field calculation
    abstract interface
        function CellBxyz_T(Model,Grid,magFlux,i,j,k) result(Bxyz)
            Import :: rp, Model_T, Grid_T, NDIM
            type(Model_T), intent(in) :: Model
            type(Grid_T) , intent(in) :: Grid
            real(rp)  , intent(in) :: magFlux(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM)
            integer, intent(in) :: i,j,k
            real(rp), dimension(NDIM) :: Bxyz
        end function CellBxyz_T
    end interface
    !This points to chosen routine
    procedure(CellBxyz_T), pointer :: CellBxyz => gamCellBxyz

    contains

    function GetCellDT(Model,Gr,State,i,j,k) result(dtijk)
        integer, intent(in) :: i,j,k
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp) :: dtijk
        call CellDT(Model,Gr,State,i,j,k,dtijk)
    end function GetCellDT

    function GetCellDT_Sonic(Model,Gr,State,i,j,k) result(dtijk)
        integer, intent(in) :: i,j,k
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp) :: dtijk
        call CellDT_Sonic(Model,Gr,State,i,j,k,dtijk)
    end function GetCellDT_Sonic

    subroutine CellDT(Model,Gr,State,i,j,k,dtijk)
        integer, intent(in) :: i,j,k
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp),intent(out) :: dtijk
        real(rp) ::  ke,e,P, dt
        real(rp) :: Vx,Vy,Vz,rho,Bx,By,Bz
        real(rp) :: dl,MagB,MagV,Vfl,Valf,Cs,vCFL,Diff,Vdiff
        !Get three speeds, fluid/sound/alfven
        rho = State%Gas(i,j,k,DEN,BLK)
        Vx  = State%Gas(i,j,k,MOMX,BLK)/rho
        Vy  = State%Gas(i,j,k,MOMY,BLK)/rho
        Vz  = State%Gas(i,j,k,MOMZ,BLK)/rho
        MagV = sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
        
        ke = 0.5*rho*(MagV**2.0)
        e = State%Gas(i,j,k,ENERGY,BLK) - ke
        P = (Model%gamma-1)*e
        !Handle multifluid case for sound speed/max flow
        if (Model%doMultiF) then
            Cs = MultiFCs(Model,State%Gas(i,j,k,:,:))
            MagV = MultiFSpeed(Model,State%Gas(i,j,k,:,:))
        else
            Cs = sqrt(Model%gamma*P/rho)
        endif
        Vfl   = MagV
        Valf  = 0.0
        VDiff = 0.0
        if (Model%doMHD) then
            Bx = State%Bxyz(i,j,k,XDIR)
            By = State%Bxyz(i,j,k,YDIR)
            Bz = State%Bxyz(i,j,k,ZDIR)
            if (Model%doBackground) then
                Bx = Bx + Gr%B0(i,j,k,XDIR)
                By = By + Gr%B0(i,j,k,YDIR)
                Bz = Bz + Gr%B0(i,j,k,ZDIR)
            endif
            !Use B to calculate Alfven speed for CFL speed
            MagB = sqrt(Bx**2.0+By**2.0+Bz**2.0)
            Valf = MagB/sqrt(rho)
            !Boris correct Alfven speed
            if (Model%doBoris) then
                Valf = Model%Ca*Valf/sqrt(Model%Ca*Model%Ca + Valf*Valf)
            endif
            
            if(Model%doResistive) then
               ! Asume t ~ x^2/(2Diff)
               Diff = EdgeScalar2CC(Model,Gr,State%Deta,i,j,k)
               Vdiff = 2.0d0*Diff/minval((/Gr%di(i,j,k),Gr%dj(i,j,k),Gr%dk(i,j,k)/))
            end if
        endif
        vCFL = Vfl + sqrt(Cs**2.0 + Valf**2.0) + Vdiff
        
        !Use min length for timestep calculation
        dl = minval((/Gr%di(i,j,k),Gr%dj(i,j,k),Gr%dk(i,j,k)/))
        dtijk = Model%CFL*dl/vCFL
    end subroutine CellDT
    
    subroutine CellDT_Sonic(Model,Gr,State,i,j,k,dtijk)
        integer, intent(in) :: i,j,k
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp),intent(out) :: dtijk
        real(rp) ::  ke,e,P, dt
        real(rp) :: Vx,Vy,Vz,rho
        real(rp) :: dl,MagV,Cs,vCFL
        !Get three speeds, fluid/sound/alfven
        rho = State%Gas(i,j,k,DEN,BLK)
        Vx  = State%Gas(i,j,k,MOMX,BLK)/rho
        Vy  = State%Gas(i,j,k,MOMY,BLK)/rho
        Vz  = State%Gas(i,j,k,MOMZ,BLK)/rho
        MagV = sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
        
        ke = 0.5*rho*(MagV**2.0)
        e = State%Gas(i,j,k,ENERGY,BLK) - ke
        P = (Model%gamma-1)*e
        !Handle multifluid case for sound speed/max flow
        if (Model%doMultiF) then
            Cs = MultiFCs(Model,State%Gas(i,j,k,:,:))
            MagV = MultiFSpeed(Model,State%Gas(i,j,k,:,:))
        else
            Cs = sqrt(Model%gamma*P/rho)
        endif
        vCFL = Cs
        
        !Use min length for timestep calculation
        dl = minval((/Gr%di(i,j,k),Gr%dj(i,j,k),Gr%dk(i,j,k)/))
        dtijk = Model%CFL*dl/vCFL
 
    end subroutine CellDT_Sonic

    !Given GasIC_T function, initialize State variable for a given Grid
    !sOpt is optional species designation
    subroutine GasIC2State(Model,Grid,State,W,sOpt)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        procedure(GasIC_T), pointer, intent(in) :: W
        integer, intent(in), optional :: sOpt

        real(rp) :: dV, D,P,Mx,My,Mz,IntE,KinE
        real(rp) :: ijkW(NVAR)
        real(rp), dimension(8,NDIM) :: xyzC
        integer :: i,j,k,s

        if (present(sOpt)) then
            s = sOpt
        else
            s = BLK
        endif

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(xyzC,dV,ijkW,D,Mx,My,Mz,P,KinE,IntE)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=Grid%is,Grid%ie
                    !Get cell corners
                    call cellCoords(Model,Grid,i,j,k,xyzC)
                    dV = Grid%volume(i,j,k)
                    ijkW = GaussianVolumeIntegral(xyzC,W)/dV
                    
                    D = ijkW(DEN)
                    Mx = D*ijkW(VELX)
                    My = D*ijkW(VELY)
                    Mz = D*ijkW(VELZ)
                    P = ijkW(PRESSURE)
                    KinE = 0.5*(Mx**2.0 + My**2.0 + Mz**2.0)/D
                    IntE = P/(Model%gamma-1.0)

                    State%Gas(i,j,k,DEN   ,s) = D
                    State%Gas(i,j,k,MOMX  ,s) = Mx
                    State%Gas(i,j,k,MOMY  ,s) = My
                    State%Gas(i,j,k,MOMZ  ,s) = Mz
                    State%Gas(i,j,k,ENERGY,s) = IntE + KinE
                enddo
            enddo
        enddo

    end subroutine GasIC2State

    !Given vector potential routine, initialize face fluxes

    subroutine VectorPot2Flux(Model,Grid,State,Axyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        procedure(VectorField_T), pointer, intent(in) :: Axyz

        integer :: i,j,k
        real(rp), dimension(NDIM) :: e0,eI,eJ,eK,A
        real(rp), allocatable, dimension(:,:,:,:) :: lA

        call allocGridVec(Model,Grid,lA,.false.,NDIM)
        
        !Loop over grid, calculate edge-integral of A
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,e0,eI,eJ,eK)
        do k=Grid%ksg,Grid%keg-1
            do j=Grid%jsg,Grid%jeg-1
                do i=Grid%isg,Grid%ieg-1
                    !Get corner points to integrate along all 3 edges
                    e0 = Grid%xyz(i  ,j  ,k  ,:)
                    eI = Grid%xyz(i+1,j  ,k  ,:)
                    eJ = Grid%xyz(i  ,j+1,k  ,:)
                    eK = Grid%xyz(i  ,j  ,k+1,:)

                    lA(i,j,k,IDIR) = dot_product(GaussianEdgeIntegral(Model,eI,e0,Axyz),eI-e0)
                    lA(i,j,k,JDIR) = dot_product(GaussianEdgeIntegral(Model,eJ,e0,Axyz),eJ-e0)
                    lA(i,j,k,KDIR) = dot_product(GaussianEdgeIntegral(Model,eK,e0,Axyz),eK-e0)

                enddo
            enddo
        enddo
        
        !Now use edge-integrals to set fluxes
        !$OMP PARALLEL DO default(shared) collapse(2)
        do k=Grid%ksg,Grid%keg-1
            do j=Grid%jsg,Grid%jeg-1
                do i=Grid%isg,Grid%ieg-1
                    State%magFlux(i,j,k,IDIR) =   lA(i,j,k,JDIR) - lA(i,j,k+1,JDIR) + lA(i,j+1,k,KDIR) - lA(i,j,k,KDIR)
                    State%magFlux(i,j,k,JDIR) = -(lA(i,j,k,IDIR) - lA(i,j,k+1,IDIR) + lA(i+1,j,k,KDIR) - lA(i,j,k,KDIR))
                    State%magFlux(i,j,k,KDIR) =   lA(i,j,k,IDIR) - lA(i,j+1,k,IDIR) + lA(i+1,j,k,JDIR) - lA(i,j,k,JDIR)
                enddo
            enddo
        enddo
        deallocate(lA)
    end subroutine VectorPot2Flux

    !Turn Bxyz into face fluxes on grid
    subroutine VectorField2Flux(Model,Grid,State,Bxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        procedure(VectorField_T), pointer, intent(in) :: Bxyz

        integer :: i,j,k
        real(rp), dimension(NDIM) :: f0,f1,f2,f3
        real(rp) :: mFlx

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(f0,f1,f2,f3,mFlx)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg
                    
                    !Calculate face fields
                    !Need specific ordering of corner points

                    !I face
                    call faceCoords(Model,Grid,i,j,k,IDIR,f0,f1,f2,f3)
                    call GaussianFaceFlux(Model,f0,f1,f2,f3,Bxyz,mFlx)
                    State%magFlux(i,j,k,IDIR) = mFlx

                    !J face
                    call faceCoords(Model,Grid,i,j,k,JDIR,f0,f1,f2,f3)
                    call GaussianFaceFlux(Model,f0,f1,f2,f3,Bxyz,mFlx)
                    State%magFlux(i,j,k,JDIR) = mFlx

                    !K face
                    call faceCoords(Model,Grid,i,j,k,KDIR,f0,f1,f2,f3)
                    call GaussianFaceFlux(Model,f0,f1,f2,f3,Bxyz,mFlx)
                    State%magFlux(i,j,k,KDIR) = mFlx

                enddo
            enddo
        enddo

    end subroutine VectorField2Flux

    !Converts face-centered fluxes to cell-centered field components
    subroutine bFlux2Fld(Model,Grid,magFlux,magFld)

        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), dimension(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM), intent(in)  :: magFlux
        real(rp), dimension(Grid%isg:Grid%ieg  ,Grid%jsg:Grid%jeg  ,Grid%ksg:Grid%keg  ,1:NDIM), intent(out) :: magFld
        integer :: i,j,k

        !$OMP PARALLEL DO default(shared) collapse(2)
        do k=Grid%ksg, Grid%keg
            do j=Grid%jsg, Grid%jeg
                do i=Grid%isg, Grid%ieg
                    magFld(i,j,k,:) = CellBxyz(Model,Grid,magFlux,i,j,k)
                enddo
            enddo
        enddo

        if (Model%doRing) then
            call RingFlux2FieldFix(Model,Grid,magFld)
        endif
    end subroutine bFlux2Fld

! compute the current from curl B using the cell centered fleid components
! the currents are computed at cell centers, which is slightly different compared to the getParallelCurrent() function in LFM
! This function calculates J from B and cell geometry information. It is based on equation (5) in Slava's notes that computes the average value of
! J for a cell using line integrals of B around the cell edges. There's the algorithm:
!
!   1. Compute B dot dl along cell edges - put in the bint(i,j,k,DIR) array
!   2. Compute the integral of B dot dl (equals to J dot dS) for i,j,k cell faces - put in the js(i,j,k,DIR) array
!   3. Use the CellBxyz function to convert J dot dS to J at cell centers - put in the curnt(i,j,k,DIR) array
!   4. Enforce periodicity, smoothing if needed

    !jSclO = optional current scaling
    subroutine bFld2Jxyz(Model,Grid,Bxyz,Jxyz,jSclO)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), intent(in)  :: Bxyz(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), intent(inout) :: Jxyz(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), intent(in), optional :: jSclO

        integer :: i,j,k,Nk,k0
        real(rp), allocatable, dimension(:,:,:,:) :: bInt, JdS
        real(rp), dimension(NDIM) :: J0,J2,J3
        real(rp), dimension(NDIM) :: bAvg,dxyz
        real(rp) :: jScl

        if (present(jSclO)) then
            jScl = jSclO
        else
            jScl = 1.0
        endif

        call allocGridVec(Model,Grid,bInt,.false.,NDIM)
        call allocGridVec(Model,Grid,JdS ,.true. ,NDIM) !Treat like flux-sized to use CellBxyz routine

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=Grid%ksg, Grid%keg
            do j=Grid%jsg, Grid%jeg
                do i=Grid%isg, Grid%ieg
                     Jxyz(i,j,k,:) = 0.0
                enddo
            enddo
        enddo

        ! integrate B dot dl along cell edges in the i,j,k directions, respectively
        ! the integral is done in a second-order accuracy by linearly interpolating the cell-centered magnetic fields to the cell edges
        ! since the MagFld() array has all the ghost cell information, the line-integral is relatively straightforward, but the calculation
        ! in the first/last computation cell (center) is probably not meaningful.
        ! due to interpolation, calculate all values one cell away from the outer boundary
        ! loop up to the outermost cells and down-select valid edges inside the loop
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(bAvg,dxyz,i,j,k)
        do k=Grid%ksg+1, Grid%keg
            do j=Grid%jsg+1, Grid%jeg
                do i=Grid%isg+1, Grid%ieg
                    ! lower bounds are easy, upper bounds are hard
                    ! since we are calculating up to one cell from the boundary, check if one edge away is valid
                    if(isGoodEdgeIJK(Model,Grid,i+1,j,k,IDIR)) then
                        ! integrating in the i-direction along i-edges, the effective index range is isg+1:ieg-1, jsg+1:jeg, ksg+1:keg
                        bAvg = 0.25*( Bxyz(i,j,k,:) + Bxyz(i,j-1,k,:) + Bxyz(i,j,k-1,:) + Bxyz(i,j-1,k-1,:) )
                        dxyz = Grid%xyz(i+1,j,k,:) - Grid%xyz(i,j,k,:)
                        bInt(i,j,k,IDIR) = dot_product(bAvg,dxyz)
                    endif

                    if(isGoodEdgeIJK(Model,Grid,i,j+1,k,JDIR)) then
                        ! integrating in the j-direction along j-edges, the effective index range is isg+1:ieg, jsg+1:jeg-1, ksg+1:keg
                        bAvg = 0.25*( Bxyz(i,j,k,:) + Bxyz(i,j,k-1,:) + Bxyz(i-1,j,k,:) + Bxyz(i-1,j,k-1,:) )
                        dxyz = Grid%xyz(i,j+1,k,:) - Grid%xyz(i,j,k,:)
                        bInt(i,j,k,JDIR) = dot_product(bAvg,dxyz)
                    endif

                    if(isGoodEdgeIJK(Model,Grid,i,j,k+1,KDIR)) then
                        ! integrating in the k-direction along k-edges, the effective index range is isg+1:ieg, jsg+1:jeg, ksg+1:keg-1
                        bAvg = 0.25*(Bxyz(i,j,k,:) + Bxyz(i-1,j,k,:) + Bxyz(i,j-1,k,:) + Bxyz(i-1,j-1,k,:) )
                        dxyz = Grid%xyz(i,j,k+1,:) - Grid%xyz(i,j,k,:)
                        bInt(i,j,k,KDIR) = dot_product(bAvg,dxyz)
                    endif
                enddo
            enddo
        enddo

        ! line integrals of bint to get J*dS on cell faces
        ! calculate values on all faces one cell away from the outer boundary
        ! loop up to the outermost cells and down-select valid faces inside the loop
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=Grid%ksg+1, Grid%keg
            do j=Grid%jsg+1, Grid%jeg
                do i=Grid%isg+1, Grid%ieg
                    ! lower bounds are easy, upper bounds are hard
                    ! since we are calculating up to one cell from the boundary, check if one face away is valid
                    if(isGoodFaceIJK(Model,Grid,i,j+1,k+1,IDIR)) JdS(i,j,k,IDIR) = bInt(i,j,k,JDIR) + bInt(i,j+1,k,KDIR) - bInt(i,j,k+1,JDIR) - bInt(i,j,k,KDIR)
                    if(isGoodFaceIJK(Model,Grid,i+1,j,k+1,JDIR)) JdS(i,j,k,JDIR) = bInt(i,j,k,KDIR) + bInt(i,j,k+1,IDIR) - bInt(i+1,j,k,KDIR) - bInt(i,j,k,IDIR)
                    if(isGoodFaceIJK(Model,Grid,i+1,j+1,k,KDIR)) JdS(i,j,k,KDIR) = bInt(i,j,k,IDIR) + bInt(i+1,j,k,JDIR) - bInt(i,j+1,k,IDIR) - bInt(i,j,k,JDIR)
                enddo
            enddo
        enddo
        ! Jds ends up being defined isg+1,ieg-1

        ! now from J-flux to J-field using the CellBxyz function (and include scaling)
        ! calculate values at cell centers one cell away from the outer boundary
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=Grid%ksg+1, Grid%keg-1
            do j=Grid%jsg+1, Grid%jeg-1
                do i=Grid%isg+1, Grid%ieg-1
                    Jxyz(i,j,k,:) = jScl*CellBxyz(Model,Grid,JdS,i,j,k)
                enddo
            enddo
        enddo
        
        !Below here can do things related to ring-avg/periodicity
        if (Model%doRing) then
            select case (Model%Ring%GridID)
            case("lfm")
                
            !Start by cleaning up about axis
                Nk = Grid%ke-Grid%ks+1
                !$OMP PARALLEL DO default(shared) &
                !$OMP private(i,k,J0,J3,k0)
                do i=Grid%isg+1,Grid%ieg-1
                    !Positive axis
                    if (Model%Ring%doS) then
                        !Calculate J @ pole by averaging @ ring 3
                        J0 = sum(Jxyz(i,Grid%js+2,Grid%ks:Grid%ke,XDIR:ZDIR),dim=1)/Nk
                        !Now loop around ring and linearly interpolate first and second j
                        !Weighting, 0 (pole) -> 0.5 (1cc) -> 1.5 (2cc) -> 2.5 (3cc)
                        !J1 = 0.8*J0 + 0.2*J3
                        !J2 = 0.4*J0 + 0.6*J3
                        do k=Grid%ks,Grid%ke
                            J3 = Jxyz(i,Grid%js+2,k,:)
                            Jxyz(i,Grid%js  ,k,:) = 0.8*J0 + 0.2*J3
                            Jxyz(i,Grid%js+1,k,:) = 0.4*J0 + 0.6*J3
                        enddo

                        !Now fill in ghosts
                        do k=Grid%ks,Grid%ke
                            k0 = mod(k+Nk/2-1,Nk) + 1
                            Jxyz(i,Grid%js-1,k0,:) = Jxyz(i,Grid%js,  k,:)
                            Jxyz(i,Grid%js-2,k0,:) = Jxyz(i,Grid%js+1,k,:)
                            Jxyz(i,Grid%js-3,k0,:) = Jxyz(i,Grid%js+2,k,:)
                            Jxyz(i,Grid%js-4,k0,:) = Jxyz(i,Grid%js+3,k,:)
                        enddo
                        Jxyz(i,Grid%jsg:Grid%js+3,-3:0     ,:) = Jxyz(i,Grid%jsg:Grid%js+3,Nk-3:Nk,:)
                        Jxyz(i,Grid%jsg:Grid%js+3,Nk+1:Nk+4,:) = Jxyz(i,Grid%jsg:Grid%js+3,1:4    ,:)
                    endif

                    !Negative axis
                    if (Model%Ring%doE) then
                        J0 = sum(Jxyz(i,Grid%je-2,Grid%ks:Grid%ke,XDIR:ZDIR),dim=1)/Nk
                        !Same as above
                        do k=Grid%ks,Grid%ke
                            J3 = Jxyz(i,Grid%je-2,k,:)
                            Jxyz(i,Grid%je  ,k,:) = 0.8*J0 + 0.2*J3
                            Jxyz(i,Grid%je-1,k,:) = 0.4*J0 + 0.6*J3
                        enddo

                        do k=Grid%ks,Grid%ke
                            k0 = mod(k+Nk/2-1,Nk) + 1
                            Jxyz(i,Grid%je+1,k0,:) = Jxyz(i,Grid%je,  k,:)
                            Jxyz(i,Grid%je+2,k0,:) = Jxyz(i,Grid%je-1,k,:)
                            Jxyz(i,Grid%je+3,k0,:) = Jxyz(i,Grid%je-2,k,:)
                            Jxyz(i,Grid%je+4,k0,:) = Jxyz(i,Grid%je-3,k,:)
                        enddo
                        Jxyz(i,Grid%je+1:Grid%jeg,-3:0     ,:) = Jxyz(i,Grid%je+1:Grid%jeg,Nk-3:Nk,:)
                        Jxyz(i,Grid%je+1:Grid%jeg,Nk+1:Nk+4,:) = Jxyz(i,Grid%je+1:Grid%jeg,1:4    ,:)

                    endif
                enddo !i loop

            end select
        endif

        deallocate(bInt)
        deallocate(JdS)
            
    end subroutine bFld2Jxyz

!Converts edge-centered electric fields to cell-centered XYZ field components
! this subroutine computes cell-centered electric field (Ex,Ey,Ez) using edge-centered electric field components (Ei, Ej, Ek)
! Here is how the (x,y,z) components of the electric field are recovered from (i,j,k) components:
! We have a system of linear equations:
!
!                                e_i = e_x*i_x + e_y*i_y + e_z*i_z
!                                e_j = e_x*j_x + e_y*j_y + e_z*j_z                (1)
!                                e_k = e_x*k_x + e_y*k_y + e_z*k_z
!
! where e_i, e_j, e_k - components of the electric FIELD (not potential) in the (i,j,k) basis, while
! e_x,e_y,e_z - el. field in (x,y,z) basis. Correspondingly, i_x,i_y,i_z etc. are components 
! of the (i,j,k) in (x,y,z) basis.
!
! The above system has a unique solution if for every cell no 5 vertices are coplanar, i.e.
! the cells are just distorted cubes. The solution is given by matrix inversion:
!
!                                e_x   [m11 m12 m13    e_i
!                                e_y =  m21 m22 m23  * e_j                        (2)
!                                e_z    m31 m32 m33]   e_k
!
! where M is the inverse matrix in (1) given by:
! 
!                               [m11 m12 m13    [i_x i_y i_z    [1 0 0
!                                m21 m22 m23  *  j_x j_y j_z  =  0 1 0            (3)
!                                m31 m32 m33]    k_x k_y k_z]    0 0 1]
!
! this is equivalent to the calculation in MHDInnerBoundaryInterface.C:
!
!                                e_x = ( e dot (ijk_y cross ijk_z) )/ det T
!                                e_y = ( ijk_x dot (e cross ijk_z) )/ det T
!                                e_z = ( ijk_x dot (ijk_y cross e) )/ det T
!
! where T is the matrix determining (1) and det T - its determinant
!
!                                ijk_x = (i_x, j_x, k_x)
!                                ijk_y = (i_y, j_y, k_y)
!                                ijk_z = (i_z, j_z, k_z)
!                                e     = (e_i, e_j, e_k)
    subroutine Eijk2xyz(Model,Grid,Eijk,Exyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), dimension(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM), intent(in) ::  Eijk
        real(rp), dimension(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM), intent(out) :: Exyz

        integer :: i,j,k
        Exyz(:,:,:,:) = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2)
        do k=Grid%ks, Grid%ke
            do j=Grid%js, Grid%je
                do i=Grid%is, Grid%ie
                    Exyz(i,j,k,:) = CellExyz(Model,Grid,Eijk,i,j,k)
                enddo
            enddo
        enddo
    end subroutine Eijk2xyz

    subroutine Eijk2cc(Model,Grid,Eijk,ccEijk)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), dimension(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM), intent(in) ::  Eijk
        real(rp), dimension(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM), intent(out) :: ccEijk

        integer :: i,j,k
        ccEijk(:,:,:,:) = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2)
        do k=Grid%ks, Grid%ke
            do j=Grid%js, Grid%je
                do i=Grid%is, Grid%ie
                    ccEijk(i,j,k,IDIR) = 0.25*(Eijk(i,j,k,IDIR)+Eijk(i,j,k+1,IDIR)+Eijk(i,j+1,k,IDIR)+Eijk(i,j+1,k+1,IDIR))
                    ccEijk(i,j,k,JDIR) = 0.25*(Eijk(i,j,k,JDIR)+Eijk(i,j,k+1,JDIR)+Eijk(i+1,j,k,JDIR)+Eijk(i+1,j,k+1,JDIR))
                    ccEijk(i,j,k,KDIR) = 0.25*(Eijk(i,j,k,KDIR)+Eijk(i+1,j,k,KDIR)+Eijk(i,j+1,k,KDIR)+Eijk(i+1,j+1,k,KDIR))
                enddo
            enddo
        enddo
    end subroutine Eijk2cc

    !Convert single cell-centered Eijk (field not EMF) to cc-Exyz
    function ccEijk2Exyz(Model,Grid,ccEijk,i,j,k) result (Exyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        real(rp)  , intent(in) :: ccEijk(NDIM)
        integer, intent(in) :: i,j,k
        real(rp), dimension(NDIM) :: Exyz

        real(rp) :: e_i,e_j,e_k  ! cell-centered i,j,k components of the electric field
        real(rp) :: i_x,i_y,i_z,j_x,j_y,j_z,k_x,k_y,k_z
        real(rp) :: m_11,m_12,m_13,m_21,m_22,m_23,m_31,m_32,m_33,determ
        real(rp), dimension(NDIM) :: iVec,jVec,kVec
        real(rp) :: ijkDet, Mijk(NDIM,NDIM)

        ! compute the ijk_xyz vectors
        !i,j,k vectors at cell center
        iVec = ijkVec(Model,Grid,i,j,k,IDIR)
        jVec = ijkVec(Model,Grid,i,j,k,JDIR)
        kVec = ijkVec(Model,Grid,i,j,k,KDIR)

        ! compute the inverse matrix for Exyz calculations 
        call ijkMatrix(Model,Grid,iVec,jVec,kVec,Mijk)        
        
        ! determinant of the ijk Matrix
        ijkDet = iVec(XDIR)*Mijk(1,1) + jVec(XDIR)*Mijk(1,2) + kVec(XDIR)*Mijk(1,3)

        ! compute the x,y,z components of the electric field using eqn (3)
        !K: Removing overall negative
        Exyz(XDIR) = (Mijk(1,1)*ccEijk(IDIR) + Mijk(1,2)*ccEijk(JDIR) + Mijk(1,3)*ccEijk(KDIR))/ijkDet
        Exyz(YDIR) = (Mijk(2,1)*ccEijk(IDIR) + Mijk(2,2)*ccEijk(JDIR) + Mijk(2,3)*ccEijk(KDIR))/ijkDet
        Exyz(ZDIR) = (Mijk(3,1)*ccEijk(IDIR) + Mijk(3,2)*ccEijk(JDIR) + Mijk(3,3)*ccEijk(KDIR))/ijkDet
        
    end function ccEijk2Exyz

    !Converts edge IJK *FIELDS* (not EMFs) to cc-Exyz
    function CellExyz(Model,Grid,Eijk,i,j,k) result (Exyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        real(rp)  , intent(in) :: Eijk(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        integer, intent(in) :: i,j,k

        real(rp), dimension(NDIM) :: Exyz
        real(rp), dimension(NDIM) :: ccEijk
        
        ! first compute the ijk component of electric field at cell center through linear interpolation
        ccEijk(IDIR) = 0.25*(Eijk(i,j,k,IDIR)+Eijk(i,j,k+1,IDIR)+Eijk(i,j+1,k,IDIR)+Eijk(i,j+1,k+1,IDIR))
        ccEijk(JDIR) = 0.25*(Eijk(i,j,k,JDIR)+Eijk(i,j,k+1,JDIR)+Eijk(i+1,j,k,JDIR)+Eijk(i+1,j,k+1,JDIR))
        ccEijk(KDIR) = 0.25*(Eijk(i,j,k,KDIR)+Eijk(i+1,j,k,KDIR)+Eijk(i,j+1,k,KDIR)+Eijk(i+1,j+1,k,KDIR))
       
        Exyz = ccEijk2Exyz(Model,Grid,ccEijk,i,j,k)

    end function CellExyz

    !Converts edge IJK *EMFS* to cc-Exyz
    function emf2xyz(Model,Grid,Efld,i,j,k) result (Exyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        real(rp)  , intent(in) :: Efld(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        integer, intent(in) :: i,j,k

        real(rp), dimension(NDIM) :: Exyz
        real(rp), dimension(NDIM) :: ccEijk
        
        ! first compute the ijk component of electric field at cell center through linear interpolation
        ccEijk(IDIR) = 0.25*(Efld(i  ,j  ,k,IDIR)/Grid%edge(i  ,j  ,k,IDIR) + Efld(i  ,j  ,k+1,IDIR)/Grid%edge(i  ,j  ,k+1,IDIR) &
                           + Efld(i  ,j+1,k,IDIR)/Grid%edge(i  ,j+1,k,IDIR) + Efld(i  ,j+1,k+1,IDIR)/Grid%edge(i  ,j+1,k+1,IDIR))
        ccEijk(JDIR) = 0.25*(Efld(i  ,j  ,k,JDIR)/Grid%edge(i  ,j  ,k,JDIR) + Efld(i  ,j  ,k+1,JDIR)/Grid%edge(i  ,j  ,k+1,JDIR) &
                           + Efld(i+1,j  ,k,JDIR)/Grid%edge(i+1,j  ,k,JDIR) + Efld(i+1,j  ,k+1,JDIR)/Grid%edge(i+1,j  ,k+1,JDIR))
        ccEijk(KDIR) = 0.25*(Efld(i  ,j  ,k,KDIR)/Grid%edge(i  ,j  ,k,KDIR) + Efld(i+1,j  ,k  ,KDIR)/Grid%edge(i+1,j  ,k  ,KDIR) &
                           + Efld(i  ,j+1,k,KDIR)/Grid%edge(i  ,j+1,k,KDIR) + Efld(i+1,j+1,k  ,KDIR)/Grid%edge(i+1,j+1,k  ,KDIR))
       
        Exyz = ccEijk2Exyz(Model,Grid,ccEijk,i,j,k)
    end function emf2xyz

    !Calculates fluxes in cell ijk based on Bxyz
    function CellFld2Flux(Model,Gr,Bxyz,i,j,k) result(Flx)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Gr
        real(rp), intent(in) :: Bxyz(NDIM)
        integer, intent(in) :: i,j,k

        real(rp) :: Flx(NDIM)
        Flx(IDIR) = Gr%face(i,j,k,IDIR)*dot_product(Bxyz,Gr%Tf(i,j,k,NORMX:NORMZ,IDIR))
        Flx(JDIR) = Gr%face(i,j,k,JDIR)*dot_product(Bxyz,Gr%Tf(i,j,k,NORMX:NORMZ,JDIR))
        Flx(KDIR) = Gr%face(i,j,k,KDIR)*dot_product(Bxyz,Gr%Tf(i,j,k,NORMX:NORMZ,KDIR))

    end function CellFld2Flux
    
    !Calculates the XYZ field in a cell from flux
    !Gamera vesion
    function gamCellBxyz(Model,Grid,magFlux,i,j,k) result(Bxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        real(rp)  , intent(in) :: magFlux(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM)
        integer, intent(in) :: i,j,k
        real(rp), dimension(NDIM) :: Bxyz

        real(rp) :: dV
        real(rp) :: fI,fJ,fK,fIp,fJp,fKp
        real(rp), dimension(NDIM) :: xfI,xfIp,xfJ,xfJp,xfK,xfKp
        real(rp) :: Div, xcc(NDIM)
        dV = Grid%volume(i,j,k)
        fI  = magFlux(i,j,k,IDIR)
        fJ  = magFlux(i,j,k,JDIR)
        fK  = magFlux(i,j,k,KDIR)
        fIp = magFlux(i+1,j,k,IDIR)
        fJp = magFlux(i,j+1,k,JDIR)
        fKp = magFlux(i,j,k+1,KDIR)

        xfI  = Grid%xfc(i,j,k,:,IDIR)
        xfJ  = Grid%xfc(i,j,k,:,JDIR)
        xfK  = Grid%xfc(i,j,k,:,KDIR)
        xfIp = Grid%xfc(i+1,j,k,:,IDIR)
        xfJp = Grid%xfc(i,j+1,k,:,JDIR)
        xfKp = Grid%xfc(i,j,k+1,:,KDIR)

        Div = fIp - fI + fJp - fJ + fKp - fK

        Bxyz = (fIp*xfIp + fJp*xfJp + fKp*xfKp - fI*xfI - fJ*xfJ - fK*xfK)/dV  
        Bxyz = Bxyz - Div*Grid%xyzcc(i,j,k,:)/dV

        
    end function gamCellBxyz

    !Calculates the XYZ field in a cell
    function lfmCellBxyz(Model,Grid,magFlux,i,j,k) result(Bxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        real(rp)  , intent(in) :: magFlux(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM)
        integer, intent(in) :: i,j,k
        real(rp), dimension(NDIM) :: Bxyz

        real(rp) :: dV
        real(rp) :: fI,fJ,fK,fIp,fJp,fKp
        real(rp), dimension(NDIM) :: xfI,xfIp,xfJ,xfJp,xfK,xfKp
        real(rp) :: Div, xcc(NDIM)

        dV = Grid%volume(i,j,k)
        fI  = magFlux(i,j,k,IDIR)
        fJ  = magFlux(i,j,k,JDIR)
        fK  = magFlux(i,j,k,KDIR)
        fIp = magFlux(i+1,j,k,IDIR)
        fJp = magFlux(i,j+1,k,JDIR)
        fKp = magFlux(i,j,k+1,KDIR)

        !Use coordinate centers, not face centers
        !IE, need 4-point averages
        xfI = 0.25*( Grid%xyz(i,j,k,:) + Grid%xyz(i,j+1,k,:) + Grid%xyz(i,j+1,k+1,:) + Grid%xyz(i,j,k+1,:) )
        xfJ = 0.25*( Grid%xyz(i,j,k,:) + Grid%xyz(i+1,j,k,:) + Grid%xyz(i+1,j,k+1,:) + Grid%xyz(i,j,k+1,:) )
        xfK = 0.25*( Grid%xyz(i,j,k,:) + Grid%xyz(i+1,j,k,:) + Grid%xyz(i+1,j+1,k,:) + Grid%xyz(i,j+1,k,:) )
        
        xfIp = 0.25*( Grid%xyz(i+1,j,k,:) + Grid%xyz(i+1,j+1,k,:) + Grid%xyz(i+1,j+1,k+1,:) + Grid%xyz(i+1,j,k+1,:) )
        xfJp = 0.25*( Grid%xyz(i,j+1,k,:) + Grid%xyz(i+1,j+1,k,:) + Grid%xyz(i+1,j+1,k+1,:) + Grid%xyz(i,j+1,k+1,:) )
        xfKp = 0.25*( Grid%xyz(i,j,k+1,:) + Grid%xyz(i+1,j,k+1,:) + Grid%xyz(i+1,j+1,k+1,:) + Grid%xyz(i,j+1,k+1,:) )
    
        Bxyz = 0.5*( (fIp+fI)*(xfIp-xfI) + (fJp+fJ)*(xfJp-xfJ) + (fKp+fK)*(xfKp-xfK) )/dV

    end function lfmCellBxyz

    subroutine DivB(Model,Grid,State,totDivB,DivOut,doTotO)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        type(State_T), intent(in) :: State
        real(rp), intent(out) :: totDivB
        real(rp), intent(inout), optional :: DivOut(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg)
        logical, intent(in), optional :: doTotO

        logical :: doTot
        real(rp) :: DivBcc(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg)
        integer :: i,j,k
        real(rp) :: dFi,dFj,dFk,dFi0,dFj0,dFk0

        if (present(doTotO)) then
            doTot = doTotO
        else
            doTot = .true.
        endif

        dFi0 = 0.0
        dFj0 = 0.0
        dFk0 = 0.0
        totDivB = 0.0

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(dFi,dFj,dFk,dFi0,dFj0,dFk0) &
        !$OMP reduction(+:totDivB)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=Grid%is,Grid%ie
                    dFi0 = 0.0
                    dFj0 = 0.0
                    dFk0 = 0.0
                                  
                    dFi = State%magFlux(i,j,k,IDIR) - State%magFlux(i+1,j  ,k  ,IDIR)
                    dFj = State%magFlux(i,j,k,JDIR) - State%magFlux(i  ,j+1,k  ,JDIR)
                    dFk = State%magFlux(i,j,k,KDIR) - State%magFlux(i  ,j  ,k+1,KDIR)
                    if (Model%doBackground .and. doTot) then
                        dFi0 = Grid%bFlux0(i,j,k,IDIR) - Grid%bFlux0(i+1,j  ,k  ,IDIR)
                        dFj0 = Grid%bFlux0(i,j,k,JDIR) - Grid%bFlux0(i  ,j+1,k  ,JDIR)
                        dFk0 = Grid%bFlux0(i,j,k,KDIR) - Grid%bFlux0(i  ,j  ,k+1,KDIR)
                    endif
                    DivBcc(i,j,k) = abs(dFi+dFj+dFk+dFi0+dFj0+dFk0)

                    totDivB = totDivB + DivBcc(i,j,k)

                enddo
            enddo
        enddo

        if (present(DivOut)) then
            DivOut = DivBcc
        endif
    end subroutine DivB

    !Calculates cell-centered value of a quantity defined at all cell edges
    !NOTE: Assuming GridVec size w/ +1 (see allocGridVec)
    function EdgeScalar2CC(Model,Grid,Qedg,i,j,k) result(Qcc)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        real(rp)  , intent(in) :: Qedg(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM)
        integer   , intent(in) :: i,j,k
        real(rp) :: Qcc

        !Get 12-pt average for cell-centered value
        Qcc = Qedg(i  ,j  ,k  ,IDIR) + Qedg(i+1,j  ,k  ,JDIR) + Qedg(i  ,j+1,k  ,IDIR) + Qedg(i  ,j  ,k  ,JDIR)  + & !bottom k face
              Qedg(i  ,j  ,k+1,IDIR) + Qedg(i+1,j  ,k+1,JDIR) + Qedg(i  ,j+1,k+1,IDIR) + Qedg(i  ,j  ,k+1,JDIR)  + & !top k face
              Qedg(i  ,j  ,k  ,KDIR) + Qedg(i+1,j  ,k  ,KDIR) + Qedg(i+1,j+1,k  ,KDIR) + Qedg(i  ,j+1,k  ,KDIR)
        Qcc = Qcc/12.0
        
    end function EdgeScalar2CC

    !Fill in the final rank of magflux by utilizing divB=0
    subroutine FixGhostFace(Model,Grid,magFlux)
        type(Model_T), intent(in)    :: Model
        type(Grid_T) , intent(in)    :: Grid
        real(rp)     , intent(inout) :: magFlux(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM)

        integer :: i,j,k
        real(rp) :: dFi,dFj,dFk,dFi0,dFj0,dFk0

        !For cell i,j,k
        !dFi  = magFlux(i,j,k,IDIR) - magFlux(i+1,j  ,k  ,IDIR)
        !dFj  = magFlux(i,j,k,JDIR) - magFlux(i  ,j+1,k  ,JDIR)
        !dFk  = magFlux(i,j,k,KDIR) - magFlux(i  ,j  ,k+1,KDIR)
        !dFi0 = Grid%bFlux0  (i,j,k,IDIR) - Grid%bFlux0  (i+1,j  ,k  ,IDIR)
        !dFj0 = Grid%bFlux0  (i,j,k,JDIR) - Grid%bFlux0  (i  ,j+1,k  ,JDIR)
        !dFk0 = Grid%bFlux0  (i,j,k,KDIR) - Grid%bFlux0  (i  ,j  ,k+1,KDIR)
        !div(dB) + div(B0) = 0 = dFi + dFj + dFk + dFi0 + dFj0 + dFk0
        !Note: This only matters in the small region where div(B0) != 0 but div(dB+B0) = 0

        !i-dir
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,dFi,dFj,dFk,dFi0,dFj0,dFk0)
        do k=Grid%ksg, Grid%keg
            do j=Grid%jsg, Grid%jeg
                !Want ieg+1,j,k I-face
                i = Grid%ieg
                dFi  = magFlux(i,j,k,IDIR) - magFlux(i+1,j  ,k  ,IDIR)
                dFj  = magFlux(i,j,k,JDIR) - magFlux(i  ,j+1,k  ,JDIR)
                dFk  = magFlux(i,j,k,KDIR) - magFlux(i  ,j  ,k+1,KDIR)
                if (Model%doBackground) then
                    dFi0 = Grid%bFlux0(i,j,k,IDIR) - Grid%bFlux0(i+1,j  ,k  ,IDIR)
                    dFj0 = Grid%bFlux0(i,j,k,JDIR) - Grid%bFlux0(i  ,j+1,k  ,JDIR)
                    dFk0 = Grid%bFlux0(i,j,k,KDIR) - Grid%bFlux0(i  ,j  ,k+1,KDIR)
                else
                    dFi0 = 0.0
                    dFj0 = 0.0
                    dFk0 = 0.0
                endif
                magFlux(i+1,j,k,IDIR) = magFlux(i,j,k,IDIR) + dFj + dFk &
                                        + dFi0 + dFj0 + dFk0
            enddo
        enddo

        !j-dir
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,dFi,dFj,dFk,dFi0,dFj0,dFk0)
        do k=Grid%ksg, Grid%keg
            do i=Grid%isg, Grid%ieg
                j = Grid%jeg
                dFi  = magFlux(i,j,k,IDIR) - magFlux(i+1,j  ,k  ,IDIR)
                dFj  = magFlux(i,j,k,JDIR) - magFlux(i  ,j+1,k  ,JDIR)
                dFk  = magFlux(i,j,k,KDIR) - magFlux(i  ,j  ,k+1,KDIR)
                if (Model%doBackground) then
                    dFi0 = Grid%bFlux0(i,j,k,IDIR) - Grid%bFlux0(i+1,j  ,k  ,IDIR)
                    dFj0 = Grid%bFlux0(i,j,k,JDIR) - Grid%bFlux0(i  ,j+1,k  ,JDIR)
                    dFk0 = Grid%bFlux0(i,j,k,KDIR) - Grid%bFlux0(i  ,j  ,k+1,KDIR)
                else
                    dFi0 = 0.0
                    dFj0 = 0.0
                    dFk0 = 0.0
                endif
                magFlux(i,j+1,k,JDIR) = magFlux(i,j,k,JDIR) + dFi + dFk &
                                        + dFi0 + dFj0 + dFk0
            enddo
        enddo   

        !k-dir
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,dFi,dFj,dFk,dFi0,dFj0,dFk0)
        do j=Grid%jsg, Grid%jeg
            do i=Grid%isg, Grid%ieg
                k = Grid%keg
                dFi  = magFlux(i,j,k,IDIR) - magFlux(i+1,j  ,k  ,IDIR)
                dFj  = magFlux(i,j,k,JDIR) - magFlux(i  ,j+1,k  ,JDIR)
                dFk  = magFlux(i,j,k,KDIR) - magFlux(i  ,j  ,k+1,KDIR)
                if (Model%doBackground) then
                    dFi0 = Grid%bFlux0(i,j,k,IDIR) - Grid%bFlux0(i+1,j  ,k  ,IDIR)
                    dFj0 = Grid%bFlux0(i,j,k,JDIR) - Grid%bFlux0(i  ,j+1,k  ,JDIR)
                    dFk0 = Grid%bFlux0(i,j,k,KDIR) - Grid%bFlux0(i  ,j  ,k+1,KDIR)
                else
                    dFi0 = 0.0
                    dFj0 = 0.0
                    dFk0 = 0.0
                endif
                magFlux(i,j,k+1,KDIR) = magFlux(i,j,k,KDIR) + dFi + dFk &
                                        + dFi0 + dFj0 + dFk0
            enddo
        enddo

    end subroutine FixGhostFace
    
    subroutine EnergyPartition(Model,Gr,State,gIntE,gKinE,gMagP)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State

        real(rp), intent(inout) :: gKinE,gMagP,gIntE

        integer :: i,j,k
        real(rp) :: Vx,Vy,Vz,D,Bx,By,Bz
        real(rp) :: KinE, IntE, MagB, Vijk

        !Reset energies
        gKinE=0.0
        gMagP=0.0
        gIntE=0.0

        !Loop over active cells
        do k=Gr%ksg, Gr%keg
            do j=Gr%jsg, Gr%jeg
                do i=Gr%isg, Gr%ieg
                    Vijk = Gr%volume(i,j,k)
                    D   = State%Gas(i,j,k,DEN,BLK)
                    Vx  = State%Gas(i,j,k,MOMX,BLK)/D
                    Vy  = State%Gas(i,j,k,MOMY,BLK)/D
                    Vz  = State%Gas(i,j,k,MOMZ,BLK)/D

                    KinE = 0.5*D*(Vx**2.0+Vy**2.0+Vz**2.0)
                    IntE = State%Gas(i,j,k,ENERGY,BLK) - KinE
                    MagB = 0.0

                    if (Model%doMHD) then
                        Bx = State%Bxyz(i,j,k,XDIR)
                        By = State%Bxyz(i,j,k,YDIR)
                        Bz = State%Bxyz(i,j,k,ZDIR)
                        if (Model%doBackground) then
                            Bx = Bx + Gr%B0(i,j,k,XDIR)
                            By = By + Gr%B0(i,j,k,YDIR)
                            Bz = Bz + Gr%B0(i,j,k,ZDIR)
                        endif
    
                        MagB = sqrt(Bx**2.0+By**2.0+Bz**2.0)
                    endif
                    !Accumulate
                    gIntE = gIntE + IntE*Vijk
                    gKinE = gKinE + KinE*Vijk
                    gMagP = gMagP + 0.5*MagB*MagB*Vijk
                enddo
            enddo
        enddo

    end subroutine EnergyPartition

    subroutine getJxyz(Model,Gr,State,Jxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp), dimension(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM), intent(inout) :: Jxyz
        
        real (rp), dimension(:,:,:,:), allocatable :: VecA  !! dummy to hold reduced bfld
        integer :: i,j,k

        call allocGridVec(Model,Gr,VecA)

        !!! Get the current
        if (Model%isMagsphere) then
            !Subtract dipole before calculating current
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k)
            do k=Gr%ksg,Gr%keg
                do j=Gr%jsg,Gr%jeg
                    do i=Gr%isg,Gr%ieg
                        if (Model%doBackground) then
                            VecA(i,j,k,:) = State%Bxyz(i,j,k,:) + Gr%B0(i,j,k,:) - MagsphereDipole(Gr%xyzcc(i,j,k,:),Model%MagM0)
                        else
                            VecA(i,j,k,:) = State%Bxyz(i,j,k,:) - MagsphereDipole(Gr%xyzcc(i,j,k,:),Model%MagM0)
                        endif
                    enddo
                enddo
            enddo
        else
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k)
            do k=Gr%ksg,Gr%keg
                do j=Gr%jsg,Gr%jeg
                    do i=Gr%isg,Gr%ieg
                        if (Model%doBackground) then
                            VecA(i,j,k,:) = State%Bxyz(i,j,k,:) + Gr%B0(i,j,k,:)
                        else
                            VecA(i,j,k,:) = State%Bxyz(i,j,k,:)
                        endif
                    enddo
                enddo
            enddo
        endif

        call bFld2Jxyz(Model,Gr,VecA,Jxyz)

        if ( Model%do25d ) then
            do k=1,4
                Jxyz(:,:,-4+k,:) = Jxyz(:,:,1,:)
                JXyz(:,:,k+1 ,:) = Jxyz(:,:,1,:)
            enddo
        endif
    end subroutine getJxyz

    function isGoodFaceIJK(Model,Grid,i,j,k,d) result(isGood)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        integer      , intent(in) :: i,j,k,d

        logical :: isGood

        isGood = .false.
        select case(d)
            case(IDIR)
                isGood =     (i >= Grid%isg) .and. (i <= Grid%ieg+1) &
                       .and. (j >= Grid%jsg) .and. (j <= Grid%jeg  ) &
                       .and. (k >= Grid%ksg) .and. (k <= Grid%keg  )

            case(JDIR)
                isGood =     (i >= Grid%isg) .and. (i <= Grid%ieg  ) &
                       .and. (j >= Grid%jsg) .and. (j <= Grid%jeg+1) &
                       .and. (k >= Grid%ksg) .and. (k <= Grid%keg  )
             case(KDIR)
                isGood =     (i >= Grid%isg) .and. (i <= Grid%ieg  ) &
                       .and. (j >= Grid%jsg) .and. (j <= Grid%jeg  ) &
                       .and. (k >= Grid%ksg) .and. (k <= Grid%keg+1)
        end select

    end function isGoodFaceIJK

    function isGoodEdgeIJK(Model,Grid,i,j,k,d) result(isGood)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        integer      , intent(in) :: i,j,k,d

        logical :: isGood

        isGood = .false.
        select case(d)
            case(IDIR)
                isGood =     (i >= Grid%isg) .and. (i <= Grid%ieg  ) &
                       .and. (j >= Grid%jsg) .and. (j <= Grid%jeg+1) &
                       .and. (k >= Grid%ksg) .and. (k <= Grid%keg+1)

            case(JDIR)
                isGood =     (i >= Grid%isg) .and. (i <= Grid%ieg+1) &
                       .and. (j >= Grid%jsg) .and. (j <= Grid%jeg  ) &
                       .and. (k >= Grid%ksg) .and. (k <= Grid%keg+1)
             case(KDIR)
                isGood =     (i >= Grid%isg) .and. (i <= Grid%ieg+1) &
                       .and. (j >= Grid%jsg) .and. (j <= Grid%jeg+1) &
                       .and. (k >= Grid%ksg) .and. (k <= Grid%keg  )
        end select

    end function isGoodEdgeIJK


end module gridutils
