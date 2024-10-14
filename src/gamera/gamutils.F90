module gamutils
    use gamtypes
    use xml_input
    use math

    implicit none
    
    contains

!Front end routines to brickette loaders
!NOTE: xxxLoadBlock assumes that memory being pulled from is aligned (ie cell-centered data
!      xxxLoadBlockI does not assume alignment, is for interface/edge variables

    !Loads a brickette of size (iMax,recLen) 
    !Into a brickette of size (vecLen,recLen)
    subroutine recLoadBlock(Model,Gr,Qb,Q,iS,j,k,iMax,d)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: iS,j,k,iMax,d
        !Remap Q to 1D from (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), intent(in ) :: Q (Gr%Ni*Gr%Nj*Gr%Nk) 
        real(rp), intent(out) :: Qb(vecLen,-Nr2:Nr2-1)

        call gamLoadBlock(Model,Gr,Qb,Q,iS,j,k,iMax,d,Nr2)
    end subroutine recLoadBlock

    !Same as above but for limiter stencil instead of reconstruction stencil
    subroutine limLoadBlock(Model,Gr,Qb,Q,iS,j,k,iMax,d)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: iS,j,k,iMax,d
        !Remap Q to 1D from (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), intent(in ) :: Q (Gr%Ni*Gr%Nj*Gr%Nk) 
        real(rp), intent(out) :: Qb(vecLen,-Nl2:Nl2-1)

        call gamLoadBlock(Model,Gr,Qb,Q,iS,j,k,iMax,d,Nl2)
    end subroutine limLoadBlock

    subroutine recLoadBlockI(Model,Gr,Qb,Q,iS,j,k,iMax,d)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: iS,j,k,iMax,d
        !Remap Q to 1D from (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), intent(in ) :: Q ( (Gr%Ni+1)*(Gr%Nj+1)*(Gr%Nk+1) ) 
        real(rp), intent(out) :: Qb(vecLen,-Nr2:Nr2-1)

        call gamLoadBlockI(Model,Gr,Qb,Q,iS,j,k,iMax,d,Nr2)
    end subroutine recLoadBlockI

    subroutine limLoadBlockI(Model,Gr,Qb,Q,iS,j,k,iMax,d)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: iS,j,k,iMax,d
        !Remap Q to 1D from (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), intent(in ) :: Q ( (Gr%Ni+1)*(Gr%Nj+1)*(Gr%Nk+1) ) 
        real(rp), intent(out) :: Qb(vecLen,-Nl2:Nl2-1)

        call gamLoadBlockI(Model,Gr,Qb,Q,iS,j,k,iMax,d,Nl2)
    end subroutine limLoadBlockI

    !----
    !Core block routines
    subroutine gamLoadBlock(Model,Gr,Qb,Q,iS,j,k,iMax,d,Nx2)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: iS,j,k,iMax,d,Nx2
        !Remap Q to 1D from (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), intent(in ) :: Q(Gr%Ni*Gr%Nj*Gr%Nk) 
        real(rp), intent(out) :: Qb(vecLen,-Nx2:Nx2-1)
        
        integer :: l,n,l0,nSi,nSj,nSk,nSt
        !DIR$ ASSUME_ALIGNED Qb: ALIGN
        !DIR$ ASSUME_ALIGNED Q: ALIGN
            
        !Mapping of iS,j,k to 1D (l0)
        l0 = ijk2n(Gr,iS,j,k)

        !Strides dep. on direction
        nSi = 1; nSj = Gr%Ni; nSk = Gr%Ni*Gr%Nj

        !Set nSt equal to nSi/nSj/nSk for given direction
        nSt = nSi*ijkD(d,IDIR) + nSj*ijkD(d,JDIR) + nSk*ijkD(d,KDIR)
        Qb = 0.0
        do n = -Nx2,Nx2-1
            l = l0 + nSt*n !Jump to next 1:iMax in direction d
            Qb(1:iMax,n) = Q(l:l+iMax-1)

        enddo

        !In comments: long-form loop using Q(isg:,jsg:,ksg:)
        !Does same thing, slower but more readable
        ! do n = -Nr2,Nr2-1
        !     do i=1,iMax
        !         iG = i+iS-1 !Global index into Q
        !         ip = iG + n*ijkD(d,IDIR)
        !         jp = j  + n*ijkD(d,JDIR)
        !         kp = k  + n*ijkD(d,KDIR)
        !         Qb(i,n) = Q(ip,jp,kp)
        !     enddo
        ! enddo   

    end subroutine gamLoadBlock

    !Same as LoadBlock but for interface-centered variables
    !Not necessarily aligned
    subroutine gamLoadBlockI(Model,Gr,Qb,Q,iS,j,k,iMax,d,Nx2)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: iS,j,k,iMax,d,Nx2
        !Remap Q to 1D from (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), intent(in) :: Q( (Gr%Ni+1)*(Gr%Nj+1)*(Gr%Nk+1) ) 
        real(rp), intent(out) :: Qb(vecLen,-Nx2:Nx2-1)
        
        integer :: l,n,l0,nSi,nSj,nSk,nSt
            
        !Mapping of iS,j,k to 1D (l0)
        l0 = ijk2nI(Gr,iS,j,k)

        !Strides dep. on direction
        nSi = 1; nSj = Gr%Ni+1; nSk = (Gr%Ni+1)*(Gr%Nj+1)

        !Set nSt equal to nSi/nSj/nSk for given direction
        nSt = nSi*ijkD(d,IDIR) + nSj*ijkD(d,JDIR) + nSk*ijkD(d,KDIR)
        Qb = 0.0
        do n = -Nx2,Nx2-1
            l = l0 + nSt*n !Jump to next 1:iMax in direction d
            Qb(1:iMax,n) = Q(l:l+iMax-1)

        enddo
    end subroutine gamLoadBlockI

    !Converts i,j,k index (isg:ieg,jsg:jeg,ksg:keg)
    !To 1D index, n= 1:Ni*Nj*Nk)

    function ijk2n(Gr,i,j,k) result (n)
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: i,j,k
        integer :: n

        integer :: ip,jp,kp

        !Convert to ip,jp,kp (1:Ni,j,k)
        ip = (i-Gr%isg)+1
        jp = (j-Gr%jsg)+1
        kp = (k-Gr%ksg)+1

        !Convert to n
        n = ip + (jp-1)*Gr%Ni + (kp-1)*(Gr%Ni*Gr%Nj)

    end function ijk2n

    function ijk2nI(Gr,i,j,k) result (n)
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: i,j,k
        integer :: n

        integer :: ip,jp,kp

        !Convert to ip,jp,kp (1:Ni,j,k)
        ip = (i-Gr%isg)+1
        jp = (j-Gr%jsg)+1
        kp = (k-Gr%ksg)+1

        !Convert to n
        n = ip + (jp-1)*(Gr%Ni+1) + (kp-1)*( (Gr%Ni+1)*(Gr%Nj+1) )

    end function ijk2nI

    subroutine allocState(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        type(State_T), intent(inout) :: State

        if ( .not. allocated(State%Gas) ) then
            allocate( State%Gas(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NVAR,0:Model%nSpc) )
            State%Gas = 0.0
        endif

        if ( .not. allocated(State%magFlux) ) then
            call allocGridVec(Model,Grid,State%magFlux,.true.,NDIM)
            State%magFlux = 0.0
        endif
        
        if ( .not. allocated(State%Efld) ) then
            call allocGridVec(Model,Grid,State%Efld,.false.,NDIM)
            State%Efld = 0.0
        endif

        if ( Model%doResistive .and. .not. allocated(State%Deta) ) then
            call allocGridVec(Model,Grid,State%Deta,.true.,NDIM)
            State%Deta = 0.0
        endif

        if ( .not. allocated(State%Bxyz) ) then
            call allocGridVec(Model,Grid,State%Bxyz,.false.,NDIM)
            State%Bxyz = 0.0
        endif

        State%time = Model%t

    end subroutine allocState

    subroutine deallocState(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        type(State_T), intent(inout) :: State
        
        if ( allocated(State%Gas)     ) deallocate(State%Gas)
        if ( allocated(State%magFlux) ) deallocate(State%magFlux)
        if ( allocated(State%Efld)    ) deallocate(State%Efld)
        if ( allocated(State%Bxyz)    ) deallocate(State%Bxyz)
        if ( allocated(State%Deta)    ) deallocate(State%Deta)

    end subroutine deallocState
    
    !Creates corner array
    !Note: Assumes is/ie,etc are already set
    subroutine allocGrid(Model,Grid)
        type(Model_T), intent(in)    :: Model
        type(Grid_T),  intent(inout) :: Grid

        if ( allocated(Grid%xyz) ) then
            deallocate(Grid%xyz)
        endif

        !Create big memory
        allocate(Grid%xyz(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,XDIR:ZDIR))
        Grid%xyz = 0.0

    end subroutine allocGrid

    !Allocates space for a grid-sized variable (w/ ghosts)
    !If doP1, create all dims+1
    subroutine allocGridVar(Model,Grid,gV,doP1)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: gV
        logical, optional, intent(in) :: doP1

        logical :: doPlus
        if ( .not. allocated(gV) ) then
            !Choose size, either all xlg:xhg or xlg:xhg+1
            if (present(doP1)) then
                doPlus = doP1
            else
                doPlus = .false.
            endif

            if (doPlus) then
                allocate( gV(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1) )
            else
                allocate( gV(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg) )
            endif
            !Go ahead and wipe newly allocated data
            gV = 0.0
            
        endif

    end subroutine allocGridVar

    !Allocates space for a grid-sized vector (w/ ghosts)
    !If doP1, create all dims+1
    !numDim is the number of vector components (ie, for state vector or spatial vector)
    subroutine allocGridVec(Model,Grid,gV,doP1,numDim)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: gV
        logical, optional, intent(in) :: doP1
        integer, optional, intent(in) :: numDim

        logical :: doPlus
        integer :: nD

        if ( .not. allocated(gV) ) then
            !Choose size, either all xlg:xhg or xlg:xhg+1
            if (present(doP1)) then
                doPlus = doP1
            else
                doPlus = .false.
            endif
            !Choose number of elements in vector, default to NDIM
            if (present(numDim)) then
                nD = numDim
            else
                nD = NDIM
            endif

            if (doPlus) then
                allocate( gV(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:nD ) )
            else
                allocate( gV(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:nD) )
            endif
            
            !Go ahead and wipe newly allocated data
            gV = 0.0
        endif       

    end subroutine allocGridVec

    !Convert conserved to primitive variables
    subroutine CellC2P(Model,Con,Prim)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in)  :: Con(NVAR)
        real(rp), intent(out) :: Prim(NVAR)

        real(rp) :: rho,Vx,Vy,Vz,E,KinE

        rho = Con(DEN)
        Vx  = Con(MOMX)/rho
        Vy  = Con(MOMY)/rho
        Vz  = Con(MOMZ)/rho
        E   = Con(ENERGY)

        KinE = 0.5*rho*(Vx**2.0+Vy**2.0+Vz**2.0)

        Prim(DEN)      = max(rho,dFloor)
        Prim(VELX)     = Vx
        Prim(VELY)     = Vy
        Prim(VELZ)     = Vz
        Prim(PRESSURE) = max((Model%gamma-1)*(E-KinE),pFloor)
    end subroutine CellC2P

    !Calculate sound speed (Cs) from conserved variables
    subroutine CellPress2Cs(Model,Con,Cs)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in)  :: Con(NVAR)
        real(rp), intent(out) :: Cs

        real(rp) :: rho,Vx,Vy,Vz,E,KinE,P

        rho = max(Con(DEN),dFloor)
        Vx  = Con(MOMX)/rho
        Vy  = Con(MOMY)/rho
        Vz  = Con(MOMZ)/rho
        E   = Con(ENERGY)

        KinE = 0.5*rho*(Vx**2.0+Vy**2.0+Vz**2.0)

        P = max((Model%gamma-1)*(E-KinE),pFloor)
        Cs = sqrt(Model%gamma*P/rho)
    end subroutine CellPress2Cs

    subroutine Con2Fast(Model,gCon,B,Cx)
        type(Model_T)  , intent(in)  :: Model
        real(rp)       , intent(in)  :: gCon(NVAR),B(NDIM)
        real(rp)       , intent(out) :: Cx
        real(rp) :: gW(NVAR)
        real(rp) :: D,Cs,Va,MagB
        call CellC2P(Model,gCon,gW) !to primitive
        !Get sound speed
        call CellPress2Cs(Model,gCon,Cs)
        !Get alfven speed
        D = gCon(DEN)
        MagB = norm2(B)
        Va = MagB/sqrt(D)
        if (Model%doBoris) then
            Va = Model%Ca*Va/sqrt(Model%Ca**2.0 + Va**2.0)
        endif
        Cx = sqrt(Cs**2.0 + Va**2.0)
    end subroutine Con2Fast
    
    subroutine CellP2C(Model,Prim,Con)
        type(Model_T), intent(in) :: Model
        real(rp), intent(out) :: Con(NVAR)
        real(rp), intent(in)  :: Prim(NVAR)

        real(rp) :: rho,Vx,Vy,Vz,KinE,IntE

        rho  = max(Prim(DEN),dFloor)
        Vx   = Prim(VELX)
        Vy   = Prim(VELY)
        Vz   = Prim(VELZ)

        IntE = max(Prim(PRESSURE),pFloor)/(Model%gamma-1)
        KinE = 0.5*rho*(Vx**2.0+Vy**2.0+Vz**2.0)

        Con(DEN)    = rho
        Con(MOMX)   = rho*Vx
        Con(MOMY)   = rho*Vy
        Con(MOMZ)   = rho*Vz
        Con(ENERGY) = IntE+KinE
        
    end subroutine CellP2C

    !Set floors (density/pressure) from XML
    !Check gamera/sim/xFloor and gamera/floors/xFloor
    subroutine SetFloors(Model,iXML)
        type(Model_T)    , intent(in) :: Model
        type(XML_Input_T), intent(in) :: iXML

        call SetFloorsWDefs(Model,iXML,TINY,TINY)
    end subroutine SetFloors

    subroutine SetFloorsWDefs(Model,iXML,dFDef,pFDef)
        type(Model_T)    , intent(in) :: Model
        type(XML_Input_T), intent(in) :: iXML
        real(rp)         , intent(in) :: dFDef,pFDef

        if (iXML%Exists("sim/dFloor"   )) call iXML%Set_Val(dFloor,"sim/dFloor"   ,dFDef)
        if (iXML%Exists("floors/dFloor")) call iXML%Set_Val(dFloor,"floors/dFloor",dFDef)

        if (iXML%Exists("sim/pFloor"   )) call iXML%Set_Val(pFloor,"sim/pFloor"   ,pFDef)
        if (iXML%Exists("floors/pFloor")) call iXML%Set_Val(pFloor,"floors/pFloor",pFDef)

    end subroutine SetFloorsWDefs

    !Calculate semi-relativistic transform
    !m = Lam*(rho x u)
    subroutine Mom2Rel(Model,D,B,Lam)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: D, B(NDIM)
        real(rp), intent(out) :: Lam(NDIM,NDIM)

        real(rp), dimension(NDIM) :: bhat
        real(rp) :: alpha

        
        if (.not. Model%doBoris) then
            Lam = Eye33
            return
        endif
        !Calc ratio of magnetic/plasma mass, va2/ca2
        alpha = dot_product(B,B)/max(D,dFloor)/(Model%Ca**2.0)
        bhat = normVec(B)

        Lam = Eye33 + alpha*(Eye33 - Dyad(bhat,bhat))
    end subroutine Mom2Rel

    !Calculate inverse of semi-rel transform
    !rho x u = Laminv * m
    subroutine Rel2Mom(Model,D,B,Laminv)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: D, B(NDIM)
        real(rp), intent(out) :: Laminv(NDIM,NDIM)

        real(rp), dimension(NDIM) :: bhat
        real(rp) :: alpha,GamA

        if (.not. Model%doBoris) then
            Laminv = Eye33
            return
        endif
        !Calc ratio of magnetic/plasma mass, va2/ca2
        alpha = dot_product(B,B)/max(D,dFloor)/(Model%Ca**2.0)
        bhat = normVec(B)

        GamA = 1.0/sqrt(1.0 + alpha)
        LamInv = (GamA**2.0)*( Eye33 + alpha*Dyad(bhat,bhat))
    end subroutine Rel2mom
end module gamutils
