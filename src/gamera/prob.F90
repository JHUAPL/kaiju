
!Contains initial conditions for various base problems
module prob
    use gamtypes
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use background

    implicit none

    !Constants for all problems

    !EM unit constants
    real(rp), parameter :: emScl = 1
    
    contains

    !Set initState function pointer based on icStr (problem id)
    subroutine setIC_T(initState,icStr,userInitFunc)
        procedure(StateIC_T), pointer, intent(out) :: initState
        character(len=*), intent(in) :: icStr
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc

        initState => NULL()

        select case (icStr)

        case ("OT2D")
            initState => initOT2D
        case ("BW")
            initState => initBW
        case ("KH")
            initState => initKH
        case ("LOOP2D")
            initState => initLoop
        case ("SHEET")
            initState => initSheet
        case ("RINGBW") 
            initState => ringBW
        case ("RINGLOOP") 
            initState => ringLoop
        case ("SOD")
            initState => initSod
        case ("KH_McNally")
            initState => initKH_McNally
        case ("ROTOR")
            initState => initRotor
        case ("ADV1D")
            initState => initADV1D
        case ("ADV2D")
            initState => initADV2D
        case ("IMPLOSION") 
            initState => initImplosion
        case ("GEM") 
            initState => initGEM
        case ("ALFVEN") 
            initState => initALFVEN
        case ("MULTIFLUIDBW")
            initState => initMultiFBW
        case ("USER","user")
            initState => userInitFunc

        case default
            write(*,*) 'Unknown problem id, exiting ...'
            write(*,*) icStr
            stop
        end select

    end subroutine

    subroutine initOT2D(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real(rp) :: xc,yc,zc
        real(rp) :: Rho0,P0,B0
        real(rp) :: Vx,Vy,Vz, KinE
        integer :: i,j,k

        write(*,*) 'Initializing OT-2D'

        !Set OT defaults, then check input deck
        Rho0 = 25.0/(36.0*pi)
        P0 = 5/(12.0*pi)
        B0 = 1/sqrt(4.0*pi)

        call inpXML%Set_Val(Rho0,"prob/d0",Rho0)
        call inpXML%Set_Val(P0  ,"prob/P0",P0)
        call inpXML%Set_Val(B0  ,"prob/B0",B0)

        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)
                    
                    State%Gas(i,j,k,DEN,BLK) = Rho0
                    Vx = -sin(2*pi*yc)
                    Vy = sin(2*pi*xc)
                    Vz = 0.0
                    KinE = 0.5*Rho0*(Vx**2.0+Vy**2.0+Vz**2.0)
                    State%Gas(i,j,k,MOMX,BLK) = Vx*Rho0
                    State%Gas(i,j,k,MOMY,BLK) = Vy*Rho0
                    State%Gas(i,j,k,MOMZ,BLK) = Vz*Rho0
                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P0/(Model%gamma-1)
                    !Get faces  
                    State%magFlux(i,j,k,XDIR) = (-B0*sin(2*pi*yc))*Grid%Face(i,j,k,XDIR)
                    State%magFlux(i,j,k,YDIR) = (B0*sin(4*pi*xc))*Grid%Face(i,j,k,YDIR)
                    State%magFlux(i,j,k,ZDIR) = 0.0
                enddo
            enddo
        enddo
        
    end subroutine

    !Initialize field loop advection w/ singularity
    subroutine ringLoop(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        if (.not. Model%doRing) then
            write(*,*) 'Need ring-avg'
            stop
        endif
        call initLoop(Model,Grid,State,inpXML)

        allocate(cylindricalPoleBC_T    :: Grid%externalBCs(INI )%p)
        allocate(zeroGradientOuterIBC_T :: Grid%externalBCs(OUTI)%p)
        allocate(periodicInnerJBC_T     :: Grid%externalBCs(INJ )%p)
        allocate(periodicOuterJBC_T     :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T     :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T     :: Grid%externalBCs(OUTK)%p)

    end subroutine ringLoop

    !Initialize blast wave w/ cylindrical singularity
    subroutine ringBW(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        if (.not. Model%doRing) then
            write(*,*) 'Need ring-avg'
            stop
        endif

        call initBW(Model,Grid,State,inpXML)

        allocate(cylindricalPoleBC_T    :: Grid%externalBCs(INI )%p)
        allocate(zeroGradientOuterIBC_T :: Grid%externalBCs(OUTI)%p)
        
    end subroutine ringBW

    !Initialize blast wave w/ multifluid setup
    subroutine initMultiFBW(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real(rp) :: dScl,pScl

        if (.not. Model%doMultiF) then
            write(*,*) 'Need multifluid'
            stop
        endif

        call initBW(Model,Grid,State,inpXML)

        !Setup fluid species
        dScl = 0.8
        pScl = 0.2
        State%Gas(:,:,:,DEN   ,1) = dScl    *State%Gas(:,:,:,DEN   ,BLK)
        State%Gas(:,:,:,DEN   ,2) = (1-dScl)*State%Gas(:,:,:,DEN   ,BLK)
        State%Gas(:,:,:,ENERGY,1) = pScl    *State%Gas(:,:,:,ENERGY,BLK)
        State%Gas(:,:,:,ENERGY,2) = (1-pScl)*State%Gas(:,:,:,ENERGY,BLK)

    end subroutine initMultiFBW

    !Initialize all blast wave variants here
    !MHD/Hydro, 2D/3D
    subroutine initBW(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real(rp) :: xc,yc,zc,x0,y0,z0,R
        real(rp) :: B0_BW, P0, pRat, dRat, D0, rC
        real(rp) :: Vx0,Vy0,Vz0, KinE, IntE, Rho
        integer :: i,j,k
        logical :: doB0, doBz
        real(rp) :: bScl,bSclz

        procedure(VectorField_T), pointer :: Axyz
        procedure(GasIC_T), pointer :: Wxyz

        write(*,*) 'Initializing Blast Wave ...'

        !Get problem parameters from input deck
        call inpXML%Set_Val(x0 ,"prob/x0" ,0.0_rp)
        call inpXML%Set_Val(y0 ,"prob/y0" ,0.0_rp)
        call inpXML%Set_Val(z0 ,"prob/z0" ,0.0_rp)
        call inpXML%Set_Val(Vx0,"prob/Vx0",0.0_rp)
        call inpXML%Set_Val(Vy0,"prob/Vy0",0.0_rp)
        call inpXML%Set_Val(Vz0,"prob/Vz0",0.0_rp)

        call inpXML%Set_Val(B0_BW,"prob/B0",1.0_rp)
        B0_BW = emScl*B0_BW

        call inpXML%Set_Val(P0  ,"prob/P0"  ,0.1_rp  )
        call inpXML%Set_Val(pRat,"prob/pRat",100.0_rp)
        call inpXML%Set_Val(dRat,"prob/dRat",1.0_rp  )
        call inpXML%Set_Val(D0  ,"prob/D0"  ,1.0_rp  )
        call inpXML%Set_Val(rC  ,"prob/rC"  ,0.1_rp  )

        !Do Z field or not
        call inpXML%Set_Val(doBz,"prob/doBz",.false.)
        !Do field via background or raw
        call inpXML%Set_Val(doB0,"prob/doBack",.false.)

        if (doBz) then
            bScl  = B0_BW/sqrt(3.0)
            bSclz = bScl
        else
            bScl  = B0_BW/sqrt(2.0)
            bSclz = 0.0
        endif

        call inpXML%Set_Val(bSclz,"prob/Bz0",bSclz) !Allow problem file to overwrite Z component

        !Initialize State variable analytic function
        Wxyz => GasIC_BW
        call GasIC2State(Model,Grid,State,Wxyz)

        !Add B field if nec., use either background or vector potential
        if (Model%doMHD) then
            if (doB0) then
                !Incorporate field via background
                Model%doBackground = .true.
                Model%B0 => BlastB0
            else
                Axyz => VectorPot_BW
                call VectorPot2Flux(Model,Grid,State,Axyz)
            endif
        endif


        !Local functions for initBW  
        contains

            !Background field for blast wave
            subroutine BlastB0(x,y,z,Ax,Ay,Az)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: Ax,Ay,Az

                Ax = bScl
                Ay = bScl
                Az = bSclz

            end subroutine BlastB0

            subroutine VectorPot_BW(x,y,z,Ax,Ay,Az)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: Ax,Ay,Az
        
                Ax = 0.0
                Ay = bSclz*x
                Az = bScl*(y - x)
            end subroutine VectorPot_BW

            subroutine GasIC_BW(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P

                real(rp) :: R

                R = sqrt( (x-x0)**2.0 + (y-y0)**2.0 + (z-z0)**2.0)
                !R = sqrt( (x-x0)**2.0 + (y-y0)**2.0 )
                D = D0
                P = P0
                Vx = Vx0
                Vy = Vy0
                Vz = Vz0
                if (R <= rC) then
                    D = dRat*D
                    P = pRat*P
                endif
            end subroutine GasIC_BW
    end subroutine initBW

    subroutine initSheet(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        procedure(GasIC_T), pointer :: Wxyz

        integer :: i,j,k,d
        real(rp) :: A,beta
        real(rp), dimension(NDIM) :: Bxyz, fN, xFC
        real(rp) :: fA

        !Get problem parameters from input deck
        call inpXML%Set_Val(A   ,"prob/A"   ,0.1_rp)
        call inpXML%Set_Val(beta,"prob/beta",0.1_rp)

        Wxyz => GasIC
        call GasIC2State(Model,Grid,State,Wxyz)
        State%magFlux(:,:,:,:) = 0.0

        
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    do d=1,NDIM
                        fA = Grid%face(i,j,k,d)
                        fN = Grid%Tf(i,j,k,NORMX:NORMZ,d)
                        call faceCenter(Grid,i,j,k,xFC(XDIR),xFC(YDIR),xFC(ZDIR),d)
                        if (abs(xFC(XDIR))>0.25) then
                            Bxyz = [0,1,0]
                        else
                            Bxyz = [0,-1,0]
                        endif
                        State%magFlux(i,j,k,d) = fA*dot_product(fN,Bxyz)
                    enddo
                enddo
            enddo
        enddo

        !Local functions
        contains
            subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P

                D = 1.0
                P = 0.5*beta

                Vx = A*sin(2*PI*y)
                Vy = 0.0
                Vz = 0.0

            end subroutine GasIC
    end subroutine initSheet

    !Initialize hydro/MHD KH
    subroutine initKH(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real(rp) :: xc,yc,zc, Rho,Vx,Vy,KinE,Theta
        real(rp) :: yCrit, B0, P0, Amp, D0,D1,V0,V1

        integer :: i,j,k

        write(*,*) 'Initializing Kelvin-Helmholtz ...'
        call inpXML%Set_Val(Model%gamma,"physics/gamma",1.4_rp)
        call inpXML%Set_Val(yCrit,"prob/yCrit",0.25_rp)
        
        call inpXML%Set_Val(P0 ,"prob/P0" ,2.5_rp )
        call inpXML%Set_Val(Amp,"prob/Amp",0.01_rp)
        call inpXML%Set_Val(D0 ,"prob/D0" ,1.0_rp )
        call inpXML%Set_Val(D1 ,"prob/D1" ,2.0_rp )
        call inpXML%Set_Val(V0 ,"prob/V0" ,-0.5_rp)
        call inpXML%Set_Val(V1 ,"prob/V1" , 0.5_rp)
        call inpXML%Set_Val(B0,"prob/B0",0.2_rp)
        call inpXML%Set_Val(Theta,"prob/theta",0.0_rp)
        B0 = emScl*B0

        if (.not. Model%doMHD) then
            B0 = 0.0
        endif

        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)

                    if (abs(yc) <= yCrit) then
                        Rho = D1
                        Vx  = V1
                    else
                        Rho = D0
                        Vx  = V0
                    endif
                    Vy = 0.0

                    !Add perturbations
                    Vx = Vx + genRand(-Amp,Amp)
                    Vy = Vy + genRand(-Amp,Amp)

                    !Add monochromatic perturbation
                    !Vy = Amp_KH*sin(2*pi*xc)

                    State%Gas(i,j,k,DEN,BLK) = Rho
                    State%Gas(i,j,k,MOMX,BLK) = Rho*Vx
                    State%Gas(i,j,k,MOMY,BLK) = Rho*Vy
                    State%Gas(i,j,k,MOMZ,BLK) = 0.0

                    KinE = 0.5*Rho*(Vx**2.0+Vy**2.0)

                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P0/(Model%gamma-1)

                    !Initialize fields
                    State%magFlux(i,j,k,IDIR) = B0*cos(Theta*pi/180.0)*Grid%Face(i,j,k,IDIR)
                    State%magFlux(i,j,k,JDIR) = 0.0
                    State%magFlux(i,j,k,KDIR) = B0*sin(Theta*pi/180.0)*Grid%Face(i,j,k,KDIR)
                enddo
            enddo
        enddo        

        call WipeBCs(Model,Grid)
        allocate(periodicInnerIBC_T      :: Grid%externalBCs(INI )%p)
        allocate(periodicOuterIBC_T      :: Grid%externalBCs(OUTI)%p)
        allocate(periodicInnerJBC_T      :: Grid%externalBCs(INJ )%p)
        allocate(periodicOuterJBC_T      :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T      :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T      :: Grid%externalBCs(OUTK)%p)

    end subroutine initKH

    !Initialize field loop advection
    subroutine initLoop(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        procedure(VectorField_T), pointer :: Axyz
        real(rp) :: xc,yc,zc, Rho,Vx,Vy,Vz,Vz0,P,KinE
        real(rp) :: Rho0,P0,V0,alpha
        real(rp) :: A0_FL, R0_FL, x0_FL, y0_FL

        integer :: i,j,k

        write(*,*) 'Initializing Field Loop Advection ...'
        call inpXML%Set_Val(Rho0,"prob/d0" ,1.0_rp)
        call inpXML%Set_Val(P0  ,"prob/P0" ,1.0_rp)
        call inpXML%Set_Val(V0  ,"prob/V0" ,1.0_rp)
        call inpXML%Set_Val(Vz0 ,"prob/Vz0",0.0_rp)

        call inpXML%Set_Val(alpha,"prob/alpha",60.0_rp)
        call inpXML%Set_Val(A0_FL,"prob/A0"   ,1.0e-3_rp)
        call inpXML%Set_Val(R0_FL,"prob/R0"   ,0.3_rp)
        call inpXML%Set_Val(x0_FL,"prob/x0"   ,0.0_rp)
        call inpXML%Set_Val(y0_FL,"prob/y0"   ,0.0_rp)

        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)
                    Rho = Rho0
                    P = P0
                    Vx = V0*sin(alpha*pi/180.0)
                    Vy = V0*cos(alpha*pi/180.0)
                    Vz = Vz0
                    State%Gas(i,j,k,DEN,BLK) = Rho
                    State%Gas(i,j,k,MOMX,BLK) = Rho*Vx
                    State%Gas(i,j,k,MOMY,BLK) = Rho*Vy
                    State%Gas(i,j,k,MOMZ,BLK) = Rho*Vz

                    KinE = 0.5*Rho*(Vx**2.0+Vy**2.0+Vz**2.0)

                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P/(Model%gamma-1)

                enddo
            enddo
        enddo
        !Initialize fields
        Axyz => VectorPot_Loop2D
        call VectorPot2Flux(Model,Grid,State,Axyz)

        !Local functions
        contains
            subroutine VectorPot_Loop2D(x,y,z,Ax,Ay,Az)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: Ax,Ay,Az
        
                real(rp) :: r
                Ax = 0.0
                Ay = 0.0
                r = sqrt( (x-x0_FL)**2.0 + (y-y0_FL)**2.0 ) ! + z**2.0)
                Az = max( 0.0 , A0_FL*(R0_FL-r) )
        
            end subroutine VectorPot_Loop2D
    end subroutine initLoop

    ! This is the APJ KH test problem by NcNally et al. [2012], doi:10.1088/0067-0049/201/2/18
    ! The simulation grid is fixed to be [0 1]x[0 1] 
    subroutine initKH_McNally(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real (rp) :: xc,yc,zc, Rho,Vx,Vy,KinE,Theta
        real (rp) :: yCrit, B0, P0, Amp, rho1, rho2, rhom, u1, u2, um, L

        integer :: i,j,k

        write(*,*) 'Initializing McNally Kelvin-Helmholtz Test ...'
        call inpXML%Set_Val(Model%gamma,"physics/gamma",1.4_rp)
        call inpXML%Set_Val(yCrit,"prob/yCrit",0.25_rp)

        call inpXML%Set_Val(P0 ,"prob/P0" ,2.5_rp )
        call inpXML%Set_Val(Amp,"prob/Amp",0.01_rp)
        call inpXML%Set_Val(rho1 ,"prob/D0" ,1.0_rp )
        call inpXML%Set_Val(rho2 ,"prob/D1" ,2.0_rp )
        call inpXML%Set_Val(u1 ,"prob/V0" , 0.5_rp)
        call inpXML%Set_Val(u2 ,"prob/V1" ,-0.5_rp)
        call inpXML%Set_Val(L,"prob/L",0.025_rp)

        call inpXML%Set_Val(B0,"prob/B0",0.2_rp)
        call inpXML%Set_Val(Theta,"prob/theta",0.0_rp)
        B0 = emScl*B0

        rhom = (rho1-rho2)/2
        um = (u1-u2)/2

        if (.not. Model%doMHD) then
            B0 = 0.0
        endif

        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)

                    if (yc < 1.0/4.0) then
                        Rho = rho1 - rhom*exp((yc-1.0/4.0)/L)                        
                        Vx  = u1   -   um*exp((yc-1.0/4.0)/L)
                    else if ((yc >= 1.0/4.0).and.(yc < 2.0/4.0)) then
                        Rho = rho2 + rhom*exp((-yc+1.0/4.0)/L)
                        Vx  = u2   +   um*exp((-yc+1.0/4.0)/L)
                    else if ((yc >= 2.0/4.0).and.(yc < 3.0/4.0)) then
                        Rho = rho2 + rhom*exp((yc-3.0/4.0)/L)
                        Vx  = u2   +   um*exp((yc-3.0/4.0)/L)
                    else
                        Rho = rho1 - rhom*exp(-(yc-3.0/4.0)/L)
                        Vx  = u1   -   um*exp(-(yc-3.0/4.0)/L)
                    end if

                    Vy = Amp*sin(4*pi*xc)

                    !Add perturbations
                    !Vx = Vx + genRand(-Amp,Amp)
                    !Vy = Vy + genRand(-Amp,Amp)

                    !Add monochromatic perturbation
                    !Vy = Amp_KH*sin(2*pi*xc)

                    State%Gas(i,j,k,DEN,BLK) = Rho
                    State%Gas(i,j,k,MOMX,BLK) = Rho*Vx
                    State%Gas(i,j,k,MOMY,BLK) = Rho*Vy
                    State%Gas(i,j,k,MOMZ,BLK) = 0.0

                    KinE = 0.5*Rho*(Vx**2.0+Vy**2.0)

                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P0/(Model%gamma-1)

                    !Initialize fields
                    State%magFlux(i,j,k,IDIR) = B0*cos(Theta*pi/180.0)*Grid%Face(i,j,k,IDIR)
                    State%magFlux(i,j,k,JDIR) = 0.0
                    State%magFlux(i,j,k,KDIR) = B0*sin(Theta*pi/180.0)*Grid%Face(i,j,k,KDIR)
                enddo
            enddo
        enddo        

        call WipeBCs(Model,Grid)
        allocate(periodicInnerIBC_T      :: Grid%externalBCs(INI )%p)
        allocate(periodicOuterIBC_T      :: Grid%externalBCs(OUTI)%p)
        allocate(periodicInnerJBC_T      :: Grid%externalBCs(INJ )%p)
        allocate(periodicOuterJBC_T      :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T      :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T      :: Grid%externalBCs(OUTK)%p)

    end subroutine initKH_McNally

    subroutine initRotor(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        procedure(VectorField_T), pointer :: Axyz
        real (rp) :: xc,yc,zc, Rho,Vx,Vy,P,KinE
        real (rp) :: yCrit, B0, P0, r

        integer :: i,j,k

        write(*,*) 'Initializing MHD Rotor Test ...'
        call inpXML%Set_Val(Model%gamma,"physics/gamma",1.4_rp)
        call inpXML%Set_Val(B0,"prob/B0",1.41_rp)
        B0 = emScl*B0
        call inpXML%Set_Val(P0,"prob/P0",0.5_rp)

        if (.not. Model%doMHD) then
            B0 = 0.0
        endif

        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)

                    P = P0
                    r = sqrt((xc-0.5)**2+(yc-0.5)**2)
                    if (r<0.1) then
                        Rho = 10.0
                        Vx = -(10.0*yc - 5.0)
                        Vy =   10.0*xc - 5.0
                    elseif (r<0.115) then
                        Rho = 1.0 + 9.0*(23.0-200.0*r)/3.0
                        Vx  = -(10.0*yc-5.0)*(23.0-200.0*r)/3.0
                        Vy  =  (10.0*yc-5.0)*(23.0-200.0*r)/3.0
                    else
                        Rho = 1.0
                        Vx = 0.0
                        Vy = 0.0
                    endif

                    State%Gas(i,j,k,DEN,BLK) = Rho
                    State%Gas(i,j,k,MOMX,BLK) = Rho*Vx
                    State%Gas(i,j,k,MOMY,BLK) = Rho*Vy
                    State%Gas(i,j,k,MOMZ,BLK) = 0.0

                    KinE = 0.5*Rho*(Vx**2.0+Vy**2.0)

                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P/(Model%gamma-1)

                enddo
            enddo
        enddo  

        !Initialize fields
        Axyz => VectorPot_Loop2D
        call VectorPot2Flux(Model,Grid,State,Axyz)

        !Local functions
        contains
            subroutine VectorPot_Loop2D(x,y,z,Ax,Ay,Az)
                real (rp), intent(in) :: x,y,z
                real (rp), intent(out) :: Ax,Ay,Az
        
                real (rp) :: r
                Ax = 0.0
                Ay = 0.0
                Az = x*B0
        
            end subroutine VectorPot_Loop2D


    end subroutine initRotor

    subroutine initADV1D(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real (rp) :: xc,yc,zc, Rho,Rho0,Rho1,Rho2,Rho3,Rho4,Vx,Vy,P,KinE
        real (rp) :: P0,B0
        integer :: i,j,k

        write(*,*) 'Initializing Advection Test (Four Shapes) ...'
        call inpXML%Set_Val(Model%gamma,"physics/gamma",1.4_rp)
        call inpXML%Set_Val(P0 ,"prob/P0" ,1.0_rp )

        if (.not. Model%doMHD) then
            B0 = 0.0
        endif

        ! This is a standard advection test problem defined at [-1 1], with advection speed of 1
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)

                    Rho0 = 0.1
                    Rho1 = exp(-((xc+0.7)/0.05)**2)
                    if (abs(xc+0.3) <= 0.1) then
                        Rho2 = 1.0
                    else
                        Rho2 = 0.0
                    endif
                    if (abs(xc-0.1)<=0.1) then
                        Rho3 = 1.0-abs(10.0*(xc-0.1))
                    else
                        Rho3 = 0.0
                    endif
                    if ((0.1**2 - (xc-0.5)**2) >=0) then
                        Rho4 = 10.0*sqrt(0.1**2 - (xc-0.5)**2)
                    else
                        Rho4 = 0.0
                    endif

                    Rho = Rho0 + Rho1 + Rho2 + Rho3 + Rho4
                    Vx = 1.0
                    Vy = 0.0
                    P = 1.0

                    State%Gas(i,j,k,DEN,BLK) = Rho
                    State%Gas(i,j,k,MOMX,BLK) = Rho*Vx
                    State%Gas(i,j,k,MOMY,BLK) = Rho*Vy
                    State%Gas(i,j,k,MOMZ,BLK) = 0.0

                    KinE = 0.5*Rho*(Vx**2.0+Vy**2.0)

                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P/(Model%gamma-1)

                    !Initialize fluxes, this is a hydro advection problem so all the fiels are zero
                    State%magFlux(i,j,k,IDIR) = 0.0
                    State%magFlux(i,j,k,JDIR) = 0.0
                    State%magFlux(i,j,k,KDIR) = 0.0
                enddo
            enddo
        enddo  

    end subroutine initADV1D

    !Initialize field loop advection
    subroutine initADV2D(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        procedure(VectorField_T), pointer :: Axyz
        real (rp) :: xc,yc,zc, r, Rho,Vx,Vy,P,KinE
        real (rp) :: Rho0,P0,V0,alpha
        real (rp) :: A0_FL, R0_FL, x0_FL, y0_FL

        integer :: i,j,k

        write(*,*) 'Initializing Field Loop Advection ...'
        call inpXML%Set_Val(Rho0 ,"prob/d0",1.0_rp)
        call inpXML%Set_Val(P0   ,"prob/P0",1.0_rp)
        call inpXML%Set_Val(V0   ,"prob/V0",1.0_rp)
        call inpXML%Set_Val(alpha,"prob/alpha",60.0_rp)

        call inpXML%Set_Val(A0_FL,"prob/A0",1.0e-3_rp)
        call inpXML%Set_Val(R0_FL,"prob/R0",0.3_rp)
        call inpXML%Set_Val(x0_FL,"prob/x0",0.0_rp)
        call inpXML%Set_Val(y0_FL,"prob/y0",0.0_rp)

        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)
                    Rho = Rho0
                    P = P0
                    Vx = V0*sin(alpha*pi/180.0)
                    Vy = V0*cos(alpha*pi/180.0)

                    r = sqrt((xc-x0_FL)**2 + (yc-y0_FL)**2)

                    if (r<R0_FL) then
                        Rho = Rho0*2.0
                    endif

                    State%Gas(i,j,k,DEN,BLK) = Rho
                    State%Gas(i,j,k,MOMX,BLK) = Rho*Vx
                    State%Gas(i,j,k,MOMY,BLK) = Rho*Vy
                    State%Gas(i,j,k,MOMZ,BLK) = 0.0

                    KinE = 0.5*Rho*(Vx**2.0+Vy**2.0)

                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P/(Model%gamma-1)

                    !Initialize fluxes, this is a hydro advection problem so all the fiels are zero
                    State%magFlux(i,j,k,IDIR) = 0.0
                    State%magFlux(i,j,k,JDIR) = 0.0
                    State%magFlux(i,j,k,KDIR) = 0.0

                enddo
            enddo
        enddo

    end subroutine initADV2D

    !Initialize field loop advection
    subroutine initSod(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real (rp) :: xc,yc,zc, r, Rho,Vx,Vy,P,KinE

        integer :: i,j,k

        write(*,*) 'Initializing SOD ...'

!% rho(xc<0)=1;
!% rho(xc>0)=0.125;
!% p(xc<0)=1;
!% p(xc>0)=0.1;
!% bi(:)=0.75;
!% bj(xj<0)=1;
!% bj(xj>0)=-1;

        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)

                    if (xc<0.0) then
                        Rho = 1.0
                        P   = 1.0
                    else
                        Rho = 0.125
                        P   = 0.1
                    endif

                    State%Gas(i,j,k,DEN,BLK) = Rho
                    State%Gas(i,j,k,MOMX,BLK) = Rho*Vx
                    State%Gas(i,j,k,MOMY,BLK) = Rho*Vy
                    State%Gas(i,j,k,MOMZ,BLK) = 0.0

                    KinE = 0.5*Rho*(Vx**2.0+Vy**2.0)

                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P/(Model%gamma-1)

                    !Initialize fluxes, this is a hydro advection problem so all the fiels are zero
                    State%magFlux(i,j,k,IDIR) = 0.0
                    State%magFlux(i,j,k,JDIR) = 0.0
                    State%magFlux(i,j,k,KDIR) = 0.0

                enddo
            enddo
        enddo

        allocate(zeroGradientInnerIBC_T :: Grid%externalBCs(INI )%p)
        allocate(zeroGradientOuterIBC_T :: Grid%externalBCs(OUTI)%p)

    end subroutine initSod

    ! Liska, R., & Wendroff, B., "Comparison of Several difference schemes on 1D and 2D Test problems for the Euler equations", 
    ! The test is described in section 4.7 of LW, although it was presented in an earlier paper by Hui et al. (JCP, 153, 596, 1999).
    ! the domain is square: [0 0.3]x[0 0.3], for (xc+yc)<0.15, rho = p = 1; otherwise rho = 0.125, p = 0.14.
    ! reflectiing boundaries in both x and y directions;
    subroutine initImplosion(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real (rp) :: xc,yc,zc, Rho,Vx,Vy,P,KinE
        real (rp) :: yCrit, B0, P0, Amp, rho1, rho2, rhom, u1, u2, um, L

        integer :: i,j,k

        write(*,*) 'Initializing Implosion Test ...'
        call inpXML%Set_Val(Model%gamma,"physics/gamma",1.4_rp)

        if (.not. Model%doMHD) then
            B0 = 0.0
        endif

        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !Find centers
                    call cellCenter(Grid,i,j,k,xc,yc,zc)

                        Vx  = 0.0
                        Vy  = 0.0

                    if ((xc+yc) > 0.15) then
                        Rho = 1.0                        
                        P   = 1.0
                    else 
                        Rho = 0.125
                        P = 0.14
                    endif

                    State%Gas(i,j,k,DEN,BLK) = Rho
                    State%Gas(i,j,k,MOMX,BLK) = Rho*Vx
                    State%Gas(i,j,k,MOMY,BLK) = Rho*Vy
                    State%Gas(i,j,k,MOMZ,BLK) = 0.0

                    KinE = 0.5*Rho*(Vx**2.0+Vy**2.0)

                    State%Gas(i,j,k,ENERGY,BLK) = KinE + P/(Model%gamma-1)

                    !Initialize fields
                    State%magFlux(i,j,k,IDIR) = 0.0
                    State%magFlux(i,j,k,JDIR) = 0.0
                    State%magFlux(i,j,k,KDIR) = 0.0
                enddo
            enddo
        enddo   

        allocate(cartesianReflectingInnerIBC_T :: Grid%externalBCs(INI )%p)
        allocate(cartesianReflectingOuterIBC_T :: Grid%externalBCs(OUTI)%p)
        allocate(cartesianReflectingInnerJBC_T :: Grid%externalBCs(INJ )%p)
        allocate(cartesianReflectingOuterJBC_T :: Grid%externalBCs(OUTJ)%p)
        allocate(zeroGradientInnerKBC_T        :: Grid%externalBCs(INK )%p)
        allocate(zeroGradientOuterKBC_T        :: Grid%externalBCs(OUTK)%p)

    end subroutine initImplosion

!------------------------------
!GEM reconnection problem 
    subroutine initGEM(Model,Grid,State,inpXML)

        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real(rp) :: xc,yc,zc,lam
        real(rp) :: B0,Lx,Ly,Psi0,n0,n_inf
        logical :: doHarris
        integer :: i,j,k
        procedure(VectorField_T), pointer :: Axyz
        procedure(GasIC_T), pointer :: Wxyz

        write(*,*) 'Initializing  GEM ...'

        !Define things up here once
        Lx = 25.6
        Ly = 12.8
        Psi0 = 0.1

        lam = 1.
        B0 = 1.
        n0 = 1.0
        n_inf = 0.2

        !Get problem parameters from input deck
        call inpXML%Set_Val(doHarris,"prob/doHarris",.false.)
        
        State%magFlux = 0.0

        !Initialize State variable analytic function
        Wxyz => GasIC_GEM
        call GasIC2State(Model,Grid,State,Wxyz)

        !Add B flux - first do vector potential for the perturbation island
        !Axyz => VectorPot_GEM
        !call VectorPot2Flux(Model,Grid,State,Axyz)

        ! then add the i-face Harris flux to the island
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%is,Grid%ie+1
                    !add face magFluxes bi(if_act,:,:) = bi(if_act,:,:)+tanh(yi(if_act,:,:)./lam);
                    call cellCenter(Grid,i,j,k,xc,yc,zc)
                    ! here the field should be evaluated at i-face centers, since it's along the i-dir, so yi = yc
                    if (doHarris) then
                       State%magFlux(i,j,k,IDIR) = Bx_GEM(Grid%xfc(i,j,k,YDIR,IDIR))*Grid%Face(i,j,k,IDIR)
                    else
                       State%magFlux(i,j,k,IDIR) = State%magFlux(i,j,k,IDIR)+tanh(yc/lam)*Grid%Face(i,j,k,IDIR)
                    endif
                enddo
            enddo
        enddo        

        call WipeBCs(Model,Grid)
        allocate(periodicInnerIBC_T      :: Grid%externalBCs(INI )%p)
        allocate(periodicOuterIBC_T      :: Grid%externalBCs(OUTI)%p)
        allocate(zeroGradientInnerJBC_T  :: Grid%externalBCs(INJ )%p)
        allocate(zeroGradientOuterJBC_T  :: Grid%externalBCs(OUTJ)%p)
        allocate(zeroExtensionInnerKBC_T :: Grid%externalBCs(INK )%p)
        allocate(zeroExtensionOuterKBC_T :: Grid%externalBCs(OUTK)%p)

        !Local functions for initBW  
        contains
            subroutine VectorPot_GEM(x,y,z,Ax,Ay,Az)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: Ax,Ay,Az


                Ax = 0.0
                Ay = 0.0
                Az = -Psi0*cos(2*pi*x/Lx)*cos(pi*y/Ly);

            end subroutine VectorPot_GEM

            function Bx_GEM(y)
                real(rp), intent(in) :: y
                real(rp) :: Bx_GEM

                Bx_GEM = B0*tanh(y/lam)

            end function Bx_GEM

            subroutine GasIC_GEM(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P

                real(rp) :: k0

                k0 = 0.5*(B0*B0)/n0

                D = n_inf+n0/cosh(y/lam)**2.0
                !P = 0.5*( n_inf+n0/cosh(y/lam)**2.0 )
                !P = 0.5*(n0+n_inf) - tanh(y/lam)**2./2
                P = k0*D
                Vx = 0.0
                Vy = 0.0
                Vz = 0.0
                
            end subroutine GasIC_GEM
    end subroutine initGEM

    subroutine initALFVEN(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        real(rp) :: xc,yc,zc
        real(rp) :: Rho0,P0,B_para,V_para,alp
        real(rp) :: Vx,Vy,Vz, KinE
        integer :: i,j,k
        procedure(VectorField_T), pointer :: Axyz
        procedure(GasIC_T), pointer :: Wxyz

        write(*,*) 'Initializing Nonlinear Circular Alfven'

        !Set ALFVEN defaults, then check input deck
        Rho0 = 1.0
        P0 = 0.1
        B_para = 1.0
        V_para = 0.0
        alp = pi/3.0
        call inpXML%Set_Val(Rho0,"prob/d0",Rho0)
        call inpXML%Set_Val(P0  ,"prob/P0",P0)
        call inpXML%Set_Val(B_para,"prob/B_para",B_para)
        call inpXML%Set_Val(V_para,"prob/V_para",V_para)
        call inpXML%Set_Val(alp,"prob/alp",alp)

        !Initiate gas state variables
        Wxyz => GasIC_ALF
        call GasIC2State(Model, Grid, State, Wxyz)
        !Initiate MagFlux using vector potential functions
        Axyz => VectorPot_ALF
        call VectorPot2Flux(Model,Grid,State,Axyz)

        !Local functions for initALFVEN
        contains
            subroutine VectorPot_ALF(x,y,z,Ax,Ay,Az)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: Ax,Ay,Az

                Ax =      B_para*sin(alp)*z
                Ay = -1.0*B_para*cos(alp)*z+0.1*sin(2*pi*(x*cos(alp)+y*sin(alp)))/(2*pi*cos(alp))
                Az =                        0.1*cos(2*pi*(x*cos(alp)+y*sin(alp)))/(2*pi)
            end subroutine VectorPot_ALF

            subroutine GasIC_ALF(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P

                D = Rho0
                P = P0
                Vx = V_para*cos(alp) - 0.1*sin(2*pi*(x*cos(alp)+y*sin(alp)))*sin(alp)
                Vy = V_para*sin(alp) + 0.1*sin(2*pi*(x*cos(alp)+y*sin(alp)))*cos(alp)
                Vz =                   0.1*cos(2*pi*(x*cos(alp)+y*sin(alp)))

            end subroutine GasIC_ALF
    end subroutine initALFVEN
end module prob
