!Setup for Neptune magnetosphere (no MI coupling)

module usergamic
    use gamtypes
    use volttypes, only: IONO_EX,IONO_EY,IONO_EZ
    use gamutils
    use xml_input
    use bcs
    use msphutils
    use wind
    
    implicit none

    !Various global would go here
    real(rp) :: Rho0,P0,Om0,RhoX
    real(rp) :: ThD,CosThD,SinThD
    
    integer :: iSh = 2

    !These values are taken from msphutils
    !integer, parameter :: PsiShells = 5
    !integer, parameter :: PsiSt = -3

    ! type for remix BC
    type, extends(innerIBC_T) :: CoroInnerBC_T

        !Main electric field structures
        real(rp), allocatable, dimension(:,:,:,:) :: inEijk,inExyz

        contains

        procedure :: doInit => InitCoroInner
        procedure :: doBC => CoroInner

    end type CoroInnerBC_T

    contains

    subroutine initUser(Model,Grid,State,inpXML)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        character(len=strLen) :: pID !Planet ID string

        pID = "Neptune"
    !Get defaults from input deck
        call inpXML%Set_Val(Rho0 ,"prob/Rho0",0.0001_rp)
        call inpXML%Set_Val(RhoX ,"prob/RhoX",10.0_rp)
        call inpXML%Set_Val(P0   ,"prob/P0"  ,0.001_rp)
        call inpXML%Set_Val(ThD  ,"prob/ThD" ,46.0_rp)

        CosThD = cos(ThD*deg2rad)
        SinThD = sin(ThD*deg2rad)

        !Set magnetosphere parameters
        call setMagsphere(Model,inpXML,pID)

        Om0 = 2*PI*gT0/(16.0*60*60)

    !Set BCs
        ! deallocate default BCs
        deallocate(Grid%ExternalBCs(INI )%p)
        deallocate(Grid%ExternalBCs(OUTI)%p)
        deallocate(Grid%ExternalBCs(INJ )%p)
        deallocate(Grid%ExternalBCs(OUTJ)%p)
        deallocate(Grid%ExternalBCs(INK )%p)
        deallocate(Grid%ExternalBCs(OUTK)%p)

        !Set BCs (spherical, RPT)
        allocate(CoroInnerBC_T      :: Grid%externalBCs(INI )%p)
        allocate(WindBC_T           :: Grid%externalBCs(OUTI)%p)
        allocate(lfmInBC_T          :: Grid%externalBCs(INJ )%p)
        allocate(lfmOutBC_T         :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(OUTK)%p)
    !Set ICs
        
        call GasIC2State(Model,Grid,State,GasIC)
        call genTiltDipole(Model,Grid,State)

    !Set bounds
        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

    !Set hack functions
        if ( (Model%Ri == Model%NumRi) .or. (Model%Ri == 1) ) then
            Model%HackE => EFix
        endif

        if (Model%Ri == 1) then
            Model%HackFlux => InnerFlux
        endif

        Model%HackStep => PerStep

        contains

        subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
            real(rp), intent(in) :: x,y,z
            real(rp), intent(out) :: D,Vx,Vy,Vz,P

            real(rp) :: r,M
            r = sqrt(x**2.0+y**2.0+z**2.0)
            M = RampDown(r,5.0_rp,10.0_rp)
            
            D = M*RhoX*Rho0 + Rho0
            P = P0
            Vx = 0.0
            Vy = 0.0
            Vz = 0.0

        end subroutine GasIC

    end subroutine initUser

    !Create a cut dipole
    subroutine genTiltDipole(Model,Grid,State)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State

        real(rp), dimension(:,:,:,:), allocatable :: bFlux0, bFluxT
        procedure(VectorField_T), pointer :: Axyz

        !Generate some holders
        call allocGridVec(Model,Grid,bFlux0,.true.,NDIM)
        call allocGridVec(Model,Grid,bFluxT,.true.,NDIM)
        bFlux0 = 0.0
        bFluxT = 0.0

        !Calculate flux from B0_xyz components
        Model%doBackground = .true.

        !Put z-axis dipole into background
        Model%B0 => Dipole
        Axyz     => Dipole
        call AddB0(Model,Grid,Model%B0)
        call VectorField2Flux(Model,Grid,State,Axyz)
        !Save flux from z-axis dipole
        bFlux0(:,:,:,:) = State%magFlux(:,:,:,:) !bFlux0 = B0

        !Now use vector potential to set total initial field
        Axyz => TiltDipoleVP
        call VectorPot2Flux(Model,Grid,State,Axyz)
        bFluxT(:,:,:,:) = State%magFlux(:,:,:,:)

        !Subtract off background face-fluxes to get dB
        State%magFlux = bFluxT - bFlux0

    end subroutine genTiltDipole

    subroutine TiltDipoleVP(x,y,z,Ax,Ay,Az)
        
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Ax,Ay,Az

        real(rp), dimension(NDIM) :: A,m,r,rhat

        m = [SinThD,0.0_rp,CosThD]

        r = [x,y,z]
        rhat = r/norm2(r)

        A = M0*cross(m,rhat)/(dot_product(r,r))
        Ax = A(XDIR)
        Ay = A(YDIR)
        Az = A(ZDIR)
    end subroutine TiltDipoleVP

!---
    subroutine TiltDipole(x,y,z,t,Ax,Ay,Az)
        
        real(rp), intent(in) :: x,y,z,t
        real(rp), intent(out) :: Ax,Ay,Az

        real(rp), dimension(NDIM) :: m,r,A
        real(rp) :: rad

        r = [x,y,z]
        rad = norm2(r)
        m = [SinThD*cos(Om0*t),SinThD*sin(Om0*t),CosThD]
        

        A = 3*dot_product(m,r)*r/rad**5.0 - m/rad**3.0

        Ax = A(XDIR)
        Ay = A(YDIR)
        Az = A(ZDIR)

    end subroutine TiltDipole
!---
    !Fixes electric field before application
    subroutine EFix(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        type(State_T), intent(inout) :: State

        !Fix inner shells
        SELECT type(iiBC=>Gr%externalBCs(INI)%p)
            TYPE IS (CoroInnerBC_T)
                
                if (Model%Ri == 1) then
                    call IonEFix(Model,Gr,State,iiBC%inEijk)
                    

                endif
        END SELECT

        !Fix outer shells
        SELECT type(pWind=>Gr%externalBCs(OUTI)%p)
            TYPE IS (WindBC_T)
                if (Model%Ri == Model%NumRi) then
                    call WindEFix(pWind,Model,Gr,State)
                endif
        END SELECT

    end subroutine EFix

!---
    !Ensure no flux through degenerate faces
    subroutine InnerFlux(Model,Gr,gFlx,mFlx)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(inout) :: gFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(inout), optional :: mFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM,1:NDIM)

        integer :: i
        if (Model%Ri == 1) then
            !gFlx(Gr%is,:,:,DEN,IDIR,BLK) = 0.0
            gFlx(Gr%is:Gr%is+iSh,:,:,1:NVAR,1:NDIM,:) = 0.0
            if (present(mFlx)) then
                mFlx(Gr%is:Gr%is+iSh,:,:,1:NDIM,1:NDIM) = 0.0
            endif
        endif

    end subroutine InnerFlux
!---
    subroutine CoroEFix(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State

        integer :: i,j,k,iG
        real(rp), dimension(NDIM) :: x0,xE,ehat,Exyz

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,iG,x0,xE,ehat,Exyz) 
        do i=Gr%is,Gr%is+iSh+1
            do k=Gr%ks,Gr%keg-1
                do j=Gr%js,Gr%jeg-1

                    x0 = Gr%xyz(i,j,k,:)

                    !IDIR
                    xE = Gr%xyz(i+1,j,k,:)
                    ehat = normVec(xE-x0)
                    Exyz = CorotationE(0.5*(x0+xE),Model%t)
                    State%Efld(i,j,k,IDIR) = Gr%edge(i,j,k,IDIR)*dot_product(ehat,Exyz)

                    !JDIR
                    xE = Gr%xyz(i,j+1,k,:)
                    ehat = normVec(xE-x0)
                    Exyz = CorotationE(0.5*(x0+xE),Model%t)
                    State%Efld(i,j,k,JDIR) = Gr%edge(i,j,k,JDIR)*dot_product(ehat,Exyz)

                    !KDIR
                    xE = Gr%xyz(i,j,k+1,:)
                    ehat = normVec(xE-x0)
                    Exyz = CorotationE(0.5*(x0+xE),Model%t)
                    State%Efld(i,j,k,KDIR) = Gr%edge(i,j,k,KDIR)*dot_product(ehat,Exyz)
                enddo
            enddo
        enddo

    end subroutine CoroEFix

    subroutine PerStep(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        type(State_T), intent(inout) :: State

        !Fix inner shells
        SELECT type(iiBC=>Gr%externalBCs(INI)%p)
            TYPE IS (CoroInnerBC_T)
                if (Model%Ri == 1) call SetInnerFields(Model,Gr,iiBC%inEijk,iiBC%inExyz)
        END SELECT  

        
    end subroutine PerStep

    subroutine SetInnerFields(Model,Gr,inEijk,inExyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        real(rp), intent(inout) :: inEijk(1:PsiSh+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,1:NDIM)
        real(rp), intent(inout) :: inExyz(1:PsiSh  ,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:Model%nvInExyz)

        integer :: i,j,k,iG
        real(rp), dimension(NDIM) :: x0,xE,ehat,Exyz
        !Calculate electric field

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,iG,x0,xE,ehat,Exyz) 
        do i=1,PsiSh+1
            do k=Gr%ks,Gr%keg-1
                do j=Gr%js,Gr%jeg-1

                    iG = i+PsiSt-1 !Global i index for Grid access
                    x0 = Gr%xyz(iG,j,k,:)

                    !IDIR
                    xE = Gr%xyz(iG+1,j,k,:)
                    ehat = normVec(xE-x0)
                    Exyz = CorotationE(0.5*(x0+xE),Model%t)
                    inEijk(i,j,k,IDIR) = dot_product(ehat,Exyz)

                    !JDIR
                    xE = Gr%xyz(iG,j+1,k,:)
                    ehat = normVec(xE-x0)
                    Exyz = CorotationE(0.5*(x0+xE),Model%t)
                    inEijk(i,j,k,JDIR) = dot_product(ehat,Exyz)

                    !KDIR
                    xE = Gr%xyz(iG,j,k+1,:)
                    ehat = normVec(xE-x0)
                    Exyz = CorotationE(0.5*(x0+xE),Model%t)
                    inEijk(i,j,k,JDIR) = dot_product(ehat,Exyz)
                enddo
            enddo
        enddo

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,iG,x0)        
        do i=1,PsiSh
            do k=Gr%ksg,Gr%keg
                do j=Gr%jsg,Gr%jeg
                    iG = i+PsiSt-1 !Global i index into grid arrays
                    x0 = Gr%xyzcc(iG,j,k,:)
                    Exyz = CorotationE(x0,Model%t)
                    inExyz(iG,j,k,IONO_EX) = Exyz(IONO_EX)
                    inExyz(iG,j,k,IONO_EY) = Exyz(IONO_EY)
                    inExyz(iG,j,k,IONO_EZ) = Exyz(IONO_EZ)
                enddo
            enddo
        enddo

    end subroutine SetInnerFields

    function CorotationE(xyz,t) result(Exyz)
        real(rp), intent(in) :: xyz(NDIM), t
        real(rp) :: Exyz(NDIM)

        real(rp) :: x,y,z,phi,Xp,r,r3,r5

        x = xyz(XDIR)
        y = xyz(YDIR)
        z = xyz(ZDIR)

        phi = Om0*t
        Xp = x*cos(phi) + y*sin(phi)
        r = norm2(xyz)
        r3 = r**3.0
        r5 = r**5.0

        Exyz(XDIR) = (3.0/r5)*z*SinThD*x*Xp - (r*r-3*z*z)*CosThD*x/r5
        Exyz(YDIR) = (3.0/r5)*z*SinThD*y*Xp - (r*r-3*z*z)*CosThD*y/r5
        Exyz(ZDIR) = -(Xp/r3)*SinThD*( (3*x*x+3*y*y)/(r*r) - 1 ) &
                   & -(3.0/r5)*CosThD*z*(x*x+y*y)

    end function CorotationE
!---
    !Initialization for Coro Inner BC
    subroutine InitCoroInner(bc,Model,Grid,State,xmlInp)
        class(CoroInnerBC_T), intent(inout) :: bc
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        !Are we on the inner (REMIX) boundary
        if (Model%Ri == 1) then

            !Create holders for coupling electric field
            allocate(bc%inExyz(1:PsiSh  ,Grid%jsg:Grid%jeg  ,Grid%ksg:Grid%keg  ,1:Model%nvInExyz))
            allocate(bc%inEijk(1:PsiSh+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM))
            bc%inExyz = 0.0
            bc%inEijk = 0.0
        endif
    end subroutine InitCoroInner

!---
    !Inner-I BC for inner corotation boundary
    subroutine CoroInner(bc,Model,Grid,State)
        class(CoroInnerBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k,ip,jp,kp,ig,n,np
        real(rp) :: rc,xc,yc,zc,Vr,Rin
        real(rp), dimension(NDIM) :: Exyz,Veb,dB,Bd,rHatG,rHatP,Vxyz,Vmir
        real(rp), dimension(NVAR) :: pW,pCon,gW,gCon

        !Get inner radius and low-latitude
        Rin = norm2(Grid%xyz(Grid%is,Grid%js,Grid%ks,:))
        
        !i-boundaries (IN)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,ip,jp,kp,ig,n,np) &
        !$OMP private(rc,xc,yc,zc,Vr) &
        !$OMP private(Exyz,Veb,Bd,dB,rHatG,rHatP,Vxyz,Vmir) &
        !$OMP private(pW,pCon,gW,gCon)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                !Map to active ip,jp,kp (i=Grid%is => ip=Grid%is)
                call lfmIJK(Model,Grid,Grid%is,j,k,ip,jp,kp)

                !Loop inward over ghosts
                do n=1,Model%Ng
                    ig = Grid%is-n
                    ip = Grid%is+n-1
                    !Map n=[1,ng] -> inExyz
                    !ASSUMING PsiSt=-3 if you're nudging, so n=[1,ng]->[4,...,1]
                    np = Model%nG-n+1 !Mapping into 1,4 of inExyz

                !-------
                !Get geometry for this ghost and matching physical

                    !call cellCenter(Grid,ig,jp,kp,xc,yc,zc)
                    !NOTE: Using j/k instead of jp/kp to deal with double-corner sign flip
                    call cellCenter(Grid,ig,j ,k ,xc,yc,zc)
                    rHatG = normVec([xc,yc,zc])

                    call cellCenter(Grid,ip,jp,kp,xc,yc,zc)
                    rHatP = normVec([xc,yc,zc])

                !-------
                !Get E, dipole/perturbation and ExB velocity/Mirror velocity

                    !Get velocity from i-reflected active cell
                    Vmir = State%Gas(ip,jp,kp,MOMX:MOMZ,BLK)/max(State%Gas(ip,jp,kp,DEN,BLK),dFloor)
                    Exyz(IONO_EX) = bc%inExyz(np,jp,kp,IONO_EX)
                    Exyz(IONO_EY) = bc%inExyz(np,jp,kp,IONO_EY)
                    Exyz(IONO_EZ) = bc%inExyz(np,jp,kp,IONO_EZ)
                    call cellCenter(Grid,ig,j ,k ,xc,yc,zc)
                    call TiltDipole(xc,yc,zc,Model%t,Bd(XDIR),Bd(YDIR),Bd(ZDIR))
                    dB = State%Bxyz(ip,jp,kp,:)
                    !Using ExB everywhere
                    Veb = cross(Exyz,Bd)/dot_product(Bd,Bd)
                    Vxyz = Veb - rHatP*dot_product(rHatP,Veb)

                !-------
                !Set ghost hydro quantities
                    !Let density float
                    call SphereWall(Model,State%Gas(ig,j,k,:,:),State%Gas(ip,jp,kp,:,:),Vxyz)
                !-------
                !Now handle magnetic quantities
                    !Mirror fluxes to minimize gradient
                    State%Bxyz(ig,j,k,:) = State%Bxyz(ip,jp,kp,:)
                    State%magFlux(ig,j,k,IDIR) = State%magFlux(ip,jp,kp,IDIR)
                    State%magFlux(ig,j,k,JDIR) = State%magFlux(ip,jp,kp,JDIR)
                    State%magFlux(ig,j,k,KDIR) = State%magFlux(ip,jp,kp,KDIR)

                enddo !n
            enddo
        enddo

    end subroutine CoroInner

end module usergamic
