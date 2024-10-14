!Toy gamera/remix magnetosphere

module uservoltic
    use gamtypes
    use gamdebug
    use cmidefs
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use background
    use msphutils
    use msphingest
    use wind
    use multifluid

    implicit none
!Earth scaling
!x0 = 6.38e6 m [1 RE]
!rho0 = 1.67e-21 kg/m^3 [1 particle/cc]
!v0 = 10e5 m/s [100 km/s]
!B0 = 4.58 nT
!P0 = 1.67e-11
!t0 = 63.8 second
!M0 = -0.31*1.0e+5/B0

    !Various module variables
    real(rp), private :: Rho0,P0

    ! type for remix BC
    type, extends(innerIBC_T) :: IonInnerBC_T
        !Main electric field structures
        real(rp), allocatable, dimension(:,:,:,:) :: inEijk,inExyz
        real(rp) :: dtCpl !Coupling timescale [code time]
        logical :: doIonPush
        integer :: nIonP

        contains

        procedure :: doInit => InitIonInner
        procedure :: doBC => IonInner

    end type IonInnerBC_T

    contains

    subroutine initUser(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML
        procedure(GasIC_T), pointer :: Wxyz
        procedure(VectorField_T), pointer :: Axyz
        procedure(HackE_T), pointer :: eHack
        procedure(HackStep_T), pointer :: tsHack

        real(rp) :: M0g
        integer :: s

        !Set user hack functions
        !NOTE: Need silly double value for GNU
        eHack  => NULL()
        tsHack => NULL()
        Model%HackE => eHack
        Model%HackStep => tsHack

        !Get defaults from input deck

        !Density for magnetosphere/wind
        call inpXML%Set_Val(Rho0 ,"prob/Rho0",0.2_rp)
        call inpXML%Set_Val(P0   ,"prob/P0"  ,0.001_rp)
        
        !Set magnetosphere parameters
        call setMagsphere(Model,inpXML)
        P0 = P0/Model%Units%gP0 !Scale to magsphere units

        ! deallocate default BCs
        call WipeBCs(Model,Grid)

        !Set BCs (spherical, RPT)
        allocate(IonInnerBC_T       :: Grid%externalBCs(1)%p)
        allocate(WindBC_T           :: Grid%externalBCs(2)%p)
        allocate(lfmInBC_T          :: Grid%externalBCs(3)%p)
        allocate(lfmOutBC_T         :: Grid%externalBCs(4)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(5)%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(6)%p)

        !Setup fields
        !Use cutoff dipole
        call genCutDipole(Model,Grid,State,inpXML)

        !Map IC to grid
        Wxyz => GasIC
        if (Model%doMultiF) then
            !Set fluid 1
            call GasIC2State(Model,Grid,State,Wxyz,1)
            !Set fluid 2 to plasmasphere
            Wxyz => PSphereIC
            call GasIC2State(Model,Grid,State,Wxyz,2)
            !Null remaining fluids
            do s=3,Model%nSpc
                State%Gas(:,:,:,DEN      ,s) = dFloor
                State%Gas(:,:,:,MOMX:MOMZ,s) = 0.0
                State%Gas(:,:,:,ENERGY   ,s) = pFloor
            enddo
            !Accumulate
            call State2Bulk(Model,Grid,State)
        else
            call GasIC2State(Model,Grid,State,Wxyz)
        endif

        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

        !Trap for unsupported cases
        if (.not. ( Grid%hasLowerBC(KDIR) .and. Grid%hasUpperBC(KDIR) ) ) then
            write(*,*) 'K-decomposition not yet supported for magnetosphere, bailing ...'
            stop
        endif
        
    !Set user hack functions
    !NOTE: Need silly double value for GNU

        !For everybody
        eHack  => EFix
        Model%HackE => eHack
        tsHack => PerStep
        Model%HackStep => tsHack

        !Local functions
        !NOTE: Don't put BCs here as they won't be visible after the initialization call
        contains
            subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P

                real(rp) :: r,lambda,cL
                real(rp) :: phi,L,Lx,Ly
                real(rp) :: M

                r = sqrt(x**2.0+y**2.0+z**2.0)
                lambda = asin(z/r)
                cL = cos(lambda)
                L = r/(cL*cL)

                phi = atan2(y,x)
                Lx = L*cos(phi)
                Ly = L*sin(phi)

                P = P0
                D = Rho0
                
                Vx = 0.0
                Vy = 0.0
                Vz = 0.0
            end subroutine GasIC

            subroutine PSphereIC(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P
                real(rp) :: r,lambda,cL,L

                r = sqrt(x**2.0+y**2.0+z**2.0)
                lambda = asin(z/r)
                cL = cos(lambda)
                L = r/(cL*cL)
                D = max(dFloor,psphD(L))
                P = P0
                Vx = 0.0
                Vy = 0.0
                Vz = 0.0

            end subroutine PSphereIC

    end subroutine initUser

    !Routines to do every timestep
    subroutine PerStep(Model,Gr,State)
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(inout) :: Gr
        class(State_T), intent(inout) :: State

        integer :: nbc

        !Call ingestion function
        if (Model%doSource) then
            call MagsphereIngest(Model,Gr,State)
        endif

        !Do heavier nudging to first shell
        if ( Gr%hasLowerBC(IDIR) ) then
            nbc = FindBC(Model,Gr,INI)
            SELECT type(iiBC=>Gr%externalBCs(nbc)%p)
                TYPE IS (IonInnerBC_T)
                    if (iiBC%doIonPush) call PushIon(iiBC,Model,Gr,State)
                    
                CLASS DEFAULT
                    write(*,*) 'Could not find Ion Inner BC in PerStep'
                    stop
            END SELECT
        endif

        !Do some nudging at the outermost cells to hit solar wind
        if (Gr%hasUpperBC(IDIR)) then
            nbc = FindBC(Model,Gr,OUTI)

            SELECT type(pWind=>Gr%externalBCs(nbc)%p)
                TYPE IS (WindBC_T)
                    call NudgeSW(pWind,Model,Gr,State)
            END SELECT
        endif

    end subroutine PerStep

    !Fixes electric field before application
    subroutine EFix(Model,Gr,State)
        class(Model_T), intent(in)    :: Model
        class(Grid_T) , intent(inout) :: Gr
        class(State_T), intent(inout) :: State

        integer :: nbc

        !Fix inner shells
        if (Gr%hasLowerBC(IDIR)) then
            nbc = FindBC(Model,Gr,INI)
            SELECT type(iiBC=>Gr%externalBCs(nbc)%p)
                TYPE IS (IonInnerBC_T)
                    call IonEFix(Model,Gr,State,iiBC%inEijk)
                CLASS DEFAULT
                    write(*,*) 'Could not find Ion Inner BC in EFix'
                    stop
            END SELECT
        endif

        !Fix outer shells
        if (Gr%hasUpperBC(IDIR)) then
            nbc = FindBC(Model,Gr,OUTI)
            SELECT type(pWind=>Gr%externalBCs(nbc)%p)
                TYPE IS (WindBC_T)
                   call WindEFix(pWind,Model,Gr,State)
            CLASS DEFAULT
                write(*,*) 'Could not find Wind BC in EFix'
                stop
            END SELECT
        endif

        call FixEFieldLFM(Model,Gr,State%Efld)
        
    end subroutine EFix

    !Ensure no flux through degenerate faces
    subroutine IonFlux(Model,Gr,gFlx,mFlx)
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(in) :: Gr
        real(rp), intent(inout) :: gFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(inout), optional :: mFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM,1:NDIM)

        integer :: n,s,j,k
        real(rp) :: igFlx(NVAR), imFlx(NDIM)

        !This is inner-most I tile
        if ( Gr%hasLowerBC(IDIR) .and. (.not. Model%doMultiF) ) then

            if (Model%doRing) then
                if (Model%Ring%doE) then
                    igFlx = sum(gFlx(Gr%is,Gr%je,Gr%ks:Gr%ke,1:NVAR,IDIR,BLK),dim=1)/Model%Ring%Np
                    imFlx = sum(mFlx(Gr%is,Gr%je,Gr%ks:Gr%ke,1:NDIM,IDIR    ),dim=1)/Model%Ring%Np
                    do k=Gr%ks,Gr%ke
                        gFlx(Gr%is,Gr%je,k,:,IDIR,BLK) = igFlx
                        mFlx(Gr%is,Gr%je,k,:,IDIR    ) = imFlx
                    enddo
                endif !doE

                if (Model%Ring%doS) then
                    igFlx = sum(gFlx(Gr%is,Gr%js,Gr%ks:Gr%ke,1:NVAR,IDIR,BLK),dim=1)/Model%Ring%Np
                    imFlx = sum(mFlx(Gr%is,Gr%js,Gr%ks:Gr%ke,1:NDIM,IDIR    ),dim=1)/Model%Ring%Np
                    do k=Gr%ks,Gr%ke
                        gFlx(Gr%is,Gr%js,k,:,IDIR,BLK) = igFlx
                        mFlx(Gr%is,Gr%js,k,:,IDIR    ) = imFlx
                    enddo
                endif !doE

            endif !doRing
            
        endif !Inner i-tile and not MF

    end subroutine IonFlux

    !Initialization for Ion Inner BC
    subroutine InitIonInner(bc,Model,Grid,State,xmlInp)
        class(IonInnerBC_T), intent(inout) :: bc
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(in) :: Grid
        class(State_T), intent(in) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        integer :: PsiShells
        real(rp) :: dtXML
        procedure(HackE_T), pointer :: eHack

        !Are we on the inner (REMIX) boundary
        if (Grid%hasLowerBC(IDIR)) then
            PsiShells = PsiSh !Coming from cmidefs

            !Create holders for coupling electric field
            allocate(bc%inExyz(1:PsiShells  ,Grid%jsg:Grid%jeg  ,Grid%ksg:Grid%keg  ,1:NDIM))
            allocate(bc%inEijk(1:PsiShells+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM))
            bc%inExyz = 0.0
            bc%inEijk = 0.0
            eHack  => EFix
            Model%HackE => eHack
            Model%HackFlux => IonFlux
            !Set value of coupling timescale
            call xmlInp%Set_Val(dtXML,"/Kaiju/voltron/coupling/dt",5.0)
            bc%dtCpl = dtXML/Model%Units%gT0
            !Get knobs for pushing
            call xmlInp%Set_Val(bc%doIonPush,"ibc/doIonPush",.true.)
            call xmlInp%Set_Val(bc%nIonP,"ibc/nIonP",1)
        endif

    end subroutine InitIonInner

    !Inner-I BC for ionosphere
    subroutine IonInner(bc,Model,Grid,State)
        class(IonInnerBC_T), intent(inout) :: bc
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(in) :: Grid
        class(State_T), intent(inout) :: State

        real(rp) :: Rin,llBC,dA,Rion,MagB0
        real(rp), dimension(NDIM) :: Bd,Exyz,Veb,rHat
        integer :: i,j,k,ig,jg,kg,ip,jp,kp,idip,n,np,d
        !integer :: ig,ip,idip,j,k,jp,kp,n,np,d
        integer, dimension(NDIM) :: dApm

        
        !Are we on the inner (REMIX) boundary
        if (.not. Grid%hasLowerBC(IDIR)) return

        Rion = RadIonosphere()
        !Get inner radius and low-latitude
        Rin = norm2(Grid%xyz(Grid%is,Grid%js,Grid%ks,:))
        llBC = 90.0 - rad2deg*asin(sqrt(Rion/Rin)) !co-lat -> lat

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(ig,ip,idip,j,k,jp,kp,n,np,d) &
        !$OMP private(Bd,Exyz,Veb,rHat,dA,dApm,MagB0)
        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1

                !Loop inward over ghosts
                do n=1,Model%Ng
                    ig = Grid%is-n
                    ip = Grid%is+n-1

                    !Map n=[1,ng] -> inExyz
                    !ASSUMING PsiSt=-3 if you're nudging, so n=[1,ng]->[4,...,1]
                    np = Model%nG-n+1 !Mapping into 1,4 of inExyz

                    !Do cell-centered stuff
                    if (isCellCenterG(Model,Grid,ig,j,k)) then
                        !Map to active ip,jp,kp (i=Grid%is => ip=Grid%is)
                        call lfmIJKcc(Model,Grid,Grid%is-1,j,k,idip,jp,kp)

                        !Get dipole value
                        Bd = VecDipole(Grid%xyzcc(idip,jp,kp,:)) !For direction to keep stencil regular
                        MagB0 = norm2( VecDipole(Grid%xyzcc(ig,jp,kp,:)) )

                        !Get remix field
                        Exyz = bc%inExyz(np,jp,kp,:)

                        !ExB velocity
                        !Veb = cross(Exyz,Bd)/dot_product(Bd,Bd)
                        Veb = cross(Exyz,normVec(Bd))/MagB0 !Use is-1 direction, but is-n field strength

                        !Remove radial component of velocity
                        rHat = normVec(Grid%xyzcc(idip,jp,kp,:))
                        Veb = Vec2Perp(Veb,rHat)

                        !Now do spherical wall BC
                        call SphereWall(Model,State%Gas(ig,j,k,:,:),State%Gas(ip,jp,kp,:,:),Veb,Bd)

                    endif !Cell-centered

                !Now do face fluxes
                    dApm(IDIR:KDIR) = 1 !Use this to hold coefficients for singularity geometry

                    if ( (Model%Ring%doS) .and. (j < Grid%js) ) then
                        dApm(JDIR:KDIR) = -1
                    endif
                    if ( (Model%Ring%doE) .and. (j >= Grid%je+1) ) then
                        dApm(JDIR:KDIR) = -1
                    endif
                    
                    !Loop over face directions
                    do d=IDIR,KDIR
                        call lfmIJKfc(Model,Grid,d,ig,j,k,ip,jp,kp)
                        dA = Grid%face(ig,j,k,d)/max(Grid%face(Grid%is,jp,kp,d),TINY)

                        if (isGoodFaceIJK(Model,Grid,ig,j,k,d)) then
                            if ( isLowLat(Grid%xfc(ig,j,k,:,d),llBC) ) then
                                State%magFlux(ig,j,k,d) = 0.0
                            else
                                State%magFlux(ig,j,k,d) = dApm(d)*dA*State%magFlux(Grid%is,jp,kp,d)
                            endif
                        endif
                    enddo !dir

                enddo !n loop (ig)
            enddo !j loop
        enddo !k loop

        !Now that every magflux is set go back and calculate Bxyz
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,ig,j,k,ip,jp,kp)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                !Loop inward over ghosts
                do n=1,Model%Ng
                    ig = Grid%is-n
                    call lfmIJKcc(Model,Grid,ig,j,k,ip,jp,kp)
                    !Use conjugate cell to set Bxyz in this cell
                    State%Bxyz(ig,j,k,:) = CellBxyz(Model,Grid,State%magFlux,ip,jp,kp)
                enddo
            enddo
        enddo

        !Now fix singular regions
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,i,k,ig,jg,kg,ip,jp,kp)
        do k=Grid%ksg,Grid%keg
            do i=Grid%isg,Grid%is-1
                do n=1,Model%Ng
                    if (Model%Ring%doS) then
                        ig = i
                        kg = k
                        jg = Grid%js-n
                        call lfmIJKcc(Model,Grid,ig,jg,kg,ip,jp,kp)
                        State%Bxyz(ig,jg,kg,:) = State%Bxyz(ip,jp,kp,:)
                    endif

                    if (Model%Ring%doE) then
                        ig = i
                        kg = k
                        jg = Grid%je+n
                        call lfmIJKcc(Model,Grid,ig,jg,kg,ip,jp,kp)
                        State%Bxyz(ig,jg,kg,:) = State%Bxyz(ip,jp,kp,:)
                    endif
                enddo !n
            enddo
        enddo !k

        contains 

            function isLowLat(xyz,llBC)
                real(rp), dimension(NDIM), intent(in) :: xyz
                real(rp), intent(in) :: llBC
                logical :: isLowLat

                real(rp) :: invlat

                invlat = rad2deg*InvLatitude(xyz)
                isLowLat = (invlat<=llBC)
            end function isLowLat

    end subroutine IonInner

    !Push velocity of first active cell
    subroutine PushIon(bc,Model,Grid,State)
        class(IonInnerBC_T), intent(inout) :: bc
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(in) :: Grid
        class(State_T), intent(inout) :: State

        integer :: i,j,k,PsiShells,dN
        real(rp) :: dt
        real(rp), dimension(NVAR) :: pW,pCon
        real(rp), dimension(NDIM) :: vMHD,xcc,Bd,Exyz,rHat,Veb,dV,B

        if ( (.not. bc%doIonPush) .or. (.not. Grid%hasLowerBC(IDIR)) ) return

        if (Model%doMultiF) then
            write(*,*) 'PushIon not implemented for MF yet ...'
            !stop
        endif

        PsiShells = PsiSh !Coming from cmidefs

        !Loop over active
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,pW,pCon) &
        !$OMP private(vMHD,xcc,Bd,Exyz,rHat,Veb,dV,dt,dN,B)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=Grid%is,Grid%is+bc%nIonP-1
                    
                !Get MHD info
                    pCon = State%Gas(i,j,k,:,BLK)
                    call CellC2P(Model,pCon,pW)
                    vMHD = pW(VELX:VELZ)
                !Get EB info
                    xcc = Grid%xyzcc(i,j,k,:)
                    Bd = VecDipole(xcc)
                    Exyz = bc%inExyz(PsiShells,j,k,:)
                    rHat = normVec(xcc)
                    if (Model%doBackground) then
                        B = State%Bxyz(i,j,k,:) + Grid%B0(i,j,k,:)
                    else
                        B = State%Bxyz(i,j,k,:)
                    endif
        
                    !Get ExB velocity w/o radial component
                    Veb = cross(Exyz,B)/dot_product(B,B)
                    Veb = Vec2Perp(Veb,rHat)
                    
                !Setup push and finish up
                    dN = i - Grid%is + 1
                    dt = Model%dt/(dN*bc%dtCpl)
                    dt = min(dt,1.0)
                    dV = Veb - vMHD
                    
                    pW(VELX:VELZ) = pW(VELX:VELZ) + dt*dV
                    call CellP2C(Model,pW,pCon)
                    State%Gas(i,j,k,:,BLK) = pCon

                enddo
            enddo
        enddo

    end subroutine PushIon

end module uservoltic
