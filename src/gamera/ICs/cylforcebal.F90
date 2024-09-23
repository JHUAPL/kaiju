!Cylindrical force balance test, see Skinner+ 2010 Sec 10.1
!Typical domain
!R = [1,2], phi = [0,pi/4], T=10

module usergamic
    use gamtypes
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use background

    implicit none

    !Various globals go here
    real(rp) :: P0,Rho0,B0,pMin,pMax,Om0

    type, extends(innerIBC_T) :: cylfbIBC_T
        contains
        procedure :: doBC => cylfb_ibcI
    end type cylfbIBC_T

    type, extends(outerIBC_T) :: cylfbOBC_T
        contains
        procedure :: doBC => cylfb_obcI
    end type cylfbOBC_T

    contains

    subroutine initUser(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        logical :: doRotation

        procedure(GasIC_T), pointer :: Wxyz
        procedure(VectorField_T), pointer :: Bxyz

        doRotation = .false.

        !Get defaults from input deck
        call inpXML%Set_Val(P0  ,"prob/P0"  ,1.0_rp)
        call inpXML%Set_Val(B0  ,"prob/B0"  ,1.0_rp)
        call inpXML%Set_Val(Om0 ,"prob/Om0" ,PI*0.25_rp)
        call inpXML%Set_Val(Rho0,"prob/Rho0",1.0_rp)
        call inpXML%Set_Val(doRotation,"prob/doRotation",doRotation)

        pMin = 0.0
        pMax = 2*PI
        if (.not. doRotation) then
            Om0 = 0
        endif
        
        ! deallocate default BCs
        call WipeBCs(Model,Grid)
        allocate(cylfbIBC_T         :: Grid%externalBCs(INI )%p)
        allocate(cylfbOBC_T         :: Grid%externalBCs(OUTI)%p)
        allocate(periodicInnerJBC_T :: Grid%externalBCs(INJ )%p)
        allocate(periodicOuterJBC_T :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(OUTK)%p)


        !Map IC to grid
        Wxyz => GasIC
        call GasIC2State(Model,Grid,State,Wxyz)

        Bxyz => MagIC
        call VectorField2Flux(Model,Grid,State,Bxyz)

        Model%doGrav = .true.
        if (.not. Model%doGrav) then
        	write(*,*) "This problem requires gravity, bailing ..."
            stop
        endif
        Model%doSphGrav = .false.
        Model%Phi => GravPot

        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

    end subroutine initUser

    subroutine GravPot(x,y,z,pot)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: pot

        real(rp) :: R,Om1,Om2

        R = sqrt( x**2.0 + y**2.0)
        Om1 = -B0*B0/(2*Rho0*R*R)
        Om2 = (Om0*R)**2.0 / 2

        pot = Om1 + Om2
        
    end subroutine GravPot

    subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: D,Vx,Vy,Vz,P

        real(rp) :: R,phi,psi

        R = sqrt( x**2.0 + y**2.0)
        phi = atan2(y,x)
        psi = 2*PI*(phi-pMin)/(pMax-pMin)

        D = Rho0*(1 + sin(psi)**2.0)
        P = P0 + B0*B0*(1 + sin(psi)**2.0)/(2*R*R)

        !Want Vphi = Om0 R
        Vx = -Om0*R*sin(phi)
        Vy =  Om0*R*cos(phi)
        Vz = 0.0

    end subroutine GasIC

    !B = B0 cos(psi)/R R-hat
    subroutine MagIC(x,y,z,Bx,By,Bz)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Bx,By,Bz

        real(rp) :: R,phi,psi

        R = sqrt( x**2.0 + y**2.0)
        phi = atan2(y,x)
        psi = 2*PI*(phi-pMin)/(pMax-pMin)

        Bx = B0*cos(psi)*cos(phi)/R
        By = B0*cos(psi)*sin(phi)/R
        Bz = 0.0

    end subroutine MagIC

    subroutine cylfb_ibcI(bc,Model,Grid,State)
        class(cylfbIBC_T), intent(inout) :: bc
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(in) :: Grid
        class(State_T), intent(inout) :: State

        integer :: df,n,ig,j,k
        real(rp) :: x,y,z,Bx,By,Bz,D,P,Vx,Vy,Vz
        real(rp), dimension(NVAR) :: pW,pCon
        real(rp), dimension(NDIM) :: xcc,xfc

        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                do n=1,Model%Ng
                    ig = Grid%is-n
                    if (isCellCenterG(Model,Grid,ig,j,k)) then
                        !Get cell center
                        xcc = Grid%xyzcc(ig,j,k,:)
                        x = xcc(XDIR)
                        y = xcc(YDIR)
                        z = xcc(ZDIR)

                        call MagIC(x,y,z,Bx,By,Bz)
                        call GasIC(x,y,z,D,Vx,Vy,Vz,P)
                        pW(VELX:VELZ) = [Vx,Vy,Vz]
                        pW(DEN) = D
                        pW(PRESSURE) = P
                        call CellP2C(Model,pW,pCon)

                        State%Gas (ig,j,k,:,BLK)  = pCon
                        State%Bxyz(ig,j,k,:)    = [Bx,By,Bz]
                    endif
                    !Always do flux
                    do df=1,3
                        !Do each face
                        xfc = Grid%xfc(ig,j,k,:,df)
                        x = xfc(XDIR)
                        y = xfc(YDIR)
                        z = xfc(ZDIR)
                        call MagIC(x,y,z,Bx,By,Bz)
                        State%magFlux(ig,j,k,df) = Project2Face(Model,Grid,[Bx,By,Bz],df,ig,j,k)

                    enddo
                enddo
            enddo
        enddo

    end subroutine cylfb_ibcI

    subroutine cylfb_obcI(bc,Model,Grid,State)
        class(cylfbOBC_T), intent(inout) :: bc
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(in) :: Grid
        class(State_T), intent(inout) :: State

        integer :: df,n,ig,j,k
        real(rp) :: x,y,z,Bx,By,Bz,D,P,Vx,Vy,Vz
        real(rp), dimension(NVAR) :: pW,pCon
        real(rp), dimension(NDIM) :: xcc,xfc

        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                do n=1,Model%Ng
                    ig = Grid%ie+n
                    if (isCellCenterG(Model,Grid,ig,j,k)) then
                        !Get cell center
                        xcc = Grid%xyzcc(ig,j,k,:)
                        x = xcc(XDIR)
                        y = xcc(YDIR)
                        z = xcc(ZDIR)

                        call MagIC(x,y,z,Bx,By,Bz)
                        call GasIC(x,y,z,D,Vx,Vy,Vz,P)
                        pW(VELX:VELZ) = [Vx,Vy,Vz]
                        pW(DEN) = D
                        pW(PRESSURE) = P
                        call CellP2C(Model,pW,pCon)

                        State%Gas (ig,j,k,:,BLK)  = pCon
                        State%Bxyz(ig,j,k,:)    = [Bx,By,Bz]
                    endif
                    !Always do flux
                    do df=1,3
                        !Do each face
                        xfc = Grid%xfc(ig,j,k,:,df)
                        x = xfc(XDIR)
                        y = xfc(YDIR)
                        z = xfc(ZDIR)
                        call MagIC(x,y,z,Bx,By,Bz)
                        State%magFlux(ig,j,k,df) = Project2Face(Model,Grid,[Bx,By,Bz],df,ig,j,k)

                    enddo
                enddo
            enddo
        enddo

    end subroutine cylfb_obcI
end module usergamic
