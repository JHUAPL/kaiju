!Simple test of BW on LFM grid

module usergamic
    use gamtypes
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use background

    implicit none

    !Various global would go here
    real(rp) :: B0,D0,P0,x0,y0,z0,rC,pRat
    real(rp) :: bMag
    logical :: doDip = .false.

    type, extends(innerIBC_T) :: bwiIBC_T
        contains
        procedure :: doBC => bw_ibcI
    end type bwiIBC_T

    type, extends(outerIBC_T) :: bwoIBC_T
        contains
        procedure :: doBC => bw_obcI
    end type bwoIBC_T

    contains

    subroutine initUser(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        procedure(GasIC_T), pointer :: Wxyz
        procedure(VectorField_T), pointer :: Axyz


        !Get defaults from input deck
        call inpXML%Set_Val(x0 ,"prob/x0" ,0.0_rp)
        call inpXML%Set_Val(y0 ,"prob/y0" ,0.0_rp)
        call inpXML%Set_Val(z0 ,"prob/z0" ,0.0_rp)
        call inpXML%Set_Val(P0  ,"prob/P0"  ,0.1_rp  )
        call inpXML%Set_Val(pRat,"prob/pRat",10.0_rp)        
        call inpXML%Set_Val(D0  ,"prob/D0"  ,1.0_rp  )
        call inpXML%Set_Val(rC  ,"prob/rC"  ,0.1_rp  )
        call inpXML%Set_Val(B0  ,"prob/B0",1.0_rp)
        call inpXML%Set_Val(doDip,"prob/doDip",.false.)

        bMag = B0/sqrt(2.0)

        ! deallocate default BCs
        call WipeBCs(Model,Grid)

        !Set BCs (spherical, RPT)
        allocate(bwiIBC_T           :: Grid%externalBCs(INI )%p)
        allocate(bwoIBC_T           :: Grid%externalBCs(OUTI)%p)
        allocate(lfmInBC_T          :: Grid%externalBCs(INJ )%p)
        allocate(lfmOutBC_T         :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(OUTK)%p)


        !Map IC to grid
        Wxyz => GasIC
        call GasIC2State(Model,Grid,State,Wxyz)

        !Map vector potential to initial field
        if (doDip) then
            Axyz => VP_Dipole
        else
            Axyz => VectorPot_BW
        endif

        call VectorPot2Flux(Model,Grid,State,Axyz)

        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke
        
        Model%HackPredictor => PredFix

        !Local functions
        !NOTE: Don't put BCs here as they won't be visible after the initialization call

        contains
            subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P

                real(rp) :: R

                R = sqrt( (x-x0)**2.0 + (y-y0)**2.0 + (z-z0)**2.0)
                D = D0
                P = P0
                Vx = 0.0
                Vy = 0.0
                Vz = 0.0
                if (R <= rC) then
                    P = pRat*P
                endif

            end subroutine GasIC

            subroutine VectorPot_BW(x,y,z,Ax,Ay,Az)
                
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: Ax,Ay,Az
        
                Ax = 0.0
                Ay = 0.0
                Az = bMag*(y - x)
            end subroutine VectorPot_BW

            subroutine VP_Dipole(x,y,z,Ax,Ay,Az)
                
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: Ax,Ay,Az

                real(rp), dimension(NDIM) :: A,m,r,rhat
                m = [0,0,1]
                r = [x,y,z]
                rhat = r/norm2(r)

                A = B0*cross(m,rhat)/(dot_product(r,r))
                Ax = A(XDIR)
                Ay = A(YDIR)
                Az = A(ZDIR)
            end subroutine VP_Dipole

    end subroutine initUser


    !Inner-I BC 
    subroutine bw_ibcI(bc,Model,Grid,State)
        class(bwiIBC_T), intent(inout) :: bc
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(in) :: Grid
        class(State_T), intent(inout) :: State

        integer :: n,ig,ip,j,k,jp,kp,d
        integer, dimension(NDIM) :: dApm
        real(rp) :: dA

        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                
                !Loop inward over ghosts
                do n=1,Model%Ng
                    ig = Grid%is-n
                    ip = Grid%is+n-1
                    !Do cell-centered stuff
                    if (isCellCenterG(Model,Grid,ig,j,k)) then
                        !Fix gas vars
                        State%Gas(ig,j,k,DEN,:)  = D0
                        State%Gas(ig,j,k,MOMX:MOMZ,:)  = 0.0
                        State%Gas(ig,j,k,ENERGY,:) = P0/(Model%gamma-1.0) !Just internal energy
                    endif

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

                        dA = Grid%face(ig,j,k,d)/Grid%face(Grid%is,jp,kp,d)
                        State%magFlux(ig,j,k,d) = dApm(d)*dA*State%magFlux(Grid%is,jp,kp,d)

                    enddo

                enddo !n loop (ig)
            enddo !j loop
        enddo !k loop

    end subroutine bw_ibcI
    
    !Outer-I BC
    subroutine bw_obcI(bc,Model,Grid,State)
        class(bwoIBC_T), intent(inout) :: bc
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(in) :: Grid
        class(State_T), intent(inout) :: State

        integer :: n,ig,di,ip,j,k,jp,kp,d
        integer, dimension(NDIM) :: dApm
        real(rp) :: dA

        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                
                !Loop outward over ghosts
                do n=1,Model%Ng
                    ig = Grid%ie+n
                    ip = Grid%ie-n+1

                    !Do cell-centered stuff
                    if (isCellCenterG(Model,Grid,ig,j,k)) then
                        !Fix gas vars
                        State%Gas(ig,j,k,DEN,:)  = D0
                        State%Gas(ig,j,k,MOMX:MOMZ,:)  = 0.0
                        State%Gas(ig,j,k,ENERGY,:) = P0/(Model%gamma-1.0) !Just internal energy
                    endif

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
                        if (d == IDIR) then
                            di = +1
                        else
                            di = 0
                        endif

                        call lfmIJKfc(Model,Grid,d,ig+di,j,k,ip,jp,kp)

                        dA = Grid%face(ig+di,j,k,d)/Grid%face(Grid%ie+di,jp,kp,d)
                        State%magFlux(ig+di,j,k,d) = dApm(d)*dA*State%magFlux(Grid%ie+di,jp,kp,d)

                    enddo

                enddo !n loop (ig)
            enddo !j loop
        enddo !k loop

    end subroutine bw_obcI

    !Fixes cell-centered fields in the predictor
    subroutine PredFix(Model,Gr,State)
        class(Model_T), intent(in) :: Model
        class(Grid_T), intent(inout) :: Gr
        class(State_T), intent(inout) :: State

        integer :: n,ip,ig,ix,jp,kp,j,k

        if (Gr%hasLowerBC(IDIR)) then
            do k=Gr%ksg,Gr%keg
                do j=Gr%jsg,Gr%jeg
                    do n=1,Model%Ng
                        ip = Gr%is
                        ig = Gr%is-n

                        call lfmIJKcc(Model,Gr,ig,j,k,ix,jp,kp)
                        State%Bxyz(ig,j,k,:) = State%Bxyz(ip,jp,kp,:)

                    enddo !n loop
                enddo !j loop
            enddo !k loop
        endif

    end subroutine PredFix

    ! !Fixes electric field before application
    ! subroutine EFix(Model,Gr,State)
    !     type(Model_T), intent(in) :: Model
    !     type(Grid_T), intent(inout) :: Gr
    !     type(State_T), intent(inout) :: State


    !     if (Gr%hasLowerBC(IDIR)) then
    !         !Zero out E fields on inner shell (too geometric anyways)
    !         State%Efld(Gr%is:Gr%is,:,:,:) = 0.0
    !     endif
    ! end subroutine EFix
end module usergamic
