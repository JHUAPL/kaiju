module usergamic
    use kdefs
    use gamtypes
    use gambctypes
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use ioH5
    use helioutils

    implicit none

    enum, bind(C)
       ! variables passed via innerbc file
       ! Br, Vr, Rho, Temperature, Br @ kface, Vr @ kface, et @ kedge
       enumerator :: BRIN=1,VRIN,RHOIN,TIN,BRKFIN,VRKFIN,ETKEDGE
    endenum 


    integer, private, parameter :: NVARSIN=7 ! SHOULD be the same as the number of vars in the above enumerator
    real(rp), dimension(:,:,:,:), allocatable :: ibcVars
    real(rp), dimension(:,:), allocatable :: ibcEt

    !Various global would go here
    real (rp) :: Rho0, P0, Vslow,Vfast, wScl, Cs0, B0, MJD_c
 
    !Scaling from innerbc to CGS
    !Elena: May be move these conversion factors to kdefs?
    real (rp) :: km2cm = 1.e5
    real (rp) :: nT2Gs = 1.e-5

    ! things we keep reusing
    real(rp), dimension(NDIM) :: xyz,xyz0,rHat,phiHat
    real(rp) :: Rfactor

    ! global grid
    integer :: gNkp

    ! FIXME
    ! setting it here temporarily. eventually need to read from HDF or something
    ! note also that this is slightly incorrect, since Rbc below is used as the radius
    ! of the center of the first ghost cell.

    !for inner heliosphere
    real(rp) :: Rbc = 21.5

    real(rp) :: Tsolar ! Solar rotation period, defined in apps/helioutils.F90

    !DeltaT sets a rotation of the wsa map relative to +X at simulation time t=0 
    !DeltaT = 0 then the map joint is at +X at t=0
    !DeltaT = Tsolar/2 then the map joint is at -X at t=0 (facing Earth)
    real(rp) :: DeltaT 

    character(len=strLen) :: wsaFile
 
    ! use this to fix the Efield at the inner boundary
    real(rp), allocatable :: inEijk(:,:,:,:)

    ! type for solar wind BC
    type, extends(innerIBC_T) :: SWInnerBC_T

        !Main electric field structures
        real(rp), allocatable, dimension(:,:,:,:) :: inEijk,inExyz

        contains

!        procedure :: doInit => InitIonInner
          ! TODO: this shoudl be made generic (wsa, mas, etc.) How
        procedure :: doBC => wsaBC
    end type SWInnerBC_T

    contains

    !> Initialize the run
    !> Set initial conditions
    subroutine initUser(Model,Grid,State,inpXML)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        procedure(GasIC_T), pointer :: Wxyz
        procedure(HackStep_T), pointer :: tsHack
        procedure(HackE_T), pointer :: eHack

        integer :: i,j,k,nvar,nr,d
        integer :: n1, n2
        integer :: kg, ke, kb, jg
        real(rp) :: a
        real(rp) :: Br_left
        real(rp) :: ibcVarsStaticBRKFN !For brkfn
        real(rp) :: ibcVarsStatic !For br only

        if (.not.allocated(inEijk)) allocate(inEijk(1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM))
        inEijk = 0.0

        ! set units and other things, like Tsolar
        call setHeliosphere(Model,inpXML,Tsolar)

        ! grab inner 
        call inpXML%Set_Val(wsaFile,"prob/wsaFile","innerbc.h5" )

        ! compute global Nkp
        gNkp = Grid%Nkp*Grid%NumRk

        ! initial conditions
        ! TODO: change using norm. units in Model, set in helioutils
        Cs0   = 0.267  ! 40 km/s
        Vslow = 1.33   ! 200 km/s
        Vfast = 5.33   ! 800 km/s
        !for inner helio
        B0    = 2.0    ! 200 nT
        Rho0  = 1.0    ! 200/cc
        P0    = 1.0e-4*Rho0*Cs0**2.0/Model%gamma

        !TO DO Elena: Move to wsa.xml
        DeltaT = Tsolar/2.


        ![OHelio] for 1-10 au helio
        !Cs0 = 0.78    !27 km/s
        !B0 = 1.       ! 5nT
        !Rho0 = 1.     ! 10/cc
        !Vslow = 8.5   !300 km/s 

        ! deallocate default BCs
        ! required because defaults are triply-periodic
        ! need to wipe
        call WipeBCs(Model,Grid)

        !Set BCs
        allocate(SWInnerBC_T        :: Grid%externalBCs(1)%p)
        allocate(helioOuterIBC_T    :: Grid%externalBCs(2)%p)
        allocate(helioInnerJBC_T    :: Grid%externalBCs(3)%p)
        allocate(helioOuterJBC_T    :: Grid%externalBCs(4)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(5)%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(6)%p)

        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

        !Fixing Efield at the inner boundary
        eHack  => EFix
        Model%HackE => eHack
   
        ! Write MJD_center_of_WSA_map to the root of H5 output 
        Model%HackIO_0 => writeMJDcH5Root

        ! everybody reads WSA data
        call readIBC(wsaFile,Model)

        !Set MJD0 which corresponds to sim time 0
        !MJD0 is MJD_center_of_WSA_map - Tsolar_synodic/2; Tsolar_synodic = 27.28
        Model%MJD0 = MJD_c - Tsolar_synodic/2.

        !if not restart set State%time according the tSpin
        State%time  = Model%t

        !Map IC to grid
        Wxyz => GasIC
        call GasIC2State(Model,Grid,State,Wxyz)

        ! NOTE, not filling ghost cells here
        ! relying on J, K boundary conditions
        ! first zero everything out
        State%magFlux = 0.0

        do k=Grid%ks,Grid%ke
           kg = k+Grid%ijkShift(KDIR)
           !initial spin of a map by Tsolar/2
           call mapKinit(kg,ke,kb,a)

           do j=Grid%js,Grid%je
              jg = j+Grid%ijkShift(JDIR)

              do i=Grid%isg,Grid%ieg+1
                 ibcVarsStatic = 1._rp

                 if ( (jg>=Grid%js).and.(jg<=size(ibcVars,2)) ) then
                 ! interpolate BRIN linearly from rotating to inertial frame
                    ibcVarsStatic = a*ibcVars(Model%Ng,jg,kb,BRIN)+(1-a)*ibcVars(Model%Ng,jg,ke,BRIN)
                 endif

                 ! FIXME: calc Rbc appropriately above!
                 Rfactor = Rbc/norm2(Grid%xfc(Grid%is,j,k,:,IDIR))
                 ! note scaling br to the first active face
                 ! this is opposite to what we did on the python side
                 ! not elegant, but otherwise we'd have to store both br and br_iface in teh innerbc.h5 file
                 !
                 ! note also that ibcVars(Model%Ng,:,:,:) corresponds to the center of the first ghost cell (just below the boundary)
                 !State%magFlux(i,j,k,IDIR) = ibcVars(Model%Ng,j+Grid%ijkShift(JDIR),k+Grid%ijkShift(KDIR),BRIN)*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
                 State%magFlux(i,j,k,IDIR) = ibcVarsStatic*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
              enddo
           enddo
        enddo

        contains
            !> Set initial conditions
            subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
                real (rp), intent(in) :: x,y,z
                real (rp), intent(out) :: D,Vx,Vy,Vz,P
                real (rp) :: r, theta, phi
                real (rp) :: r_unit(NDIM)

                D = Rho0
                P = P0

                !Calculate radial hat vector
                r = sqrt(x**2.0 + y**2.0 + z**2.0)
                theta = acos(z/R)
                phi = atan2(y,x)
                r_unit = [x,y,z]/r

                !Set primitives (already have D/P)
                Vx = r_unit(XDIR)*Vslow
                Vy = r_unit(YDIR)*Vslow
                Vz = r_unit(ZDIR)*Vslow
            end subroutine GasIC

            !> Rotates initial WSA map so that at simulation time 0
            !> the start date of the map faces Earth (-X axis)
            subroutine mapKinit(k,ke,kb,a)
                ! find the lower and upper k-index of the rotating cell on the static grid
                integer, intent(in) :: k
                integer, intent(out) :: ke,kb  ! upper (ke) and lower (kb) indices
                real(rp), intent(out):: a      ! interp coefficient

                real(rp) :: kprime

                kprime = modulo(k-gNkp*(State%time+DeltaT)/Tsolar,real(gNkp,kind=rp))

                if ((kprime.ge.0).and.(kprime.lt.1)) then
                    kb=gNkp
                    ke=1
                else
                    kb=floor(kprime)
                    ke=ceiling(kprime)
                endif
                a = ke-kprime

            end subroutine mapKinit

    end subroutine initUser

    !> Setting inner-I boundary conditions in four ghost i-layers of cells
    !> Linear interpolation from rotating to inertial grid
    !> Setting cell-centered variables and magnetic fluxes at cell faces
    !> Calculate user defined Efield at the inner boundary
    subroutine wsaBC(bc,Model,Grid,State)
      class(SWInnerBC_T), intent(inout) :: bc
      class(Model_T), intent(in) :: Model
      class(Grid_T), intent(in) :: Grid
      class(State_T), intent(inout) :: State
      procedure(HackE_T), pointer :: eHack

      ! local variables
      integer :: i,j,k, kb, ke
      integer :: kg, jg, ig ! global indices
      integer :: var ! ibcVar variable number
      real(rp) :: a
      real(rp) :: ibcVarsStatic(NVARSIN)
      real(rp) :: R, Theta, Phi
      real(rp) :: Theta_kf, R_kf ! kface
      real(rp), dimension(NVAR) :: conVar, pVar
      real(rp) :: Br_left

      !i-boundaries (IN)
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,k,jg,kg,ke,kb,a,var,xyz,R,Theta,Phi,rHat,phiHat) &
      !$OMP private(ibcVarsStatic,pVar,conVar,xyz0,R_kf,Theta_kf)
      do k=Grid%ksg,Grid%keg+1  ! note, going all the way to last face for mag fluxes
         kg = k+Grid%ijkShift(KDIR)
         ! map rotating to static grid
         call mapK(kg,ke,kb,a)

         do j=Grid%js,Grid%je+1 ! note, going all the way to last face for mag fluxes
            jg = j+Grid%ijkShift(JDIR)

            do ig=1,Model%Ng
               i=ig-Model%Ng
               ! note, ibcVars are not defined in the global j-corners or k-corners
               ! access to k-corners is fine because mapK function will always map into active domain (FIXME: check for k=nk+1 face!)
               ! access to j-corners gets into unallocated space for ibcVars. Use this trick instead:
               ! set everything arbitrarily to 1. in the global corners (on low and high j boundaries)
               ! we then apply the j-boundary after i boundary anyway, so the corners will be overwritten
               ibcVarsStatic = 1._rp  

               ! otherwise
               ! interpolate linearly from rotating to inertial frame 
               if ( (jg>=Grid%js).and.(jg<=size(ibcVars,2)) ) then
                  do var=1,NVARSIN
                     ibcVarsStatic(var) = a*ibcVars(ig,jg,kb,var)+(1-a)*ibcVars(ig,jg,ke,var)
                  end do
                  ! use Br value from the left (kb) cell of the rotating grid  to calculate Ej
                  ! Ensures constant Ej between Br updates and linearly changing Br in time. 
                  Br_left = ibcVars(ig,jg,kb,BRIN)
               end if

               ! do cell centered things for cell-centers only
               if ( (j/=Grid%jeg+1).and.(k/=Grid%keg+1) ) then
                  ! various geometrical quantities for the cell center
                  xyz  = Grid%xyzcc(i,j,k,:)
                  R    = norm2(xyz)
                  Theta = acos(xyz(3)/R)
                  Phi = atan2(xyz(2),xyz(1))
                  rHat = xyz/R
                  phiHat = [-sin(phi),cos(phi),0._rp]

                  ! NOTE, WSA data were already scaled appropriately in the python code

                  !Set primitives
                  pVar(VELX:VELZ) = rHat*ibcVarsStatic(VRIN)
                  ! note conversion to my units with B0^2/4pi in the denominator
                  pVar(PRESSURE)  = ibcVarsStatic(RHOIN)*Model%Units%gD0*Kbltz*ibcVarsStatic(TIN)/(Model%Units%gP0)
                  pVar(DEN)       = ibcVarsStatic(RHOIN)

                  !Swap prim->con in ghost variables
                  call CellP2C(Model,pVar,conVar)
                  State%Gas(i,j,k,:,BLK) = conVar

                  ! note, no need to set cc Bxyz because we run flux2field through ghosts
               end if

               ! also need to calculate theta at k-face for k-flux
               ! although we're assuming theta and R don't change from cell center to k-face center,
               ! xyzcc used under if statement above is only defined for cell centers so need to define it here
               xyz0 = Grid%xfc(i,j,k,:,KDIR) !just reusing a temp var here.
               R_kf = norm2(xyz0)
               Theta_kf = acos(xyz0(ZDIR)/R_kf)

               ! Set magnetic fluxes
               ! note scaling for Bi. See note above in InitUser
               Rfactor = Rbc/norm2(Grid%xfc(Grid%is,j,k,:,IDIR))
               State%magFlux(i,j,k,IDIR) = ibcVarsStatic(BRIN)*Rfactor**2*Grid%face(Grid%is,j,k,IDIR)
               State%magFlux(i,j,k,JDIR) = 0.0
               State%magFlux(i,j,k,KDIR) = -2*PI/Tsolar*R_kf*sin(Theta_kf)/ibcVarsStatic(VRKFIN)*ibcVarsStatic(BRKFIN)*Grid%face(i,j,k,KDIR)

               ! Calculate user`s electric fields at the inner boundary
               if (ig == Model%Ng) then
                    ! calculate E_theta at edges using Br_left value defined above
                    inEijk(Grid%is,j,k,JDIR) = -2*PI/Tsolar*Rbc*sin(Theta_kf)*Br_left*Rfactor**2
               end if

            end do
               
         end do
      end do
     
      contains 

      !> Find the lower and upper k-index for interpolation of  the rotating cell on the static grid
      subroutine mapK(k,ke,kb,a)
        integer, intent(in) :: k
        integer, intent(out) :: ke,kb  ! upper (ke) and lower (kb) indices
        real(rp), intent(out):: a      ! interp coefficient

        real(rp) :: kprime

        kprime = modulo(k-gNkp*(State%time+DeltaT)/Tsolar,real(gNkp,kind=rp)) 

        if ((kprime.ge.0).and.(kprime.lt.1)) then
           kb=gNkp 
           ke=1
        else
           kb=floor(kprime)
           ke=ceiling(kprime)
        endif
        a = ke-kprime

      end subroutine mapK

    end subroutine wsaBC

    !> Set user defined tangential Efield to enforce rotating WSA map
    subroutine EFix(Model,Gr,State)
      class(Model_T), intent(in) :: Model
      class(Grid_T), intent(inout) :: Gr
      class(State_T), intent(inout) :: State

      if (Gr%hasLowerBC(IDIR)) then
          State%Efld(Gr%is,:,:,JDIR) = inEijk(1,:,:,JDIR)*Gr%edge(Gr%is  ,:,:,JDIR)
          State%Efld(Gr%is,:,:,KDIR) = 0.
      endif

    end subroutine EFix

    !> Read in innerbc.h5 file with inner boundary conditions
    subroutine readIBC(ibcH5,Model)
      character(len=*), intent(in) :: ibcH5
      type(Model_T), intent(in) :: Model
      logical :: fExist
      integer :: i,nvar,dims(3), dimsE(2)
      integer, parameter :: MAXIOVAR = 50
      type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
    
      character(len=10) :: Grp = "Step#0"

      !Reset IO chain
      call ClearIO(IOVars)

      inquire(file=ibcH5,exist=fExist)
      if (.not. fExist) then
         !Error out and leave
         write(*,*) 'Unable to open innerbc file, exiting'
         stop
      endif

      !Setup input chain
      call AddInVar(IOVars,"vr",       vSclO=km2cm/Model%Units%gv0)
      call AddInVar(IOVars,"vr_kface", vSclO=km2cm/Model%Units%gv0)
      call AddInVar(IOVars,"rho",      vSclO=1/Model%Units%gD0)
      call AddInVar(IOVars,"temp")     !Temp in K
      call AddInVar(IOVars,"br",       vSclO=nT2Gs/Model%Units%gB0)
      call AddInVar(IOVars,"br_kface", vSclO=nT2Gs/Model%Units%gB0)
      call AddInVar(IOVars,"et_kedge", vSclO=1.e-3/(Model%Units%gv0*1.e-2*Model%Units%gB0*1.e-4))
      call AddInVar(IOVars,"MJD")

      call ReadVars(IOVars,.false.,ibcH5,Grp) !Don't use io precision

      dims=IOVars(1)%dims(1:3) ! i,j,k
      dimsE(:) = dims(2:3) ! (Nj,Nk) for E_theta 

      if (.not.allocated(ibcVars)) allocate(ibcVars(dims(1),dims(2),dims(3),NVARSIN))
      if (.not.allocated(ibcEt))allocate(ibcEt(dimsE(1),dimsE(2)))

      ! Zero out
      ibcVars = 0.
      ibcEt = 0.


      do i=1,NVARSIN
         select case (i)
         case (BRIN)
            nvar= FindIO(IOVars,"br")
         case (VRIN)
            nvar= FindIO(IOVars,"vr")
         case (RHOIN)
            nvar= FindIO(IOVars,"rho")
         case (TIN)
            nvar= FindIO(IOVars,"temp")
         case (BRKFIN)
            nvar= FindIO(IOVars,"br_kface")
         case (VRKFIN)
            nvar= FindIO(IOVars,"vr_kface")
         case (ETKEDGE)
            nvar= FindIO(IOVars,"et_kedge")
         end select

            if (i <= NVARSIN-1) then
               ibcVars(:,:,:,i) = reshape(IOVars(nvar)%data*IOVars(nvar)%scale,dims)
            else if (i == ETKEDGE) then
               ibcEt = reshape(IOVars(nvar)%data*IOVars(nvar)%scale,dimsE)
            end if
      end do

      !reading modified julian date from innerbc
      MJD_c = GetIOReal(IOVars,"MJD")
   
    end subroutine readIBC

    !> Writing Modified Julian Date at the center of WSA map to the root of h5 file
    !> MJDc is required to set the GAMERA Helio coordinate system
    subroutine writeMJDcH5Root(Model,Grid,IOVars)
        class(Model_T), intent(in)    :: Model
        class(Grid_T) , intent(in)    :: Grid
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars

        call AddOutVar(IOVars,"MJDc", MJD_c)

    end subroutine writeMJDcH5Root

end module usergamic
