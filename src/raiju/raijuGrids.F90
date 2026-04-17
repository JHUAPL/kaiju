

module raijugrids
    use ioh5
    use xml_input
    use shellgrid

    use raijudefs
    use raijutypes
    use raijuRecon
    use raijuSpeciesHelper
    use planethelper

    implicit none

    interface wrapJcc
        module procedure wrapJcc_real, wrapJcc_int, wrapJcc_logical
    end interface

    interface wrapJcorners
        module procedure wrapJcorners_real, wrapJcorners_int, wrapJcorners_logical
    end interface

    contains

!------
! Spatial grid setup stuff
!------

    subroutine raijuGenUniSphGrid(Model, Grid, iXML)
        type(raijuModel_T), intent(in   ) :: Model
        type(raijuGrid_T) , intent(inout) :: Grid
        type(XML_Input_T), intent(in)    :: iXML

        real(rp), dimension(:), allocatable :: Theta
        real(rp), dimension(:), allocatable :: Phi
        real(rp) :: dTheta, dPhi, thetaL, thetaU
        integer :: Nt,Np,Ng
        integer, dimension(4) :: nGhosts
        integer :: i


        call iXML%Set_Val(thetaL , "grid/ThetaL", 15.0)
            !! Lower colat boundary [deg], ~15 Re in dipole
        call iXML%Set_Val(thetaU , "grid/ThetaU", 45.0)
            !! Upper colat boundary [deg], 2 Re in dipole
        call iXML%Set_Val(Nt, "grid/Nt", 31 )  ! 1 deg resolution
        call iXML%Set_Val(Np, "grid/Np", 361)  ! 1 deg resolution
        call iXML%Set_Val(Ng, "grid/Ng", 4  )  ! Number of ghosts, in every direction for now

        ! Turn degrees into radians
        thetaL = thetaL*deg2rad
        thetaU = thetaU*deg2rad

        ! Express numGhosts in the way GenShellGrid expects
        nGhosts = Ng

        ! Allocate arrays
        allocate(Theta(Nt+1))
        allocate(Phi(Np+1))

        ! Create uniform grids
        dTheta = (thetaU-thetaL)/Nt
        dPhi = 2.0*PI/Np

        do i=1,Nt+1
            Theta(i) = thetaL + (i-1)*dTheta
        enddo

        do i=1,Np+1
            Phi(i) = (i-1)*dPhi
        enddo
        ! Catch for slight overshoot
        if ((Phi(Np+1) > 2*PI) .and. (Phi(Np+1) - 2*PI < TINY) ) then
            Phi(Np+1) = 2.0*PI
        endif
        
        call GenShellGrid(Grid%shGrid,Theta,Phi,RAI_SG_NAME,nGhosts=nGhosts,radO=Model%planet%Ri_m/Model%planet%Rp_m)

    end subroutine raijuGenUniSphGrid



    subroutine raijuGenWarpSphGrid_Fok2021(Model, Grid, iXML)
        !! Currently just implementation of Fok+2021 scaling of L = 1/cos^n(lat)
        type(raijuModel_T), intent(in   ) :: Model
        type(raijuGrid_T) , intent(inout) :: Grid
        type(XML_Input_T), intent(in)    :: iXML

        real(rp), dimension(:), allocatable :: Theta
        real(rp), dimension(:), allocatable :: Phi
        real(rp) :: dPhi, thetaL, thetaU
        real(rp) :: sL, sU, ds, s
        integer :: Nt,Np,Ng
        integer :: n
        integer, dimension(4) :: nGhosts
        integer :: i

        call iXML%Set_Val(thetaL , "grid/ThetaL", 15.0)
            !! Lower colat boundary [deg], ~15 Re in dipole
        call iXML%Set_Val(thetaU , "grid/ThetaU", 45.0)
            !! Upper colat boundary [deg], 2 Re in dipole
        call iXML%Set_Val(n , "grid/FokN", 5 )
        call iXML%Set_Val(Nt, "grid/Nt", 71 )
        call iXML%Set_Val(Np, "grid/Np", 361)  ! 1 deg resolution
        call iXML%Set_Val(Ng, "grid/Ng", 4  )  ! Number of ghosts, in every direction for now

        ! Turn degrees into radians
        thetaL = thetaL*deg2rad
        thetaU = thetaU*deg2rad

        ! Get mapping of high/low theta boundaries
        sL = sin(thetaL)**(-1.0*n)
        sU = sin(thetaU)**(-1.0*n)
        ds = (sU - sL) / Nt
        ! Note: ds only equals dipole dL if n=2, otherwise it doesn't map exactly to equatorial spacing
        ! e.g. Equatorial spacing is only constant (assuming dipole) when n=2
        dPhi = 2.0*PI/Np

        ! Express numGhosts in the way GenShellGrid expects
        nGhosts = Ng

        ! Allocate arrays
        allocate(Theta(Nt+1))
        allocate(Phi(Np+1))

        do i=1,Nt+1
            s = sL + (i-1)*ds 
            Theta(i) = asin(s**(-1.0/n))
        enddo

        do i=1,Np+1
            Phi(i) = (i-1)*dPhi
        enddo
        ! Catch for slight overshoot
        if ((Phi(Np+1) > 2*PI) .and. (Phi(Np+1) - 2*PI < TINY) ) then
            Phi(Np+1) = 2.0*PI
        endif

        call GenShellGrid(Grid%shGrid,Theta,Phi,RAI_SG_NAME,nGhosts=nGhosts,radO=Model%planet%Ri_m/Model%planet%Rp_m)        
        
    end subroutine raijuGenWarpSphGrid_Fok2021


    subroutine raijuGenWarpSphGrid_Shafee2008(Model, Grid, iXML)
        !! Modified version of Shafee+ 2008 (10.1086/593148)
        type(raijuModel_T), intent(in   ) :: Model
        type(raijuGrid_T) , intent(inout) :: Grid
        type(XML_Input_T), intent(in)    :: iXML

        real(rp), dimension(:), allocatable :: Theta
        real(rp), dimension(:), allocatable :: Phi
        real(rp) :: dPhi, thetaL, thetaU, x, dX
        integer :: Nt,Np,Ng
        integer, dimension(4) :: nGhosts
        integer :: i
        integer :: nPow
            !! Power applied to non-linear scaling term
        real(rp) :: hWgt
            !! Weight between linear and non-linear term. 1=linear, 0=non-linear
        real(rp) :: x_low, x_high
            !! Bounds for generating x values between 0 and 1
            !! FIXME: replace with theta_center and x_scale, and calculate our x_low and x_high from that
        real(rp) :: a1, a2
            !! Coefficients to map from (x_low, x_high) to (thetaL, thetaU), using the scaling function
            !! If x_low=0 and x_high=1, a1 = thetaL and a2 = thetaU
            !! But this makes the resolution focus symmetric between thetaLa nd thetaU which is not necessarily what we want

        call iXML%Set_Val(thetaL , "grid/ThetaL", 15.0)
            !! Lower colat boundary [deg], ~15 Re in dipole
        call iXML%Set_Val(thetaU , "grid/ThetaU", 45.0)
            !! Upper colat boundary [deg], 2 Re in dipole
        call iXML%Set_Val(Nt, "grid/Nt", 71 )
        call iXML%Set_Val(Np, "grid/Np", 361)  ! 1 deg resolution
        call iXML%Set_Val(Ng, "grid/Ng", 4  )  ! Number of ghosts, in every direction for now
        call iXML%Set_Val(nPow,"grid/nPow", 7)
        call iXML%Set_Val(hWgt,"grid/h"   , 0.7)
        call iXML%Set_Val(x_low ,"grid/xLow", 0.1)
        call iXML%Set_Val(x_high,"grid/xHigh", 0.95)

        ! Turn degrees into radians
        thetaL = thetaL*deg2rad
        thetaU = thetaU*deg2rad

        if (x_low < 0 ) then
            write(*,*) "ERROR in raijuGenWarpSphGrid_Shafee2008:"
            write(*,*) "xLow must be greater than zero"
            stop
        endif
        if (x_high > 1) then
            write(*,*) "ERROR in raijuGenWarpSphGrid_Shafee2008:"
            write(*,*) "xHigh must be less than one"
            stop
        endif
        if (x_low > x_high) then
            write(*,*) "ERROR in raijuGenWarpSphGrid_Shafee2008:"
            write(*,*) "xLow must be less than xHigh"
            stop
        endif

        call calcCoeffs(thetaL, thetaU, x_low, x_high, hWgt, nPow, a1, a2)
        dX = (x_high - x_low) / Nt

        dPhi = 2.0*PI/Np

        ! Express numGhosts in the way GenShellGrid expects
        nGhosts = Ng

        ! Allocate arrays
        allocate(Theta(Nt+1))
        allocate(Phi(Np+1))

        do i=1,Nt+1
            x = x_low + (i-1)*dX
            Theta(i) = a1 + (a2-a1)*scaleFunc(x, hWgt, nPow)
        enddo

        do i=1,Np+1
            Phi(i) = (i-1)*dPhi
        enddo
        ! Catch for slight overshoot
        if ((Phi(Np+1) > 2*PI) .and. (Phi(Np+1) - 2*PI < TINY) ) then
            Phi(Np+1) = 2.0*PI
        endif

        call GenShellGrid(Grid%shGrid,Theta,Phi,RAI_SG_NAME,nGhosts=nGhosts,radO=Model%planet%Ri_m/Model%planet%Rp_m)

        contains

        ! Job security
        subroutine calcCoeffs(thL, thU, xL, xU, wgt, n, coeff1, coeff2)
            !! Calculates coefficients to use in final theta grid generation
            real(rp), intent(in) :: thL, thU, xL, xU
            real(rp), intent(in) :: wgt
            integer, intent(in) :: n
            real(rp), intent(out) :: coeff1, coeff2

            real(rp) :: f_l, f_u

            ! thetaL = c1 + (c2-c1)*scaleFunc(xL)
            ! thetaU = c1 + (c2-c1)*scaleFunc(xU)
            ! Solve for c1 and c2

            f_l = scaleFunc(xL, wgt, n)
            f_u = scaleFunc(xU, wgt, n)

            coeff1 = (thL*f_u - thU*f_l) / (f_u - f_l)
            coeff2 = (thU - a1)/f_u + a1

        end subroutine calcCoeffs

        function scaleFunc(x, wgt, n) result(val)
            !! for x=[0,1], maps with linear and non-linear scaling between [0,1]
            real(rp), intent(in) :: x
            real(rp), intent(in) :: wgt
            integer, intent(in) :: n
            
            real(rp) :: val
            val = 0.5*(wgt*(2*x-1) + (1-wgt)*(2*x-1)**n + 1)
        end function scaleFunc

    end subroutine raijuGenWarpSphGrid_Shafee2008


    subroutine raijuGenGridFromShGrid(raijuGrid, shGrid, iXML, opts, doVerboseO, okO, errMsgO)
        type(ShellGrid_T)  , intent(inout) :: raijuGrid
        type(ShellGrid_T), intent(in) :: shGrid
        type(XML_Input_T), intent(in) :: iXML
        type(raijuOptions_T), intent(in) :: opts
        logical, optional, intent(in) :: doVerboseO
        logical, optional, intent(out) :: okO
        character(len=*), optional, intent(out) :: errMsgO
        
        real(rp) :: thetaL, thetaU
        integer :: Ng
        integer, dimension(4) :: nGhosts
        integer :: iStart, iEnd
        logical :: doVerbose

        if(present(doVerboseO)) then
            doVerbose=doVerboseO
        else
            doVerbose=.false.
        endif

        if (present(okO)) okO = .false.
        if (present(errMsgO)) errMsgO = ""

        if(opts%thetaL%isSet()) then
            thetaL = opts%thetaL%get()
        else
            call iXML%Set_Val(thetaL , "grid/ThetaL", 15.0)
                !! Lower colat boundary [deg], ~15 Re in dipole
        endif
        if(opts%thetaU%isSet()) then
            thetaU = opts%thetaU%get()
        else
        call iXML%Set_Val(thetaU , "grid/ThetaU", 45.0)
            !! Upper colat boundary [deg], 2 Re in dipole
        endif

        call iXML%Set_Val(Ng, "grid/Ng", 4  )  ! Number of ghosts, in every direction for now
        nGhosts = Ng
        
        ! Convert to radians
        thetaL = thetaL*deg2rad
        thetaU = thetaU*deg2rad

        associate (Nt=>shGrid%Nt, Np=>shGrid%Np)
        ! Do some checks first
        if (thetaL < shGrid%th(1)) then
            call failGridGen("Error in raijuGenGridFromShGrid: lower theta bound outside of parent grid")
            return
        endif
        if (thetaU > shGrid%th(Nt+1)) then
            call failGridGen("Error in raijuGenGridFromShGrid: upper theta bound outside of parent grid")
            return
        endif

        ! Determine what our start and end i bounds should be within parent grid
        iStart = shGrid%is
        do while (shGrid%th(iStart + 1) < thetaL)
            iStart = iStart + 1
        enddo
        iEnd = iStart
        do while (iEnd < Nt .and. shGrid%th(iEnd+1) < thetaU - TINY)
            iEnd = iEnd + 1
        enddo

        if (iStart >= iEnd) then
            call failGridGen("Error in raijuGenGridFromShGrid: no valid active theta band remains after applying requested theta bounds")
            return
        endif

        if (iEnd + Ng + 1 > Nt + 1) then
            call failGridGen("Error in raijuGenGridFromShGrid: requested low-lat coverage cannot be satisfied with current parent grid and ghost count")
            return
        endif

        if (shGrid%th(iEnd+Ng+1) > PI/2.0_rp - TINY) then
            call failGridGen("Error in raijuGenGridFromShGrid: requested low-lat coverage would force RAIJU ghost cells to reach or cross the equator")
            return
        endif

        ! Now we can make our grid
        call GenChildShellGrid(shGrid, raijuGrid, RAI_SG_NAME, nGhosts, sub_is=iStart, sub_ie=iEnd, ghostsOnParentGridO=.true.)
        call validateRaijuChildShGrid(raijuGrid, thetaL=thetaL, thetaU=thetaU, okO=okO, errMsgO=errMsgO)
        if (present(okO)) then
            if (.not. okO) return
        else if (present(errMsgO)) then
            if (len_trim(errMsgO) > 0) return
        endif
        if (doVerbose) then
            write(*,'(a,f5.2,a,f5.2,a,f5.2,a,f5.2,a)')"RAIJU active theta range [deg] (rad) = [",raijuGrid%th(raijuGrid%is )*rad2deg,",",raijuGrid%th(raijuGrid%ie +1)*rad2deg,"] (",raijuGrid%th(raijuGrid%is ),",",raijuGrid%th(raijuGrid%ie +1),")"
            write(*,'(a,f5.2,a,f5.2,a,f5.2,a,f5.2,a)')"RAIJU full   theta range [deg] (rad) = [",raijuGrid%th(raijuGrid%isg)*rad2deg,",",raijuGrid%th(raijuGrid%ieg+1)*rad2deg,"] (",raijuGrid%th(raijuGrid%isg),",",raijuGrid%th(raijuGrid%ieg+1),")"
        endif
        if (present(okO)) okO = .true.
        end associate

        contains

        subroutine failGridGen(msg)
            character(len=*), intent(in) :: msg

            if (present(okO) .or. present(errMsgO)) then
                if (present(okO)) okO = .false.
                if (present(errMsgO)) errMsgO = trim(msg)
            else
                write(*,*) trim(msg)
                write(*,*) " thetaL = ",thetaL,", thetaU = ",thetaU
                write(*,*) " parent theta range = ",shGrid%th(1), shGrid%th(shGrid%Nt+1)
                stop
            endif
        end subroutine failGridGen

    end subroutine raijuGenGridFromShGrid


    subroutine validateRaijuChildShGrid(raijuGrid, thetaL, thetaU, okO, errMsgO)
        type(ShellGrid_T), intent(in) :: raijuGrid
        real(rp), optional, intent(in) :: thetaL, thetaU
        logical, optional, intent(out) :: okO
        character(len=*), optional, intent(out) :: errMsgO

        if (present(okO)) okO = .false.
        if (present(errMsgO)) errMsgO = ""

        if (.not. raijuGrid%isChild) then
            call failValidation("ERROR in validateRaijuChildShGrid: expected RAIJU ShellGrid to be a child grid")
            return
        endif
        if (.not. raijuGrid%isPeriodic) then
            call failValidation("ERROR in validateRaijuChildShGrid: expected RAIJU ShellGrid to wrap full phi")
            return
        endif
        if (any(raijuGrid%th(:) > PI/2.0_rp - TINY)) then
            call failValidation("ERROR in validateRaijuChildShGrid: RAIJU ghost cells reach or cross the equator")
            return
        endif
        if (present(thetaL)) then
            if (raijuGrid%th(raijuGrid%is) > thetaL + TINY .or. raijuGrid%th(raijuGrid%is+1) < thetaL - TINY) then
                call failValidation("ERROR in validateRaijuChildShGrid: first active theta cell does not bound requested thetaL")
                return
            endif
        endif
        if (present(thetaU)) then
            if (raijuGrid%th(raijuGrid%ie) > thetaU + TINY .or. raijuGrid%th(raijuGrid%ie+1) < thetaU - TINY) then
                call failValidation("ERROR in validateRaijuChildShGrid: last active theta cell does not bound requested thetaU")
                return
            endif
        endif

        if (present(okO)) okO = .true.

        contains

        subroutine failValidation(msg)
            character(len=*), intent(in) :: msg

            if (present(okO) .or. present(errMsgO)) then
                if (present(okO)) okO = .false.
                if (present(errMsgO)) errMsgO = trim(msg)
            else
                write(*,*) trim(msg)
                stop
            endif
        end subroutine failValidation
    end subroutine validateRaijuChildShGrid


    subroutine finalizeLLGrid(Grid, planet)
        !! Use a fully-created shell grid to allocate and populate the rest of the grid parameters
        !! NOTE: All lengths are either in units of meters or plentary radii, as specified below
        !!  no units of ionospheric radius (even though many of these parameters "live" on the ionospheric grid)
        type(raijuGrid_T), intent(inout) :: Grid
        type(planet_T), intent(in) :: planet

        integer :: i,j
        real(rp) :: L, thc_extra
        real(rp), dimension(3) :: xyz
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg) :: cosThc, BMagTh
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1) :: cosTh

        associate(shGr=>Grid%shGrid)
            ! First allocate remaining arrays
            allocate( Grid%thRp(shGr%isg:shGr%ieg+1) )
            allocate( Grid%thcRp(shGr%isg:shGr%ieg  ) )
            allocate( Grid%delTh(shGr%isg:shGr%ieg+1) )
            allocate( Grid%delPh(shGr%jsg:shGr%jeg+1) )
            allocate( Grid%areaCC  (shGr%isg:shGr%ieg  ,shGr%jsg:shGr%jeg) )
            allocate( Grid%areaFace(shGr%isg:shGr%ieg+1,shGr%jsg:shGr%jeg+1, 2) )
            allocate( Grid%lenFace (shGr%isg:shGr%ieg+1,shGr%jsg:shGr%jeg+1, 2) )
            allocate( Grid%Bmag  (shGr%isg:shGr%ieg,shGr%jsg:shGr%jeg) )
            allocate( Grid%cosdip(shGr%isg:shGr%ieg,shGr%jsg:shGr%jeg) )
            allocate( Grid%BrFace(shGr%isg:shGr%ieg+1,shGr%jsg:shGr%jeg+1, 2) )
            allocate( Grid%Brcc(shGr%isg:shGr%ieg,shGr%jsg:shGr%jeg) )
            

            ! For each theta value (defined on spherical grid in the ionosphere), 
            ! calculate corresponding theta of field line mapped to planet's surface
            xyz = 0.0
            do i=shGr%isg, shGr%ieg+1
                ! Don't need longitude information so just assume we are at Y=0
                xyz(1) = planet%ri_m/planet%rp_m*sin(shGr%th(i))
                xyz(3) = planet%ri_m/planet%rp_m*cos(shGr%th(i))
                L = DipoleL(xyz)
                Grid%thRp(i) = abs(asin(sqrt(1.0_rp/L)))
            enddo
            do i=shGr%isg, shGr%ieg
                xyz(1) = planet%ri_m/planet%rp_m*sin(shGr%thc(i))
                xyz(3) = planet%ri_m/planet%rp_m*cos(shGr%thc(i))
                L = DipoleL(xyz)
                Grid%thcRp(i) = abs(asin(sqrt(1.0_rp/L)))
            enddo

            ! Calculate theta/phi delta between cells [radians]
            ! First and last elements should be zero since there's no cells below/above isg/ieg
            Grid%delTh(shGr%isg) = 0.0
            Grid%delTh(shGr%ieg+1) = 0.0  
            do i=shGr%isg+1,shGr%ieg
                ! Define delta between th(i-1) and th(i)
                Grid%delTh(i) = abs(shGr%thc(i) - shGr%thc(i-1))
            enddo

            ! [radians]
            Grid%delPh(shGr%jsg) = 0.0
            Grid%delPh(shGr%jeg+1) = 0.0  
            do j=shGr%jsg+1,shGr%jeg
                ! Define delta between ph(i-1) and ph(i)
                Grid%delPh(j) = abs(shGr%phc(j) - shGr%phc(j-1))
            enddo

            ! Calc areas at cell centers [Rp^2]
            do i=shGr%isg,shGr%ieg
                do j=shGr%jsg,shGr%jeg
                    ! r^2 * (cos(th1)-cos(th2)) * dPh
                    Grid%areaCC(i,j) = (planet%ri_m/planet%rp_m)**2 &
                                     * ( cos(shGr%th(i)) - cos(shGr%th(i+1)) ) &
                                     * (shGr%ph(j+1) - shGr%ph(j))
                enddo
            enddo

            ! Area at faces [Rp^2]
            ! Kinda overkill, but just in case: use 8-centered reconstruction for each face
            ! Only do for faces of active cells, because we don't need them in ghosts
            Grid%areaFace = 0.0
            do i=shGr%is,shGr%ie+1
                do j=shGr%js,shGr%je+1
                    ! Theta dir
                    Grid%areaFace(i,j,RAI_TH) = Central8(Grid%areaCC(i-4:i+3, j))
                    ! Phi dir
                    Grid%areaFace(i,j,RAI_PH) = Central8(Grid%areaCC(i, j-4:j+3))
                enddo
            enddo


            ! Faces in theta dir are i +/- 1/2 faces, each of which spans the phi direction
            do j=shGr%jsg,shGr%jeg
                do i=shGr%isg,shGr%ieg+1
                    Grid%lenFace(i, j, RAI_TH) = (planet%ri_m/planet%rp_m) &
                                               * sin(shGr%th(i)) &
                                               * (shGr%ph(j+1) - shGr%ph(j))
                enddo
            enddo
            do j=shGr%jsg,shGr%jeg+1
                do i=shGr%isg,shGr%ieg
                    Grid%lenFace(i,j, RAI_PH) = (planet%ri_m/planet%rp_m) &
                                              * (shGr%th(i+1) - shGr%th(i))
                enddo
            enddo


            ! Calc ionospheric Bmag [nT] and cos of dip angle at cell centers
            cosThc = cos(shGr%thc)
            BMagTh = planet%magMoment*G2nT &
                    /(planet%ri_m/planet%rp_m)**3.0 &
                    * sqrt(1.0+3.0*cosThc**2.0)  ! [nT]
            do j=shGr%jsg,shGr%jeg
                Grid%Bmag(:,j) = BMagTh
                Grid%cosdip(:,j) = 2.0*cosThc/sqrt(1.0 + 3.0*cosThc**2.0) 
            enddo

            
            ! Radial component of Bmag [nT] at cell edges, used for velocity calculation
            Grid%BrFace = 0.0
            ! Theta-dir faces have theta locations at cell interfaces
            cosTh = cos(shGr%th)  ! Ni+1
            do j=shGr%jsg,shGr%jeg+1
                Grid%BrFace(:,j,RAI_TH) = planet%magMoment*G2nT &
                                        /(planet%ri_m/planet%rp_m)**3.0 &
                                        * 2*cosTh  ! [nT]
            enddo
            ! Phi-dir faces have theta locations at cell centers
            cosThc = cos(shGr%thc)  ! Ni
            do j=shGr%jsg,shGr%jeg+1
                Grid%BrFace(shGr%isg:shGr%ieg,j,RAI_PH) &
                                        = planet%magMoment*G2nT &
                                        /(planet%ri_m/planet%rp_m)**3.0 &
                                        * 2*cosThc  ! [nT]
            enddo
            ! There's a row of ieg+1 phi interfaces we are never gonna use but might as well mark them anyways
            ! So we don't get nans when dividing by BrFace when we calculate velocities
            ! Linearly extrapolate from last theta edge + difference of theta edge to theta center
            thc_extra = 2.0*shGr%th(shGr%ieg+1) - shGr%thc(shGr%ieg)
            Grid%BrFace(shGr%ieg+1,:,RAI_PH) &
                                       = planet%magMoment*G2nT &
                                       /(planet%ri_m/planet%rp_m)**3.0 &
                                       * 2*cos(thc_extra)  ! [nT]
            ! Also do Br for cell centers, can re-use cosTh
            do j=shGr%jsg,shGr%jeg
                Grid%Brcc(:,j) = planet%magMoment*G2nT &
                               /(planet%ri_m/planet%rp_m)**3.0 &
                               * 2*cosThc  ! [nT]
            enddo

        end associate

    end subroutine finalizeLLGrid


!------
! Grid helpers
!------

    function checkResGrid(sh, shR) result(isSame)
        type(ShellGrid_T), intent(in) :: sh
            !! Raiju's currently loaded grid
        type(ShellGrid_T), intent(in) :: shR
            !! ShellGrid read from restart file
        logical :: isSame
            !! Whether or not we things these are equivalent

        if ( (sh%Nt .ne. shR%Nt) &
        .or. (sh%Np .ne. shR%Np) &
        .or. (abs(sh%minGTheta - shR%minGTheta) > TINY) &
        .or. (abs(sh%maxGTheta - shR%maxGTheta) > TINY) ) then
            isSame = .false.
        else
            isSame = .true.
        endif

    end function checkResGrid

!------
! Lambda grid stuff
!------

!! Maybe should leave just spatial grid stuff in raijuGrids and move lambda stuff to a lambdaUtils

    subroutine initLambdaGrid(Model, Grid, configfname)
        type(raijuModel_T)  , intent(inout) :: Model
        type(raijuGrid_T), intent(inout) :: Grid
        character(len=strLen), intent(in) :: configfname

        ! First read in species information from config file
        !if (.not. Model%isRestart) then
        !    call populateSpeciesFromConfig(Model, Grid, configfname)
        !else
        !    call populateSpeciesFromConfig(Model, Grid, Model%ResF)
        !endif
        call populateSpeciesFromConfig(Model, Grid, configfname)

        ! Now prepare our alamc grid, tell each Species what its range is in alamc array
        call initAlamc(Grid)
        
    end subroutine initLambdaGrid



!------
! Spatial grid helpers
!------

    function wrapJIdx(sh, j, locEnum) result (jWrap)
        !! Takes j index and returns with wrapping applied
        !! Only meant for ghost cells, assumes jsg <= j <= jeg (+1 for corner centered)
        type(ShellGrid_T), intent(in) :: sh
        integer, intent(in) :: j
        integer, intent(in) :: locEnum
            !! SHGR_CC,SHGR_CORNER,SHGR_FACE_THETA,SHGR_FACE_PHI

        integer :: jWrap, jExtra

        if (locEnum == SHGR_CC .or. locEnum == SHGR_FACE_THETA) then
            jExtra = 0
        else if (locEnum == SHGR_CORNER .or. locEnum == SHGR_FACE_PHI) then
            jExtra = 1
        else
            write(*,*) "ERROR in wrapJIdx: locEnum provided is invalid"
            write(*,*) "locEnum:",locEnum
            stop
        endif

        if (j < sh%jsg .or. j > sh%jeg+jExtra) then
            write(*,*) "ERROR in wrapJIdx: j is outside of grid"
            write(*,*) "jsg,jeg+extra,j",sh%jsg, sh%jeg+jExtra, j
            stop
        endif

        if (j < sh%js) then
            jWrap = j + sh%Np
        elseif (j > sh%je) then
            jWrap = j - sh%Np
        else
            jWrap = j
        endif
    end function wrapJIdx


!------
! Spatial grid overloads
!------
    subroutine wrapJcc_real(sh, Q)
        !! Populates overlapping j ghosts with active cell-centered data
        !! Assumes that all values within (:, js:je) are valid to wrap
        !! And assumes that j ghosts overlap with last active cells on the other side
        type(ShellGrid_T), intent(in) :: sh
        real(rp), dimension(sh%isg:sh%ieg, sh%jsg:sh%jeg), intent(inout) :: Q

        ! Starting ghost cells
        Q(:, sh%jsg:sh%js-1) = Q(:, sh%je-sh%Nge+1:sh%je)
        ! Ending ghosts cells
        Q(:, sh%je+1:sh%jeg) = Q(:, sh%js:sh%js+sh%Ngw-1)

    end subroutine wrapJcc_real


    subroutine wrapJcc_int(sh, Q)
        !! Populates overlapping j ghosts with active cell-centered data
        !! Assumes that all values within (:, js:je) are valid to wrap
        !! And assumes that j ghosts overlap with last active cells on the other side
        type(ShellGrid_T), intent(in) :: sh
        integer, dimension(sh%isg:sh%ieg, sh%jsg:sh%jeg), intent(inout) :: Q

        ! Starting ghost cells
        Q(:, sh%jsg:sh%js-1) = Q(:, sh%je-sh%Ngw+1:sh%je)
        ! Ending ghosts cells
        Q(:, sh%je+1:sh%jeg) = Q(:, sh%js:sh%js+sh%Nge-1)

    end subroutine wrapJcc_int


    subroutine wrapJcc_logical(sh, Q)
        !! Populates overlapping j ghosts with active cell-centered data
        !! Assumes that all values within (:, js:je) are valid to wrap
        !! And assumes that j ghosts overlap with last active cells on the other side
        type(ShellGrid_T), intent(in) :: sh
        logical, dimension(sh%isg:sh%ieg, sh%jsg:sh%jeg), intent(inout) :: Q

        ! Starting ghost cells
        Q(:, sh%jsg:sh%js-1) = Q(:, sh%je-sh%Nge+1:sh%je)
        ! Ending ghosts cells
        Q(:, sh%je+1:sh%jeg) = Q(:, sh%js:sh%js+sh%Ngw-1)

    end subroutine wrapJcc_logical


    subroutine wrapJcorners_real(sh, Q)
        !! Populates overlapping j ghosts with active corner data
        !! Assumes that all values within (:, js:je) are valid to wrap
        !! And assumes that j ghosts overlap with last active cells on the other side
        type(ShellGrid_T), intent(in) :: sh
        real(rp), dimension(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1), intent(inout) :: Q

        ! Starting ghost cells
        Q(:, sh%jsg:sh%js) = Q(:, sh%je-sh%Nge+1:sh%je+1)
        ! Ending ghosts cells
        Q(:, sh%je+1:sh%jeg+1) = Q(:, sh%js:sh%js+sh%Ngw)

    end subroutine wrapJcorners_real


    subroutine wrapJcorners_int(sh, Q)
        !! Populates overlapping j ghosts with active corner data
        !! Assumes that all values within (:, js:je) are valid to wrap
        !! And assumes that j ghosts overlap with last active cells on the other side
        type(ShellGrid_T), intent(in) :: sh
        integer, dimension(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1), intent(inout) :: Q

        ! Starting ghost cells
        Q(:, sh%jsg:sh%js) = Q(:, sh%je-sh%Nge+1:sh%je+1)
        ! Ending ghosts cells
        Q(:, sh%je+1:sh%jeg+1) = Q(:, sh%js:sh%js+sh%Ngw)

    end subroutine wrapJcorners_int


    subroutine wrapJcorners_logical(sh, Q)
        !! Populates overlapping j ghosts with active corner data
        !! Assumes that all values within (:, js:je) are valid to wrap
        !! And assumes that j ghosts overlap with last active cells on the other side
        type(ShellGrid_T), intent(in) :: sh
        logical, dimension(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1), intent(inout) :: Q

        ! Starting ghost cells
        Q(:, sh%jsg:sh%js) = Q(:, sh%je-sh%Nge+1:sh%je+1)
        ! Ending ghosts cells
        Q(:, sh%je+1:sh%jeg+1) = Q(:, sh%js:sh%js+sh%Ngw)

    end subroutine wrapJcorners_logical

end module raijugrids
