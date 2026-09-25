module streamutils
    use chmpdefs
    use ebtypes
    use math
    use ebinterp
    use xml_input

    implicit none

    integer , parameter, private :: NumRKF45 = 6
    integer , parameter, private :: NumRKF23 = 3
    integer , parameter, private :: NumBS23  = 4
    real(rp), parameter, private :: StreamTol = 1.0e-3 !Units of planetary radius
    real(rp), dimension(NumRKF45), parameter, private :: & !Coefficients for RKF45
                    RKF45_LO = [25.0/216,0.0,1408.0/2565 ,2197.0/4104  ,-1.0/5 ,0.0   ], &
                    RKF45_HO = [16.0/216,0.0,6656.0/12825,28561.0/56430,-9.0/50,2.0/55]

    !Holder for data defining point on grid
    type GridPoint_T
        real(rp) :: xyz(NDIM)
        real(rp) :: t,dl !Time/lengthscale of current cell
        integer :: ijkG(NDIM) !Guess for cell location
    end type GridPoint_T

    !Generic one-step streamline routine
    !NOTE: TO handle FSAL methods, possible we get first B and output B at final position
    abstract interface
        subroutine OneStep_T(gpt,Model,ebState,eps,h,dx,iB,oB)
            import :: rp,NDIM,GridPoint_T,chmpModel_T,ebState_T
            type(GridPoint_T), intent(inout) :: gpt
            type(chmpModel_T), intent(in)    :: Model
            type(ebState_T)  , intent(in)    :: ebState
            real(rp), intent(in) :: eps
            real(rp), intent(inout) :: h,dx(NDIM)
            real(rp), intent(in) , dimension(NDIM), optional :: iB
            real(rp), intent(out), dimension(NDIM), optional :: oB
        end subroutine OneStep_T
    end interface

    procedure(OneStep_T), pointer :: StreamStep=>Step_MAGE

    !Some knobs for tracing cut-off
    real(rp), private :: bMinC = 0.0_rp !Min allowable field strength (in chimp units)
    contains

    subroutine setStreamlineKnobs(Model,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: inpXML

        character(len=strLen) :: sStr
        real(rp) :: bMin_nT

        if (Model%isMAGE) then
            call inpXML%Set_Val(sStr,'streamline/steptype',"MAGE")
        else
            call inpXML%Set_Val(sStr,'streamline/steptype',"RK4L")
        endif

        select case(trim(toUpper(sStr)))
            case("MAGE")
                StreamStep=>Step_MAGE
            case("RK4L")
                StreamStep=>Step_RK4L
        end select

        if (Model%isMAGE) then
            call inpXML%Set_Val(bMin_nT,"streamline/bMin_C",TINY)
        else
            bMin_nT = 0.0 !Don't do this for non-mage case
        endif

        !Convert bmin from nT to chimp eb coordinates
        bMinC = bMin_nT/oBScl

    end subroutine setStreamlineKnobs

    !Combines speedy method away from the axis w/ more careful stepping nearby
    subroutine Step_MAGE(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB

        integer :: Nr,j0
        logical :: isAxS,isAxE

        if ( (.not. present(iB)) .or. (.not. present(oB)) ) then
            write(*,*) "Non-optional error in step_mage"
        endif

        !Test field magnitude
        if (norm2(iB) <= bMinC) then
            !Get outta here
            dx = 0.0
            h  = 0.0
            if (present(oB)) oB = 0.0
            return
        endif

        j0 = gpt%ijkG(JDIR)
        ! !Get approx number of rings, 4/8/12/16 (DQOH)
        ! Nr = nint( 4*( log(1.0*ebState%ebGr%Nkp/64.0)/log(2.0) + 1) )
        Nr = 2 
        isAxS = (j0 < ebState%ebGr%js+Nr)
        isAxE = (j0 > ebState%ebGr%je-Nr)

        if (isAxS .or. isAxE) then
            call Step_RKF45(gpt,Model,ebState,eps,h,dx,iB,oB)
        else
            call Step_RK4L (gpt,Model,ebState,eps,h,dx,iB,oB)
        endif
    end subroutine Step_MAGE

    subroutine Step_FE(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB

        real(rp), dimension(NDIM) :: x0
        logical :: isGood

        x0 = gpt%xyz
        if (present(iB)) then
            dx = h*normVec(iB)
            isGood = .true.
        else
            dx = h*FastHat(x0,gPt%t,Model,ebState,isGood,gpt%ijkG)
        endif

        if (present(oB)) then
            oB = FastMag(x0+dx,gPt%t,Model,ebState,isGood,gpt%ijkG)
        endif

    end subroutine Step_FE

    !RK4 w/ linearized substeps
    subroutine Step_RK4L(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB

        real(rp), dimension(NDIM) :: Jb,Jb2,Jb3,F1,F2,F3,F4
        real(rp), dimension(NDIM) :: x0,B,B1,B2,B3,B4
        real(rp), dimension(NDIM,NDIM) :: JacB
        real(rp) :: dsmag,dsOld,ds,sig,MagB,MagJb,Lb
        real(rp) :: ds2,ds3,ds4
        logical :: isGood

        x0 = gpt%xyz
        !Need Jacobian, using streamlined routine
        JacB = FastJacB(x0,gPt%t,Model,ebState,gPt%ijkG)

        if (present(iB)) then
            B = iB
        else
            B = FastMag(x0,gPt%t,Model,ebState,isGood,gPt%ijkG)
        endif

        MagB = norm2(B)
        if (MagB<TINY .or. (.not. CheckMagmag(B))) then
            !Get outta here
            dx = 0.0
            h  = 0.0
            if (present(oB)) oB = 0.0
            return
        endif

        MagJb = norm2(JacB)
        Lb = MagB/max(MagJb,TINY)

    !Get step length
        !h is signed step, dsmag is absolute
        sig = sign(1.0_rp,h)
        dsOld = abs(h)
        !Pick ds = eps*Lb
        dsmag = eps*Lb
        !Now constrain ds to be:
        ! <= n x dsOld, dl
        ! >= eps*dl
        call ClampValue(dsmag,eps*gpt%dl,min(2*dsOld,gpt%dl))

        h = sig*dsmag !Step length
    !Do step
        ds = h/max(MagB,TINY) !Streamline units
        !Get powers of jacobian
        Jb  = matmul(JacB,B  )
        Jb2 = matmul(JacB,Jb )
        Jb3 = matmul(JacB,Jb2)
        !Get powers of ds
        ds2 = ds*ds
        ds3 = ds2*ds
        ds4 = ds3*ds
        !Calculate steps
        F1 = ds*B
        F2 = F1 + (ds2/2)*Jb
        F3 = F2 + (ds3/4)*Jb2
        F4 = ds*B + ds*ds*Jb + (ds**3.0/2.0)*Jb2 + (ds**4.0/4.0)*Jb3
        !F4 = 2*F3 - F1 + (ds4/4)*Jb3
        
        !Construct field values for testing
        B1 = F1/ds
        B2 = F2/ds
        B3 = F3/ds
        B4 = F4/ds

        isGood = CheckMagmag(B1) .and. CheckMagmag(B2) .and. &
                 CheckMagmag(B3) .and. CheckMagmag(B4)
        if (.not. isGood) then
            dx = 0.0
            h  = 0.0
            if (present(oB)) oB = 0.0
            !write(*,*) "Bad value in RK4L!"
            return
        endif

    !Advance
        dx = (F1+2*F2+2*F3+F4)/6.0
    !Finish up
        if (present(oB)) then
            oB = FastMag(x0+dx,gPt%t,Model,ebState,isGood,gpt%ijkG)
        endif

    end subroutine Step_RK4L

    !RK23 (RKF23) w/ linearized field
    subroutine StepRKF23Lin(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB

        real(rp), dimension(NDIM) :: x0,x2,x3,k1,k2,k3,dxLO,dxHO,B0
        real(rp), dimension(NDIM,NDIM) :: JacB0

        real(rp) :: ddx,sScl,absh
        logical :: isGood

        x0 = gpt%xyz
        !Need Jacobian, using streamlined routine
        JacB0 = FastJacB(x0,gPt%t,Model,ebState,gPt%ijkG)

        if (present(iB)) then
            B0 = iB
            isGood = .true.
        else
            B0 = FastMag(x0,gPt%t,Model,ebState,isGood,gPt%ijkG)
        endif
        !Get first hat
        k1 = h*normVec(B0)

        if (isGood) then
            x2 = x0 + k1
            k2 = h*LinHat(x0,B0,JacB0,x2)

            x3 = x0 + (k1+k2)/4.0
            k3 = h*LinHat(x0,B0,JacB0,x3)
        else
            !Bad step
            dx = 0.0
            h = 0.0
            if (present(oB)) oB = 0.0
            return
        endif

        !If we're still here, the step was good so let's calculate
        dxLO = 0.5*(k1+k2)
        dxHO = (k1+k2+4.0*k3)/6.0

        ddx = max( norm2(dxLO-dxHO),TINY )
        sScl = 0.7937*( (StreamTol*gpt%dl/ddx)**(1.0/3.0) ) !Relative
        
        !Now calculate new step
        absh = 0.8*abs(h)*sScl
        !Clamp min/max step
        call ClampValue(absh,eps*gpt%dl,min(0.5*gpt%dl,2*abs(h)))
        h = sign(absh,h)

        !Using HO step
        dx = dxHO
        if (present(oB)) then
            oB = FastMag(x0+dx,gPt%t,Model,ebState,isGood,gPt%ijkG)
        endif

        contains

            function LinHat(x0,B0,Jb0,x) result(bhat)
                real(rp), dimension(NDIM), intent(in) :: x0,B0,x
                real(rp), dimension(NDIM,NDIM), intent(in) :: Jb0

                real(rp), dimension(NDIM) :: bhat
                real(rp), dimension(NDIM) :: dx,B

                dx = x-x0
                B  = B0 + matmul(Jb0,dx)
                bhat = normVec(B)
                 
            end function LinHat
    end subroutine StepRKF23Lin

    subroutine Step_RKF45(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB

        real(rp), dimension(NDIM) :: x0,x2,x3,x4,x5,x6,k1,k2,k3,k4,k5,k6
        real(rp), dimension(NDIM) :: dxLO,dxHO
        real(rp) :: ddx,sScl,absh
        logical :: isGoods(NumRKF45),isG

        isGoods = .false.
        x0 = gpt%xyz

        if (present(iB)) then
            k1 = h*normVec(iB)
            isGoods(1) = .true.
        else
            k1 = h*FastHat(x0,gPt%t,Model,ebState,isGoods(1),gPt%ijkG)
        endif

        if (all(isGoods(1:1))) then
        	x2 = x0 + k1/4.0
        	k2 = h*FastHat(x2,gPt%t,Model,ebState,isGoods(2),gPt%ijkG)
        endif
        if (all(isGoods(1:2))) then
        	x3 = x0 + (3.0/32)*k1 + (9.0/32)*k2
        	k3 = h*FastHat(x3,gPt%t,Model,ebState,isGoods(3),gPt%ijkG)
        endif
        if (all(isGoods(1:3))) then
        	x4 = x0 + (1932.0/2197)*k1 - (7200.0/2197)*k2 + (7296.0/2197)*k3
        	k4 = h*FastHat(x4,gPt%t,Model,ebState,isGoods(4),gPt%ijkG)
        endif
        if (all(isGoods(1:4))) then
        	x5 = x0 + (439.0/216)*k1 - 8.0*k2 + (3680.0/513)*k3 - (845.0/4104)*k4
        	k5 = h*FastHat(x5,gPt%t,Model,ebState,isGoods(5),gPt%ijkG)
        endif
        if (all(isGoods(1:5))) then
        	x6 = x0 - (8.0/27)*k1 + 2.0*k2 - (3544.0/2565)*k3 + (1859.0/4104)*k4 - (11.0/40)*k5
        	k6 = h*FastHat(x6,gPt%t,Model,ebState,isGoods(6),gPt%ijkG)
        endif

        if (.not. all(isGoods)) then
            !At least one step was bad
            dx = 0.0
            h = 0.0
            if (present(oB)) oB = 0.0
            return
        endif

        !If we're still here, the step was good so let's calculate
        dxLO = RKF45_LO(1)*k1 + RKF45_LO(2)*k2 + RKF45_LO(3)*k3 + RKF45_LO(4)*k4 + RKF45_LO(5)*k5 + RKF45_LO(6)*k6
        dxHO = RKF45_HO(1)*k1 + RKF45_HO(2)*k2 + RKF45_HO(3)*k3 + RKF45_HO(4)*k4 + RKF45_HO(5)*k5 + RKF45_HO(6)*k6

        ddx = max( norm2(dxLO-dxHO),TINY )
        !sScl = 0.84*( (StreamTol/ddx)**(0.25) )
        sScl = 0.84*( (StreamTol*gpt%dl/ddx)**(0.25) ) !Relative

        !Now calculate new step
        absh = 0.95*abs(h)*sScl !Optimal value according to math
        !Clamp min/max step based on fraction of cell size and 3x old value
        call ClampValue(absh,eps*gpt%dl,min(3*abs(h),1.5*gpt%dl))
        h = sign(absh,h)

        !Using HO step (but see some comments about using LO for stiff problems)
        dx = dxHO
        if (present(oB)) then
            oB = FastMag(x0+dx,gPt%t,Model,ebState,isG,gPt%ijkG)
        endif

    end subroutine Step_RKF45

    subroutine Step_BS23(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB


        real(rp), dimension(NDIM) :: x0,x2,x3,x4,k1,k2,k3,k4
        real(rp), dimension(NDIM) :: dxLO,dxHO,B
        real(rp) :: ddx,sScl,absh
        logical :: isGoods(NumBS23),isG
        real(rp), parameter :: eStpT = 1.0e-3
        real(rp), parameter :: dlMax = 2.0 !Max cell length

        write(*,*) 'This method needs to be debugged!'
        isGoods = .false.
        x0 = gpt%xyz

        if (present(iB)) then
            k1 = h*normVec(iB)
            isGoods(1) = .true.
        else
            k1 = h*FastHat(x0,gPt%t,Model,ebState,isGoods(1),gPt%ijkG)
        endif
        if (all(isGoods(1:1))) then
            x2 = x0 + k1/2.0
            k2 = h*FastHat(x2,gPt%t,Model,ebState,isGoods(2),gPt%ijkG)
        endif
        if (all(isGoods(1:2))) then
            x3 = x0 + (3.0/4)*k2
            k3 = h*FastHat(x3,gPt%t,Model,ebState,isGoods(3),gPt%ijkG)
        endif
        
        if (all(isGoods(1:3))) then
            dxHO = (2*k1+3*k2+4*k3)/9.0
            x4 = x0 + dxHO !Actual solution
            B = FastMag(x4,gPt%t,Model,ebState,isGoods(4),gPt%ijkG)
            if (present(oB)) then
                oB = B
            endif
            k4 = h*normVec(B)
        endif
        !Check for failure
        if (.not. all(isGoods)) then
            !At least one step was bad
            dx = 0.0
            h = 0.0
            if (present(oB)) oB = 0.0
            return
        endif
    !Now finish up
        dxLO = (7*k1+6*k2+8*k3+3*k4)/24.0
        dx = dxHO
        
    ! !Kutta-Merson style
    !     ddx = norm2(dxHO-dxLO)/gpt%dl
    !     if (ddx >= eStpT) then
    !         h = h/2.0 !Reduce step
    !     else if (ddx <= eStpT/64.0) then
    !         h = 2.0*h !Increase step
    !     endif
    !     absh = abs(h)
    !     call ClampValue(absh,eps*gpt%dl,eMax*gpt%dl)
        

    !Embedded opt style
        ! ddx = max( norm2(dxLO-dxHO),TINY )
        ! !sScl = 0.9*0.7*( (StreamTol/ddx)**(0.5) )
        ! !sScl = 0.8*(StreamTol/ddx)**(0.25)

        ! !Now calculate new step
        ! absh = abs(h)*sScl !Optimal value according to math
        ! !Clamp min/max step based on fraction of cell size
        ! call ClampValue(absh,eps*gpt%dl,eMax*gpt%dl)

        h = sign(absh,h)

    end subroutine Step_BS23

    !Return magnetic field unit vector at given point/time
    function FastHat(xyz,t,Model,ebState,isIn,ijkG) result(bhat)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out) :: isIn
        integer, intent(inout) :: ijkG(NDIM)
        real(rp) :: bhat(NDIM)

        real(rp), dimension(NDIM) :: B

        B = FastMag(xyz,t,Model,ebState,isIn,ijkG)
        bhat = normVec(B)
        !Test field magnitude
        if (norm2(B) <= bMinC) then
            !Field too weak
            bhat = 0.0
            isIn = .false.
        endif
        
    end function FastHat

    !Check if a given field is sufficiently strong (above bMinC)
    function CheckMagmag(B) result(isGood)
        real(rp), intent(in) :: B(NDIM)
        logical :: isGood
        real(rp) :: bmag

        bmag = norm2(B)
        if (bmag <= bMinC) then
            isGood = .false.
        else
            isGood = .true.
        endif

    end function CheckMagmag

    function FastMag(xyz,t,Model,ebState,isIn,ijkG) result(B)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out) :: isIn
        integer, intent(inout) :: ijkG(NDIM)
        real(rp), dimension(NDIM) :: B
        integer :: ijkGuess(NDIM)

        ! The locator initializes its output before inspecting the guess.
        ! Keep the previous cell separate so the output cannot overwrite it.
        ijkGuess = ijkG
        call locate(xyz,ijkG,Model,ebState%ebGr,isIn,ijkGuess)
        B = fldInterp(xyz,t,Model,ebState,BFLD,isIn,ijkG)
    end function FastMag

    function FastJacB(xyz,t,Model,ebState,ijk) result(JacB)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        integer, intent(in) :: ijk(NDIM)
        real(rp), dimension(NDIM,NDIM) :: JacB

        type(gcFields_T) :: gcFields
        real(rp), dimension(NDIM) :: E,B

        call ebFields(xyz,t,Model,ebState,E,B,gcFields=gcFields)
        JacB = gcFields%JacB
    end function FastJacB

    ! !Return magnetic field vector at given point/time
    ! !NOTE: This reproduces a lot of ebinterp code but is designed to be streamlined/faster
    ! function FastMag(xyz,t,Model,ebState,isIn,ijkG) result(B)
    !     real(rp), intent(in) :: xyz(NDIM),t
    !     type(chmpModel_T), intent(in) :: Model
    !     type(ebState_T), intent(in)   :: ebState
    !     logical, intent(out) :: isIn
    !     integer, intent(inout) :: ijkG(NDIM)
    !     real(rp), dimension(NDIM) :: B

    !     real(rp), dimension(NDIM) :: B0
    !     real(rp), dimension(Nw,Nw,Nw) :: W
    !     real(rp), dimension(Nw,Nw,Nw,NDIM) :: Qb !Buffer
    !     integer :: n,i0,j0,k0

    !     real(rp), dimension(NDIM) :: Bold
    !     integer :: i1,i2,i3

    !     !Bold = fldInterp(xyz,t,Model,ebState,BFLD,isIn,ijkG)
    !     B = fldInterp(xyz,t,Model,ebState,BFLD,isIn,ijkG)

    !     return

    !     !Locate w/ guess
    !     !call locate(xyz,ijkG,Model,ebState%ebGr,isIn,ijkG)
    !     if (.not. isIn) then
    !         B = Model%B0(xyz)
    !         return
    !     endif

    ! !If still here, do work. NOTE: Assuming static field to avoid extra work
    !     B0 = Model%B0(xyz)
    !     call GetWeights(xyz,ijkG,W,Model,ebState%ebGr)
    !     i0 = ijkG(IDIR) ; j0 = ijkG(JDIR) ; k0 = ijkG(KDIR)

    !     Qb = ebState%eb1%dB(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,XDIR:ZDIR)

    !     ! do n=1,NDIM
    !     !     !NOTE: Avoiding doing below to avoid array temporary creation
    !     !     !B(n) = sum( W(:,:,:)*Qb(:,:,:,n) ) + B0(n)
    !     !     B(n) = B0(n) + FastColon(W,Qb(:,:,:,n))
    !     ! enddo

    !     !Bad!
    !     ! do n=1,NDIM
    !     !     B(n) = B0(n) + sum(W*Qb(:,:,:,n))
    !     ! enddo

    !     do n=1,NDIM
    !         !V(n) = V0(n) + w1*sum(W*v1b(:,:,:,n)) + w2*sum(W*v2b(:,:,:,n))
    !         B(n) = B0(n)
    !         do i3=1,Nw
    !             do i2=1,Nw
    !                 do i1=1,Nw
    !                     B(n) = B(n) + W(i1,i2,i3)*Qb(i1,i2,i3,n)
    !                 enddo
    !             enddo
    !         enddo
    !     enddo

    !     if (norm2(Bold-B) > 1.0e-5) then
    !         !$OMP CRITICAL
    !         write(*,*) '---'
    !         write(*,*) 'dB = ', norm2(Bold-B),Bold,B
    !         write(*,*) 'dB = ', Bold-B
    !         write(*,*) 'W = ', W
    !         do n=1,NDIM
    !             write(*,*) 'Qbn = ', Qb(:,:,:,n)
    !         enddo
    !         write(*,*) '---'
    !         !$OMP END CRITICAL
    !     endif
    ! end function FastMag

    ! !Fast tensor contraction, may want to toy w/ this for vectorization
    ! function FastColon(A,B) result(ab)
    !     real(rp), dimension(NDIM,NDIM), intent(in) :: A,B
    !     real(rp) :: ab

    !     ab = sum(A*B)
    ! end function FastColon

    !Return jacobian of B at given point/time
    !NOTE: This assumes ijk is CORRECT and xyz is indomain and eb is static
    !NOTE: This reproduces a lot of ebinterp code but is designed to be streamlined/faster
    ! recursive function FastJacDB(xyz,t,Model,ebState,ijk) result(JacDB)
    !     real(rp), intent(in) :: xyz(NDIM),t
    !     type(chmpModel_T), intent(in) :: Model
    !     type(ebState_T), intent(in)   :: ebState
    !     integer, intent(in) :: ijk(NDIM)
    !     real(rp), dimension(NDIM,NDIM) :: JacDB

    !     real(rp), dimension(NDIM,NDIM) :: Tix
    !     real(rp), dimension(NDIM)      :: ezp,wE,wZ,wP,wEp,wZp,wPp
    !     real(rp), dimension(Nw,Nw,Nw)  :: eW,zW,pW !Interpolation weights
    !     real(rp), dimension(Nw,Nw,Nw,NDIM) :: dB  !Interpolation stencils
    !     real(rp), dimension(Nw,Nw,Nw)  :: dBn !Holder
    !     integer :: i,j,k,n,m, i0,j0,k0

    !     !Variables for axis handling
    !     integer :: ip,jp,kp,ijkAx(NDIM)
    !     real(rp) :: wAx
    !     real(rp), dimension(NDIM) :: Xp,Xm,Xc
    !     logical :: isAxisS,isAxisE
    !     real(rp), dimension(NDIM,NDIM) :: pJacDB,mJacDB
    !     integer, parameter :: dAxI=1

    !     i0 = ijk(IDIR) ; j0 = ijk(JDIR) ; k0 = ijk(KDIR)

    ! !Handle axis
    !     isAxisS = .false.
    !     isAxisE = .false.
    !     if ( j0 <= ebState%ebGr%js+dAxI ) then
    !         isAxisS = .true.
    !     endif
    !     if ( j0 >= ebState%ebGr%je-dAxI ) then
    !         isAxisE = .true.
    !     endif

    !     if (isAxisS .or. isAxisE) then
    !         associate( ebGr=>ebState%ebGr )
    !         !Do different calculation at axis
    !         Xc = ebGr%xyzcc(i0,j0,k0,XDIR:ZDIR)
    !         if (isAxisS) then
    !         !Positive displacement
    !             ip = i0;jp = ebGr%js+dAxI+1;kp = k0
    !             Xp = ebGr%xyzcc(ip,jp,kp,XDIR:ZDIR)
    !             ijkAx = [ip,jp,kp]
    !             pJacDB = FastJacDB(Xp,t,Model,ebState,ijkAx)
    !         !Negative displacement
    !             call ijk2Active(Model,ebGr,i0,ebGr%js-2-dAxI,k0,ip,jp,kp)
    !             Xm = ebGr%xyzcc(ip,jp,kp,XDIR:ZDIR)
    !             ijkAx = [ip,jp,kp]
    !             mJacDB = FastJacDB(Xm,t,Model,ebState,ijkAx)
    !         !Weight for P (closer)
    !             wAx = norm2(Xm-Xc)/norm2(Xp-Xm)
    !         else !Negative axis
    !         !Positive displacement, flipping positive (to closer point)   
    !             ip = i0;jp = ebGr%je-1-dAxI;kp = k0
    !             Xp = ebGr%xyzcc(ip,jp,kp,XDIR:ZDIR)
    !             ijkAx = [ip,jp,kp]
    !             pJacDB = FastJacDB(Xp,t,Model,ebState,ijkAx)
    !         !Negative displacement
    !             call ijk2Active(Model,ebGr,i0,ebGr%je+2+dAxI,k0,ip,jp,kp)
    !             Xm = ebGr%xyzcc(ip,jp,kp,XDIR:ZDIR)
    !             ijkAx = [ip,jp,kp]
    !             mJacDB = FastJacDB(Xm,t,Model,ebState,ijkAx)
    !         !Weight for P (closer)
    !             wAx = norm2(Xm-Xc)/norm2(Xp-Xm)   
    !         endif

    !         !Interpolate dJac across axis (both cases)
    !         JacDB = wAx*pJacDB + (1-wAx)*mJacDB

    !         return !We're done here
    !         end associate
    !     endif !Axis

    ! !Pull stencil
    !     dB = ebState%eb1%dB(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,XDIR:ZDIR)
    ! !Get mapping and weights
    !     !Map to ezp
    !     ezp = Map2ezp(xyz,ijk,Model,ebState%ebGr)

    !     !Get 1D spatial weights
    !     wE = Wgt1D(ezp(IDIR))
    !     wZ = Wgt1D(ezp(JDIR))
    !     wP = Wgt1D(ezp(KDIR))

    ! !Now get JacDB
    !     !---------
    !     !Get 1D weights for Jacobians            
    !     wEp = Wgt1Dp(ezp(IDIR))
    !     wZp = Wgt1Dp(ezp(JDIR))
    !     wPp = Wgt1Dp(ezp(KDIR))

    !     !Turn 1D weights into 3D weights
    !     do k=1,Nw
    !         do j=1,Nw
    !             do i=1,Nw
    !                 !Partial derivatives of weights wrt eta,zeta,psi
    !                 eW(i,j,k) = wEp(i)*wZ (j)*wP (k)
    !                 zW(i,j,k) = wE (i)*wZp(j)*wP (k)
    !                 pW(i,j,k) = wE (i)*wZ (j)*wPp(k)
    !             enddo
    !         enddo
    !     enddo

    !     !Calculate Jacobians
    !     !JacA(i,j) = d B_Xi / dXj
    !     !Tix(i0,j0,k0,ezp,xyz) = ezp derivs wrt xyz

    !     !Pull metric terms
    !     Tix = ebState%ebGr%Tix(i0,j0,k0,:,:)

    !     !Do main calculation
    !     do m=1,NDIM !Derivative direction (x,y,z)
    !         do n=1,NDIM !Vector component
    !             dBn = dB(:,:,:,n)
    !             JacDB(n,m) =   Tix(IDIR,m)*FastColon(eW,dBn)  &
    !                          + Tix(JDIR,m)*FastColon(zW,dBn)  &
    !                          + Tix(KDIR,m)*FastColon(pW,dBn)
    !         enddo
    !     enddo

    ! end function FastJacDB

    function getDiag(ebGr,ijk) result (dl)
        type(ebGrid_T), intent(in)   :: ebGr
        integer, intent(in) :: ijk(NDIM)
        real(rp) :: dl
        integer :: i,j,k
        i = ijk(IDIR) ; j = ijk(JDIR) ; k = ijk(KDIR)

        dl = norm2( ebGr%xyz(i+1,j+1,k+1,:)-ebGr%xyz(i,j,k,:) )

    end function getDiag

    subroutine cleanLine(fL)
        type(magLine_T), intent(inout) :: fL

        fL%nMax = MaxFL
        fL%x0 = 0.0
        fL%Nm = 0
        fL%Np = 0
        fL%isGood = .false.

        if (allocated(fL%xyz )) deallocate(fL%xyz )
        if (allocated(fL%ijk )) deallocate(fL%ijk )
        if (allocated(fL%Bxyz)) deallocate(fL%Bxyz)
        if (allocated(fL%magB)) deallocate(fL%magB)
        if (allocated(fL%Gas )) deallocate(fL%Gas )

    end subroutine cleanLine


    subroutine deepCopyLine(flIn, flOut)
        type(magLine_T), intent(in) :: flIn
        type(magLine_T), intent(out) :: flOut

        integer :: N1,N2,Nv,Ns
        
        N1 = flIn%Nm
        N2 = flIn%Np
        Nv = size(flIn%Gas,dim=2)
        Ns = ubound(flIn%Gas,dim=3)

        call cleanLine(flOut)

        flOut%nMax = flIn%nMax
        flOut%x0 = flIn%x0
        flOut%Nm = flIn%Nm
        flOut%Np = flIn%Np
        flOut%isGood = flIn%isGood

        allocate(flOut%xyz (-N1:+N2,NDIM))
        allocate(flOut%ijk (-N1:+N2,NDIM))
        allocate(flOut%Bxyz(-N1:+N2,NDIM))
        allocate(flOut%magB(-N1:+N2))
        allocate(flOut%Gas(-N1:+N2,Nv,0:Ns))

        flOut%xyz = flIn%xyz
        flOut%ijk = flIn%ijk
        flOut%Bxyz = flIn%Bxyz
        flOut%magB = flIn%magB
        flOut%Gas  = flIn%Gas

    end subroutine deepCopyLine

end module streamutils
