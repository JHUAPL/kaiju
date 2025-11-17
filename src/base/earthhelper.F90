  ! various help functions for earth

module earthhelper
    use kdefs
    use planethelper

    implicit none
  

    !Magic numbers for Gallagher 2D
        ! the basic structure of these coefficients is as follows: 
        ! xn,yn are for nightside, and xd,yd dayside
        ! for the first index it is for kp=1,3 and 5
        ! for the second index:
        !       the numbers 1-6 are the parameters for the linear fit to
        !       the outer, middle and inner regions of the curve
        !       x numbers 7-8 are the x-centers where the curves are blended 
        !       y numbers 7-8 are the slopes of the blend

    real(rp), dimension(8), parameter, private :: xn1  = [ 8.35, 0.0 , 8.0 , 0.0 , 1.0 , 1.75, 5.85, 1.1]
    real(rp), dimension(8), parameter, private :: xn2  = [ 7.5 , 1.8 , 8.0 , 0.0 , 1.0 , 1.75, 4.1 , 1.1]
    real(rp), dimension(8), parameter, private :: xn3  = [ 7.5 , 1.8 , 8.0 , 0.0 , 1.0 , 1.75, 3.2 , 1.1]
    real(rp), dimension(8), parameter, private :: yn1  = [ 0.0 , 2.2 , 0.0 , 4.8 , 5.3 ,-1.0 , 0.1 , 0.1]
    real(rp), dimension(8), parameter, private :: yn2  = [-1.0 , 1.0 , 0.0 , 4.8 , 5.3 ,-1.0 , 0.1 , 0.1]
    real(rp), dimension(8), parameter, private :: yn3  = [-1.0 , 1.0 , 0.0 , 4.8 , 5.3 ,-1.0 , 0.1 , 0.1]
    real(rp), dimension(8), parameter, private :: xd1  = [ 8.0 , 3.4 , 8.0 , 0.0 , 1.0 , 2.1 , 5.5 , 1.1]
    real(rp), dimension(8), parameter, private :: xd2  = [ 8.0 , 3.4 , 8.0 , 0.0 , 1.0 , 2.1,  4.5 , 1.1]
    real(rp), dimension(8), parameter, private :: xd3  = [ 8.0 , 3.4 , 8.0 , 0.0 , 1.0 , 2.1 , 4.0 , 1.1]
    real(rp), dimension(8), parameter, private :: yd1  = [ 0.3 , 1.8 , 0.0 , 4.9 , 5.8 ,-1.0 , 0.1 , 0.1]
    real(rp), dimension(8), parameter, private :: yd2  = [ 0.3 , 1.8 , 0.0 , 4.9 , 5.8 ,-1.0 , 0.1 , 0.1]
    real(rp), dimension(8), parameter, private :: yd3  = [ 0.3 , 1.8 , 0.0 , 4.9 , 5.8 ,-1.0 , 0.1 , 0.1]

    real(rp), parameter, private :: GLMin = 1.0
    real(rp), parameter, private :: GLMax = 8.0

    real(rp), private :: kpDefault = 1.0


    !Toy code for putting in a quiet time RC
    !See Liemohn 2003
    real(rp), private :: QTRC_P0  = 1.0 !Scaling parameter for pressure
    real(rp), private :: QTRC_Lpk = 4.0 !Location of peak
    real(rp), private :: QTRC_dL  = 0.625 !Width

    contains

!Adapted by Kareem from FRT's original code
    !----------------------------------------
    ! Simple fit to Gallagher et al plasmasphere model, based on figure 2 from:
    ! Side note: In the figure from the paper, the x-axis is reversed
    ! Gallagher, D., Craven, P., & Comfort, R. (2000). Global core plasma model, JGR, 105, 18819–18833.
    ! PLEASE NOTE: It returns log10(n/cc) as a function of x,y in the equatorial plane and kp (1,3,5)
    ! it blends the dayside and nightside values from the model as a linear interpolation in x 
    ! to produce a crude 2d fit
    !
    ! v 1.1 1/16 frt
    !----------------------------------------

    !Routine to change kp-default on the fly if necessary
    subroutine SetKp0(kp)
        real(rp), intent(in) :: kp

        !Can only be 1,5
        kpDefault = GallagherKp(kp)

    end subroutine SetKp0

    ! Routine to enforce kp to Gallagher's valid range [1,5]
    pure function GallagherKp(kp) result(validKp)
    real(rp), intent(in) :: kp
    real(rp) :: validKp

    if (kp < 1.0_rp) then
        validKp = 1.0_rp
    else if (kp > 5.0_rp) then
        validKp = 5.0_rp
    else
        validKp = kp
    endif
    end function GallagherKp

    !Return 2D gallagher density afa r,phi (rad)
    function GallagherRP(r,phi,kpO) result(D)
        real(rp), intent(in) :: r,phi
        real(rp), intent(in), optional :: kpO
        real(rp) :: D
        real(rp) :: x,y
        x = r*cos(phi)
        y = r*sin(phi)
        if (present(kpO)) then
            D = GallagherXY(x,y,kpO)
        else
            D = GallagherXY(x,y)
        endif

    end function GallagherRP

    !Return 2D gallagher density
    function GallagherXY(x,y,kpO) result(D)
        real(rp), intent(in) :: x,y
        real(rp), intent(in), optional :: kpO
        real(rp) :: D

        real(rp), dimension(8) :: xfn,yfn,xfd,yfd
        real(rp) :: nout,nin,ninner,r,nd,nn,nfin
        integer :: kp

        D = 0.0

        if (present(kpO)) then
            !Can only be 1,5
            kp = nint(GallagherKp(kpO))
        else
            kp = nint(kpDefault)
        endif

        select case(kp)
            case(1)
                xfn = xn1
                yfn = yn1
                xfd = xd1
                yfd = yd1
            case(2)
                xfn = 0.5*( xn1 + xn2 ) 
                yfn = 0.5*( yn1 + yn2 ) 
                xfd = 0.5*( xd1 + xd2 ) 
                yfd = 0.5*( yd1 + yd2 ) 

            case(3)
                xfn = xn2
                yfn = yn2
                xfd = xd2
                yfd = yd2

            case(4)
                xfn = 0.5*( xn2 + xn3 ) 
                yfn = 0.5*( yn2 + yn3 ) 
                xfd = 0.5*( xd2 + xd3 ) 
                yfd = 0.5*( yd2 + yd3 ) 

            case(5)
                xfn = xn3
                yfn = yn3
                xfd = xd3
                yfd = yd3
        end select

        r = sqrt(x**2.0 + y**2.0)
        if ( (r>=GLMin) .and. (r<=GLMax) ) then
        !Night side
            nout   = linfunc(r,xfn(1),xfn(2),yfn(1),yfn(2))
            nin    = linfunc(r,xfn(3),xfn(4),yfn(3),yfn(4))
            ninner = linfunc(r,xfn(5),xfn(6),yfn(5),yfn(6))

            !Blend arrays
            nn = blend(r,nout,nin,xfn(7),yfn(7))
            !Now add inner spike
            nn = blend(r,nn,ninner,xfn(8),yfn(8))
        !Day side
            nout   = linfunc(r,xfd(1),xfd(2),yfd(1),yfd(2))
            nin    = linfunc(r,xfd(3),xfd(4),yfd(3),yfd(4))
            ninner = linfunc(r,xfd(5),xfd(6),yfd(5),yfd(6))

            !Blend arrays
            nd = blend(r,nout,nin,xfd(7),yfd(7))
            !Now add inner spike
            nd = blend(r,nd,ninner,xfd(8),yfd(8))

        !Mix day/night
            nfin = ((nd*(r+x)+(r-x)*nn)/(2*r))
            D = 10.0**nfin
        endif

        contains

        function linfunc(x,xi1,xi2,yi1,yi2)
            real(rp),intent(in) :: x,xi1,xi2,yi1,yi2
            real(rp) :: linfunc

            linfunc = yi1*(x-xi2)/(xi1-xi2)+yi2*(x-xi1)/(xi2-xi1)
        end function linfunc

        function tanhout(x,x0,h)
            real(rp), intent(in) :: x,x0,h
            real(rp) :: tanhout

            tanhout = 0.5*(tanh((x-x0)/h)+1.)
        end function tanhout

        function tanhin(x,x0,h)
            real(rp), intent(in) :: x,x0,h
            real(rp) :: tanhin
            tanhin = 0.5*(-tanh((x-x0)/h)+1.)
        end function tanhin

        ! blends values from functions inner and outer using tanh 
        ! centred at x0 with rate h
        function blend(x,youter,yinner,x0,h)
            real(rp), intent(in) :: x,youter,yinner,x0,h
            real(rp) :: blend

            blend = youter*0.5*(tanh((x-x0)/h)+1.)+yinner*0.5*(1.-tanh((x-x0)/h))
        end function blend

    end function GallagherXY

    function pwolf(r)
          implicit none
          real(rp), intent(in) :: r
          real(rp) :: pwolf
          ! fit coefficients
          real(rp), parameter :: rmax = 3.2015
          real(rp), parameter :: a1 = 1.1315
          real(rp), parameter :: b1 = 0.943
          real(rp), parameter :: c1 = 3.0
          real(rp), parameter :: a2 = 2.31
          real(rp), parameter :: b2 = 0.38
          real(rp), parameter :: d1 = 0.141
          real(rp), parameter :: d2 = 0.0412
          real(rp), parameter :: d3 = 3.0
          real(rp), parameter :: d4 = 0.3
                
          ! inner ring current pressure
          if(r <rmax)then
               pwolf = 10**(a1-b1*(r-c1)**2)
          else
               pwolf = 10**(a2-b2*r)
          end if
    ! outer pressure
          pwolf = pwolf + 10**(d1-d2*r)/(1.0+exp((d3-r)/d4))
    end function pwolf

    function psk(x,y,kpO)

    ! from the lui et al paper jgr 99, Jan 1, 1994, pp155 eqn 10.
    ! and from spence and kilveson jgr 98 sept 93 pp15490
    ! based on a spence and kivelson equation
    ! depending on the value of the kp input, it will use a
    ! pressure curve to match a kp 0, 3, or 5 level of disturbance
        INTEGER, intent(in), optional :: kpO
        REAL(rp), intent(in) :: x,y
        REAL(rp) :: p1,p2,a,b,rho,phi
        REAL(rp) :: factor = 1.0
        REAL(rp) :: Psk
        INTEGER :: kp

        if (present(kpO)) then
            kp = kpO
        else
            kp = kpDefault
        endif

        rho = sqrt(x*x+y*y); phi = - atan2(y,x)

        if (kp >= 5) then
           p1 = 1.9239e-7; p2 = 3.4201e-7; a = -0.8492; b = -2.2153
        else if (kp >= 2) then
           p1 = 3.1486e-7; p2 = 1.0802e-7; a = -0.9672; b = -1.9637
        else
           p1 = 2.9676e-7; p2 = 0.2958e-7; a = -0.9432; b = -1.7005
        endif
        psk = factor*(p1*exp(a*abs(rho))+p2*abs(rho)**(b))
        !Scale from P => nPa
        psk = (1.0e+9)*psk

    end function psk

    !Gallagher plasmasphere density model
    !Yoinked from Slava's code
    function psphD(L) result(D)
        ! approx based on gallagher et al, figure 1
        ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 105, NO. A8, PAGES 18,819-18,833, AUGUST 1, 2000
        ! returns values in ples/cc

        implicit none
        real(rp),intent(in) :: L
        real(rp)  :: D
        real(rp), parameter :: L0 = 4.5
        real(rp), parameter :: alpha = 10.
        real(rp), parameter :: a1 = -0.25
        real(rp), parameter :: b1 = 2.4
        real(rp), parameter :: a2 = -0.5
        real(rp), parameter :: b2 = 4.5
        real ::f,q

        f = 0.5*(1.0+tanh(alpha*(L-L0)))
        q = f*(a1*L + b1) + (1.0-f)*(a2*L + b2)
        D = 10.**q
    end function psphD
    
!---------
!Code for quiet time RC
    !Set parameters based on desired dst, Lpeak and dL
    subroutine SetQTRC(dst0,LpkO,dLO,doVerbO)
        real(rp), intent(in) :: dst0
        real(rp), intent(in), optional :: LpkO,dLO
        logical , intent(in), optional :: doVerbO

        logical :: doVerb
        real(rp) :: K,dB
    
        if (present(LpkO)) then
            QTRC_Lpk = LpkO
        endif

        if (present(dLO)) then
            QTRC_dL = dLO
        endif

        if (present(doVerbO)) then
            doVerb = doVerbO
        else
            doVerb = .false.
        endif

        if (dst0 >= 0) then
            QTRC_P0 = 0.0
            return
        else
            QTRC_P0 = 1.0
        endif

        !Get energy content w/ P0=1
        K = KIn()
        dB = -4.2*(1.0e-30)*K !Dst from DPS

        !Rescale to get target dst0
        QTRC_P0 = dst0/dB 
        
        !Calculate new energy content
        K = KIn()

        if (doVerb) then
            write(*,*) '---------------'
            write(*,*) 'Adding quiet-time ring current'
            write(*,*) 'Target Dst [nT]    = ', dst0
            write(*,*) 'RC Energy  [keV]   = ', K
            write(*,*) 'L-Peak     [Re]    = ', QTRC_Lpk
            write(*,*) 'dL         [Re]    = ', QTRC_dL
            write(*,*) 'P-Peak     [nPa]   = ', QTRC_P0
            write(*,*) '---------------'
        endif

        contains
            !Integrate total ring current energy content (keV)
            function KIn()
                real(rp) :: KIn

                real(rp) :: LIn,LOut,dl,K
                real(rp) :: Li,Lip,Lc,Pc,dV
                integer :: n,Nl
                !Now just doing lazy numerical integral
                LIn  = 1.5
                LOut = 10.0
                Nl = 100

                dl = (LOut-LIn)/Nl

                KIn = 0.0 !Cumulative energy
                do n=1,Nl
                    Li  = (n-1)*dl + LIn
                    Lip = Li+dl
                    Lc = Li+0.5*dl

                    Pc = P_QTRC(Lc)
                    dV = MagDV(Li,Lip)
                    KIn = KIn + Pc*dV
                enddo
                !K has units nPa*m3
                !Convert to keV
                KIn = ((1.0e-9)/kev2J)*KIn           
            end function KIn

    end subroutine SetQTRC

    !Pressure profile, See Liemohn 2003
    function P_QTRC(L) result(P)
        real(rp), intent(in) :: L
        real(rp) :: P
        real(rp) :: Lm,A,x,eA
        !Set values
        Lm = QTRC_Lpk - QTRC_dL
        
        A = QTRC_dL*exp(-1.0 + QTRC_Lpk/QTRC_dL)

        if (L>=Lm) then
            x = L-Lm
            eA = -(x-QTRC_Lpk)/QTRC_dL
            P = x*exp(eA)/A
        else
            P = 0.0
        endif

        P = QTRC_P0*P

    end function P_QTRC

    !Just return peak L value
    function LPk_QTRC() result(LPk)
        real(rp) :: LPk

        LPk = QTRC_Lpk
    end function LPk_QTRC
    
end module earthhelper
