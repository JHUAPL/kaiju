!Various low-level math routines
module math
    use kdefs

    implicit none

    !Matrix for 5x5 smoothing (unscaled)
    real(rp), private, parameter, dimension(5,5) ::  &
             SmoothOp5x5 = reshape( [ 1, 4, 7, 4, 1, &
                                      4,16,26,16, 4, &
                                      7,26,41,26, 7, &
                                      4,16,26,16, 4, &
                                      1, 4, 7, 4, 1  ], [5,5] )

    !Marix for 3x3 smoothing (scaled)
    real(rp), parameter, dimension(3,3) ::  &
              SmoothOpTSC = reshape( [ 0.0625,0.1250,0.0625, &
                                       0.1250,0.2500,0.1250, &
                                       0.0625,0.1250,0.0625  ], [3,3] )

    !3x3 identity tensor
    real(rp), parameter, dimension(3,3) ::  &
              Eye33 = reshape( [ 1.0,0.0,0.0, &
                                 0.0,1.0,0.0, &
                                 0.0,0.0,1.0  ], [3,3] )
    contains

    !Generates a random number between vMin/vMax
    function genRand(vMin, vMax, seedIN)
        real(rp), intent(in) :: vMin, vMax
        integer, intent(in), optional :: seedIN
        integer, allocatable :: seed(:)
        real(rp), allocatable :: rseed(:)
        integer :: n, i

        real(rp) :: genRand
        real(rp) :: r

        ! usfull for testing and comparsions of runs
        if (present(seedIN)) then
            call random_seed(size=n)
            if (.not. allocated(seed)) &
                allocate (seed(n), rseed(n))
            seed(:) = seedIN
            call random_seed(put=seed)
            call random_number(rseed)
            seed = 1e6*rseed
            call random_seed(put=seed)
        end if

        call random_number(r)
        genRand = vMin + (vMax - vMin)*r

    end function genRand

    !Generates log-spaced random numbers between vMin/vMax
    function genRandLog(vMin, vMax)
        real(rp), intent(in) :: vMin, vMax
        real(rp) :: genRandLog

        real(rp) :: lMin, lMax, lR

        !Find bounds in log-space
        lMin = log10(vMin)
        lMax = log10(vMax)

        !Generate uniform random # in the log interval
        lR = genRand(lMin, lMax)

        !Un-log number
        genRandLog = 10**lR

    end function genRandLog

    !Mollifier (smooth bump function, compactly supported)
    !Exp[-1/(1-x2/L2)], or 0
    function Mollify(r, rScl) result(M)
        real(rp), intent(in) :: r, rScl
        real(rp) :: M

        real(rp) :: eArg, a, unity
        unity = 1.0_rp

        a = abs(r/rScl)
        if ((a + TINY) <= unity) then
            eArg = unity - a**2.0
            M = exp(unity)*exp(-unity/eArg)
        else
            M = 0.0
        endif

    end function Mollify

    !Smooth rampdown function
    !M = 1, r<=rC
    !M = 0, r>=rC+lC
    function RampDown(r, rC, lC) result(M)
        real(rp), intent(in) :: r, rC, lC
        real(rp) :: M

        real(rp) :: rScl, eArg
        if (r <= rC) then
            M = 1.0
        else if (r >= rC + lC) then
            M = 0.0
        else
            rScl = (r - rC)/lC
            eArg = -1.0/(1.0-rScl*rScl)
            M = exp(eArg + 1.0)
        endif
    end function RampDown

    function LinRampDown(r, rC, lC) result(M)
        real(rp), intent(in) :: r, rC, lC
        real(rp) :: M

        real(rp) :: rScl

        M = min(1.0, max((rC + lC - r)/(lC), 0.0))

    end function LinRampDown

    function LinRampUp(r, rC, lC) result(M)
        real(rp), intent(in) :: r, rC, lC
        real(rp) :: M

        real(rp) :: rScl

        rScl = (r - rC)/lC
        if (rScl <= 0.0) then
            M = 0.0
        else if (rScl >= 1.0) then
            M = 1.0
        else
            M = rScl
        endif

    end function LinRampUp

    function CubicRampDown(r, rC, lC) result(M)
        real(rp), intent(in) :: r, rC, lC
        real(rp) :: M

        real(rp) :: rScl

        rScl = (r - rC)/lC
        if (rScl <= 0.0) then
            M = 1.0
        else if (rScl >= 1.0) then
            M = 0.0
        else
            M = 1.0-3*(rScl**2.0) + 2*(rScl**3.0)
        endif
    end function CubicRampDown

    function PenticRampDown(r, rC, lC) result(M)
        real(rp), intent(in) :: r, rC, lC
        real(rp) :: M, Mc

        real(rp) :: rScl

        rScl = (r - rC)/lC
        if (rScl <= 0.0) then
            M = 1.0
        else if (rScl >= 1.0) then
            M = 0.0
        else
            Mc = 6*(rScl**5.0) - 15*(rScl**4.0) + 10*(rScl**3.0)
            M = 1.0-Mc
        endif

    end function PenticRampDown

    function RampUp(r, rC, lC) result(M)
        real(rp), intent(in) :: r, rC, lC
        real(rp) :: M

        real(rp) :: rScl, eArg
        rScl = (r - rC)/lC
        if (rScl <= 0.0) then
            M = 0.0
        else if (rScl >= 1.0) then
            M = 1.0
        else
            eArg = -1.0/(rScl*rScl)
            M = exp(1.0)*exp(eArg)
        endif

    end function RampUp

    !Generic min mod function
    function minmod(x, y) result(d)
        real(rp), intent(in) :: x, y
        real(rp) :: d
        real(rp), parameter :: ONE = 1.0
        d = 0.5*(sign(ONE, x) + sign(ONE, y))*min(abs(x), abs(y))
    end function minmod
!--------------------------------------
!Vector/tensor stuff

    !Calculate vector dot tensor, = V \cdot T
    function VdT(V, T)
        real(rp), intent(in) :: V(NDIM), T(NDIM, NDIM)
        real(rp), dimension(NDIM) :: VdT

        VdT(XDIR) = dot_product(V, T(:, XDIR))
        VdT(YDIR) = dot_product(V, T(:, YDIR))
        VdT(ZDIR) = dot_product(V, T(:, ZDIR))
    end function VdT

    !Turn vector into matrix xA so that, AxB = xA*B (matrix-mult)
    function xMat(A) result(xA)
        real(rp), dimension(NDIM), intent(in) :: A
        real(rp), dimension(NDIM, NDIM) :: xA

        xA(XDIR, :) = [0.0_rp, -A(ZDIR), A(YDIR)]
        xA(YDIR, :) = [A(ZDIR), 0.0_rp, -A(XDIR)]
        xA(ZDIR, :) = [-A(YDIR), A(XDIR), 0.0_rp]

    end function xMat

    !Calculate dyadic tensor from A,B
    function Dyad(A, B) result(AB)
        real(rp), dimension(NDIM), intent(in) :: A, B
        real(rp), dimension(NDIM, NDIM) :: AB

        AB(XDIR, :) = A(XDIR)*B
        AB(YDIR, :) = A(YDIR)*B
        AB(ZDIR, :) = A(ZDIR)*B

    end function Dyad

    !Calculate ab \cross^{\cdot} M
    !a,b vectors, M tensor
    ! ab dCross M = a \cross (b \cdot M)
    function dCross(a, b, M)
        real(rp), intent(in), dimension(NDIM) :: a, b
        real(rp), intent(in) :: M(NDIM, NDIM)
        real(rp) :: dCross(NDIM)

        dCross = cross(a, VdT(b, M))

    end function dCross

    !A:B, tensor contraction
    function tColon(A, B)
        real(rp), intent(in), dimension(NDIM, NDIM) :: A, B
        real(rp) :: tColon

        tColon = sum(A*B)

    end function tColon

!--------------------------------------
!Geometry stuff
    !atan2 w/ trap for pi/-pi jump
    elemental function katan2(a,b) result(c)
        real(rp), intent(in) :: a,b
        real(rp) :: c

        !Ensure arctan is between 0,2pi
        c = modulo(atan2(a,b) + 2*PI,2*PI)
    end function katan2

    function cross(a, b)
        !DEC$ ATTRIBUTES NOINLINE :: cross
        real(rp), dimension(NDIM), intent(in) :: a, b
        real(rp), dimension(NDIM) :: cross

        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)

    end function cross

    !Clamp a in [aL,aH]
    subroutine ClampValue(a,aL,aH)
        real(rp), intent(inout) :: a
        real(rp), intent(in) :: aL,aH
        if (a<aL) a = aL
        if (a>aH) a = aH
    end subroutine ClampValue
    
    !Circular mean, in: rad / out: rad
    function CircMean(alpha) result(alphabar)
        real(rp), intent(in) :: alpha(:)
        real(rp) :: alphabar
        integer :: N
        real(rp) :: X,Y
        
        N = size(alpha)
        Y = sum(sin(alpha))/N
        X = sum(cos(alpha))/N

        alphabar = modulo( atan2(Y,X),2*PI )
    end function CircMean
    
    !Circular mean, in: deg / out: deg
    function CircMeanDeg(alphadeg) result(alphabar)
        real(rp), intent(in) :: alphadeg(:)
        real(rp) :: alphabar
        
        alphabar = CircMean(deg2rad*alphadeg)*rad2deg

    end function CircMeanDeg

    function ArithMean(alpha) result(alphabar)
        real(rp), intent(in) :: alpha(:)
        real(rp) :: alphabar
        integer :: N
        real(rp) :: X,Y
        
        N = size(alpha)
        alphabar = sum(alpha)/N

    end function ArithMean

    !Weighted means
    function wArithMean(alpha,wgt) result(alphabar)
        real(rp), intent(in) :: alpha(:),wgt(:)
        real(rp) :: alphabar
        integer :: N
        real(rp) :: X,Y
        
        if (size(alpha) /= size(wgt)) then
            write(*,*) "wArithMean size error ..."
            stop
        endif
        alphabar = dot_product(alpha,wgt)/sum(wgt)
    end function wArithMean

    function normVec(a)
        real(rp), dimension(NDIM) :: a, normVec
        real(rp) :: normA

        normA = sqrt(dot_product(a, a))

        if (normA <= TINY) then
            normVec(:) = 0.0
        else
            normVec = a/normA
        endif
    end function normVec

    !Calculates angle between two vectors
    function angVec(a, b) result(alpha)
        real(rp), dimension(3), intent(in) :: a, b
        real(rp) :: alpha, nProd

        nProd = max(norm2(a), TINY)*max(norm2(b), TINY)
        alpha = acos(dot_product(a, b)/(nProd))

    end function

    !Simple function for print statements
    function ang2deg(a)
        real(rp), intent(in) :: a
        real(rp) :: ang2deg
        ang2deg = a*180.0/PI
    end function
    
    !Return // (to nhat) component of vector A
    !Assuming nhat is already unit vector
    function Vec2Para(A, nhat) result(Apara)
        real(rp), dimension(NDIM), intent(in) :: A, nhat
        real(rp), dimension(NDIM) :: Apara

        Apara = dot_product(A, nhat)*nhat
    end function Vec2Para

    function Vec2Perp(A, nhat) result(Aperp)
        real(rp), dimension(NDIM), intent(in) :: A, nhat
        real(rp), dimension(NDIM) :: Aperp

        Aperp = A - Vec2Para(A, nhat)
    end function Vec2Perp

    !Check if 2D point xp is in 2D cell bounded by xCs
    !Assuing CCW orientation, 0->1->2->3->0
    function inCell2D(xp, xCs) result(inCell)
        real(rp), intent(in) :: xp(2), xCs(4, 2)
        logical :: inCell
        logical :: inT1, inT2

        inT1 = inTri(xp, xCs(1, :), xCs(2, :), xCs(3, :))
        inT2 = inTri(xp, xCs(3, :), xCs(4, :), xCs(1, :))
        inCell = inT1 .or. inT2

    end function inCell2D

    !Check if 2D point xp is in triangle formed by p0,p1,p2
    function inTri(xp, p0, p1, p2)
        real(rp), dimension(2), intent(in) :: xp, p0, p1, p2
        logical :: inTri

        real(rp) :: p0x, p0y, p1x, p1y, p2x, p2y, x, y
        real(rp) :: A, s, t

        !Pull variables
        p0x = p0(1)
        p0y = p0(2)
        p1x = p1(1)
        p1y = p1(2)
        p2x = p2(1)
        p2y = p2(2)
        x = xp(1)
        y = xp(2)

        !Do test
        A = 0.5*(-p1y*p2x+p0y*(p2x-p1x) + p0x*(p1y-p2y) + p1x*p2y)
        s = (1/(2*A))*(p0y*p2x-p0x*p2y+(p2y-p0y)*x + (p0x-p2x)*y)
        t = (1/(2*A))*(p0x*p1y-p0y*p1x+(p0y-p1y)*x + (p1x-p0x)*y)

        InTri = (s >= 0) .and. (t >= 0) .and. ((1 - s - t) >= 0)
        !write(*,*) 's,t,1-s-t = ',s,t,1-s-t
    end function inTri

    !Convert XYZ to spherical (RTP) vector at given point
    function xyz2rtp(xyz,Axyz) result(Artp)
        real(rp), dimension(NDIM), intent(in) :: xyz,Axyz
        real(rp), dimension(NDIM) :: Artp

        real(rp) :: Ax,Ay,Az,Ar,At,Ap,rad,theta,phi
        Ax = Axyz(XDIR)
        Ay = Axyz(YDIR)
        Az = Axyz(ZDIR)

        rad = norm2(xyz)
        if (rad < TINY) then
            Artp = 0.0
            return
        endif

        theta = acos(xyz(ZDIR)/rad)
        phi = atan2(xyz(YDIR),xyz(XDIR))

        Ar =  Ax*sin(theta)*cos(phi) + Ay*sin(theta)*sin(phi) + Az*cos(theta)
        At =  Ax*cos(theta)*cos(phi) + Ay*cos(theta)*sin(phi) - Az*sin(theta)
        Ap = -Ax           *sin(phi) + Ay           *cos(phi) 
        Artp = [Ar,At,Ap]
    end function xyz2rtp
    
    !Convert RTP spherical to XYZ vector at given point
    function rtp2xyz(xyz,Artp) result(Axyz)
        real(rp), dimension(NDIM), intent(in) :: xyz,Artp
        real(rp), dimension(NDIM) :: Axyz

        real(rp) :: Ax,Ay,Az,Ar,At,Ap,rad,theta,phi
        Ar = Artp(1)
        At = Artp(2)
        Ap = Artp(3)

        rad = norm2(xyz)
        if (rad < TINY) then
            Axyz = 0.0
            return
        endif
        theta = acos(xyz(ZDIR)/rad)
        phi = atan2(xyz(YDIR),xyz(XDIR))

        Ax = Ar*sin(theta)*cos(phi) + At*cos(theta)*cos(phi) - Ap*sin(phi)
        Ay = Ar*sin(theta)*sin(phi) + At*cos(theta)*sin(phi) + Ap*cos(phi)
        Az = Ar*cos(theta)          - At*sin(theta)          !No phi contribution

        Axyz = [Ax,Ay,Az]
    end function rtp2xyz
!--------------------------------------
!Linear algebra stuff
    !Calculate determinant of 3X3 matrix
    function DetJ(J) result(dJ)
        real(rp), intent(in) :: J(NDIM, NDIM)
        real(rp) :: dJ

        dJ =   J(1, 1)*(J(2, 2)*J(3, 3) - J(2, 3)*J(3, 2)) &
             - J(1, 2)*(J(2, 1)*J(3, 3) - J(2, 3)*J(3, 1)) &
             + J(1, 3)*(J(2, 1)*J(3, 2) - J(2, 2)*J(3, 1))

    end function DetJ

    !Invert 3x3 matrix
    !Code adapted from http://web.hku.hk/~gdli/UsefulFiles/matrix/m33inv_f90.txt

    subroutine matinv3(A, B)
        real(rp), intent(in)  :: A(3, 3)
        real(rp), intent(out) :: B(3, 3)

        real(rp) :: dA, cofac(3, 3)

        dA = DetJ(A)

        cofac(1, 1) = +(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))
        cofac(1, 2) = -(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))
        cofac(1, 3) = +(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
        cofac(2, 1) = -(A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))
        cofac(2, 2) = +(A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))
        cofac(2, 3) = -(A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))
        cofac(3, 1) = +(A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))
        cofac(3, 2) = -(A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))
        cofac(3, 3) = +(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))
        if (dA < TINY) then
            B = 0
            write (*, *) 'Warning! Inverting singular matrix, A = ', A
        else
            B = transpose(cofac)/dA
        endif
    end subroutine matinv3

    subroutine matinv4(A, B)
        !! Performs a direct calculation of the inverse of a 4×4 matrix.
        !! grabbed by vgm from the internet
        !! FIXME: check coeffs
        real(rp), intent(in) :: A(4, 4) !! Matrix
        real(rp), intent(out):: B(4, 4) !! Inverse matrix
        real(rp)             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))) &
                 - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))  &
                 + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))  &
                 - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

        ! Calculate the inverse of the matrix
        B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
        B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
        B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
        B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
        B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
        B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
        B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
        B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
        B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
        B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
        B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
        B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    end subroutine matinv4

!--------------------------------------
!Simple smoothing stuff
    !Use 5x5 smoothing window, isGO is optional logical for good cells
    function SmoothOperator55(Q,isGO) result(Qs)
        real(rp), dimension(5,5), intent(in)  :: Q
        logical , dimension(5,5), intent(in), optional  :: isGO
        real(rp) :: Qs
        logical, dimension(5,5) :: isG
        real(rp) :: qAvg,qScl

        if (present(isGO)) then
            isG = isGO
        else
            isG = .true.
        endif

        if (any(isG)) then
            qAvg = sum(SmoothOp5x5*Q,mask=isG)
            qScl = sum(SmoothOp5x5  ,mask=isG)
            Qs = qAvg/qScl
        else
            Qs = 0.0
        endif
        
    end function SmoothOperator55

    !Use 3x3 smoothing window, isGO is optional logical for good cells
    function SmoothOperator33(Q,isGO) result(Qs)
        real(rp), dimension(3,3), intent(in)  :: Q
        logical , dimension(3,3), intent(in), optional  :: isGO
        real(rp) :: Qs
        logical, dimension(3,3) :: isG
        real(rp) :: qAvg,qScl

        if (present(isGO)) then
            isG = isGO
        else
            isG = .true.
        endif

        if (any(isG)) then
            qAvg = sum(SmoothOpTSC*Q,mask=isG)
            qScl = sum(SmoothOpTSC  ,mask=isG)
            Qs = qAvg/qScl
        else
            Qs = 0.0
        endif
        
    end function SmoothOperator33

    !Soft floor (smooth max)
    elemental function SoftFloor(X,XMin) result(sX)
        real(rp), intent(in) :: X, XMin
        real(rp), parameter  :: f = 0.1_rp

        real(rp) :: sX

        sX = 0.5_rp * ( X + XMin + &
             sqrt( (X - XMin)**2.0_rp + (f*XMin)**2.0_rp ) )
    end function SoftFloor

    function isAscending(A) result(isAsc)
      real(rp), dimension(:), intent(in) :: A
      logical :: isAsc
      integer :: N

      N = size(A)

      if ( all(A(2:N)-A(1:N-1)>0) ) then
         isAsc = .True. 
      else
         isAsc = .False. 
      end if
    end function isAscending

    !> Simple corner-to-center function where everyone can get to it easily
    function toCenter2D(v) result (vcc)
        !! 2D cell center of a value
        real(rp), dimension(2,2), intent(in) :: v

        real(rp) :: vcc

        vcc = 0.25*(v(1,1) + v(2,1) + v(1,2) + v(2,2))
    end function toCenter2D
end module math
