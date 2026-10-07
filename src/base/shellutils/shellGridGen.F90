!> Helper module for generating theta and phi arrays to use within a ShellGrid
module shellGridGen
    use math
    use kdefs

    implicit none

    contains

    subroutine genThetaPhi_uniform(Nt, Np, thetaL, thetaU, theta, phi)
        integer, intent(in) :: Nt
            !! Number of active cells in theta direction
        integer, intent(in) :: Np
            !! Number of active cells in phi direction
        real(rp), intent(in) :: thetaL
            !! [deg] Lower theta bound (0 = pole, 90 = equator)
        real(rp), intent(in) :: thetaU
            !! [deg] Upper theta bound (0 = pole, 90 = equator)
        real(rp), dimension(Nt+1), intent(inout) :: theta
            !! Populated active theta grid we return
        real(rp), dimension(Np+1), intent(inout) :: phi
            !! Populated active phi grid we return

        real(rp) :: thetaL_rad, thetaU_rad
        real(rp) :: dTheta, dPhi
        integer :: i

        ! Get theta bounds in radians
        thetaL_rad = thetaL*deg2rad
        thetaU_rad = thetaU*deg2rad

        dTheta = (thetaU_rad-thetaL_rad)/Nt
        dPhi = 2.0*PI/Np

        do i=1,Nt+1
            theta(i) = thetaL_rad + (i-1)*dTheta
        enddo

        do i=1,Np+1
            phi(i) = (i-1)*dPhi
        enddo
        ! Catch for slight overshoot
        if ((phi(Np+1) > 2*PI) .and. (phi(Np+1) - 2*PI < TINY) ) then
            phi(Np+1) = 2.0*PI
        endif

    end subroutine genThetaPhi_Uniform


    subroutine genThetaPhi_Shafee2008(Nt, Np, thetaL, thetaU, nPow, hWgt, xLow, xHigh, theta, phi)
        !! Modified version of Shafee+ 2008 (10.1086/593148)
        integer, intent(in) :: Nt
            !! Number of active cells in theta direction
        integer, intent(in) :: Np
            !! Number of active cells in phi direction
        real(rp), intent(in) :: thetaL
            !! [deg] Lower theta bound (0 = pole, 90 = equator)
        real(rp), intent(in) :: thetaU
            !! [deg] Upper theta bound (0 = pole, 90 = equator)
        integer :: nPow
            !! Power applied to non-linear scaling term
        real(rp) :: hWgt
            !! Weight between linear and non-linear term. 1=linear, 0=non-linear
        real(rp) :: xLow, xHigh
            !! Bounds for generating x values between 0 and 1
            !! FIXME: replace with theta_center and x_scale, and calculate our x_low and x_high from that
        real(rp), dimension(Nt+1), intent(inout) :: theta
        real(rp), dimension(Np+1), intent(inout) :: phi

        
        real(rp) :: a1, a2
            !! Coefficients to map from (x_low, x_high) to (thetaL, thetaU), using the scaling function
            !! If x_low=0 and x_high=1, a1 = thetaL and a2 = thetaU
            !! But this makes the resolution focus symmetric between thetaLa nd thetaU which is not necessarily what we want
        real(rp) :: thetaL_rad, thetaU_rad
        real(rp) :: x, dX, dPhi
        integer :: i

        ! Turn degrees into radians
        thetaL_rad = thetaL*deg2rad
        thetaU_rad = thetaU*deg2rad

        if (xLow < 0 ) then
            write(*,*) "ERROR in genThetaPhi_Shafee2008:"
            write(*,*) "xLow must be greater than zero"
            stop
        endif
        if (xHigh > 1) then
            write(*,*) "ERROR in genThetaPhi_Shafee2008:"
            write(*,*) "xHigh must be less than one"
            stop
        endif
        if (xLow > xHigh) then
            write(*,*) "ERROR in genThetaPhi_Shafee2008:"
            write(*,*) "xLow must be less than xHigh"
            stop
        endif

        call calcCoeffs(thetaL_rad, thetaU_rad, xLow, xHigh, hWgt, nPow, a1, a2)
        dX = (xHigh - xLow) / Nt

        dPhi = 2.0*PI/Np


        do i=1,Nt+1
            x = xLow + (i-1)*dX
            theta(i) = a1 + (a2-a1)*scaleFunc(x, hWgt, nPow)
        enddo

        ! Pin the prescribed bounds before callers mirror a hemispheric grid;
        ! cancellation in the affine map can otherwise put them just outside.
        theta(1)    = thetaL_rad
        theta(Nt+1) = thetaU_rad

        do i=1,Np+1
            phi(i) = (i-1)*dPhi
        enddo
        ! Catch for slight overshoot
        if ((phi(Np+1) > 2*PI) .and. (phi(Np+1) - 2*PI < TINY) ) then
            phi(Np+1) = 2.0*PI
        endif

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

    end subroutine genThetaPhi_Shafee2008


end module shellGridGen
