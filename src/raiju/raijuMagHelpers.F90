module raijuMagHelpers
    !! Helper routines to calculate certain misc magnetosphere-related quantities

    use raijudefs
    use raijutypes

    implicit none

    contains

    function DriftPeriod(D, P, L, Bmin, Rc, Rp_m) result(tauD)
        real(rp), intent(in) :: D
            !! number density [#/cc]
        real(rp), intent(in) :: P
            !! [nPa]
        real(rp), intent(in) :: L
            !! [Rp]
        real(rp), intent(in) :: Bmin
            !! [nT]
        real(rp), intent(in) :: Rc
            !! [Rp] Curvature radius
        real(rp), intent(in) :: Rp_m
            !! [m] planetary radius 

        real(rp) :: tauD
            !! [s]
        real(rp) :: kT, Vd
            !! [keV], [m/s]

        kT = DP2kT(D, P)
        ! V_drift = 156*K/Rc/B [km/s], K[keV], Rc[Rp], B[nT]
        Vd = 1.e3*156.0*kT/Rc/Bmin  ! [m/s]
        tauD = 2.0*PI*L*Rp_m/Vd
    end function DriftPeriod


    function DipoleDriftPeriod(D, P, L) result(tauD)
        real(rp), intent(in) :: D
            !! number density [#/cc]
        real(rp), intent(in) :: P
            !! [nPa]
        real(rp), intent(in) :: L
            !! [Rp]

        real(rp) :: tauD
            !! [s]
        real(rp) :: kT

        kT = DP2kT(D, P)
        ! Using Td = 700/K/L hours
        tauD = 3600.0*700/kT/L  ! [s]
    end function DipoleDriftPeriod

end module raijuMagHelpers