module raijuEleLossHelper
    !! Helper functions for common electron loss calculations

    use kdefs
    use math

    use raijudefs
    use raijutypes

    

    implicit none

    contains

    function CalcTau_StrongScattering(Model, Grid, State, i, j, k) result(tau)
        !! Calculates strong scattering rate, according to Schulz 1998
        !! tau ~ [2*FTV*Bfp/(1-eta)](gamma*m0/p)
        !! FTV = flux tube volume, Bfp = B-field at foot point, eta - back-scattering rate
        !! eta is backscatter rate at alitude h, here eta=2/3.
        !! gamma = m/m0 is relativisitc factor, p is particle momentum.
        !!       = mc2/m0c2 = (m0c2+K)/m0c2 = 1+K/mec2 ! mec2=0.511 is me*c^2 in MeV
        !! m = m0/sqrt(1-v^2/c^2)
        !! V = c*1/sqrt(1-1/gammar2)
        type(raijuModel_T ), intent(in) :: Model
        type(raijuGrid_T  ), intent(in) :: Grid
        type(raijuState_T ), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau

        real(rp) :: KE, gammar, V, vm
        real(rp) :: eta_scatter = 2./3.

        tau = HUGE

        if (State%bvol_cc(i,j) < TINY) then
            write(*,*) "ERROR in CalcTau_StrongScattering: bVol for requested point is below TINY!"
            write(*,*)State%bvol_cc(i,j)
            call State%myRaiApp%WriteBlackBox()
            write(*,*)"Goodbye"
            stop
        else
            vm = State%bvol_cc(i,j)**(-2./3.)
        endif

        KE = abs(Grid%alamc(k))*vm * 1.0D-3  ! Energy [keV]
        gammar = 1.0 + (KE*1.0D-3)/mec2  ! Gamma with 1 + MeV/MeV
        V = (vc_cgs*1e-2)*sqrt(1.0 - 1.0/gammar**2)/Model%planet%rp_m  ! [Rp/s]
        tau = 2.0*State%bvol_cc(i,j)*Grid%Bmag(i,j)/(1.0 - eta_scatter) / V*gammar  ! [Rp/nT * nT / (Rp/s) = s]

        call ClampValue(tau, TINY, HUGE)

    end function CalcTau_StrongScattering


    function CalcTau_WeakScattering(Model, Grid, State, i, j, k) result(tau)
        !! Calculates weak scattering rate, according to Chen et al. 2005
        !! This is the version only dependent on E and L, not MLT
        !! I think the MLT dependence tries to account for the plasmasphere, but we have hiss so we will use that directly
        !! eq. 10: l_o = min(0.08(E, MEV)^(-1.32), 0.4*10^(2L - 6 + 0.4*log_2(E, MEV))) [1/days]
        type(raijuModel_T ), intent(in) :: Model
        type(raijuGrid_T  ), intent(in) :: Grid
        type(raijuState_T ), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau

        real(rp) :: KE, L, vm
        real(rp) :: A, B, l_0, tau_SS

        if (State%bvol_cc(i,j) < TINY) then
            write(*,*) "ERROR in CalcTau_WeakScattering: bVol for requested point is below TINY!"
            write(*,*)State%bvol_cc(i,j)
            call State%myRaiApp%WriteBlackBox()
            write(*,*)"Goodbye"
            stop
        else
            vm = State%bvol_cc(i,j)**(-2./3.)
        endif

        KE = abs(Grid%alamc(k))*vm * 1.0D-6  ! Energy [MeV]

        L = sqrt(State%xyzMincc(i,j,XDIR)**2 + State%xyzMincc(i,j,YDIR)**2)  ! Eq projection of bmin point, kinda L shell [Rp]
        A = 0.08_rp*KE**(-1.32_rp)
        B = 0.4*10**(2.0_rp*L - 6.0_rp + 0.4_rp*log10(KE)/log10(2.0_rp))
        l_0 = min(A, B)  ! [1/days]
        tau_SS = CalcTau_StrongScattering(Model, Grid, State, i, j, k)/86400.0_rp  ! [s --> days]
        
        ! Smooth transition between weak and strong scattering
        tau = (1.0_rp + l_0*tau_SS) / l_0  ! [days]
        tau = tau*86400.0_rp  ! [s]
        call ClampValue(tau, TINY, HUGE)  ! Just in case something weird happened

    end function CalcTau_WeakScattering


    function CalcTau_Hiss(MLTin, Lin, Ein, Kpin) result(tau)
        ! Empirical lifetime against plasmaspheric hiss pitch angle diffusion, based on Orlova et al. 2015JA021878.
        ! Improvements relative to 2014GL060100: 1. Hiss wave intensity distribution model is based on new data 
        ! (O14 was based on single-component E field in CRRES data. O16 used Spasojevic+2015 model based on EMFISIS B data on VAP); 
        ! 2. Wave spectrum is assumed differently (O14 assume Gaussian spectrum based on CRRES data).
        ! Electron lifetime tau(L,E,MLT,Kp) = tau_av(L,E)/g(MLT)/h(Kp), 
        ! where 1.5<L<5.5, E=log10(Ek[MeV]) for 1 KeV < Ek < 10 MeV.
        ! log10(tau_av(L,E)) = a1+a2*L+a3*E+...+a20*E^3, when E >= f(L).
        ! f(L) = 0.1328*L^2-2.1463*L+3.7857.
        ! g(MLT) = 10^g0(MLT)/G0
        ! h(Kp) = 10^h0(Kp)/H0
        !   G0 = int_0^24(10^g0(MLT))dMLT / 24 = 782.3.
        !   g0(MLT) = b2*MLT^2 + b1*MLT + b0
        !   H0 = 1315.
        !   h0(Kp) = c2*Kp^2 + c1*Kp + c0
        real(rp), intent(in) :: MLTin, Lin
            !! L in Re
        real(rp), intent(in) :: Ein, Kpin
            !! E in MeV

        real(rp) :: tau_av, tau, rateK
            !! tau in s, rate in 1/s
        real(rp) :: MLT, L, E, Kp
            !! Params we clamp if needed
        real(rp) :: L2, L3, L4, fL, E2, E3, E4, E5, LE
            !! Polynomials
        real(rp), dimension(20) :: le_pol
            !! Container for polynomials
        real(rp) :: b0, b1, b2, G0, g0_MLT, g_MLT, c0, c1, c2, H0, h0_Kp, h_Kp
            !! Coefficients for g_MLT and h_Kp
        real(rp), dimension(20), parameter :: a1_20 = [77.323, -92.641, -55.754, 44.497, 48.981, 8.9067, -10.704, &
                                                         -15.711, -3.3326, 1.5189, 1.294, 2.2546, 0.31889, -0.85916, & 
                                                         -0.22182, 0.034318, 0.097248, -0.12192, -0.062765, 0.0063218]
            !! Coefficients for polynomials

        rateK = 0.0

        MLT = MLTin
        Kp = Kpin
        L = Lin
        call ClampValue(L,1.5_rp,5.5_rp)
        L2 = L*L
        fL = 0.1328*L2 - 2.1463*L + 3.7857
        E = log10(Ein)
        call ClampValue(E,max(-3.0_rp,fL),1.0_rp) ! 1 keV to 1 Mev

        E2 = E*E
        E3 = E2*E
        E4 = E3*E
        E5 = E4*E
        L3 = L2*L
        L4 = L3*L
        LE = L*E
        le_pol = (/1.D0,L,E,L2,LE,E2,L3,L2*E,L*E2,E3,L4,L3*E,L2*E2,L*E3,E4,L*E4,L2*E3,L4*E,L2*L3,E5/)

        b0 = 2.080
        b1 = 0.1773
        b2 = -0.007338
        G0 = 782.3
        g0_MLT = b2*MLT*MLT + b1*MLT + b0
        g_MLT = 10**g0_MLT/G0
        c0 = 2.598
        c1 = 0.2321
        c2 = -0.01414
        H0 = 1315.0
        Kp = min(Kp,5.0_rp) ! 0<Kp<5, 1.2% data in Kp bin of 4.3-7.6
        h0_Kp = c2*Kp*Kp + c1*Kp + c0
        h_Kp = 10**h0_Kp/H0

        tau_av = 10.0**(dot_product(a1_20,le_pol))*86400.D0 ! seconds
        tau = tau_av/g_MLT/h_Kp
        call ClampValue(tau, TINY, HUGE)
        !rateK = 1.0/tau

        !write(*,*) Ein, L, fL, E, tau, rateK, '--'
        !stop

    end function CalcTau_Hiss    

end module raijuEleLossHelper