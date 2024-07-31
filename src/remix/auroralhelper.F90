!Module to hold routines to help calculate auroral conductance

module auroralhelper
    use kdefs
    use mixdefs
    use mixtypes
    
    implicit none

    real(rp), private, parameter :: ET_E0 = 35.0 !E field magnitude [mV/m] for ET amplification
    real(rp), private, parameter :: KMax = 100.0 ![keV], max energy to allow for Kap23 conductance formula
    real(rp), private, parameter :: KMin = 1.0   ![keV], min energy to allow for Kap23 conductance formula

    contains

!Functions for auroral conductance from precipitation values
    !Returns Robinson's SigP from eavg [keV] and eflux [ergs/cm^2]
    elemental function SigmaP_Robinson(eavg,eflux) result(SigP)
        real(rp), intent(in) :: eavg,eflux
        real(rp) :: SigP
        SigP = 40.0*eavg*sqrt(eflux)/(16.0 + eavg**2.0)
    end function SigmaP_Robinson

    !Returns Robinson's SigH from eavg [keV] and eflux [ergs/cm^2]
    !NOTE: Extra correction from Fedder
    elemental function SigmaH_Robinson(eavg,eflux) result(SigH)
        real(rp), intent(in) :: eavg,eflux
        real(rp) :: SigH
        real(rp) :: SigP

        SigP = SigmaP_Robinson(eavg,eflux)
        ! Option to use Fedder correction to curb values for high eavg: SigH = 0.45*SigP*(eavg**0.85)/(1.0 + 0.0025*eavg**2.0)
        ! By default, use Kaeppler+ 15 correction to SigH/SigP from Robinson 
        SigH = 0.57*SigP*(eavg**0.53)
    end function SigmaH_Robinson


!New routines based on Kaeppler++ 2023, 10.1029/2023JA031764
    elemental function SigmaP_Kaep23(eavg,eflux,dE) result(SigP)
        real(rp), intent(in) :: eavg,eflux,dE
        real(rp) :: SigP

        real(rp), dimension(4), parameter :: A_M = [21.27, 1.32, 7.51, 2.09] !Mono
        real(rp), dimension(4), parameter :: A_D = [17.38, 1.67, 6.77, 2.37] !Diffuse

        !Choose mono vs. diffuse based on dE
        if (dE > mixeTINY) then
            SigP = SigXK(eavg,eflux,A_M)
        else
            SigP = SigXK(eavg,eflux,A_D)
        endif
    end function SigmaP_Kaep23

    elemental function SigmaH_Kaep23(eavg,eflux,dE) result(SigH)
        real(rp), intent(in) :: eavg,eflux,dE
        real(rp) :: SigH

        real(rp), dimension(4), parameter :: A_M = [26.65, 1.47,10.24, 1.74] !Mono
        real(rp), dimension(4), parameter :: A_D = [23.72, 1.24, 6.55, 1.45] !Diffuse
        !Choose mono vs. diffuse based on dE
        if (dE > mixeTINY) then
            SigH = SigXK(eavg,eflux,A_M)
        else
            SigH = SigXK(eavg,eflux,A_D)
        endif

    end function SigmaH_Kaep23

    !Calculate Kaeppler++ 23 format fit given coefficients A
    pure function SigXK(eavgin,eflux,A) result(Sig_X)
        real(rp), intent(in) :: eavgin,eflux
        real(rp), intent(in) :: A(4)
        real(rp) :: Sig_X
        real(rp) :: eavg

        eavg = eavgin
        !Clamp to validity range (wonky to use other routines b/c pure)
        eavg = min(eavg,KMax)
        eavg = max(eavg,KMin)

        Sig_X = sqrt(eflux)*A(1)*(eavg**A(2))/(A(3) + eavg**A(4))
    end function SigXK

!---
!Functions for electrojet turbulence correction (useful for keeping CPCP down)

    !Calculate magnitude of E field (unsigned b/c mobius strip)
    !pot is kV, MagE is mV/m
    subroutine CalcMagE(Gr,pot,MagE)
        type(mixGrid_T), intent(in) :: Gr
        real(rp), dimension(Gr%Np,Gr%Nt), intent(in)  :: pot
        real(rp), dimension(Gr%Np,Gr%Nt), intent(out) :: MagE

        integer :: i,j
        real(rp) :: R0,Et,Ep,dth,dph,th,EAvg


        MagE = 0.0

        !Start by calculating kV/Re
        R0 = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,dth,dph,th,Et,Ep)
        do i=1,Gr%Np
            do j=1+1,Gr%Nt-1
                !NOTE: only care about magnitude
                dth = Gr%dt(i,j)
                dph = Gr%dp(i,j)
                th  = Gr%t (i,j)

                !Get theta deriv
                if (j == 1) then
                    Et = ( pot(i,j+1)-pot(i,j) )/(R0*dth)
                else if (j == Gr%Nt) then
                    Et = ( pot(i,j-1) - pot(i,j) )/(R0*dth)
                else
                    Et = ( pot(i,j+1) - pot(i,j-1) )/(2*R0*dth)
                endif

                !Get phi deriv
                if (i == 1) then
                    Ep = ( pot(i+1,j) - pot(i,j) )/(R0*sin(th)*dph)
                else if (i == Gr%Np) then
                    Ep = ( pot(i,j) - pot(i-1,j) )/(R0*sin(th)*dph)
                else
                    Ep = ( pot(i+1,j) - pot(i-1,j) )/(2*R0*sin(th)*dph)
                endif

                MagE(i,j) = sqrt( Et**2.0 + Ep**2.0 )
            enddo
        enddo

        EAvg = sum(MagE(:,2))/Gr%Np
        MagE(:,1) = EAvg

        EAvg = sum(MagE(:,Gr%Nt-1))/Gr%Np
        MagE(:,Gr%Nt) = EAvg
        
        !Have MagE in kV/Re, convert to mV/m
        MagE = (1.0e+3)*(1.0e+3)*MagE/REarth

    end subroutine CalcMagE

    elemental function SigET_P(E)
        real(rp), intent(in) :: E
        real(rp) :: SigET_P
        if (E>ET_E0) then
            SigET_P = 1.0 + 0.01*(E-ET_E0) + 1.3*(1.0e-5)*(E-ET_E0)**2.0
            SigET_P = max(SigET_P,1.0)
        else
            SigET_P = 1.0
        endif

    end function SigET_P

    elemental function SigET_H(E)
        real(rp), intent(in) :: E
        real(rp) :: SigET_H

        if (E>ET_E0) then
            SigET_H = 1.0 + 0.01172*(E-ET_E0) - 1.207*(1.0e-5)*(E-ET_E0)**2.0
            SigET_H = max(SigET_H,1.0)
        else
            SigET_H = 1.0
        endif

    end function SigET_H

end module auroralhelper
