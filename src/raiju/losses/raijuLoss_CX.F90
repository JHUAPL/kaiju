module raijuLoss_CX
    !! Charge exchange losses for RAIJU

    use raijudefs
    use raijutypes
    use raijuSpeciesHelper

    implicit none

    type, extends(baseRaijuLoss_T) :: raiLoss_CX_T
        logical :: reqsGood = .false.
        real(rp), dimension(2) :: kevBnd_H = [0.005_rp, 250.0_rp]
        real(rp), dimension(2) :: kevBnd_O = [0.025_rp, 600.0_rp]
        contains

        procedure :: doInit => CXLossInit
        !procedure :: doUpdate
        procedure :: isValidSpc => CXLossValidSpc
        procedure :: calcTau => CXLossCalcTau

    end type raiLoss_CX_T

    contains

    subroutine CXLossInit(this, Model, Grid, xmlInp)
        class(raiLoss_CX_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(XML_Input_T) , intent(in) :: xmlInp

        if (trim(toUpper(Model%planet%name)) .ne. "EARTH") then
            write(*,*)"WARNING in raijuLoss_CX: Not simulating Earth, idk what neutrals you want me to use"
            write(*,*)"No CX happening"
        else
            this%reqsGood = .true.
        endif
    end subroutine CXLossInit


    function CXLossValidSpc(this, spc) result(isValid)
        class(raiLoss_CX_T), intent(in) :: this
        type(raijuSpecies_T), intent(in) :: spc
        logical :: isValid

        isValid = .false.
        if (.not. this%reqsGood) return
        if (spc%flav == F_PSPH) return

        if ( (spc%spcType .eq. RAIJUHPLUS) &
        .or. (spc%spcType .eq. RAIJUOPLUS) ) then
            isValid = .true.
        endif

    end function CXLossValidSpc
 

    function CXLossCalcTau(this, Model, Grid, State, i, j, k) result(tau)
        class(raiLoss_CX_T), intent(in) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau

        real(rp) :: energy, rLoc, Ngeo, cxSig, M, V
        real(rp) :: tauEq, tauBA, vm

        associate(spc => Grid%spc(Grid%k2spc(k)))

        tau = HUGE
        if (.not. this%reqsGood) then
            return
        endif

        if (State%bvol_cc(i,j) < TINY) then
            write(*,*) "ERROR in CXLossCalcTau: bVol for requested point is below TINY!"
            write(*,*)State%bvol_cc(i,j)
            call State%myRaiApp%WriteBlackBox()
            write(*,*)"Goodbye"
            stop
        else
            vm = State%bvol_cc(i,j)**(-2./3.)
        endif

        ! Neutral
        rLoc = norm2(State%xyzMincc(i,j,:))  ! [Rp]
        Ngeo = OstgaardGeocorona(rLoc)  ! [#/cc]

        ! Ion
        energy = abs(Grid%alamc(k))*vm * 1D-3  ! keV
        if (spc%spcType .eq. RAIJUHPLUS) then
            call ClampValue(energy,this%kevBnd_H(1),this%kevBnd_H(2))
        elseif (spc%spcType .eq. RAIJUOPLUS) then
            call ClampValue(energy,this%kevBnd_O(1),this%kevBnd_O(2))
        endif
        cxSig = CXSigma(energy, spc%spcType)
        M = spc%amu * dalton  ! [kg]
        V = sqrt(2*(energy*kev2J)/M)*100.0  ! [m/s -> cm/s]
        
        ! Timescale
        tauEq = 1.0/(Ngeo*V*cxSig)  ! [s]
        ! Bounce-averaged approximation from Smith & Bewtra 1976
        ! NOTE: Technically, not derived from Ostgaard distribution. Fix me soon
        ! Also, angle is mirror lat, not pitch angle. So 45 deg is probably a bad estimate
        tauBA = tauEq * cos(45*PI/180.0)**3.5

        tau = max(tauBA, TINY)

        end associate

    end function CXLossCalcTau


!------
! Helpers
!------

    function OstgaardGeocorona(L) result(Ngeo)
        !! Geocoronal density afa L [#/cc], Taken from Ostgaard 2003 
        real(rp), intent(in) :: L
            !! L shell (radial dist?) [Rp]
        real(rp) :: Ngeo
            !! Density [#/cc]

        Ngeo = 10000.0*exp(-L/1.02) + 70.0*exp(-L/8.2)
    end function OstgaardGeocorona


    !Charge exchange cross-section for K [keV] and species ispc
    !Sig in cm2
    !Using Lindsay & Stebbings 2005
    function CXSigma(K,ispc) result(Sig)
        real(rp), intent(in) :: K
            !! Energy [keV]
        integer, intent(in) :: ispc
            !! RAIJU enum for species type
        real(rp) :: Sig
            !! Cross-section [cm^2]
        real(rp) :: Sig0, KSig,a1,a2,a3,B1,B2

        Sig0 = 1.0e-16

        select case(ispc)
            case(RAIJUHPLUS)
            !Charge exchange cross-section for H+/H
                !Cap for validity of CX cross-section
                KSig = K
                call ClampValue(KSig,0.005_rp,250.0_rp)

                a1 = 4.15
                a2 = 0.531
                a3 = 67.3

                B1 = (a1-a2*log(KSig))**2.0
                B2 = 1.0-exp(-a3/KSig) 
                Sig =  Sig0*B1*(B2**(4.5))
            case(RAIJUOPLUS)
            !Charge exchange cross-section for O+/H
                !Cap for validity of CX cross-section
                KSig = K
                call ClampValue(KSig,0.025_rp,600.0_rp)
                a1 = 3.13
                a2 = 0.170
                a3 = 87.5

                B1 = (a1-a2*log(KSig))**2.0
                B2 = 1.0-exp(-a3/KSig) 
                Sig =  Sig0*B1*(B2**(0.8))
            case default
                Sig = 0.0
        end select

    END FUNCTION CXSigma


!    subroutine setCXValidSpecies(this, Model, Grid)
!        class(raiLoss_CX_T), intent(inout) :: this
!        type(raijuModel_T), intent(in) :: Model
!        type(raijuGrid_T) , intent(in) :: Grid
!
!        integer :: s
!
!        allocate(lp%isValidSpc(Grid%nSpc))
!        this%isValidSpc = .false.
!
!        if (trim(toUpper(Model%planet%name)) .ne. "EARTH") then
!            return
!        endif
!
!        do s=1,Grid%nSpc
!            if ( (SpcType(Grid%spc(s)) .eq. RAIJUHPLUS) &
!            .or. (SpcType(Grid%spc(s)) .eq. RAIJUOPLUS) ) then
!                this%isValidSpc(s) = .true.
!            endif
!        enddo
!
!    end subroutine setCXValidSpecies

end module raijuLoss_CX