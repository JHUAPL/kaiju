module raijuLoss_CC
    !! Coulomb collision losses for RAIJU

    use raijudefs
    use raijutypes
    use raijuSpeciesHelper

    implicit none

    type, extends(baseRaijuLoss_T) :: raiLoss_CC_T
        logical :: reqsGood = .false.

        ! Idea: Have a local 3D copy of our taus and heat fluxes
        ! On step update, zero out locals and populate tau array. calcTau just returns already calculated value for this time step
        ! Additiona subroutine "AccumHeatFlux" will accum in local head flux array
        ! Then IO can reach in here for these outputs if desired

        contains

        procedure :: doInit => CCLossInit
        !procedure :: doUpdate
        procedure :: isValidSpc => CCLossValidSpc
        procedure :: calcTau => CCLossCalcTau

    end type raiLoss_CC_T

    contains


    subroutine CCLossInit(this, Model, Grid, xmlInp)
        class(raiLoss_CC_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(XML_Input_T) , intent(in) :: xmlInp

        ! Technically yes, but actually no
        !this%isPrecip = .true.

        ! Nothing to initialize, just make sure plasmasphere actually exists
        if (.not. Model%doPlasmasphere) then
            write(*,*)"WARNING in raijuLoss_CC.F90 init: need plasmasphere to do Coulomb collision losses."
            write(*,*)"i.e. Not cool enough for CC, so CC is off."
        else if (trim(toUpper(Model%planet%name)) .ne. "EARTH") then
            write(*,*)"WARNING in raijuLoss_CC: Not simulating EARTH, my CC model was fit to Earth"
            write(*,*)"If you want me to do CC, do something about it"
        else
            this%reqsGood = .true.
        endif
    end subroutine CCLossInit


    function CCLossValidSpc(this, spc) result(isValid)
        class(raiLoss_CC_T), intent(in) :: this
        type(raijuSpecies_T), intent(in) :: spc
        logical :: isValid

        isValid = .false.
        if (.not. this%reqsGood) return
        if (spc%flav == F_PSPH) return

        if ( (spc%spcType .eq. RAIJUHPLUS) ) then
            isValid = .true.
        endif
    end function


    function CCLossCalcTau(this, Model, Grid, State, i, j, k) result(tau)
        class(raiLoss_CC_T), intent(in) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau

        integer :: psphIdx
        real(rp) :: vm

        tau = HUGE
        if (.not. this%reqsGood) then
            return
        endif

        if (State%bvol_cc(i,j) < TINY) then
            write(*,*) "ERROR in CCLossCalcTau: bVol for requested point is below TINY!"
            write(*,*)State%bvol_cc(i,j)
            write(*,*)"i,j,k=",i,j,k
            call State%myRaiApp%WriteBlackBox()
            write(*,*)"Goodbye"
            stop
        else
            vm = State%bvol_cc(i,j)**(-2./3.)
        endif

        associate(spc => Grid%spc(Grid%k2spc(k)))
        psphIdx = spcIdx(Grid, F_PSPH)
        tau = CCTau(spc%spcType, Grid%alamc(k), &
                    vm, State%Den(psphIdx)%data(i,j))
        end associate

    end function CCLossCalcTau



    ! Simple Coulomb collision losses, using fit to Ebihara+ 98 Fig #5
    !  Note: assumes psph e- temp is 1 eV. Derived from Fok+ 1991
    function CCTau(spcType,alam,vm,Dpp) result(tau)
        integer, intent(in) :: spcType
        real(rp), intent(in) :: alam,vm,Dpp !Dpp is plasmasphere density in #/cc
        real(rp) :: tau

        real(rp), parameter :: a3 = -0.113288
        real(rp), parameter :: a2 = +0.659057
        real(rp), parameter :: a1 = +0.319542
        real(rp), parameter :: a0 = +2.16253        
        real(rp), parameter :: day2s = 24.0*60.0*60

        real(rp) :: K,x,y,nTau

        tau = HUGE

        K = abs( alam*vm*(1.0e-3) )!Energy [keV]
        x = log10(K)
        

        if (Dpp < TINY) return
        if (spcType == RAIJUHPLUS) then
            y = a3*(x**3.0) + a2*(x**2.0) + a1*x + a0
            nTau = 10.0**y !Normalized lifetime, days/cm3
            tau = nTau*day2s/Dpp !Lifetime, [s]
            
            tau = max(tau, TINY)
        endif

    end function CCTau


end module raijuLoss_CC