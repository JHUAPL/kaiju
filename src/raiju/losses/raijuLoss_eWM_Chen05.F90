module raijuLoss_eWM_Chen05
    !! Electron wave model blending Chen 2005 and Hiss 2014
    !! Based on Chen et al. 2005 and Orlova et al. 2016

    use kdefs
    use math  ! good idea
    
    use raijudefs
    use raijutypes
    use raijuSpeciesHelper
    use raijuEleLossHelper


    implicit none

    type, extends(baseRaijuLoss_T) :: raiLoss_eWM_C05_T
        logical :: reqsGood = .false.

        ! -- Model -- !
        real(rp) :: NpsphHigh  = 100.0  ! [#/cc]
        real(rp) :: NpsphLow   = 10.0   ! [#/cc]
        
        type(TimeSeries_T) :: KpTS
            !! Kp data from wind file

        real(rp) :: nowT
            !! Current time as of most recent update
        real(rp) :: nowKp
            !! Current kp as of most recent update

        contains

        procedure :: doInit     => eWM_C05_LossInit
        procedure :: isValidSpc => eWM_C05_LossValidSpc
        procedure :: doUpdate   => eWM_C05_LossUpdate
        procedure :: calcTau    => eWM_C05_LossCalcTau

    end type raiLoss_eWM_C05_T

    contains

    subroutine eWM_C05_LossInit(this, Model, Grid, xmlInp)
        class(raiLoss_eWM_C05_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(XML_Input_T) , intent(in) :: xmlInp

        ! Only valid for Earth
        if (trim(toUpper(Model%planet%name)) .ne. "EARTH") then
            write(*,*)"WARNING in raijuLoss_eWM_BW: Not simulating Earth, idk what waves are like elswhere"
            write(*,*)"No electron wave model happening"
            return
        else
            this%reqsGood = .true.
        endif
        this%isPrecip = .true.

        ! We need Kp from wind file
        call xmlInp%Set_Val(this%KpTS%wID,"/Kaiju/Gamera/wind/tsfile","NONE")
        call this%KpTS%initTS("Kp",doLoudO=.false.)

    end subroutine eWM_C05_LossInit


    function eWM_C05_LossValidSpc(this, spc) result(isValid)
        class(raiLoss_eWM_C05_T), intent(in) :: this
        type(raijuSpecies_T), intent(in) :: spc
        logical :: isValid
        
        isValid = .false.
        if (.not. this%reqsGood) return

        if ( (spc%spcType .eq. RAIJUELE) ) then
            isValid = .true.
        endif
    end function eWM_C05_LossValidSpc


    subroutine eWM_C05_LossUpdate(this, Model, Grid, State)
        class(raiLoss_eWM_C05_T), intent(inout) :: this
        type(raijuModel_T ), intent(in) :: Model
        type(raijuGrid_T  ), intent(in) :: Grid
        type(raijuState_T ), intent(in) :: State

        this%nowT = State%t
        this%nowKp = evalTimeSeriesAt(this%KpTS,this%nowT)
    end subroutine eWM_C05_LossUpdate


    function eWM_C05_LossCalcTau(this, Model, Grid, State, i, j, k) result(tau)
        class(raiLoss_eWM_C05_T), intent(in) :: this
        type(raijuModel_T ), intent(in) :: Model
        type(raijuGrid_T  ), intent(in) :: Grid
        type(raijuState_T ), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau

        integer :: psphIdx
        real(rp) :: NpsphPnt
            !! [#/cc] plasmasphere density at point i,j
        real(rp) :: wgtHiss
            !! Weighting factor for hiss contribution
        real(rp) :: L, MLT, E, Kp
        real(rp) :: tauHiss, tauChen

        !psphIdx = spcIdx(Grid, F_PSPH)
        !NpsphPnt = State%Den(psphIdx)%data(i,j)
        !wgtHiss = log(NpsphPnt/this%NpsphLow) / log(this%NpsphHigh/this%NpsphLow)
        !call ClampValue(wgtHiss, 0.0_rp, 1.0_rp)
            !! 1 => Psphere Hiss, 0 => Chen

        !L = sqrt(State%xyzMincc(i,j,XDIR)**2 + State%xyzMincc(i,j,YDIR)**2)  ! [Re]
        !MLT = atan2(State%xyzMincc(i,j,YDIR),State%xyzMincc(i,j,XDIR))/pi*12.D0+12.D0
        !E = abs(Grid%alamc(k)) * State%bvol_cc(i,j)**(-2./3.) * 1.0D-6  ! [MeV]
        !tauHiss = CalcTau_Hiss(MLT, L, E, this%nowKp)  ! [s]

        tauChen = CalcTau_WeakScattering(Model, Grid, State, i, j, k)  ! [s]
        
        !! 1/tau = w_h/tau_h + (1-w_h)/tau_c
        !tau = tauHiss*tauChen / (wgtHiss*tauChen + (1.0_rp - wgtHiss)*tauHiss)  ! [s]

        tau = tauChen  ! test
        call ClampValue(tau, TINY, HUGE)  ! [s]
            !! Clamp to avoid division by zero

    end function eWM_C05_LossCalcTau

end module raijuLoss_eWM_Chen05