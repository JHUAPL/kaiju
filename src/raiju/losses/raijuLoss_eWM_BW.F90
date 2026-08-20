module raijuLoss_eWM_BW
    !! Bao-Wang electron Wave Model - based loss and precipitation model

    use kdefs, only : mec2, vc_cgs
    use raijudefs
    use raijutypes
    use raijuSpeciesHelper
    use raijuLoss_SS

    use arrayutil

    implicit none

    integer, parameter, private :: MAXIOVAR = 50

    type, extends(baseRaijuLoss_T) :: raiLoss_eWM_BW_T
        logical :: reqsGood = .false.

        type(raiLoss_SS_T) :: loss_SS
            !! Our personal instance of the strong scattering implementation to apply where needed

        ! -- Model -- !
        real(rp) :: NpsphHigh  = 100.0  ! [#/cc]
        real(rp) :: NpsphLow   = 10.0   ! [#/cc]
        real(rp) :: ChorusLMax = 7.0  ! [Re]
        real(rp) :: PsheetLMin = 8.0  ! [Re]
        real(rp) :: ChorusEMin = 1.0  ! [keV]
        
        type(TimeSeries_T) :: KpTS
            !! Kp data from wind file

        ! -- Grid -- !
        ! Chorus info
        integer :: Nkp, Nmlt, Nl, Ne
            !! Number of bins for Kp, MLT, L shell, and Energy
        real(rp), allocatable, dimension(:) :: Kp1D
            !! 1D array of Kp dimension for Tau4D
        real(rp), allocatable, dimension(:) :: MLT1D
            !! 1D array of MLT dimension for Tau4D
        real(rp), allocatable, dimension(:) :: L1D
            !! 1D array of L shell dimension for Tau4D [Re]
        real(rp), allocatable, dimension(:) :: Energy1D
            !! 1D array of energy dimension for Tau4D [MeV]
        real(rp), allocatable, dimension(:,:,:,:) :: Tau4D
            !! Tau(Kp, MLT, L, E) table electron scattering table [seconds]

        ! -- State -- !
        real(rp), allocatable, dimension(:,:) :: wPS
        real(rp), allocatable, dimension(:,:) :: wHISS
        real(rp), allocatable, dimension(:,:) :: wCHORUS

        real(rp), allocatable, dimension(:,:,:) :: tauTotal
            !! (i,j,k) Final taus considering all processes and weights
        

        contains
        
        procedure :: doInit => eWM_BW_LossInit
        procedure :: doUpdate => eWM_BW_DoUpdate
        procedure :: isValidSpc => eWM_BW_LossValidSpc
        procedure :: calcTau    => eWM_BW_CalcTau
        procedure :: doOutput => eWM_BW_DoOutput

    end type raiLoss_eWM_BW_T

    contains

    subroutine eWM_BW_LossInit(this, Model, Grid, xmlInp)
        class(raiLoss_eWM_BW_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(XML_Input_T) , intent(in) :: xmlInp
        
        type(raijuSpecies_T) :: spc_ele
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        type(IOVar_T) :: tauVAR

        ! Only valid for Earth
        if (trim(toUpper(Model%planet%name)) .ne. "EARTH") then
            write(*,*)"WARNING in raijuLoss_eWM_BW: Not simulating Earth, idk what waves are like elswhere"
            write(*,*)"No electron wave model happening"
            return
        else
            this%reqsGood = .true.
        endif
        this%isPrecip = .true.

        ! If we are still here then we are gonna be used, set everything up

        spc_ele = Grid%spc(spcIdx(Grid, F_HOTE))

        ! Init our instance of RAIJU strong scattering model
        call this%loss_SS%doInit(Model, Grid, xmlInp)

        ! We need Kp from wind file, try to get that first
        call xmlInp%Set_Val(this%KpTS%wID,"/Kaiju/Gamera/wind/tsfile","NONE")
        call this%KpTS%initTS("Kp",doLoudO=.false.)

        call ClearIO(IOVars)

        if(.not. ioExist(Model%configFName, "Tau", "waveModel")) then
            write(*,*) "This config file not structured for RAIJU, use genRAIJU.py. Good day."
            stop
        endif

        !Chorus wave
        call AddInVar(IOVars,"Kp")  
        call AddInVar(IOVars,"MLT") 
        call AddInVar(IOVars,"L") 
        call AddInVar(IOVars,"Ek") 
        call AddInVar(IOVars,"Tau")
        call ReadVars(IOVars,.false.,Model%configFname, "waveModel")

        this%Nkp  = IOVars(FindIO(IOVars, "Kp" ))%N
        this%Nmlt = IOVars(FindIO(IOVars, "MLT"))%N
        this%Nl   = IOVars(FindIO(IOVars, "L"  ))%N
        this%Ne   = IOVars(FindIO(IOVars, "Ek" ))%N

        tauVAR = IOVars(FindIO(IOVars, "Tau"))

        ! Do some dimension checks
        if (tauVAR%Nr /= 4) then
            write(*,*) "tauDim:",tauVAR%Nr
            write(*,*) 'Currently only support tau model files in the form tau(Kp,MLT,L,Ek)'
            stop
        endif

        if (this%Nkp /=  tauVAR%dims(1) .or. this%Nmlt /=  tauVAR%dims(2) .or. this%Nl /=  tauVAR%dims(3) .or. this%Ne /=  tauVAR%dims(4)) then
            write(*,*) "tauDims:",tauVAR%dims,"Nk:",this%Nkp,"Nm:",this%Nmlt,"Nl:",this%Nl,"Ne:",this%Ne
            write(*,*) 'Dimensions of tau arrays are not compatible'
            stop
        endif

        ! If still here we're good to allocate and store data

        allocate(this%Kp1D    (this%Nkp ))
        allocate(this%MLT1D   (this%Nmlt))
        allocate(this%L1D     (this%Nl  ))
        allocate(this%Energy1D(this%Ne  ))
        allocate(this%Tau4D(this%Nkp,this%Nmlt,this%Nl,this%Ne))

        call IOArray1DFill(IOVars,"Kp" , this%Kp1D    )
        call IOArray1DFill(IOVars,"MLT", this%MLT1D   )
        call IOArray1DFill(IOVars,"L"  , this%L1D     )
        call IOArray1DFill(IOVars,"Ek" , this%Energy1D)
        call IOArray4DFill(IOVars,"Tau", this%Tau4D   )

        call ClearIO(IOVars)
        !Array order check: array is in acsending order
           !Chorus
        if(this%Kp1D(1) > this%Kp1D(this%Nkp)) then
            write(*,*) "Kp: ",this%Kp1D
            write(*,*) "reorder wave model so Kp is in ascending order"
            stop
        end if

        if(this%MLT1D(1) > this%MLT1D(this%Nmlt)) then
            write(*,*) "MLT: ",this%MLT1D
            write(*,*) "reorder wave model so MLT is in ascending order"
            stop
        end if

        if(this%L1D(1) > this%L1D(this%Nl)) then
            write(*,*) "L: ",this%L1D
            write(*,*) "reorder wave model so L shell is in ascending order"
            stop
        end if

        if(this%Energy1D(1) > this%Energy1D(this%Ne)) then
            write(*,*) "Ek: ",this%Energy1D
            write(*,*) "reorder wave model so Ek is in ascending order"
            stop
        end if

        
        allocate(this%wPS    (Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg))
        allocate(this%wHISS  (Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg))
        allocate(this%wCHORUS(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg))

        allocate(this%tauTotal(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, spc_ele%kStart:spc_ele%kEnd))

        
    end subroutine eWM_BW_LossInit


    subroutine eWM_BW_DoUpdate(this, Model, Grid, State)
        !! We are called at the start of every step, once coupling info, boundary conditions, etc. have beenr resolved
        !! We calculate all of the weighting, loss taus, for all electrons, and store it here
        !! When our calcTau gets called for an i,j,k, we'll just report the right value from our State
        class(raiLoss_eWM_BW_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State

        integer :: psphIdx, eleIdx
        logical :: isGood
        integer :: i,j,k
        real(rp) :: NpsphPnt
            !! Density [#/cc] of plasmasphere at point i,j
        real(rp) :: L, MLT, E, Kp, Etemp
            !! L shell and MLT of given point, channel energy in MeV, current Kp, and temporary channel energy in keV
        real(rp) :: wLBlend, wNBlend, tScl
            !! L-weighting of blending between IMAG and PS. 0=PS
            !! Density-weighting between Chorus and Hiss
            !! Scale factor in lifetime between sub-1keV and 1keV+ energy
        real(rp) :: tauPS, tauHiss, tauChorus

        ! Zero everyone out in prep for new values
        call fillArray(this%wPS     , 0.0_rp)
        call fillArray(this%wHISS   , 0.0_rp)
        call fillArray(this%wCHORUS , 0.0_rp)
        call fillArray(this%tauTotal, HUGE  )

        call this%loss_SS%doUpdate(Model, Grid, State)  ! Just in case its doUpdate actually does something in the future

        ! Things that aren't i,j,k-dependent
        psphIdx = spcIdx(Grid, F_PSPH)
        eleIdx = spcIdx(Grid, F_HOTE)
        Kp = evalTimeSeriesAt(this%KpTS,State%t)

        associate(sh=>Grid%shGrid, spc=>Grid%spc(eleIdx))

            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,k,isGood,L,MLT,E,Etemp,tScl,NpsphPnt,wNBlend,wLBlend,tauPS,tauHiss,tauChorus)
            do j=sh%jsg,sh%jeg
                do i=sh%isg,sh%ieg
                    isGood = State%active(i,j) == RAIJUACTIVE
                    if (.not. isGood) then
                        cycle
                    endif

                    L = State%radMax_cc(i,j)  ! [Re]
                    MLT = atan2(State%xyzMincc(i,j,2),State%xyzMincc(i,j,1))/pi*12.D0+12.D0

                    NpsphPnt = State%Den(psphIdx)%data(i,j)  ! [#/cc]
                    
                    ! Calculate blending  
                    wNBlend = dlog(NpsphPnt/this%NpsphLow) / dlog(this%NpsphHigh/this%NpsphLow)
                    call ClampValue(wNBlend, 0.0_rp, 1.0_rp)
                        !! 1 => Psphere Hiss, 0 => Other
                    wLBlend = RampDown(L, this%ChorusLMax, this%PsheetLMin - this%ChorusLMax)
                        !! 1 => Chorus, 0 => PS
                    
                    ! If psphere density is high enough, we always apply Hiss regardless of L
                        ! This means Hiss model needs to appropriately handle any L value its given
                    ! If psphere density is low, we choose between Chorus and plasma sheet strong-scattering based on L
                    ! Also note: the weighting is not k-dependent, but if its not too costly to do for each k then ¯\_(ツ)_/¯
                    
                    ! Calculate weights
                    this%wHISS(i,j) = wNBlend
                    this%wCHORUS(i,j) = (1 - wNBlend)*wLBlend
                    this%wPS(i,j) = (1 - wNBlend)*(1-wLBlend)

                    ! Now calculate loss taus and accumulate
                
                    do k=spc%kStart,spc%kEnd        
                        
                        Etemp = abs(Grid%alamc(k) * State%bvol_cc(i,j)**(-2.0/3.0)) * 1.0E-3  ! [KeV]
                         !Scale up lifetime for sub-1keV energy in Chorus
                         !Calculate E [MeV] to evaluate wave model
                        if (Etemp >= this%ChorusEMin) then  ! [KeV]
                            E = Etemp*1.0E-3_rp !Energy [MeV] 
                            tScl = 1.0_rp !No change needed 
                        else
                            !Define a scaling factor to multiply tau (lifetime)
                            !Lower energy particles are slower which increases lifetime
                            E = this%ChorusEmin*1.0E-3_rp !Energy [MeV]
                            tScl = sqrt(this%ChorusEmin/Etemp)
                        endif

                        if (this%wHISS(i,j) > TINY) then
                            tauHiss = CalcTau_Hiss(MLT, L, E, Kp) * tScl
                        else
                            tauHiss = HUGE
                        endif

                        if (this%wPS(i,j) > TINY) then
                            tauPS = this%loss_SS%calcTau(Model, Grid, State, i, j, k)
                        else
                            tauPS = HUGE
                        endif

                        if (this%wCHORUS(i,j) > TINY) then
                            !! Implement chorus tau calculation here
                            tauChorus = CalcTau_Chorus(this, MLT, L, E, Kp) * tScl
                        else
                            tauChorus = HUGE
                        endif

                        this%tauTotal(i,j,k) = this%wHISS(i,j)*tauHiss + this%wPS(i,j)*tauPS  + this%wCHORUS(i,j)*tauCHORUS 
                    enddo
                enddo
            enddo

        end associate

    end subroutine eWM_BW_DoUpdate


    function eWM_BW_CalcTau(this, Model, Grid, State, i, j, k) result(tau)
        !! We are not doing any heavy calculations since we did it in DoUpdate
        !! Just give back thr right value after some checks
        class(raiLoss_eWM_BW_T), intent(in) :: this
        type(raijuModel_T ), intent(in) :: Model
        type(raijuGrid_T  ), intent(in) :: Grid
        type(raijuState_T ), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau

        associate(spc => Grid%spc(Grid%k2spc(k)))

        tau = HUGE
        if (.not. this%reqsGood) then
            return
        endif

        if (k < spc%kStart .or. k > spc%kEnd) then
            return
        endif

        ! If still here then the requester is right to be calling us, give back our tau
        tau = this%tauTotal(i,j,k)

        end associate

    end function eWM_BW_CalcTau


    function eWM_BW_LossValidSpc(this, spc) result(isValid)
        class(raiLoss_eWM_BW_T), intent(in) :: this
        type(raijuSpecies_T), intent(in) :: spc
        logical :: isValid
        
        isValid = .false.
        if (.not. this%reqsGood) return

        if ( (spc%spcType .eq. RAIJUELE) ) then
            isValid = .true.
        endif
    end function eWM_BW_LossValidSpc


    subroutine eWM_BW_DoOutput(this, Model, Grid, State, gStr, doGhostsO)
        class(raiLoss_eWM_BW_T), intent(in) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        character(len=strLen), intent(in) :: gStr
            !! Group ("Step#X") to append our data to in a sub-group ("eWM")
        logical, intent(in), optional :: doGhostsO


        character(len=strLen) :: gStrWM
        integer :: is, ie, js, je, ks, ke
        logical :: doGhosts
        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .false.
        endif

        if (doGhosts) then
            is = Grid%shGrid%isg
            ie = Grid%shGrid%ieg
            js = Grid%shGrid%jsg
            je = Grid%shGrid%jeg
        else
            is = Grid%shGrid%is
            ie = Grid%shGrid%ie
            js = Grid%shGrid%js
            je = Grid%shGrid%je
        endif

        write(gStrWM, "(A, A)") trim(gStr), "/eWM"
        !Reset IO chain
        call ClearIO(IOVars)

        ! Fat output for precip info
        associate(spc=>Grid%spc(spcIdx(Grid, F_HOTE)))
        call AddOutVar(IOVars, "wPS"    , this%wPS    (is:ie,js:je), dStr="0-1 weighting for plasma sheet strong scattering losses")
        call AddOutVar(IOVars, "wHISS"  , this%wHISS  (is:ie,js:je), dStr="0-1 weighting for plasmasphere hiss losses")
        call AddOutVar(IOVars, "wCHORUS", this%wCHORUS(is:ie,js:je), dStr="0-1 weighting for inner mag chorus wave losses")
        call AddOutVar(IOVars, "tauTotal", this%tauTotal(is:ie,js:je,spc%kStart:spc%kEnd), dStr="[s] combined lifetimes of all contributing loss processes")
        call WriteVars(IOVars,.true.,Model%raijuH5, gStrWM)
        end associate
    end subroutine eWM_BW_DoOutput

    function CalcTau_Chorus(this,MLTin,Lin,Ekin,Kpin) result(tau)
    ! linearly interpolate tau from EWMTauInput to current MLT,L,Kp,Ek value
        class(raiLoss_eWM_BW_T), intent(in) :: this        
        real(rp), intent(in) :: MLTin,Lin,Ekin,Kpin
        real(rp) :: tau
        real(rp) :: tauKMLE(2,2,2,2),tauMLE(2,2,2),tauLE(2,2),tauE(2)! tauKMLE(1,2,2,2) means tauKlMuLuEu, l:lower bound, u: upper bound in the NN method
        real(rp) :: dK,wK,dM,wM,dL,wL,dE,wE
        integer :: kL,kU,mL,mU,lL,lU,eL,eU

        associate(Nm=>this%Nmlt,Nl=>this%Nl,Nk=>this%Nkp,Ne=>this%Ne,&
                  Kpi=>this%Kp1D,MLTi=>this%MLT1D,Li=>this%L1D,Eki=>this%Energy1D,&
                  taui=>this%Tau4D)
        ! Find the nearest neighbors in Kp
        if (Kpin >= maxval(Kpi)) then
            kL = Nk !use Kp maximum 
            kU = Nk
        else if (Kpin <= minval(Kpi)) then
            kL = 1  !use Kp minimum
            kU = 1
        else
            kL = maxloc(Kpi,dim=1,mask=(Kpi<Kpin))
            kU = kL+1
        endif
        ! Find the nearest neighbours in MLT
        if ((MLTin >= maxval(MLTi)) .or. (MLTin <= minval(MLTi)))  then ! maxval of MLT is 24, minval of MLT is 0
            mL = 1 !use MLT = 0
            mU = 1
        else
            mL = maxloc(MLTi,dim=1,mask=(MLTi<MLTin))
            mU = mL+1
        endif
        ! Find the nearest neighbours in L
        if (Lin >= maxval(Li)) then
            lL = Nl !use L maximum
            lU = Nl
        else if (Lin <= minval(Li)) then
            lL = 1 ! use L minimum
            lU = 1
        else
            lL = maxloc(Li,dim=1,mask=(Li<Lin))
            lU = lL+1
        endif
         ! Find the nearest neighbours in Ek
        if (Ekin < minval(Eki)) then
            tau = HUGE ! For low energies, assign a huge lifetime.
            return
        else if (Ekin >= maxval(Eki)) then
            eL = Ne !use Ek maximum
            eU = Ne
        else
            eL = maxloc(Eki,dim=1,mask=(Eki<Ekin))
            eU = eL + 1
        endif
        !linear interpolation in Kp
        if (kL == kU) then
            tauMLE(1,1,1) = log10(taui(kL,mL,lL,eL))!Interpolation in log10(taui) space
            tauMLE(1,1,2) = log10(taui(kL,mL,lL,eU))
            tauMLE(1,2,1) = log10(taui(kL,mL,lU,eL))
            tauMLE(1,2,2) = log10(taui(kL,mL,lU,eU))
            tauMLE(2,1,1) = log10(taui(kL,mU,lL,eL))
            tauMLE(2,1,2) = log10(taui(kL,mU,lL,eU))
            tauMLE(2,2,1) = log10(taui(kL,mU,lU,eL))
            tauMLE(2,2,2) = log10(taui(kL,mU,lU,eU))
        else
            dK = Kpi(kU)-Kpi(kL)
            wK = (Kpin-Kpi(kL))/dK
            tauKMLE(1,1,1,1) = log10(taui(kL,mL,lL,eL))
            tauKMLE(2,1,1,1) = log10(taui(kU,mL,lL,eL))
            tauMLE(1,1,1) = tauKMLE(1,1,1,1) + wK*(tauKMLE(2,1,1,1)-tauKMLE(1,1,1,1))
            tauKMLE(1,1,1,2) = log10(taui(kL,mL,lL,eU))
            tauKMLE(2,1,1,2) = log10(taui(kU,mL,lL,eU))
            tauMLE(1,1,2) = tauKMLE(1,1,1,2) + wK*(tauKMLE(2,1,1,2)-tauKMLE(1,1,1,2))
            tauKMLE(1,1,2,1) = log10(taui(kL,mL,lU,eL))
            tauKMLE(2,1,2,1) = log10(taui(kU,mL,lU,eL))
            tauMLE(1,2,1) = tauKMLE(1,1,2,1) + wK*(tauKMLE(2,1,2,1)-tauKMLE(1,1,2,1))
            tauKMLE(1,1,2,2) = log10(taui(kL,mL,lU,eU))
            tauKMLE(2,1,2,2) = log10(taui(kU,mL,lU,eU))
            tauMLE(1,2,2) = tauKMLE(1,1,2,2) + wK*(tauKMLE(2,1,2,2)-tauKMLE(1,1,2,2))
            tauKMLE(1,2,1,1) = log10(taui(kL,mU,lL,eL))
            tauKMLE(2,2,1,1) = log10(taui(kU,mU,lL,eL))
            tauMLE(2,1,1) = tauKMLE(1,2,1,1) + wK*(tauKMLE(2,2,1,1)-tauKMLE(1,2,1,1))
            tauKMLE(1,2,1,2) = log10(taui(kL,mU,lL,eU))
            tauKMLE(2,2,1,2) = log10(taui(kU,mU,lL,eU))
            tauMLE(2,1,2) = tauKMLE(1,2,1,2) + wK*(tauKMLE(2,2,1,2)-tauKMLE(1,2,1,2))
            tauKMLE(1,2,2,1) = log10(taui(kL,mU,lU,eL))
            tauKMLE(2,2,2,1) = log10(taui(kU,mU,lU,eL))
            tauMLE(2,2,1) = tauKMLE(1,2,2,1) + wK*(tauKMLE(2,2,2,1)-tauKMLE(1,2,2,1))
            tauKMLE(1,2,2,2) = log10(taui(kL,mU,lU,eU))
            tauKMLE(2,2,2,2) = log10(taui(kU,mU,lU,eU))
            tauMLE(2,2,2) = tauKMLE(1,2,2,2) + wK*(tauKMLE(2,2,2,2)-tauKMLE(1,2,2,2))
        end if
        ! linear interpolation in mlt
        if (mL == mU) then
            tauLE(1,1) = tauMLE(2,1,1)
            tauLE(1,2) = tauMLE(2,1,2)
            tauLE(2,1) = tauMLE(2,2,1)
            tauLE(2,2) = tauMLE(2,2,2)
        else
            dM = MLTi(mU)-MLTi(mL)
            wM = (MLTin-MLTi(mL))/dM
            tauLE(1,1) = tauMLE(1,1,1) + wM*(tauMLE(2,1,1)-tauMLE(1,1,1))
            tauLE(1,2) = tauMLE(1,1,2) + wM*(tauMLE(2,1,2)-tauMLE(1,1,2))
            tauLE(2,1) = tauMLE(1,2,1) + wM*(tauMLE(2,2,1)-tauMLE(1,2,1))
            tauLE(2,2) = tauMLE(1,2,2) + wM*(tauMLE(2,2,2)-tauMLE(1,2,2))
        end if
        ! linear interpolation in L
        if (lL == lU) then
            tauE(1) = tauLE(2,1)
            tauE(2) = tauLE(2,2)
        else
            dL = Li(lU)-Li(lL)
            wL = (Lin-Li(lL))/dL
            tauE(1) = tauLE(1,1)+ wL*(tauLE(2,1)-tauLE(1,1))
            tauE(2) = tauLE(1,2)+ wL*(tauLE(2,2)-tauLE(1,2))
        end if
        ! linear interpolation in Ek
        if (eL == eU) then
            tau = tauE(1)
        else
            dE = Eki(eU)-Eki(eL)
            wE = (Ekin-Eki(eL))/dE
            tau = tauE(1) + wE*(tauE(2)-tauE(1))
        end if
        tau = 10.0**tau !convert back to tau in seconds
        end associate

    END FUNCTION CalcTau_Chorus

end module raijuLoss_eWM_BW
