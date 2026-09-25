!Various useful empirical inner magnetosphere models
module imaghelper
    use kdefs
    use kronos
    use earthhelper
    use math
    use geopack
    use raijudefs, ONLY: def_tiote
    
    implicit none

    !Stuff to calculate TM03 plasma sheet model
    !TM03 plasma sheet model, adapted from Matlab code by Meg Noah adapted from TM03

    type TM03_T
    	!Vectors should be in GSM
    	!#/cc,km/s,nPa,nT
        real(rp) :: dSW,vSW,PdynSW,BxSW,BySW,BzSW,MJD0,kT
        logical :: isInit = .false. !Has it been initialized
        type(TimeSeries_T) :: tsD,tsVx,tsVy,tsVz,tsBx,tsBy,tsBz,tsMJD,tsTemp
    end type TM03_T

	type(TM03_T), private :: TM03
    !Coefficients for TM03
    real(rp), dimension(11), parameter, private :: AN = [-0.15930, 0.60805,0.50555,0.07959,0.27462,0.03611, &
    											         -0.03419,-0.79347,1.16222,0.47559,0.71166]

    real(rp), dimension(14), parameter, private :: AP = [ 0.05701,0.52401,0.09075,0.52720,0.07814,-4.42226,-1.53325, &
    											         -1.21666,2.53964,0.31988,0.75433,1.04862,-0.07390,1.01540]
  
    contains

    !Initialize TM03 from solar wind file (wID = filename)
    subroutine InitTM03(wID,t0)
    	character(len=*), intent(in) :: wID
    	real(rp), intent(in) :: t0

    	if (TM03%isInit) return !Already initialized

    	!write(*,*) "Initializing TM03 model ..."

    	!Get time series for TM03 object
    	call SetupTS("D" ,wID,TM03%tsD )
    	call SetupTS("Vx",wID,TM03%tsVx)
    	call SetupTS("Vy",wID,TM03%tsVy)
    	call SetupTS("Vz",wID,TM03%tsVz)
    	call SetupTS("Bx",wID,TM03%tsBx)
    	call SetupTS("By",wID,TM03%tsBy)
    	call SetupTS("Bz",wID,TM03%tsBz)

    	call SetupTS("MJD",wID,TM03%tsMJD)
    	call SetupTS("Temp",wID,TM03%tsTemp)

    	TM03%isInit = .true.
    	call UpdateTM03(t0)

    	contains
    		subroutine SetupTS(vID,fID,tsQ)
    			character(len=*)  , intent(in)    :: vID,fID
    			type(TimeSeries_T), intent(inout) :: tsQ

    			tsQ%wID = trim(fID)
    			call tsQ%initTS(trim(vID),doLoudO=.false.)
    		end subroutine SetupTS
    end subroutine InitTM03

    !Just return current MJD of TM model
    function TM03_MJD()
    	real(rp) :: TM03_MJD

    	TM03_MJD = TM03%MJD0

    end function TM03_MJD

    !Update time for TM03
    subroutine UpdateTM03(t0)
    	real(rp), intent(in) :: t0
    	real(rp) :: D,Vx,Vy,Vz,Bx,By,Bz,T

    	if (.not. TM03%isInit) then
    		write(*,*) "Attempting to update uninitialized TM03 model, bailing ..."
    		stop
    	endif
    	D  = evalTimeSeriesAt(TM03%tsD,t0)        !#/cc
    	Vx = evalTimeSeriesAt(TM03%tsVx,t0)*1.0e-3 !km/s
    	Vy = evalTimeSeriesAt(TM03%tsVy,t0)*1.0e-3 !km/s
    	Vz = evalTimeSeriesAt(TM03%tsVz,t0)*1.0e-3 !km/s
    	Bx = evalTimeSeriesAt(TM03%tsBx,t0)        !nT
    	By = evalTimeSeriesAt(TM03%tsBy,t0)        !nT
    	Bz = evalTimeSeriesAt(TM03%tsBz,t0)        !nT
    	T  = evalTimeSeriesAt(TM03%tsTemp,t0)      !K

    	TM03%MJD0 = evalTimeSeriesAt(TM03%tsMJD,t0)
        !Start storing, want GSM
        TM03%dSW  = D
        TM03%vSW  = norm2([Vx,Vy,Vz])
        TM03%PdynSW = PV2PDyn(TM03%dSW,TM03%vSW)

        !Do B field transform
        call mjdRecalc(TM03%MJD0)
        call SM2GSW(Bx,By,Bz,TM03%BxSW,TM03%BySW,TM03%BzSW)

        TM03%kT = T/(kev2J*J2K) !K=>keV
    end subroutine UpdateTM03

	!Vectors should properly be in GSM
	!sScl is safety factor, >1 means less likely to be called out
	function inShueMP(xyz,sSclO)
		real(rp), intent(in) :: xyz(NDIM)
		logical :: inShueMP
		real(rp), intent(in), optional :: sSclO
		real(rp) :: R,R0,RM,sScl
		
		if (present(sSclO)) then
			sScl = sSclO
		else
			sScl = 1.0
		endif
		R = norm2(xyz)
		call ShueHelper(xyz,R0,RM)
		if (R<RM*sScl) then
			inShueMP = .true.
		else
			inShueMP = .false.
		endif

	end function inShueMP

	function inShueMP_SM(xyzSM,sSclO)
		real(rp), intent(in) :: xyzSM(NDIM)
		logical :: inShueMP_SM
		real(rp), intent(in), optional :: sSclO

		real(rp), dimension(NDIM) :: xyzGSM
		call SM2GSW(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGSM(XDIR),xyzGSM(YDIR),xyzGSM(ZDIR))
		if (present(sSclO)) then
			inShueMP_SM = inShueMP(xyzGSM,sSclO)
		else
			inShueMP_SM = inShueMP(xyzGSM)
		endif
		
	end function inShueMP_SM

	!Map xyz coordinates to terminator (dusk/dawn) along Shue MP contours
	subroutine ShueMP2Terminator(xyz,xyzD)
		real(rp), intent(in)  :: xyz (NDIM)
		real(rp), intent(out) :: xyzD(NDIM)

		real(rp) :: R,theta,CosT,R0,alpha,rScl
		xyzD = 0.0

		R = norm2(xyz)
		theta = acos(xyz(XDIR)/R)
		CosT = cos(theta)
		alpha = (0.58-0.007*TM03%BzSW)*(1.0+0.024*log(TM03%PdynSW))

		!Solve for R0
		rScl = (2.0/max(1.0+CosT,TINY))**(alpha)
		R0 = R/rScl

		if (xyz(YDIR)>0) then
			xyzD = [0.0_rp,+R0*(2.0**alpha),xyz(ZDIR)] !Dusk
		else
			xyzD = [0.0_rp,-R0*(2.0**alpha),xyz(ZDIR)] !Dawn
		endif
	end subroutine ShueMP2Terminator

	!Map xyz coordinates to dusk along Shue MP contours
	subroutine ShueMP2Dusk(xyz,xyzD)
		real(rp), intent(in)  :: xyz (NDIM)
		real(rp), intent(out) :: xyzD(NDIM)

		real(rp) :: R,theta,CosT,R0,alpha,rScl
		xyzD = 0.0

		R = norm2(xyz)
		theta = acos(xyz(XDIR)/R)
		CosT = cos(theta)
		alpha = (0.58-0.007*TM03%BzSW)*(1.0+0.024*log(TM03%PdynSW))

		!Solve for R0
		rScl = (2.0/max(1.0+CosT,TINY))**(alpha)
		R0 = R/rScl

		xyzD = [0.0_rp,R0*(2.0**alpha),xyz(ZDIR)]

	end subroutine ShueMP2Dusk

	subroutine ShueHelper(xyz,R0,RM)
		real(rp), intent(in ) :: xyz (NDIM)
		real(rp), intent(out) :: R0,RM
		real(rp) :: R,theta,CosT,alpha

		R = norm2(xyz)
		theta = acos(xyz(XDIR)/R)
		CosT = cos(theta)
		
		R0 = (10.22+1.29*tanh(0.184*(TM03%BzSW+8.14)))*(TM03%PdynSW)**(-1.0/6.6)
		alpha = (0.58-0.007*TM03%BzSW)*(1.0+0.024*log(TM03%PdynSW))
		RM = R0*(2.0/max(1.0+CosT,TINY))**(alpha)

	end subroutine ShueHelper

	!Return density [#/cc] and pressure [nPa] from TM03
	!Assuming XYZ are in GSM
	subroutine EvalTM03(xyzIN,D,P,isIn)
		real(rp), intent(in)    :: xyzIN(NDIM)
		real(rp), intent(out)   :: D,P
		logical , intent(out)   :: isIn

		real(rp) :: R,clockang,Fst,Pst,rho,rho10,Nst,BNst,BSst,phi
		real(rp) :: xyz(NDIM)

		D = 0.0
		P = 0.0
		isIn = .false.

		rho = norm2([xyzIN(XDIR),xyzIN(YDIR)])

		if (xyzIN(XDIR)>0) then
			!Map back to Shue MP
			call ShueMP2Terminator(xyzIN,xyz)
		else
			xyz = xyzIN
		endif

		isIn = inShueMP(xyz) .and. (rho>10) .and. (xyz(XDIR)<=0) &
			   .and. (xyz(XDIR)>-50)
		if (.not. isIn) return

		!----
		!If still here, get to work
		clockang = katan2(TM03%BySW,TM03%BzSW)
		Fst = norm2([TM03%BySW,TM03%BzSW])*sin(clockang/2)/5.0
		Pst = TM03%PdynSW/3.0
		
		rho10 = rho/10.0
		Nst = TM03%dSW/10.0
		BNst = TM03%BzSW/5.0
		if (BNst<0) BNst = 0.0
		BSst = (-TM03%BzSW/5.0)*(TM03%vSW/500.0)
		if (BSst<0) BSst = 0.0

		phi = katan2(xyz(YDIR),-xyz(XDIR))

		!Get TS D [#/cc], P [nPa]
		D = (AN(1) + AN(2)*(Nst**AN(10)) + AN(3)*BNst + AN(4)*BSst)*(RHO10**AN(8)) + &
		    (AN(5)*(Nst**AN(11)) + AN(6)*BNst + AN(7)*BSst)*(RHO10**AN(9))*(sin(phi)**2.0)


		P =  AP(1)*(RHO10**AP(6)) + &
		     AP(2)*(Pst**AP(11))*(RHO10**AP(7)) + &
		     AP(3)*(Fst**AP(12))*(RHO10**AP(8)) + &
		    (AP(4)*(Pst**AP(13))*exp(-AP( 9)*RHO10) + &
		     AP(5)*(Fst**AP(14))*exp(-AP(10)*RHO10))*(sin(phi)**2.0)

	end subroutine EvalTM03

	!Evaluate for SM coordinates incoming
	subroutine EvalTM03_SM(xyzSM,D,P,isIn)
		real(rp), intent(in)    :: xyzSM(NDIM)
		real(rp), intent(out)   :: D,P
		logical , intent(out)   :: isIn

		real(rp), dimension(NDIM) :: xyzGSM

		call SM2GSW(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGSM(XDIR),xyzGSM(YDIR),xyzGSM(ZDIR))

		call EvalTM03(xyzGSM,D,P,isIn)

	end subroutine EvalTM03_SM

	!Returns electron temperature in keV based on DGSR2016 plasma sheet model
	!10.1002/2016JA022947
	!Adapted from code written by Sina and provided by Shanshan
	function Tps_dgsr2016(xgsm,ygsm,zgsm) result(Tps)
		real(rp), intent(in) :: xgsm,ygsm,zgsm
		real(rp) :: Tps
		real(rp) :: A1,A2,A3,A4,A5,A6,A7,A8,A9,phi0
		real(rp) :: rho,phi,rho_s,phi_s,n_sw_s,v_sw_s,imf_bs_s,imf_bn_s,Q
		A1 = -0.0215
		A2 = -0.426 
		A3 =  1.47  
		A4 =  0.587 
		A5 = -0.538 
		A6 = -0.489 
		A7 =  0.32  
		A8 =  0.36  
		A9 =  2.31  
		rho = sqrt(xgsm**2.0 + ygsm**2.0 + zgsm**2.0)
		phi0 = katan2(ygsm,xgsm) !Azimuth relative to +X
		!Note: Need funny phi w/ midnight = 0deg, dusk = +90deg, dawn = -90deg
		phi = PI - phi0
		!Before scaling, do some clamping
		call ClampValue(rho,6.0_rp,10.0_rp) !Radial bounds for DGSR
		call ClampValue(phi,-PI/2,+PI/2) !Force on nightside
		rho_s = rho / 10.0
		phi_s = 2.0*phi/PI
		n_sw_s = TM03%dSW/10.0
		v_sw_s = TM03%vSW/400.0 
		if (TM03%BzSW<0) then
			imf_bs_s = -TM03%BzSW/2.0
			imf_bn_s = 0.0
		else
			imf_bs_s = 0.0
			imf_bn_s = TM03%BzSW/2.0
		endif
 		Q = A1 + A2*phi_s + A3*v_sw_s + (A4+A5*phi_s*phi_s*rho_s)*(imf_bs_s**A7) &
 		    + A6*rho_s*(imf_bn_s**A8)
 		if (Q > TINY) then
 			Tps = Q**A9 !Temperature in keV
 		else
 			Tps = TINY
 		endif
 	end function Tps_dgsr2016

 	!Get empirical Ti over Te from xyzGSM
 	!Using TM03 and DSGR16 respectively, creating MLT profile at R=10
 	!NOTE: This routine clamps to nightside
 	function GSM2TioTe_Empirical(xyzGSM) result(TioTe)
 		real(rp), intent(in)    :: xyzGSM(NDIM)
		real(rp) :: TioTe,ionD,ionP,ionT,eleT,phi
		real(rp) :: R0 = 10.0,rEps=0.1,rho
		real(rp), dimension(NDIM) :: xyzGSM_TM03,xyzGSM_DSGR16
		logical :: isIn
		phi = katan2(xyzGSM(YDIR),xyzGSM(XDIR)) !0,2pi
		!Ensure we're on the nightside
		call ClampValue(phi,PI/2+TINY,3*PI/2-TINY)
		rho = norm2(xyzGSM)
		!Now construct a point to evaluate TM03
		if (rho < (R0+rEps)) then
			xyzGSM_TM03   = (R0+rEps)*[cos(phi),sin(phi),0.0_rp]
		else
			xyzGSM_TM03   = xyzGSM
		endif
		!Evaluate DSGR16
		if (rho > (R0-rEps)) then
			xyzGSM_DSGR16 = (R0-rEps)*[cos(phi),sin(phi),0.0_rp]
		else
			xyzGSM_DSGR16 = xyzGSM
		endif
		call EvalTM03(xyzGSM_TM03,ionD,ionP,isIn)
		if (.not. isIn) then
			TioTe = def_tiote
		else
			ionT = DP2kT(ionD,ionP)
			eleT = Tps_dgsr2016(xyzGSM_DSGR16(XDIR),xyzGSM_DSGR16(YDIR),xyzGSM_DSGR16(ZDIR))
			
			if ( (ionT<TINY) .or. (eleT<TINY) ) then
				TioTe = def_tiote
			else
				TioTe = ionT/eleT
			endif
		endif
		!Clamp to physically reasonable values
		call ClampValue(TioTe,2.0_rp,7.0_rp) !See Runov++ 2015
 	end function GSM2TioTe_Empirical
 	
 	!Map SM to TioTe w/ MLT interpolation on the dayside
 	function TioTe_Empirical(xyzSM) result(TioTe)
 		real(rp), intent(in)    :: xyzSM(NDIM)
 		real(rp) :: TioTe,rho,phi
		real(rp) :: TioTe_Dawn,TioTe_Dusk
		real(rp), dimension(NDIM) :: xyzGSM, xyzDawn, xyzDusk
 		call SM2GSW(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGSM(XDIR),xyzGSM(YDIR),xyzGSM(ZDIR))
 		phi = atan2(xyzGSM(YDIR),xyzGSM(XDIR)) !-pi,+pi
 		rho = norm2(xyzGSM)
 		if ( xyzGSM(XDIR) < 0 ) then
 			!On nightside, just evaluate directly
 			TioTe = GSM2TioTe_Empirical(xyzGSM)
 		else
 			!On dayside, evaluate at dawn/dusk
 			xyzDawn = [0.0_rp,-rho,0.0_rp]
 			xyzDusk = [0.0_rp,+rho,0.0_rp]
 			Tiote_Dawn = GSM2TioTe_Empirical(xyzDawn)
 			Tiote_Dusk = GSM2TioTe_Empirical(xyzDusk)
 			TioTe = Tiote_Dawn + LinRampUp(phi,-PI/2,PI)*(Tiote_Dusk-Tiote_Dawn)
 		endif
 			
 	end function TioTe_Empirical

	!> Plasmasphere temperature, use Genestreti+ 2016 linear fit
    !> Density [#/cc] => Temperature [keV]
	function PsphTemp_Genestreti(npsph) result(Tkev)
		real(rp), intent(in) :: npsph
			!! Density [#/cc]
		real(rp) :: Tkev
			!! Temperature [keV]
		real(rp), parameter :: A = -0.15,B=1.4
		!T [eV] = B* n^A, n [#/cc], A = -0.15, B = 1.4
		!Floor at 0.01/cc
		Tkev = (1.0e-3)*B*(max(npsph,0.01)**A)
	end function PsphTemp_Genestreti

 end module imaghelper
