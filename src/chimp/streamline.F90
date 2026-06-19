module streamline
    use chmpdefs
    use ebtypes
    use math
    use ebinterp
    use planethelper
    use streamutils
    use imaghelper

    implicit none

    real(rp), private :: ShueScl = 2.0  !Safety factor for Shue MP
    real(rp), private :: rShue   = 6.0  !Radius to start checking Shue
    integer , private :: NpChk   = 10   !Cadence for Shue checking
    logical , private :: doShueG = .false. !Global doShue, can be overriden by optional arguments to routines

    contains

    subroutine setShue(Model,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: inpXML

        if (Model%isMAGE .and. (trim(toUpper(Model%uID)) == "EARTHCODE")) then
            !This is for Earth and we're running in tandem w/ mage
            call inpXML%Set_Val(ShueScl,'streamshue/ShueScl' ,ShueScl)
            call inpXML%Set_Val(rShue  ,'streamshue/rShue'   ,rShue  )
            call inpXML%Set_Val(NpChk  ,'streamshue/NpChk'   ,NpChk  )
            call inpXML%Set_Val(doShueG,'streamshue/doShue' ,.false.)
            if (doShueG) write(*,*) "Initializing SHUE-MP checking ..."
        else
            rShue = HUGE
            doShueG = .false.
        endif

    end subroutine setShue

    !Trace field line w/ seed point x0 (both directions)
    !doShueO = Use Shue MP model (w/ safety factor) for short circuiting

    !doNHO = T, assume doing RCM coupling
    subroutine genLine(Model,ebState,x0,t,fL,MaxStepsO,doShueO,doNHO,doSHO)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        type(magLine_T)  , intent(inout) :: fL
        integer , intent(in), optional :: MaxStepsO
        logical , intent(in), optional :: doShueO,doNHO,doSHO

        real(rp), allocatable :: xyzn(:,:,:), nBxyz(:,:,:), nGas(:,:,:,:)
        integer , allocatable :: ijkn(:,:,:)

        logical :: inDom,doShue,doNH,doSH
        integer :: N1,N2,i,MaxN,n0

        !Start by empting line
        call cleanLine(fL)

        !Do some pre-checking
        fL%x0 = x0
        inDom = inDomain(x0,Model,ebState%ebGr)
        if (.not. inDom) then
            !Bad tube, seed point isn't in domain
            fL%isGood = .false.
            allocate(fL%xyz(0:0,NDIM))
            fL%xyz(0,:) = x0
            return

        endif

        !Okay, we're still here so let's do this thing
        !Start by dealing w/ optional and defaults
        if (present(MaxStepsO)) then
            MaxN = min(MaxStepsO, MaxFL)
        else
            MaxN = MaxFL
        endif
        if ((present(doNHO))) then
            doNH = doNHO
        else
            doNH = .false.
        endif

        if ((present(doSHO))) then
            doSH = doSHO
        else
            doSH = .false.
        endif

        !If optional argument, do whatever it says. Otherwise default to global value
        if (present(doShueO)) then
            doShue = doShueO
        else
            doShue = doShueG
        endif

        !Allocate temp arrays to hold information along each direction
        !Trailing dimension (1-2) represents -/+ directions
        allocate(xyzn (0:MaxFL,NDIM,2))
        allocate(ijkn (0:MaxFL,NDIM,2))
        allocate(nBxyz(0:MaxFL,NDIM,2))
        allocate(nGas (0:MaxFL,NVARMHD,0:Model%nSpc,2))

    !Get traces and store in temp arrays
        n0 = 1 !Default to forward side
        !Get trace in negative direction
        if (doSH) then
            !Starting in southern hemisphere so don't go in negative direction
            N1 = 0
            n0 = 2 !Can't use forward side later
        else
            call genTrace(Model,ebState,x0,t,xyzn(:,:,1),ijkn(:,:,1),nBxyz(:,:,1),nGas(:,:,:,1),N1,-1,MaxN,doShue)
        endif

        !Get trace in positive direction
         if (doNH) then
            !We're starting at northern hemisphere so don't need to go in positive direction
            N2 = 0
            n0 = 1
        else
            call genTrace(Model,ebState,x0,t,xyzn(:,:,2),ijkn(:,:,2),nBxyz(:,:,2),nGas(:,:,:,2),N2,+1,MaxN,doShue)
        endif

    !Now create field line and scrape values out of temp arrays
        fL%nMax = MaxN
        fL%isGood = .true.
        fL%x0 = x0
        fL%Nm = N1
        fL%Np = N2
    
        !Do allocations
        allocate(fL%xyz (-N1:+N2,NDIM))
        allocate(fL%ijk (-N1:+N2,NDIM))
        allocate(fL%Bxyz(-N1:+N2,NDIM))
        allocate(fL%magB(-N1:+N2))
        allocate(fL%Gas (-N1:+N2,NVARMHD,0:Model%nSpc) ) !0 = BULK

        !Load field lines
        fL%xyz (0,:)   = x0
        fL%ijk (0,:)   = ijkn (0,:,n0)
        fL%Bxyz(0,:)   = nBxyz(0,:,n0)
        fL%magB(0)     = norm2(fL%Bxyz(0,:))
        fL%Gas (0,:,:) = nGas (0,:,:,n0)

        do i=1,N1 !Negative direction
            fL%xyz (-i,:)   = xyzn (i,:,1)
            fL%ijk (-i,:)   = ijkn (i,:,1)
            fL%Bxyz(-i,:)   = nBxyz(i,:,1)
            fL%magB(-i)     = norm2(fL%Bxyz(-i,:))
            fL%Gas (-i,:,:) = nGas(i,:,:,1)
        enddo

        do i=1,N2 !positive direction
            fL%xyz (+i,:)   = xyzn (i,:,2)
            fL%ijk (+i,:)   = ijkn (i,:,2)
            fL%Bxyz(+i,:)   = nBxyz(i,:,2)
            fL%magB(+i)     = norm2(fL%Bxyz(+i,:))
            fL%Gas (+i,:,:) = nGas(i,:,:,2)
        enddo

    end subroutine genLine

    !Gathers field line topology information
    subroutine SliceFL(Model,ebState,x0,t,ebTrc)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        type(ebTrc_T)    , intent(inout) :: ebTrc

        type(magLine_T) :: bTrc

        !Initialize the values
        ebTrc%OCb = 0.0
        ebTrc%dvB = 0.0
        ebTrc%bD  = 0.0
        ebTrc%bP  = 0.0
        ebTrc%bS  = 0.0
        ebTrc%bMin = 0.0
        ebTrc%stdP = 0.0
        ebTrc%stdD = 0.0

        ebTrc%MagEQ(:) = 0.0
        ebTrc%xEPm (:) = 0.0
        ebTrc%xEPp (:) = 0.0

        if (.not. inDomain(x0,Model,ebState%ebGr)) return
        !Trace field line
        call genLine(Model,ebState,x0,t,bTrc)

        !Get diagnostics
        ebTrc%OCb = 1.0*FLTop(Model,ebState%ebGr,bTrc)
        if (ebTrc%OCb > 0) then
            !Get flux-tube integrals
            call FLThermo(Model,ebState%ebGr,bTrc,ebTrc%bD,ebTrc%bP,ebTrc%dvB)
            call FLStdev(Model, ebState%ebGr,bTrc,ebTrc%bD,ebTrc%bP,ebTrc%stdD,ebTrc%stdP)
            !call FLStdev(Model, ebState%ebGr,bTrc,ebTrc%stdD,ebTrc%stdP)
            ebTrc%bS   = FLEntropy(Model,ebState%ebGr,bTrc)

            !Get magnetic equator info
            call FLEq(Model,bTrc,ebTrc%MagEQ,ebTrc%bMin)

            !Get endpoints info
            associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
            ebTrc%xEPm = bTrc%xyz(-Nm,:)
            ebTrc%xEPp = bTrc%xyz(+Np,:)
            end associate

        endif

    end subroutine SliceFL

!---------------------------------
!Field line diagnostics
    !Flux tube volume
    function FLVol(Model,ebGr,bTrc) result(dvB)
        type(chmpModel_T)  , intent(in) :: Model
        type(ebGrid_T)     , intent(in) :: ebGr
        type(magLine_T)    , intent(in) :: bTrc
        real(rp) :: dvB

        integer :: OCb,k
        real(rp) :: Phi0,dl,bAvg,dA

        Phi0 = 1.0 !Fixed flux, could be B0*A0
        dvB = 0.0

        OCb = FLTop(Model,ebGr,bTrc)
        if (OCb<2) return
        
        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
        do k=-Nm,Np-1
            dl   = norm2(bTrc%xyz (k+1,:)-bTrc%xyz (k,:))
            bAvg = 0.5* (bTrc%magB(k+1) + bTrc%magB(k))

            dA = Phi0/bAvg
            dvB = dvB + dA*dl
        enddo

        end associate

    end function FLVol

    !Calculate arc length of field line
    function FLArc(Model,ebGr,bTrc) result(L)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        real(rp) :: L
        real(rp) :: dL
        integer :: k

        L = 0.0
        if (.not. bTrc%isGood) return

        do k=-bTrc%Nm,bTrc%Np-1
            dL = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            L = L + dL
        enddo

    end function FLArc

    !Calculate Alfven crossing time on line
    function FLAlfvenX(Model,ebGr,bTrc) result(dtX)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        real(rp) :: dtX

        integer :: k
        real(rp) :: dL,eD,bMag,Va

        dtX = 0.0
        if (.not. bTrc%isGood) return

        do k=-bTrc%Nm,bTrc%Np-1
            dL = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            dL = dL*L0*1.0e-5 !Corner units to km
            !Get egde-centered quantities, use BULK
            eD   = 0.5*( bTrc%Gas(k+1,DEN,BLK) + bTrc%Gas(k,DEN,BLK) )
            bMag = 0.5* (bTrc%magB(k+1) + bTrc%magB(k))

            !Convert B to nT, eD in #/cc
            bMag = oBScl*bMag
            Va = 22.0*bMag/sqrt(eD) !Alfven speed in km/s, NRL formulary
            dtX = dtX + dL/Va
        enddo

    end function FLAlfvenX

    !Calculate FAC in the ionosphere from the MHD inner boundary
    function FLIonFAC(Model,ebGr,ebState,bTrc,Rion,i0) result(J11i)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(ebState_T)  , intent(in) :: ebState
        type(magLine_T  ), intent(in) :: bTrc
        real(rp)         , intent(in) :: Rion
        integer          , intent(in) :: i0

        real(rp) :: J11i
        real(rp) :: t,J11
        real(rp), dimension(NDIM) :: B0,bhat,Jxyz,xyzMHD

        J11i = 0.0 !Placeholder
        if (.not. Model%doJ) return !Nothing to be done

        call DipoleShift2Shell(Model,ebGr,i0,bTrc%x0,xyzMHD)
        t = ebState%eb1%time

        Jxyz = jInterp(xyzMHD,t,Model,ebState) !Evaluate J @ xyzMHD
        B0 = Model%B0(xyzMHD)
        bhat = normVec(B0)

        J11 = dot_product(Jxyz,bhat) !A/m2 J11 @ xyzMHD
        !Now get scaling factor
        J11i = Bi2M(xyzMHD,Rion)*J11*(1.0e+6) !Rescale to ionosphere, conserve J/B, and convert to muA/m2
    end function FLIonFAC

    !Return Bxyz and Jxyz at i0 shell
    !NOTE: Jxyz includes Bi2m scaling factor already
    subroutine FLIonBJ(Model,ebGr,ebState,bTrc,Rion,i0,Bxyz,Jxyz)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(ebState_T)  , intent(in) :: ebState
        type(magLine_T  ), intent(in) :: bTrc
        real(rp)         , intent(in) :: Rion
        integer          , intent(in) :: i0
        real(rp), dimension(NDIM), intent(out) :: Bxyz,Jxyz

        real(rp) :: t,xyzMHD(NDIM)
        integer :: ijk(NDIM)
        logical :: isIn

        Bxyz = 0.0
        Jxyz = 0.0
        if (.not. Model%doJ) return !Nothing to be done

        call DipoleShift2Shell(Model,ebGr,i0,bTrc%x0,xyzMHD)
        t = ebState%eb1%time
        call locate(xyzMHD,ijk,Model,ebGr,isIn)
        !Evaluate B & J @ xyzMHD
        Bxyz = fldInterp(xyzMHD,t,Model,ebState,BFLD,isIn,ijk)
        Jxyz = jInterp(xyzMHD,t,Model,ebState,ijk) 

        Bxyz = Bxyz*oBScl !EB => nT
        !Rescale to ionosphere, conserve J/B, and convert to muA/m2
        Jxyz = (1.0e+6)*Bi2M(xyzMHD,Rion)*Jxyz

    end subroutine FLIonBJ

    function Bi2M(xyz,Rion)
        real(rp), intent(in) :: xyz(NDIM), Rion 
        real(rp) :: Bi2M

        real(rp) :: R0,Rx,zor2,mZ,pZ
        R0 = norm2(xyz)
        Rx = R0/Rion
        zor2 = (xyz(ZDIR)/R0)**2.0
        mZ = (1.0 - zor2)/Rx
        pZ = 1.0 + 3.0*zor2
        Bi2M = (Rx**3.0)*sqrt( 1.0 + 3.0*(1.0-mZ) )/sqrt(pZ)
    end function Bi2M

    !Get plasma information from near MHD inner boundary from ionosphere
    subroutine FLIonDP(Model,ebGr,ebState,bTrc,bD,bP,i0,sOpt) 
        type(chmpModel_T), intent(in ) :: Model
        type(ebGrid_T)   , intent(in ) :: ebGr
        type(ebState_T)  , intent(in ) :: ebState
        type(magLine_T  ), intent(in ) :: bTrc
        integer          , intent(in ) :: i0
        real(rp)         , intent(out) :: bD,bP
        integer          , intent(in ), optional :: sOpt

        integer :: s0
        real(rp) :: t
        real(rp), dimension(NDIM) :: xyzMHD
        
        bD = 0.0
        bP = 0.0
        if (present(sOpt)) then
            if (sOpt > Model%nSpc) then
                !Got nothing for ya
                return
            else
                s0 = sOpt
            endif
        else
            s0 = BLK
        endif

        call DipoleShift2Shell(Model,ebGr,i0,bTrc%x0,xyzMHD)
        t = ebState%eb1%time
        call mhdInterpDP(xyzMHD,t,Model,ebState,bD,bP,s0)

    end subroutine FLIonDP

    !Take point xyz0 and dipole shift to radius of the i0-th shell
    subroutine DipoleShift2Shell(Model,ebGr,i0,xyz0,xyzMHD)
        type(chmpModel_T), intent(in)  :: Model
        type(ebGrid_T)   , intent(in)  :: ebGr
        integer          , intent(in)  :: i0
        real(rp)         , intent(in)  :: xyz0(NDIM)
        real(rp)         , intent(out) :: xyzMHD(NDIM)

        real(rp) :: R0
        R0 = norm2(ebGr%xyzcc(ebGr%is+i0-1,ebGr%js,ebGr%ks,:))
        xyzMHD = DipoleShift(xyz0,R0) !Dipole shift to i0 shell radius

        !Check we went far out enough and die hard if we didn't
        if (inGap(xyzMHD,Model,ebGr)) then
            write(*,*) "Mind the gap!"
            stop
        endif      
    end subroutine DipoleShift2Shell

    !Averaged density/pressure along a segment of the field line
    subroutine FL2iBC(Model,ebGr,bTrc,bD,bP,i0,L,sOpt)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        real(rp), intent(out) :: bD,bP
        integer , intent(in ) :: i0
        real(rp), intent(in ) :: L
        integer , intent(in ), optional :: sOpt 

        integer  :: s0,k,kS,kE,dk
        real(rp) :: bMag,dl,dvB,dX,eD,eP

        bD = 0.0
        bP = 0.0
        if (.not. bTrc%isGood) return

        if (present(sOpt)) then
            s0 = min(sOpt,Model%nSpc)
        else
            s0 = BLK
        endif

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)

        !Figure out which direction to go
        dX  = 0.0
        dvB = 0.0
        if (Np > Nm) then
            !More points in + direction, this is a SH seed
            kS = 0
            kE = Np-1
            dk = +1
        else
            !More points in - direction, this is a NH seed
            kS = 0
            kE = -Nm+1
            dk = -1
        endif

        !Loop over edges
        do k=kS,kE,dk
            !Get edge-centered quantities
            dl = norm2(bTrc%xyz(k+dk,:) - bTrc%xyz(k,:)) !Edge length

            eD   = 0.5*( bTrc%Gas(k+dk,DEN     ,s0) + bTrc%Gas(k,DEN     ,s0) )
            eP   = 0.5*( bTrc%Gas(k+dk,PRESSURE,s0) + bTrc%Gas(k,PRESSURE,s0) )

            bMag = 0.5* (bTrc%magB(k+dk) + bTrc%magB(k))
            !Now accumulate into flux-tube integrals if in region of interest
            if (bTrc%ijk(k,IDIR) >= i0 .and. bTrc%ijk(k+dk,IDIR) >= i0) then
                bD  = bD  +  eD*dl/bMag
                bP  = bP  +  eP*dl/bMag
                dvB = dvB + dl/bMag
                dX = dX + dl
            endif
            if (dX >= L) exit !Break out of loop
        enddo
        dvB = max(dvB,TINY)
        bD  = bD/dvB  
        bP  = bP/dvB        
        
        end associate
    end subroutine FL2iBC

    !Averaged density/pressure (of species sOpt or BLK)
    subroutine FLThermo(Model,ebGr,bTrc,bD,bP,dvB,bBetaO,sOpt)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        real(rp), intent(out) :: bD,bP,dvB
        real(rp), intent(out), optional :: bBetaO
        integer , intent(in ), optional :: sOpt 

        integer :: k,s0
        real(rp) :: dvB_active,bMag,dl,eP,eD,ePb !Edge-centered values
        real(rp) :: bPb,bBeta

        !Zero out accumulators
        bD = 0.0
        bP = 0.0
        dvB = 0.0
        bPb = 0.0
        bBeta = 0.0
        dvB_active = 0.0
        
        if (.not. bTrc%isGood) return

        if (present(sOpt)) then
            s0 = min(sOpt,Model%nSpc)
        else
            s0 = BLK
        endif

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)

        !Loop over edges
        do k=-Nm,Np-1
            !Get edge-centered quantities
            dl = norm2(bTrc%xyz(k+1,:) - bTrc%xyz(k,:)) !Edge length

            eD   = 0.5*( bTrc%Gas(k+1,DEN     ,s0) + bTrc%Gas(k,DEN     ,s0) )
            eP   = 0.5*( bTrc%Gas(k+1,PRESSURE,s0) + bTrc%Gas(k,PRESSURE,s0) )

            bMag = 0.5* (bTrc%magB(k+1) + bTrc%magB(k))

            !Get edge mag pressure, bmag=>nT(oBScl)=>T
            !(NRL Plasma formulary):
            !3.98x10^6 * (B/B0)^2 = Pb [dynes/cm2] = 0.1 Pa, x10^8 0.1 Pa => nPa
            !ePb = (1.0e+8)*(3.98*1.0e+6)*(bMag*oBScl*1.0e-9)**2.0
            ePb = 1.0e+14*(bMag*oBScl*1.0e-9/0.501)**2.0 !Edge mag pressure in nPa

            !Now accumulate into flux-tube integrals
            ! Always accumulate total flux tube volume
            dvB = dvB + dl/bMag
            ! Only accumulate other quantities if in active domain, so we have more sensible flux tube averages

            if (bTrc%ijk(k,IDIR) >= ebGr%is .and. bTrc%ijk(k+1,IDIR) >= ebGr%is) then
                bD  = bD  +  eD*dl/bMag
                bP  = bP  +  eP*dl/bMag
                bPb = bPb + ePb*dl/bMag
                bBeta = bBeta + (eP/ePb)*dl/bMag
                dvB_active = dvB_active + dl/bMag
            endif
        enddo

        dvB_active = max(dvB_active, TINY)  ! Max just do prevent nans for lines that were completely under the gap+padding
        !Now turn flux-tube integrals of quantities into flux-tube averages
        bD  = bD /dvB_active  
        bP  = bP /dvB_active
        bPb = bPb/dvB_active

        !bBeta = bP/bPb
        bBeta = bBeta/dvB_active


        if (present(bBetaO)) then
            bBetaO = bBeta
        endif

        end associate
    end subroutine FLThermo


    subroutine FLStdev(Model, ebGr, bTrc, bD, bP, stdD, stdP, sOpt)
    !subroutine FLStdev(Model, ebGr, bTrc, stdD, stdP, sOpt)
        !! Computes volume-weighted standard deviation of thermo quantities
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        real(rp), intent(in)  :: bD, bP
        !    !! Field line averages calculated from FLThermo
        real(rp), intent(out) :: stdD, stdP
            !! Standard deviation of pressure and density we calculate and return
        integer , intent(in ), optional :: sOpt 
            !! Which fluid to use

        integer :: k,s0
        real(rp) :: bMag,dl,eP,eD !Edge-centered values
        real(rp) :: dvB
        real(rp), dimension(NDIM) :: xyzTmp

        !Zero out accumulators
        stdP = 0.0
        stdD = 0.0
        dvB = 0.0
        xyzTmp = 0.0
        
        if (.not. bTrc%isGood) return

        if (present(sOpt)) then
            s0 = min(sOpt,Model%nSpc)
        else
            s0 = BLK
        endif

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)


        ! Now to stdev
        !Loop over edges
        do k=-Nm,Np-1
            if ( (bTrc%ijk(k,XDIR) .le. ebGr%is) .or. (bTrc%ijk(k+1,XDIR) .le. ebGr%is) ) then
                cycle
            endif
            !Get edge-centered quantities
            dl = norm2(bTrc%xyz(k+1,:) - bTrc%xyz(k,:)) !Edge length

            eD = 0.5*( bTrc%Gas(k+1,DEN     ,s0) + bTrc%Gas(k,DEN     ,s0) )
            eP = 0.5*( bTrc%Gas(k+1,PRESSURE,s0) + bTrc%Gas(k,PRESSURE,s0) )
            bMag = 0.5* (bTrc%magB(k+1) + bTrc%magB(k))

            dvB = dvB + dl/bMag
            ! Perform volume-weighted summations
            stdP = stdP + (eP - bP)**2 * dl/bMag
            stdD = stdD + (eD - bD)**2 * dl/bMag
            ! Keep track of flux tube volume to normalize by later
        enddo

        dvB = max(dvB, TINY)
        ! Turn summations into standard deviations
        stdP = sqrt(stdP / dvB)
        stdD = sqrt(stdD / dvB)

        end associate

    end subroutine FLStdev
    

    !Flux tube entropy (from BLK)
    function FLEntropy(Model,ebGr,bTrc,GamO) result(S)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        real(rp), optional :: GamO

        real(rp) :: S

        integer :: OCb,k
        real(rp) :: Phi0,dl,bAvg,pAvg,dV,Gam

        if (present(GamO)) then
            Gam = GamO
        else
            Gam = 5.0/3.0
        endif

        Phi0 = 1.0
        S = 0.0

        OCb = FLTop(Model,ebGr,bTrc)
        if ( OCb<2 .or. (.not. Model%doMHD) ) return

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
        do k=-Nm,Np-1
            dl = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            bAvg = 0.5* (bTrc%magB(k+1) + bTrc%magB(k))
            dV = dl*Phi0/bAvg
            pAvg = 0.5*( bTrc%Gas(k+1,PRESSURE,BLK) + bTrc%Gas(k,PRESSURE,BLK) )
            S = S + dV*(pAvg**(1.0/Gam))
        enddo

        end associate

    end function FLEntropy

    function FLTop(Model,ebGr,bTrc) result(OCb)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        integer :: OCb

        logical :: isCP,isCM,isFin,isStart
        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)

        !Test topology
        !OCB =  0 (solar wind), 1 (half-closed), 2 (both ends closed)
        !OCB = -1 (plasmoid/timeout)

        if (.not. bTrc%isGood) then
            OCb = -1
            return
        endif

        isCP  = isClosed(bTrc%xyz(+Np,:),Model)
        isCM  = isClosed(bTrc%xyz(-Nm,:),Model)

        
        !isFin = (Np<MaxFL-1) .and. (Nm<MaxFL-1) !Check if finished
        isFin = (Np<bTrc%nMax-1) .and. (Nm<bTrc%nMax-1) !Check if finished
        isStart = (Np>0) .and. (Nm>0) !Check if both sides went somewhere

        OCb = 0

        if ( isFin ) then
        !Both sides ended normally
            if (isCP .or. isCM) then
                !At least one side is closed
                if (isCP .and. isCM) then
                    OCb = 2
                else
                    OCb = 1
                endif
            else
                !Neither side closed and both sides ended normally
                !IMF
                OCb = 0
            endif
        else
        !At least one side ended badly
            OCb = -1 !Just set timeout flag
            !if ( (Np<MaxFL-1) ) OCb = OCb-1
            !if ( (Nm<MaxFL-1) ) OCb = OCb-1
        endif
        end associate
    end function FLTop

    !Get conjugate point from field line
    subroutine FLConj(Model,ebGr,bTrc,xyzC)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        real(rp), intent(out) :: xyzC(NDIM)

        real(rp), dimension(NDIM) :: xP,xM,x0

        if (.not. bTrc%isGood) then
            xyzC = 0.0
            return
        endif

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
        !Get endpoints of field line
        xP = bTrc%xyz(+Np,:)
        xM = bTrc%xyz(-Nm,:)
        x0 = bTrc%xyz(0  ,:) !Seed point

        !Use seed point to decide which hemisphere we're in
        if (x0(ZDIR) > 0) then
            !We're in NH, find SH conjgate
            xyzC = xM
        else
            !We're in SH, find NH conjugate
            xyzC = xP
        endif

        end associate
    end subroutine FLConj
    
    !Get minimum field strength and location
    subroutine FLEq(Model,bTrc,xeq,Beq)
        type(chmpModel_T), intent(in) :: Model
        type(magLine_T  ), intent(in) :: bTrc
        real(rp), intent(out) :: xeq(NDIM),Beq

        integer :: i0,iMin,OCb
        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)

        if (.not. bTrc%isGood) then
            !Bad field line
            xeq = 0.0
            Beq = 0.0
            return
        endif

        !Find minimum field and where it occurs
        Beq = minval(bTrc%magB) !Min field strength
        i0  = minloc(bTrc%magB,dim=1) !Note this is between 1:N b/c fortran
        !Need to offset i0, iMin = i0-Nm-1
        iMin = i0-Nm-1
        xeq = bTrc%xyz(iMin,:)

        end associate
    end subroutine FLEq

    !Get curvature radius at equator and ExB velocity [km/s]
    subroutine FLCurvRadius(Model,ebGr,ebState,bTrc,rCurv,vEB)
        type(chmpModel_T), intent(in)  :: Model
        type(ebGrid_T)   , intent(in)  :: ebGr
        type(ebState_T)  , intent(in)  :: ebState
        type(magLine_T)  , intent(in)  :: bTrc
        real(rp)         , intent(out) :: rCurv
        real(rp)         , intent(out) :: vEB

        type(gcFields_T) :: gcFields
        real(rp), dimension(NDIM,NDIM) :: Jacbhat
        real(rp), dimension(NDIM) :: xeq,E,B,vExB
        real(rp) :: t,Beq

        rCurv = 0.0
        if (.not. bTrc%isGood) return

        !Get equator
        call FLEq(Model,bTrc,xeq,Beq)

        !Now get field information there
        t = ebState%eb1%time
        call ebFields(xeq,t,Model,ebState,E,B,vExB=vExB,gcFields=gcFields)

        rCurv = getRCurv(B,gcFields%JacB)
        vEB = norm2(vExB)*oVScl
    end subroutine FLCurvRadius

!---------------------------------
!Projection routines
    !Project to SM EQ (Z=0)
    subroutine getEquatorProjection(Model,ebState,x0,t,xe)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in)   :: Model
        type(ebState_T)  , intent(in)   :: ebState
        real(rp), intent(out) :: xe(NDIM) ! end point
        logical :: failEquator

        if (.not. inDomain(x0,Model,ebState%ebGr)) then
           xe = HUGE
           return
        endif
        
        if (x0(ZDIR)>=TINY) then
           ! assume first that northern hemisphere is always traced in -B direction
           call project(Model,ebState,x0,t,xe,-1,.true.,failEquator)
        else if (x0(ZDIR)<=-TINY) then 
           call project(Model,ebState,x0,t,xe,+1,.true.,failEquator)
        else
           ! if we're for some reason at equator don't do anything.
           xe = x0
           return
        endif

        ! at this point, we have either gotten to equator
        ! or to the domain boundary
        ! because we trace along -B for z>0 and along B for z<0
        ! can get to weird places on non-closed field lines or 
        ! even on closed if strongly tilted
        ! trap for those points here 
        if (failEquator) xe = HUGE

      end subroutine getEquatorProjection

    !Project to magnetic equator (min along field line)
    subroutine getMagEQ(Model,ebState,x0,t,xeq,Beq,tOpt)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(out) :: xeq(NDIM),Beq
        integer, intent(out), optional :: tOpt

        type(magLine_T) :: bTrc
        integer :: OCb
        
        !Initialize output
        xeq = 0.0
        Beq = 0.0
        if (present(tOpt)) tOpt = -1

        if (.not. inDomain(x0,Model,ebState%ebGr)) return

        !Trace field line
        call genLine(Model,ebState,x0,t,bTrc)

        !Get diagnostics
        call FLEq(Model,bTrc,xeq,Beq)
        OCb = FLTop(Model,ebState%ebGr,bTrc)

        if (present(tOpt)) then
            tOpt = OCb
        endif

    end subroutine getMagEQ

!---------------------------------
    !Project XYZ to lat-lon on ionosphere
    !TODO: Better merge this with voltron RCM code, too redundant now

    subroutine Map2NH(ebModel,ebState,xyz,t,x1,x2)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebState_T)  , intent(in) :: ebState
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(in) :: t
        real(rp), intent(out) :: x1,x2

        real(rp), dimension(NDIM) :: xE,xIon,xyz0
        real(rp) :: dX,rC,startEps,rEps
        logical :: isGood

        startEps = 0.05
        rEps = 0.125

        x1 = 0.0
        x2 = 0.0

        ! trap for when we're within epsilon of the inner boundary
        ! (really, it's probably only the first shell of nodes at R=Rinner_boundary that doesn't trace correctly)
        if ( (norm2(xyz)-rClosed)/rClosed < startEps ) then
           ! dipole-shift to startEps
           xyz0 = DipoleShift(xyz,norm2(xyz)+startEps)
        else
           xyz0 = xyz
        end if

        !Use one-sided projection routine from chimp
        !Trace along field line (i.e. to northern hemisphere)
        call project(ebModel,ebState,xyz0,t,xE,+1,toEquator=.false.)

        dX = norm2(xyz0-xE)
        rC = rClosed*(1.+rEps)
        isGood = (dX>TINY) .and. (norm2(xE) <= rC) .and. (xE(ZDIR) > 0)

        if (isGood) then
            !Get invariant lat/lon
            x1 = InvLatitude(xE)
            x2 = katan2(xE(YDIR),xE(XDIR))
            if (x2 < 0) x2 = x2 + 2*PI
        else
            !Set 0/0 for projection failure
            x1 = 0.0
            x2 = 0.0
        endif

    end subroutine Map2NH

!---------------------------------
!Tracing routines

    !Calculate one-sided trace (in sgn direction)
    subroutine genTrace(Model,ebState,x0,t,xyzn,ijkn,bxyzn,gasn,Np,sgn,MaxStepsO,doShueO)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in) :: x0(NDIM),t
        real(rp), intent(inout) :: xyzn(0:MaxFL,NDIM),bxyzn(0:MaxFL,NDIM),gasn(0:MaxFL,NVARMHD,0:Model%nSpc)
        integer , intent(inout) :: ijkn(0:MaxFL,NDIM)
        integer , intent(out) :: Np
        integer , intent(in) :: sgn
        integer , intent(in), optional :: MaxStepsO
        logical , intent(in), optional :: doShueO

        integer :: MaxSteps
        type(GridPoint_T) :: gPt
        real(rp) :: h
        real(rp), dimension(NDIM) :: B,oB,dx
        real(rp), dimension(NVARMHD,0:Model%nSpc) :: Q
        logical :: inDom,doShue
        
        if (present(MaxStepsO)) then
            MaxSteps = min(MaxStepsO,MaxFL)
        else
            MaxSteps = MaxFL
        endif
        if (present(doShueO)) then
            doShue = doShueO
        else
            doShue = doShueG
        endif

    !Initialize
        gPt%xyz = x0
        gPt%t   = t
        Np = 0

        xyzn  = 0.0
        gasn  = 0.0
        bxyzn = 0.0
        ijkn  = 0
        
        inDom = inDomain(gPt%xyz,Model,ebState%ebGr)

    !Start main loop
        do while (inDom .and. Np <= MaxSteps)
            if (Np == 0) then
                !First time don't have guess
                call locate(gPt%xyz,gPt%ijkG,Model,ebState%ebGr,inDom)
                gPt%dl = getDiag(ebState%ebGr,gPt%ijkG)
                !Pick first h
                h = sgn*Model%epsds*gPt%dl
                !Get first B field
                B = FastMag(gPt%xyz,gPt%t,Model,ebState,inDom,gPt%ijkG)
            else
                !Locate w/ guess (last known cell)
                !call locate(gPt%xyz,gPt%ijkG,Model,ebState%ebGr,inDom,gPt%ijkG)
                !NOTE: Don't need to re-locate b/c last ijkG is correct from FSAL
                gPt%dl = getDiag(ebState%ebGr,gPt%ijkG)
            endif

        !Get values at this point on the tube
            !Already have B field from FSAL (first same as last) of previous iteration
            !B = fldInterp(gPt%xyz,gPt%t,Model,ebState,BFLD,inDom,gPt%ijkG)
            xyzn (Np,:) = gPt%xyz
            ijkn (Np,:) = gPt%ijkG
            bxyzn(Np,:) = B

            if (Model%doMHD) then
                Q = mhdInterpMF(gPt%xyz,gPt%t,Model,ebState,gPt%ijkG)
                gasn(Np,1:NVARMHD,0:Model%nSpc) = Q
            endif

        !Now do step
            !Call streamstep providing first B and getting last B
            call StreamStep(gPt,Model,ebState,Model%epsds,h,dx,B,oB)
            B = oB !Setup B field for next iteration
            gPt%xyz = gPt%xyz + dx
            !Verify sign of h to be sure
            h = sign(h,sgn*1.0_rp)

        !Check if we're done
            inDom = inDomain(gPt%xyz,Model,ebState%ebGr) .and. (norm2(dx)>TINY)
            if ( doShue .and. (modulo(Np,NpChk) == 0) .and. (norm2(gPt%xyz)>=rShue) ) then
                inDom = inDom .and. inShueMP_SM(gPt%xyz,ShueScl)
            endif

            if (inDom) Np = Np + 1
        enddo

        ! check if exceeded tube bounds
        if(Np > MaxSteps) then
            Np = MaxSteps
            !Removing warning that is triggered too often and probably not beneficial
            !!$OMP CRITICAL
            !write(*,*) ANSIRED
            !write(*,*) "<WARNING! genTrace hit max tube size!>"
            !write(*,*) "Seed: ", x0
            !write(*,*) "End : ", xyzn(Np,:)
            !write(*,'(a)',advance="no") ANSIRESET, ''
            !!$OMP END CRITICAL
        endif                            
    end subroutine genTrace

    !Slimmed down projection to northern hemisphere for MAGE
    !RinO is optional cut-off inner radius when in northern hemisphere
    !epsO is optional epsilon (otherwise use Model default)
    subroutine mageproject(Model,ebState,x0,t,xyz,Np,isG,epsO,MaxStepsO,doShueO)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in)  :: x0(NDIM),t
        real(rp), intent(out) :: xyz(NDIM)
        integer , intent(out) :: Np
        logical,  intent(out) :: isG
        real(rp), intent(in), optional :: epsO
        integer , intent(in), optional :: MaxStepsO
        logical , intent(in), optional :: doShueO

        type(GridPoint_T) :: gPt
        integer :: sgn,MaxSteps
        real(rp) :: eps,h
        real(rp), dimension(NDIM) :: dx,B,oB
        logical :: inDom,isSC,isDone,doShue

        if (present(epsO)) then
            eps = epsO
        else
            eps = Model%epsds
        endif

        if (present(MaxStepsO)) then
            MaxSteps = MaxStepsO
        else
            MaxSteps = MaxFL
        endif

        if (present(doShueO)) then
            doShue = doShueO
        else
            doShue = doShueG
        endif

        sgn = +1 !Step towards NH
        isG = .false.
    !Initialize
        gPt%xyz = x0
        gPt%t   = t
        Np = 0
        call CheckDone(Model,ebState,gPt%xyz,inDom,isSC,isDone,Np)

        do while ( (.not. isDone) .and. (Np <= MaxSteps) )
        !Locate and get fields
            !Get location in ijk using old ijk as guess if possible
            if (Np == 0) then
                !First step, no guess yet
                call locate(gPt%xyz,gPt%ijkG,Model,ebState%ebGr,inDom)
                gPt%dl = getDiag(ebState%ebGr,gPt%ijkG)
                !Pick first h
                h = sgn*eps*gPt%dl
                !Get first B field
                B = FastMag(gPt%xyz,gPt%t,Model,ebState,inDom,gPt%ijkG)
            else
                !Otherwise don't need to locate b/c last ijkG is correct
                gPt%dl = getDiag(ebState%ebGr,gPt%ijkG)                
            endif

        !Now do step
            call StreamStep(gPt,Model,ebState,eps,h,dx,B,oB)
            B = oB !Setup B field for next iteration
            gPt%xyz = gPt%xyz + dx
            !Verify sign of h to be sure
            h = sign(h,sgn*1.0_rp)

        !Check if we're done
            if (norm2(dx)<TINY) then
                !Multistep integrator had one out value
                isDone = .true.
                inDom  = .false.
                isSC   = .false.
            else
                call CheckDone(Model,ebState,gPt%xyz,inDom,isSC,isDone,Np)
            endif
            Np = Np + 1

        enddo !Main stepping loop

        xyz = gPt%xyz
    !Finished loop somehow, decide if it was good
        
        if (isSC .or. (Np>=MaxSteps)) then
            !Got short circuited or timed out, not good
            isG = .false.
            return
        endif

        if (inDom) then
            !This shouldn't happen
            !$OMP CRITICAL
            write(*,*) 'How did you get here? Bad thing in mageproject'
            !$OMP END CRITICAL
            stop
        else
            !Left domain, but not sure if left from inner boundary
            if (isClosed(xyz,Model)) then
                isG = .true.
            else
                isG = .false.
            endif
            return
        endif

        contains
             
            subroutine CheckDone(Model,ebState,xyz,inDom,isSC,isDone,Np)
                type(chmpModel_T), intent(in) :: Model
                type(ebState_T), intent(in)   :: ebState
                real(rp), intent(inout) :: xyz(NDIM)
                logical, intent(out) :: inDom,isSC,isDone
                integer, intent(in) :: Np
                logical :: inMP

                inDom = inDomain(xyz,Model,ebState%ebGr)
                if ( doShue .and. (modulo(Np,NpChk) == 0) .and. (norm2(gPt%xyz)>=rShue) ) then
                    inMP  = inShueMP_SM(xyz,ShueScl)
                else
                    inMP = .true.
                endif
                
                if (inDom) then
                    !In domain, check if in Shue MP w/ safety factor
                    if (inMP) then
                        isSC   = .false.
                        isDone = .false.
                    else
                        !Not in shue mp
                        isSC = .true.
                        isDone = .true.
                    endif
                else
                    !Straight up out of domain
                    isSC = .false.
                    isDone = .true.
                endif

            end subroutine CheckDone

    end subroutine mageproject

!====      
    !Calculate one-sided projection (in sgn direction)
    subroutine project(Model,ebState,x0,t,xn,sgn,toEquator,failEquator)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in) :: x0(NDIM),t
        integer :: Np
        integer, intent(in) :: sgn
        real(rp), intent(inout) :: xn(NDIM)
        logical, optional, intent(in) :: toEquator
        logical, optional, intent(out) :: failEquator

        logical :: inDom
        real(rp), dimension(NDIM) :: B,E,dx
        real(rp), dimension(NDIM) :: Jb,Jb2,Jb3,F1,F2,F3,F4
        real(rp), dimension(NDIM,NDIM) :: JacB
        real(rp) :: ds,dsmag,dl,MagB,MagJb,dzSgn
        real(rp), dimension(NVARMHD) :: Q
        integer, dimension(NDIM) :: ijk,ijkG
        type(gcFields_T) :: gcF

        !Prime pump
        call locate(x0,ijk,Model,ebState%ebGr,inDom)
        B = fldInterp(x0,t,Model,ebState,BFLD,inDom,ijk)
        
        !Prep for iteration
        Np = 0
        Xn = x0
        dl = getDiag(ebState%ebGr,ijk)
        ds = sgn*Model%epsds*dl

        ijkG = ijk

        if (present(toEquator)) then
           if ((toEquator).and.(.not.(present(failEquator)))) then
              write(*,*) 'Project operator called incorrectly.'
              stop
           endif
        end if

        if (present(failEquator)) failEquator = .true.

        !write(*,*) 'sgn/ds/X0 = ', sgn,ds,x0
        do while (inDom .and. Np <= MaxFL)
        !Locate and get fields
            !Get location in ijk using old ijk as guess
            call locate(Xn,ijk,Model,ebState%ebGr,inDom,ijkG)
            call ebFields(Xn,t,Model,ebState,E,B,ijk,gcFields=gcF)
            MagB = norm2(B)

            ! get the jacobian and new ds
            JacB = gcF%JacB
            MagJb = norm2(JacB)
            dl = getDiag(ebState%ebGr,ijk)

            if (MagJb <= TINY) then
                !Field is constant-ish, use local grid size
                dsmag = dl
            else
                dsmag = MagB/MagJb
            endif
            ds = sgn*Model%epsds*min(dl,dsmag)
        !Update position
            !Convert ds to streamline units
            ds = ds/max(MagB,TINY)
            
            !Get powers of jacobian
            Jb  = matmul(JacB,B  )
            Jb2 = matmul(JacB,Jb )
            Jb3 = matmul(JacB,Jb2)

            !Calculate steps
            F1 = ds*B
            F2 = F1 + (ds*ds/2)*Jb
            F3 = F2 + (ds*ds*ds/4)*Jb2
            F4 = ds*B + ds*ds*Jb + (ds**3.0/2.0)*Jb2 + (ds**4.0/4.0)*Jb3

            !Advance
            dx = (F1+2*F2+2*F3+F4)/6.0
            Xn = Xn + dx
        !Prep for next step, test exit criteria
            inDom = inDomain(xn,Model,ebState%ebGr)
            if (inDom) then
                Np = Np+1

                if (toEquator) then
                    !Check for x'ing

                   ! (trap for a rare weird case)
                   ! if we happened to be just at the equator before making the step, don't do anything
                   ! this can happen if we start tracing exactly from z=0
                   ! the calling getEquatorProjection function should capture this case
                   ! but still add this trap here, just in case
                    dzSgn = Xn(ZDIR)*( Xn(ZDIR)-dx(ZDIR) )
                    if ((dzSgn < 0).or.(Xn(ZDIR)-dx(ZDIR).eq.0.)) then
                        ! interpolate exactly to equator
                        Xn = Xn-dx*abs(Xn(ZDIR))/abs(dx(ZDIR))

                        failEquator = .false.
                        return
                    endif !dzSgn
                 end if
            endif !inDom
 
        enddo
        !write(*,*) 'Found sign/points/distance = ', sgn,Np,norm2(x0-Xn)

    end subroutine project
    
    !Simple routine to print tube information to console
    subroutine PrintTube(ebModel,ebGr,bTrc)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc

        integer :: k
        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
        write(*,*) "---"
        write(*,*) "Bad Tube! Nm/Np = ", Nm,Np
        do k=-Nm,Np
            write(*,*) "k / xyz = ", k,bTrc%xyz (k,:)
        enddo
        end associate
    end subroutine PrintTube

end module streamline
