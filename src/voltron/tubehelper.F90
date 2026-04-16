
!Helper routines to handle tube tracing and data extraction
module tubehelper
    use volttypes
    use voltcpltypes
    use streamline
    use chmpdbz
    use shellUtils
    use shellGridIO
    use msphutils, only : RadIonosphere
    use raijudefs, only : def_tiote

    implicit none

    contains

    subroutine init_TubeShell(sh, tubeShell, TioTeO)
        type(ShellGrid_T), intent(in) :: sh
            !! Voltron shellgrid
        type(TubeShell_T), intent(inout) :: tubeShell
            !! IMAGTubeShell object we are initializing
        !logical, dimension(:,:), intent(in), optional :: maskO
        real(rp), intent(in), optional :: TioteO
            !! Default Ti/Te ratio

        integer :: i
        real(rp) :: tiote

        if (present(TioteO)) then
            tiote = TioteO
        else 
            tiote = def_tiote
        endif


        do i=1,NDIM
            call initShellVar(sh, SHGR_CORNER, tubeShell%xyz0(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%X_bmin(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%X_maxRad(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%Jxyz_iono(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%Bxyz_iono(i))
        enddo
        call initShellVar(sh, SHGR_CORNER, tubeShell%lat0   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%lon0   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%latc   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%lonc   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%J11i   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%invlat )
        call initShellVar(sh, SHGR_CORNER, tubeShell%topo   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%bmin   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%bVol   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%Lb     )
        call initShellVar(sh, SHGR_CORNER, tubeShell%Tb     )
        call initShellVar(sh, SHGR_CORNER, tubeShell%wMAG   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%rCurv  )
        call initShellVar(sh, SHGR_CORNER, tubeShell%avgBeta)
        do i=0,MAXTUBEFLUIDS
            call initShellVar(sh, SHGR_CORNER, tubeShell%avgP(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%avgN(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%stdP(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%stdN(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%bcN (i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%bcP (i))
        enddo
        call initShellVar(sh, SHGR_CORNER, tubeShell%losscone )
        call initShellVar(sh, SHGR_CORNER, tubeShell%lossconec)
        call initShellVar(sh, SHGR_CORNER, tubeShell%TioTe0   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%nTrc     )

        tubeShell%TioTe0%data = tiote


    end subroutine init_TubeShell

    subroutine FakeTube(P,xyz0,bTube)
        !! Given a seed xyz0 generate a fake tube that seems plausible
        type(planet_T), intent(in)  :: P
        real(rp)      , intent(in)  :: xyz0(NDIM)
        type(Tube_T)  , intent(out) :: bTube

        call FreshenTube(bTube)
        call DipoleTube(P,xyz0,bTube)
        !Add fake plasma stuff
        !TODO: Add calls to gallagher/TM03/etc

    end subroutine FakeTube

    
    subroutine DipoleTube(P,xyz0,bTube)
        !! Given a seed xyz0 generate an empty tube w/ dipole quantities
        type(planet_T), intent(in)  :: P
        real(rp)      , intent(in)  :: xyz0(NDIM)
        type(Tube_T)  , intent(out) :: bTube

        real(rp) :: L,bIon

        call FreshenTube(bTube)
        !Start w/ recording seed information and such
        bTube%xyz0(:) = xyz0
        call xyz2LL(bTube%xyz0,bTube%lat0,bTube%lon0)
        bTube%invlat = InvLatitude(bTube%xyz0)

        !Topology stuff
        bTube%topo = TUBE_CLOSED !Always closed
        bTube%latc = -bTube%lat0
        bTube%lonc =  bTube%lon0
        
        !Magnetic stuff
        L = DipoleL(xyz0)
        bTube%bmin = abs(P%magMoment)*G2T*1.0e+9/L**3.0 !Min B in nT
        bTube%X_bmin = L*[cos(bTube%lon0),sin(bTube%lon0),0.0_rp]
        bTube%X_maxRad = bTube%X_bmin
        bTube%bVol = DipFTV_L(L,P%magMoment) !Rx/nT
        bTube%Lb = L !Just lazily using L shell
        bTube%rCurv = L/3.0
        bTube%wMAG = 1.0
        
        !Loss cone stuff
        bIon = norm2(DipoleB0(xyz0))*oBScl !Ionospheric field strength [nT]
        bTube%losscone = asin(sqrt(bTube%bmin/bIon))
        bTube%lossconec = bTube%losscone !Assuming symmetry

        !Null plasma stuff
        bTube%Tb = 0.0
        bTube%avgBeta = 0.0
        bTube%avgP(:) = 0.0
        bTube%stdP(:) = 0.0
        bTube%avgN(:) = 0.0
        bTube%stdN(:) = 0.0

        bTube%TioTe0 = 1.0 !Meh
        
    end subroutine DipoleTube

    subroutine Line2Tube(ebTrcApp,P,bTrc,bTube)
        !! Given a traced field line, populate tube object
        type(ebTrcApp_T), intent(in) :: ebTrcApp
        type(planet_T)  , intent(in)  :: P
        type(magLine_T) , intent(in)  :: bTrc
        type(Tube_T)    , intent(out) :: bTube

        integer , parameter :: i0=2
        real(rp), parameter :: L_ibc = 2.0 !Length to integrate along field line for ibc value
        integer :: OCb, s
        real(rp) :: bD,bP,bD_rc,bP_rc,dvB,bBeta,stdD,stdP,bmin
        real(rp) :: bIon,bIonC,rCurv,VebMKS,VaMKS,TiEV,CsMKS

        real(rp), dimension(NDIM) :: bEq,xyzC,xyzIonC

        call FreshenTube(bTube)

        !Start w/ recording seed information and such
        bTube%xyz0(:) = bTrc%xyz(0,:)

        call xyz2LL(bTube%xyz0,bTube%lat0,bTube%lon0)

        bTube%invlat = InvLatitude(bTube%xyz0)

        associate(ebModel=>ebTrcApp%ebModel,ebGr=>ebTrcApp%ebState%ebGr,ebState=>ebTrcApp%ebState)

        !Record number of points
        bTube%nTrc = bTrc%Np + bTrc%Nm !Line array is -Nm:+Np w/ 0=seed
        bTube%J11i = 0.0

        !Figure out topology
        bTube%topo = VoltTopo(ebModel,ebGr,bTrc)
        if (bTube%topo == TUBE_DISASTER) then
            !$OMP CRITICAL
            call PrintTube(ebModel,ebGr,bTrc)
            !$OMP END CRITICAL
            stop
        endif

    !If we're still here, let's do some stuff
        !Start w/ stuff for both closed and open
        bTube%TioTe0 = def_tiote !For now

        !Get J11 at radius of shell i0
        bTube%J11i = FLIonFAC(ebModel,ebGr,ebState,bTrc,RadIonosphere(),i0)
        call FLIonBJ(ebModel,ebGr,ebState,bTrc,RadIonosphere(),i0,bTube%Bxyz_iono,bTube%Jxyz_iono)

        if (bTube%topo == TUBE_CLOSED) then
            do s=0,ebModel%nSpc
                call FLThermo(ebModel,ebGr,bTrc,bD,bP,dvB,bBeta,s)
                call FLStdev (ebModel,ebGr,bTrc,bD,bP,stdD,stdP,s)
                bTube%avgP(s) = bP
                bTube%stdP(s) = stdP
                bTube%avgN(s) = bD
                bTube%stdN(s) = stdD
                if (s == BLK) then
                    bTube%avgBeta = bBeta
                endif
                !call FL2iBC(ebModel,ebGr,bTrc,bD,bP,i0,L_ibc,s)
                call FLIonDP(ebModel,ebGr,ebState,bTrc,bD,bP,i0,s)
                bTube%bcN (s) = bD
                bTube%bcP (s) = bP
            enddo !spc
        else
            !Open or undefined
            bTube%avgBeta = 0.0
            do s=0,ebModel%nSpc
                call FLIonDP(ebModel,ebGr,ebState,bTrc,bD,bP,i0,s)
                bTube%avgN(s) = bD
                bTube%avgP(s) = bP
                bTube%stdN(s) = 0.0
                bTube%stdP(s) = 0.0

                !call FL2iBC(ebModel,ebGr,bTrc,bD,bP,i0,L_ibc,s)
                bTube%bcN (s) = bD
                bTube%bcP (s) = bP
            enddo !spc
        endif

        if (bTube%topo == TUBE_CLOSED) then
        !Do closed field line stuff
            bTube%bVol = dvB/oBScl !Rx/EB => Rx/nT

        !Minimal surface (bEq in Rp, bMin in EB)
            call FLEq(ebModel,bTrc,bEq,bmin)
            bmin = bmin*oBScl !EB => nT
            bTube%bmin = bmin
            bTube%X_bmin = bEq

            bTube%X_maxRad = FLMaxRadX(bTrc)
            
        !Line stuff
            bTube%Lb = FLArc(ebModel,ebGr,bTrc) !Rx
            !TODO: Check units of AlfvenX calculation
            bTube%Tb = FLAlfvenX(ebModel,ebGr,bTrc)
        !Conjugate stuff
            call FLConj(ebModel,ebGr,bTrc,xyzC)
            !Shift conjugate point to match radius of seed point
            xyzIonC = DipoleShift(xyzC, norm2(bTube%xyz0))

            call xyz2LL(xyzIonC,bTube%latc,bTube%lonc)
        !Loss cone stuff
            bIon  = norm2(DipoleB0(bTube%xyz0))*oBScl !Ionospheric field strength [nT]
            bIonC = norm2(DipoleB0(xyzIonC   ))*oBScl !Ionospheric field strength [nT]
            bTube%losscone   = asin(sqrt(bmin/bIon))
            bTube%lossconec  = asin(sqrt(bmin/bIonC)) !Keeping these separate for IGRF later

        !Get energy partition and curvature stuff
            !Get curvature radius and ExB velocity [km/s]
            call FLCurvRadius(ebModel,ebGr,ebState,bTrc,rCurv,VebMKS)
            bTube%rCurv = rCurv
            bD = bTube%avgN(BLK)
            bP = bTube%avgP(BLK)
            bD_rc = bTube%avgN(RCFLUID)
            bP_rc = bTube%avgP(RCFLUID)

            !VaMKS = flux tube arc length [km] / Alfven crossing time [s]
            VaMKS = (bTube%Lb*P%rp_m*1.0e-3)/bTube%Tb 
            !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
            TiEV = (1.0e+3)*DP2kT(bD,bP) !Temp in eV
            !TiEV = (1.0e+3)*DP2kT(bD_rc,bP_rc) !Temp in eV
            CsMKS = 9.79*sqrt((5.0/3)*TiEV)
            bTube%wMAG = VaMKS/( sqrt(VaMKS**2.0 + CsMKS**2.0) + VebMKS)
        endif

        end associate !ebTrc stuff

    end subroutine Line2Tube

    function FLMaxRadX(bTrc) result(xMaxRad)
        type(magLine_T), intent(in) :: bTrc
        real(rp) :: xMaxRad(NDIM)

        integer :: k
        real(rp) :: rNow, rMax

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
            xMaxRad = bTrc%xyz(0,:)
            rMax = norm2(xMaxRad)

            do k=-Nm,Np
                rNow = norm2(bTrc%xyz(k,:))
                if (rNow > rMax) then
                    rMax = rNow
                    xMaxRad = bTrc%xyz(k,:)
                endif
            enddo
        end associate

    end function FLMaxRadX

    !NOTE: This function should only be called for field lines traced from the Voltron grid
    !If you're calling this in any other situation, you'll probably get what you deserve
    function VoltTopo(ebModel,ebGr,bTrc) result(topo)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        integer :: topo

        integer, parameter :: i0=3 !Closed field line criteria
        integer :: iP,iM
        real(rp), dimension(NDIM) :: xyzP,xyzM
        real(rp) :: R0,zz
        logical :: isIonP,isIonM

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
            !Get end points
            xyzP = bTrc%xyz(+Np,:)
            xyzM = bTrc%xyz(-Nm,:)
            zz = xyzP(ZDIR)*xyzM(ZDIR) !<0 for cross-hemisphere

            !Get i cell of last point for each side
            iP = bTrc%ijk(+Np,IDIR)
            iM = bTrc%ijk(-Nm,IDIR)

            !Assuming one end point is attached to voltron grid
            isIonP = .false. !xyzP is @ ion
            isIonM = .false. !xyzM is @ ion

            if (Np == 0) isIonP = .true.
            if (Nm == 0) isIonM = .true.

            if ( (.not. isIonP) .and. (.not. isIonM) ) then
                !What the shit?!
                topo = TUBE_DISASTER
                return
            endif
            R0 = norm2(ebGr%xyzcc(i0,ebGr%js,ebGr%ks,:))
            if (isIonP) then
                if ((norm2(xyzM) <= R0) .and. (zz <= 0)) then
                    !Both sides are closed and cross equator
                    topo = TUBE_CLOSED
                elseif (isOut(ebModel,ebGr,iM)) then
                    topo = TUBE_OPEN
                else
                    topo = TUBE_UNDEF
                endif
            elseif (isIonM) then
                !Same logic as above
                if ((norm2(xyzP) <= R0) .and. (zz <= 0)) then
                    topo = TUBE_CLOSED
                elseif (isOut(ebModel,ebGr,iP)) then
                    topo = TUBE_OPEN
                else
                    topo = TUBE_UNDEF
                endif
            endif

        end associate

        contains
            function isOut(ebModel,ebGr,iX)
                type(chmpModel_T), intent(in) :: ebModel
                type(ebGrid_T)   , intent(in) :: ebGr
                integer          , intent(in) :: iX

                logical :: isOut
                if ( abs(iX-ebGr%ie) <= 4+1 ) then
                    !This cell is pretty close to the outer boundary
                    isOut = .true.
                else
                    !Meh, could be out-er
                    isOut = .false.
                endif
            end function isOut

    end function VoltTopo


    subroutine xyz2LL(xyz,lat,lon)
        real(rp), intent(in)  :: xyz(NDIM)
        real(rp), intent(out) :: lat,lon

        lon = katan2(xyz(YDIR),xyz(XDIR) )
        lat = asin  (xyz(ZDIR)/norm2(xyz))

    end subroutine xyz2LL

    !Just zero out values
    subroutine FreshenTube(bTube)
        type(Tube_T), intent(out) :: bTube
        bTube%xyz0    = 0.0
        bTube%lat0    = 0.0
        bTube%lon0    = 0.0
        bTube%invlat  = 0.0
        bTube%latc    = 0.0
        bTube%lonc    = 0.0
        
        bTube%topo = TUBE_UNDEF
        bTube%bmin = 0.0
        bTube%X_bmin = 0.0
        bTube%X_maxRad = 0.0
        bTube%Jxyz_iono = 0.0
        bTube%Bxyz_iono = 0.0
        bTube%bVol = 0.0
        bTube%Lb = 0.0
        bTube%Tb = 0.0
        bTube%TioTe0 = def_tiote

        bTube%wMAG = 0.0
        bTube%rCurv = 0.0
        bTube%avgBeta = 0.0

        bTube%avgP(:) = 0.0
        bTube%stdP(:) = 0.0
        bTube%avgN(:) = 0.0
        bTube%stdN(:) = 0.0

        bTube%bcN(:)  = 0.0
        bTube%bcP(:)  = 0.0

        bTube%losscone  = 0.0
        bTube%lossconec = 0.0

        bTube%TioTe0 = 0.0

        bTube%nTrc = 0

        !Zero pot
        bTube%pot   = 0.0
        bTube%crpot = 0.0
        bTube%potc  = 0.0

    end subroutine FreshenTube


    subroutine tubes2Shell(shGr, ijTubes, tubeShell)
        type(ShellGrid_T), intent(in) :: shGr
            !! ShellGrid that ijTubes and tubeShell live on
        type(Tube_T), dimension(shGr%is:shGr%ie+1,shGr%js:shGr%je+1), intent(in) :: ijTubes
            !! Nt+1,Np+1 array of Tube_T's
        type(TubeShell_T), intent(inout) :: tubeShell
            !! ShellGridVar versions of 2D Tube_T data we are populating

        integer :: i,j,k
        logical, dimension(shGr%isg:shGr%ieg+1,shGr%jsg:shGr%jeg+1) :: topoMask

        ! Copy tubes to active tubeShell domain
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k)
        do j=shGr%js,shGr%je+1
            do i=shGr%is,shGr%ie+1
                do k=1,NDIM
                    tubeShell%xyz0  (k)%data(i,j) = ijTubes(i,j)%xyz0(k)
                    tubeShell%X_bmin(k)%data(i,j) = ijTubes(i,j)%X_bmin(k)
                    tubeShell%X_maxRad(k)%data(i,j) = ijTubes(i,j)%X_maxRad(k)
                    tubeShell%Jxyz_iono(k)%data(i,j) = ijTubes(i,j)%Jxyz_iono(k)
                    tubeShell%Bxyz_iono(k)%data(i,j) = ijTubes(i,j)%Bxyz_iono(k)
                enddo
                tubeShell%lat0   %data(i,j) = ijTubes(i,j)%lat0
                tubeShell%lon0   %data(i,j) = ijTubes(i,j)%lon0
                tubeShell%invlat %data(i,j) = ijTubes(i,j)%invlat
                tubeShell%latc   %data(i,j) = ijTubes(i,j)%latc
                tubeShell%lonc   %data(i,j) = ijTubes(i,j)%lonc
                tubeShell%J11i   %data(i,j) = ijTubes(i,j)%J11i
                tubeShell%topo   %data(i,j) = ijTubes(i,j)%topo
                tubeShell%bmin   %data(i,j) = ijTubes(i,j)%bmin
                tubeShell%bVol   %data(i,j) = ijTubes(i,j)%bVol
                tubeShell%Lb     %data(i,j) = ijTubes(i,j)%Lb
                tubeShell%Tb     %data(i,j) = ijTubes(i,j)%Tb
                tubeShell%wMAG   %data(i,j) = ijTubes(i,j)%wMAG
                tubeShell%rCurv  %data(i,j) = ijTubes(i,j)%rCurv
                tubeShell%avgBeta%data(i,j) = ijTubes(i,j)%avgBeta
                do k=0,MAXTUBEFLUIDS
                    tubeShell%avgP(k)%data(i,j) = ijTubes(i,j)%avgP(k)
                    tubeShell%avgN(k)%data(i,j) = ijTubes(i,j)%avgN(k)
                    tubeShell%stdP(k)%data(i,j) = ijTubes(i,j)%stdP(k)
                    tubeShell%stdN(k)%data(i,j) = ijTubes(i,j)%stdN(k)
                    tubeShell%bcN (k)%data(i,j) = ijTubes(i,j)%bcN (k)
                    tubeShell%bcP (k)%data(i,j) = ijTubes(i,j)%bcP (k)
                enddo
                tubeShell%losscone %data(i,j) = ijTubes(i,j)%losscone
                tubeShell%lossconec%data(i,j) = ijTubes(i,j)%lossconec
                tubeShell%TioTe0   %data(i,j) = ijTubes(i,j)%TioTe0
                tubeShell%nTrc     %data(i,j) = ijTubes(i,j)%nTrc
            enddo
        enddo

        ! Now wrap everyone so ghosts are filled in
        do k=1,NDIM
            call wrapJ_SGV(shGr, tubeShell%xyz0(k))
            call wrapJ_SGV(shGr, tubeShell%X_bmin(k))
            call wrapJ_SGV(shGr, tubeShell%X_maxRad(k))
            call wrapJ_SGV(shGr, tubeShell%Jxyz_iono(k))
            call wrapJ_SGV(shGr, tubeShell%Bxyz_iono(k))
        enddo
        call wrapJ_SGV(shGr, tubeShell%lat0   )
        call wrapJ_SGV(shGr, tubeShell%lon0   )
        call wrapJ_SGV(shGr, tubeShell%invlat )
        call wrapJ_SGV(shGr, tubeShell%latc   )
        call wrapJ_SGV(shGr, tubeShell%lonc   )
        call wrapJ_SGV(shGr, tubeShell%J11i   )
        call wrapJ_SGV(shGr, tubeShell%topo   )
        call wrapJ_SGV(shGr, tubeShell%bmin   )
        call wrapJ_SGV(shGr, tubeShell%bVol   )
        call wrapJ_SGV(shGr, tubeShell%Lb     )
        call wrapJ_SGV(shGr, tubeShell%Tb     )
        call wrapJ_SGV(shGr, tubeShell%wMAG   )
        call wrapJ_SGV(shGr, tubeShell%rCurv  )
        call wrapJ_SGV(shGr, tubeShell%avgBeta)
        do k=0,MAXTUBEFLUIDS
            call wrapJ_SGV(shGr, tubeShell%avgP(k))
            call wrapJ_SGV(shGr, tubeShell%avgN(k))
            call wrapJ_SGV(shGr, tubeShell%stdP(k))
            call wrapJ_SGV(shGr, tubeShell%stdN(k))
            call wrapJ_SGV(shGr, tubeShell%bcN (k))
            call wrapJ_SGV(shGr, tubeShell%bcP (k))
        enddo
        call wrapJ_SGV(shGr, tubeShell%losscone )
        call wrapJ_SGV(shGr, tubeShell%lossconec)
        call wrapJ_SGV(shGr, tubeShell%TioTe0   )
        call wrapJ_SGV(shGr, tubeShell%nTrc     )


        ! Then apply topo mask to vars that only make sense for closed tubes
        where (abs(tubeShell%topo%data - TUBE_CLOSED) < TINY)
            topoMask = .true.
        elsewhere
            topoMask = .false.
        end where

        do k=1,NDIM
            ! xyz0 is good
            tubeShell%X_bmin(k)%mask(:,:) = topoMask
            tubeShell%X_maxRad(k)%mask(:,:) = topoMask
        enddo
        ! lat0,lon0,invlat are good
        tubeShell%latc   %mask(:,:) = topoMask
        tubeShell%lonc   %mask(:,:) = topoMask
        tubeShell%bmin   %mask(:,:) = topoMask
        tubeShell%bVol   %mask(:,:) = topoMask
        tubeShell%Lb     %mask(:,:) = topoMask
        tubeShell%Tb     %mask(:,:) = topoMask
        tubeShell%wMAG   %mask(:,:) = topoMask
        tubeShell%rCurv  %mask(:,:) = topoMask
        tubeShell%avgBeta%mask(:,:) = topoMask
        do k=0,MAXTUBEFLUIDS
            tubeShell%avgP(k)%mask(:,:) = topoMask
            tubeShell%avgN(k)%mask(:,:) = topoMask       
            tubeShell%stdP(k)%mask(:,:) = topoMask
            tubeShell%stdN(k)%mask(:,:) = topoMask
        enddo
        tubeShell%TioTe0%mask(:,:) = topoMask
        ! idk about losscone, lossconec
        ! nTrc is good
    end subroutine tubes2Shell

    subroutine writeTubeShell(shGr, tubeShell, outF, gStrO, doGhostsO,doFatO)
        type(ShellGrid_T), intent(in) :: shGr
        type(TubeShell_T), intent(in) :: tubeShell
        character(len=strLen), intent(in) :: outF
        character(len=strLen), intent(in), optional :: gStrO
        !! Group name to use instead of 'TubeShell' at root level
        logical, intent(in), optional :: doGhostsO
        logical, intent(in), optional :: doFatO

        integer :: is,ie,js,je
        character(len=strLen) :: gStrOut
        logical :: doWriteMask, isStepOutput, doGhosts,doFat
        integer, dimension(4) :: outSGVBnds_corner
        type(IOVAR_T), dimension(80) :: IOVars

        if(present(gStrO)) then
            gStrOut = trim(gStrO)
        else
            gStrOut = "/TubeShell"
        endif

        isStepOutput = (gStrOut(1:5) == "Step#")

        if (isStepOutput) then
            doWriteMask = .false.  ! Don't add all the masks to every Step, thats a lot of unnecessary data
        else
            doWriteMask = .true.
        endif

        if (isStepOutput .and. present(doGhostsO)) then
            doGhosts = doGhostsO
        else
            doGhosts = .true.  ! Default to saving ghosts
        endif
        
        if (doGhosts) then
            is = shGr%isg
            ie = shGr%ieg
            js = shGr%jsg
            je = shGr%jeg
        else
            is = shGr%is
            ie = shGr%ie
            js = shGr%js
            je = shGr%je
        endif

        if (present(doFatO)) then
            doFat = doFatO
        else
            doFat = .true.
        endif
        outSGVBnds_corner = (/is,ie+1,js,je+1/)

        call AddOutSGV(IOVars,"X_bmin"   ,tubeShell%X_bmin   , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="Rp")
        call AddOutSGV(IOVars,"X_maxRad" ,tubeShell%X_maxRad , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="Rp")
        call AddOutSGV(IOVars,"latc"     ,tubeShell%latc     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="deg")
        call AddOutSGV(IOVars,"lonc"     ,tubeShell%lonc     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="deg")
        call AddOutSGV(IOVars,"J11i"     ,tubeShell%J11i     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="muA/m^2")
        
        call AddOutSGV(IOVars,"topo"     ,tubeShell%topo     , outBndsO=outSGVBnds_corner, doWriteMaskO=.true.     , uStr="enum")
        call AddOutSGV(IOVars,"bmin"     ,tubeShell%bmin     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="nT")
        call AddOutSGV(IOVars,"Jxyz_iono",tubeShell%Jxyz_iono, outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="muA/m^2")
        call AddOutSGV(IOVars,"Bxyz_iono",tubeShell%Bxyz_iono, outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="nT")
        call AddOutSGV(IOVars,"bVol"     ,tubeShell%bVol     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="Rp/nT")
        call AddOutSGV(IOVars,"Lb"       ,tubeShell%Lb       , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="Rp")
        call AddOutSGV(IOVars,"Tb"       ,tubeShell%Tb       , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="s")
        call AddOutSGV(IOVars,"wMAG"     ,tubeShell%wMAG     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="normalized")
        call AddOutSGV(IOVars,"rCurv"    ,tubeShell%rCurv    , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="Rp")
        call AddOutSGV(IOVars,"avgBeta"  ,tubeShell%avgBeta  , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="normalized")
        call AddOutSGV(IOVars,"avgP"     ,tubeShell%avgP     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="nPa")
        call AddOutSGV(IOVars,"avgN"     ,tubeShell%avgN     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="#/cc")
        call AddOutSGV(IOVars,"stdP"     ,tubeShell%stdP     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="nPa")
        call AddOutSGV(IOVars,"stdN"     ,tubeShell%stdN     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="#/cc")
        call AddOutSGV(IOVars, "bcP"     ,tubeShell%bcP      , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="nPa")
        call AddOutSGV(IOVars, "bcN"     ,tubeShell%bcN      , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="#/cc")
        
        call AddOutSGV(IOVars,"losscone" ,tubeShell%losscone , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="rad")
        call AddOutSGV(IOVars,"lossconec",tubeShell%lossconec, outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="rad")
        call AddOutSGV(IOVars,"Tiote0"   ,tubeShell%TioTe0   , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="normalized")
        call AddOutSGV(IOVars,"nTrc"     ,tubeShell%nTrc     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="counts")

        if (doFat) then
            call AddOutSGV(IOVars,"xyz0"     ,tubeShell%xyz0     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="Rp")
            call AddOutSGV(IOVars,"lat0"     ,tubeShell%lat0     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="deg")
            call AddOutSGV(IOVars,"lon0"     ,tubeShell%lon0     , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="deg")
            call AddOutSGV(IOVars,"invlat"   ,tubeShell%invlat   , outBndsO=outSGVBnds_corner, doWriteMaskO=doWriteMask, uStr="deg")
        endif
        if (present(gStrO)) then
            call WriteVars(IOVars, .false., outF, trim(gStrO))
        else
            call WriteVars(IOVars, .false., outF, '/TubeShell')
        endif

    end subroutine writeTubeShell


    subroutine readTubeShellRestart(tubeShell, ResF, gStrO)
        type(TubeShell_T), intent(inout) :: tubeShell
        character(len=strLen), intent(in) :: ResF
        character(len=strLen), intent(in), optional :: gStrO

        character(len=strLen) :: gStr
        
        if (present(gStrO)) then
            gStr = gStrO
        else
            gStr = "/TubeShell"
        endif

        call ReadInSGV(tubeShell%xyz0     , ResF, "xyz0"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%X_bmin   , ResF, "X_bmin"   , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%X_maxRad , ResF, "X_maxRad" , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%lat0     , ResF, "lat0"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%lon0     , ResF, "lon0"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%latc     , ResF, "latc"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%lonc     , ResF, "lonc"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%J11i     , ResF, "J11i"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%invlat   , ResF, "invlat"   , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%topo     , ResF, "topo"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%bmin     , ResF, "bmin"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%Jxyz_iono, ResF, "Jxyz_iono", gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%Bxyz_iono, ResF, "Bxyz_iono", gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%bVol     , ResF, "bVol"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%Lb       , ResF, "Lb"       , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%Tb       , ResF, "Tb"       , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%wMAG     , ResF, "wMAG"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%rCurv    , ResF, "rCurv"    , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%avgBeta  , ResF, "avgBeta"  , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%avgP     , ResF, "avgP"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%avgN     , ResF, "avgN"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%stdP     , ResF, "stdP"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%stdN     , ResF, "stdN"     , gStrO=gStr, doIOpO=.false.)
        if (ioExist(ResF,"bcP",gStr)) then
            call ReadInSGV(tubeShell%bcP      , ResF, "bcP"      , gStrO=gStr, doIOpO=.false.)
            call ReadInSGV(tubeShell%bcN      , ResF, "bcN"      , gStrO=gStr, doIOpO=.false.)
        else
            !Make do with what you got
            call ReadInSGV(tubeShell%bcP      , ResF, "avgP"     , gStrO=gStr, doIOpO=.false.)
            call ReadInSGV(tubeShell%bcN      , ResF, "avgN"     , gStrO=gStr, doIOpO=.false.)
        endif            
        call ReadInSGV(tubeShell%losscone , ResF, "losscone" , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%lossconec, ResF, "lossconec", gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%TioTe0   , ResF, "Tiote0"   , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%nTrc     , ResF, "nTrc"     , gStrO=gStr, doIOpO=.false.)
        
    end subroutine readTubeShellRestart
end module tubehelper
