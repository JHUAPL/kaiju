
!Helper routines to handle tube tracing and data extraction
module tubehelper
    use volttypes
    use voltcpltypes
    use streamline
    use chmpdbz
    use shellUtils
    use shellGridIO

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
            tiote = 4.0
        endif


        do i=1,NDIM
            call initShellVar(sh, SHGR_CORNER, tubeShell%xyz0(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%X_bmin(i))
        enddo
        call initShellVar(sh, SHGR_CORNER, tubeShell%lat0   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%lon0   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%latc   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%lonc   )
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

        !Figure out topology
        bTube%topo = VoltTopo(ebModel,ebGr,bTrc)
        if (bTube%topo == TUBE_DISASTER) return


    !If we're still here, let's do some stuff
        !Start w/ stuff for both closed and open
        bTube%TioTe0 = 4.0 !For now

        !Right now doing species loop for open or closed
        !TODO: Rewrite open to get near-earth average instead of field line average

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
        enddo

        if (bTube%topo == TUBE_CLOSED) then
        !Do closed field line stuff
            bTube%bVol = dvB/oBScl !Rx/EB => Rx/nT

        !Minimal surface (bEq in Rp, bMin in EB)
            call FLEq(ebModel,bTrc,bEq,bmin)
            bmin = bmin*oBScl !EB => nT
            bTube%bmin = bmin
            bTube%X_bmin = bEq

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

    function VoltTopo(ebModel,ebGr,bTrc) result(topo)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebGrid_T)   , intent(in) :: ebGr
        type(magLine_T  ), intent(in) :: bTrc
        integer :: topo

        integer :: iP,iM
        real(rp), dimension(NDIM) :: xyzP,xyzM
        real(rp) :: zz
        logical :: isFin,isCP,isCM,isOP,isOM

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
            !Get end points
            xyzP = bTrc%xyz(+Np,:)
            xyzM = bTrc%xyz(-Nm,:)

            !Check if finished
            isFin = (Np<bTrc%nMax-1) .and. (Nm<bTrc%nMax-1) !Check if finished
            if (.not. isFin) then
                topo = TUBE_UNDEF !Didn't finish, dunno what's going on
                return
            endif

            isCP  = isClosed(xyzP,ebModel)
            isCM  = isClosed(xyzM,ebModel)
            zz = xyzP(ZDIR)*xyzM(ZDIR)

            if (isCP .and. isCM .and. (zz<=0)) then
                topo = TUBE_CLOSED
                return
            endif

            if ( (.not. isCP) .and. (.not. isCM) ) then
                !what the shit?!
                topo = TUBE_DISASTER
                return
            endif

        !If still here, one side is closed but not the other
            !Get i cell of last point for each side
            iP = bTrc%ijk(+Np,IDIR)
            iM = bTrc%ijk(-Nm,IDIR)

            isOP = isOut(ebModel,ebGr,iP)
            isOM = isOut(ebModel,ebGr,iM)

            if (isCP .and. isOM) then
                topo = TUBE_OPEN
            elseif (isOP .and. isCM) then
                topo = TUBE_OPEN
            else
                topo = TUBE_DISASTER
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
        bTube%bVol = 0.0
        bTube%Lb = 0.0
        bTube%Tb = 0.0

        bTube%wMAG = 0.0
        bTube%rCurv = 0.0
        bTube%avgBeta = 0.0

        bTube%avgP(:) = 0.0
        bTube%stdP(:) = 0.0
        bTube%avgN(:) = 0.0
        bTube%stdN(:) = 0.0

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
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k)
        do j=shGr%js,shGr%je+1
            do i=shGr%is,shGr%ie+1
                do k=1,NDIM
                    tubeShell%xyz0  (k)%data(i,j) = ijTubes(i,j)%xyz0(k)
                    tubeShell%X_bmin(k)%data(i,j) = ijTubes(i,j)%X_bmin(k)
                enddo
                tubeShell%lat0   %data(i,j) = ijTubes(i,j)%lat0
                tubeShell%lon0   %data(i,j) = ijTubes(i,j)%lon0
                tubeShell%invlat %data(i,j) = ijTubes(i,j)%invlat
                tubeShell%latc   %data(i,j) = ijTubes(i,j)%latc
                tubeShell%lonc   %data(i,j) = ijTubes(i,j)%lonc
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
        enddo
        call wrapJ_SGV(shGr, tubeShell%lat0   )
        call wrapJ_SGV(shGr, tubeShell%lon0   )
        call wrapJ_SGV(shGr, tubeShell%invlat )
        call wrapJ_SGV(shGr, tubeShell%latc   )
        call wrapJ_SGV(shGr, tubeShell%lonc   )
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
            tubeShell%stdP(k)%mask(:,:) = topoMask
            tubeShell%avgN(k)%mask(:,:) = topoMask
            tubeShell%stdN(k)%mask(:,:) = topoMask
        enddo
        tubeShell%TioTe0%mask(:,:) = topoMask
        ! idk about losscone, lossconec
        ! nTrc is good
    end subroutine tubes2Shell


    subroutine writeTubeShellRestart(tubeShell, ResF, gStrO)
        type(TubeShell_T), intent(in) :: tubeShell
        character(len=strLen), intent(in) :: ResF
        character(len=strLen), intent(in), optional :: gStrO
            !! Group name to use instead of 'TubeShell' at root level

        type(IOVAR_T), dimension(50) :: IOVars
        
        call AddOutSGV(IOVars,"xyz0"     ,tubeShell%xyz0     , doWriteMaskO=.true., uStr="Rp")
        call AddOutSGV(IOVars,"X_bmin"   ,tubeShell%X_bmin   , doWriteMaskO=.true., uStr="Rp")
        call AddOutSGV(IOVars,"lat0"     ,tubeShell%lat0     , doWriteMaskO=.true., uStr="deg")
        call AddOutSGV(IOVars,"lon0"     ,tubeShell%lon0     , doWriteMaskO=.true., uStr="deg")
        call AddOutSGV(IOVars,"latc"     ,tubeShell%latc     , doWriteMaskO=.true., uStr="deg")
        call AddOutSGV(IOVars,"lonc"     ,tubeShell%lonc     , doWriteMaskO=.true., uStr="deg")
        call AddOutSGV(IOVars,"invlat"   ,tubeShell%invlat   , doWriteMaskO=.true., uStr="deg")
        call AddOutSGV(IOVars,"topo"     ,tubeShell%topo     , doWriteMaskO=.true., uStr="enum")
        call AddOutSGV(IOVars,"bmin"     ,tubeShell%bmin     , doWriteMaskO=.true., uStr="nT")
        call AddOutSGV(IOVars,"bVol"     ,tubeShell%bVol     , doWriteMaskO=.true., uStr="Rp/nT")
        call AddOutSGV(IOVars,"Lb"       ,tubeShell%Lb       , doWriteMaskO=.true., uStr="Rp")
        call AddOutSGV(IOVars,"Tb"       ,tubeShell%Tb       , doWriteMaskO=.true., uStr="s")
        call AddOutSGV(IOVars,"wMAG"     ,tubeShell%wMAG     , doWriteMaskO=.true., uStr="normalized")
        call AddOutSGV(IOVars,"rCurv"    ,tubeShell%rCurv    , doWriteMaskO=.true., uStr="Rp")
        call AddOutSGV(IOVars,"avgBeta"  ,tubeShell%avgBeta  , doWriteMaskO=.true., uStr="normalized")
        call AddOutSGV(IOVars,"avgP"     ,tubeShell%avgP     , doWriteMaskO=.true., uStr="nPa")
        call AddOutSGV(IOVars,"avgN"     ,tubeShell%avgN     , doWriteMaskO=.true., uStr="#/cc")
        call AddOutSGV(IOVars,"stdP"     ,tubeShell%stdP     , doWriteMaskO=.true., uStr="nPa")
        call AddOutSGV(IOVars,"stdN"     ,tubeShell%stdN     , doWriteMaskO=.true., uStr="#/cc")
        call AddOutSGV(IOVars,"losscone" ,tubeShell%losscone , doWriteMaskO=.true., uStr="rad")
        call AddOutSGV(IOVars,"lossconec",tubeShell%lossconec, doWriteMaskO=.true., uStr="rad")
        call AddOutSGV(IOVars,"Tiote0"   ,tubeShell%TioTe0   , doWriteMaskO=.true., uStr="normalized")
        call AddOutSGV(IOVars,"nTrc"     ,tubeShell%nTrc     , doWriteMaskO=.true., uStr="counts")
        
        
        if (present(gStrO)) then
            call WriteVars(IOVars, .false., ResF, trim(gStrO))
        else
            call WriteVars(IOVars, .false., ResF, '/TubeShell')
        endif

    end subroutine writeTubeShellRestart

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
        call ReadInSGV(tubeShell%lat0     , ResF, "lat0"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%lon0     , ResF, "lon0"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%latc     , ResF, "latc"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%lonc     , ResF, "lonc"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%invlat   , ResF, "invlat"   , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%topo     , ResF, "topo"     , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%bmin     , ResF, "bmin"     , gStrO=gStr, doIOpO=.false.)
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
        call ReadInSGV(tubeShell%losscone , ResF, "losscone" , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%lossconec, ResF, "lossconec", gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%TioTe0   , ResF, "Tiote0"   , gStrO=gStr, doIOpO=.false.)
        call ReadInSGV(tubeShell%nTrc     , ResF, "nTrc"     , gStrO=gStr, doIOpO=.false.)
    end subroutine readTubeShellRestart
end module tubehelper
