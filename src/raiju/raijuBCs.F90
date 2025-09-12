module raijuBCs

    use planethelper
    use raijudefs
    use raijutypes
    use raijuetautils
    use raijudomain

    implicit none

    contains

    subroutine applyRaijuBCs(Model, Grid, State, doWholeDomainO)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        logical, optional, intent(in) :: doWholeDomainO
            !! Whether we should apply certain BCs (moments2eta mapping) to entire active domain or not

        logical :: doWholeDomain
        integer :: s
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: doMomentIngest
            !! Whether we do moments to eta mapping at given grid cell

        if (present(doWholeDomainO)) then
            doWholeDomain = doWholeDomainO
        else
            doWholeDomain = .false.
        endif

        ! Now that topo is set, we can calculate active domain
        call setActiveDomain(Model, Grid, State)

        call calcMomentIngestionLocs(Model, Grid, State, doWholeDomain, doMomentIngest)
        call applyMomentIngestion(Model, Grid, State, doMomentIngest)


        if (Model%doActiveShell ) then
            ! Do first round of determining active shells for each k
            ! We do this everywhere
            do s=1,Grid%nSpc
                call setActiveShellsByContribution(Grid%shGrid, Grid%spc(s), State%bVol, &
                    State%eta(:,:,Grid%spc(s)%kStart:Grid%spc(s)%kEnd), &
                    State%activeShells(:,Grid%spc(s)%kStart:Grid%spc(s)%kEnd), &  ! This is what we're writing to
                    worthyFracO=Model%worthyFrac)
            enddo
        endif

    end subroutine applyRaijuBCs


    subroutine calcMomentIngestionLocs(Model, Grid, State, doWholeDomain, doMomentIngest)
        !! Determine which locations will have etas set by MHD moments
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        logical, intent(in) :: doWholeDomain
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: doMomentIngest
        
        integer :: i,j

        doMomentIngest = .false.
        ! Determine where to do BCs
        if(doWholeDomain) then
            where (State%active .ne. RAIJUINACTIVE)
                doMomentIngest = .true.
            end where
        else

            associate(sh=>Grid%shGrid)

            do j=sh%jsg,sh%jeg
                do i=sh%isg,sh%ie  ! Note: not evaluating low lat ghost cells, we don't wanna reset them with incoming plasma for now
                    
                    ! All buffer cells get set to moments
                    if (State%active(i,j) .eq. RAIJUBUFFER) then
                        doMomentIngest(i,j) = .true.
                    ! If a cell is currently active but it wasn't previously, we want to give it the freshest moment information
                    elseif (State%active(i,j) .eq. RAIJUACTIVE .and. State%active_last(i,j) .ne. RAIJUACTIVE) then
                        doMomentIngest(i,j) = .true.
                    ! If we are a not-invalid i-dir ghost cell we get reset no matter what
                    ! Note: This is only for high-lat boundary. Low lat does not get reset for now
                    elseif ((State%active(i,j) .ne. RAIJUINACTIVE) .and. (i < sh%is ) ) then
                        doMomentIngest(i,j) = .true.
                    endif

                enddo
            enddo

            end associate

        endif

    end subroutine calcMomentIngestionLocs


    subroutine applyMomentIngestion(Model, Grid, State, doMomentIngest)
        !! Based on doMomentIngest mask, map moments to eta channels for all species
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: doMomentIngest

        integer :: i,j,fm,fIdx
        integer :: psphIdx, eleIdx
            !! Index of plasmasphere species
        integer :: s
        real(rp) :: kT,vm
        real(rp) :: etaBelow
            !! Amount of eta below lowest lambda bound (every i,j gets a copy)
        real(rp) :: tmp_D, tmp_P
        real(rp) :: tmp_kti, tmp_kte
        real(rp) :: eMin

        psphIdx = spcIdx(Grid, F_PSPH)
        eleIdx  = spcIdx(Grid, F_HOTE)
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,s,fIdx,fm,vm,kT,etaBelow,tmp_kti,tmp_kte,eMin,tmp_D,tmp_P)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (State%active(i,j) .eq. RAIJUINACTIVE) then
                    State%eta(i,j,:) = 0.0
                endif
                if(.not. doMomentIngest(i,j)) then
                    cycle  ! This cycle should be okay because its inside the second loop
                endif

                vm = State%bvol_cc(i,j)**(-2./3.)

                ! Before we map to any RAIJU species, zero them out in case we want to accumulate
                do s=1,Model%nSpc
                    if (Grid%spc(s)%isMappedTo) then
                        State%eta(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd) = 0.0
                    endif

                    !!!!!!!!
                    !! TODO: Implement proper electron mapping
                    !!!!!!!!
                    if (Grid%spc(s)%flav .eq. F_HOTE) then
                        State%eta(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd) = 0.0
                    endif
                enddo
                
                ! Now go ahead and do mapping
                do fm=1,Model%nFluidIn
                    fIdx = Model%fluidInMaps(fm)%idx_mhd
                    s = spcIdx(Grid, Model%fluidInMaps(fm)%flav)
                    
                    ! Calculate moment we want to use for this point
                    call getDomWeightedMoments(i, j, fIdx, tmp_D, tmp_P)

                    kT = DP2kT(tmp_D, tmp_P)  ! [keV]
                    !!!!!!!!!!!!!
                    !! TODO: Implement proper Te map calculation
                    !!!!!!!!!!!!!
                    if (Grid%spc(s)%spcType .eq. RAIJUHPLUS) then
                        tmp_kti = kT / (1.0 + 1.0/Model%tiote)
                        tmp_kte = kT / (1.0 + Model%tiote)

                        eMin = Grid%spc(s)%alami(Grid%spc(s)%kStart)*vm*1e-3  ! [kev] Lowest energy for this channel

                        if (tmp_kti < eMin) then
                            ! Probably plasmasphere, dump it all into there and don't map to hot channels
                            etaBelow = tmp_D*State%bvol_cc(i,j)*sclEta  ! Basically the inverse of etak2Den in raijuEtautils
                        else

                            call DkT2SpcEta(Model,Grid%spc(s), &
                                State%eta(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd), &
                                tmp_D, tmp_kti, &
                                vm, doAccumulateO=.true., etaBelowO=etaBelow)
                            call DkT2SpcEta(Model,Grid%spc(eleIdx), &
                                State%eta(i,j,Grid%spc(eleIdx)%kStart:Grid%spc(eleIdx)%kEnd), &
                                tmp_D, tmp_kte, &
                                vm, doAccumulateO=.true.)
                        endif


                    else
                        call DkT2SpcEta(Model,Grid%spc(s), &
                            State%eta(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd), &
                            tmp_D, kT, &
                            vm, doAccumulateO=.true., etaBelowO=etaBelow)
                    endif

                    ! etaBelow has the amount of eta that is below the lowest lambda channel bound
                    !! TODO: Check to see if we are missing too much pressure
                    if (Model%doPsphEvol .and. Model%doExcesstoPsph .and. Model%fluidInMaps(fm)%doExcessToPsph) then
                        State%eta(i,j,Grid%spc(psphIdx)%kStart) = State%eta(i,j,Grid%spc(psphIdx)%kStart) + etaBelow
                    endif
                    
                enddo  ! f
            enddo  ! j
        enddo  ! i


        contains

        subroutine getDomWeightedMoments(i, j, fIdx, D, P)
            integer, intent(in) :: i, j, fIdx
            real(rp), intent(out) :: D, P

            real(rp), dimension(5,5) :: den2D, press2D, wgt2D
            logical, dimension(5,5) :: isG2D

            if (abs(State%domWeights(i,j) - 1.0_rp) < TINY) then
                D = State%Davg(i, j, fIdx)
                P = State%Pavg(i, j, fIdx)
            else
                ! Do weighted average by neighbors
                den2D   = State%Davg      (i-2:i+2, j-2:j+2, fIdx)
                press2D = State%Pavg      (i-2:i+2, j-2:j+2, fIdx)
                wgt2D   = State%domWeights(i-2:i+2, j-2:j+2)
                isG2D   = State%active    (i-2:i+2, j-2:j+2) .ne. RAIJUINACTIVE
                D = sum(den2D  * wgt2D, mask=isG2D)/sum(wgt2D, mask=isG2D)
                P = sum(press2D* wgt2D, mask=isG2D)/sum(wgt2D, mask=isG2D)
            endif
        end subroutine getDomWeightedMoments

    end subroutine applyMomentIngestion

!------
! Active I Shell calculations
!------

    subroutine setActiveShellsByContribution(shGrid, spc, bVol, etas, activeShellsOut, nSpacesO, worthyFracO)
        !! For each lambda channel, calculate active shells based on how much they contribute to the total pressure and/or density
        type(ShellGrid_T), intent(in) :: shGrid
        type(raijuSpecies_T), intent(in) :: spc
        real(rp), dimension(shGrid%isg:shGrid%ieg,shGrid%jsg:shGrid%jeg), intent(in) :: bVol
        real(rp), dimension(shGrid%isg:shGrid%ieg,shGrid%jsg:shGrid%jeg,spc%kStart:spc%kEnd), intent(in) :: etas
        logical, dimension(shGrid%isg:shGrid%ieg, spc%kStart:spc%kEnd), intent(inout) :: activeShellsOut
        integer, optional, intent(in) :: nSpacesO
            !! Number of i spaces between last good value and active i for species
        real(rp), optional, intent(in) :: worthyFracO

        integer :: i,j,k, iL, iU
        integer :: nSpaces
        real(rp) :: worthyFrac
        real(rp) :: alamc, kDen, kPress
        real(rp), dimension(shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg) :: spcDen, spcPress
        logical, dimension(shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg) :: as2D

        if(present(nSpacesO)) then
            nSpaces = nSpacesO
        else
            nSpaces = nSpacesDef
        endif

        if(present(worthyFracO)) then
            worthyFrac = worthyFracO
        else
            worthyFrac = fracWorthyDef
        endif
        
        as2D = .false.
        activeShellsOut = .false.

        ! Start by getting total density and pressure for each i,j point
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j)
        do j=shGrid%jsg,shGrid%jeg
            do i=shGrid%isg,shGrid%ieg
                spcDen  (i,j) = SpcEta2Den  (spc, etas(i,j,:), bVol(i,j))
                spcPress(i,j) = spcEta2Press(spc, etas(i,j,:), bVol(i,j))
            enddo
        enddo

        ! Then calculate active shells for each k
        do k=spc%kStart,spc%kEnd

            ! Setup for this k
            as2D = .false.
            alamc = 0.5*abs(spc%alami(k) + spc%alami(k+1))

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,kDen,kPress)
            do j=shGrid%jsg,shGrid%jeg
                do i=shGrid%isg,shGrid%ieg

                    kDen   = etak2Den  (etas(i,j,k), bVol(i,j))
                    kPress = etak2Press(etas(i,j,k), alamc, bVol(i,j))
            
                    if (kDen/spcDen(i,j) > worthyFrac .or. kPress/spcPress(i,j) > worthyFrac) then
                        as2D(i,j) = .true.
                    endif
                enddo
            enddo

            ! Now collapse j
            !! Re-use as2D, first j element 
            do i=shGrid%isg,shGrid%ieg
                as2D(i,shGrid%jsg) = any(as2D(i,:))
            enddo

            ! Finally, assign to activeShells, including buffer spaces
            do i=shGrid%isg,shGrid%ieg
                iL = max(i-nSpaces, shGrid%isg)
                iU = min(i+nSpaces, shGrid%ieg)
                activeShellsOut(i,k) = any(as2D(iL:iU, shGrid%jsg))
            enddo

        enddo
    end subroutine setActiveShellsByContribution


end module raijuBCs
