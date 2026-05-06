module raijuRecon
    !! Reconstruction routines
    use raijudefs
    use raijuTypes

    implicit none

    integer, parameter :: maxOrderSupported = 8
    real(rp), dimension(raiRecLen), parameter :: interpWgt_4c = [ 0, 0,  -1,  7,  7,  -1, 0, 0]/60.0_rp
    real(rp), dimension(raiRecLen), parameter :: interpWgt_6c = [ 0, 1,  -8, 37, 37,  -8, 1, 0]/60.0_rp
    real(rp), dimension(raiRecLen), parameter :: interpWgt_8c = [-3,29,-139,533,533,-139,29,-3]/840.0_rp

    real(rp), dimension(7), parameter :: C7Up = [-3,25,-101,319,214,-38,4]/420.0_rp
    real(rp), dimension(5), parameter :: C5Up = [2,-13,47,27,-3]/60.0_rp
    real(rp), dimension(3), parameter :: C3Up = [-1,5,2]/6.0_rp

    contains


!------
! Reconstruction stencil routines
!------

    !> 8th order central interpolation for a single interface
    function Central8(Qcc) result(Qi)
        real(rp), dimension(raiRecLen), intent(in) :: Qcc
        
        real(rp) :: Qi

        Qi = dot_product(interpWgt_8c,Qcc(1:raiRecLen))

    end function Central8
    

    !> 6th order central interpolation for a single interface
    function Central6(Qcc) result(Qi)
        !! Note: still takes 8-element stencil
        real(rp), dimension(raiRecLen), intent(in) :: Qcc
        
        real(rp) :: Qi

        Qi = dot_product(interpWgt_6c,Qcc(1:raiRecLen))

    end function Central6


    !> 4th order central interpolation for a single interface
    function Central4(Qcc) result(Qi)
        !! Note: still takes 8-element stencil
        real(rp), dimension(raiRecLen), intent(in) :: Qcc
        
        real(rp) :: Qi

        Qi = dot_product(interpWgt_4c,Qcc(1:raiRecLen))

    end function Central4


    function Upwind7(q1, q2, q3, q4, q5, q6, q7) result(qi)
        real(rp), intent(in) :: q1, q2, q3, q4, q5, q6, q7
        real(rp) :: qi

        qi = q1*C7Up(1) + q2*C7Up(2) + q3*C7Up(3) + q4*C7Up(4) &
           + q5*C7Up(5) + q6*C7Up(6) + q7*C7Up(7)
    end function Upwind7


    function Upwind5(q1, q2, q3, q4, q5) result(qi)
        real(rp), intent(in) :: q1, q2, q3, q4, q5
        real(rp) :: qi

        qi = q1*C5Up(1) + q2*C5Up(2) + q3*C5Up(3) + q4*C5Up(4) &
           + q5*C5Up(5)
    end function Upwind5

    function Upwind3(q1, q2, q3) result(qi)
        real(rp), intent(in) :: q1, q2, q3
        real(rp) :: qi

        qi = q1*C3Up(1) + q2*C3Up(2) + q3*C3Up(3)
    end function Upwind3


!------
! PDM implementations
!------
    subroutine raijuPDM(qm, q0, qp, Qi, pdmb, Qpdm)
        !! Applies PDM limiting on a single interface value from one direction
        real(rp), intent(in) :: qm, q0, qp
            !! Cell-centered values at donor cell (q0), cell behind donor (qm) and in front of donor (qp)
        real(rp), intent(in) :: Qi
            !! Reconstructed interface value, between q0 and qo
        real(rp), intent(in) :: pdmb
            !! PDM limiter value to use
        real(rp), intent(out) :: Qpdm
            !! PDM-ed value we return

        !! Using Left notation, but real direction is determined by ordering of arguments given to us
        real(rp) :: minQ, maxQ, Qn
        real(rp) :: dQL, dQi, sL, si, qL
        real(rp) :: dQnL

        ! Determine which is the intermediate value between Qi, Qi+1/2, Qi+1
        minQ = min(q0, qp)
        maxQ = max(q0, qp)
        Qn   = max(minQ, min(Qi, maxQ))

        ! Calculate differences to either side of interface, and at interface
        dQL = pdmb*( q0 - qm )
        dQi = pdmb*( qp - q0 )

        ! Signs of differences
        sL = sign(1.0_rp, dQL)
        si = sign(1.0_rp, dQi)

        ! Will be 0 if local extrema, or 2 if monotonic
        qL = abs(sL + si)

        ! Slope between median Qn and L cell centers
        ! Note, may be 0 in the case where Qi was the median value
        dQnL = Qn - q0

        ! Finally, calculate L reconstructed value
        Qpdm = Qn - si*max(0.0, abs(dQnL) - 0.5*qL*abs(dQL))

    end subroutine raijuPDM


    subroutine raijuPDMLR(Qcc, Qi, pdmb, QpdmL, QpdmR)
        !! Calculate PDM-limited value of Q at interface
        !! Using notation where interface is i+1/2
        !! Qcc is 4 elements long, i-1, i, i+1, i+2
        !!   Note, this is independent from indexing of arrays when they're handed to this subroutine
        !! Left is < i+1/2, Right is > i+1/2
        real(rp), dimension(4), intent(in) :: Qcc
            !! Cell-centered values of quantity. First two are left of interface, last 2 are right of interface
        real(rp), intent(in) :: Qi
            !! High-order reconstruction of quantity at interface
        real(rp), intent(in) :: pdmb
            !! PDM limiter value
        real(rp), intent(out) :: QpdmL, QpdmR
            !! PDM-limited quantity to return

        real(rp) :: minQ, maxQ, Qn
        real(rp) :: dQL, dQi, dQR, sL, si, sR, qL, qR
        real(rp) :: dQnL, dQnR

        ! Determine which is the intermediate value between Qi, Qi+1/2, Qi+1
        minQ = min(Qcc(2), Qcc(3))
        maxQ = max(Qcc(2), Qcc(3))
        Qn   = max(minQ, min(Qi, maxQ))

        ! Calculate differences to either side of interface, and at interface
        dQL = pdmb*( Qcc(2) - Qcc(1) )
        dQi = pdmb*( Qcc(3) - Qcc(2) )
        dQR = pdmb*( Qcc(4) - Qcc(3) )

        ! Signs of differences
        sL = sign(1.0_rp, dQL)
        si = sign(1.0_rp, dQi)
        sR = sign(1.0_rp, dQR)

        ! Will be 0 if local extrema, or 2 if monotonic
        qL = abs(sL + si)
        qR = abs(si + sR)

        ! Slope between median Qn and L/R cell centers
        ! Note, may be 0 in the case where Qi or Qi+1 was the median value
        dQnL = Qn     - Qcc(2)
        dQnR = Qcc(3) - Qn

        ! Finally, calculate L/R reconstructed values
        QpdmL = Qn - si*max(0.0, abs(dQnL) - 0.5*qL*abs(dQL))
        QpdmR = Qn + si*max(0.0, abs(dQnR) - 0.5*qR*abs(dQR))

    end subroutine raijuPDMLR

!------
! Interface reconstruction
!------
! ReconFaces: responsible for doing full face reconstruction procedure, including pdm, for all active cells
    !  -> reconFaceLR: is given 8-cell stencil, decides which order to use, applies pdm, returns L/R values

    subroutine ReconFaceLR(isG, Qcc, areaCC, areaFace, BrCC, BrFace, maxOrder, pdmb, QfaceL, QfaceR, QreconLO, QreconRO)
        !! Receives an 8-element stencil and returns reconstructed L/R values of a single face
        logical , dimension(raiRecLen), intent(in) :: isG
        real(rp), dimension(raiRecLen), intent(in) :: Qcc
            !! Cell centers on either side of interface
            !! (Interface is between cells 4 and 5)
        real(rp), dimension(raiRecLen), intent(in) :: areaCC
            !! Cell-centered areas
        real(rp), intent(in) :: areaFace
            !! Estimated cell area at the face we are reconstructing at
        real(rp), dimension(raiRecLen), intent(in) :: BrCC
            !! Cell-centered radial magnetic field
        real(rp), intent(in) :: BrFace
            !! Estimated radial magnetic field at the face we are reconstructing at
        integer , intent(in) :: maxOrder
            !! Maximum order to use for reconstruction
        real(rp), intent(in) :: pdmb
            !! PDM limiter value to pass along to PDM functions
        real(rp), intent(inout) :: QfaceL, QfaceR
            !! L/R face values we return
        real(rp), optional, intent(inout) :: QreconLO, QreconRO
            !! We will optionally return the reconstructed states if these variables are provided to us

        real(rp), dimension(raiRecLen) :: QccA
            !! Area-weighted cell-centered quantity
        real(rp) :: QreconL, QreconR
            !! Face values after reconstruction but before limiting

        logical :: doUpwind
            !! Whether we will do upwind schemes or centered schemes

        !! If we are given an odd max order, we assume we should use an odd order even if it must be reduced
        !! Same for even max order

        QfaceL = 0.0
        QfaceR = 0.0
        QreconL = 0.0
        QreconR = 0.0

        doUpwind = mod(maxOrder,2)==1

        QccA = Qcc*areaCC*BrCC  ! [Q * Rp^2 * nT]
        if (doUpwind) then
            ! Starting at highest available order, we first see if caller wants us to use this order or higher 
            ! and if cells are good enough to use this order.
            ! If not, we keep going till we land at an acceptable order
            if (maxOrder == 7 .and. all(isG)) then
                !! Entirety of stencil must be good because we need to do both sides of the interface
                QreconL = Upwind7(QccA(1), QccA(2), QccA(3), QccA(4), QccA(5), QccA(6), QccA(7))/areaFace/BrFace
                QreconR = Upwind7(QccA(8), QccA(7), QccA(6), QccA(5), QccA(4), QccA(3), QccA(2))/areaFace/BrFace
                call raijuPDM(Qcc(3), Qcc(4), Qcc(5), QreconL, pdmb, QfaceL)
                call raijuPDM(Qcc(6), Qcc(5), Qcc(4), QreconR, pdmb, QfaceR)
            else if (maxOrder >= 5 .and. all(isG(2:7))) then
                QreconL = Upwind5(QccA(2), QccA(3), QccA(4), QccA(5), QccA(6))/areaFace/BrFace
                QreconR = Upwind5(QccA(7), QccA(6), QccA(5), QccA(4), QccA(3))/areaFace/BrFace
                call raijuPDM(Qcc(3), Qcc(4), Qcc(5), QreconL, pdmb, QfaceL)
                call raijuPDM(Qcc(6), Qcc(5), Qcc(4), QreconR, pdmb, QfaceR)
            else if (maxOrder >= 3 .and. all(isG(3:6))) then
                QreconL = Upwind3(QccA(3), QccA(4), QccA(5))/areaFace/BrFace
                QreconR = Upwind3(QccA(6), QccA(5), QccA(4))/areaFace/BrFace
                call raijuPDM(Qcc(3), Qcc(4), Qcc(5), QreconL, pdmb, QfaceL)
                call raijuPDM(Qcc(6), Qcc(5), Qcc(4), QreconR, pdmb, QfaceR)
            else if (all(isG(4:5))) then
                QreconL = QccA(4)/areaFace/BrFace
                QreconR = QccA(5)/areaFace/BrFace
            endif
        else
            ! Even orders
            ! Note: Because directionality doesn't matter for centered, we don't need to make sure we order things correctly
            ! and we can just do a dot product on the whole 8-cell stencil. Any bad points in lower-order stencils will be zero-ed out
            ! Also, left and right states are equal, so we only save the reconstructed value to QreconL,
            !  and give it to a different version of PDM that will split up the single value into L/R states
            if(maxOrder == 8 .and. all(isG)) then
                QreconL = Central8(QccA)/areaFace/BrFace
                call raijuPDMLR(Qcc(3:6), QreconL, pdmb, QfaceL, QfaceR)
            else if(maxOrder >= 6 .and. all(isG(2:7))) then
                QreconL = Central6(QccA)/areaFace/BrFace
                call raijuPDMLR(Qcc(3:6), QreconL, pdmb, QfaceL, QfaceR)
            else if(maxOrder >= 4 .and. all(isG(3:6))) then
                QreconL = Central4(QccA)/areaFace/BrFace
                call raijuPDMLR(Qcc(3:6), QreconL, pdmb, QfaceL, QfaceR)
            else if(maxOrder >= 2 .and. all(isG(4:5))) then
                QreconL = 0.5*(QccA(4) + QccA(5))/areaFace/BrFace
                QfaceL = QreconL
                QfaceR = QreconL
            endif
        endif


        if (present(QreconLO)) then
            QreconLO = QreconL
        endif

        if (present(QreconRO)) then
            QreconRO = QreconR
        endif

    end subroutine ReconFaceLR


    subroutine ReconFaces(Model, Grid, isG, Qcc, QfaceL, QfaceR, QreconLO, QreconRO, doOMPO)
        !! Performs full face reconstruction procedure, including pdm, on Qcc for all active cells
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isG
                           !! Whether a cell is safe to use in reconstruction
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: Qcc
                           !! Cell-centered variable to interpolate to faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: QfaceL, QfaceR
                           !! Left/Right face-interpolated values
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), optional, intent(inout) :: QreconLO, QreconRO
        logical, optional, intent(in) :: doOMPO

        integer :: i,j
        logical :: doOMP

        if (present(doOMPO)) then
            doOMP = doOMPO
        else
            doOMP = .false.
        endif
        QfaceL = 0.0
        QfaceR = 0.0

        !! Note: We only populate from is:ie+1, js:je+1 because that's the only place we'll use face values
        !!  But full size includes ghost faces for output purposes

        
        ! I just wanna start by saying I'm not a fan of this
        ! But I think its important to make it relatively easy to debug
        ! And if we do it this way then we don't need to do if conditions on every loop
        if (present(QreconLO) .and. present(QreconRO)) then
            QreconLO = 0.0
            QreconRO = 0.0
            !!$OMP PARALLEL DO default(shared) &
            !!$OMP schedule(dynamic) &
            !!$OMP private(i,j) &
            !!$OMP IF(doOMP)
            do j=Grid%shGrid%js,Grid%shGrid%je+1
                do i=Grid%shGrid%is,Grid%shGrid%ie+1
                    ! Theta dir
                    !ReconFaceLR(isG, Qcc, areaCC, areaFace, maxOrder, pdmb, QfaceL, QfaceR, QreconRO, QreconLO)
                    call ReconFaceLR(isG       (i-4:i+3,j), Qcc(i-4:i+3,j)           , &
                                    Grid%areaCC(i-4:i+3,j), Grid%areaFace(i,j,RAI_TH), &
                                    Grid%Brcc  (i-4:i+3,j), Grid%BrFace  (i,j,RAI_TH), &
                                    Model%maxOrder        , Model%PDMB               , &
                                    QfaceL  (i,j,RAI_TH)  , QfaceR  (i,j,RAI_TH)     , &
                                    QreconLO(i,j,RAI_TH)  , QreconRO(i,j,RAI_TH)     )
                    ! Phi dir
                    call ReconFaceLR(isG       (i,j-4:j+3), Qcc(i,j-4:j+3)          , &
                                    Grid%areaCC(i,j-4:j+3), Grid%areaFace(i,j,RAI_PH), &
                                    Grid%Brcc  (i,j-4:j+3), Grid%BrFace  (i,j,RAI_PH), &
                                    Model%maxOrder        , Model%PDMB               , &
                                    QfaceL  (i,j,RAI_PH)  , QfaceR  (i,j,RAI_PH)     , &
                                    QreconLO(i,j,RAI_PH)  , QreconRO(i,j,RAI_PH)     )
                enddo
            enddo

        else
            !!$OMP PARALLEL DO default(shared) &
            !!$OMP schedule(dynamic) &
            !!$OMP private(i,j) &
            !!$OMP IF(doOMP)
            do j=Grid%shGrid%js,Grid%shGrid%je+1
                do i=Grid%shGrid%is,Grid%shGrid%ie+1
                    ! Theta dir
                    !ReconFaceLR(isG, Qcc, areaCC, areaFace, maxOrder, pdmb, QfaceL, QfaceR)
                    call ReconFaceLR(isG       (i-4:i+3,j), Qcc(i-4:i+3,j)           , &
                                    Grid%areaCC(i-4:i+3,j), Grid%areaFace(i,j,RAI_TH), &
                                    Grid%Brcc  (i-4:i+3,j), Grid%BrFace  (i,j,RAI_TH), &
                                    Model%maxOrder        , Model%PDMB               , &
                                    QfaceL(i,j,RAI_TH)    , QfaceR(i,j,RAI_TH)       )
                    ! Phi dir
                    call ReconFaceLR(isG       (i,j-4:j+3), Qcc(i,j-4:j+3)           , &
                                    Grid%areaCC(i,j-4:j+3), Grid%areaFace(i,j,RAI_PH), &
                                    Grid%Brcc  (i,j-4:j+3), Grid%BrFace  (i,j,RAI_PH), &
                                    Model%maxOrder        , Model%PDMB               , &
                                    QfaceL(i,j,RAI_PH)    , QfaceR(i,j,RAI_PH)       )
                enddo
            enddo
        endif

    end subroutine ReconFaces

!------
! Flux calculation
!------

    subroutine calcBoundaryFluxes(sh, active, Qflux)
        !! Sets fluxes at invalid/buffer boundaries
        type(ShellGrid_T), intent(in) :: sh
        integer , dimension(sh%isg:sh%ieg, &
                            sh%jsg:sh%jeg)     , intent(in) :: active
        real(rp), dimension(sh%isg:sh%ieg+1, &
                            sh%jsg:sh%jeg+1, 2), intent(inout) :: Qflux

        integer :: i,j

        do j=sh%js,sh%je+1
            do i=sh%is,sh%ie+1
                
                ! Theta dir
                if (  (active(i-1,j       ) .eq. RAIJUINACTIVE) &
                .and. (active(i  ,j       ) .eq. RAIJUBUFFER  ) &
                .and. (Qflux (i+1,j,RAI_TH) < 0.0 ) ) then
                ! If current cell is good, lower i is bad, only copy i+1 flux if we have outflow
                    Qflux(i,j,RAI_TH) = Qflux(i+1,j,RAI_TH)
                else if (  (active(i-1,j       ) .eq. RAIJUBUFFER  ) &
                     .and. (active(i  ,j       ) .eq. RAIJUINACTIVE) &
                     .and. (Qflux (i-1,j,RAI_TH) > 0.0 ) ) then
                ! If current cell is bad, lower i is good, only copy i-1 flux if we have outflow
                    Qflux(i,j,RAI_TH) = Qflux(i-1,j,RAI_TH)
                endif

                ! Psi dir
                if (  (active(i,j-1       ) .eq. RAIJUINACTIVE) &
                .and. (active(i,j         ) .eq. RAIJUBUFFER  ) &
                .and. (Qflux (i,j+1,RAI_PH) < 0.0 ) ) then
                ! If current cell is good, lower j is bad, only copy j+1 flux if we have outflow
                    Qflux(i,j,RAI_PH) = Qflux(i,j+1,RAI_PH)
                else if (  (active(i,j-1       ) .eq. RAIJUBUFFER  ) &
                     .and. (active(i,j         ) .eq. RAIJUINACTIVE) &
                     .and. (Qflux (i,j-1,RAI_PH) > 0.0 ) ) then
               ! If current cell is bad, lower j is good, only copy j-1 flux if we have outflow
                    Qflux(i,j,RAI_PH) = Qflux(i,j-1,RAI_PH)
               endif

            enddo
        enddo

    end subroutine calcBoundaryFluxes


    subroutine calcFluxes(Model, Grid, State, k, isGoodRecon, Qcc, Qflux)
        !! Takes cell-centered quantity (Qcc) and uses state information to calculate flux of Q through cell faces (Qflux)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isGoodRecon
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: Qcc
            !! Cell-centered quantity (eta)
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: Qflux
            !! Flux of Q through faces
        
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: QfaceL, QfaceR
        real(rp) :: QfluxL, QfluxR
        integer :: i,j
        
        QfaceL = 0.0
        QfaceR = 0.0
        QfluxL = 0.0
        QfluxR = 0.0
        Qflux  = 0.0
        
        ! ReconFaces(Model, Grid, isG, Qcc, QfaceL, QfaceR, QreconLO, QreconRO)
        if (Model%doOutput_debug) then
            call ReconFaces(Model, Grid, isGoodRecon, Qcc, QfaceL, QfaceR, State%etaFaceReconL(:,:,k,:), State%etaFaceReconR(:,:,k,:))
        else
            call ReconFaces(Model, Grid, isGoodRecon, Qcc, QfaceL, QfaceR)
        endif

        if (Model%doUseVelLRs) then
            do j=Grid%shGrid%js,Grid%shGrid%je+1
                do i=Grid%shGrid%is,Grid%shGrid%ie+1
                    ! Only allow velocities passing through the interface
                    QfluxL = QfaceL(i,j,RAI_TH)*max(0.0, State%iVelL(i,j,k,RAI_TH))
                    QfluxR = QfaceR(i,j,RAI_TH)*min(0.0, State%iVelR(i,j,k,RAI_TH))
                    Qflux(i,j,RAI_TH) = QfluxL + QfluxR

                    QfluxL = QfaceL(i,j,RAI_PH)*max(0.0, State%iVelL(i,j,k,RAI_PH))
                    QfluxR = QfaceR(i,j,RAI_PH)*min(0.0, State%iVelR(i,j,k,RAI_PH))
                    Qflux(i,j,RAI_PH) = QfluxL + QfluxR
                enddo
            enddo
        else
            ! Since we know the exact velocity at the interface, 
            ! our final flux is QL*v if stuff is leaving the cell, or QR*v if stuff is entering the cell
            do j=Grid%shGrid%js,Grid%shGrid%je+1
                do i=Grid%shGrid%is,Grid%shGrid%ie+1
                    if (State%iVel(i,j,k,RAI_TH) > 0.0) then
                        Qflux(i,j,RAI_TH) = QfaceL(i,j,RAI_TH)*State%iVel(i,j,k,RAI_TH)
                    else
                        Qflux(i,j,RAI_TH) = QfaceR(i,j,RAI_TH)*State%iVel(i,j,k,RAI_TH)
                    endif

                    if (State%iVel(i,j,k,RAI_PH) > 0.0) then
                        Qflux(i,j,RAI_PH) = QfaceL(i,j,RAI_PH)*State%iVel(i,j,k,RAI_PH)
                    else
                        Qflux(i,j,RAI_PH) = QfaceR(i,j,RAI_PH)*State%iVel(i,j,k,RAI_PH)
                    endif
                enddo
            enddo
        endif
        !Qflux = merge(QfaceL*State%iVel(:,:,k,:), QfaceR*State%iVel(:,:,k,:), State%iVel(:,:,k,:) > 0.0)  ! [Q * m/s]
        Qflux = Qflux * Grid%BrFace / Model%planet%rp_m  ! [Q * nT * Rp/s]

        ! Thus far we have ignored fluxes of faces at invalid/buffer boundary
        !  (ReconFaces set them to zero)
        ! Now that valid faces are decided, we can apply our invalid/buffer boundary conditions
        !call calcBoundaryFluxes(Grid%shGrid, State%active, Qflux)

        if (Model%doOutput_debug) then
            State%etaFacePDML(:,:,k,:) = QfaceL
            State%etaFacePDMR(:,:,k,:) = QfaceR
            State%etaFlux    (:,:,k,:) = Qflux
        endif
        
    end subroutine calcFluxes

end module raijuRecon