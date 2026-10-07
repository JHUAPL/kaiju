!Routines for voltron to calculate MI coupling quantities

module cmiutils
    use gamtypes
    use cmidefs
    use volttypes, only: IONO_EX,IONO_EY,IONO_EZ
    use metric
    use ringutils

    implicit none

    logical, parameter, private :: doRingFAC = .true. !Do some ring-processing on currents before sending to remix

    contains

    !-----------------------------
    !Convert electic potential from ionosphere to E fields for inner BCs
    subroutine Ion2MHD(Model,Grid,gPsi,inEijk,inExyz,pSclO)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), intent(inout) :: gPsi  (1:PsiSh+1,Grid%js:Grid%je+1  ,Grid%ks:Grid%ke+1)
        real(rp), intent(inout) :: inEijk(1:PsiSh+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM)
        real(rp), intent(inout) :: inExyz(1:PsiSh  ,Grid%jsg:Grid%jeg  ,Grid%ksg:Grid%keg,1:Model%nvInExyz)
        real(rp), intent(in), optional :: pSclO
        integer :: i,j,k,iG
        integer :: NumP
        real(rp) :: ijkDet, Mijk(NDIM,NDIM)
        real(rp), dimension(NDIM) :: iCC,jCC,kCC,ccEijk
        real(rp) :: ionScl

        !Set Eijk fields
        inEijk = 0.0
        inExyz = 0.0

        !Set scaling if present
        ionScl = 1.0
        if (present(pSclO)) then
            ionScl = pSclO
        endif

        NumP = Grid%ke-Grid%ks+1

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,iG) 
        do i=1,PsiSh+1
            do k=Grid%ks,Grid%ke+1
                do j=Grid%js,Grid%je+1

                    iG = i+PsiSt-1 !Global i index for Grid access

                    if (i <= PsiSh) then
                        inEijk(i,j,k,IDIR) = -( gPsi(i+1,j,k) - gPsi(i,j,k) )/Grid%edge(iG,j,k,IDIR)
                    endif
                    if (j <= Grid%je) then
                        inEijk(i,j,k,JDIR) = -( gPsi(i,j+1,k) - gPsi(i,j,k) )/Grid%edge(iG,j,k,JDIR)
                    endif
                    if (k <= Grid%ke) then
                        inEijk(i,j,k,KDIR) = -( gPsi(i,j,k+1) - gPsi(i,j,k) )/Grid%edge(iG,j,k,KDIR)
                    endif
                    !Scale field
                    inEijk(i,j,k,:) = inEijk(i,j,k,:)/ionScl

                enddo
            enddo
        enddo !Shell loop

        !$OMP PARALLEL DO default(shared)
        do i=1,PsiSh+1
            !Correct for pole
            !K/NOTE: Not killing Ek b/c small value properly cancels js+1,Ek in CT
            if (Model%Ring%doS) then
                !inEijk(i,Grid%js,:,KDIR) = 0.0
                inEijk(i,Grid%js,:,IDIR) = sum(inEijk(i,Grid%js,Grid%ks:Grid%ke,IDIR))/NumP
            endif
            if (Model%Ring%doE) then
                inEijk(i,Grid%je+1,:,JDIR) = 0.0
                !inEijk(i,Grid%je+1,:,KDIR) = 0.0
                inEijk(i,Grid%je+1,:,IDIR) = sum(inEijk(i,Grid%je+1,Grid%ks:Grid%ke,IDIR))/NumP
            endif
            !Enforce periodicity constraints (needed for differencing in next step)
            inEijk(i,:,Grid%ke+1,IDIR:KDIR) = inEijk(i,:,Grid%ks,IDIR:KDIR)

        enddo !Shell loop

        !Fill ghosts of inEijk
        call WrapBallEC(Model,Grid,inEijk)

        !Now turn nc-IJK FIELDS into cc-XYZ fields
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,iG,ccEijk,iCC,jCC,kCC,Mijk,ijkDet)
        do i=1,PsiSh
            do k=Grid%ks,Grid%ke
                do j=Grid%js,Grid%je
                    iG = i+PsiSt-1 !Global i index into grid arrays
                    !Get cell-centered XYZ fields
                  
                    ccEijk(IDIR) = 0.25*(inEijk(i,j,k,IDIR)+inEijk(i,j,k+1,IDIR)+inEijk(i,j+1,k,IDIR)+inEijk(i,j+1,k+1,IDIR))
                    ccEijk(JDIR) = 0.25*(inEijk(i,j,k,JDIR)+inEijk(i,j,k+1,JDIR)+inEijk(i+1,j,k,JDIR)+inEijk(i+1,j,k+1,JDIR))
                    ccEijk(KDIR) = 0.25*(inEijk(i,j,k,KDIR)+inEijk(i+1,j,k,KDIR)+inEijk(i,j+1,k,KDIR)+inEijk(i+1,j+1,k,KDIR))

                    !Get cell-centered ijk vectors @ cell-center
                    iCC = ijkVec(Model,Grid,iG,j,k,IDIR)
                    jCC = ijkVec(Model,Grid,iG,j,k,JDIR)
                    kCC = ijkVec(Model,Grid,iG,j,k,KDIR)

                    !Get inverse matrix
                    call ijkMatrix(Model,Grid,iCC,jCC,kCC,Mijk)

                    !Determinant of ijk matrix
                    ijkDet = iCC(XDIR)*Mijk(1,1) + jCC(XDIR)*Mijk(1,2) + kCC(XDIR)*Mijk(1,3)

                    inExyz(i,j,k,IONO_EX) = (Mijk(1,1)*ccEijk(IDIR) + Mijk(1,2)*ccEijk(JDIR) + Mijk(1,3)*ccEijk(KDIR))/ijkDet
                    inExyz(i,j,k,IONO_EY) = (Mijk(2,1)*ccEijk(IDIR) + Mijk(2,2)*ccEijk(JDIR) + Mijk(2,3)*ccEijk(KDIR))/ijkDet
                    inExyz(i,j,k,IONO_EZ) = (Mijk(3,1)*ccEijk(IDIR) + Mijk(3,2)*ccEijk(JDIR) + Mijk(3,3)*ccEijk(KDIR))/ijkDet

                enddo
            enddo
                !Recalculate inner-most ring (degenerate metric)
            if (Model%doRing) then
                if (Model%Ring%doS) then
                    call FixRAVec_S(Model,Grid,inExyz(i,Grid%js:Grid%js+1,Grid%ks:Grid%ke,IONO_EX:IONO_EZ))
                endif
                if (Model%Ring%doE) then
                    call FixRAVec_E(Model,Grid,inExyz(i,Grid%je-1:Grid%je,Grid%ks:Grid%ke,IONO_EX:IONO_EZ))
                endif
            endif

        enddo !Shell loop
        
        !Fill ghosts of inExyz
        call WrapBallCC(Model,Grid,inExyz)

    end subroutine Ion2MHD

    subroutine WrapBallCC(Model,Grid,inExyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Grid
        real(rp)  , intent(inout) :: inExyz(1:PsiSh,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:Model%nvInExyz)

        integer :: i,j,k,ip,jp,kp
        logical :: isInI,isInJ,isInK,isAct

        !Fill all ghosts, lazily loop through everything
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,ip,jp,kp,isInI,isInJ,isInK,isAct)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do i=1,PsiSh
                    !Check if loop iteration is in active
                    isInI = .true.
                    isInJ = (j>=Grid%js) .and. (j<=Grid%je)
                    isInK = (k>=Grid%ks) .and. (k<=Grid%ke)
                    isAct = isInI .and. isInJ .and. isInK
                    if (isAct) then
                        cycle
                    else
                        !This is outside active, map to conjugate
                        call lfmIJKcc(Model,Grid,i,j,k,ip,jp,kp)
                        inExyz(i,j,k,:) = inExyz(ip,jp,kp,:)
                    endif
                enddo
            enddo
        enddo
    end subroutine WrapBallCC

    subroutine WrapBallEC(Model,Grid,inEijk)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp) ,  intent(inout) :: inEijk(1:PsiSh+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM)

        integer :: i,j,k,ip,jp,kp,d,dijk

        !Just lazily loop over everything
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,ip,jp,kp,d,dijk) 
        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                do i=1,PsiSh+1
                    do d=1,NDIM
                        call lfmIJKec(Model,Grid,d,i,j,k,ip,jp,kp)
                        !If any index is different then do something
                        dijk = abs(i-ip) + abs(j-jp) + abs(k-kp)
                        
                        if ( dijk > 0 ) then
                            inEijk(i,j,k,d) = inEijk(ip,jp,kp,d)
                        endif
                    enddo !d
                enddo !i
            enddo !j
        enddo !k
        
    end subroutine WrapBallEC
    !-----------------------------
    !Calculate currents on inner-most active radial shells
    !NOTE: Only using perturbation field as we are assuming B0 is curl-free in inner region
    !See gridutils/bFld2Jxyz for full commented calculation
    !TODO: Better wrap this in routines to avoid replicated code
    !NOTE: Assuming this is being called on Voltron w/ full J/K bounds 
    subroutine GetShellJ(Model,Grid,Bxyz,Jxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), intent(in)  :: Bxyz(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), intent(inout) :: Jxyz(1:JpSh,Grid%js:Grid%je,Grid%ks:Grid%ke,1:NDIM)

        real(rp), dimension(:,:,:,:), allocatable :: bInt,JdS
        integer :: iG,j,k,n,NumP
        real(rp), dimension(NDIM) :: dl,bEdg,xcc,JdV
        real(rp) :: fI,fJ,fK,fIp,fJp,fKp,dV,Div
        real(rp), dimension(NDIM) :: xfI,xfIp,xfJ,xfJp,xfK,xfKp

    !bInt = edge-integrated magnetic field
        allocate(bInt(1:JpSh+1,Grid%js:Grid%je+1,Grid%ks:Grid%ke+1,1:NDIM))
        bInt = 0.0
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(n,j,k,iG,dl,bEdg)
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do n=1,JpSh+1
                    iG = JpSt+n-1

                    !Calculate I-dir
                    dl = Grid%xyz(iG+1,j,k,:) - Grid%xyz(iG,j,k,:)
                    bEdg = 0.25*( Bxyz(iG  ,j  ,k  ,:) + Bxyz(iG  ,j-1,k  ,:) &
                                + Bxyz(iG  ,j  ,k-1,:) + Bxyz(iG  ,j-1,k-1,:) )
                    bInt(n,j,k,IDIR) = dot_product(bEdg,dl)  
                    
                    !Calculate J-dir
                    dl = Grid%xyz(iG,j+1,k,:) - Grid%xyz(iG,j,k,:)
                    bEdg = 0.25*( Bxyz(iG  ,j  ,k  ,:) + Bxyz(iG-1,j  ,k  ,:) &
                                + Bxyz(iG  ,j  ,k-1,:) + Bxyz(iG-1,j  ,k-1,:) )
                    bInt(n,j,k,JDIR) = dot_product(bEdg,dl)  

                    !Calculate K-dir
                    dl = Grid%xyz(iG,j,k+1,:) - Grid%xyz(iG,j,k,:)
                    bEdg = 0.25*( Bxyz(iG  ,j  ,k  ,:) + Bxyz(iG-1,j  ,k  ,:) &
                                + Bxyz(iG  ,j-1,k  ,:) + Bxyz(iG-1,j-1,k  ,:) )
                    bInt(n,j,k,KDIR) = dot_product(bEdg,dl)  

                enddo
            enddo
        enddo

        !Now enforce matching edges
        do n=1,JpSh+1
            !K-seam
            bInt(n,:,Grid%ke+1,IDIR) = bInt(n,:,Grid%ks,IDIR)
            bInt(n,:,Grid%ke+1,JDIR) = bInt(n,:,Grid%ks,JDIR)
            !Ring (js)
            bInt(n,Grid%js  ,Grid%ks:Grid%ke+1,KDIR) = sum(bInt(n,Grid%js  ,Grid%ks:Grid%ke,KDIR))/Grid%Nkp
            bInt(n,Grid%js  ,Grid%ks:Grid%ke+1,IDIR) = sum(bInt(n,Grid%js  ,Grid%ks:Grid%ke,IDIR))/Grid%Nkp

            !Ring (je+1)
            bInt(n,Grid%je+1,Grid%ks:Grid%ke+1,KDIR) = sum(bInt(n,Grid%je+1,Grid%ks:Grid%ke,KDIR))/Grid%Nkp
            bInt(n,Grid%je+1,Grid%ks:Grid%ke+1,IDIR) = sum(bInt(n,Grid%je+1,Grid%ks:Grid%ke,IDIR))/Grid%Nkp

        enddo

    !JdS = face-integrated current flux
        allocate(JdS(1:JpSh+1,Grid%js:Grid%je+1,Grid%ks:Grid%ke+1,1:NDIM))
        JdS = 0.0
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(n,j,k,iG)
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do n=1,JpSh+1
                    iG = JpSt+n-1
                    !Do i-face
                    if ( (k<=Grid%ke) .and. (j<=Grid%je) ) then
                        JdS(n,j,k,IDIR) = bInt(n,j,k,JDIR) + bInt(n,j+1,k,KDIR) - bInt(n,j,k+1,JDIR) - bInt(n,j,k,KDIR)
                    endif
                    !Do j-face
                    if ( (n<=JpSh) .and. (k<=Grid%ke) ) then
                        JdS(n,j,k,JDIR) = bInt(n,j,k,KDIR) + bInt(n,j,k+1,IDIR) - bInt(n+1,j,k,KDIR) - bInt(n,j,k,IDIR)
                    endif
                    !Do k-face
                    if ( (n<=JpSh) .and. (j<=Grid%je) ) then
                        JdS(n,j,k,KDIR) = bInt(n,j,k,IDIR) + bInt(n+1,j,k,JDIR) - bInt(n,j+1,k,IDIR) - bInt(n,j,k,JDIR)
                    endif
                enddo
            enddo
        enddo

        !Now handle ring stuff
        do n=1,JpSh+1
            !Zero out degenerate faces
            JdS(n,Grid%js  ,:,JDIR) = 0.0
            JdS(n,Grid%je+1,:,JDIR) = 0.0
            !Match k fluxes at seam
            JdS(n,:,Grid%ke+1,KDIR) = JdS(n,:,Grid%ks,KDIR)
        enddo

    !Now handle Jxyz
        Jxyz = 0.0
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(n,j,k,iG,dV,fI,fJ,fK,fIp,fJp,fKp) &
        !$OMP private(xfI,xfJ,xfK,xfIp,xfJp,xfKp,Div,xcc,JdV)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do n=1,JpSh
                    iG = JpSt+n-1

                    !Grab values
                    dV = Grid%volume(iG,j,k)
                    fI  = JdS(n  ,j  ,k  ,IDIR)
                    fJ  = JdS(n  ,j  ,k  ,JDIR)
                    fK  = JdS(n  ,j  ,k  ,KDIR)
                    fIp = JdS(n+1,j  ,k  ,IDIR)
                    fJp = JdS(n  ,j+1,k  ,JDIR)
                    fKp = JdS(n  ,j  ,k+1,KDIR)

                    xfI  = Grid%xfc(iG  ,j  ,k  ,:,IDIR)
                    xfJ  = Grid%xfc(iG  ,j  ,k  ,:,JDIR)
                    xfK  = Grid%xfc(iG  ,j  ,k  ,:,KDIR)
                    xfIp = Grid%xfc(iG+1,j  ,k  ,:,IDIR)
                    xfJp = Grid%xfc(iG  ,j+1,k  ,:,JDIR)
                    xfKp = Grid%xfc(iG  ,j  ,k+1,:,KDIR)

                    Div = fIp - fI + fJp - fJ + fKp - fK
                    xcc = Grid%xyzcc(iG,j,k,:)
                    JdV = (fIp*xfIp + fJp*xfJp + fKp*xfKp - fI*xfI - fJ*xfJ - fK*xfK) - Div*xcc

                    Jxyz(n,j,k,XDIR:ZDIR) = JdV/dV

                enddo
            enddo
        enddo

        !Do some smoothing on currents
        if (doRingFAC .and. Model%doRing) then
            
            NumP = Model%Ring%Np
            do n=1,JpSh
                if (Model%Ring%doS) then
                    call FixRAVec_S(Model,Grid,Jxyz(n,Grid%js:Grid%js+1,Grid%ks:Grid%ke,1:NDIM))
                endif

                if (Model%Ring%doE) then
                    call FixRAVec_E(Model,Grid,Jxyz(n,Grid%je-1:Grid%je,Grid%ks:Grid%ke,1:NDIM))
                endif

            enddo !JpSh loop

        endif !doRingFAC

    end subroutine GetShellJ

    !Pack from corners to 1D array
    subroutine SquishCorners(Q,Qs)
        real(rp), dimension(2,2,2), intent(in)  :: Q
        real(rp), dimension(8)    , intent(out) :: Qs

        integer :: k,j,i,n

        !Doing this tedious to avoid temporary arrays
        n = 1
        do k=1,2
            do j=1,2
                do i=1,2
                    Qs(n) = Q(i,j,k)
                    n = n+1
                enddo
            enddo
        enddo

    end subroutine SquishCorners
end module cmiutils
