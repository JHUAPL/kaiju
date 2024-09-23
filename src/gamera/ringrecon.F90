!-----
!Ring reconstruction routines
module ringrecon

    use gamtypes
    use math
    use ringutils
    use multifluid

    implicit none

    integer, parameter, private :: RingNg = 2 !Number of ghost zones for ring reconstruction
    
    !RingLR_T
    !Generic reconstruction routine for ring-avg
    !5-point stencil -> L/R
    abstract interface
        subroutine RingLR_T(fm2,fm1,f,fp1,fp2,vm,vp)
            Import :: rp
            real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
            real(rp), intent(out) :: vm,vp
        end subroutine RingLR_T
    end interface

    !Set choice of ring reconstruction here
    procedure(RingLR_T), pointer :: RingLR

    contains
    
    !Reconstruct chunked ring values
    !m indices are over ring
    !n indices are over chunks
    !isGO (optional), which cells are good
    subroutine ReconstructRing(Model,rW,Nc,isGO)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: rW(Model%Ring%Np)
        integer, intent(in) :: Nc !Number of chunks in ring
        logical, intent(in), optional :: isGO(Model%Ring%Np)

        !Hold expanded ring (ie ghosts)
        real(rp) :: chW(1-RingNg:Nc+RingNg)
        real(rp) :: fm2,fm1,f,fp1,fp2,fL,fR
        real(rp) :: a,b,c,fI
        integer :: dJ,m,n,nv, mS,mE,ngood
        logical :: isG(Model%Ring%Np)

        !DIR$ ASSUME_ALIGNED rW: ALIGN

        if (present(isGO)) then
            isG = isGO
        else
            isG = .true.
        endif

        associate(Np=>Model%Ring%Np)
        dJ = Np/Nc

        !Fill inner region and add ghosts
        chW(1:Nc) = rW(1:Np:dJ)
        do n=1,RingNg
            chW(1-n)  = chW(Nc-n+1)
            chW(Nc+n) = chW(n)
        enddo

        !Loop over each chunk, create interpolant for chunk
        do n=1,Nc
            !Indices of this chunk within 1,Np array
            mS = 1 + (n-1)*dJ
            mE = mS + dJ - 1
            if (Model%doMultiF) then
                ngood = count( isG(mS:mE) ) !Number of good elements in chunk
                if (ngood == 0) cycle !Nothing to do
            endif

            !Grab stencil for interval LR's
            fm2 = chW(n-2)
            fm1 = chW(n-1)
            f   = chW(n  )
            fp1 = chW(n+1)
            fp2 = chW(n+2)

            call RingLR(fm2,fm1,f,fp1,fp2,fL,fR)
            
            !Calculate coefficients for parabolic interpolant
            a = 3*(fL + fR - 2*f)
            b = 2*(3*f - fR - 2*fL)
            c = fL

            !Reconstruct within this chunk
            do m=1,dJ
                fI = (a/3.0)*(3*m*m-3*m+1)/(dJ*dJ) + 0.5*b*(2*m-1)/dJ + c
                rW(mS+m-1) = fI
            enddo
        enddo !Chunks

        end associate

    end subroutine ReconstructRing

    subroutine WgtReconstructRing(Model,rW,xW,Nc,isGO)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: rW(Model%Ring%Np)
        real(rp), intent(in) :: xW(Model%Ring%Np)
        integer, intent(in) :: Nc !Number of chunks in ring
        logical, intent(in), optional :: isGO(Model%Ring%Np)

        !Hold expanded ring (ie ghosts)
        real(rp) :: chW(1-RingNg:Nc+RingNg),tMass(1-RingNg:Nc+RingNg)
        
        real(rp) :: fm2,fm1,f,fp1,fp2,fL,fR
        real(rp) :: a,b,c,fI
        integer :: dJ,m,n,nv, mS,mE,ngood
        logical :: isG(Model%Ring%Np)
        real(rp) :: dwL,dwR,dwC,dwM,min1,min2
        real(rp), allocatable, dimension(:) :: xWC,Mxc,Mxi,dMx
        real(rp) :: xWC0

        !DIR$ ASSUME_ALIGNED rW: ALIGN

        if (present(isGO)) then
            isG = isGO
        else
            isG = .true.
        endif

        associate(Np=>Model%Ring%Np)
        dJ = Np/Nc

        allocate(Mxc(1:dJ),dMx(1:dJ),xWC(1:dJ),Mxi(1:dJ+1))

        !Loop over each chunk and get total mass
        do n=1,Nc
            !Indices of this chunk within 1,Np array
            mS = 1 + (n-1)*dJ
            mE = mS + dJ - 1

            tMass(n) = sum(xW(mS:mE))
        enddo

        !Fill inner region and add ghosts
        chW(1:Nc) = rW(1:Np:dJ)
        do n=1,RingNg
            chW(1-n)  = chW(Nc-n+1)
            chW(Nc+n) = chW(n)

            tMass(1-n) = tMass(Nc-n+1)
            tMass(Nc+n)= tMass(n)
        enddo

        do n=1-RingNg,Nc+RingNg
            chW(n) = chW(n)/(tMass(n)/dJ)
        enddo

        !Loop over each chunk, create interpolant for chunk
        do n=1,Nc
            !Indices of this chunk within 1,Np array
            mS = 1 + (n-1)*dJ
            mE = mS + dJ - 1
            if (Model%doMultiF) then
                ngood = count( isG(mS:mE) ) !Number of good elements in chunk
                if (ngood == 0) cycle
            endif

            !Grab stencil for interval LR's
            fm1 = chW(n-1)
            f   = chW(n  )
            fp1 = chW(n+1)
            
            dwL = f - fm1
            dwR = fp1 - f
            dwC = ( fp1 - fm1 )/2.0

            min1 = min(2*abs(dwL),2*abs(dwR))
            min2 = min(min1,abs(dwC))
            !SIGN(A,B) returns the value of A with the sign of B
            dwM = sign(min2,dwC)

            !Create cell-centered mass-coordinates and Jacobian (dMx/dJ)
            xWC = xW(mS:mE)
            xWC0 = sum(xWC)

            Mxi(1) = 0.0
            do m=2,dJ+1
                Mxi(m) = sum(xWC(1:m-1))/xWC0
            enddo

            do m=1,dJ
                Mxc(m) = 0.5*(Mxi(m)+Mxi(m+1))
                dMx(m) = Mxi(m+1)-Mxi(m)
            enddo

            !Reconstruct within this chunk
            do m=1,dJ
                fI = f + dwM*(Mxc(m)-0.5)*dMx(m)
                rW(mS+m-1) = fI*xW(mS+m-1)
            enddo
        enddo !Chunks
        end associate
    end subroutine WgtReconstructRing

    !Lazy routine to make things equivalent to PCM
    !Both L & R are f
    subroutine PCM_IntLR(fm2,fm1,f,fp1,fp2,vm,vp)
        real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
        real(rp), intent(out) :: vm,vp

        vm = f
        vp = f

    end subroutine PCM_IntLR

    !Lazy routine to make things equivalent to PLM
    !Need that vp+vm=2*f to zero out quadratic term

    subroutine PLM_IntLR(fm2,fm1,f,fp1,fp2,vm,vp)
        real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
        real(rp), intent(out) :: vm,vp

        real(rp) :: dwL,dwR,dwC,dwM
        real(rp) :: min1,min2

        dwL = f - fm1
        dwR = fp1 - f
        dwC = ( fp1 - fm1 )/2.0

        min1 = min(2*abs(dwL),2*abs(dwR))
        min2 = min(min1,abs(dwC))
        !SIGN(A,B) returns the value of A with the sign of B
        dwM = sign(min2,dwC)

        vp = f + 0.5*dwM
        vm = f - 0.5*dwM

    end subroutine PLM_IntLR

    subroutine WENO_IntLR(fm2,fm1,f,fp1,fp2,vm,vp)
        real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
        real(rp), intent(out) :: vm,vp

        vm = Weno(fp2,fp1,f,fm1,fm2)
        vp = Weno(fm2,fm1,f,fp1,fp2)

    end subroutine WENO_IntLR

    !Calculates PPM *interval* L/R
    subroutine PPM_IntLR(fm2,fm1,f,fp1,fp2,vm,vp)
        real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
        real(rp), intent(out) :: vm,vp

        real(rp) :: dvm1,dvm2,dvp1,dvp2
        real(rp) :: SM,Sm1,Sp1,Sp2
        real(rp) :: dvc,am,ap

        dvm2 = fm1 - fm2
        dvm1 = f   - fm1
        dvp1 = fp1 - f
        dvp2 = fp2 - fp1
 
        dvc = 0.5*(dvm2 + dvm1)
        SM  = 2.0*minmod(dvm2, dvm1)
        Sm1 = minmod(dvc, SM)
 
        dvc = 0.5*(dvm1 + dvp1)
        SM  = 2.0*minmod(dvm1, dvp1)
        Sp1 = minmod(dvc, SM)
 
        dvc = 0.5*(dvp2 + dvp1)
        SM  = 2.0*minmod(dvp2, dvp1)
        Sp2 = minmod(dvc, SM)
 
        vp = 0.5*(f + fp1) - (Sp2 - Sp1)/6.0
        vm = 0.5*(f + fm1) - (Sp1 - Sm1)/6.0
 
        ap = vp - f
        am = vm - f
 
        if (ap*am >= 0.0) then
            ap = 0.0
            am = 0.0
        else 
            if (abs(ap) >= 2.0*abs(am)) then
                ap = -2.0*am
            endif

            if (abs(am) >= 2.0*abs(ap)) then
                am = -2.0*ap
            endif
        endif
 
        vp = f + ap
        vm = f + am

    end subroutine PPM_IntLR

    !WENO-5 reconstruction
    function Weno(hm2,hm1,h,hp1,hp2) result(hI)
        real(rp), intent(in) :: hm2,hm1,h,hp1,hp2
        real(rp) :: hI

        real(rp) :: h0,h1,h2, is0,is1,is2
        real(rp) :: a0,a1,a2, w0,w1,w2

        h0 =  (1.0/3)*hm2 - (7.0/6)*hm1 + (11.0/6)*h
        h1 = -(1.0/6)*hm1 + (5.0/6)*h   +  (1.0/3)*hp1
        h2 =  (1.0/3)*h   + (5.0/6)*hp1 -  (1.0/6)*hp2

        is0 = (13.0/12)*(hm2 - 2.0*hm1 +   h)**2.0 + (1.0/4)*(hm2 - 4.0*hm1 + 3.0*h)**2.0
        is1 = (13.0/12)*(hm1 - 2.0*h   + hp1)**2.0 + (1.0/4)*(hm1 - hp1)**2.0
        is2 = (13.0/12)*(h   - 2.0*hp1 + hp2)**2.0 + (1.0/4)*(3.0*h - 4.0*hp1 + hp2)**2.0

        a0 = (1.0/10)/(is0+TINY)**2.0
        a1 = (6.0/10)/(is1+TINY)**2.0
        a2 = (3.0/10)/(is2+TINY)**2.0

        w0 = a0/(a0+a1+a2)
        w1 = a1/(a0+a1+a2)
        w2 = a2/(a0+a1+a2)
 
        hI = w0*h0 + w1*h1 + w2*h2
    end function Weno
end module ringrecon