module mgmres

    use kdefs
    
    implicit none

    real(rp), private, parameter :: delta = 1.0e-3
    logical , private, parameter :: doWarn = .true.
    contains

!---------------
    !A*x (compressed row)
    subroutine ax_cr(n,nz_num,ia,ja,A,x,W)
        integer , intent(in)  :: n,nz_num
        integer , intent(in)  :: ia(n+1),ja(nz_num)
        real(rp), intent(in)  :: A(nz_num)
        real(rp), intent(in)  :: x(n)
        real(rp), intent(inout) :: W(n)

        integer :: i,k1,k2

        W(1:n) = 0.0
        !TODO: Add OMP bindings here
        do i=1,n
            k1 = ia(i)
            k2 = ia(i+1) - 1
            W(i) = W(i) + dot_product( A(k1:k2), x(ja(k1:k2)) )
        enddo
    end subroutine ax_cr
!---------------
    subroutine diagonal_pointer_cr(n,nz_num,ia,ja,ua)
        integer, intent(in) :: n,nz_num
        integer, intent(in) :: ia(n+1),ja(nz_num)
        integer, intent(out) :: ua(n)

        integer :: i,k

        ua(1:n) = -1
        do i=1,n
            do k=ia(i),ia(i+1)-1
                if ( ja(k) == i ) then
                    ua(i) = k
                endif
            enddo
        enddo

    end subroutine diagonal_pointer_cr

!---------------
    !Incomplete LU factorization
    subroutine ilu_cr(n,nz_num,ia,ja,A,ua,L)
        integer, intent(in) :: n,nz_num
        integer, intent(in) :: ia(n+1),ja(nz_num)
        integer, intent(inout) :: ua(n)
        real(rp), intent(in) :: A(nz_num)
        real(rp), intent(out) :: L(nz_num)

        integer :: i,j,jj,jrow,jw,k
        integer :: iw(n)
        real(rp) :: tl

        L(1:nz_num) = A(1:nz_num)
        do i=1,n
            iw(1:n) = -1
            do k=ia(i),ia(i+1)-1
                iw(ja(k)) = k
            enddo

            do j=ia(i),ia(i+1)-1
                jrow = ja(j)
                if ( i<= jrow ) exit
                tl = L(j)*L(ua(jrow))
                L(j) = tl
                do jj=ua(jrow)+1,ia(jrow+1)-1
                    jw = iw(ja(jj))
                    if ( jw /= -1) then
                        L(jw) = L(jw) - tl*L(jj)
                    endif
                enddo
            enddo

            ua(i) = j
            L(j) = 1.0/L(j)
        enddo
        L(ua(1:n)) = 1.0/L(ua(1:n))
    end subroutine ilu_cr

!---------------
    subroutine lus_cr(n,nz_num,ia,ja,L,ua,R,Z)
        integer, intent(in) :: n,nz_num
        integer, intent(in) :: ia(n+1),ja(nz_num),ua(n)
        real(rp), intent(inout) :: L(nz_num)
        real(rp), intent(in) :: R(n)
        real(rp), intent(out) :: Z(n)

        integer :: i,j
        real(rp), dimension(n) :: W

        !Copy in
        W(1:n) = R(1:n)

        !Solve L*W=W, L unit lower triangular
        !TODO: Add OMP bindings here
        do i=2,n
            do j=ia(i),ua(i)-1
                W(i) = W(i) - L(j)*W(ja(j))
            enddo
        enddo

        !Solve U*W=W, U upper triangular
        do i=n,1,-1
            do j=ua(i)+1,ia(i+1)-1
                W(i) = W(i) - L(j)*W(ja(j))
            enddo
            W(i) = W(i)/L(ua(i))
        enddo

        !Copy out
        Z(1:n) = W(1:n)

    end subroutine lus_cr

!---------------
    !Givens rotation
    subroutine mult_givens(c,s,k,g)
        real(rp), intent(in) :: c,s !cos/sin of rotation
        integer, intent(in) :: k
        real(rp), intent(inout) :: g(1:k+1)

        real(rp) :: g1,g2
        g1 = c*g(k) - s*g(k+1)
        g2 = s*g(k) + c*g(k+1)

        g(k  ) = g1
        g(k+1) = g2

    end subroutine mult_givens

!---------------
    !Preconditioned restarted GMRES
    subroutine pmgmres_ilu_cr(n,nz_num,ia,ja,A,x,rhs,itr_max,mr,tol_abs,tol_rel, isGoodSolveO)
        integer, intent(in) :: n,nz_num
        integer, intent(inout) :: ia(n+1),ja(nz_num)
        real(rp), intent(inout) :: A(nz_num)
        real(rp), intent(inout) :: x(n)
        real(rp), intent(in) :: rhs(n)
        integer, intent(in) :: itr_max,mr
        real(rp), intent(in) :: tol_abs,tol_rel
        logical, optional, intent(inout) :: isGoodSolveO

        real(rp) :: av,htmp,mu,rho,rho_tol
        real(rp) :: r(n)
        integer :: ua(n)
        real(rp), dimension(mr+1) :: c,g,s,y
        real(rp) :: H(mr+1,mr)
        real(rp) :: L(ia(n+1)+1)
        real(rp) :: V(n,mr+1)
        integer :: i,itr,itr_used,j,k,k_copy

        logical :: isConverged

        itr_used = 0
        call rearrange_cr(n,nz_num,ia,ja,A)
        call diagonal_pointer_cr(n,nz_num,ia,ja,ua)
        call ilu_cr(n,nz_num,ia,ja,A,ua,L)

        do itr=1,itr_max
            call ax_cr(n,nz_num,ia,ja,A,x,r)
            r(1:n) = rhs(1:n) - r(1:n)
            call lus_cr(n,nz_num,ia,ja,L,ua,r,r)
            rho = sqrt(dot_product(r,r))

            if (itr == 1) then
                rho_tol = rho*tol_rel
            endif

            V(1:n,1) = r(1:n)/rho
            g(1) = rho
            g(2:mr+1) = 0.0
            H(1:mr+1,1:mr) = 0.0

            do k=1,mr
                k_copy = k
                call ax_cr (n,nz_num,ia,ja,A,V(1:n,k),V(1:n,k+1) )
                call lus_cr(n,nz_num,ia,ja,L,ua,V(1:n,k+1),V(1:n,k+1) )
                av = sqrt( dot_product( V(1:n,k+1), V(1:n,k+1) ) )
                do j=1,k
                    H(j,k) = dot_product ( V(1:n,k+1), V(1:n,j) )
                    V(1:n,k+1) = V(1:n,k+1) - V(1:n,j) * H(j,k)
                enddo

                H(k+1,k) = sqrt ( dot_product ( V(1:n,k+1), V(1:n,k+1) ) )

                if ( ( av + delta * H(k+1,k)) == av ) then
                    do j=1,k
                        htmp = dot_product ( V(1:n,k+1), V(1:n,j) )
                        H(j,k) = H(j,k) + htmp
                        V(1:n,k+1) = V(1:n,k+1) - htmp * V(1:n,j)
                    enddo
                    H(k+1,k) = sqrt ( dot_product ( V(1:n,k+1), V(1:n,k+1) ) )
                endif

                if ( abs(H(k+1,k)) >= TINY ) then
                    V(1:n,k+1) = V(1:n,k+1)/H(k+1,k)
                endif

                if ( 1<k ) then
                    y(1:k+1) = H(1:k+1,k)
                    do j=1,k-1
                        call mult_givens(c(j),s(j),j,y)
                    enddo
                    H(1:k+1,k) = y(1:k+1)
                endif

                mu = sqrt( H(k,k)**2.0 + H(k+1,k)**2.0)
                c(k) =  H(k  ,k)/mu
                s(k) = -H(k+1,k)/mu
                H(k,k) = c(k)*H(k,k) - s(k)*H(k+1,k)
                H(k+1,k) = 0.0
                call mult_givens(c(k),s(k),k,g)
                rho = abs(g(k+1))

                itr_used = itr_used+1
                if ( rho <= rho_tol .and. rho <= tol_abs ) then
                    !Exit if done
                    exit
                endif
            enddo !K-loop
            k = k_copy - 1
            y(k+1) = g(k+1)/H(k+1,k+1)

            do i=k,1,-1
                y(i) = ( g(i) - dot_product ( H(i,i+1:k+1), y(i+1:k+1) ) ) / H(i,i)
            enddo

            do i=1,n
                x(i) = x(i) + dot_product ( V(i,1:k+1), y(1:k+1) )
            enddo

            if ( (rho<=rho_tol) .and. (rho<=tol_abs) ) then
                exit
            endif


        enddo !Iteration loop

        !Check convergence
        isConverged = (rho<=rho_tol) .and. (rho<=tol_abs)

        if ( doWarn .and. (.not. isConverged) ) then
            write(*,*) '---'
            write(*,*) 'WARNING: PMGMRES failed to converge'
            write(*,*) '  Iterations        = ', itr, itr_used, itr_max
            write(*,*) '  Final Residual    = ', rho
            write(*,*) '  Rel/Abs Tolerance = ', rho_tol,tol_abs
            
            write(*,*) '---'

            if (present(isGoodSolveO)) isGoodSolveO = .false.

        endif

        if (present(isGoodSolveO)) isGoodSolveO = .true.

    end subroutine pmgmres_ilu_cr

!---------------
    !Sort a compressed row matrix
    subroutine rearrange_cr(n,nz_num,ia,ja,A)
        integer, intent(in) :: n,nz_num
        integer, intent(in) :: ia(n+1)
        integer, intent(inout) :: ja(nz_num)
        real(rp), intent(inout) :: A(nz_num)

        integer :: i,itemp,k,l
        real(rp) :: rtemp

        do i=1,n
            do k=ia(i),ia(i+1)-2
                do l=k+1,ia(i+1)-1
                    if ( ja(l) < ja(k) ) then
                        itemp = ja(l)
                        ja(l) = ja(k)
                        ja(k) = itemp

                        rtemp = A(l)
                        A(l) = A(k)
                        A(k) = rtemp
                    endif
                enddo
            enddo
        enddo

    end subroutine rearrange_cr

end module mgmres
