! sets up and solves the stencil matrix
#ifdef USEMKL
  include 'mkl_pardiso.f90'
#endif
module mixsolver
  use mixdefs
  use mixtypes
  use mgmres
  use clocks
#ifdef USEMKL
  use mkl_pardiso
#endif

  implicit none

  logical, private :: isSolverInit=.false.

  contains
    ! running index
    integer(kind=8) function Kindex(j,i,Nj)
      integer, intent(in) :: i,j,Nj
      Kindex=(i-1)*Nj+j
    end function Kindex

    ! map 2 to -1, -2 to 1, and [-1,0,1] to themselves
    ! wt stands for wrap tripls
    integer function wt(q)
      integer,intent(in) :: q
      wt = nint(asin(sin(q*2*pi/3))/asin(sin(2*pi/3)))
    end function wt

    subroutine init_solver(P,G,S)
      type(mixParams_T), intent(in) :: P
      type(mixGrid_T), intent(in) :: G
      type(Solver_T), intent(inout) :: S
      integer :: i,j

      ! this is the number of non-zeros in the matrix. Note, this
      ! depends on the stencil and boundary conditions
      S%nnz = G%Np*(G%Nt-2)*5 + G%Np + G%Np*(G%Np+1) 

      ! allocate RHS vector
      if (.not.allocated(S%RHS)) allocate(S%RHS(G%Nt*G%Np))
      ! allocate row index
      if (.not.allocated(S%rowI)) allocate(S%rowI(G%Nt*G%Np+1))
      ! Matrix 
      if (.not.allocated(S%data)) allocate(S%data(S%nnz))
      if (.not.allocated(S%II)) allocate(S%II(S%nnz))
      if (.not.allocated(S%JJ)) allocate(S%JJ(S%nnz))
      ! Equation terms
      if (.not.allocated(S%F11)) allocate(S%F11(G%Np,G%Nt))
      if (.not.allocated(S%F22)) allocate(S%F22(G%Np,G%Nt))
      if (.not.allocated(S%F12)) allocate(S%F12(G%Np,G%Nt))

      ! low lat boudnary array
      if (.not.allocated(S%LLBC)) allocate(S%LLBC(G%Np))
      ! solution storate
      if (.not. allocated(S%solution)) allocate(S%solution(G%Np*G%Nt))
      ! assign 0s for solver initial guess
      S%solution = 0.0_rp

    end subroutine init_solver

    subroutine set_solver_terms(P,G,St,S)
      type(mixState_T), intent(in) :: St
      type(mixParams_T), intent(in) :: P  ! passing parameters just in case. Not used.
      type(mixGrid_T), intent(in) :: G
      type(Solver_T), intent(inout) :: S
      
      S%F11 = sin(G%t)*St%Vars(:,:,SIGMAP)/G%cosd**2
      S%F22 = St%Vars(:,:,SIGMAP)
      S%F12 =-St%Vars(:,:,SIGMAH)/G%cosd ! (assuming F21=-F12)
    end subroutine set_solver_terms

    subroutine set_solver_matrix_and_rhs(P,G,St,S)
      type(mixState_T), intent(in) :: St
      type(mixParams_T), intent(in) :: P  ! passing parameters just in case. Not used.
      type(mixGrid_T), intent(in) :: G
      type(Solver_T), intent(inout) :: S
      integer :: i,j,jm1,jp1,jj,a1,a2,q,u
      integer(kind=8) :: count,nextRow,nextRowI
      real(rp) :: dF12t,dF12p
      real(rp) :: d(-2:2),c(-2:2)

      ! init to zero. This is important, since we're only filling in
      ! non-zero elements in the matrix construction
      S%RHS = 0.0D0 
      S%data=0.0D0
      S%II = 0.0D0
      S%JJ = 0.0D0
      
      count=1
      nextRowI=1
      ! make a note about the fact that the order is corrects for the CSR format

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Pole boundary 
      do j=1,G%Np
         S%data(count) = 1.0D0
         S%II(count)   = Kindex(j,1,G%Np)
         S%JJ(count)   = Kindex(j,1,G%Np)
         count=count+1

         do jj=1,G%Np   ! with (Np,Nt) definition of the grid, this is a natural alignment
            S%data(count) = -G%dp(jj,2)/(2*pi)
            S%II(count) = Kindex(j,1,G%Np)
            S%JJ(count) = Kindex(jj,2,G%Np)
            count=count+1
         enddo
         
         ! These are points on the pole boundary Thus, each row of the
         ! matrix has (Np+1) entries: for the point itself + Np
         ! entries for averaging the points at i=2 (see above).
         ! Therefore, the rows start with indices 1, Np+2, 2Np+3, etc.
         S%rowI(Kindex(j,1,G%Np)) = nextRowI       !1+(j-1)*(G%Np+1)
         nextRowI=nextRowI+G%Np+1
      enddo
      ! note, not setting RHS because it's initializaed to zero anyway
      ! end pole boundary
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! inner block
      do i=2,G%Nt-1  ! excluding pole and low lat boundaries
         do j=1,G%Np
            jm1 = merge(G%Np,j-1,j.eq.1)   ! maps j-1=0 to j-1=Np, otherwise returns j-1
            jp1 = merge(1,j+1,j.eq.G%Np)   ! similar for j+1

            ! the above functions involve if statements, which I don't
            ! want to pack inside inner loop. Do this instead: these
            ! wonderful functions of my invention maps j-1=0 to Np,
            ! otherwise j-1 maps to itself. Similarly, j+1=Np is
            ! mapped to j+1=1, otherwise nothing's done
!            jm1 = modulo(j-1,G%Np)+G%Np*(1-int(ceiling(real(j-1)/G%Np)))  
!            jp1 = modulo(j+1,G%Np)+G%Np*(1-int(ceiling(real(modulo(j+1,G%Np))/G%Np)))

            ! derivatives for off diagonal conductance terms
            dF12p = G%fp(j,i)*( G%dp(jm1,i)/G%dp(j,i)*S%F12(jp1,i) + G%dpdp(j,i)*S%F12(j,i) - G%dp(j,i)/G%dp(jm1,i)*S%F12(jm1,i) )
            dF12t = G%ft(j,i)*( G%dt(j,i-1)/G%dt(j,i)*S%F12(j,i+1) + G%dtdt(j,i)*S%F12(j,i) - G%dt(j,i)/G%dt(j,i-1)*S%F12(j,i-1) )

            ! for periodic boundaries jm1>jp1 maybe the case, which
            ! breaks the order of things. So, let's figure it out

            a1 = (jp1-j)/abs(jp1-j)
            a2 = (jm1-j)/abs(jm1-j)
            q  = nint(0.5*(a1+a2))
            
            d(-2) = G%ft(j,i)*(S%F11(j,i)+S%F11(j,i-1))/G%dt(j,i-1)+&
                 dF12p*G%ft(j,i)*G%dt(j,i)/G%dt(j,i-1)
            c(-2) = Kindex(j,i-1,G%Np)

            d(wt(-q-1)) = G%fp(j,i)/sin(G%t(j,i))**2*(S%F22(j,i)+S%F22(jm1,i))/G%dp(jm1,i)-&
                 dF12t*G%fp(j,i)*G%dp(j,i)/G%dp(jm1,i)
            c(wt(-q-1)) = Kindex(jm1,i,G%Np)

            d(-q) = -G%ft(j,i)*( (S%F11(j,i)+S%F11(j,i+1))/G%dt(j,i)+(S%F11(j,i)+S%F11(j,i-1))/G%dt(j,i-1) ) - &
                 G%fp(j,i)/sin(G%t(j,i))**2*( (S%F22(j,i)+S%F22(jp1,i))/G%dp(j,i)+(S%F22(j,i)+S%F22(jm1,i))/G%dp(jm1,i) ) + &
                 dF12t*G%fp(j,i)*G%dpdp(j,i)-&
                 dF12p*G%ft(j,i)*G%dtdt(j,i) 
            c(-q) = Kindex(j,i,G%Np)

            d(wt(-q+1)) = G%fp(j,i)/sin(G%t(j,i))**2*(S%F22(j,i)+S%F22(jp1,i))/G%dp(j,i)+&
                 dF12t*G%fp(j,i)*G%dp(jm1,i)/G%dp(j,i)
            c(wt(-q+1)) = Kindex(jp1,i,G%Np)
            
            d(2) = G%ft(j,i)*(S%F11(j,i)+S%F11(j,i+1))/G%dt(j,i)-&
                 dF12p*G%ft(j,i)*G%dt(j,i-1)/G%dt(j,i)
            c(2) = Kindex(j,i+1,G%Np)

            do u=-2,2
               s%data(count) = d(u)
               S%JJ(count)   = c(u)
               S%II(count)   = Kindex(j,i,G%Np)
               count=count+1
            enddo

            ! At this point, coming out of the pole boundary condition, we have: S%rowI(Kindex(j,i,G%Np)-1) = G%Np*(G%Np+1)+1
!            S%rowI(Kindex(j,i,G%Np)) = S%rowI(Kindex(j,i,G%Np)-1)+G%Np+1 + 5*(Kindex(j,i,G%Np)-1)
!            S%rowI(Kindex(j,i,G%Np)) = G%Np*(G%Np+1)+1 + 5*(Kindex(j,i,G%Np)-1)
            S%rowI(Kindex(j,i,G%Np)) = nextRowI
            nextRowI = nextRowI+5
            S%RHS(Kindex(j,i,G%Np)) = St%Vars(j,i,FAC)*G%cosd(j,i)
         enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! low lat boundary
      do j=1,G%Np
         S%data(count)  = 1.0D0
         S%II(count)    = Kindex(j,G%Nt,G%Np)
         S%JJ(count)    = Kindex(j,G%Nt,G%Np)
         count=count+1

         S%rowI(Kindex(j,G%Nt,G%Np)) = nextRowI 
         nextRowI=nextRowI+1
         S%RHS(Kindex(j,G%Nt,G%Np)) = S%LLBC(j)
      enddo
      S%rowI(Kindex(G%Np,G%Nt,G%Np)+1) = nextRowI   ! the dummy last element as required by the CSR storage
    end subroutine set_solver_matrix_and_rhs
    ! end low lat boundary
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine run_solver(P,G,St,S,isGoodSolveO)
      type(mixParams_T), intent(in) :: P
      type(mixGrid_T),   intent(in) :: G
      type(mixState_T),  intent(in) :: St
      type(Solver_T),    intent(inout) :: S
      logical, optional, intent(inout) :: isGoodSolveO

      ! low latitude boundary condition
      S%LLBC = P%llbc_value

      ! this uses conductdance values to set up the Poisson equation
      call set_solver_terms(P,G,St,S)
      ! this computes the matrix coefs and RHS
      call set_solver_matrix_and_rhs(P,G,St,S)

      ! reset initial solution guess to 0 each time
      S%solution = 0.0_rp

      ! actual solve
      ! note initial condition is set to 0.
      ! there is a version of this in mix.F90 that stores the old solution
      ! and initializes the next one from that but I found that doesn't
      ! affect the convergence speed
      if(norm2(S%RHS) < 1e-6_rp .and. norm2(S%solution) < 1e-6_rp) then
          ! dealing with the edge case of guessing 0 for a rhs of 0
          S%solution = 0.01_rp
      endif

#ifdef USEMKL
      !Call MKL-Pardiso solver
      call MKLSolve(P,G,St,S,isGoodSolveO=isGoodSolveO)
#else
      ! could drop tolerances by a factor of 100 each (from default 1.d-6)
      ! -- still same result
      call pmgmres_ilu_cr (&
           G%Np*G%Nt,&
           S%nnz,&
           S%rowI,&
           S%JJ,&
           S%data,&
           S%solution,&
           S%RHS,&
           P%maxitr,&
           P%mr,&
           P%tol_abs,&
           P%tol_rel,isGoodSolveO=isGoodSolveO)
#endif
    end subroutine run_solver

#ifdef USEMKL
    subroutine MKLSolve(P,G,St,S,isGoodSolveO)
      type(mixParams_T), intent(in) :: P
      type(mixGrid_T),   intent(in) :: G
      type(mixState_T),  intent(in) :: St
      type(Solver_T),    intent(inout) :: S
      logical, optional, intent(inout) :: isGoodSolveO

      integer :: Npt
      TYPE(MKL_PARDISO_HANDLE) :: pt(64)
      integer :: error,nrhs,msglvl,mtype,iparm(64)
      integer, dimension(:), allocatable :: perm
      integer :: maxfct,mnum,phase !pardiso control params
      ! see description of parameters here: https://software.intel.com/en-us/node/470284#E44B4021-701A-48DA-BA29-70CFA20766AA

      if (P%mklmsglvl > 0) then
        write(*,*) 'Running MKL-Pardiso'
      endif

      mtype = 11 ! real nonsymmetric matrix

      call Tic("MKL-Init")
      !Always initialize/tear-down
      call pardisoinit(pt,mtype,iparm)

      Npt = G%Np*G%Nt
      allocate(perm(Npt))

      iparm(1)  = 0  ! all default parameters
      iparm(27) = 0  ! Matrix checker
      !Options added by K
      iparm(2)  = 3
      iparm(24) = 1
      !iparm(28) =0; iparm (6) =0;
      
      maxfct = 1
      mnum   = 1
      nrhs   = 1
      msglvl = P%mklmsglvl  ! no verbosity
      call Toc("MKL-Init")

      phase = 13 !Analysis, numerical factorization, solve, iterative refinement
      call Tic("MKL-Solve")
      !Call solver
      call pardiso(pt,maxfct,mnum,mtype,phase,Npt, &
                   S%data,S%rowI,S%JJ,perm, &
                   nrhs,iparm,msglvl, &
                   S%RHS,S%solution,error)

      if (error < 0) then
        write(*,*) 'MKL-Pardiso error, you should deal with that'
        if (present(isGoodSolveO)) then
          isGoodSolveO = .false.
          return
        else
          stop
        endif
      endif
      call Toc("MKL-Solve")

      call Tic("MKL-Fin")
      !Release memory
      phase =-1  ! release
      call pardiso(pt,maxfct,mnum,mtype,phase,Npt, &
                   S%data,S%rowI,S%JJ,perm, &
                   nrhs,iparm,msglvl, &
                   S%RHS,S%solution,error)
      call Toc("MKL-Fin")

      if (present(isGoodSolveO)) isGoodSolveO = .true.

    end subroutine MKLSolve

#endif

end module mixsolver
