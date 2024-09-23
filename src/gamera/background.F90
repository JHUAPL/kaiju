!Definitions and routines for background field incorporation

module background

    use gamtypes
    use gamutils
    use math
    use gridutils
    
#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    contains


    !Adds background field data to the Grid data structure
    subroutine AddB0(Model,Grid,B0)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        procedure(VectorField_T), pointer, intent(in) :: B0

        integer :: i,j,k
        real(rp), dimension(NDIM) :: f0,f1,f2,f3,eInt,fInt,fInt2,fIntX
        real(rp), dimension(NDIM) :: e1,e2
        real(rp), dimension(8,NDIM) :: xyzC
        real(rp), dimension(:,:,:,:,:), allocatable :: faceStress !(i,j,k,IJKDIR,XYZDIR)
        real(rp) :: ijkB(NVAR)
        real(rp) :: MagP

        procedure(GasIC_T), pointer :: Wxyz !Lazy wrapper for volume integral

        
        if (.not. Model%doBackground .or. .not. Grid%doB0Init) then
            !Do nothing if incorrectly configured
            return
        endif
        

        !Create background field data structures
        allocate(Grid%fcB0 (Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM,NDIM))
        allocate(Grid%edgB0(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,2   ,NDIM))
        allocate(Grid%B0(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM))
        allocate(Grid%bFlux0(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM))
        allocate(Grid%dpB0(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM))
        allocate(faceStress(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM,NDIM))
        

        Grid%fcB0 = 0.0
        Grid%edgB0 = 0.0
        Grid%B0 = 0.0
        Grid%dpB0 = 0.0
        Grid%bFlux0 = 0.0

        Wxyz => B0Wrap
        
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(f0,f1,f2,f3,e1,e2,eInt,xyzC,ijkB,fInt,fInt2,fIntX)
        do k=Grid%ksg, Grid%keg
           do j=Grid%jsg, Grid%jeg
               do i=Grid%isg, Grid%ieg
                    !Calculate face fields
                    !Need specific ordering of corner points

                !I face
                    call faceCoords(Model,Grid,i,j,k,IDIR,f0,f1,f2,f3)
                    call GaussianFaceFlux(Model,f0,f1,f2,f3,B0,Grid%bFlux0(i,j,k,IDIR))

                    call GaussianFaceIntegral(Model,f0,f1,f2,f3,B0,fInt,fInt2,fIntX)
                    Grid%fcB0(i,j,k,:,IDIR) = fInt/Grid%Face(i,j,k,IDIR)

                    !Get face stress
                    call GaussianFaceStress(Model,f0,f1,f2,f3,B0,faceStress(i,j,k,IDIR,:))

                !J face
                    call faceCoords(Model,Grid,i,j,k,JDIR,f0,f1,f2,f3)
                    call GaussianFaceFlux(Model,f0,f1,f2,f3,B0,Grid%bFlux0(i,j,k,JDIR))

                    call GaussianFaceIntegral(Model,f0,f1,f2,f3,B0,fInt,fInt2,fIntX)
                    Grid%fcB0(i,j,k,:,JDIR) = fInt/Grid%Face(i,j,k,JDIR)

                    !Get face stress
                    call GaussianFaceStress(Model,f0,f1,f2,f3,B0,faceStress(i,j,k,JDIR,:))

                !K face
                    call faceCoords(Model,Grid,i,j,k,KDIR,f0,f1,f2,f3)
                    call GaussianFaceFlux(Model,f0,f1,f2,f3,B0,Grid%bFlux0(i,j,k,KDIR))
                    call GaussianFaceIntegral(Model,f0,f1,f2,f3,B0,fInt,fInt2,fIntX)
                    Grid%fcB0(i,j,k,:,KDIR) = fInt/Grid%Face(i,j,k,KDIR)

                    !Get face stress
                    call GaussianFaceStress(Model,f0,f1,f2,f3,B0,faceStress(i,j,k,KDIR,:))

                    !Calculate edge integrals and mapping to 1/2 system using velocity mapping
                    !I edge
                    call edgeCoords(Model,Grid,i,j,k,IDIR,e1,e2)
                    
                    eInt = GaussianEdgeIntegral(Model,e1,e2,B0)

                    Grid%edgB0(i,j,k,1,IDIR) = eInt(XDIR)*Grid%Te(i,j,k,TAN1X,IDIR) + &
                                               eInt(YDIR)*Grid%Te(i,j,k,TAN1Y,IDIR) + &
                                               eInt(ZDIR)*Grid%Te(i,j,k,TAN1Z,IDIR) 
                    Grid%edgB0(i,j,k,2,IDIR) = eInt(XDIR)*Grid%Te(i,j,k,TAN2X,IDIR) + &
                                               eInt(YDIR)*Grid%Te(i,j,k,TAN2Y,IDIR) + &
                                               eInt(ZDIR)*Grid%Te(i,j,k,TAN2Z,IDIR) 
                    !J edge
                    call edgeCoords(Model,Grid,i,j,k,JDIR,e1,e2)
                    eInt = GaussianEdgeIntegral(Model,e1,e2,B0)

                    Grid%edgB0(i,j,k,1,JDIR) = eInt(XDIR)*Grid%Te(i,j,k,TAN1X,JDIR) + &
                                               eInt(YDIR)*Grid%Te(i,j,k,TAN1Y,JDIR) + &
                                               eInt(ZDIR)*Grid%Te(i,j,k,TAN1Z,JDIR) 
                    Grid%edgB0(i,j,k,2,JDIR) = eInt(XDIR)*Grid%Te(i,j,k,TAN2X,JDIR) + &
                                               eInt(YDIR)*Grid%Te(i,j,k,TAN2Y,JDIR) + &
                                               eInt(ZDIR)*Grid%Te(i,j,k,TAN2Z,JDIR) 

                    !K edge
                    call edgeCoords(Model,Grid,i,j,k,KDIR,e1,e2)
                    eInt = GaussianEdgeIntegral(Model,e1,e2,B0)

                    Grid%edgB0(i,j,k,1,KDIR) = eInt(XDIR)*Grid%Te(i,j,k,TAN1X,KDIR) + &
                                               eInt(YDIR)*Grid%Te(i,j,k,TAN1Y,KDIR) + &
                                               eInt(ZDIR)*Grid%Te(i,j,k,TAN1Z,KDIR) 
                    Grid%edgB0(i,j,k,2,KDIR) = eInt(XDIR)*Grid%Te(i,j,k,TAN2X,KDIR) + &
                                               eInt(YDIR)*Grid%Te(i,j,k,TAN2Y,KDIR) + &
                                               eInt(ZDIR)*Grid%Te(i,j,k,TAN2Z,KDIR) 

                    !Add cell center XYZ fields
                    !Get cell corners
                    call cellCoords(Model,Grid,i,j,k,xyzC)

                    ijkB = GaussianVolumeIntegral(xyzC,Wxyz)/Grid%volume(i,j,k)
                    Grid%B0(i,j,k,:) = ijkB(1:3)
                enddo
            enddo
        enddo       

        !Turn face stresses into volume integrated source term
        !$OMP PARALLEL DO default(shared) collapse(2) 
        do k=Grid%ks, Grid%ke
            do j=Grid%js, Grid%je
                do i=Grid%is, Grid%ie
                    Grid%dpB0(i,j,k,:) = -(faceStress(i+1,j,k,IDIR,:) &
                                          -faceStress(i  ,j,k,IDIR,:) &
                                          +faceStress(i,j+1,k,JDIR,:) &
                                          -faceStress(i,j  ,k,JDIR,:) &
                                          +faceStress(i,j,k+1,KDIR,:) &
                                          -faceStress(i,j,k  ,KDIR,:) )/Grid%volume(i,j,k)
               enddo
            enddo
        enddo
                
        Grid%doB0Init = .false. !Don't do again
        
        contains
            !Wrapper (to look like GasIC_T)
            subroutine B0Wrap(x,y,z,Bx,By,Bz,Aa,Ab)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: Bx,By,Bz,Aa,Ab

                Aa = 0.0
                Ab = 0.0
                call B0(x,y,z,Bx,By,Bz)
            end subroutine B0Wrap
    end subroutine AddB0

    !Calculates accelerations from gravitational potential on grid
    subroutine AddGrav(Model,Grid,Phi)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        procedure(ScalarFun_T), pointer, intent(in) :: Phi

        integer :: i,j,k
        real(rp), dimension(NDIM) :: nI,nJ,nK,nIp,nJp,nKp
        real(rp), dimension(NDIM) :: g,rHat
        real(rp) :: PhiI,PhiJ,PhiK,PhiIp,PhiJp,PhiKp
        real(rp) ::  daI, daJ, daK, daIp, daJp, daKp
        real(rp) :: dV
        if (.not. Model%doGrav .or. .not. Grid%doG0Init) then
            !Do nothing if incorrectly configured
            return
        endif

        !Create and zero out forces
        allocate(Grid%gxyz(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM))
        Grid%gxyz = 0.0


        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(nI,nJ,nK,nIp,nJp,nKp) &
        !$OMP private(PhiI,PhiJ,PhiK,PhiIp,PhiJp,PhiKp) &
        !$OMP private(daI,daJ,daK,daIp,daJp,daKp) &
        !$OMP private(dV,g,rHat)
        do k=Grid%ks, Grid%ke
            do j=Grid%js, Grid%je
                do i=Grid%is, Grid%ie
                    !Get volumes/face center points
                    dV = Grid%volume(i,j,k)

                    !Evaluate phi at each face center
                    call PWrap(Grid%xfc(i,j,k  ,:,IDIR),  PhiI  )
                    call PWrap(Grid%xfc(i,j,k  ,:,JDIR),  PhiJ  )
                    call PWrap(Grid%xfc(i,j,k  ,:,KDIR),  PhiK  )
                    call PWrap(Grid%xfc(i+1,j,k,:,IDIR),  PhiIp )
                    call PWrap(Grid%xfc(i,j+1,k,:,JDIR),  PhiJp )
                    call PWrap(Grid%xfc(i,j,k+1,:,KDIR),  PhiKp )

                    !Get face normals at each face center
                    nI  = Grid%Tf(i,j,k  ,NORMX:NORMZ,IDIR)
                    nJ  = Grid%Tf(i,j,k  ,NORMX:NORMZ,JDIR)
                    nK  = Grid%Tf(i,j,k  ,NORMX:NORMZ,KDIR)
                    nIp = Grid%Tf(i+1,j,k,NORMX:NORMZ,IDIR)
                    nJp = Grid%Tf(i,j+1,k,NORMX:NORMZ,JDIR)
                    nKp = Grid%Tf(i,j,k+1,NORMX:NORMZ,KDIR)

                    !Get face areas
                    daI  = Grid%face(i,j,k  ,IDIR)
                    daJ  = Grid%face(i,j,k  ,JDIR)
                    daK  = Grid%face(i,j,k  ,KDIR)
                    daIp = Grid%face(i+1,j,k,IDIR)
                    daJp = Grid%face(i,j+1,k,JDIR)
                    daKp = Grid%face(i,j,k+1,KDIR)

                    !Finally calculate g = -grad Phi
                    g = -(  PhiIp*daIp*nIp - PhiI*daI*nI &
                          + PhiJp*daJp*nJp - PhiJ*daJ*nJ &
                          + PhiKp*daKp*nKp - PhiK*daK*nK )/dV
                    rHat = normVec(Grid%xyzcc(i,j,k,:))
                    if (Model%doSphGrav) then
                        !Only radial component
                        Grid%gxyz(i,j,k,:) = Vec2Para(g,rHat)
                    else
                        Grid%gxyz(i,j,k,:) = g
                    endif
                enddo
            enddo
        enddo

        !Don't initialize again
        Grid%doG0Init = .false.

        contains
            !Wrapper (to take vector)
            subroutine PWrap(r,pot)
                real(rp), intent(in) :: r(NDIM)
                real(rp), intent(out) :: pot

                call Phi(r(XDIR),r(YDIR),r(ZDIR),pot)
            end subroutine PWrap
    end subroutine AddGrav

end module background

