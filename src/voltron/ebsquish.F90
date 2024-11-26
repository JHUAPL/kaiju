!Routines to handle CHIMP 2D<->3D mapping for deep coupling
module ebsquish
    use kdefs
    use gamtypes
    use ebtypes
    use volttypes
    use streamline
    use planethelper
    use imaghelper
    use clocks, only: Tic,Toc
    
    implicit none

    !Projection type
    abstract interface
        subroutine Projection_T(ebApp,xyz,t,x1,x2,isGood,MaxStep)
            Import :: rp,NDIM,ebTrcApp_T
            type(ebTrcApp_T), intent(in) :: ebApp
            real(rp), intent(in)  :: xyz(NDIM), t
            real(rp), intent(out) :: x1,x2
            logical , intent(out) :: isGood
            integer , intent(in), optional :: MaxStep
        end subroutine Projection_T
    end interface

    real(rp), parameter, private :: startEps  = 0.05
    real(rp), parameter, private :: rEps      = 0.125
    real(rp), parameter, private :: ShueScl   = 2.0 !Safety factor for Shue MP
    contains

    !Adjust the starting indices of the squish blocks according to the passed in array values
    subroutine LoadBalanceBlocks(vApp, balanceVals)
        class(voltApp_T), intent(inout) :: vApp
        real(rp), optional, intent(in) :: balanceVals(vApp%ebTrcApp%ebSquish%numSquishBlocks)

        real(rp), dimension(:), allocatable, save :: assignedPercents
        real(rp) :: totalPercent, balanceValsSum, defaultBalance
        integer :: b

        associate(ebSquish=>vApp%ebTrcApp%ebSquish,ebGr=>vApp%ebTrcApp%ebState%ebGr)

        if(.not. allocated(assignedPercents)) allocate(assignedPercents(ebSquish%numSquishBlocks))
        if(ebSquish%numSquishBlocks /= SIZE(assignedPercents)) then
            deallocate(assignedPercents)
            allocate(assignedPercents(ebSquish%numSquishBlocks))
        endif

        if(.not. present(balanceVals)) then
            ! make all blocks equal
            assignedPercents = 1.0_rp/ebSquish%numSquishBlocks
        else
            ! use the mean of balances as the default
            defaultBalance = sum(balanceVals)/size(balanceVals)
            ! calculate the relative size of each block
            balanceValsSum = 0.0_rp
            do b=1,ebSquish%numSquishBlocks
                if(balanceVals(b) <= 0) then
                    ! we can't have a block with size zero
                    balanceValsSum = balanceValsSum + defaultBalance
                else
                    ! this input is valid
                    balanceValsSum = balanceValsSum + balanceVals(b)
                endif
            enddo

            assignedPercents = 0.0_rp
            do b=1,ebSquish%numSquishBlocks
                if(balanceVals(b) <= 0) then
                    ! using default of 1
                    assignedPercents(b) = defaultBalance / balanceValsSum
                else
                    assignedPercents(b) = balanceVals(b) / balanceValsSum
                endif
            enddo
        endif

        ! use the calculated percentage of each block to assign starting indices
        totalPercent = 0.0_rp
        do b=1,ebSquish%numSquishBlocks
            ebSquish%blockStartIndices(b) = GetIndexAtPercent(vApp, totalPercent)
            totalPercent = totalPercent + assignedPercents(b)
        enddo

        end associate

    end subroutine

    !Find i-index of outer boundary of coupling domain
    function ShellBoundary(gModel,Gr,R) result(iMax)
        type(Model_T), intent(in) :: gModel
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(in) :: R
        integer :: iMax
        integer :: i

        !Just using dayside line
        !Avoiding findloc because it's not well supported by GNU
        do i=Gr%is,Gr%ie
            iMax = i
            if (Gr%xyz(i,Gr%js,Gr%ks,XDIR)>=R) exit
        enddo

    end function ShellBoundary

    !Squish functions below
    !Do 3D->2D mapping from Gamera volume to 2D inner mag model grid

    ! Helper subroutine to perform all squish blocks
    subroutine Squish(vApp)
        class(voltApp_T), intent(inout) :: vApp

        call Tic("Squish",.true.)
        do while(SquishBlocksRemain(vApp))
            call DoSquishBlock(vApp)
        enddo
        call Toc("Squish",.true.)

    end subroutine Squish

    function SquishBlocksRemain(vApp)
        class(voltApp_T), intent(in) :: vApp
        logical :: SquishBlocksRemain

        integer :: blockCount

        associate(ebSquish=>vApp%ebTrcApp%ebSquish)

        if(ebSquish%myNumBlocks < 0) then
            ! do all blocks
            SquishBlocksRemain = ebSquish%curSquishBlock <= ebSquish%numSquishBlocks
        else
            ! do some blocks
            blockCount = ebSquish%curSquishBlock - ebSquish%myFirstBlock ! blocks completed
            SquishBlocksRemain = blockCount < ebSquish%myNumBlocks
        endif

        end associate

    end function

    !Setup squishy data
    subroutine SquishStart(vApp)
        class(voltApp_T), intent(inout) :: vApp

        call UpdateTM03(vApp%time)
        call Tic("Squish",.true.)
        associate(ebGr=>vApp%ebTrcApp%ebState%ebGr, &                  
                  xyzSquish=>vApp%chmp2mhd%xyzSquish,isGood=>vApp%chmp2mhd%isGood, &
                  ebSquish=>vApp%ebTrcApp%ebSquish)

        ebSquish%Rinner = norm2(ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,XDIR:ZDIR))

        xyzSquish = 0.0
        isGood = .false.

        !Force iDeep to be even
        vApp%iDeep = 2*ceiling(vApp%iDeep/2.0)

        ebSquish%curSquishBlock = ebSquish%myFirstBlock

        end associate
        call Toc("Squish",.true.)

    end subroutine SquishStart

    ! perform bulk of squish opration
    subroutine DoSquishBlock(vApp)
        class(voltApp_T), intent(inout) :: vApp

        integer :: i,j,k,nSkp
        integer :: ksB,keB
        real(rp) :: t,x1,x2
        real(rp), dimension(NDIM) :: xyz
        logical :: isGoodIJK

        procedure(Projection_T), pointer :: ProjectXYZ

        associate(ebModel=>vApp%ebTrcApp%ebModel,ebGr=>vApp%ebTrcApp%ebState%ebGr,ebState=>vApp%ebTrcApp%ebState, &
                  xyzSquish=>vApp%chmp2mhd%xyzSquish,isGood=>vApp%chmp2mhd%isGood, &
                  ebSquish=>vApp%ebTrcApp%ebSquish)

        ProjectXYZ => NULL()

        !Use possibly looser tolerance for squish
        call SetTrcEpsilon(ebModel,vApp%chmp2mhd%epsSquish)

        select case (vApp%prType)
        case(LPPROJ)
            ProjectXYZ => Proj2LP
        case(LLPROJ)
            ProjectXYZ => Proj2LL
        case DEFAULT
            write(*,*) 'Unkown projection type, bailing ...'
            stop
        end select

        t = ebState%eb1%time

        if (vApp%doQkSquish) then
            nSkp = vApp%qkSquishStride !Stride through grid for projections
        else
            nSkp = 1
        endif

        call GetSquishBds(vApp,ksB,keB)
        call Tic("SQ-Project")

        if (ebModel%doDip) write(*,*) "Using fake projection for testing!"

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,xyz,x1,x2,isGoodIJK)
        do k=ksB,keB,nSkp
            do j=ebGr%js,ebGr%je+1,nSkp
                do i=ebGr%is,vApp%iDeep+1,nSkp
                    xyz = ebGr%xyz(i,j,k,XDIR:ZDIR)
                    if (norm2(xyz) <= vApp%rTrc) then
                        !Do projection
                        call ProjectXYZ(vApp%ebTrcApp,xyz,t,x1,x2,isGoodIJK,vApp%nTrc)
                    else
                        !Set null projection because outside radius
                        x1 = 0.0
                        x2 = 0.0
                        isGoodIJK = .false.
                    endif
                    
                    xyzSquish(i,j,k,1) = x1
                    xyzSquish(i,j,k,2) = x2

                    isGood(i,j,k) = isGoodIJK
                enddo
            enddo
        enddo

        !Reset tolerance
        call SetTrcEpsilon(ebModel,vApp%chmp2mhd%epsds0)

        ebSquish%curSquishBlock = ebSquish%curSquishBlock + 1

        end associate
        call Toc("SQ-Project")
    end subroutine DoSquishBlock

    !Get the k index pct % of the way through the squish data
    !where this will obey all rules about quick squish stride,
    !using the correct start index, etc..
    function GetIndexAtPercent(vApp, pct) result(pctIndex)
        class(voltApp_T), intent(in) :: vApp
        real(rp), intent(in) :: pct
        integer :: pctIndex

        integer :: Nk,nSkp,offset
        real(rp) :: correctPct

        associate(ebSquish=>vApp%ebTrcApp%ebSquish,ebGr=>vApp%ebTrcApp%ebState%ebGr)

        correctPct = MIN(1.0_rp, MAX(0.0_rp, pct)) ! ensure pct is between 0 and 1

        Nk = ebGr%ke+1 - ebGr%ks + 1
        if (vApp%doQkSquish) then
            nSkp = vApp%qkSquishStride !Stride through grid for projections
        else
            nSkp = 1
        endif

        offset = Nk*correctPct
        offset = offset - MOD(offset,nSkp) ! ensure index is aligned with quick squish
        pctIndex = ebGr%ks + offset
        end associate

    end function

    !Get squish bounds for block n (out of Nblk)
    !ksB and keB are the bounds of this block
    !n in [1,Nblk] (now more in line with Fortran indexing)
    subroutine GetSquishBds(vApp,ksB,keB,blockNum)
        class(voltApp_T), intent(in) :: vApp
        integer, intent(out) :: ksB,keB
        integer, optional, intent(in) :: blockNum

        integer :: curB

        associate(ebSquish=>vApp%ebTrcApp%ebSquish,ebGr=>vApp%ebTrcApp%ebState%ebGr)

        if(present(blockNum)) then
            curB = blockNum
        else
            curB = ebSquish%curSquishBlock
        endif

        if(curB<ebSquish%numSquishBlocks) then
            ksB = ebSquish%blockStartIndices(curB)
            keB = ebSquish%blockStartIndices(curB+1)
        else
            ksB = ebSquish%blockStartIndices(curB)
            keB = ebGr%ke+1 !Make sure last block finishes everything
        endif

        end associate

    end subroutine GetSquishBds

    ! Perform final operations on squishy data
    subroutine SquishEnd(vApp)
        class(voltApp_T), intent(inout) :: vApp

        integer :: i,Nk

        call Tic("Squish",.true.)
        associate(ebModel=>vApp%ebTrcApp%ebModel,ebGr=>vApp%ebTrcApp%ebState%ebGr, &
                  xyzSquish=>vApp%chmp2mhd%xyzSquish,isGood=>vApp%chmp2mhd%isGood, &
                  ebSquish=>vApp%ebTrcApp%ebSquish)

        if (vApp%doQkSquish) then
            call Tic("SQ-FillSkip")
            call FillSkips(ebModel,ebGr,vApp%iDeep,xyzSquish,isGood,vApp%qkSquishStride)
            call Toc("SQ-FillSkip")
        endif

        Nk = ebGr%ke-ebGr%ks+1

        !Ensure same value for degenerate axis points
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i)
        do i=ebGr%is,vApp%iDeep+1
            !x1 (regular average)
            xyzSquish(i,ebGr%js  ,:,1)=ArithMean(xyzSquish(i,ebGr%js  ,ebGr%ks:ebGr%ke,1))
            xyzSquish(i,ebGr%je+1,:,1)=ArithMean(xyzSquish(i,ebGr%je+1,ebGr%ks:ebGr%ke,1))

            !x2 (circular average)
            xyzSquish(i,ebGr%js  ,:,2)=CircMean(xyzSquish(i,ebGr%js  ,ebGr%ks:ebGr%ke,2))
            xyzSquish(i,ebGr%je+1,:,2)=CircMean(xyzSquish(i,ebGr%je+1,ebGr%ks:ebGr%ke,2))
        enddo

        vApp%chmp2mhd%iMax = vApp%iDeep
        
        end associate
        call Toc("Squish",.true.)

    end subroutine SquishEnd

    !Linearly interpolate between the stride 2 projections
    recursive subroutine FillSkips(ebModel,ebGr,iDeep,xyzSquish,isGood,nSkp)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebGrid_T)   , intent(in) :: ebGr
        integer          , intent(in) :: iDeep
        real(rp), intent(inout) :: xyzSquish(ebGr%is:ebGr%ie+1,ebGr%js:ebGr%je+1,ebGr%ks:ebGr%ke+1,2)
        logical , intent(inout) :: isGood   (ebGr%is:ebGr%ie+1,ebGr%js:ebGr%je+1,ebGr%ks:ebGr%ke+1)
        integer          , intent(in) :: nSkp

        integer :: i,j,k
        real(rp), dimension(2) :: X1s,X2s
        

        integer :: hSkp
        hSkp = nSkp/2


        !Start w/ i edge sweep
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,X1s,X2s)
        do k=ebGr%ks,ebGr%ke+1,nSkp
            do j=ebGr%js,ebGr%je+1,nSkp
                do i=ebGr%is+hSkp,iDeep,nSkp
                    if ( isGood(i-hSkp,j,k) .and. isGood(i+hSkp,j,k) ) then
                        X1s = [xyzSquish(i-hSkp,j,k,1),xyzSquish(i+hSkp,j,k,1)]
                        X2s = [xyzSquish(i-hSkp,j,k,2),xyzSquish(i+hSkp,j,k,2)]
                        xyzSquish(i,j,k,1) = ArithMean(X1s)
                        xyzSquish(i,j,k,2) = CircMean (X2s)

                        isGood(i,j,k) = .true.
                    endif
                enddo
            enddo
        enddo

        !Next do denser j sweep
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,X1s,X2s)
        do k=ebGr%ks,ebGr%ke+1,nSkp
            do j=ebGr%js+hSkp,ebGr%je,nSkp
                do i=ebGr%is,iDeep+1,hSkp
                    if ( isGood(i,j-hSkp,k) .and. isGood(i,j+hSkp,k) ) then
                        X1s = [xyzSquish(i,j-hSkp,k,1),xyzSquish(i,j+hSkp,k,1)]
                        X2s = [xyzSquish(i,j-hSkp,k,2),xyzSquish(i,j+hSkp,k,2)]
                        xyzSquish(i,j,k,1) = ArithMean(X1s)
                        xyzSquish(i,j,k,2) = CircMean (X2s)

                        isGood(i,j,k) = .true.
                    endif
                enddo
            enddo
        enddo

        !Finally, finish with denser-er k sweep
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,X1s,X2s)
        do k=ebGr%ks+hSkp,ebGr%ke,nSkp
            do j=ebGr%js,ebGr%je+1,hSkp
                do i=ebGr%is,iDeep+1,hSkp
                    if ( isGood(i,j,k-hSkp) .and. isGood(i,j,k+hSkp) ) then
                        X1s = [xyzSquish(i,j,k-hSkp,1),xyzSquish(i,j,k+hSkp,1)]
                        X2s = [xyzSquish(i,j,k-hSkp,2),xyzSquish(i,j,k+hSkp,2)]
                        xyzSquish(i,j,k,1) = ArithMean(X1s)
                        xyzSquish(i,j,k,2) = CircMean (X2s)

                        isGood(i,j,k) = .true.
                    endif

                enddo
            enddo
        enddo

        if(nSkp > 2) call FillSkips(ebModel,ebGr,iDeep,xyzSquish,isGood,nSkp/2) ! recurse with half skip size

    end subroutine FillSkips

    !Project XYZ to R,phi at Z=0 plane
    subroutine Proj2LP(ebApp,xyz,t,x1,x2,isGood,MaxStep)
        type(ebTrcApp_T)  , intent(in) :: ebApp
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(in) :: t
        real(rp), intent(out) :: x1,x2
        logical , intent(out) :: isGood
        integer , intent(in), optional :: MaxStep

        real(rp) :: L,phi,z
        real(rp), dimension(NDIM) :: xyzSeed,xy0

        ! trap for when we're within epsilon of the inner boundary
        ! (really, it's probably only the first shell of nodes at R=Rinner_boundary that doesn't trace correctly)
        if ( (norm2(xyz)-ebApp%ebSquish%Rinner)/ebApp%ebSquish%Rinner < startEps ) then
           ! dipole-shift to startEps
           xyzSeed = DipoleShift(xyz,norm2(xyz)+startEps)
        else
           xyzSeed = xyz
        end if

        ! note, this assumes internally tracing along -B for z>0 and along B for z<0
        ! if that fails, it returns xy0=HUGE
        call getEquatorProjection(ebApp%ebModel,ebApp%ebState,xyzSeed,t,xy0)

        !Map projection to L,phi
        L = sqrt(xy0(XDIR)**2.0+xy0(YDIR)**2.0)
        phi = atan2(xy0(YDIR),xy0(XDIR))
        if (phi<0) phi = phi+2*PI
        
        x1 = L
        x2 = phi
        if (x1 > 0) then
            isGood = .true.
        else
            isGood = .false.
        endif

    end subroutine Proj2LP

    !Project XYZ to lat-lon on ionosphere
    subroutine Proj2LL(ebApp,xyz,t,x1,x2,isGood,MaxStep)
        type(ebTrcApp_T), intent(in) :: ebApp
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(in) :: t
        real(rp), intent(out) :: x1,x2
        logical , intent(out) :: isGood
        integer , intent(in), optional :: MaxStep

        real(rp), dimension(NDIM) :: xE,xyz0
        integer :: Np
        logical :: isGP
        real(rp) :: dX,rC
        real(rp) :: x11,x22

        x1 = 0.0
        x2 = 0.0
        xE = 0.0
        
        if (ebApp%ebModel%doDip) then
            isGood = .true.
            x1 = InvLatitude(xyz)
            x2 = katan2(xyz(YDIR),xyz(XDIR))
            return
        endif

        !Do quick short-cut to save us some effort
        isGood = inShueMP_SM(xyz,ShueScl)
        if (.not. isGood) return

        ! trap for when we're within epsilon of the inner boundary
        ! (really, it's probably only the first shell of nodes at R=Rinner_boundary that doesn't trace correctly)
        if ( (norm2(xyz)-ebApp%ebSquish%Rinner)/ebApp%ebSquish%Rinner < startEps ) then
           ! dipole-shift to startEps
           xyz0 = DipoleShift(xyz,norm2(xyz)+startEps)
        else
           xyz0 = xyz
        end if
        if (present(MaxStep)) then
            call mageproject(ebApp%ebModel,ebApp%ebState,xyz0,t,xE,Np,isGP,MaxStepsO=MaxStep)
        else
            call mageproject(ebApp%ebModel,ebApp%ebState,xyz0,t,xE,Np,isGP)
        endif
        dX = norm2(xyz0-xE)
        rC = ebApp%ebSquish%Rinner*(1.+rEps)

        isGood = isGP .and. (dX>TINY) .and. (xE(ZDIR)>0.0) .and. (norm2(xE)<rC)

        if (isGood) then
            !Get invariant lat/lon
            x1 = InvLatitude(xE)
            x2 = atan2(xE(YDIR),xE(XDIR))
            if (x2 < 0) x2 = x2 + 2*PI
        else
            x1 = 0.0
            x2 = 0.0
        endif

    end subroutine Proj2LL

end module ebsquish
