!Routines related to multifluid simulations

module multifluid
    use gamtypes
    use xml_input
    use gamutils
    
    implicit none

    !Information for a given fluid species
    type FLUID_T
        character(len=strLen) :: idStr="NONE"
        real(rp) :: dVac = TINY !Threshold for "vacuum" for this species
    end type FLUID_T

    !Array of fluid descriptors (1-nSpc)
    type(FLUID_T), dimension(:), allocatable :: Spcs

    contains

    subroutine InitMultiF(Model,xmlInp)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: xmlInp

        integer :: n
        real(rp) :: gDFloor !Global density floor
        character(len=strLen) :: sID,sTag

        !Need floors here, so may end up doing it twice
        call SetFloors(Model,xmlInp)
        gDFloor = dFloor

        !Find number of species
        call xmlInp%Set_Val(Model%nSpc,'multifluid/nSpc',1)

        !Allocate number of descriptors
        allocate(Spcs(Model%nSpc))

        !Loop through species and get information
        do n=1,Model%nSpc
            write(sID ,'(A,I0)') "fluid" , n
            write(sTag,'(A,I0,A)') "fluid" , n,"/id"

            call xmlInp%Set_Val(Spcs(n)%idStr,sTag,sID)
            write(sTag,'(A,I0,A)') "fluid" , n,"/dVac"            
            call xmlInp%Set_Val(Spcs(n)%dVac,sTag,gDFloor)
        enddo

    end subroutine InitMultiF

    !Convert multi-species conserved variables into bulk flow quantities
    !(and enforce Vperp equality w/ provided B)

    subroutine MultiF2Bulk(Model,U,B)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U(NVAR,BLK:Model%nSpc)
        real(rp), intent(in) :: B(NDIM)

        logical , dimension(Model%nSpc) :: isGood
        real(rp), dimension(NDIM) :: Vblk,bhat,Vperp,Vs
        real(rp), dimension(NVAR) :: pW,pCon
        integer :: n
        
        !Which fluids in this cell are good
        isGood = ( U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        
        if (.not. any(isGood)) then
        !Well this is awkward, let's see if we can figure something out
            !Fix up the first fluid
            pW(DEN) = max(Spcs(1)%dVac,dFloor)+TINY
            pW(VELX:VELZ) = 0.0
            pW(PRESSURE) = pFloor
            call CellP2C(Model,pW,U(:,1))
            !Nuke the rest
            if (Model%nSpc>1) then
                do n=2,Model%nSpc
                    pW(DEN) = dFloor
                    pW(VELX:VELZ) = 0.0
                    pW(PRESSURE) = pFloor
                    call CellP2C(Model,pW,U(:,n))
                enddo
            endif
            isGood = ( U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        endif

        !Sum conserved quantities to get bulk
        do n=1,NVAR
            U(n,BLK) = sum(U(n,1:Model%nSpc),mask=isGood)
        enddo

        !Now get perp velocity from bulk
        pCon = U(:,BLK)
        call CellC2P(Model,pCon,pW)
        Vblk  = pW(VELX:VELZ)
        bhat  = normVec(B)
        Vperp = Vec2Perp(Vblk,bhat) !Bulk perp velocity

        do n=1,Model%nSpc
            pCon = U(:,n)
            call CellC2P(Model,pCon,pW)

            if (isGood(n)) then
                !Enforce vperp agreement
                Vs = pW(VELX:VELZ)
                pW(VELX:VELZ) = Vperp + dot_product(Vs,bhat)*bhat
            else
                !no good, either set velocity to bulk/perp/0
                !Should vacuum move at the bulk speed?
                pW(VELX:VELZ) = 0.0
                !pW(VELX:VELZ) = Vblk
                pW(PRESSURE)  = pFloor
            endif

            call CellP2C(Model,pW,pCon)
            U(:,n) = pCon
        enddo

        !Finish up and get outta here
        do n=1,NVAR
            U(n,BLK) = sum(U(n,1:Model%nSpc),mask=isGood)
        enddo

    end subroutine MultiF2Bulk
    
    !Convert full Gas state to bulk, ie include loop
    subroutine State2Bulk(Model,Grid,State,doGhostsO)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        logical, optional, intent(in) :: doGhostsO

        integer :: i,j,k
        integer :: iMin,iMax,jMin,jMax,kMin,kMax
        real(rp), dimension(NDIM) :: B

        logical :: doGhosts = .true.

        if (present(doGhostsO)) then
            doGhosts = doGhostsO
        endif

        if (doGhosts) then
            iMin = Grid%isg
            iMax = Grid%ieg
            jMin = Grid%jsg
            jMax = Grid%jeg
            kMin = Grid%ksg
            kMax = Grid%keg
        else
            iMin = Grid%is
            iMax = Grid%ie
            jMin = Grid%js
            jMax = Grid%je
            kMin = Grid%ks
            kMax = Grid%ke
        endif

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,B)
        do k=kMin,kMax
            do j=jMin,jMax
                do i=iMin,iMax
                    B = State%Bxyz(i,j,k,:)
                    if (Model%doBackground) then
                        B = B + Grid%B0(i,j,k,:)
                    endif
                    call MultiF2Bulk(Model,State%Gas(i,j,k,:,:),B)
                enddo
            enddo
        enddo
    end subroutine State2Bulk
    
    function MultiFCs(Model,U) result(MaxCs)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: U(NVAR,BLK:Model%nSpc)
        real(rp) :: MaxCs

        integer :: n
        real(rp) :: Csn
        logical, dimension(Model%nSpc) :: isGood
        
        !Which fluids in this cell are good
        isGood = ( U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        MaxCs = 0.0

        do n=1,Model%nSpc
            if (isGood(n)) then
                call CellPress2Cs(Model,U(:,n),Csn)
                MaxCs = max(MaxCs,Csn)
            endif
        enddo

    end function MultiFCs

    function MultiFSpeed(Model,U) result(MaxV)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: U(NVAR,BLK:Model%nSpc)
        real(rp) :: MaxV

        integer :: n
        logical, dimension(Model%nSpc) :: isGood
        real(rp), dimension(NVAR) :: pW,pCon
        real(rp) :: Vn
        
        !Which fluids in this cell are good
        isGood = ( U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
        MaxV = 0.0

        do n=1,Model%nSpc
            if (isGood(n)) then
                pCon = U(:,n)
                call CellC2P(Model,pCon,pW)
                Vn = norm2(pW(VELX:VELZ))
                MaxV = max(MaxV,Vn)
            endif
        enddo

    end function MultiFSpeed

    !Answers the age old question, is fluid s0 good?
    function isGoodFluid(Model,U,s0)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: U(NVAR,BLK:Model%nSpc)
        integer , intent(in) :: s0
        logical :: isGoodFluid

        if ( (s0 > Model%nSpc) .or. (s0 < 0) ) then
            isGoodFluid = .false.
            return
        endif
        if (s0 == BLK) then
            !Bulk is always good, otherwise you fucked up
            isGoodFluid = .true.
            return
        endif
        isGoodFluid = U(DEN,s0) >= Spcs(s0)%dVac
    end function isGoodFluid
    
end module multifluid
