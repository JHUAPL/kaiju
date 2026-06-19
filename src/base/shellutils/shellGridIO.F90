module shellGridIO

    use kdefs
    use ioH5
    use shellGrid

    implicit none

    integer, parameter, private :: MAXIOVAR = 50

    ! Overloader to write a shellGrid var to file
    interface AddOutSGV
        module procedure AddOutSGV_1D, AddOutSGV_0D
    end interface
    ! Overloader to read a shellGrid var to file
    interface ReadInSGV
        module procedure ReadInSGV_1D, ReadInSGV_0D
    end interface
    
    contains

    subroutine writeShellGrid(sg, outH5, gStrO)
        !! Writes ShellGrid data to file outH5
        !! Note: This does not save every variable,
        !!  just ones needed for restarts and generally useful output
        type(ShellGrid_T), intent(in) :: sg
            !! Shell Grid object to write
        character(len=*), intent(in) :: outH5
            !! Output file name
        character(len=*), optional, intent(in) :: gStrO
            !! Optional group to write to, default = /ShellGrid

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

        call ClearIO(IOVars)

        ! Attrs
        call AddOutVar(IOVars, "name", trim(sg%name))
        call AddOutVar(IOVars, "radius", sg%radius,uStr="Rp")
        call AddOutVar(IOVars, "Nt", sg%Nt)
        call AddOutVar(IOVars, "Np", sg%Np)
        call AddOutVar(IOVars, "nGhosts_n", sg%Ngn)
        call AddOutVar(IOVars, "nGhosts_s", sg%Ngs)
        call AddOutVar(IOVars, "nGhosts_e", sg%Nge)
        call AddOutVar(IOVars, "nGhosts_w", sg%Ngw)
        if (sg%isChild) then
            call AddOutVar(IOVars, "isChild", 1)
        else
            call AddOutVar(IOVars, "isChild", 0)
        endif
        call AddOutVar(IOVars, "parentName", trim(sg%parentName))
        call AddOutVar(IOVars, "bndis", sg%bndis)
        call AddOutVar(IOVars, "bndie", sg%bndie)
        call AddOutVar(IOVars, "bndjs", sg%bndjs)
        call AddOutVar(IOVars, "bndje", sg%bndje)

        ! Arrays
        call AddOutVar(IOVars, "theta", sg%th, uStr="radians")
        call AddOutVar(IOVars, "phi"  , sg%ph, uStr="radians")


        if (present(gStrO)) then
            call WriteVars(IOVars, .false., outH5, trim(gStrO))
        else
            call WriteVars(IOVars, .false., outH5, '/ShellGrid')
        endif

    end subroutine writeShellGrid


    subroutine GenShellGridFromFile(sg, sgName, inH5, gStrO)
        !! Generates ShellGrid object using previously output SG data
        type(ShellGrid_T), intent(inout) :: sg
            !! Shell Grid object to write
        character(len=*), intent(in) :: sgName
            !! Name to assign to generated ShellGrid since we can't get strings from h5 files
        character(len=*), intent(in) :: inH5
            !! Output file name
        character(len=*), optional, intent(in) :: gStrO
            !! Optional group to read from, default = /ShellGrid

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        character(len=strLen) :: gStr
        integer :: Nt, Np
        real(rp), dimension(:), allocatable :: theta, phi, theta_active, phi_active
        character(len=strLen) :: name
        integer, dimension(4) :: nGhosts 
        real(rp) :: radius
        logical :: isChild
        character(len=strLen) :: parentName
        integer :: sub_is, sub_ie, sub_js, sub_je

        call ClearIO(IOVars)

        if(present(gStrO)) then
            gStr = trim(gStrO)
        else
            gStr = "/ShellGrid"
        endif
        
        ! No name reads yet cause we can't read strings. 
        ! Also makes it tricky to properly do child grids, so that's not an option right now
        

        call AddInVar(IOVars, "name",vTypeO=IOSTR)
        call AddInVar(IOVars, "radius")
        call AddInVar(IOVars, "Nt")
        call AddInVar(IOVars, "Np")
        call AddInVar(IOVars, "nGhosts_n")
        call AddInVar(IOVars, "nGhosts_s")
        call AddInVar(IOVars, "nGhosts_e")
        call AddInVar(IOVars, "nGhosts_w")
        call AddInVar(IOVars, "theta")
        call AddInVar(IOVars, "phi")
        call AddInVar(IOVars, "isChild")
        call AddInVar(IOVars, "parentName",vTypeO=IOSTR)
        call AddInVar(IOVars, "bndis")
        call AddInVar(IOVars, "bndie")
        call AddInVar(IOVars, "bndjs")
        call AddInVar(IOVars, "bndje")

        call ReadVars(IOVars, .false., inH5, trim(gStr))
        
        !isChild = GetIOInt(IOVars, 'isChild') .eq. 1
        !if (isChild) then
        !    write(*,*) "ERROR: Reading child ShellGrid from file currently not supported."
        !    write(*,*) "  Gotta make some decisions regarding how this should be handled."
        !    write(*,*) "  Goodbye."
        !    stop
        !endif


        Nt = GetIOInt(IOVars, 'Nt')
        Np = GetIOInt(IOVars, 'Np')
        nGhosts(NORTH) = GetIOInt(IOVars, 'nGhosts_n')
        nGhosts(SOUTH) = GetIOInt(IOVars, 'nGhosts_s')
        nGhosts(EAST)  = GetIOInt(IOVars, 'nGhosts_e')
        nGhosts(WEST)  = GetIOInt(IOVars, 'nGhosts_w')
        radius = GetIOReal(IOVars, 'radius')
        
        allocate(theta(Nt + nGhosts(NORTH) + nGhosts(SOUTH) + 1))
        allocate(phi  (Np + nGhosts(WEST)  + nGhosts(EAST)  + 1))
        allocate(theta_active(Nt + 1))
        allocate(phi_active  (Np + 1))
        call IOArray1DFill(IOVars, 'theta', theta)
        call IOArray1DFill(IOVars, 'phi'  , phi)
        theta_active(:) = theta(1 + nGhosts(NORTH) : Nt + nGhosts(NORTH) + 1)
        phi_active(:)   = phi  (1 + nGhosts(WEST)  : Np + nGhosts(WEST)  + 1)

        call GenShellGrid(sg, theta_active, phi_active, trim(sgName), nGhosts=nGhosts, radO=radius)

        ! If we are actually a child grid, sneak in afterward and set that stuff up
        isChild = GetIOInt(IOVars, 'isChild') .eq. 1
        if (isChild) then
            sg%isChild = .true.
            sg%parentName = trim(GetIOStr(IOVars, 'parentName'))
            sg%bndis = GetIOInt(IOVars, 'bndis')
            sg%bndie = GetIOInt(IOVars, 'bndie')
            sg%bndjs = GetIOInt(IOVars, 'bndjs')
            sg%bndje = GetIOInt(IOVars, 'bndje')
        endif
        
    end subroutine GenShellGridFromFile


!------
! ShellGridVar write overloads
!------

    subroutine AddOutSGV_0D(IOVars, idStr, sgv, uStr, dStr, outBndsO, doWriteMaskO)
        type(IOVar_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        type(ShellGridVar_T), intent(in) :: sgv
        character(len=*), intent(in), optional :: uStr,dStr
        integer, dimension(4), intent(in), optional :: outBndsO
        logical, intent(in), optional :: doWriteMaskO

        integer :: is, ie, js, je
        character(len=strLen) :: idStr_mask
        character(len=strLen) :: mask_desc
        logical :: doWriteMask
        real(rp), dimension(:,:), allocatable :: Q_mask

        if(present(outBndsO)) then
            is = outBndsO(1); ie = outBndsO(2); js = outBndsO(3); je = outBndsO(4)
        else
            is = sgv%isv    ; ie = sgv%iev    ; js = sgv%jsv    ; je = sgv%jev
        endif

        if(present(doWriteMaskO)) then
            doWriteMask = doWriteMaskO
        else
            doWriteMask = .false.
        endif

        call AddOutVar(IOVars, idStr, sgv%data(is:ie,js:je), uStr=uStr, dStr=dStr)
        if(doWriteMask) then
            write(idStr_mask, "(A,A)") trim(idStr), "_mask"
            write(mask_desc , "(A,A)") "Mask array for variable ",trim(idStr)
            allocate( Q_mask(sgv%isv:sgv%iev, sgv%jsv:sgv%jev) )
            Q_mask = merge(1_rp, 0_rp, sgv%mask)
            call AddOutVar(IOVars, idStr_mask, Q_mask(is:ie,js:je), dStr=mask_desc)
        endif
        
    end subroutine AddOutSGV_0D


    subroutine AddOutSGV_1D(IOVars, idStr, sgv, uStr, dStr, outBndsO, doWriteMaskO)
        type(IOVar_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        type(ShellGridVar_T), dimension(:), intent(in) :: sgv
        character(len=*), intent(in), optional :: uStr,dStr
        integer, dimension(4), intent(in), optional :: outBndsO
        logical, intent(in), optional :: doWriteMaskO

        integer :: is, ie, js, je, k
        integer, dimension(1) :: sgv_shape
        character(len=strLen) :: idStr_mask
        character(len=strLen) :: mask_desc
        logical :: doWriteMask
        real(rp), dimension(:,:,:), allocatable :: Q
        real(rp), dimension(:,:,:), allocatable :: Q_mask

        if(present(doWriteMaskO)) then
            doWriteMask = doWriteMaskO
        else
            doWriteMask = .false.
        endif

        if(present(outBndsO)) then
            is = outBndsO(1); ie = outBndsO(2); js = outBndsO(3); je = outBndsO(4)
        else
            is = sgv(1)%isv ; ie = sgv(1)%iev ; js = sgv(1)%jsv ; je = sgv(1)%jev
        endif

        sgv_shape = shape(sgv)
        allocate(Q(sgv(1)%isv:sgv(1)%iev, sgv(1)%jsv:sgv(1)%jev, sgv_shape(1)))
        do k=1,sgv_shape(1)
            ! We are assuming all sgv's in array have identical bounds. Make sure that's true
            if ((sgv(k)%isv .ne. sgv(1)%isv) .or. (sgv(k)%iev .ne. sgv(1)%iev) .or. (sgv(k)%jsv .ne. sgv(1)%jsv) .or. (sgv(k)%jev .ne. sgv(1)%jev)) then
                write(*,*) "ERROR writing SGV_1D with id=",trim(idStr)
                write(*,*) "Dims not equal, dying..."
                stop
            endif
            Q(:,:,k) = sgv(k)%data
        enddo
        call AddOutVar(IOVars, idStr, Q(is:ie,js:je,:), uStr=uStr, dStr=dStr)

        if (doWriteMask) then
            write(idStr_mask, "(A,A)") trim(idStr), "_mask"
            write(mask_desc , "(A,A)") "Mask array for variable ",trim(idStr)
            allocate( Q_mask(sgv(1)%isv:sgv(1)%iev, sgv(1)%jsv:sgv(1)%jev, sgv_shape(1)) )
            do k=1,sgv_shape(1)
                Q_mask(:,:,k) = merge(1_rp, 0_rp, sgv(k)%mask)
            enddo
            call AddOutVar(IOVars, idStr_mask, Q_mask(is:ie,js:je,:), dStr=mask_desc)
        endif
        
    end subroutine AddOutSGV_1D

    subroutine ReadInSGV_0D(sgv, baseStr, idStr, gStrO, doIOpO)
        type(ShellGridVar_T), intent(inout) :: sgv
        character(len=*), intent(in) :: baseStr
        character(len=*), intent(in) :: idStr
        character(len=*), intent(in), optional :: gStrO
        logical, intent(in), optional :: doIOpO
        
        type(IOVAR_T), dimension(5) :: IOVars
        logical :: doGroup, maskExists
        character(len=strLen) :: gStr
        logical :: doIOp
        character(len=strLen) :: idStr_mask
        logical :: doReadMask
        real(rp), dimension(:,:), allocatable :: Q

        if(present(gStrO)) then
            doGroup = .true.
        else
            doGroup = .false.
        endif

        doIOp = .false.
        if (present(doIOpO)) doIOp = doIOpO
        allocate(Q(sgv%isv:sgv%iev, sgv%jsv:sgv%jev))
        write(idStr_mask, "(A,A)") trim(idStr), "_mask"

        call ClearIO(IOVars)
        call AddInVar(IOVars, idStr)

        if (doGroup) then
            maskExists = ioExist(baseStr, idStr_mask, gStrO)
        else
            maskExists = ioExist(baseStr, idStr_mask)
        endif

        if (maskExists) then
            doReadMask = .true.
            call AddInVar(IOVars, idStr_mask)
        else
            write(*,*)"ReadInSGV_0D: Did not find mask variable for id=",trim(idStr)
            doReadMask = .false.
        endif
        if(doGroup) then
            call ReadVars(IOVars, doIOp, baseStr, gStrO)
        else
            call ReadVars(IOVars, doIOp, baseStr)
        endif

        call IOArray2DFill(IOVars, idStr, sgv%data)
        if (doReadMask) then
            call IOArray2DFill(IOVars, idStr_mask, Q)
            sgv%mask = merge(.true., .false., Q > 0.5)
        endif
    end subroutine ReadInSGV_0D


    subroutine ReadInSGV_1D(sgv, baseStr, idStr, gStrO, doIOpO)
        type(ShellGridVar_T), dimension(:), intent(inout) :: sgv
        character(len=*), intent(in) :: baseStr
        character(len=*), intent(in) :: idStr
        character(len=*), intent(in), optional :: gStrO
        logical, intent(in), optional :: doIOpO
        
        integer :: k
        type(IOVAR_T), dimension(5) :: IOVars
        logical :: doGroup, maskExists
        character(len=strLen) :: gStr
        logical :: doIOp
        character(len=strLen) :: idStr_mask
        integer, dimension(1) :: sgv_shape
        logical :: doReadMask
        real(rp), dimension(:,:,:), allocatable :: Q

        if(present(gStrO)) then
            doGroup = .true.
        else
            doGroup = .false.
        endif

        doIOp = .false.
        if (present(doIOpO)) doIOp = doIOpO
        sgv_shape = shape(sgv)

        allocate(Q(sgv(1)%isv:sgv(1)%iev, sgv(1)%jsv:sgv(1)%jev, sgv_shape(1)))
        write(idStr_mask, "(A,A)") trim(idStr), "_mask"

        call ClearIO(IOVars)
        call AddInVar(IOVars, idStr)

        if (doGroup) then
            maskExists = ioExist(baseStr, idStr_mask, gStrO)
        else
            maskExists = ioExist(baseStr, idStr_mask)
        endif

        if (maskExists) then
            doReadMask = .true.
            call AddInVar(IOVars, idStr_mask)
        else
            write(*,*)"ReadInSGV_1D: Did not find mask variable for id=",trim(idStr)
            doReadMask = .false.
        endif
        if(doGroup) then
            call ReadVars(IOVars, doIOp, baseStr, gStrO)
        else
            call ReadVars(IOVars, doIOp, baseStr)
        endif

        call IOArray3DFill(IOVars, idStr, Q)
        do k=1,sgv_shape(1)
            ! We are assuming all sgv's in array have identical bounds. Make sure that's true
            if ((sgv(k)%isv .ne. sgv(1)%isv) .or. (sgv(k)%iev .ne. sgv(1)%iev) .or. (sgv(k)%jsv .ne. sgv(1)%jsv) .or. (sgv(k)%jev .ne. sgv(1)%jev)) then
                write(*,*) "ERROR reading SGV_1D with id=",trim(idStr)
                write(*,*) "Dims not equal, dying..."
                stop
            endif
            sgv(k)%data = Q(:,:,k)
        enddo
        
        if (doReadMask) then
            !allocate(Q_mask(sgv%isv:sgv%iev, sgv%jsv:sgv%jev))
            call IOArray3DFill(IOVars, idStr_mask, Q)
            do k=1,sgv_shape(1)
                sgv(k)%mask = merge(.true., .false., Q(:,:,k) > 0.5)
            enddo
        endif
    end subroutine ReadInSGV_1D


end module shellGridIO
