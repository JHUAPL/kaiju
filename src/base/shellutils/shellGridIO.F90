module shellGridIO

    use kdefs
    use ioH5
    use shellGrid

    implicit none

    integer, parameter, private :: MAXIOVAR = 50
    
    contains

    subroutine writeShellGrid(sg, outH5, gStrO)
        !! Writes ShellGrid data to file outH5
        !! Note: This does not save every variable,
        !!  just ones needed for restarts and generally useful output
        type(ShellGrid_T), intent(in) :: sg
            !! Shell Grid object to write
        character(len=strLen), intent(in) :: outH5
            !! Output file name
        character(len=strLen), optional, intent(in) :: gStrO
            !! Optional group to write to, default = /ShellGrid

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

        call ClearIO(IOVars)

        ! Attrs
        call AddOutVar(IOVars, "name", sg%name)
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
            call WriteVars(IOVars, .false., outH5, gStrO)
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
        character(len=strLen), intent(in) :: inH5
            !! Output file name
        character(len=strLen), optional, intent(in) :: gStrO
            !! Optional group to read from, default = /ShellGrid

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        integer :: Nt, Np
        real(rp), dimension(:), allocatable :: theta, phi, theta_active, phi_active
        character(len=strLen) :: name
        integer, dimension(4) :: nGhosts 
        real(rp) :: radius
        logical :: isChild
        character(len=strLen) :: parentName
        integer :: sub_is, sub_ie, sub_js, sub_je

        call ClearIO(IOVars)

        ! No name reads yet cause we can't read strings. 
        ! Also makes it tricky to properly do child grids, so that's not an option right now

        !call AddInVar(IOVars, "name")
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
        !call AddInVar(IOVars, "parentName")
        call AddInVar(IOVars, "bndis")
        call AddInVar(IOVars, "bndie")
        call AddInVar(IOVars, "bndjs")
        call AddInVar(IOVars, "bndje")

        if (present(gStrO)) then
            call ReadVars(IOVars, .false., inH5, gStrO)
        else
            call ReadVars(IOVars, .false., inH5, '/ShellGrid')
        endif

        
        isChild = GetIOInt(IOVars, 'isChild') .eq. 1
        if (isChild) then
            write(*,*) "ERROR: Reading child ShellGrid from file currently not supported."
            write(*,*) "  Gotta make some decisions regarding how this should be handled."
            write(*,*) "  Goodbye."
            stop
        endif


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

        call GenShellGrid(sg, theta_active, phi_active, sgName, nGhosts=nGhosts, radO=radius)

    end subroutine GenShellGridFromFile

end module shellGridIO