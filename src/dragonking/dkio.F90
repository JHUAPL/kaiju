! Native DragonKing restart and file output on the DK shell grid.
! Mirrors the structure of remix's mixio 
! Keeps DK self-sufficient across restarts for EMA persistence
module dkio
    use kdefs
    use dktypes
    use ioH5
    use files, only: CheckAndKill, MapSymLink

    implicit none
    private

    integer, parameter :: MAXDKIOVAR = 32
    type(IOVAR_T), dimension(MAXDKIOVAR), private :: IOVars

    ! Cached output file + run tag, set once in initDKIO
    character(len=strLen), private :: h5File = "", h5RunID = ""

    ! Dataset names shared by restart and file output so they round-trip.
    character(len=*), parameter :: DK_SIGP = "Pedersen conductance"
    character(len=*), parameter :: DK_SIGH = "Hall conductance"
    character(len=*), parameter :: DK_FAC  = "Field-aligned current"
    character(len=*), parameter :: DK_EAVG = "Average energy"
    character(len=*), parameter :: DK_NFLX = "Number flux"
    character(len=*), parameter :: DK_POT  = "Potential"
    character(len=*), parameter :: DK_ATYP = "Auroral type"
    character(len=*), parameter :: DK_DELE = "Mono potential drop"
    character(len=*), parameter :: DK_ZEAV = "MHD average energy"
    character(len=*), parameter :: DK_ZNFL = "MHD number flux"
    character(len=*), parameter :: DK_BETA = "RCM beta"

    ! Number of dkVariables_T arrays in the restart.
    integer, parameter, private :: NDKRESVAR = 6

    public :: initDKIO, writeDKRestart, readDKRestart, writeDKFileOutput

    contains

    ! Only SigmaP/SigmaH are genuinely stateful
    ! Others are recomputed but restart state is read instead of zeros.
    subroutine dkResVarRef(dkApp, n, name, uStr, dat)
        class(dkApp_T), target, intent(inout) :: dkApp
        integer,               intent(in)  :: n
        character(len=strLen), intent(out) :: name, uStr
        real(rp), pointer,     intent(out) :: dat(:,:)

        select case (n)
        ! --- conductance ---
        case ( 1); name = DK_SIGP; uStr = "S";           dat => dkApp%vars%SigmaP%data
        case ( 2); name = DK_SIGH; uStr = "S";           dat => dkApp%vars%SigmaH%data
        ! --- precipitation outputs ---
        case ( 3); name = DK_EAVG; uStr = "keV";         dat => dkApp%vars%Avg_Eng%data
        case ( 4); name = DK_NFLX; uStr = "1/cm2/s";     dat => dkApp%vars%Num_Flux%data
        case ( 5); name = DK_ATYP; uStr = "";            dat => dkApp%vars%Aur_Type%data
        case ( 6); name = DK_DELE; uStr = "kV";          dat => dkApp%vars%DeltaE%data
        case default
            write(*,*) "dkResVarRef: index out of range: ", n
            stop
        end select
    end subroutine dkResVarRef

    ! Initialize DK IO: cache the output filename/run tag
    subroutine initDKIO(dkApp, RunID, isRestart, nRes)
        class(dkApp_T),   intent(inout) :: dkApp
        character(len=*), intent(in)    :: RunID
        logical,          intent(in)    :: isRestart
        integer,          intent(in)    :: nRes

        logical :: fExist
        integer :: Np, Nt

        h5File  = trim(RunID) // ".dk.h5"
        h5RunID = trim(RunID)
        dkApp%IO%RunID = trim(RunID)

        inquire(file=h5File, exist=fExist)

        ! Restore DK state
        if (isRestart) call readDKRestart(dkApp, RunID, nRes)

        ! On restart with the output file already present, keep accumulating
        ! Step#N groups into it (grid already at root) and bail.
        if (isRestart .and. fExist) return

        ! Fresh run (or restart whose output file is missing): start clean and
        ! write the static grid once at root.
        call CheckAndKill(h5File)
        Np = dkApp%grid%Np
        Nt = dkApp%grid%Nt
        call ClearIO(IOVars)
        call AddOutVar(IOVars, "X", dkApp%grid%x(1:Nt,1:Np), uStr="unit sphere")
        call AddOutVar(IOVars, "Y", dkApp%grid%y(1:Nt,1:Np), uStr="unit sphere")
        call AddOutVar(IOVars, "UnitsID", "DragonKing")
        call WriteVars(IOVars, .false., h5File)
    end subroutine initDKIO

    ! Persist the DK conductance EMA state to <RunID>.dk.Res.<nRes>.h5.
    subroutine writeDKRestart(dkApp, nRes)
        class(dkApp_T), target, intent(inout) :: dkApp
        integer,        intent(in)    :: nRes

        character(len=strLen) :: h5Str, lnStr, nStr
        character(len=strLen) :: vName, vUnit
        real(rp), pointer     :: dat(:,:)
        integer               :: n

        if (nRes == -1) then
            write(*,*) 'nRes is -1 in writeDKRestart, now die'
            stop
        endif
        write(nStr,'(I0.5)') nRes
        write(h5Str,'(A,A,A,A)') trim(h5RunID), '.dk.Res.', trim(nStr), '.h5'

        call CheckAndKill(h5Str)

        call ClearIO(IOVars)
        do n=1,NDKRESVAR
            call dkResVarRef(dkApp, n, vName, vUnit, dat)
            call AddOutVar(IOVars, trim(vName), dat, uStr=trim(vUnit))
        enddo
        call AddOutVar(IOVars, "nRes", nRes)
        call AddOutVar(IOVars, "time", dkApp%IO%time)
        call AddOutVar(IOVars, "MJD",  dkApp%IO%mjd)
        call AddOutVar(IOVars, "tilt", dkApp%vars%tilt)
        call WriteVars(IOVars, .false., h5Str)

        ! Rolling "latest" symlink so a restart with nRes=-1 picks up the most
        ! recent DK restart (mirrors mix/gamera/voltron; respects doSymLinks).
        write(lnStr,'(A,A,A,A)') trim(h5RunID), '.dk.Res.', 'XXXXX', '.h5'
        call MapSymLink(h5Str, lnStr)
    end subroutine writeDKRestart


    ! Restore the dkResVarRef variables from <resId>.dk.Res.<nRes>.h5. If the file
    ! is absent leave vars as initialized and warn; then cold-start as before.
    subroutine readDKRestart(dkApp, resId, nRes)
        class(dkApp_T),   target, intent(inout) :: dkApp
        character(len=*), intent(in)    :: resId
        integer,          intent(in)    :: nRes

        character(len=strLen) :: h5Str, nStr
        character(len=strLen) :: vName, vUnit
        real(rp), pointer     :: dat(:,:)
        logical :: fExist
        integer :: n0, n, nMiss

        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write(nStr,'(I0.5)') nRes
        endif
        write(h5Str,'(A,A,A,A)') trim(resId), '.dk.Res.', trim(nStr), '.h5'

        inquire(file=trim(h5Str), exist=fExist)
        if (.not. fExist) then
            write(*,*) "DK restart file not found; conductance EMA will cold-start: ", trim(h5Str)
            return
        endif
        write(*,*) "DK restarting from: ", trim(h5Str)

        call ClearIO(IOVars)
        do n=1,NDKRESVAR
            call dkResVarRef(dkApp, n, vName, vUnit, dat)
            call AddInVar(IOVars, trim(vName))
        enddo
        call ReadVars(IOVars, .false., h5Str)

        ! Restart written before a variable was added keeps its initialized value
        nMiss = 0
        do n=1,NDKRESVAR
            call dkResVarRef(dkApp, n, vName, vUnit, dat)
            n0 = FindIO(IOVars, trim(vName))
            if (n0 > 0) then
                if (IOVars(n0)%isDone) then
                    call IOArray2DFill(IOVars, trim(vName), dat)
                else
                    nMiss = nMiss + 1
                endif
            endif
        enddo
        if (nMiss > 0) write(*,*) "DK restart: ", nMiss, &
            " of ", NDKRESVAR, " variables absent (older restart file); they will be recomputed."

        write(*,*) "Done readDKRestart"
    end subroutine readDKRestart


    ! Append the DK fields for this step as a Step#<nStep> group in the single
    ! <RunID>.dk.h5 file. Grid X/Y live at root (written once in initDKIO).
    subroutine writeDKFileOutput(dkApp, nStep)
        class(dkApp_T), intent(inout) :: dkApp
        integer,        intent(in)    :: nStep

        character(len=strLen) :: gStr

        ! DK vars are ShellGridVar_T on the ghost-free DK grid: %data is (Nt,Np) (theta,phi).
        write(gStr,'(A,I0)') "Step#", nStep

        call ClearIO(IOVars)
        call AddOutVar(IOVars, DK_SIGP, dkApp%vars%SigmaP%data,   uStr="S")
        call AddOutVar(IOVars, DK_SIGH, dkApp%vars%SigmaH%data,   uStr="S")
        call AddOutVar(IOVars, DK_FAC,  dkApp%vars%FAC%data,      uStr="microA/m2")
        call AddOutVar(IOVars, DK_EAVG, dkApp%vars%Avg_Eng%data,  uStr="keV")
        call AddOutVar(IOVars, DK_NFLX, dkApp%vars%Num_Flux%data, uStr="1/cm2/s")
        call AddOutVar(IOVars, DK_POT,  dkApp%vars%Pot%data,      uStr="kV")
        call AddOutVar(IOVars, DK_ATYP, dkApp%vars%Aur_Type%data, uStr="")
        call AddOutVar(IOVars, DK_DELE, dkApp%vars%DeltaE%data,   uStr="kV")
        call AddOutVar(IOVars, DK_BETA, dkApp%aurora%beta_RCM,    uStr="")
        call AddOutVar(IOVars, DK_ZEAV, dkApp%vars%Z_Eavg%data,   uStr="keV")
        call AddOutVar(IOVars, DK_ZNFL, dkApp%vars%Z_Nflux%data,  uStr="1/cm2/s")
        call AddOutVar(IOVars, "time", dkApp%IO%time)
        call AddOutVar(IOVars, "MJD",  dkApp%IO%mjd)
        call WriteVars(IOVars, .false., h5File, gStr)
    end subroutine writeDKFileOutput

end module dkio
