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

    public :: initDKIO, writeDKRestart, readDKRestart, writeDKFileOutput

    contains

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

        ! Restore the EMA history (SigmaP/SigmaH) before the first DKcompute.
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
        call AddOutVar(IOVars, "X", dkApp%grid%x(1:Np,1:Nt), uStr="Ri")
        call AddOutVar(IOVars, "Y", dkApp%grid%y(1:Np,1:Nt), uStr="Ri")
        call AddOutVar(IOVars, "UnitsID", "DragonKing")
        call WriteVars(IOVars, .false., h5File)
    end subroutine initDKIO

    ! Persist the DK conductance EMA state to <RunID>.dk.Res.<nRes>.h5.
    subroutine writeDKRestart(dkApp, nRes)
        class(dkApp_T), intent(inout) :: dkApp
        integer,        intent(in)    :: nRes

        character(len=strLen) :: h5Str, lnStr, nStr

        if (nRes == -1) then
            write(*,*) 'nRes is -1 in writeDKRestart, now die'
            stop
        endif
        write(nStr,'(I0.5)') nRes
        write(h5Str,'(A,A,A,A)') trim(h5RunID), '.dk.Res.', trim(nStr), '.h5'

        call CheckAndKill(h5Str)

        call ClearIO(IOVars)
        call AddOutVar(IOVars, "Pedersen conductance", dkApp%vars%SigmaP, uStr="S")
        call AddOutVar(IOVars, "Hall conductance", dkApp%vars%SigmaH, uStr="S")
        call AddOutVar(IOVars, "nRes", nRes)
        call AddOutVar(IOVars, "time", dkApp%IO%time)
        call AddOutVar(IOVars, "MJD",  dkApp%IO%mjd)
        call WriteVars(IOVars, .false., h5Str)

        ! Rolling "latest" symlink so a restart with nRes=-1 picks up the most
        ! recent DK restart (mirrors mix/gamera/voltron; respects doSymLinks).
        write(lnStr,'(A,A,A,A)') trim(h5RunID), '.dk.Res.', 'XXXXX', '.h5'
        call MapSymLink(h5Str, lnStr)
    end subroutine writeDKRestart


    ! Restore SigmaP/SigmaH from <resId>.dk.Res.<nRes>.h5. If the file is
    ! absent (e.g. restarting from a run that predates native DK restart),
    ! leave vars as initialized and warn; the EMA then cold-starts as before.
    subroutine readDKRestart(dkApp, resId, nRes)
        class(dkApp_T),   intent(inout) :: dkApp
        character(len=*), intent(in)    :: resId
        integer,          intent(in)    :: nRes

        character(len=strLen) :: h5Str, nStr
        logical :: fExist
        integer :: n0

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
        call AddInVar(IOVars, "Pedersen conductance")
        call AddInVar(IOVars, "Hall conductance")
        call ReadVars(IOVars, .false., h5Str)

        n0 = FindIO(IOVars, "Pedersen conductance")
        if (IOVars(n0)%isDone) call IOArray2DFill(IOVars, "Pedersen conductance", dkApp%vars%SigmaP)
        n0 = FindIO(IOVars, "Hall conductance")
        if (IOVars(n0)%isDone) call IOArray2DFill(IOVars, "Hall conductance", dkApp%vars%SigmaH)

        write(*,*) "Done readDKRestart"
    end subroutine readDKRestart


    ! Append the DK fields for this step as a Step#<nStep> group in the single
    ! <RunID>.dk.h5 file. Grid X/Y live at root (written once in initDKIO).
    subroutine writeDKFileOutput(dkApp, nStep)
        class(dkApp_T), intent(inout) :: dkApp
        integer,        intent(in)    :: nStep

        character(len=strLen) :: gStr
        integer :: Np, Nt

        ! DK runs on the active grid, so vars are all (Np,Nt); the explicit
        ! (1:Np,1:Nt) slices below are just defensive shape pins.
        Np = dkApp%grid%Np
        Nt = dkApp%grid%Nt
        write(gStr,'(A,I0)') "Step#", nStep

        call ClearIO(IOVars)
        call AddOutVar(IOVars, "Pedersen conductance",   dkApp%vars%SigmaP(1:Np,1:Nt),   uStr="S")
        call AddOutVar(IOVars, "Hall conductance",       dkApp%vars%SigmaH(1:Np,1:Nt),   uStr="S")
        call AddOutVar(IOVars, "Field-aligned current",  dkApp%vars%FAC(1:Np,1:Nt),      uStr="microA/m2")
        call AddOutVar(IOVars, "Average energy",         dkApp%vars%Avg_Eng(1:Np,1:Nt),  uStr="keV")
        call AddOutVar(IOVars, "Number flux",            dkApp%vars%Num_Flux(1:Np,1:Nt), uStr="1/cm2/s")
        call AddOutVar(IOVars, "Potential",              dkApp%vars%Pot(1:Np,1:Nt),      uStr="kV")
        call AddOutVar(IOVars, "Auroral type",           dkApp%vars%Aur_Type(1:Np,1:Nt), uStr="")
        call AddOutVar(IOVars, "time", dkApp%IO%time)
        call AddOutVar(IOVars, "MJD",  dkApp%IO%mjd)
        call WriteVars(IOVars, .false., h5File, gStr)
    end subroutine writeDKFileOutput

end module dkio
