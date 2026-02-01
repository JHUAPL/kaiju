!Standard initialization routine for eb fields
!Must specify initFields / updateFields
module userebic
    use chmpdefs
    use chmpunits
    use chmpdbz
    use ebtypes
    use xml_input
    use ioH5
    use chmpfields
    use ebinit
    use wpihelper
    
    implicit none

    contains

    !Standard EB init routine
    !Reads from gamera HDF5 data
    subroutine initFields(Model,ebState,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(inout)   :: ebState
        type(XML_Input_T), intent(inout) :: inpXML

        character(len=strLen) :: ebFile
        logical :: fExist
        integer :: Ni,Nj,Nk
        integer :: i1,i2

        logical :: doCalcLpp = .false.

        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        !Get info for input data
        call inpXML%Set_Val(ebTab%bStr,"fields/ebfile","ebdata.h5")
        call inpXML%Set_Val(ebTab%isMPI,"fields/isMPI",.false.)
        if (ebTab%isMPI) then
            call inpXML%Set_Val(ebTab%Ri,"parallel/Ri",1)
            call inpXML%Set_Val(ebTab%Rj,"parallel/Rj",1)
            call inpXML%Set_Val(ebTab%Rk,"parallel/Rk",1)
            call inpXML%Set_Val(Model%doOldNaming,"parallel/doOldNaming",.false.)
        else
            ebTab%Ri = 1
            ebTab%Rj = 1
            ebTab%Rk = 1
        endif
        ebFile = genName(ebTab%bStr,ebTab%Ri,ebTab%Rj,ebTab%Rk,1,1,1,Model%doOldNaming) !Get first (or only) eb file
        
        !Check file
        inquire(file=ebFile,exist=fExist)
        if (.not. fExist) then
            write(*,*) 'Error opening ebfile, exiting ...'
            write(*,*) 'ebfile = ', trim(ebFile)
            write(*,*) ''
            stop
        endif

        !Get time series/individual grid data from H5 file
        call InitIOTab(ebTab,inpXML,ebFile,inTScl,oTScl)

        ! Get number of MHD fluids
        call setNSpc(Model, ebTab)

        ! Check if Src data in file and save result to Model
        call CheckHasSrcData(Model, ebTab)

        !Start by getting grid from H5 file
        call rdGrid(Model,ebGr,ebTab,inpXML)

        !Allocate eb data (includes space for ghosts)
        call allocEB(Model,ebGr,eb1)
        call allocEB(Model,ebGr,eb2)

        !Do various types of initialization
        !ie, static, numb0, etc

        !Initialize and allocate plasmapause parameters
        if (Model%doPP) call initPP(ebState,inpXML)

        !Initialize wave particle interactions if present
        if ( Model%doWPI ) then
            call initWPI(ebState%ebWmodel,ebState%ebWave,inpXML)
        endif

    !Numerical background, put B0 field on grid
        if (Model%doNumB0) then
            write(*,*) '<Setting Numerical B0>'
        !Set eb1/eb2
            !Zero out E fields
            eb1%E = 0.0
            eb2%E = 0.0
            !Set dB's to B0 @ cell centers
            eb1%dB = ebGr%B0cc
            eb2%dB = ebGr%B0cc
            call ebGhosts(Model,ebGr,eb1)
            call ebGhosts(Model,ebGr,eb2)
            
        !Zero out background terms
            ebGr%B0cc = 0.0
            Model%B0 => NullB0
            Model%JacB0 => NullJacB0
            ebState%doStatic = .true.
            return
        endif
    !Pure background, only use analytic form
        if (Model%doPureB0) then
            !Zero out all fields
            write(*,*) '<Setting Pure B0>'
            eb1%E  = 0.0
            eb1%dB = 0.0
            eb2%E  = 0.0
            eb2%dB = 0.0
            ebState%doStatic = .true.
            return
        endif
        
        !Initialize eb data, find time slices
        call GetTabSlc(ebTab,Model%T0,i1,i2)
        call readEB(Model,ebState,ebGr,ebTab,eb1,ebTab%gStrs(i1),doCalcLpp)
        call readEB(Model,ebState,ebGr,ebTab,eb2,ebTab%gStrs(i2),doCalcLpp)
        if (ebTab%N == 1) ebState%doStatic = .true.
        end associate
    end subroutine initFields

    !Standard eb update
    subroutine updateFields(Model,ebState,t)
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T), intent(inout)   :: ebState 
        real(rp), intent(in) :: t

        logical :: skipUpdate
        integer :: i1,i2
        !Check if doStatic should be flipped
        if (chkStatic .and. (t >= ebState%eb1%time) ) then            
            ebState%doStatic = .false.
            chkStatic = .false.
        endif
        !write(*,*) 'Updating to ', t
        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        skipUpdate = (t >= ebState%eb1%time .and. t <= ebState%eb2%time) .or. ebState%doStatic
        if (t < ebState%eb1%time) then
            !Set to static until time is within interval
            ebState%doStatic = .true.
            skipUpdate = .true.
            chkStatic = .true.
        endif

        if (skipUpdate) return

        !Test for out of time
        if (t>=maxval(ebTab%times(:))) then
            write(*,*) 'Out of time data, switching to static fields ...'
            ebState%doStatic = .true.

            !Copy eb2->eb1
            eb1%gStr = eb2%gStr
            eb1%time = eb2%time
            eb1%dB   = eb2%dB
            eb1%E    = eb2%E
            if (Model%doMHD) eb1%W = eb2%W
            if (Model%doJ) eb1%Jxyz = eb2%Jxyz
            !Bail out
            return        
        endif

        !Do work if still here
        !Go ahead and reread both (lazy way of avoiding corner cases)
        call GetTabSlc(ebState%ebTab,t,i1,i2)

        ! Put in a trap so copies 2 to 1 and doesnt reread step
        if ( trim(eb2%gStr) == trim(ebTab%gStrs(i1)) ) then
            !Copy eb2->eb1
            eb1%gStr = eb2%gStr
            eb1%time = eb2%time
            eb1%dB   = eb2%dB
            eb1%E    = eb2%E
            if (Model%doMHD) eb1%W = eb2%W
            if (Model%doJ) eb1%Jxyz = eb2%Jxyz
            write(*,'(5a)') '<Copying eb1 from ', trim(ebTab%bStr), '/', trim(ebTab%gStrs(i1)), '>'
        else
            !Read eb1
            call readEB(Model,ebState,ebGr,ebTab,eb1,ebTab%gStrs(i1))
        end if

        !Read eb2
        call readEB(Model,ebState,ebGr,ebTab,eb2,ebTab%gStrs(i2))

        end associate

    end subroutine updateFields


    subroutine setNSpc(Model, ebTab)
        !! Set the number of species based on value in gamera h5 files
        !! Steps have the nSpc attribute we can read from
        type(chmpModel_T), intent(inout) :: Model
        type(ioTab_T), intent(in) :: ebTab
        
        character(len=strLen) :: ebFile
        type(IOVAR_T), dimension(2) :: IOVars

        ebFile = genName(ebTab%bStr,ebTab%Ri,ebTab%Rj,ebTab%Rk,1,1,1,Model%doOldNaming)

        call ClearIO(IOVars)
        call AddInVar(IOVars, 'nSpc', vTypeO=IOINT) 
        call ReadVars(IOVars, .false., ebFile,gStrO=ebTab%gStrs(1))
        
        Model%nSpc = GetIOInt(IOVars, 'nSpc')
        
    end subroutine setNSpc


    subroutine CheckHasSrcData(Model, ebTab)
        type(chmpModel_T), intent(inout) :: Model
        type(ioTab_T), intent(in) :: ebTab
        
        character(len=strLen) :: ebFile
        type(IOVAR_T), dimension(2) :: IOVars

        ebFile = genName(ebTab%bStr,ebTab%Ri,ebTab%Rj,ebTab%Rk,1,1,1,Model%doOldNaming)
        Model%doReadSrc = (ioExist(ebFile, "SrcX1", gStrO=ebTab%gStrs(1)))

    end subroutine CheckHasSrcData

end module userebic
