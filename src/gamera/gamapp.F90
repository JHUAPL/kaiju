! implementation of subroutines in GameraApp type

module gamapp
    use gamtypes
    use step
    use init
    use mhdgroup
    use output

    implicit none

    contains

    subroutine initGamIO(Model, xmlInp)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: xmlInp

        logical :: doFatIO

        !Output/Restart (IOCLOCK)
        call Model%IO%init(xmlInp,Model%t,Model%ts,Model%isRestart)
        call xmlInp%Set_Val(Model%doDivB ,'output/DivB'    ,.true. )
        call xmlInp%Set_Val(doFatIO      ,'output/doFatIO' ,.false.)
        if (doFatIO) then
            call SetFatIO()
        endif
    end subroutine

    subroutine initGamera(gameraApp, xmlInp)
        class(gamApp_T), intent(inout) :: gameraApp
        type(XML_Input_T), intent(inout) :: xmlInp

        character(len=strLen) :: kaijuRoot

        gameraApp%Model%isLoud = .true.

        call xmlInp%SetRootStr('Kaiju/Gamera')
        call xmlInp%SetVerbose(.true.)

        ! try to verify that the XML file has "Kaiju" as a root element
        kaijuRoot = ""
        call xmlInp%Get_Key_Val("/Gamera/sim/H5Grid",kaijuRoot,.false.)
        if(len(trim(kaijuRoot)) /= 0) then
            write(*,*) "The input XML appears to be of an old style."
            write(*,*) "As of June 12th, 2021 it needs a root element of <Kaiju>."
            write(*,*) "Please modify your XML config file by adding this line at the top:"
            write(*,*) "<Kaiju>"
            write(*,*) "and this line at the bottom:"
            write(*,*) "</Kaiju>"
            write(*,*) "OR (preferred) convert your configuration to an INI file and use"
            write(*,*) " the XMLGenerator.py script to create conforming XML files."
            write(*,*) "Please refer to the python script or"
            write(*,*) " the [Generating XML Files] wiki page for additional info."
            stop
        endif

        ! read debug flags
        call xmlInp%Set_Val(writeGhosts,"debug/writeGhosts",.false.)
        call xmlInp%Set_Val(writeMagFlux,"debug/writeMagFlux",.false.)

        !Initialize Grid/State/Model (Hatch Gamera)
        !Will enforce 1st BCs, caculate 1st timestep, set oldState
        call Hatch(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Solver,xmlInp,gameraApp%gOptions%userInitFunc)
        call cleanClocks()

        if (.not. gameraApp%Model%isSub) then
            if (.not. gameraApp%Model%isRestart) call fOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
            call consoleOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        
    end subroutine initGamera

    subroutine stepGamera(gameraApp)
        class(gamApp_T), intent(inout) :: gameraApp

        !update the state variables to the next timestep
        call UpdateStateData(gameraApp)

        !Calculate new timestep
        call Tic("DT")
        gameraApp%Model%dt = CalcDT(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call Toc("DT")

        !Enforce BCs
        call Tic("BCs")
        call EnforceBCs(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        !Update Bxyz's
        call bFlux2Fld (gameraApp%Model,gameraApp%Grid,gameraApp%State%magFlux,gameraApp%State%Bxyz) 
        call Toc("BCs")
        
    end subroutine stepGamera

    subroutine UpdateStateData(gameraApp)
        class(gamApp_T), intent(inout) :: gameraApp

        call Tic("Gamera",.true.)
        !Advance system
        call AdvanceMHD(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Solver,gameraApp%Model%dt)
        call Toc("Gamera",.true.)

        !Call user-defined per-step function
        !NOTE: Do this before updating time
        if (associated(gameraApp%Model%HackStep)) then
            call Tic("HackStep")
            call gameraApp%Model%HackStep(gameraApp%Model,gameraApp%Grid,gameraApp%State)
            call Toc("HackStep")
        endif

        !Update info
        gameraApp%Model%ts = gameraApp%Model%ts+ 1
        gameraApp%Model%t  = gameraApp%Model%t + gameraApp%Model%dt

    end subroutine UpdateStateData

end module

