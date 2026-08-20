program raijuOWDx
    !! One way driving of RAIJU using GAMERA and REMIX output files

    use kdefs
    use planethelper
    use xml_input
    use clocks

    ! RAIJU stuff
    use raijudefs
    use raijustarter
    use raijuBCs
    use raijuPreAdvancer
    use raijuAdvancer

    ! Chimp stuff
    use ebtypes
    use starter
    use chmpfields
    use chopio
    use chmpunits

    ! Remix, actually
    !use calcdbtypes
    !use calcdbio
    use remixReader

    ! Voltron stuff
    use volttypes
    use voltAppHelper, only : initVoltState, calcCorotPotential
    use voltCplHelper, only : genVoltTubes
    use raijuuseric
    use raijuCplHelper

    implicit none

    type(voltApp_T)   :: vApp
        !! Kinda hacky, just create a voltApp object and initialize only the parts we're gonna use ourselves
    !Holder for remix data
    type(rmReader_T) :: rmReader
    type(raijuCoupler_T) :: raiCplApp

    character(len=strLen) :: XMLStr, gStr, ftag
    type(XML_Input_T) :: inpXML    
    
    character(len=strLen) :: FLH5, resId, swF
    
    logical :: doChmpOut,doFLOut, doColdStart

    real(rp) :: mjd0, Rin

    call initClocks()
    call Tic("Omega")

    ! Init xml
    call getIDeckStr(XMLStr)

    ! Init voltron stuff
    call raiOWD_initVoltron(vApp, XMLStr)

    ! Init CHIMP
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/Chimp",.true.)
    call goApe(vApp%ebTrcApp%ebModel,vApp%ebTrcApp%ebState,iXML=inpXML)
    
    allocate(imagOptions_T :: raiCplApp%opt)
    ! set all opt parameters for good measure
    call inpXML%Set_Val(raiCplApp%opt%swF, '/Kaiju/GAMERA/wind/tsfile', 'bcwind.h5')
    call inpXML%Set_Val(raiCplApp%opt%doColdStart,'/Kaiju/RAIJU/driver/doColdStart',.false.)
    associate(ebGr=>vApp%ebTrcApp%ebState%ebGr, ebState=>vApp%ebTrcApp%ebState)
        raiCplApp%opt%mjd0 = T2MJD(-1.0*ebState%ebTab%times(1)/inTscl, ebState%ebTab%MJDs(1))
        raicplApp%opt%mhdRin  = norm2(ebGr%xyzcc(ebGr%is ,ebGr%js,ebGr%ks,:))
        raicplApp%opt%mhdRinG = norm2(ebGr%xyzcc(ebGr%isg,ebGr%js,ebGr%ks,:))
        raiCplApp%opt%voltGrid = vApp%shGrid
    end associate

    ! Now we do raiju init
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/RAIJU",.true.)
    call inpXML%Set_Val(doChmpOut,'driver/doChmpOut',.false.)
    call inpXML%Set_Val(doFLOut,'driver/doFLOut',.false.)
    call raiCplApp%InitModel(inpXML)
    call raiCplApp%InitIO(inpXML)

    associate(raiApp=>raiCplApp%raiApp, ebModel=>vApp%ebTrcApp%ebModel, ebState=>vApp%ebTrcApp%ebState)
        
    call inpXML%Set_Val(raiApp%Model%doClockConsoleOut,'driver/doClockOut',.false.)

    if (raiApp%Model%isRestart) then
        call inpXML%Set_Val(resId,'restart/resId',raiApp%Model%RunID)
        call raiCplApp%ReadRestart(resId, raiApp%Model%nResIn)
        raiApp%State%isFirstCpl = .false.
    endif

    ! Init Remix reader
    call initRM("msphere", inpXML, rmReader)
    rmReader%time = raiApp%State%t
    vApp%time     = raiApp%State%t
    ! Also tell chimp when we're starting
    ebModel%t = inTscl*raiApp%State%t

    !call outputRMSG(rmReader,"rmReader.h5", .true.)
    
    ! Init outputs
    ebModel%doEBOut = doChmpOut
    if (ebModel%doEBOut) then
        call initEB3Dio(ebModel,ebState,inpXML)
    endif

    FLH5   = trim(raiApp%Model%RunID) // ".fl.h5" !RCM field lines
    call CheckAndKill(FLH5)

    ! Ready to loop
    do while ( raiApp%State%t < (raiApp%Model%tFin + 0.5) )

        call Tic("Output")
        ! Output if ready
        if (raiApp%State%IO%doConsole(raiApp%State%t)) then
            call raiCplApp%WriteConsoleOutput()
        endif

        if (raiApp%State%IO%doRestart(raiApp%State%t)) then
            call Tic("RAIJU Restart")
            call raiCplApp%WriteRestart(raiApp%State%IO%nRes)
            call Toc("RAIJU Restart")
        endif
        
        if (raiApp%State%IO%doOutput(raiApp%State%t)) then
            call Tic("RAIJU Output")
            call raiCplApp%WriteFileOutput(raiApp%State%IO%nOut)
            call Toc("RAIJU Output")

            if (ebModel%doEBOut) then
                ! Write eb at the same output cadence as imag
                write(gStr,'(A,I0)') "Step#", ebModel%nOut
                call Tic("CHIMP EB3D")
                call writeEB3D(ebModel,ebState,gStr)
                ebModel%nOut = ebModel%nOut + 1
                ebModel%tOut = inTscl*raiApp%State%IO%tOut  ! Idk if we need to set this since chimp isn't controlling its own output
                call Toc("CHIMP EB3D")
            endif

            if (doFLOut) then
                write(*,*)"No FL out"
            endif

            !write(gStr,'(A,I0)') "Step#", raiApp%State%IO%nOut-1  ! nOut got advanced by raijuOutput above
            !call Tic("Remix rmReader")
            !call outputRMSG(rmReader,"rmReader.h5",.false., gStr)
            !call Toc("Remix rmReader")
        endif
        call Toc("Output")

        ! Update other models
        call Tic("CHIMP update")
        call updateFields(ebModel, ebState, ebModel%t)
        call Toc("CHIMP update")

        call Tic("REMIX update")
        call updateRM(rmReader, rmReader%time)
        call Toc("REMIX update")

        ! Populate RAIJU's coupling variables with updated model info
        call Tic("fromV packing")
        call packRaijuCoupler_OWD(raiCplApp, vApp, rmReader)
        call Toc("fromV packing")

        call Tic("fromV to State")
        ! Now put coupler info into RAIJU's State
        call raiCplApp%toIMAG(vApp)
        call Toc("fromV to State")

        ! Step RAIJU
        call Tic("RAIJU Advance")
        call raiCplApp%AdvanceModel(raiApp%Model%dt)
        call Toc("RAIJU Advance")


        ! Advance model times
        raiApp%State%t  = raiApp%State%t  + raiApp%Model%dt
        raiApp%State%ts = raiApp%State%ts + 1
        raiApp%State%mjd = T2MJD(raiApp%State%t,raiCplApp%opt%mjd0)

        ebModel%t     = inTscl*raiApp%State%t
        ebModel%ts    = ebModel%ts + 1

        rmReader%time  = raiApp%State%t
        vApp%time      = raiApp%State%t

        if (raiApp%Model%doClockConsoleOut) then
            call printClocks()
            call cleanClocks()
        endif
    enddo

    end associate

    call Toc("Omega")
    write(*,*)"Main Objective complete!"

    contains

    subroutine raiOWD_initVoltron(vApp, xmlStr, nResO)
        type(voltApp_T), intent(inout) :: vApp
        character(len=*), intent(in) :: xmlStr
        integer, optional :: nResO

        type(XML_Input_T) :: xmlInp
        integer :: nRes
        character(len=strLen) :: resID, nStr, ResF

        if (present(nResO)) then
            nRes = nResO
        else
            nRes = 0
        endif

        xmlInp = New_XML_Input(trim(XMLStr),"Kaiju/Voltron",.true.)

        call xmlInp%Set_Val(resID,"/Kaiju/gamera/restart/resID","msphere")
        call xmlInp%Set_Val(nRes,"/Kaiju/gamera/restart/nRes" ,-1)
        !Get number string
        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') nRes
        endif
        write (ResF, '(A,A,A,A)') trim(resID), ".volt.Res.", trim(nStr), ".h5"
        call GenShellGridFromFile(vApp%shGrid, "VOLTRON", ResF)
        call getPlanetParams(vApp%planet, xmlInp)

        call initVoltState(vApp)

    end subroutine raiOWD_initVoltron

    subroutine packRaijuCoupler_OWD(raiCpl, vApp, rmReader)
        class(raijuCoupler_T), intent(inout) :: raiCpl
        type(voltApp_T), intent(inout) :: vApp
        type(rmReader_T) :: rmReader

        ! Update coupling time
        raiCpl%tLastUpdate = raiCpl%raiApp%State%t

        ! Draw tubes from volt grid
        call genVoltTubes(vApp)
        ! Set potential
        call calcCorotPotential(vApp%planet, vApp%shGrid, vApp%State%potential_corot,doGeoCorotO=vApp%doGeoCorot)
        call InterpShellVar_TSC_SG(rmReader%shGr, rmReader%nsPot(1), vApp%shGrid, vApp%State%potential_total)
        vApp%State%potential_total%data = vApp%State%potential_total%data + vApp%State%potential_corot%data
        vApp%State%potential_total%mask = .true.
        
        call packRaijuCoupler_RT(raiCpl, vApp)

    end subroutine packRaijuCoupler_OWD

end program raijuOWDx
