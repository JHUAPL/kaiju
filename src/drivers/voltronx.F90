!Driver for Gamera coupled with Voltron (remix only)

program voltronx
    use clocks
    use gamapp
    use voltapp
    use voltio
    use uservoltic

    implicit none

    type(voltApp_T) :: vApp
    real(rp) :: nextDT
    logical :: doResetClocks = .false.

    call initClocks()

    vApp%vOptions%gamUserInitFunc => initUser

    call initVoltron(vApp)

    do while (vApp%time < vApp%tFin)
        !Start root timer
        call Tic("Omega", .true.)
        
        !Advance voltron models one coupling step
        nextDT = min(vApp%tFin-vApp%time, vApp%IO%nextIOTime()-vApp%time)

        call Tic("StepVoltron")
        call stepVoltron(vApp, nextDT)
        call Toc("StepVoltron")
        
        !IO checks
        call Tic("IO", .true.)
        !Console output
        if (vApp%IO%doConsole(vApp%time)) then
            !Using console output from Gamera
            call consoleOutputV(vApp,vApp%gApp)
            !Timing info
            if (vApp%IO%doTimerOut) call printClocks()
            doResetClocks = .true.
        elseif (vApp%IO%doTimer(vApp%time)) then
            if (vApp%IO%doTimerOut) call printClocks()
            doResetClocks = .true.
        endif
        
        !Data output
        if (vApp%IO%doOutput(vApp%time)) then
            call fOutputV(vApp,vApp%gApp)
        endif
        !Restart output
        if (vApp%IO%doRestart(vApp%time)) then
            call resOutputV(vApp,vApp%gApp)
        endif
        !Reset clocks last
        if (doResetClocks) then
            call cleanClocks()
            doResetClocks = .false.
        endif

        call Toc("IO", .true.)

        call Toc("Omega", .true.)

    end do
    write(*,*) "Fin"
    
end program voltronx

