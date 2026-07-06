module voltCplHelper

    use kdefs
    use voltTypes
    use voltCplTypes
    use tubehelper
    use streamline,   only: FLIonFAC, FLIonDP
    use msphutils,    only: RadIonosphere
    use gdefs,        only: BLK

    implicit none

    contains

    subroutine genVoltTubes(vApp)
        !! Update Tube_Ts on Voltron's grid
        class(voltApp_T), intent(inout) :: vApp

        ! calculate field lines
        call calcTubes(vApp)

        ! Now pack into tubeShell
        call tubes2Shell(vApp%shGrid, vApp%State%ijTubes, vApp%State%tubeShell)

    end subroutine genVoltTubes

    subroutine genVoltTubesSpinup(vApp)
        !! Spinup stand-in for genVoltTubes: fill J11i/bcN/bcP by seed projection; deep products stay zero
        class(voltApp_T), intent(inout) :: vApp

        type(magLine_T) :: bTrc
        real(rp) :: mhd_Rin, seedR, Rion, eqR, bD, bP
        integer  :: i, j, s
        integer, parameter :: i0 = 2   ! sampling shell index (matches Line2Tube)

        associate(sh      => vApp%shGrid,                 &
                  ebModel => vApp%ebTrcApp%ebModel,       &
                  ebGr    => vApp%ebTrcApp%ebState%ebGr,  &
                  ebState => vApp%ebTrcApp%ebState,       &
                  ts      => vApp%State%tubeShell)

            ! Same no-data cutoff as calcTubes: dipole apex below the MHD inner corner
            mhd_Rin = norm2(ebGr%xyz(ebGr%is+2, ebGr%js, ebGr%ks, :))
            seedR   = sh%radius
            Rion    = RadIonosphere()

            ts%J11i%data = 0.0_rp
            do s = 0, ebModel%nSpc
                ts%bcN(s)%data = 0.0_rp
                ts%bcP(s)%data = 0.0_rp
            enddo

            !$OMP PARALLEL DO default(shared) schedule(dynamic) &
            !$OMP private(i,j,s,eqR,bTrc,bD,bP)
            do j = sh%js, sh%je+1
                do i = sh%is, sh%ie+1
                    eqR = DipColat2L(sh%thRp(i))
                    if (eqR <= mhd_Rin) cycle   ! dipole region: no MHD sample (as DipoleTube)

                    ! Seed-projection sample: FLIonFAC/FLIonDP use only the seed point
                    bTrc%x0 = seedR*[ sin(sh%th(i))*cos(sh%ph(j)), &
                                      sin(sh%th(i))*sin(sh%ph(j)), &
                                      cos(sh%th(i)) ]
                    ts%J11i%data(i,j) = FLIonFAC(ebModel, ebGr, ebState, bTrc, Rion, i0)
                    do s = 0, ebModel%nSpc
                        call FLIonDP(ebModel, ebGr, ebState, bTrc, bD, bP, i0, s)
                        ts%bcN(s)%data(i,j) = bD
                        ts%bcP(s)%data(i,j) = bP
                    enddo
                enddo
            enddo
            !$OMP END PARALLEL DO

            if (mod(nint(vApp%time), 60) == 0) then
                write(*,'(a,f9.1,a,es10.3,a,es10.3)') ' <spinupTubes> t=', vApp%time, &
                    '  max|J11i|=', maxval(abs(ts%J11i%data)), '  max bcN=', maxval(ts%bcN(BLK)%data)
            endif

        end associate
    end subroutine genVoltTubesSpinup

    subroutine calcTubes(vApp)
        !! Calculate field lines in ijTubes
        class(voltApp_T), intent(inout) :: vApp

        integer :: i,j,iL,iR,i1,i2
        real(rp) :: seedR, eqR, mhd_Rin,TbFill
        real(rp), dimension(NDIM) :: xyz0
        logical :: doSH,doNH
        logical, dimension(:,:), allocatable :: isDipT
        type(magLine_T) :: magLine
        ! Dipole region filling stuff
        logical :: foundDipRun, haveLeft, haveRight
        integer :: nDipRuns
        real(rp) :: TbLeft, TbRight

        associate(sh=>vApp%shGrid, ebApp=>vApp%ebTrcApp, ebGr=>vApp%ebTrcApp%ebState%ebGr)
            mhd_Rin = norm2(ebGr%xyz(ebGr%is+2,ebGr%js,ebGr%ks,:))
            seedR = sh%radius  ! Ionosphere radius in Rp

            allocate(isDipT(sh%is:sh%ie+1,sh%js:sh%je+1))

            ! Do field line tracing, populate ijTubes
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,eqR,magLine,doSH,doNH,xyz0)
            do j=sh%js,sh%je+1
                do i=sh%is,sh%ie+1
                    !Calculate seed point

                    xyz0 = seedR*[sin(sh%th(i))*cos(sh%ph(j)), &
                                  sin(sh%th(i))*sin(sh%ph(j)), &
                                  cos(sh%th(i))]
                    eqR = DipColat2L(sh%thRp(i))  ! Function assumes colat coming from 1 Rp, make sure we use the right theta value
                    if (eqR .le. mhd_Rin) then
                        !No MHD to tube from
                        call DipoleTube(vApp%planet,xyz0,vApp%State%ijTubes(i,j))
                        isDipT(i,j) = .true.
                    else
                        isDipT(i,j) = .false.
                        xyz0 = DipoleShift(xyz0, seedR + TINY)
                        if (xyz0(ZDIR) < 0) then
                            doNH = .false.
                            doSH = .true.
                        else
                            doNH = .true.
                            doSH = .false.
                        endif

                        call CleanLine(magLine)
                        !Note: Not using volt time b/c chimp wants time in its units
                        call genLine(ebApp%ebModel,ebApp%ebState,xyz0,ebApp%ebState%eb1%time, magLine,&
                                     doNHO=doNH,doSHO=doSH)
                        call Line2Tube(ebApp,vApp%planet,magLine,vApp%State%ijTubes(i,j))
                    endif

                enddo
            enddo

        ! Go back and patch Tb for dipole-only low-lat values
            do j=sh%js,sh%je+1
                i = sh%is
                foundDipRun = .false.
                nDipRuns = 0

                do while (i <= sh%ie+1)

                    if (.not. isDipT(i,j)) then
                        i = i + 1
                        cycle
                    endif

                    foundDipRun = .true.
                    nDipRuns = nDipRuns + 1

                    i1 = i
                    do while (i <= sh%ie+1 .and. isDipT(i,j))
                        i = i + 1
                    enddo
                    i2 = i - 1

                    haveLeft  = (i1 > sh%is)
                    haveRight = (i2 < sh%ie+1)

                    if (.not. haveLeft .and. .not. haveRight) then
                        write(*,*) "ERROR in calcTubes: dipole segment spans full i-range at j=", j
                        stop
                    endif

                    if (haveLeft)  TbLeft  = vApp%State%ijTubes(i1-1,j)%Tb
                    if (haveRight) TbRight = vApp%State%ijTubes(i2+1,j)%Tb

                    if (haveLeft .and. haveRight) then
                        TbFill = 0.5_rp*(TbLeft + TbRight)
                    else if (haveLeft) then
                        TbFill = TbLeft
                    else
                        TbFill = TbRight
                    endif

                    if (TbFill <= 0.0_rp) then
                        write(*,*) "ERROR in calcTubes: invalid TbFill for dipole segment at j=", j, " i1/i2=", i1, i2
                        write(*,*) "  TbLeft/TbRight=", TbLeft, TbRight
                        stop
                    endif

                    do i=i1,i2
                        vApp%State%ijTubes(i,j)%Tb = TbFill
                    enddo

                enddo

                if (.not. foundDipRun) then
                    write(*,*) "ERROR in calcTubes: no dipole tube segment found at j=", j
                    stop
                endif

                if (nDipRuns > 1) then
                    write(*,*) "WARNING in calcTubes: multiple dipole segments found at j=", j, " count=", nDipRuns
                endif
            enddo

        end associate
        
    end subroutine calcTubes

end module voltCplHelper
