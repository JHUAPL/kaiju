! functions to implement voltron coupling between gamera and other models

module gamCoupleHelper
    use gamtypes
    use volttypes
    use step
    use init
    use mhdgroup
    use output
    use cmiutils
    use uservoltic ! required to have IonInnerBC_T defined


    implicit none

    contains

    subroutine getCPCP(mhdvarsin,cpcp)
        real(rp), dimension(:,:,:,:,:),intent(in) :: mhdvarsin
        real(rp), intent(out) :: cpcp(2)
        cpcp(NORTH) = maxval(mhdvarsin(1,:,:,MHDPSI,NORTH))-minval(mhdvarsin(1,:,:,MHDPSI,NORTH))
        cpcp(SOUTH) = maxval(mhdvarsin(1,:,:,MHDPSI,SOUTH))-minval(mhdvarsin(1,:,:,MHDPSI,SOUTH))

    end subroutine getCPCP

    subroutine writeCouplerFileOutput(App, nStep)
        class(gamCoupler_T), intent(inout) :: App
        integer, intent(in) :: nStep

        integer :: i,j,k
        real(rp) :: cpcp(2)
        real(rp), dimension(:,:,:,:), allocatable, save :: inEijk,inExyz,Veb
        real(rp), dimension(:,:,:), allocatable, save :: psi
        real(rp), dimension(NDIM) :: Exyz,Bdip,xcc
        character(len=strLen) :: gStr
        type(IOVAR_T), dimension(50) :: IOVars

        associate(Gr => App%Grid)

        if (.not. allocated(inEijk)) allocate(inEijk(PsiSh+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,1:NDIM))
        if (.not. allocated(inExyz)) allocate(inExyz(PsiSh  ,Gr%jsg:Gr%jeg  ,Gr%ksg:Gr%keg  ,1:App%Model%nvInExyz))
        if (.not. allocated(psi))    allocate(psi(PsiSh,Gr%js:Gr%je,Gr%ks:Gr%ke)) !Cell-centered potential
        if (.not. allocated(Veb))    allocate(Veb(PsiSh,Gr%js:Gr%je,Gr%ks:Gr%ke,1:NDIM))

        call getCPCP(App%mixOutput,cpcp)
        call Ion2MHD(App%Model,App%Grid,App%gPsi,inEijk,inExyz,App%rm2g)

        !Subtract dipole before calculating current
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xcc,Bdip,Exyz)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=1,PsiSh
                    psi(i,j,k) = 0.125*( App%gPsi(i+1,j  ,k) + App%gPsi(i+1,j  ,k+1) &
                                       + App%gPsi(i  ,j+1,k) + App%gPsi(i  ,j+1,k+1) &
                                       + App%gPsi(i+1,j+1,k) + App%gPsi(i+1,j+1,k+1) &
                                       + App%gPsi(i  ,j  ,k) + App%gPsi(i  ,j  ,k+1) )
                    xcc = Gr%xyzcc(i,j,k,:)
                    Bdip = MagsphereDipole(xcc,App%Model%MagM0)
                    Exyz(IONO_EX) = inExyz(i,j,k,IONO_EX)
                    Exyz(IONO_EY) = inExyz(i,j,k,IONO_EY)
                    Exyz(IONO_EZ) = inExyz(i,j,k,IONO_EZ)
                    Veb(i,j,k,:) = cross(Exyz,Bdip)/dot_product(Bdip,Bdip)
                enddo
            enddo
        enddo

        call ClearIO(IOVars)
        call AddOutVar(IOVars,"Ex",inExyz(:,Gr%js:Gr%je,Gr%ks:Gr%ke,IONO_EX))
        call AddOutVar(IOVars,"Ey",inExyz(:,Gr%js:Gr%je,Gr%ks:Gr%ke,IONO_EY))
        call AddOutVar(IOVars,"Ez",inExyz(:,Gr%js:Gr%je,Gr%ks:Gr%ke,IONO_EZ))

        call AddOutVar(IOVars,"Vx",Veb(:,:,:,XDIR))
        call AddOutVar(IOVars,"Vy",Veb(:,:,:,YDIR))
        call AddOutVar(IOVars,"Vz",Veb(:,:,:,ZDIR))

        call AddOutVar(IOVars,"psi",psi)

        call AddOutVar(IOVars,"cpcpN",cpcp(1))
        call AddOutVar(IOVars,"cpcpS",cpcp(2))

        write(gStr,'(A,I0)') "Step#", nStep
        call WriteVars(IOVars,.true., App%vh5File,gStr)

        end associate

    end subroutine

    subroutine writeGamCouplerRestart(App, nRes)
        class(gamCoupler_T), intent(inout) :: App
        integer, intent(in) :: nRes

        character(len=strLen) :: ResF,lnResF
        integer, parameter :: MAXGCIOVAR = 20
        type(IOVAR_T), dimension(MAXGCIOVAR) :: IOVars

        !Restart Filename
        write (ResF, '(A,A,I0.5,A)') trim(App%Model%RunID), ".gamCpl.Res.", nRes, ".h5"

        !Reset IO chain
        call ClearIO(IOVars)

        !Remix Coupling Variables
        call AddOutVar(IOVars,"mixOutput",App%mixOutput)
        call AddOutVar(IOVars,"gPsi"     ,App%gPsi)
        call AddOutVar(IOVars,"Gas0"     ,App%Grid%Gas0)
        SELECT type(iiBC=>App%Grid%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                    call AddOutVar(IOVars,"inEijk",iiBC%inEijk)
                    call AddOutVar(IOVars,"inExyz",iiBC%inExyz)
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in writeGamCouplerRestart'
                stop
        END SELECT

        !Write out, force real precision
        call WriteVars(IOVars,.false.,ResF)

        !Create link to latest restart
        write (lnResF, '(A,A,A,A)') trim(App%Model%RunID), ".gamCpl.Res.", "XXXXX", ".h5"
        call MapSymLink(ResF,lnResF)

    end subroutine

    subroutine readGamCouplerRestart(App, resId, nRes)
        class(gamCoupler_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        character(len=strLen) :: ResF,nStr
        logical :: fExist
        integer, parameter :: MAXGCIOVAR = 20
        type(IOVAR_T), dimension(MAXGCIOVAR) :: IOVars
        integer :: n0

        !Restart Filename
        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') nRes
        endif

        write (ResF, '(A,A,A,A)') trim(resID), ".gamCpl.Res.", trim(nStr), ".h5"

        inquire(file=ResF,exist=fExist)
        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open input restart file for mix2gam, exiting'
            stop
        endif

        !Reset IO chain
        call ClearIO(IOVars)

        !Read Remix Coupling Variables
        call AddInVar(IOVars,"mixOutput")
        call AddInVar(IOVars,"gPsi")
        call AddInVar(IOVars,"Gas0")
        call AddInVar(IOVars,"inEijk")
        call AddInVar(IOVars,"inExyz")

        !Get data
        call ReadVars(IOVars,.false.,ResF)

        call IOArray5DFill(IOVars,"mixOutput",App%mixOutput)
        call IOArray3DFill(IOVars,"gPsi",App%gPsi)

        n0 = FindIO(IOVars,"Gas0")
        if (IOVars(n0)%isDone) then
            ! restart is newer and has full data
            call IOArray4DFill(IOVars,"Gas0",App%Grid%Gas0)
            SELECT type(iiBC=>App%Grid%externalBCs(INI)%p)
                TYPE IS (IonInnerBC_T)
                        call IOArray4DFill(IOVars,"inEijk",iiBC%inEijk)
                        call IOArray4DFill(IOVars,"inExyz",iiBC%inExyz)
                CLASS DEFAULT
                    write(*,*) 'Could not find Ion Inner BC in readGamCouplerRestart'
                    stop
            END SELECT
        endif

    end subroutine

end module
