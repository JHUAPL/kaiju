module chmpfields
    use chmpdefs
    use chmpunits
    use ebtypes
    use ebutils
    use plasmaputils
    use ioH5
    use math
    
    implicit none

    logical, parameter :: doRering = .true.
    integer, parameter :: NumRR = 1

    integer, parameter :: EBIOVARS = 15
    type(IOVAR_T), dimension(EBIOVARS) :: ebIOs

    logical :: chkStatic = .false. !Whether to check for static->time interval switch

    contains


    !Reads specific slice given by group string
    !File: bStr, Group: gStr
    !NOTE: Calculating convective electric field to avoid diffusive terms
    subroutine readEB(Model,ebState,ebGr,ebTab,ebF,gStr,doCalcLppO)
        type(chmpModel_T), intent(in)     :: Model
        type(  ebState_T), intent(inout)  :: ebState
        type(   ebGrid_T), intent(in)     :: ebGr
        type(    ioTab_T), intent(in)     :: ebTab
        type(  ebField_T), intent(inout)  :: ebF
        character(len=strLen), intent(in) :: gStr
        logical, optional, intent(in) :: doCalcLppO

        character(len=strLen) :: ebFile

        real(rp), dimension(:,:,:), allocatable :: Bx,By,Bz,Jx,Jy,Jz
        real(rp), dimension(:,:,:,:), allocatable :: Vx,Vy,Vz,D,P
        integer :: Nip,Njp,Nkp
        integer :: i,j,k,s,is,ie,js,je,ks,ke
        integer :: dN(NDIM)
        integer :: nioD,nioP,nioJx,nioJy,nioJz
        real(rp), dimension(NDIM) :: B0xyz, Vxyz,Bxyz,Exyz,xcc,Jxyz
        real(rp) :: jScl
        logical :: doCalcLpp
        character(len=strLen) :: dID,pID,vxID,vyID,vzID

        jScl = 1.0
        if (present(doCalcLppO)) then
            doCalcLpp = doCalcLppO
        else
            doCalcLpp = .true.
        endif


        write(*,'(5a)') '<Reading eb from ', trim(ebTab%bStr), '/', trim(gStr), '>'

        !------------
        !Init holders (global # of active cells)
        Nip = ebGr%Nip
        Njp = ebGr%Njp
        Nkp = ebGr%Nkp

        allocate(Bx(Nip,Njp,Nkp))
        allocate(By(Nip,Njp,Nkp))
        allocate(Bz(Nip,Njp,Nkp))

        allocate(Vx(Nip,Njp,Nkp,0:Model%nSpc))
        allocate(Vy(Nip,Njp,Nkp,0:Model%nSpc))
        allocate(Vz(Nip,Njp,Nkp,0:Model%nSpc))

        if (Model%doMHD) then
            allocate(D (Nip,Njp,Nkp,0:Model%nSpc))
            allocate(P (Nip,Njp,Nkp,0:Model%nSpc))
        endif

        if (Model%doJ) then
            allocate(Jx(Nip,Njp,Nkp))
            allocate(Jy(Nip,Njp,Nkp))
            allocate(Jz(Nip,Njp,Nkp))
        endif

        !------------
        !Get data from files
        dN = [ebTab%dNi,ebTab%dNj,ebTab%dNk] !Individual blocks
        do k=1,ebTab%Rk
            do j=1,ebTab%Rj
                do i=1,ebTab%Ri
                    ebFile = genName(ebTab%bStr,ebTab%Ri,ebTab%Rj,ebTab%Rk,i,j,k,Model%doOldNaming)

                    ! Get piece from file
                    call ClearIO(ebIOs)
                    ! Start with getting non-fluid and bulk variables
                    call AddInVar(ebIOs,"Bx")
                    call AddInVar(ebIOs,"By")
                    call AddInVar(ebIOs,"Bz")
                    call AddInVar(ebIOs,"Vx")
                    call AddInVar(ebIOs,"Vy")
                    call AddInVar(ebIOs,"Vz")
                    
                    if (Model%doMHD) then
                        call AddInVar(ebIOs,"D")
                        call AddInVar(ebIOs,"P")
                    endif
                    if (Model%doJ) then
                        call AddInVar(ebIOs,"Jx")
                        call AddInVar(ebIOs,"Jy")
                        call AddInVar(ebIOs,"Jz")
                    endif                            
                    call ReadVars(ebIOs,.true.,ebFile,gStr)

                    !Push piece to grid
                    is = (i-1)*ebTab%dNi + 1
                    js = (j-1)*ebTab%dNj + 1
                    ks = (k-1)*ebTab%dNk + 1
                    ie = is + ebTab%dNi - 1
                    je = js + ebTab%dNj - 1
                    ke = ks + ebTab%dNk - 1

                    Bx(is:ie,js:je,ks:ke) = inBScl*reshape(ebIOs(1)%data,dN)
                    By(is:ie,js:je,ks:ke) = inBScl*reshape(ebIOs(2)%data,dN)
                    Bz(is:ie,js:je,ks:ke) = inBScl*reshape(ebIOs(3)%data,dN)
                    Vx(is:ie,js:je,ks:ke,BLK) = inVScl*reshape(ebIOs(4)%data,dN)
                    Vy(is:ie,js:je,ks:ke,BLK) = inVScl*reshape(ebIOs(5)%data,dN)
                    Vz(is:ie,js:je,ks:ke,BLK) = inVScl*reshape(ebIOs(6)%data,dN)
                    if (Model%doMHD) then
                        nioD = FindIO(ebIOs,"D")
                        nioP = FindIO(ebIOs,"P")

                        !Convert incoming data to [#/cc] and [nPa] using defined scaling params
                        D(is:ie,js:je,ks:ke,BLK) = inDScl*reshape(ebIOs(nioD)%data,dN)
                        P(is:ie,js:je,ks:ke,BLK) = inPScl*reshape(ebIOs(nioP)%data,dN)
                    endif
                    if (Model%doJ) then
                        nioJx = FindIO(ebIOs,"Jx")
                        nioJy = FindIO(ebIOs,"Jy")
                        nioJz = FindIO(ebIOs,"Jz")
                        Jx(is:ie,js:je,ks:ke) = reshape(ebIOs(nioJx)%data,dN)
                        Jy(is:ie,js:je,ks:ke) = reshape(ebIOs(nioJy)%data,dN)
                        Jz(is:ie,js:je,ks:ke) = reshape(ebIOs(nioJz)%data,dN)
                    endif

                    if (Model%nSpc > 0) then
                        do s=1,Model%nSpc
                            write(dID  ,'(A,I0)') "D"  , s
                            write(pID  ,'(A,I0)') "P"  , s
                            write(vxID ,'(A,I0)') "Vx" , s
                            write(vyID ,'(A,I0)') "Vy" , s
                            write(vzID ,'(A,I0)') "Vz" , s
                            call ClearIO(ebIOs)
                            call AddInVar(ebIOs,dID )
                            call AddInVar(ebIOs,pID )
                            call AddInVar(ebIOs,vxID)
                            call AddInVar(ebIOs,vyID)
                            call AddInVar(ebIOs,vzID)
                            call ReadVars(ebIOs,.true.,ebFile,gStr)
                            D (is:ie,js:je,ks:ke,s) = inDScl*reshape(ebIOs(FindIO(ebIOs,dID ))%data,dN)
                            P (is:ie,js:je,ks:ke,s) = inPScl*reshape(ebIOs(FindIO(ebIOs,pID ))%data,dN)
                            Vx(is:ie,js:je,ks:ke,s) = inVScl*reshape(ebIOs(FindIO(ebIOs,vxID))%data,dN)
                            Vy(is:ie,js:je,ks:ke,s) = inVScl*reshape(ebIOs(FindIO(ebIOs,vyID))%data,dN)
                            Vz(is:ie,js:je,ks:ke,s) = inVScl*reshape(ebIOs(FindIO(ebIOs,vzID))%data,dN)
                        enddo
                    endif

                enddo
            enddo
        enddo

        !Do any rescaling (that needs info from file) necessary
        if (Model%doJ) then
            !Guarantee that Jxyz is SI current density [A/m2]
            nioJx = FindIO(ebIOs,"Jx")

            jScl = 1.0
            select case (trim(toUpper(ebIOs(nioJx)%unitStr)))
            case("NA/M2")
                !nA/m2 current density, standard for mspheres
                !Convert to typical SI A/m2
                jScl = (1.0e-9)
            case("CODE")
               if ( trim(toUpper(Model%uID))=="EARTH") then
                  ! note, using an ugly way to access B and L scaling by using in2G and L0 (globals from chmpunits)
                  ! FIXME: make chimp unit treatment more elegant, like gamera
                  jScl = in2G*1.e-5/(L0*Mu0)*1.0e-9   ! 1.e-5 to convert from G to nT and 1.e-9 to convert the result from nA to A
               else
                  !Not (yet) supported units
                  write(*,*) "------------------------"
                  write(*,*) "Error, units of current density are not [nA/m2] !"
                  write(*,*) "Units: ", trim(toUpper(ebIOs(nioJx)%unitStr))
                  write(*,*) "This is likely because the GAMERA simulation is too old."
                  write(*,*) "Either regenerate the MHD data or add the proper unit scaling."
                  write(*,*) "Womp womp womp ..."
                  write(*,*) "------------------------"
                  stop
             end if
            end select
        endif

        !Get time data from last file
        call ClearIO(ebIOs)
        call AddInVar(ebIOs,"time",vTypeO=IOREAL)
        call ReadVars(ebIOs,.true.,ebFile,gStr)
        i = FindIO(ebIOs,"time")
        ebF%time = inTScl*ebIOs(i)%data(1)

        !inclde step info
        ebF%gStr = gStr

        !------------
        !Pull data into ebField
        ebF%dB = 0.0
        ebF%E  = 0.0
        if (Model%doMHD) then
            ebF%W = 0.0
        endif


        !FIXME: Lazily assuming is=1,ie=Nxp
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,s,B0xyz,Bxyz,Vxyz,Exyz,Jxyz)    
        do k=1,Nkp
            do j=1,Njp
                do i=1,Nip
                    !Get residual field
                    B0xyz = ebGr%B0cc(i,j,k,:) !Background field

                    Bxyz = [Bx(i,j,k),By(i,j,k),Bz(i,j,k)]
                    Vxyz = [Vx(i,j,k,BLK),Vy(i,j,k,BLK),Vz(i,j,k,BLK)]

                    ebF%dB(i,j,k,:) = Bxyz-B0xyz

                    !Get electric field from V/B
                    Exyz = -cross(Vxyz,Bxyz)
                    ebF%E(i,j,k,:) = Exyz
                    if (Model%doMHD) then
                        ebF%W(i,j,k,DEN      ,BLK) = D(i,j,k,BLK)
                        ebF%W(i,j,k,VELX:VELZ,BLK) = Vxyz
                        ebF%W(i,j,k,PRESSURE ,BLK) = P(i,j,k,BLK)
                    endif
                    if (Model%doJ) then
                        Jxyz = [Jx(i,j,k),Jy(i,j,k),Jz(i,j,k)]
                        ebF%Jxyz(i,j,k,:) = jScl*Jxyz
                    endif

                    if (Model%nSpc>0) then
                        do s=1,Model%nSpc
                            Vxyz = [Vx(i,j,k,s),Vy(i,j,k,s),Vz(i,j,k,s)]
                            ebF%W(i,j,k,DEN      ,s) = D(i,j,k,s)
                            ebF%W(i,j,k,VELX:VELZ,s) = Vxyz
                            ebF%W(i,j,k,PRESSURE ,s) = P(i,j,k,s)
                        enddo
                    endif
                    
                enddo
            enddo
        enddo
        
        !------------
        !Fill in ghosts
        call ebGhosts(Model,ebGr,ebF)

        !------------
        !Calculate Lpp(MLT) if needed
        if (Model%doPP .and. doCalcLpp) call calcLppMLT(Model,ebState,ebF%time,ebF%Lpp)
        
        !------------
        !Clean up
        call ClearIO(ebIOs)
        deallocate(Bx,By,Bz,Vx,Vy,Vz)
        if (Model%doMHD) deallocate(D,P)

    end subroutine readEB

    !Fill in ghost cells for EB data
    subroutine ebGhosts(Model,ebGr,ebF)
        type(chmpModel_T), intent(in)    :: Model
        type(   ebGrid_T), intent(in)    :: ebGr
        type(  ebField_T), intent(inout) :: ebF

        integer :: n,d,i,j,k,ip,jp,kp,s
        integer :: Nk,Nk2
        logical :: iGh,jGh,kGh

        select case(ebGr%GrID)
        case(LFMGRID,EGGGRID)
            if (doRering) then
                !Average ring
                !NumRR = # of rings to average over
                Nk = ebGr%Nkp
                do n=0,NumRR-1
                    !$OMP PARALLEL DO default(shared) &
                    !$OMP private(i,d)
                    do i=ebGr%is,ebGr%ie
                        !Fields
                        do d=1,NDIM
                            ebF%dB(i,ebGr%js+n,ebGr%ks:ebGr%ke,d) = sum(ebF%dB(i,ebGr%js+n,ebGr%ks:ebGr%ke,d))/Nk
                            ebF%E (i,ebGr%js+n,ebGr%ks:ebGr%ke,d) = sum(ebF%E (i,ebGr%js+n,ebGr%ks:ebGr%ke,d))/Nk

                            ebF%dB(i,ebGr%je-n,ebGr%ks:ebGr%ke,d) = sum(ebF%dB(i,ebGr%je-n,ebGr%ks:ebGr%ke,d))/Nk
                            ebF%E (i,ebGr%je-n,ebGr%ks:ebGr%ke,d) = sum(ebF%E (i,ebGr%je-n,ebGr%ks:ebGr%ke,d))/Nk
                        enddo
                        !MHD variables
                        if (Model%doMHD) then
                            do s=0,Model%nSpc
                                do d=1,NVARMHD
                                    ebF%W (i,ebGr%js+n,ebGr%ks:ebGr%ke,d,s) = sum(ebF%W (i,ebGr%js+n,ebGr%ks:ebGr%ke,d,s))/Nk
                                    ebF%W (i,ebGr%je-n,ebGr%ks:ebGr%ke,d,s) = sum(ebF%W (i,ebGr%je-n,ebGr%ks:ebGr%ke,d,s))/Nk
                                enddo
                            enddo !spc
                        endif
                    enddo !i-loop
                enddo
            endif

            !Now do ghosts, lazily loop through *EVERYTHING*
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,ip,jp,kp,iGh,jGh,kGh)
            do k=ebGr%ksg,ebGr%keg
                do j=ebGr%jsg,ebGr%jeg
                    do i=ebGr%isg,ebGr%ieg
                        iGh = (i<ebGr%is) .or. (i>ebGr%ie)
                        jGh = (j<ebGr%js) .or. (j>ebGr%je)
                        kGh = (k<ebGr%ks) .or. (k>ebGr%ke)
                        if (iGh .or. jGh .or. kGh) then
                            call ijk2Active(Model,ebGr,i,j,k,ip,jp,kp)
                            ebF%dB(i,j,k,:) = ebF%dB(ip,jp,kp,:)
                            ebF%E (i,j,k,:) = ebF%E (ip,jp,kp,:)
                            if (Model%doMHD) then
                                ebF%W (i,j,k,:,:) = ebF%W (ip,jp,kp,:,:)
                            endif
                        endif
                    enddo
                enddo
            enddo
                    
        case(CARTGRID)
            write(*,*) 'Cartesian not implemented'
            stop
        end select

    end subroutine ebGhosts

end module chmpfields
