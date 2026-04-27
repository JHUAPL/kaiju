module sliceio

    use chmpdefs
    use chmpunits
    use ebtypes
    use streamline
    use ioH5
    use xml_input
    use files
    use plasmaputils
    use parintime

    use volttypes

    implicit none

    character(len=strLen) :: ebOutF
    integer, parameter :: MAXEBVS = 80
    !Parameters for output EB grid (2D equatorial)
    !Fixme: Clean up and generalize
    integer  :: Nx1 = 128, Nx2 = 256
    real(rp) :: xSun = 12.5,xTail=-20.0,yM=20.0 !Default Cartesian slice
    real(rp) :: dx0=0.05
    
    real(rp), dimension(:,:), allocatable :: xxi,yyi,xxc,yyc
    real(rp), dimension(:,:,:), allocatable :: B02D
    logical , private :: doXY = .true. !Do XY or XZ slice
    real(rp), private :: z0=0.0 !Default XY slice
    real(rp), private :: y0=0.0 !Default XZ slice
    real(rp), private :: dSlc=0.05 !Spacing for slice averaging

    contains

    subroutine initEBio(Model,ebState,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(XML_Input_T), intent(in) :: inpXML

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars
        real(rp) :: dx1,dx2,Rin,Rout,Pin,Pout,x1,x2,xMax
        integer :: i,j,ijk(NDIM)
        integer :: iS,iE,n,Npow
        real(rp) :: xcc(NDIM)
        real(rp), dimension(:,:,:), allocatable :: ijkXY
        character(len=strLen) :: idStr
        logical :: doLog

        !write(ebOutF,'(2a)') trim(adjustl(Model%RunID)),'.eb.h5'
        ! Check for time parallelism
        call InitParInTime(Model,inpXML,"eb",ebOutF)

        associate( ebGr=>ebState%ebGr )

        call CheckAndKill(ebOutF)

        !Figure out grid slice type
        call inpXML%Set_Val(doXY,'slice/doXY',.true.)
        call inpXML%Set_Val(idStr,'slice/grType',"XY")
        call inpXML%Set_Val(Npow,'slice/Npow',0) !Number of times to double
        call inpXML%Set_Val(dSlc,'slice/dSlc',dSlc)

        if (doXY) then
            call inpXML%Set_Val(z0,'slice/z0',z0)
        else
            call inpXML%Set_Val(y0,'slice/y0',y0)
        endif

        select case(trim(toUpper(idStr)))
            case("XY")
                call inpXML%Set_Val(xSun,'slice/xSun',xSun)
                call inpXML%Set_Val(xTail,'slice/xTail',xTail)
                call inpXML%Set_Val(yM,'slice/yM',yM)
                Nx1 = nint((xSun-xTail)/dx0)
                Nx2 = nint(2*yM/dx0)

                call inpXML%Set_Val(Nx1,'slice/Nx1',Nx1)
                call inpXML%Set_Val(Nx2,'slice/Nx2',Nx2)
            case("RP")
                Nx1 = ebState%ebGr%Nip
                Nx2 = 2*ebState%ebGr%Njp
                call inpXML%Set_Val(Nx1,'slice/Nx1',Nx1)
                call inpXML%Set_Val(Nx2,'slice/Nx2',Nx2)
                !Do polar grid, take Rin/Rout from sunward line
                Rin  = ebGr%xyz(ebGr%is+1,ebGr%js,ebGr%ks,XDIR)
                Rout = ebGr%xyz(ebGr%ie-1,ebGr%js,ebGr%ks,XDIR)

                call inpXML%Set_Val(Rin,'slice/Rin',Rin)
                call inpXML%Set_Val(Rout,'slice/Rout',Rout)
                write(*,*) 'Radial grid bounds = ', Rin,Rout
                call inpXML%Set_Val(doLog,'slice/doLog',.true.)
                !Get phi bounds (deg)
                call inpXML%Set_Val(Pin ,'slice/Pin' ,0.0  )
                call inpXML%Set_Val(Pout,'slice/Pout',360.0)

            case("LFM2D")
                !Create 2D LFM slice (w/ full 2pi)
                call inpXML%Set_Val(xSun,'slice/xSun',xSun)
                
                !Make sure xSun doesn't run over the grid
                xMax = ebGr%xyz(ebGr%ie,ebGr%js,ebGr%ks,XDIR) - TINY
                xSun = min(xSun,xMax)

                !Find first cell that's "in" according to set locator
                i = ebGr%is
                xcc = ebGr%xyz(i,ebGr%js,ebGr%ks,:)
                do while ( .not. inDomain(xcc,Model,ebGr) )
                    i = i+1
                    xcc = ebGr%xyz(i,ebGr%js,ebGr%ks,:)
                enddo
                iS = min(i,ebGr%ie)
                iE = minloc(abs(ebGr%xyz(:,ebGr%js,ebGr%ks,XDIR)-xSun),dim=1)
                write(*,*) 'i-Start/End/Total = ',iS,iE,ebGr%Nip
                
                Nx1 = iE-iS+1
                Nx2 = 2*ebState%ebGr%Njp
        end select

        !Create 2D grid, wait to allocate rest until after embiggening
        allocate(xxi(Nx1+1,Nx2+1))
        allocate(yyi(Nx1+1,Nx2+1))

        !Create grid
        select case(trim(toUpper(idStr)))
        !---------
        case("XY")
            
            dx1 = (xSun-xTail)/Nx1
            dx2 = 2*yM/Nx2
            do j=1,Nx2+1
                do i=1,Nx1+1
                    x1 = xTail + (i-1)*dx1
                    x2 = -yM + (j-1)*dx2
                    xxi(i,j) = x1
                    yyi(i,j) = x2

                enddo
            enddo
        !---------
        case("RP")
            if (doLog) then 
                dx1 = ( log10(Rout)-log10(Rin) )/Nx1
            else
                dx1 = (Rout-Rin)/Nx1
            endif
            dx2 = (Pout-Pin)/Nx2
            do j=1,Nx2+1
                do i=1,Nx1+1
                    if (doLog) then
                        x1 = 10**( log10(Rin) + (i-1)*dx1 )
                    else
                        x1 = Rin + (i-1)*dx1
                    endif

                    x2 = Pin + (j-1)*dx2 !Degrees
                    xxi(i,j) = x1*cos(x2*PI/180.0)
                    yyi(i,j) = x1*sin(x2*PI/180.0)
                enddo
            enddo
        !---------
        case("LFM2D")
            !Do upper half plane
            do j=1,(Nx2/2)+1
                do i=1,Nx1+1
                    xxi(i,j) = ebGr%xyz(i+iS-1,j,ebGr%ks,XDIR)
                    yyi(i,j) = ebGr%xyz(i+iS-1,j,ebGr%ks,YDIR)
                enddo
            enddo
            !Reflect for bottom
            do j=1,Nx2/2
                do i=1,Nx1+1
                    xxi(i,j+1+Nx2/2) =  ebGr%xyz(i+iS-1,Nx2/2 - j+1,ebGr%ks,XDIR)
                    yyi(i,j+1+Nx2/2) = -ebGr%xyz(i+iS-1,Nx2/2 - j+1,ebGr%ks,YDIR)
                enddo
            enddo
        end select

        !Embiggen grid (i.e. upscale)
        do n=1,Npow
            write(*,*) 'Embiggening grid ...'
            write(*,*) 'Before: ',Nx1,Nx2
            call Embiggen(xxi,yyi,Nx1,Nx2)
            write(*,*) 'After : ',Nx1,Nx2
        enddo

        !Allocate remaining holders
        allocate(xxc(Nx1,Nx2),yyc(Nx1,Nx2))
        !Hold mappings from cell center of 2D grid to original XYZ grid cells
        allocate(ijkXY(Nx1,Nx2,NDIM))

        !Calculate cell-centers/background field
        allocate(B02D(Nx1,Nx2,NDIM))
        do j=1,Nx2
            do i=1,Nx1
                xxc(i,j) = 0.25*( xxi(i,j) + xxi(i+1,j) + xxi(i,j+1) + xxi(i+1,j+1) )
                yyc(i,j) = 0.25*( yyi(i,j) + yyi(i+1,j) + yyi(i,j+1) + yyi(i+1,j+1) )
                if (doXY) then
                    xcc(XDIR) = xxc(i,j)
                    xcc(YDIR) = yyc(i,j)
                    xcc(ZDIR) = z0
                else
                    !Assuming XZ, x->x and y->z
                    xcc(XDIR) = xxc(i,j)
                    xcc(YDIR) = y0
                    xcc(ZDIR) = yyc(i,j)
                endif
                
                if (inDomain(xcc,Model,ebGr)) then
                    B02D(i,j,:) = oBScl*Model%B0(xcc)
                    !Get 2D->3D mapping info for debugging
                    call locate(xcc,ijk,Model,ebGr)
                    ijkXY(i,j,:) = 1.0*ijk + Map2ezp(xcc,ijk,Model,ebGr)
                else
                    B02D(i,j,:) = 0.0
                    ijkXY(i,j,:)= 0.0
                endif
            enddo
        enddo

        !Write grid
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",xxi)
        call AddOutVar(IOVars,"Y",yyi)
        call AddOutVar(IOVars,"Bx0",B02D(:,:,XDIR))
        call AddOutVar(IOVars,"By0",B02D(:,:,YDIR))
        call AddOutVar(IOVars,"Bz0",B02D(:,:,ZDIR))
        call AddOutVar(IOVars,"iG",ijkXY(:,:,IDIR))
        call AddOutVar(IOVars,"jG",ijkXY(:,:,JDIR))
        call AddOutVar(IOVars,"kG",ijkXY(:,:,KDIR))            

        call WriteVars(IOVars,.true.,ebOutF)
        call ClearIO(IOVars)

        end associate
    end subroutine initEBio

    !Do output slice for eb data
    !Scale output data using values set in chimpunits
    subroutine writeEB(Model,ebState,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars
        real(rp), dimension(:,:,:,:), allocatable :: Q
        real(rp), dimension(:,:,:), allocatable :: dB2D,E2D,J2D,Vr,Src
        real(rp), dimension(:,:), allocatable :: Lb,LbXY,dLpp,rCurv

        integer :: i,j,s
        real(rp), dimension(NDIM) :: xp,xm,dB,Ep,Em,Bp,Bm,Eeq,Beq,B
        real(rp) :: MagB,MagJ,oVGScl
        real(rp), dimension(NVARMHD,0:Model%nSpc) :: Qij
        type(gcFields_T) :: gcFieldsP,gcFieldsM,gcFieldsEq
        real(rp), dimension(NDIM,NDIM) :: jB
        real(rp), dimension(NVARVOLTSRC) :: SrcijkP,SrcijkM
        character(len=strLen) :: dID,pID,vxID,vyID,vzID,vrID

        real(rp) :: oJScl
        !Data for tracing
        type(ebTrc_T), dimension(:,:), allocatable :: ebTrcIJ

        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        allocate( dB2D(Nx1,Nx2,NDIM))
        allocate(  E2D(Nx1,Nx2,NDIM))
        allocate(  J2D(Nx1,Nx2,NDIM))
        allocate(Vr  (Nx1,Nx2,0:Model%nSpc))
        allocate(Lb  (Nx1,Nx2))
        allocate(LbXY(Nx1,Nx2))
        allocate(rCurv(Nx1,Nx2))
        if(Model%doReadSrc) allocate(Src(Nx1,Nx2,NVARVOLTSRC))


        if (Model%doTrc) then
            !Create an ebTrc for each point on slice
            allocate(ebTrcIJ(Nx1,Nx2))
        endif

        if (Model%doPP) then
            allocate(dLpp(Nx1,Nx2))
            dLpp = 10.0 ! not 0, which corresponds to PP location
        endif

        dB2D  = 0.0
        E2D   = 0.0
        J2D   = 0.0
        Vr   = 0.0
        Lb   = 0.0
        LbXY = 0.0

        if (Model%doMHD) then
            allocate(Q(Nx1,Nx2,NVARMHD,0:Model%nSpc))
            Q = 0.0
        endif
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,xp,xm,Bp,Bm,Ep,Em,dB,Qij,gcFieldsP,gcFieldsM,jB,MagB,MagJ,B)
        do j=1,Nx2
            do i=1,Nx1
                !Straddle slice plane
                !Get fields at x,y,z0
                if (doXY) then
                    xp = [xxc(i,j),yyc(i,j),z0+dSlc]
                    xm = [xxc(i,j),yyc(i,j),z0-dSlc]
                else
                    xp = [xxc(i,j),y0+dSlc,yyc(i,j)]
                    xm = [xxc(i,j),y0-dSlc,yyc(i,j)]
                endif

                if(Model%doReadSrc) then
                    call ebFields(xp,Model%t,Model,ebState,Ep,Bp,gcFields=gcFieldsP,src=SrcijkP)
                    call ebFields(xm,Model%t,Model,ebState,Em,Bm,gcFields=gcFieldsM,src=SrcijkM)
                    Src(i,j,:) = 0.5*(SrcijkP + SrcijkM)
                else
                    call ebFields(xp,Model%t,Model,ebState,Ep,Bp,gcFields=gcFieldsP)
                    call ebFields(xm,Model%t,Model,ebState,Em,Bm,gcFields=gcFieldsM)
                endif

                jB = 0.5*(gcFieldsP%JacB + gcFieldsM%JacB)

                !Background already scaled to output units
                db = oBScl*0.5*(Bp+Bm)-B02D(i,j,:)
                dB2D(i,j,:) = db
                E2D (i,j,:) = oEScl*0.5*(Em+Ep)
                B = db+B02D(i,j,:)
                MagB = norm2(B)
                MagJ = oBScl*sqrt(sum(jB**2.0))

                ! radius of curvature
                !rCurv(i,j) = getRCurv(B/oBScl,jB)

                !Get MHD vars if requested
                if (Model%doMHD) then
                    !Qij = mhdInterp(xp,Model%t,Model,ebState)
                    !Q(i,j,:) = Qij
                    Qij = 0.5*(mhdInterpMF(xp,Model%t,Model,ebState) + mhdInterpMF(xm,Model%t,Model,ebState))
                    Q(i,j,:,:) = Qij
                    do s=0,Model%nSpc
                        Vr(i,j,s) = (xxc(i,j)*Qij(VELX,s) + yyc(i,j)*Qij(VELY,s))/norm2([xxc(i,j),yyc(i,j)])
                    enddo
                endif

                !Current
                J2D(i,j,:) = Jac2Curl(jB)

                !Get B lengthscale (3D)
                Lb(i,j) = MagB/max(MagJ,TINY)

                !Get B lengthscale (2D)
                !Note, want derivatives of all 3 B components wrt X,Y (not Z)
                !MagB is unchanged (all 3 components)
                !Change MagJ
                MagJ = oBScl*sqrt(sum(jB(XDIR:ZDIR,XDIR:YDIR)**2.0))
                LbXY(i,j) = MagB/max(MagJ,TINY)

                if (Model%doTrc) then
                    !Get field line topology stuff
                    call SliceFL(Model,ebState,0.5*(xp+xm),Model%t,ebTrcIJ(i,j))

                    call ebFields(ebTrcIJ(i,j)%MagEQ,Model%t,Model,ebState,Eeq,Beq,gcFields=gcFieldsEq)
                    ! radius of curvature
                    rCurv(i,j) = getRCurv(Beq,gcFieldsEq%JacB)

                endif

                if (Model%doPP) then
                    !Get distance to plasmapause for each cell
                    call deltaLpp(ebState,0.5*(xp+xm),Model%t,dLpp(i,j))
                endif
            enddo
        enddo

        !Calculate output J (current) scaling
        !Mu0 [Tm/A] 
        ! Curl(B) * oBScl/L0 => nT/cm x (1.0e-9)/(1.0e-2) => T/m
        !oJScl => A/m^2
        oJScl = ( (1.0e-7)*oBScl/L0 )/Mu0

        call ClearIO(IOVars)

        !-------------------
        !Variables to always output
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"MJD",ioTabMJD(ebState%ebTab,Model%t))

        call AddOutVar(IOVars,"dBx" ,db2D(:,:,XDIR))
        call AddOutVar(IOVars,"dBy" ,db2D(:,:,YDIR))
        call AddOutVar(IOVars,"dBz" ,db2D(:,:,ZDIR))
        call AddOutVar(IOVars,"Lb"  ,Lb  (:,:)     )
        call AddOutVar(IOVars,"LbXY",LbXY(:,:)     )
        call AddOutVar(IOVars,"RCurvature",rCurv(:,:)    )

        if (Model%doPP) then
            call AddOutVar(IOVars,"dLppeq",dLpp(:,:))
        end if

        if (Model%doTrc) then
            !Field line tracing metrics
            call AddOutVar(IOVars,"OCb" ,ebTrcIJ(:,:)%OCb )
            call AddOutVar(IOVars,"dvB" ,ebTrcIJ(:,:)%dvB/oBScl,uStr="Rx/nT")
            call AddOutVar(IOVars,"bD"  ,ebTrcIJ(:,:)%bD,uStr="#/cc")
            call AddOutVar(IOVars,"bP"  ,ebTrcIJ(:,:)%bP,uStr="nPa")
            call AddOutVar(IOVars,"bS"  ,ebTrcIJ(:,:)%bS,uStr="Wolf^(1/gamma)")
            call AddOutVar(IOVars,"bMin",ebTrcIJ(:,:)%bMin*oBScl,uStr="nT")
            call AddOutVar(IOVars,"stdD",ebTrcIJ(:,:)%stdD)
            call AddOutVar(IOVars,"stdP",ebTrcIJ(:,:)%stdP)

            !Equator and end-points
            call AddOutVar(IOVars,"xBEQ",ebTrcIJ(:,:)%MagEQ(XDIR))
            call AddOutVar(IOVars,"yBEQ",ebTrcIJ(:,:)%MagEQ(YDIR))
            call AddOutVar(IOVars,"zBEQ",ebTrcIJ(:,:)%MagEQ(ZDIR))

            call AddOutVar(IOVars,"xP",ebTrcIJ(:,:)%xEPp(XDIR))
            call AddOutVar(IOVars,"yP",ebTrcIJ(:,:)%xEPp(YDIR))
            call AddOutVar(IOVars,"zP",ebTrcIJ(:,:)%xEPp(ZDIR))

            call AddOutVar(IOVars,"xM",ebTrcIJ(:,:)%xEPm(XDIR))
            call AddOutVar(IOVars,"yM",ebTrcIJ(:,:)%xEPm(YDIR))
            call AddOutVar(IOVars,"zM",ebTrcIJ(:,:)%xEPm(ZDIR))

        endif

        if (.not. Model%doSlim) then
            call AddOutVar(IOVars,"Ex" , E2D(:,:,XDIR),uStr="mV/m")
            call AddOutVar(IOVars,"Ey" , E2D(:,:,YDIR),uStr="mV/m")
            call AddOutVar(IOVars,"Ez" , E2D(:,:,ZDIR),uStr="mV/m")
            call AddOutVar(IOVars,"Jx" , oJScl*J2D(:,:,XDIR),uStr="A/m2")
            call AddOutVar(IOVars,"Jy" , oJScl*J2D(:,:,YDIR),uStr="A/m2")
            call AddOutVar(IOVars,"Jz" , oJScl*J2D(:,:,ZDIR),uStr="A/m2")
        endif        

        if (Model%doMHD) then
            call AddOutVar(IOVars,"Vx" , oVScl*Q(:,:,VELX,BLK)    ,uStr="km/s")
            call AddOutVar(IOVars,"Vy" , oVScl*Q(:,:,VELY,BLK)    ,uStr="km/s")
            call AddOutVar(IOVars,"Vz" , oVScl*Q(:,:,VELZ,BLK)    ,uStr="km/s")
            call AddOutVar(IOVars,"Vr" , oVScl*Vr(:,:,BLK)        ,uStr="km/s")
            call AddOutVar(IOVars,"D"  ,       Q(:,:,DEN,BLK)     ,uStr="#/cc")
            call AddOutVar(IOVars,"P"  ,       Q(:,:,PRESSURE,BLK),uStr="nPa" )
        endif

        if (Model%nSpc > 0) then
            do s=1,Model%nSpc
                write(dID  ,'(A,I0)') "D"  , s
                write(pID  ,'(A,I0)') "P"  , s
                write(vxID ,'(A,I0)') "Vx" , s
                write(vyID ,'(A,I0)') "Vy" , s
                write(vzID ,'(A,I0)') "Vz" , s
                write(vrID ,'(A,I0)') "Vr" , s
                call AddOutVar(IOVars,vxID , oVScl*Q(:,:,VELX,s)    ,uStr="km/s")
                call AddOutVar(IOVars,vyID , oVScl*Q(:,:,VELY,s)    ,uStr="km/s")
                call AddOutVar(IOVars,vzID , oVScl*Q(:,:,VELZ,s)    ,uStr="km/s")
                call AddOutVar(IOVars,vrID , oVScl*Vr(:,:,s)        ,uStr="km/s")
                call AddOutVar(IOVars,dID  ,       Q(:,:,DEN,s)     ,uStr="#/cc")
                call AddOutVar(IOVars,pID  ,       Q(:,:,PRESSURE,s),uStr="nPa" )
            enddo
        endif

        if (Model%doReadSrc) then
            call AddOutVar(IOVars, "SrcX1"    , Src(:,:,PROJLAT  ))
            call AddOutVar(IOVars, "SrcX2"    , Src(:,:,PROJLON  ))
            call AddOutVar(IOVars, "SrcIONEx" , Src(:,:,IONEX    ))
            call AddOutVar(IOVars, "SrcIONEy" , Src(:,:,IONEY    ))
            call AddOutVar(IOVars, "SrcIONEz" , Src(:,:,IONEZ    ))
            call AddOutVar(IOVars, "SrcD_RING", Src(:,:,IM_D_RING))
            call AddOutVar(IOVars, "SrcP_RING", Src(:,:,IM_P_RING))
            call AddOutVar(IOVars, "SrcD_COLD", Src(:,:,IM_D_COLD))
            call AddOutVar(IOVars, "SrcP_COLD", Src(:,:,IM_P_COLD))
            call AddOutVar(IOVars, "SrcDT"    , Src(:,:,IM_TSCL  ))
        endif

        call WriteVars(IOVars,.true.,ebOutF,gStr)
        call ClearIO(IOVars)
        deallocate(dB2D,E2D,rCurv)

        end associate
    end subroutine writeEB

    !Double grid from corners
    subroutine Embiggen(xxi,yyi,Nx1,Nx2)
        real(rp), allocatable, intent(inout) :: xxi(:,:),yyi(:,:)
        integer, intent(inout) :: Nx1,Nx2

        real(rp), dimension(:,:), allocatable :: xxO,yyO
        integer :: i,j,ip,jp,nNx1,nNx2

        write(*,*) shape(xxi)
        !Start by saving old corners
        allocate(xxO(Nx1+1,Nx2+1))
        allocate(yyO(Nx1+1,Nx2+1))

        xxO = xxi
        yyO = yyi

        !Now reallocate xxi/yyi to bigger size
        nNx1 = 2*Nx1
        nNx2 = 2*Nx2

        deallocate(xxi)
        deallocate(yyi)
        allocate(xxi(nNx1+1,nNx2+1))
        allocate(yyi(nNx1+1,nNx2+1))
        !Embed old points into new grid
        !$OMP PARALLEL DO default(shared) private(i,j,ip,jp)
        do j=1,Nx2+1
            do i=1,Nx1+1
                ip = 2*i-1
                jp = 2*j-1
                xxi(ip,jp) = xxO(i,j)
                yyi(ip,jp) = yyO(i,j)
            enddo
        enddo

        !Create I-midpoints
        !$OMP PARALLEL DO default(shared) private(i,j,ip,jp)
        do i=1,Nx1
            xxi(2*i,:) = 0.5*( xxi(2*i-1,:) + xxi(2*i+1,:) )
            yyi(2*i,:) = 0.5*( yyi(2*i-1,:) + yyi(2*i+1,:) )
        enddo
        !Create J-midpoints
        !$OMP PARALLEL DO default(shared) private(i,j,ip,jp)
        do j=1,Nx2
            xxi(:,2*j) = 0.5*( xxi(:,2*j-1) + xxi(:,2*j+1) )
            yyi(:,2*j) = 0.5*( yyi(:,2*j-1) + yyi(:,2*j+1) )
        enddo

        !$OMP PARALLEL DO default(shared) private(i,j,ip,jp)
        do j=1,Nx2
            do i=1,Nx1
                ip = 2*i
                jp = 2*j
                !Create new I-J midpoints
                xxi(ip,jp) = 0.25*( xxi(ip-1,jp-1) + xxi(ip-1,jp+1) + xxi(ip+1,jp-1) + xxi(ip+1,jp+1) )
                yyi(ip,jp) = 0.25*( yyi(ip-1,jp-1) + yyi(ip-1,jp+1) + yyi(ip+1,jp-1) + yyi(ip+1,jp+1) )
            enddo
        enddo

        !Finally, change sizing information
        Nx1 = nNx1
        Nx2 = nNx2
    end subroutine Embiggen
end module sliceio
