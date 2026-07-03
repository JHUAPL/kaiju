module chopio

    use chmpdefs
    use chmpunits
    use ebtypes
    use parintime
    use ioH5
    use xml_input
    use ebinterp
    use files
    use streamline
    
    implicit none

    character(len=strLen) :: eb3DOutF
    integer, parameter :: MAXEBVS = 50
    
    integer  :: Nx1 = 64, Nx2 = 64, Nx3 = 64
    real(rp) :: originOff = 0.0 !Default origin offset
    real(rp) :: lookDir = 0.0   !Default CONE grid tilt

    real(rp) :: xyzMax = 10.0 !Default bound
    real(rp), dimension(:,:,:), allocatable :: xxi,yyi,zzi,xxc,yyc,zzc
    integer, private :: NumB = 0

    contains

    subroutine initEB3Dio(Model,ebState,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(XML_Input_T), intent(in) :: inpXML

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars
        character(len=strLen) :: idStr
        real(rp) :: x1Max,x1Min,x2Max,x2Min,x3Max,x3Min
        real(rp) :: x1,x2,x3
        real(rp) :: dx1,dx2,dx3
        integer :: i,j,k
        integer :: iS,iE
        real(rp) :: xcc(NDIM)
        logical :: doLogR
        real(rp) :: xOff,yOff,zOff  !Cartesian offsets for CONE grid origin
        real(rp) :: theta_pitch,phi_yaw !Theta and phi for y and z rotations in 3D space
        real(rp) :: xxi_temp,yyi_temp,zzi_temp  !Temp variables for 3D rotations
        real(rp) :: xxi_tilt,yyi_tilt,zzi_tilt,xxi_pan,yyi_pan,zzi_pan  !More temp variables for 3D rotations

    !Check for time parallelism
        call InitParInTime(Model,inpXML,"eb3",eb3DOutF)
    !Setup output file
        associate( ebGr=>ebState%ebGr )
        
        call CheckAndKill(eb3DOutF)

        !Figure out information for chop region
        call inpXML%Set_Val(idStr,'chop/grType',"XYZ")
        !Get 3D bounds, use xyzMax as default max
        call inpXML%Set_Val(x1Max,'chop/x1Max',xyzMax)
        call inpXML%Set_Val(x2Max,'chop/x2Max',xyzMax)
        call inpXML%Set_Val(x3Max,'chop/x3Max',xyzMax)

        !Use negative of maxes as default mins
        call inpXML%Set_Val(x1Min,'chop/x1Min',-x1Max)
        call inpXML%Set_Val(x2Min,'chop/x2Min',-x2Max)
        call inpXML%Set_Val(x3Min,'chop/x3Min',-x3Max)

        call inpXML%Set_Val(Nx1,'chop/Nx1',Nx1)
        call inpXML%Set_Val(Nx2,'chop/Nx2',Nx2)
        call inpXML%Set_Val(Nx3,'chop/Nx3',Nx3)

        !Get origin offset for CONE chopping. Default to originOff (zero)
        call inpXML%Set_Val(xOff,'chop/xOff',originOff)
        call inpXML%Set_Val(yOff,'chop/yOff',originOff)
        call inpXML%Set_Val(zOff,'chop/zOff',originOff)

        !Get pitch and yaw angles for CONE chopping. Default to lookDir (zero)
        call inpXML%Set_Val(theta_pitch,'chop/lookPol',lookDir)
        call inpXML%Set_Val(phi_yaw,'chop/lookAzi',lookDir)

        !for LFM grid, take all i-shells within x1Max, in the Sun direction
        select case(trim(toUpper(idStr)))
        !---------
        case("LFM")
            !Find first cell that's "in" according to set locator
            i = ebGr%is
            xcc = ebGr%xyz(i,ebGr%js,ebGr%ks,:)
            do while ( .not. inDomain(xcc,Model,ebGr) )
                i = i+1
                xcc = ebGr%xyz(i,ebGr%js,ebGr%ks,:)
            enddo
            iS = i
            iE = minloc(abs(ebGr%xyz(:,ebGr%js,ebGr%ks,XDIR)-x1Max),dim=1)
            write(*,*) 'i-Start/End/Total = ',iS,iE,ebGr%Nip
            
            Nx1 = iE-iS+1
            Nx2 = ebState%ebGr%Njp
            Nx3 = ebState%ebGr%Nkp
        
        end select

        !Allocate chop grids
        allocate(xxi(Nx1+1,Nx2+1,Nx3+1))
        allocate(yyi(Nx1+1,Nx2+1,Nx3+1))
        allocate(zzi(Nx1+1,Nx2+1,Nx3+1))
        
        allocate(xxc(Nx1  ,Nx2  ,Nx3  ))
        allocate(yyc(Nx1  ,Nx2  ,Nx3  ))
        allocate(zzc(Nx1  ,Nx2  ,Nx3  ))

        !Create corners
        select case(trim(toUpper(idStr)))
        !---------
        case("XYZ")
            dx1 = (x1Max-x1Min)/Nx1
            dx2 = (x2Max-x2Min)/Nx2
            dx3 = (x3Max-x3Min)/Nx3
            !$OMP PARALLEL DO
            do k=1,Nx3+1
                do j=1,Nx2+1
                    do i=1,Nx1+1
                        xxi(i,j,k) = x1Min + (i-1)*dx1
                        yyi(i,j,k) = x2Min + (j-1)*dx2
                        zzi(i,j,k) = x3Min + (k-1)*dx3
                    enddo
                enddo
            enddo
        !---------
        ! Spherical
        case("RTP")
            !Check for log spacing in r
            call inpXML%Set_Val(doLogR,'chop/doLogR',.false.)
            if (doLogR) then
                dx1 = ( log10(x1Max)-log10(x1Min) )/Nx1
            else
                dx1 = (x1Max - x1Min)/Nx1
            endif
            dx2 = (x2Max-x2Min)*deg2rad/Nx2
            dx3 = (x3Max-x3Min)*deg2rad/Nx3
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,k,x1,x2,x3)
            do k=1, Nx3+1
                do j=1,Nx2+1    
                    do i=1,Nx1+1
                        !x1 = Rin(i,j,k) + (i-1)*dx1
                        if (doLogR) then
                            x1 = 10**( log10(x1Min) + (i-1)*dx1 )
                        else
                            x1 = x1Min + (i-1)*dx1
                        endif
                        x2 = x2Min*deg2rad + (j-1)*dx2 !Theta
                        x3 = x3Min*deg2rad + (k-1)*dx3 !Phi

                        xxi(i,j,k) = x1*cos(x3)*sin(x2)
                        yyi(i,j,k) = x1*sin(x3)*sin(x2)
                        zzi(i,j,k) = x1*cos(x2)

                    enddo
                enddo
            enddo
        !---------
        ! Conical/Solid angle (based on RTP)
        case("CONE")
            !Check for log spacing in r
            call inpXML%Set_Val(doLogR,'chop/doLogR',.false.)
            if (doLogR) then
                dx1 = ( log10(x1Max)-log10(x1Min) )/Nx1
            else
                dx1 = (x1Max - x1Min)/Nx1
            endif
            dx2 = (x2Max-x2Min)*deg2rad/Nx2
            dx3 = (x3Max-x3Min)*deg2rad/Nx3
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(i,j,k,x1,x2,x3,xxi_temp,yyi_temp,zzi_temp,xxi_tilt,yyi_tilt,zzi_tilt,xxi_pan,yyi_pan,zzi_pan)
            do k=1,Nx3+1
                do j=1,Nx2+1    
                    do i=1,Nx1+1
                        !x1 = Rin(i,j,k) + (i-1)*dx1
                        if (doLogR) then
                            x1 = 10**( log10(x1Min) + (i-1)*dx1 )
                        else
                            x1 = x1Min + (i-1)*dx1
                        endif
                        x2 = x2Min*deg2rad + (j-1)*dx2 !Theta
                        x3 = x3Min*deg2rad + (k-1)*dx3 !Phi
                        
                        ! pre-rotation variables
                        xxi_temp = x1*cos(x3)*sin(x2)
                        yyi_temp = x1*sin(x3)*sin(x2)
                        zzi_temp = x1*cos(x2)

                        ! first tilt (polar) then pan (azimuthal) (pitch then yaw)
                        ! perform rotation about y axis (theta, polar angle)
                        xxi_tilt = (cos(theta_pitch*deg2rad)*xxi_temp) + (0) + (sin(theta_pitch*deg2rad)*zzi_temp)
                        yyi_tilt = (0) + (yyi_temp) + (0)
                        zzi_tilt = (-sin(theta_pitch*deg2rad)*xxi_temp) + (0) + (cos(theta_pitch*deg2rad)*zzi_temp)

                        ! perform rotation about z axis (phi, azimuthal angle)
                        ! apply to the tilted xxi/yyi/zzi
                        xxi_pan = (cos(phi_yaw*deg2rad)*xxi_tilt) + (-sin(phi_yaw*deg2rad)*yyi_tilt) + (0)
                        yyi_pan = (sin(phi_yaw*deg2rad)*xxi_tilt) + (cos(phi_yaw*deg2rad)*yyi_tilt) + (0)
                        zzi_pan = (0) + (0) + (zzi_tilt)
                        
                        ! post-rotation, add cartesian offsets
                        xxi(i,j,k) = xxi_pan + xOff
                        yyi(i,j,k) = yyi_pan + yOff
                        zzi(i,j,k) = zzi_pan + zOff

                    enddo
                enddo
            enddo
        !---------
        ! only x1Max is used to set how far downtail you chop
        case("LFM")
            !$OMP PARALLEL DO
            do k=1,Nx3+1
                do j=1,Nx2+1
                    do i=1,Nx1+1
                        xxi(i,j,k) = ebGr%xyz(i+iS-1,j,k,XDIR)
                        yyi(i,j,k) = ebGr%xyz(i+iS-1,j,k,YDIR)
                        zzi(i,j,k) = ebGr%xyz(i+iS-1,j,k,ZDIR)
                    enddo
                enddo
            enddo

        case default
            write(*,*) 'Unknown chop type, exiting ...'
            stop

        end select
        !Calculate cell centers (just simple 8 point average)
        !$OMP PARALLEL DO
        do k=1,Nx3
            do j=1,Nx2
                do i=1,Nx1

                    xxc(i,j,k) = 0.125*( xxi(i,j,k  ) + xxi(i+1,j,k  ) + xxi(i,j+1,k  ) + xxi(i+1,j+1,k  ) + &
                                       & xxi(i,j,k+1) + xxi(i+1,j,k+1) + xxi(i,j+1,k+1) + xxi(i+1,j+1,k+1)    )
                    yyc(i,j,k) = 0.125*( yyi(i,j,k  ) + yyi(i+1,j,k  ) + yyi(i,j+1,k  ) + yyi(i+1,j+1,k  ) + &
                                       & yyi(i,j,k+1) + yyi(i+1,j,k+1) + yyi(i,j+1,k+1) + yyi(i+1,j+1,k+1)    )
                    zzc(i,j,k) = 0.125*( zzi(i,j,k  ) + zzi(i+1,j,k  ) + zzi(i,j+1,k  ) + zzi(i+1,j+1,k  ) + &
                                       & zzi(i,j,k+1) + zzi(i+1,j,k+1) + zzi(i,j+1,k+1) + zzi(i+1,j+1,k+1)    )


                enddo
            enddo
        enddo

        !Write grid
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",xxi)
        call AddOutVar(IOVars,"Y",yyi)
        call AddOutVar(IOVars,"Z",zzi)

        call WriteVars(IOVars,.true.,eb3DOutF)

        call ClearIO(IOVars)  
        end associate     
    end subroutine initEB3Dio

    !Do output slice for eb data
    !Scale output data using values set in chimpunits
    subroutine writeEB3D(Model,ebState,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars

        real(rp), dimension(:,:,:,:)  , allocatable :: B,E,J3
        real(rp), dimension(:,:,:,:,:), allocatable :: Q
        real(rp), dimension(NVARMHD,0:Model%nSpc) :: Qijk
        real(rp), dimension(NDIM) :: Bijk,Eijk,xyz
        type(gcFields_T) :: gcFields
        real(rp), dimension(NDIM,NDIM) :: jB
        integer :: i,j,k,s
        character(len=strLen) :: dID,pID,vxID,vyID,vzID

        real(rp) :: oJScl

        !Data for tracing
        type(ebTrc_T), dimension(:,:,:), allocatable :: ebTrcIJK

        allocate(B (Nx1,Nx2,Nx3,NDIM))
        allocate(E (Nx1,Nx2,Nx3,NDIM))
        allocate(J3(Nx1,Nx2,Nx3,NDIM))
        allocate(Q (Nx1,Nx2,Nx3,NVARMHD,0:Model%nSpc))

        if (Model%doTrc) then
            !Create an ebTrc for each point on slice
            allocate(ebTrcIJK(Nx1,Nx2,Nx3))
        endif

        B = 0.0
        E = 0.0
        J = 0.0
        Q = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyz,Bijk,Eijk,Qijk,gcFields,jB)
        do k=1,Nx3
            do j=1,Nx2
                do i=1,Nx1
                    xyz = [xxc(i,j,k),yyc(i,j,k),zzc(i,j,k)]

                    call ebFields(xyz,Model%t,Model,ebState,Eijk,Bijk,gcFields=gcFields)
                    B(i,j,k,:) = Bijk
                    E(i,j,k,:) = Eijk

                    jB = gcFields%JacB

                    J3(i,j,k,:) = Jac2Curl(jB)

                    if (Model%doMHD) then
                        Qijk = mhdInterpMF(xyz,Model%t,Model,ebState)
                        Q(i,j,k,:,:) = Qijk
                    endif

                    if (Model%doTrc) then
                        !Get field line topology stuff
                        call SliceFL(Model,ebState,xyz,Model%t,ebTrcIJK(i,j,k))
                    endif

                enddo
            enddo
        enddo

        !Calculate output J (current) scaling
        !Mu0 [Tm/A] 
        ! Curl(B) * oBScl/L0 => nT/cm x (1.0e-9)/(1.0e-2) => T/m
        !oJScl => A/m^2
        oJScl = ( (1.0e-7)*oBScl/L0 )/Mu0
        
        !Now write out
        call ClearIO(IOVars)

        !-------------------
        !Variables to always output
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"MJD",ioTabMJD(ebState%ebTab,Model%t))
        call AddOutVar(IOVars,"Bx"  ,oBScl*B(:,:,:,XDIR),uStr="nT")
        call AddOutVar(IOVars,"By"  ,oBScl*B(:,:,:,YDIR),uStr="nT")
        call AddOutVar(IOVars,"Bz"  ,oBScl*B(:,:,:,ZDIR),uStr="nT")

        call AddOutVar(IOVars,"Ex"  ,oEScl*E(:,:,:,XDIR),uStr="mV/m")
        call AddOutVar(IOVars,"Ey"  ,oEScl*E(:,:,:,YDIR),uStr="mV/m")
        call AddOutVar(IOVars,"Ez"  ,oEScl*E(:,:,:,ZDIR),uStr="mV/m")

        call AddOutVar(IOVars,"Jx"  ,oJScl*J3(:,:,:,XDIR),uStr="A/m2")
        call AddOutVar(IOVars,"Jy"  ,oJScl*J3(:,:,:,YDIR),uStr="A/m2")
        call AddOutVar(IOVars,"Jz"  ,oJScl*J3(:,:,:,ZDIR),uStr="A/m2")

        if (Model%doMHD) then
            call AddOutVar(IOVars,"Vx"  ,oVScl*Q(:,:,:,VELX    ,BLK),uStr="km/s")
            call AddOutVar(IOVars,"Vy"  ,oVScl*Q(:,:,:,VELY    ,BLK),uStr="km/s")
            call AddOutVar(IOVars,"Vz"  ,oVScl*Q(:,:,:,VELZ    ,BLK),uStr="km/s")
            call AddOutVar(IOVars,"D"   ,      Q(:,:,:,DEN     ,BLK),uStr="#/cc")
            call AddOutVar(IOVars,"P"   ,      Q(:,:,:,PRESSURE,BLK),uStr="nPa" )
            if (Model%nSpc > 0) then
                do s=1,Model%nSpc
                    write(dID  ,'(A,I0)') "D"  , s
                    write(pID  ,'(A,I0)') "P"  , s
                    write(vxID ,'(A,I0)') "Vx" , s
                    write(vyID ,'(A,I0)') "Vy" , s
                    write(vzID ,'(A,I0)') "Vz" , s
                    call AddOutVar(IOVars,vxID , oVScl*Q(:,:,:,VELX    ,s),uStr="km/s")
                    call AddOutVar(IOVars,vyID , oVScl*Q(:,:,:,VELY    ,s),uStr="km/s")
                    call AddOutVar(IOVars,vzID , oVScl*Q(:,:,:,VELZ    ,s),uStr="km/s")
                    call AddOutVar(IOVars,dID  ,       Q(:,:,:,DEN     ,s),uStr="#/cc")
                    call AddOutVar(IOVars,pID  ,       Q(:,:,:,PRESSURE,s),uStr="nPa" )
                enddo
            endif
        endif

        if (Model%doTrc) then
            !Field line tracing metrics
            call AddOutVar(IOVars,"OCb" ,ebTrcIJK(:,:,:)%OCb )
            call AddOutVar(IOVars,"dvB" ,ebTrcIJK(:,:,:)%dvB/oBScl,uStr="Rx/nT")
            call AddOutVar(IOVars,"bD"  ,ebTrcIJK(:,:,:)%bD,uStr="#/cc")
            call AddOutVar(IOVars,"bP"  ,ebTrcIJK(:,:,:)%bP,uStr="nPa")
            call AddOutVar(IOVars,"bS"  ,ebTrcIJK(:,:,:)%bS,uStr="Wolf^(1/gamma)")
            call AddOutVar(IOVars,"bMin",ebTrcIJK(:,:,:)%bMin*oBScl,uStr="nT")

            !Equator and end-points
            call AddOutVar(IOVars,"xBEQ",ebTrcIJK(:,:,:)%MagEQ(XDIR))
            call AddOutVar(IOVars,"yBEQ",ebTrcIJK(:,:,:)%MagEQ(YDIR))
            call AddOutVar(IOVars,"zBEQ",ebTrcIJK(:,:,:)%MagEQ(ZDIR))

            call AddOutVar(IOVars,"xP",ebTrcIJK(:,:,:)%xEPp(XDIR))
            call AddOutVar(IOVars,"yP",ebTrcIJK(:,:,:)%xEPp(YDIR))
            call AddOutVar(IOVars,"zP",ebTrcIJK(:,:,:)%xEPp(ZDIR))

            call AddOutVar(IOVars,"xM",ebTrcIJK(:,:,:)%xEPm(XDIR))
            call AddOutVar(IOVars,"yM",ebTrcIJK(:,:,:)%xEPm(YDIR))
            call AddOutVar(IOVars,"zM",ebTrcIJK(:,:,:)%xEPm(ZDIR))

        endif

        call WriteVars(IOVars,.true.,eb3DOutF,gStr)
        call ClearIO(IOVars)
        
    end subroutine writeEB3D
end module chopio
