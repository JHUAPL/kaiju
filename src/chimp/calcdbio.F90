module calcdbio

	use chmpdefs
	use chmpunits
    use parintime
	use ebtypes
	use ioH5
	use xml_input
	use files
	use clocks
    use calcdbutils
    use calcdbcore
    use ebinterp
    
	implicit none

    character(len=strLen), private :: dbOutF
    integer, parameter, private :: MAXDBVS = 80
    integer, parameter, private :: RDIR=1,TDIR=2,PDIR=3
    real(rp), private :: dzGG = 10.0 !Default height spacing [km]
    real(rp), private :: z0 = 0.0 !Starting height [km]

    contains

    subroutine initDBio(Model,ebState,gGr,inpXML,NumP)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)    :: ebState
        type(grGrid_T), intent(inout) :: gGr
        type(XML_Input_T), intent(in) :: inpXML
        integer, intent(inout) :: NumP

        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        character(len=strLen) :: cID,inH5

        integer :: i,j,k,NLat,NLon,Nz,dOut
        integer, dimension(NDIM) :: Nijk
        real(rp) :: z,R,lat,phi
        real(rp), dimension(:,:,:,:), allocatable :: SphI,SphC !Spherical coordinates
        logical :: doH5g

        associate( ebGr=>ebState%ebGr )
        !Equate dtout/dt since the difference doesn't matter here
        Model%dtOut = Model%dt

    !Read info from XML
        call inpXML%Set_Val(doH5g,'Grid/doH5g',.false.)
        if (doH5g) then
            call inpXML%Set_Val(inH5,"Grid/H5Grid","grid.h5")
            call CheckFileOrDie(inH5,"Input grid file not found")
            Nijk = GridSizeH5(inH5)
            
            NLat = Nijk(IDIR)-1
            NLon = Nijk(JDIR)-1
            Nz   = Nijk(KDIR)-1
        else
            !Generate grid
            call inpXML%Set_Val(NLat,'Grid/NLat',45) !Number of latitudinal cells
            call inpXML%Set_Val(NLon,'Grid/NLon',90) !Number of longitudinal cells
            call inpXML%Set_Val(Nz  ,'Grid/Nz'  , 2) !Number of longitudinal cells
            call inpXML%Set_Val(dzGG,'Grid/dzGG',dzGG)
            call inpXML%Set_Val(z0  ,'Grid/z0'  ,z0)


        endif


        !Stuff to read for either grid type
        call inpXML%Set_Val(gGr%rMax,'CalcDB/rMax',gGr%rMax)
        call inpXML%Set_Val(gGr%doGEO,'Grid/doGEO',gGr%doGEO) !Whether to do GEO on ground
        if (gGr%doGEO) then
            cID = "GEO"
        else
            cID = "SM"
        endif
    !Check for time parallelism
        call InitParInTime(Model,inpXML,"deltab",dbOutF)

    !Setup output file
        call CheckAndKill(dbOutF)

    !Create ground grid
        NumP = NLat*NLon*Nz

        !Store sizes
        gGr%NLat = NLat
        gGr%NLon = NLon
        gGr%Nz   = Nz

        !Calculate corner/center
        allocate(gGr%GxyzI (NLat+1,NLon+1,Nz+1,NDIM))
        allocate(gGr%GxyzC (NLat  ,NLon  ,Nz  ,NDIM))
        allocate(gGr%SMxyzC(NLat  ,NLon  ,Nz  ,NDIM))

        if (doH5g) then
        !Read grid from file
            write(*,*) "Reading grid from file ..."
            call ClearIO(IOVars)

            call AddInVar(IOVars,"X")
            call AddInVar(IOVars,"Y")
            call AddInVar(IOVars,"Z")

            call AddInVar(IOVars,"Xcc")
            call AddInVar(IOVars,"Ycc")
            call AddInVar(IOVars,"Zcc")

            call ReadVars(IOVars,.false.,inH5) !Don't use io precision

            !Now reshape and store
            call IOArray3DFill(IOVars,"X"  ,gGr%GxyzI(:,:,:,XDIR))
            call IOArray3DFill(IOVars,"Y"  ,gGr%GxyzI(:,:,:,YDIR))
            call IOArray3DFill(IOVars,"Z"  ,gGr%GxyzI(:,:,:,ZDIR))

            call IOArray3DFill(IOVars,"Xcc",gGr%GxyzC(:,:,:,XDIR))
            call IOArray3DFill(IOVars,"Ycc",gGr%GxyzC(:,:,:,YDIR))
            call IOArray3DFill(IOVars,"Zcc",gGr%GxyzC(:,:,:,ZDIR))

        else
        !Generate grid
            write(*,*) "Generating grid ..."
            !Do corners
            do k=1,Nz+1
            	z = -0.5*dzGG + (k-1)*dzGG + z0 !km above ground
            	z = (1.0e+5)*z !cm above ground
            	R = 1.0 + z/Re_cgs !Assuming Earth here

            	do j=1,NLon+1
            		phi = (j-1)*360.0/NLon 
            		do i=1,NLat+1
            			lat = -90.0 + (i-1)*180.0/NLat
            			gGr%GxyzI(i,j,k,XDIR) = R*cos(lat*PI/180.0)*cos(phi*PI/180.0)
            			gGr%GxyzI(i,j,k,YDIR) = R*cos(lat*PI/180.0)*sin(phi*PI/180.0)
            			gGr%GxyzI(i,j,k,ZDIR) = R*sin(lat*PI/180.0)

            		enddo !i
            	enddo !j
            enddo

            !Do centers
            !Not using 8-pt average due to non-uniform radius
            write(*,*) "Height range [km] : ", z0,(Nz-1)*dzGG+z0
            
            do k=1,Nz
                z = (k-1)*dzGG + z0 !km above ground
                z = (1.0e+5)*z !cm above ground
                R = 1.0 + z/Re_cgs !Assuming Earth here

            	do j=1,NLon
                    phi = (j-0.5)*360.0/NLon
            		do i=1,NLat
                        lat = -90.0 + (i-0.5)*180.0/NLat
                        gGr%GxyzC(i,j,k,XDIR) = R*cos(lat*PI/180.0)*cos(phi*PI/180.0)
                        gGr%GxyzC(i,j,k,YDIR) = R*cos(lat*PI/180.0)*sin(phi*PI/180.0)
                        gGr%GxyzC(i,j,k,ZDIR) = R*sin(lat*PI/180.0)

            		enddo
            	enddo
            enddo
        endif !Grid generation, doH5g

        !Just set SM=G for now
        gGr%SMxyzC = gGr%GxyzC

        end associate

    !Get matching spherical grids
        allocate(SphI(NLat+1,NLon+1,Nz+1,NDIM))
        allocate(SphC(NLat  ,NLon  ,Nz  ,NDIM))

        SphI(:,:,:,RDIR) = norm2(gGr%GxyzI,dim=4)
        SphC(:,:,:,RDIR) = norm2(gGr%GxyzC,dim=4)

        SphI(:,:,:,TDIR) = acos(gGr%GxyzI(:,:,:,ZDIR)/SphI(:,:,:,RDIR))
        SphC(:,:,:,TDIR) = acos(gGr%GxyzC(:,:,:,ZDIR)/SphC(:,:,:,RDIR))

        SphI(:,:,:,PDIR) = katan2(gGr%GxyzI(:,:,:,YDIR),gGr%GxyzI(:,:,:,XDIR))
        SphC(:,:,:,PDIR) = katan2(gGr%GxyzC(:,:,:,YDIR),gGr%GxyzC(:,:,:,XDIR))


        !Allocate db holders (XYZ)
        allocate(gGr%dbMAG_xyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Magnetospheric delta-B
        allocate(gGr%dbION_xyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Ionospheric delta-B
        allocate(gGr%dbFAC_xyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !FAC delta-B
        !Spherical holders
        allocate(gGr%dbMAG_rtp(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Magnetospheric delta-B
        allocate(gGr%dbION_rtp(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Ionospheric delta-B
        allocate(gGr%dbFAC_rtp(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !FAC delta-B

        gGr%dbMAG_xyz = 0.0
        gGr%dbION_xyz = 0.0
        gGr%dbFAC_xyz = 0.0
        gGr%dbMAG_rtp = 0.0
        gGr%dbION_rtp = 0.0
        gGr%dbFAC_rtp = 0.0

        !Allocate possibly time-varying coordinates for geomagnetic lat/lon and northward deflection
        allocate(gGr%smlat(gGr%NLat,gGr%NLon,gGr%Nz))
        allocate(gGr%smlon(gGr%NLat,gGr%NLon,gGr%Nz))
        allocate(gGr%dBn  (gGr%NLat,gGr%NLon,gGr%Nz))

        gGr%smlat = 0.0
        gGr%smlon = 0.0
        gGr%dBn   = 0.0
        
    !Write grid
        !Adding argument to cc-variables to make them not get squeezed out
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",gGr%GxyzI(:,:,:,XDIR),uStr="Re")
        call AddOutVar(IOVars,"Y",gGr%GxyzI(:,:,:,YDIR),uStr="Re")
        call AddOutVar(IOVars,"Z",gGr%GxyzI(:,:,:,ZDIR),uStr="Re")

        call AddOutVar(IOVars,"Xcc",gGr%GxyzC(:,:,:,XDIR),uStr="Re",doSqzO=.false.)
        call AddOutVar(IOVars,"Ycc",gGr%GxyzC(:,:,:,YDIR),uStr="Re",doSqzO=.false.)
        call AddOutVar(IOVars,"Zcc",gGr%GxyzC(:,:,:,ZDIR),uStr="Re",doSqzO=.false.)

        call AddOutVar(IOVars,"Rad"  ,SphI(:,:,:,RDIR),uStr="Re")
        call AddOutVar(IOVars,"Theta",SphI(:,:,:,TDIR),uStr="Re")
        call AddOutVar(IOVars,"Phi"  ,SphI(:,:,:,PDIR),uStr="Re")

        call AddOutVar(IOVars,"Radcc"  ,SphC(:,:,:,RDIR),uStr="Re",doSqzO=.false.)
        call AddOutVar(IOVars,"Thetacc",SphC(:,:,:,TDIR),uStr="Re",doSqzO=.false.)
        call AddOutVar(IOVars,"Phicc"  ,SphC(:,:,:,PDIR),uStr="Re",doSqzO=.false.)

        call AddOutVar(IOVars,"CoordinatesID",cID)
        call AddOutVar(IOVars,"Re",Re_km,uStr="km")

        call WriteVars(IOVars,.true.,dbOutF)
        call ClearIO(IOVars)

    end subroutine initDBio

    subroutine initRM(Model,ebState,rmState,inpXML)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(rmState_T), intent(inout) :: rmState
        type(XML_Input_T), intent(in) :: inpXML

        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        integer :: Ni,Nj
        character(len=strLen) :: rmF !Remix file
        write(rmF,'(2a)') trim(adjustl(ebState%ebTab%bStr)),'.mix.h5'

        !Decide whether to add corotation potential
        call inpXML%Set_Val(rmState%doHall ,'CalcDB/doHall' ,.true. ) 
        call inpXML%Set_Val(rmState%doPed  ,'CalcDB/doPed'  ,.true. )

        write(*,*) 'Initializing w/ ', trim(rmF)

        call ClearIO(IOVars)
        call AddInVar(IOVars,"X")
        call AddInVar(IOVars,"Y")
        call ReadVars(IOVars,.true.,rmF)

        Ni = IOVars(1)%dims(1) - 1
        Nj = IOVars(2)%dims(2) - 1

        !Initialize 4 hemispheres (N/S i1/i2)
        call initHemi(rmState%rmN1,Ni,Nj)
        call initHemi(rmState%rmN2,Ni,Nj)
        call initHemi(rmState%rmS1,Ni,Nj)
        call initHemi(rmState%rmS2,Ni,Nj)

        !Now do main rmState
        rmState%Np  = Ni !Phi cells
        rmState%Nth = Nj !Theta cells

        allocate(rmState%XY(Ni+1,Nj+1,XDIR:YDIR))

        allocate(rmState%nFac (Ni,Nj))
        allocate(rmState%nPot (Ni,Nj))
        allocate(rmState%nSigP(Ni,Nj))
        allocate(rmState%nSigH(Ni,Nj))
        allocate(rmState%sFac (Ni,Nj))
        allocate(rmState%sPot (Ni,Nj))
        allocate(rmState%sSigP(Ni,Nj))
        allocate(rmState%sSigH(Ni,Nj))

        call IOArray2DFill(IOVars,"X",rmState%XY(:,:,XDIR))
        call IOArray2DFill(IOVars,"Y",rmState%XY(:,:,YDIR))


        call hemi2rm(rmState,0.0_rp,0.0_rp) !Zero out main state arrays

        contains

            subroutine initHemi(rmHemi,Ni,Nj)
                type(rmHemi_T), intent(inout) :: rmHemi
                integer, intent(in) :: Ni,Nj

                rmHemi%nStp = -1 !Not yet set
                rmHemi%time = 0.0
                rmHemi%Np  = Ni !Phi cells
                rmHemi%Nth = Nj !Theta cells

                allocate(rmHemi%xFac (Ni,Nj))
                allocate(rmHemi%xPot (Ni,Nj))
                allocate(rmHemi%xSigP(Ni,Nj))
                allocate(rmHemi%xSigH(Ni,Nj))

                rmHemi%xFac  = 0.0
                rmHemi%xPot  = 0.0
                rmHemi%xSigP = 0.0
                rmHemi%xSigH = 0.0

            end subroutine initHemi
    end subroutine initRM

    !Update remix data
    subroutine updateRemix(Model,ebState,t,rmState)
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)   :: ebState
        real(rp), intent(in) :: t
        type(rmState_T), intent(inout)   :: rmState

        character(len=strLen) :: rmF !Remix file    
        integer :: i1,i2
        real(rp) :: w1,w2

        !TODO: Remove redundant code here
        write(rmF,'(2a)') trim(adjustl(ebState%ebTab%bStr)),'.mix.h5'
        
        call GetTabSlc(ebState%ebTab,t,i1,i2)

        !Read 4 hemispheres
        call readHemi(rmState%rmN1,rmF,ebState%ebTab,i1,NORTH)
        call readHemi(rmState%rmS1,rmF,ebState%ebTab,i1,SOUTH)

        call readHemi(rmState%rmN2,rmF,ebState%ebTab,i2,NORTH)
        call readHemi(rmState%rmS2,rmF,ebState%ebTab,i2,SOUTH)

        !Now fill in remix main state for this t
        call GetTWgts(Model,ebState,t,w1,w2)
        call hemi2rm(rmState,w1,w2)

        contains

        !Read hemisphere data
        !NOTE: nStp here refers to array index in ebTab, not necessarily Step#X
        subroutine readHemi(rmHemi,rmF,ebTab,nStp,nsID)
            type(rmHemi_T), intent(inout) :: rmHemi
            character(len=strLen), intent(in) :: rmF
            type(ioTab_T), intent(in) :: ebTab
            integer, intent(in) :: nStp,nsID

            character(len=strLen) :: hID,gStr
            type(IOVAR_T), dimension(MAXDBVS) :: IOVars

            !Check to see if we need to read
            if (rmHemi%nStp == nStp) then
                !We've already read this data
                return
            endif
            !Otherwise get the data
            !Which hemisphere?
            if (nsID == NORTH) then
                hID = "NORTH"
            else
                hID = "SOUTH"
            endif
            gStr = trim(ebState%ebTab%gStrs(nStp))

            rmHemi%time = ebTab%times(nStp)
            rmHemi%nStp = nStp

            call ClearIO(IOVars)
            call AddInVar(IOVars,"Field-aligned current " // hID)
            call AddInVar(IOVars, "Pedersen conductance " // hID)
            call AddInVar(IOVars,     "Hall conductance " // hID)
            call AddInVar(IOVars,            "Potential " // hID)
            call ReadVars(IOVars,.true.,rmF,gStr)

            !Pull from arrays
            call IOArray2DFill(IOVars,"Field-aligned current " // hID,rmHemi%xFac )
            call IOArray2DFill(IOVars, "Pedersen conductance " // hID,rmHemi%xSigP)
            call IOArray2DFill(IOVars,     "Hall conductance " // hID,rmHemi%xSigH)
            call IOArray2DFill(IOVars,            "Potential " // hID,rmHemi%xPot )

        end subroutine readHemi

    end subroutine updateRemix

    !Calculate and output delta-B data
    subroutine writeDB(Model,ebState,gGr,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(grGrid_T), intent(inout) :: gGr
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        real(rp) :: mjd
        real(rp), dimension(:,:,:,:), allocatable :: dbRTP
        real(rp), dimension(:,:,:)  , allocatable :: dbJ

        logical :: doDebugOut = .true.

        write(*,*) 'Writing ', trim(gStr)

        mjd = ioTabMJD(ebState%ebTab,Model%t)
        
        !Do conversion to spherical
        call Grxyz2rtp(gGr,gGr%dbMAG_xyz,gGr%dbMAG_rtp)
        call Grxyz2rtp(gGr,gGr%dbION_xyz,gGr%dbION_rtp)
        call Grxyz2rtp(gGr,gGr%dbFAC_xyz,gGr%dbFAC_rtp)

        !Get total perturbation
        allocate(dbRTP(gGr%NLat,gGr%NLon,gGr%Nz,NDIM))
        allocate(dbJ  (gGr%NLat,gGr%NLon,gGr%Nz))

        dbRTP = gGr%dbMAG_rtp + gGr%dbION_rtp + gGr%dbFAC_rtp

        !Write output data
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"MJD",mjd)

        !Indices
        call AddOutVar(IOVars,"SML"     ,gGr%SML      ,uStr="nT")
        call AddOutVar(IOVars,"SML_MLat",gGr%SML_MLat,uStr="deg")
        call AddOutVar(IOVars,"SML_MLon",gGr%SML_MLon,uStr="deg")

        call AddOutVar(IOVars,"SMU"     ,gGr%SMU      ,uStr="nT")
        call AddOutVar(IOVars,"SMU_MLat",gGr%SMU_MLat,uStr="deg")
        call AddOutVar(IOVars,"SMU_MLon",gGr%SMU_MLon,uStr="deg")

        call AddOutVar(IOVars,"SME",gGr%SME,uStr="nT")
        call AddOutVar(IOVars,"SMO",gGr%SMO,uStr="nT")
        call AddOutVar(IOVars,"SMR",gGr%SMR,uStr="nT")

        call AddOutVar(IOVars,"SMR_00",gGr%SMR_00,uStr="nT")
        call AddOutVar(IOVars,"SMR_06",gGr%SMR_06,uStr="nT")
        call AddOutVar(IOVars,"SMR_12",gGr%SMR_12,uStr="nT")
        call AddOutVar(IOVars,"SMR_18",gGr%SMR_18,uStr="nT")

        call AddOutVar(IOVars,"SML_00",gGr%SML_00,uStr="nT")
        call AddOutVar(IOVars,"SML_06",gGr%SML_06,uStr="nT")
        call AddOutVar(IOVars,"SML_12",gGr%SML_12,uStr="nT")
        call AddOutVar(IOVars,"SML_18",gGr%SML_18,uStr="nT")

        call AddOutVar(IOVars,"SMU_00",gGr%SMU_00,uStr="nT")
        call AddOutVar(IOVars,"SMU_06",gGr%SMU_06,uStr="nT")
        call AddOutVar(IOVars,"SMU_12",gGr%SMU_12,uStr="nT")
        call AddOutVar(IOVars,"SMU_18",gGr%SMU_18,uStr="nT")

    !Write out spherical vectors (XDIR:ZDIR = RDIR,TDIR,PDIR)
        if (Model%doFat) then
            !Magnetospheric
            call AddOutVar(IOVars,"dBrM" ,gGr%dbMAG_rtp(:,:,:,RDIR),uStr="nT",doSqzO=.false.)
            call AddOutVar(IOVars,"dBtM" ,gGr%dbMAG_rtp(:,:,:,TDIR),uStr="nT",doSqzO=.false.)
            call AddOutVar(IOVars,"dBpM" ,gGr%dbMAG_rtp(:,:,:,PDIR),uStr="nT",doSqzO=.false.)
            call CalcJdb(gGr,gGr%dbMAG_rtp,dbJ,"MAG")
            call AddOutVar(IOVars,"dbJM" ,dbJ,uStr="microA/m2",doSqzO=.false.)
            !Ionospheric
            call AddOutVar(IOVars,"dBrI" ,gGr%dbION_rtp(:,:,:,RDIR),uStr="nT",doSqzO=.false.)
            call AddOutVar(IOVars,"dBtI" ,gGr%dbION_rtp(:,:,:,TDIR),uStr="nT",doSqzO=.false.)
            call AddOutVar(IOVars,"dBpI" ,gGr%dbION_rtp(:,:,:,PDIR),uStr="nT",doSqzO=.false.)
            call CalcJdb(gGr,gGr%dbION_rtp,dbJ,"ION")
            call AddOutVar(IOVars,"dbJI" ,dbJ,uStr="microA/m2",doSqzO=.false.)
            !Field-aligned
            call AddOutVar(IOVars,"dBrF" ,gGr%dbFAC_rtp(:,:,:,RDIR),uStr="nT",doSqzO=.false.)
            call AddOutVar(IOVars,"dBtF" ,gGr%dbFAC_rtp(:,:,:,TDIR),uStr="nT",doSqzO=.false.)
            call AddOutVar(IOVars,"dBpF" ,gGr%dbFAC_rtp(:,:,:,PDIR),uStr="nT",doSqzO=.false.)
            call CalcJdb(gGr,gGr%dbFAC_rtp,dbJ,"FAC")
            call AddOutVar(IOVars,"dbJF" ,dbJ,uStr="microA/m2",doSqzO=.false.)

        endif
        
        if (.not. Model%doSlim) then
            call AddOutVar(IOVars,"smlat" ,gGr%smlat,uStr="deg",doSqzO=.false.)
            call AddOutVar(IOVars,"smlon" ,gGr%smlon,uStr="deg",doSqzO=.false.)
            call AddOutVar(IOVars,"dBn"   ,gGr%dBn  ,uStr="nT" ,doSqzO=.false.)
        endif

        call AddOutVar(IOVars,"dBr" ,dbRTP(:,:,:,RDIR),uStr="nT",doSqzO=.false.)
        call AddOutVar(IOVars,"dBt" ,dbRTP(:,:,:,TDIR),uStr="nT",doSqzO=.false.)
        call AddOutVar(IOVars,"dBp" ,dbRTP(:,:,:,PDIR),uStr="nT",doSqzO=.false.)
        
        if (Model%doFat) then
            call CalcJdb(gGr,dbRTP,dbJ,"TOT")
            call AddOutVar(IOVars,"dbJ" ,dbJ,uStr="microA/m2",doSqzO=.false.)
        endif

        call WriteVars(IOVars,.true.,dbOutF,gStr)
        call ClearIO(IOVars)
        
    end subroutine writeDB

    !Convert XYZ to RTP vectors
    subroutine Grxyz2rtp(gGr,dbXYZ,dbRTP)
        type(grGrid_T), intent(in) :: gGr
        real(rp), dimension(:,:,:,:), intent(in ) :: dbXYZ
        real(rp), dimension(:,:,:,:), intent(out) :: dbRTP

        integer :: i,j,k
        real(rp) :: rad,theta,phi,dBx,dBy,dBz
        dbRTP = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,rad,theta,phi,dBx,dBy,dBz)
        do k=1,gGr%Nz
            do j=1,gGr%NLon
                do i=1,gGr%NLat
                    rad  = norm2(gGr%GxyzC(i,j,k,:))
                    theta = acos(gGr%GxyzC(i,j,k,ZDIR)/rad)
                    phi  = atan2(gGr%GxyzC(i,j,k,YDIR),gGr%GxyzC(i,j,k,XDIR))
                    dBx = dbXYZ(i,j,k,XDIR)
                    dBy = dbXYZ(i,j,k,YDIR)
                    dBz = dbXYZ(i,j,k,ZDIR)

                    !Radial,Theta,Phi system
                    dbRTP(i,j,k,RDIR) =  dBx*sin(theta)*cos(phi) + dBy*sin(theta)*sin(phi) + dBz*cos(theta)
                    dbRTP(i,j,k,TDIR) =  dBx*cos(theta)*cos(phi) + dBy*cos(theta)*sin(phi) - dBz*sin(theta)
                    dbRTP(i,j,k,PDIR) = -dBx           *sin(phi) + dBy           *cos(phi) 

                enddo
            enddo
        enddo

    end subroutine Grxyz2rtp
    
    subroutine CalcJdb(gGr,dbRTP,dbJ,jID)
        type(grGrid_T), intent(in) :: gGr
        real(rp), intent(in)  :: dbRTP(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)
        real(rp), intent(out) :: dbJ  (gGr%NLat,gGr%NLon,gGr%Nz)
        character(len=*),intent(in) :: jID

        real(rp), dimension(gGr%NLat,gGr%NLon) :: dS !Surface area at k0 [m2]
        integer :: i,j,k,jP,jM,k0
        real(rp) :: rad,theta,thP,thM,dth,dphi,DelA,DelB
        real(rp) :: jScl,jMin,jMax,jRMS,jNet,jTot

        dbJ = 0.0
        k0 = 1
        
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,jP,jM,rad,theta,thP,thM,dth,dphi,DelA,DelB)
        do k=1,gGr%Nz
            do j=1,gGr%NLon
                do i=1+1,gGr%NLat-1
                    rad  = norm2(gGr%GxyzC(i,j,k,:))
                    theta = acos(gGr%GxyzC(i,j,k,ZDIR)/rad)
                    
                    thP = acos(gGr%GxyzC(i+1,j,k,ZDIR)/rad)
                    thM = acos(gGr%GxyzC(i-1,j,k,ZDIR)/rad)
                    dth =  thP - thM
                    dphi = 2*2*PI/gGr%NLon !Assuming uniform longitude and centered difference

                    DelA = ( sin(thP)*dbRTP(i+1,j,k,PDIR) - sin(thM)*dbRTP(i-1,j,k,PDIR) )/dth
                    if (j == 1) then
                        jM = gGr%NLon
                        jP = j+1
                    else if (j == gGr%NLon) then
                        jM = j-1
                        jP = 1
                    else
                        jP = j+1
                        jM = j-1
                    endif
                    DelB = ( dbRTP(i,jP,k,TDIR) - dbRTP(i,jM,k,TDIR) )/dphi
                    dbJ(i,j,k) = (DelA-DelB)/(rad*sin(theta))

                    if (k == k0) then
                        !Calculate dS for this patch [m2]
                        dS(i,j) = (rad**2.0)*(REarth**2.0)*sin(theta)*dphi*dth
                    endif
                enddo
            enddo
        enddo

        !This is raw curl, nT/Re
        !Convert to current, ie convert nT/Re => T/m, multiply by Mu0 = 4pi x 10^-7 Tm/A
        jScl = (1.0e-9)/(Re_cgs*1.0e-2)/Mu0 !Converts to A/m2
        dbJ = (1.0e+6)*jScl*dbJ !microA/m2

        jMin = minval(dbJ(:,:,k0))
        jMax = maxval(dbJ(:,:,k0))
        jRMS = norm2 (dbJ(:,:,k0))/sqrt(1.0*size(dbJ(:,:,k0)))

        jNet = sum(dbJ(:,:,k0)     *dS)*(1.0e-6)*(1.0e-6)
        jTot = sum(abs(dbJ(:,:,k0))*dS)*(1.0e-6)*(1.0e-6)
        write(*,'(a,a,f8.3,f8.3,f8.3,a)') trim(jID),' anomalous current [microA/m2]: ',jMin,jMax,jRMS,' (Min/Max/RMS)'
        write(*,*) "   Tot/Net anomalous [MA] = ",jTot,jNet 
    end subroutine CalcJdb

end module calcdbio
