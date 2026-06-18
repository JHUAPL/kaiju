!Main CHIMP field/grid data types

module ebtypes
    use chmpdefs
    use wpitypes
    use iotable
    
    implicit none

    !Primitive variables
    integer, parameter :: NVARMHD = 5

    !Data necessary for CHIMP 2D<->3D mapping for deep coupling
    type ebSquish_T
        ! base variables
        real(rp) :: Rinner
        ! block tracking for splitting up calculations and dividing work
        integer :: numSquishBlocks = 3 ! total blocks to divide work into
        integer :: curSquishBlock = 1
        ! start and end blocks in case work is divided across multiple ranks
        integer :: myFirstBlock = 1 ! first block for me to work on
        integer :: myNumBlocks = -1 ! how many blocks I should solve
        ! dynamic block size adjustment
        integer, dimension(:), allocatable :: blockStartIndices
    end type ebSquish_T

    !Holds single slice of field data, attached to grid data in ebState
    type ebField_T
        real(rp), dimension(:,:,:,:)  , allocatable :: dB,E !Fields
        real(rp), dimension(:,:,:,:,:), allocatable :: W !Primitive MHD variables
        !W = ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NVARMHD,0:Model%nSpc
        real(rp), dimension(:,:,:,:), allocatable :: src
            !! Source variables (NVARVOLT0 and NVARIMAG0 vars)

        real(rp), dimension(:,:,:,:), allocatable :: Jxyz !Current density, must be [A/m2]
        real(rp), dimension(:), allocatable :: Lpp !Plasmapause location,Lpp(MLT)
        real(rp) :: time !Time in code units for this slice
        character(len=strLen) :: gStr ! holds step in 'Step#N'
    end type ebField_T

    !Grid data for eb fields
    type ebGrid_T
        integer :: Ni,Nj,Nk !Grid dimensions (w/ ghosts)
        integer :: Nip,Njp,Nkp !Active cells (w/o ghosts)
        integer :: is,ie,js,je,ks,ke !Grid indices (w/o ghosts)
        integer :: isg,ieg,jsg,jeg,ksg,keg !Grid indices (w ghosts)
        integer :: GrID=CARTGRID !Grid ID (see chmpdefs)
        !4D cell corners, cell centers
        !(isg:ieg,jsg:jeg,ksg:keg,1:NDIM)
        real(rp), dimension(:,:,:,:), allocatable :: xyz,xyzcc,B0cc
        !3D cell centered dV
        real(rp), dimension(:,:,:), allocatable :: dV

        !Cell-centered Jacobians.  Txi: (xyz),ijk / Tix: (ijk),xyz
        !Ie, Txi(XDIR,JDIR) = dx/dj
        !(isg:ieg,jsg:jeg,ksg:keg,1:NDIM,1:NDIM)
        real(rp), dimension(:,:,:,:,:), allocatable :: Txi,Tix
    end type ebGrid_T

    !Main field data structure
    type ebState_T
        type(ebGrid_T)  :: ebGr
        type(ebField_T) :: eb1,eb2 !Two slices of field data
        type(ioTab_T)   :: ebTab
        type(wave_T)    :: ebWave
        type(wModel_T)  :: ebWmodel
        logical :: doStatic = .false.
    end type ebState_T
    
    !Type used for tracing field lines
    type ebTrcApp_T
        type(chmpModel_T) :: ebModel
        type(ebState_T)   :: ebState
        type(ebSquish_T)  :: ebSquish
    end type ebTrcApp_T

    !Data holder for doing field line tracing at point
    type ebTrc_T
        real(rp) :: OCb !Topology
        real(rp) :: dvB,bD,bP,bS !Flux-tube volume, averaged density/pressure, integrated entropy
        real(rp) :: stdD, stdP  ! Standard deviation of density and pressure
        real(rp), dimension(NDIM) :: MagEQ, xEPm,xEPp !xyz of equator/ -/+ field endpoints
        real(rp) :: bMin !Minimum B (@ equator)
    end type ebTrc_T

    !Fields and derivatives necessary for GC update
    type gcFields_T
        real(rp), dimension(NDIM) :: DotE,DotB
        real(rp), dimension(NDIM,NDIM) :: JacE,JacB

    end type gcFields_T

    integer, parameter :: MaxFL = MAXTUBESIZE !Reduced for multi-threading speed

    !Magnetic field line type, w/ MF plasma information
    type magLine_T
        integer :: nMax ! Max number of steps this field line is allowed to have
        real(rp), dimension(NDIM) :: x0 !Seed point
        integer :: Nm=0,Np=0 !Length of line, -Nm:+Np w/ x0=>0

        !Whether this is a degenerate line (seed point not in domain) or didn't bother to trace
        logical :: isGood = .false.

        real(rp), allocatable, dimension(:,:) :: xyz
        !xyz(-Nm:+Np,XDIR:ZDIR), w/ xyz(0,XDIR:ZDIR) = seed point

        !Localization data, ie ijk of each node of field line
        !Same size as xyz array
        integer, allocatable, dimension(:,:) :: ijk

        !Magnetic field along field line, (-Nm:+Np,XDIR:ZDIR)
        real(rp), allocatable, dimension(:,:) :: Bxyz
        real(rp), allocatable, dimension(:)   :: magB !Just vector magnitude of above

        !Plasma information (primitive variables)
        !Gas = (-Nm:+Np,NVARMHD,0:Model%nSpc) w/ 0 = BULK
        real(rp), allocatable, dimension(:,:,:) :: Gas

    end type magLine_T
        
    contains

    
end module ebtypes
