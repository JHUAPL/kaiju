!Gamera definitions/constants

module gdefs
    use kdefs
    implicit none

!Constant sizes
    integer, parameter :: NVAR=5
    integer, parameter :: BLK=0
    integer, parameter :: MAXBC=6
    !Density/pressure floors
    real(rp) :: pFloor = TINY, dFloor = TINY

!Numerical knobs
    !Interpolation
    integer, parameter :: recLen = 8 !Reconstruction stencil size
    !Limiter
    integer, parameter :: limLen = 4 !Limiter stencil size
    !Vectorization, size of work buffers
    integer, parameter :: vecLen = BRICKSIZE !Taken from preprocessor definition

    !Derived values
    integer, parameter :: Nr2 = recLen/2, Nl2 = limLen/2

!Numbered accessors
    !Transformation ordering
    enum, bind(C)
        enumerator :: NORMX=1,NORMY,NORMZ
        enumerator :: TAN1X,  TAN1Y,TAN1Z
        enumerator :: TAN2X,  TAN2Y,TAN2Z
    endenum 

    !Transformation ordering (eb-edge)
    enum, bind(C)
        enumerator :: XNQI=1,YNQI,XNQJ,YNQJ
    endenum 

    !Identity tensor, useful for i,j,k displacements
    !For i direction, displacements
    !di = ijkD(IDIR,IDIR), dj = ijkD(IDIR,JDIR), ...
    integer, dimension(NDIM,NDIM), parameter :: ijkD = reshape([1,0,0,0,1,0,0,0,1],[NDIM,NDIM])

    !Default BC ordering in Halo array, user may not abide by this
    enum, bind(C)
        enumerator :: INI=1,OUTI,INJ,OUTJ,INK,OUTK
    endenum
end module gdefs
