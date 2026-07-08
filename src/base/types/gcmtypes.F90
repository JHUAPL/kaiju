module gcmtypes
  use shellGrid
  use voltCplTypes

  implicit none

    type gcm_grid_T
        type(ShellGrid_T) :: grid
        type(ShellGridVar_T), allocatable, dimension(:) :: gcmInput,gcmOutput
        real(rp),dimension(:),allocatable :: inlat, inlon
        real(rp),dimension(:),allocatable :: glon,gclat
        integer :: nlat,nlon,nlev,nhlat,lonshift
        character(len=32), dimension(:), allocatable :: outlist,inlist
        integer, dimension(:), allocatable :: t2N, t2S
        real(rp), dimension(:,:,:), allocatable :: invar2d, outvar2d
        integer :: mix2gcm_nvar
        integer :: gcm2mix_nvar
    end type gcm_grid_T 

end module gcmtypes
