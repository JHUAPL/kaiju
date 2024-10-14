!Main BC data types

module gambctypes
    use gdefs
    use gamtypes

    implicit none

    ! derived types for BCs that apply to specific external directions
    type, extends(baseBC_T) :: innerIBC_T
        contains
        procedure :: bcDir => innerIDir
    end type
    type, extends(baseBC_T) :: outerIBC_T
        contains
        procedure :: bcDir => outerIDir
    end type
    type, extends(baseBC_T) :: innerJBC_T
        contains
        procedure :: bcDir => innerJDir
    end type
    type, extends(baseBC_T) :: outerJBC_T
        contains
        procedure :: bcDir => outerJDir
    end type
    type, extends(baseBC_T) :: innerKBC_T
        contains
        procedure :: bcDir => innerKDir
    end type
    type, extends(baseBC_T) :: outerKBC_T
        contains
        procedure :: bcDir => outerKDir
    end type

    ! derived type for BCs that can apply to multiple external directions
    ! these require extra initialization by the user, either via
    !    explicit setting of the member variable after initialziation, or
    !    sourced allocation with a default structure constructor
    type, extends(baseBC_T) :: multiBC_T
        integer :: dir = -1
        contains
        procedure :: bcDir => multiDir
    end type

    contains

    ! BC direction functions, to be used by extension BC types
    function innerIDir(bc)
        class(innerIBC_T), intent(in) :: bc
        integer :: innerIDir
        innerIDir = INI
    end function

    function outerIDir(bc)
        class(outerIBC_T), intent(in) :: bc
        integer :: outerIDir
        outerIDir = OUTI
    end function

    function innerJDir(bc)
        class(innerJBC_T), intent(in) :: bc
        integer :: innerJDir
        innerJDir = INJ
    end function

    function outerJDir(bc)
        class(outerJBC_T), intent(in) :: bc
        integer :: outerJDir
        outerJDir = OUTJ
    end function

    function innerKDir(bc)
        class(innerKBC_T), intent(in) :: bc
        integer :: innerKDir
        innerKDir = INK
    end function

    function outerKDir(bc)
        class(outerKBC_T), intent(in) :: bc
        integer :: outerKDir
        outerKDir = OUTK
    end function

    function multiDir(bc)
        class(multiBC_T), intent(in) :: bc
        integer :: multiDir
        multiDir = bc%dir
    end function

end module gambctypes
