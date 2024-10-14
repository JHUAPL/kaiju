!Various routines to handle inner magnetosphere models coupled to Voltron

module innermagsphere
    use kdefs
    use gamtypes
    use ebtypes
    use volttypes
    use gamapp
!    use sstimag
    use sstLLimag
    use rcmimag
    use msphutils, only : RadIonosphere
    use rcmXimag
    use cmiutils, only : SquishCorners
    
    implicit none

    contains

    !Figure out which inner magnetosphere model we're using and initialize it
    !subroutine InitInnerMag(vApp,isRestart,iXML)
    subroutine InitInnerMag(vApp,gApp,iXML)
        type(voltApp_T)  , intent(inout) :: vApp
        !logical, intent(in) :: isRestart
        class(gamApp_T), intent(in) :: gApp
        type(XML_Input_T), intent(inout) :: iXML
        character(len=strLen) :: imStr

        if (.not. vApp%doDeep) return !Why are you even here?

        call iXML%Set_Val(imStr,"coupling/imType","SST")

        !NOTE: Using the fact that x2 is longitude and 2P periodic for both inner mag models
        select case (trim(toUpper(imStr)))
        case("SST","TS07")
            vApp%imType = IMAGSST
            vApp%prType = LPPROJ !R-phi
            allocate(empData_T :: vApp%imagApp)
        ! case("SSTLL")  ! on lon-lat grid in the ionosphere -- like RCM
        !     vApp%imType = IMAGSSTLL
        !     vApp%prType = LLPROJ !R-phi
        !     allocate(empData_T :: vApp%imagApp)
        case("RCM")
            vApp%imType = IMAGRCM
            vApp%prType = LLPROJ !Lat-lon
            allocate(rcmIMAG_T :: vApp%imagApp)
        case("RCMX")
            vApp%imType = IMAGRCMX
            vApp%prType = LLPROJ !Lat-lon
            allocate(rcmXIMAG_T :: vApp%imagApp)
        case DEFAULT
            write(*,*) 'Unkown imType, bailing ...'
            stop
        end select

        call vApp%imagApp%doInit(iXML,gApp%Model%isRestart,vApp)

    end subroutine InitInnerMag

end module innermagsphere
