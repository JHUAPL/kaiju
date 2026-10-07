! Routines to handle ionosphere models coupled to Voltron.
! Mirrors innermagsphere.F90 in structure.

module ionosphere
    use kdefs
    use volttypes

    implicit none

    contains

    ! Allocate and initialize the ionosphere coupler based on XML config.
    subroutine InitIonosphere(vApp, iXML)
        type(voltApp_T),     intent(inout) :: vApp
        type(XML_Input_T),   intent(inout) :: iXML

        character(len=strLen) :: ionoStr

        call iXML%Set_Val(ionoStr, "coupling/ionoType", "REMIX")

        select case (trim(toUpper(ionoStr)))
        case ("REMIX")
            allocate(remixCoupler_T :: vApp%ionoApp)
            allocate(ionoOptions_T  :: vApp%ionoApp%opt)
        case DEFAULT
            write(*,*) 'Unknown ionoType: ', trim(ionoStr), ' — bailing.'
            stop
        end select

        ! Pass Voltron grid reference into options (same pattern as imagOptions_T%voltGrid)
        vApp%ionoApp%opt%voltGrid = vApp%shGrid

        call vApp%ionoApp%InitModel(iXML)

    end subroutine InitIonosphere

end module ionosphere
