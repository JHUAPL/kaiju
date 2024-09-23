module usergamic
    use gamtypes
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use background

    implicit none

    contains
    
    subroutine initUser(Model,Grid,State,xmlInp)
        class(Model_T), intent(inout) :: Model
        class(Grid_T), intent(inout) :: Grid
        class(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        write(*,*) 'Attempting to use null user routine, bailing out!'
        stop
    end subroutine initUser

end module usergamic
