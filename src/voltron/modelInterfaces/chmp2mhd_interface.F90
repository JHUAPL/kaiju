! Converting chimp data to gamera data

module chmp2mhd_interface
    use gamtypes
    use volttypes
    use gamapp
    use ebtypes
    
    implicit none

    contains

    subroutine init_chmp2Mhd(chmp2mhd, ebTrcApp, gamApp)
        type(chmp2Mhd_T), intent(inout) :: chmp2mhd
        type(ebTrcApp_T), intent(inout) :: ebTrcApp
        class(gamApp_T)  , intent(in)    :: gamApp

        associate(Gr=>gamApp%Grid)
        !Changing squish to do nodes
        allocate(chmp2mhd%xyzSquish(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,2))
        !Good projection or not
        allocate(chmp2mhd%isGood   (Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1))
        !Safe to ingest or not
        allocate(chmp2mhd%isEdible (Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1))

        chmp2mhd%xyzSquish = 0.0
        chmp2mhd%isGood    = .false.
        chmp2mhd%isEdible  = .false.

        if (ebTrcApp%ebModel%epsds > TINY) then
            chmp2mhd%epsds0 = ebTrcApp%ebModel%epsds
        else
            chmp2mhd%epsds0 = chmp2mhd%epsSquish
        endif

        end associate

    end subroutine init_chmp2Mhd

end module chmp2mhd_interface
