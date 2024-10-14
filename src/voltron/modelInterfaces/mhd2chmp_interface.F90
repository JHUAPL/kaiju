! Converting gamera data to chimp data

module mhd2chmp_interface
    use gamtypes
    use volttypes
    use math
    use clocks
    use gamapp
    use ebtypes
    use chmpunits
    use chmpfields
    use msphutils, only : RadIonosphere
    implicit none

    contains

    subroutine init_mhd2Chmp(mhd2chmp, gamApp, ebTrcApp)
        type(mhd2Chmp_T), intent(inout) :: mhd2chmp
        class(gamApp_T)  , intent(in)    :: gamApp
        type(ebTrcApp_T), intent(inout) :: ebTrcApp

        real(rp) :: rIon

        !Set lowlat BC
        !rIon = (RionE*1.0e+6)/REarth
        rIon = RadIonosphere() ! Grabs RIon in planetary radii
        associate(Gr=>gamApp%Grid)
        !mhd2chmp%Rin = norm2(Gr%xyz(Gr%is+1,Gr%js,Gr%ks,:))
        mhd2chmp%Rin = norm2(Gr%xyz(Gr%is,Gr%js,Gr%ks,:))

        end associate
        mhd2chmp%lowlatBC = 90.0 - rad2deg*asin(sqrt(rIon/mhd2chmp%Rin)) !co-lat -> lat
        mhd2chmp%lowlatBC = mhd2chmp%lowlatBC/rad2deg
        
        write(*,*) "MAXTUBESIZE = ", MAXTUBESIZE
        
    end subroutine init_mhd2Chmp

    subroutine convertGameraToChimp(mhd2chmp,gamApp,ebTrcApp)
        type(mhd2Chmp_T), intent(inout) :: mhd2chmp
        type(gamApp_T)  , intent(in)    :: gamApp
        type(ebTrcApp_T), intent(inout) :: ebTrcApp

        real(rp), dimension(NVAR) :: pW, pCon
        real(rp), dimension(NDIM) :: Bxyz,Vxyz
        integer :: i,j,k

        associate(Gr=>gamApp%Grid,State=>gamApp%State,ebGr=>ebTrcApp%ebState%ebGr,ebF=>ebTrcApp%ebState%eb1)
        
        !Scrape Gamera MHD data into CHIMP structure
        !NOTE: Chimp and Gamera both have background fields that are not equal
        !Gamera uses a cutdipole background, CHIMP uses a pure dipole background
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,pW,pCon,Bxyz,Vxyz)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie
                    !Get conserved state from Gamera, convert to prim, scale and load into Chimp
                    pCon = State%Gas(i,j,k,:,BLK)
                    call CellC2P(gamApp%Model,pCon,pW)
                    !Add Gamera field + Gamera background, scale to CHIMP
                    Bxyz = inBScl*State%Bxyz(i,j,k,:)
                    if (gamApp%Model%doBackground) then
                        Bxyz = Bxyz + inBScl*Gr%B0(i,j,k,:)
                    endif
                    Vxyz = inVScl*pW(VELX:VELZ)
                    
                    !Now store data, subtract CHIMP background
                    ebF%dB(i,j,k,:) = Bxyz - ebGr%B0cc(i,j,k,:)
                    ebF%E(i,j,k,:)  = -cross(Vxyz,Bxyz)

                    if (ebTrcApp%ebModel%doMHD) then
                        ebF%W(i,j,k,DEN)       = inDScl*pW(DEN)
                        ebF%W(i,j,k,PRESSURE)  = inPScl*pW(PRESSURE)
                        ebF%W(i,j,k,VELX:VELZ) = Vxyz !Already scaled
                    endif
                enddo
            enddo
        enddo

        !Impose some ghost conditions
        call ebGhosts(ebTrcApp%ebModel,ebGr,ebF)

        end associate
        
        ebTrcApp%ebState%eb1%time = 0.0
        ebTrcApp%ebState%eb2%time = 1.0
        ebTrcApp%ebState%doStatic = .true.
        
    end subroutine convertGameraToChimp

end module mhd2chmp_interface

