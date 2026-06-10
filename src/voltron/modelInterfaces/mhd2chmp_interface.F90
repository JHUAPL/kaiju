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
        
        !Make sure we have J set properly since we're gonna use it
        if (.not. ebTrcApp%ebModel%doJ) then
            write(*,*) "Must set chimp/fields/doJ = T"
            stop
        endif

    end subroutine init_mhd2Chmp

    subroutine convertGameraToChimp(mhd2chmp,gamApp,ebTrcApp)
        type(mhd2Chmp_T), intent(inout) :: mhd2chmp
        class(gamApp_T)  , intent(in)    :: gamApp
        type(ebTrcApp_T), intent(inout) :: ebTrcApp

        real(rp), dimension(NVAR) :: pW, pCon
        real(rp), dimension(NDIM) :: Bxyz,Vxyz
        real (rp), dimension(:,:,:,:), allocatable :: dB,Jxyz !Full-sized arrays
        real(rp) :: jScl
        
        integer :: i,j,k,n

        associate(Model=>gamApp%Model,Gr=>gamApp%Grid,State=>gamApp%State,ebGr=>ebTrcApp%ebState%ebGr,ebF=>ebTrcApp%ebState%eb1)
                
        if (ebTrcApp%ebModel%nSpc /= gamApp%Model%nSpc) then
            write(*,*) "CHIMP/GAMERA num species mismatch!"
            stop
        endif
        
    !Calculate current [A/m2]
        call allocGridVec(Model,Gr,dB  )
        call allocGridVec(Model,Gr,Jxyz)

        !Subtract dipole before calculating current
        !$OMP PARALLEL DO default(shared) collapse(2)
        do k=Gr%ksg,Gr%keg
            do j=Gr%jsg,Gr%jeg
                do i=Gr%isg,Gr%ieg
                    dB(i,j,k,:) = State%Bxyz(i,j,k,:) - MagsphereDipole(Gr%xyzcc(i,j,k,:),Model%MagM0)
                    if (Model%doBackground) then
                        dB(i,j,k,:) = dB(i,j,k,:) + Gr%B0(i,j,k,:)
                    endif
                enddo
            enddo
        enddo

        !Turn dB => Jxyz and scale to A/m2
        jScl = (Model%gamOut%jScl)*(1.0e-9) 
        call bFld2Jxyz(Model,Gr,dB,Jxyz,jSclO=jScl)

        !Scrape Gamera MHD data into CHIMP structure
        !NOTE: Chimp and Gamera both have background fields that are not equal
        !Gamera uses a cutdipole background, CHIMP uses a pure dipole background
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,pW,pCon,Bxyz,Vxyz,n)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie
                    
                !Add Gamera field + Gamera background, scale to CHIMP
                    !Use bulk velocity to get E field
                    pCon = State%Gas(i,j,k,:,BLK)
                    call CellC2P(gamApp%Model,pCon,pW)
                    Vxyz = inVScl*pW(VELX:VELZ)

                    Bxyz = inBScl*State%Bxyz(i,j,k,:)
                    if (gamApp%Model%doBackground) then
                        Bxyz = Bxyz + inBScl*Gr%B0(i,j,k,:)
                    endif
                    
                    !Now store data, subtract CHIMP background
                    ebF%dB(i,j,k,:) = Bxyz - ebGr%B0cc(i,j,k,:)
                    ebF%E(i,j,k,:)  = -cross(Vxyz,Bxyz)

                    ebF%Jxyz(i,j,k,:) = Jxyz(i,j,k,:)

                    if (.not. ebTrcApp%ebModel%doMHD) cycle
                    
                    do n=0,gamApp%Model%nSpc
                        !Get conserved state from Gamera, convert to prim, scale and load into Chimp
                        pCon = State%Gas(i,j,k,:,n)
                        call CellC2P(gamApp%Model,pCon,pW)
                        pW(DEN)      = max(pW(DEN)     ,dFloor+TINY)
                        pW(PRESSURE) = max(pW(PRESSURE),pFloor+TINY)

                        if (Model%doMultiF .and. (n > 0)) then
                            if ( inDScl*pW(DEN) < Spcs(n)%dVac ) then
                                !This is a terrible fluid
                                pW(DEN)      = 0.0
                                pW(PRESSURE) = 0.0
                            endif
                        endif !MF check
                        ebF%W(i,j,k,DEN      ,n) = inDScl*pW(DEN)
                        ebF%W(i,j,k,PRESSURE ,n) = inPScl*pW(PRESSURE)
                        ebF%W(i,j,k,VELX:VELZ,n) = inVScl*pW(VELX:VELZ)
                                                
                    enddo !species/n
                enddo !i
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

