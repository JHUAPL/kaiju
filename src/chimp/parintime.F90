module parintime
	use chmpdefs
	use xml_input
	use chmpunits

	implicit none
	integer, private :: NumB = 0
	logical, private :: doParInT = .false. !// in time

	contains

	subroutine InitParInTime(Model,inpXML,OutID,OutF)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: inpXML
        character(len=*) , intent(in) :: OutID
        character(len=strLen), intent(out) :: OutF

        integer :: NumO,n,dn,i0,i1
        real(rp), dimension(:), allocatable :: dtOuts

        !real(rp) :: dtB,T0
    !Possible // in time
        call inpXML%Set_Val(NumB,'parintime/NumB',NumB)
        if (NumB > 1) then
            doParInT = .true.
            if ( (Model%Nblk>NumB) .or. (Model%Nblk<1) ) then
                write(*,*) "This block outside of acceptable bounds"
                write(*,*) "Block = ",Model%Nblk
                write(*,*) "Bounds = ",1,NumB
                write(*,*) "Bailing ..."
                stop
            endif
            !Reset time bounds
            write(*,*) '------'
            write(*,*) 'Resetting T0/TFin = ',Model%T0*oTScl,Model%tFin*oTScl

            NumO = floor( (Model%tFin-Model%T0)/Model%dtOut ) + 1
            !To avoid dumb divisibility things, just make the entire damn array and pick out a chunk
            allocate(dtOuts(NumO))

            !Create total list of output times
            do n=1,NumO
                dtOuts(n) = Model%T0 + (n-1)*Model%dtOut
            enddo

            !Pick the range that this block is responsible for
            !Proportional split, block sizes differ by at most 1
            i0 = ((Model%Nblk-1)*NumO)/NumB + 1
            i1 = ( Model%Nblk   *NumO)/NumB
            Model%T0 = dtOuts(i0)
            Model%tFin = dtOuts(i1) + 0.01*Model%dtOut !Add a smidge to the end to make sure it gets output

            write(*,*) 'To        T0/TFin = ',Model%T0*oTScl,Model%tFin*oTScl

            !Don't bother offsetting step #, let the concatenating script handle it
            write(OutF,'(a,a,I0.4,a,a,a)') trim(adjustl(Model%RunID)),'.',Model%Nblk,'.',trim(adjustl(OutID)),'.h5'
            write(*,*) '------'            
        else
            doParInT = .false.
            NumB = 0
            write(OutF,'(4a)') trim(adjustl(Model%RunID)),'.',trim(adjustl(OutID)),'.h5'
        endif
        write(*,*) "Writing output to ", trim(OutF)
	end subroutine InitParInTime


end module parintime