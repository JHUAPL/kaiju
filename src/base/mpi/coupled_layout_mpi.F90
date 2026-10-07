! Centralized MPI world-rank layout parser for coupled Voltron runs.
module coupled_layout_mpi
    use mpi_f08
    use mpidefs
    use kdefs, ONLY: strLen
    use xml_input
    use strings

    implicit none

    integer, parameter :: maxLayoutSegments = 64

    type :: coupledLayout_T
        logical :: hasCustomLayout = .false.
        integer :: nSegments = 0
        integer, dimension(maxLayoutSegments) :: appId = 0
        integer, dimension(maxLayoutSegments) :: rankStart = 0
        integer, dimension(maxLayoutSegments) :: rankEnd = -1
        integer, dimension(maxLayoutSegments) :: rankCount = 0
        character(len=32), dimension(maxLayoutSegments) :: label = ''
    end type coupledLayout_T

contains

    function isAppIdInLayout(layout, id) result(res)
        type(coupledLayout_T), intent(in) :: layout
        integer, intent(in) :: id
        logical :: res

        integer :: seg

        res = .false.

        if (.not. layout%hasCustomLayout) then
            res = .false.
            return
        endif
        if (layout%nSegments <= 0) then
            res = .false.
            return
        endif

        do seg=1,layout%nSegments
            if (layout%appId(seg) == id) res = .true.
        enddo

    end function isAppIdInLayout

    function isAppLabelInLayout(layout, label) result(res)
        type(coupledLayout_T), intent(in) :: layout
        character(len=*), intent(in) :: label
        logical :: res

        integer :: seg

        res = .false.

        if (.not. layout%hasCustomLayout) then
            res = .false.
            return
        endif
        if (layout%nSegments <= 0) then
            res = .false.
            return
        endif

        do seg=1,layout%nSegments
            if (trim(toUpper(layout%label(seg))) == trim(toUpper(label))) res = .true.
        enddo

    end function isAppLabelInLayout

    subroutine initVoltronCoupledLayout(xmlInp, worldSize, layout)
        type(XML_Input_T), intent(inout) :: xmlInp
        integer, intent(in) :: worldSize
        type(coupledLayout_T), intent(out) :: layout

        character(len=strLen) :: layoutStr
        character(len=1) :: symbol
        integer :: i, nextRank, extRanks
        integer :: gNumI, gNumJ, gNumK, gCount
        integer :: numHelpers
        logical :: useHelpers
        logical :: sawG, sawV, sawT, sawH, sawS

        layout%hasCustomLayout = .false.
        layout%nSegments = 0
        layout%appId(:) = 0
        layout%rankStart(:) = 0
        layout%rankEnd(:) = -1
        layout%rankCount(:) = 0
        layout%label(:) = ''
        nextRank = 0
        sawG = .false.
        sawV = .false.
        sawT = .false.
        sawH = .false.
        sawS = .false.

        call xmlInp%Set_Val(layoutStr, "/Kaiju/Voltron/MPI/layout", "")
        layoutStr = trim(toUpper(layoutStr))

        if(len_trim(layoutStr) == 0) then
            layout%hasCustomLayout = .false.
            return
        endif

        layout%hasCustomLayout = .true.

        do i=1,len_trim(layoutStr)
            symbol = layoutStr(i:i)

            select case (symbol)
                case ('G')
                    if(sawG) call layoutAbort("Duplicate symbol 'G' in /Kaiju/Voltron/MPI/layout")
                    sawG = .true.
                    call xmlInp%Set_Val(gNumI, "/Kaiju/Gamera/iPdir/N", 1)
                    call xmlInp%Set_Val(gNumJ, "/Kaiju/Gamera/jPdir/N", 1)
                    call xmlInp%Set_Val(gNumK, "/Kaiju/Gamera/kPdir/N", 1)
                    gCount = gNumI*gNumJ*gNumK
                    if(gCount <= 0) call layoutAbort("Model 'G' is listed in layout but resolved rank count is <= 0")
                    call addLayoutSegment(layout, "GAMERA", gamId, gCount, nextRank)

                case ('V')
                    if(sawV) call layoutAbort("Duplicate symbol 'V' in /Kaiju/Voltron/MPI/layout")
                    sawV = .true.

                    call xmlInp%Set_Val(useHelpers, "/Kaiju/Voltron/Helpers/useHelpers", .false.)
                    call xmlInp%Set_Val(numHelpers, "/Kaiju/Voltron/Helpers/numHelpers", 0)
                    if(.not. useHelpers) numHelpers = 0
                    if(numHelpers < 0) call layoutAbort("Voltron helper count cannot be negative")

                    ! The V model always includes exactly one Voltron master rank.
                    call addLayoutSegment(layout, "VOLTRON", voltId, 1, nextRank)
                    if(numHelpers > 0) then
                        call addLayoutSegment(layout, "VOLTRON_HELPERS", helperId, numHelpers, nextRank)
                    endif

                case ('T')
                    if(sawT) call layoutAbort("Duplicate symbol 'T' in /Kaiju/Voltron/MPI/layout")
                    sawT = .true.

                    ! If/when T's XML schema changes, update this path.
                    ! For now, T rank counts must be explicitly provided here.
                    call xmlInp%Set_Val(extRanks, "/Kaiju/Voltron/MPI/T/numRanks", -1)
                    if(extRanks <= 0) then
                        call layoutAbort("Model 'T' is listed in layout but /Kaiju/Voltron/MPI/T/numRanks is missing or <= 0")
                    endif
                    call addLayoutSegment(layout, "TIEGCM", tgcmId, extRanks, nextRank)

                case ('H')
                    if(sawH) call layoutAbort("Duplicate symbol 'H' in /Kaiju/Voltron/MPI/layout")
                    sawH = .true.

                    call xmlInp%Set_Val(extRanks, "/Kaiju/Hidra/ddf", -1)
                    extRanks = 8 * extRanks
                    if(extRanks <= 0) then
                        call layoutAbort("Model 'H' is listed in layout but /Kaiju/Hidra/ddf is missing or <= 0")
                    endif
                    call addLayoutSegment(layout, "HIDRA_N", hidraNId, extRanks, nextRank)
                    call addLayoutSegment(layout, "HIDRA_S", hidraSId, extRanks, nextRank)
                
                case ('S')
                    if(sawS) call layoutAbort("Duplicate symbol 'S' in /Kaiju/Voltron/MPI/layout")
                    sawS = .true.

                    ! If/when S's XML schema changes, update this path.
                    ! For now, S rank counts must be explicitly provided here.
                    call xmlInp%Set_Val(extRanks, "/Kaiju/sami3/mpi/numWorkers", -1)
                    if(extRanks <= 0) then
                        call layoutAbort("Model 'S' is listed in layout but /Kaiju/sami3/mpi/numWorkers is missing or <= 0")
                    endif
                    call addLayoutSegment(layout, "SAMI3", sami3Id, extRanks, nextRank)

                case default
                    call layoutAbort("Unknown layout symbol '" // symbol // "' in /Kaiju/Voltron/MPI/layout")
            end select
        enddo

        if(nextRank /= worldSize) then
            write(*,'(A,I0,A,I0,A)') "MPI layout rank mismatch: layout defines ", nextRank, &
                " total ranks but MPI_COMM_WORLD has ", worldSize, "."
            call layoutAbort("Check /Kaiju/Voltron/MPI/layout and per-model rank counts.")
        endif

    end subroutine initVoltronCoupledLayout

    subroutine getLayoutAppId(layout, worldRank, appId)
        type(coupledLayout_T), intent(in) :: layout
        integer, intent(in) :: worldRank
        integer, intent(out) :: appId

        integer :: i

        appId = -1
        do i=1,layout%nSegments
            if(worldRank >= layout%rankStart(i) .and. worldRank <= layout%rankEnd(i)) then
                appId = layout%appId(i)
                return
            endif
        enddo

        call layoutAbort("Rank was not assigned by the computed MPI layout")
    end subroutine getLayoutAppId

    subroutine addLayoutSegment(layout, label, appId, rankCount, nextRank)
        type(coupledLayout_T), intent(inout) :: layout
        character(len=*), intent(in) :: label
        integer, intent(in) :: appId, rankCount
        integer, intent(inout) :: nextRank

        integer :: idx

        if(rankCount <= 0) call layoutAbort("Attempted to add layout segment with non-positive rank count")

        idx = layout%nSegments + 1
        if(idx > maxLayoutSegments) call layoutAbort("Too many layout segments; increase maxLayoutSegments")

        layout%nSegments = idx
        layout%label(idx) = label
        layout%appId(idx) = appId
        layout%rankCount(idx) = rankCount
        layout%rankStart(idx) = nextRank
        layout%rankEnd(idx) = nextRank + rankCount - 1

        nextRank = nextRank + rankCount
    end subroutine addLayoutSegment

    subroutine layoutAbort(message)
        character(len=*), intent(in) :: message
        integer :: ierr

        write(*,*) "MPI layout error: ", trim(message)
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end subroutine layoutAbort

end module coupled_layout_mpi
