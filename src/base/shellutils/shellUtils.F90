! Heap of misc. helpful functions for ShellGrids and ShellGridVars
module shellUtils

    use shellGrid

    implicit none

    contains


    subroutine ijExtra_SGV(loc, iExtra, jExtra)
        !! Determine which dimensions have extra index relative to # cells based on variable's location on grid
        integer, intent(in) :: loc
        integer, intent(out) :: iExtra
        integer, intent(out) :: jExtra
        
        select case(loc)
            case(SHGR_CC)
                iExtra = 0
                jExtra = 0
            case(SHGR_CORNER)
                iExtra = 1
                jExtra = 1
            case(SHGR_FACE_THETA)
                iExtra = 1
                jExtra = 0
            case(SHGR_FACE_PHI)
                iExtra = 0
                jExtra = 1
            case default
                write(*,*) "initShellGridVar got an invalid data location:",loc
                stop
        end select
    end subroutine ijExtra_SGV

    subroutine wrapJ_SGV(sh, sgVar)
        !! Wrap a ShellGridVar around the periodic j/phi boundary
        !! Assumes all vars within active domain are valid to wrap
        !! Note: anything in the i-direction ghost cells get wrapped too
        type(ShellGrid_T), intent(in) :: sh
        type(ShellGridVar_T), intent(inout) :: sgVar

        associate(Q => sgVar%data)

        ! Theta faces are cell-centered w.r.t. j direction
        !! TODO: Catch for zero ghosts in j direction
        if (sgVar%loc == SHGR_CC .or. sgVar%loc == SHGR_FACE_THETA) then
            ! Starting ghost cells
            Q(:, sh%jsg:sh%js-1) = Q(:, sh%je-sh%Nge+1:sh%je)
            ! Ending ghosts cells
            Q(:, sh%je+1:sh%jeg) = Q(:, sh%js:sh%js+sh%Ngw-1)
        elseif (sgVar%loc==SHGR_CORNER .or. sgVar%loc == SHGR_FACE_PHI) then
            ! Starting ghost cells
            Q(:, sh%jsg:sh%js) = Q(:, sh%je-sh%Nge+1:sh%je+1)
            ! Ending ghosts cells
            Q(:, sh%je+1:sh%jeg+1) = Q(:, sh%js:sh%js+sh%Ngw)
        endif

        end associate

    end subroutine wrapJ_SGV


    function getSGJ_periodic(sh, j) result(jWrap)
        !! Get j index with periodic boundary applied to make sure we are always in the active domain
        type(ShellGrid_T), intent(in) :: sh
        integer, intent(in) :: j

        integer :: jWrap

        if (.not. sh%isPeriodic) then
            write(*,*)"ERROR: getSGJ_periodic called for a non-periodic shell grid"
            write(*,*)"ShGr name: shGr%name"
            stop
        endif

        jWrap = j
        ! Catch for j > Np
        do while (jWrap > sh%je)
            jWrap = jWrap - sh%Np
        enddo
        ! Catch for j < Np
        do while (jWrap < sh%js)
            jWrap = jWrap + sh%Np
        enddo

    end function getSGJ_periodic


    subroutine getSGCellILoc(shGr, t, iLoc, tLocO)
        ! Gets i location of ShellGrid cell containing coordinate with theta t
        type(ShellGrid_T), intent(in) :: shGr
        real(rp), intent(in) :: t
        integer, intent(out) :: iLoc
        real(rp), optional, intent(out) :: tLocO

        real(rp) :: tLoc, dTheta

        !! Variable is defined at center w.r.t. theta direction

        ! Catch for essentially exactly on the edge
        if (abs(t - shGr%maxGTheta) < TINY) then
            iLoc = shGr%ieg
            tLoc = shGr%thc(shGr%ieg)
        else if (abs(t - shGr%minGTheta) < TINY) then
            iLoc = shGr%isg
            tLoc = shGr%thc(shGr%isg)
        ! Next handle grid overshoot
        else if ( (t>shGr%maxGTheta) ) then                
            iLoc = shGr%ieg+ceiling((t-shGr%maxGTheta)/(shGr%th(shGr%ieg+1)-shGr%th(shGr%ieg)))
            dTheta = shGr%thc(shGr%ieg) - shGr%thc(shGr%ieg-1)
            tLoc = shGr%thc(shGr%ieg) + dTheta*(iLoc - shGr%ieg)
            !write(*,*)"getSGCellILoc: Changed my mind, no going out of bounds"
            !write(*,*)"shGr%name=",shGr%name
            !write(*,*)"th=",shGr%th
            !write(*,*)"t,iLoc,tLoc=",t,iLoc,tLoc
            !stop
        else if ( (t<shGr%minGTheta) ) then
            iLoc = shGr%isg-ceiling((shGr%minGTheta-t)/(shGr%th(shGr%isg+1)-shGr%th(shGr%isg)))
            dTheta = shGr%thc(shGr%isg+1) - shGr%thc(shGr%isg)
            tLoc = shGr%thc(shGr%isg) - dTheta*(shGr%isg - iLoc)
            !write(*,*)"getSGCellILoc: Changed my mind, no going out of bounds"
            !write(*,*)"shGr%name=",shGr%name
            !write(*,*)"th=",shGr%th
            !write(*,*)"t,iLoc,tLoc=",t,iLoc,tLoc
            !stop
        else
        ! If still here then the lat bounds are okay, find closest lat cell center
            
#if defined __INTEL_COMPILER && __INTEL_COMPILER >= 1800
            iLoc = findloc(t < shGr%th,.true.,dim=1) - 1 + shGr%isg - 1 !Extra -1 to revert one cell
#else
            iLoc = minloc( abs(shGr%thc-t),dim=1 ) - 1 + shGr%isg
#endif
            tLoc = shGr%thc(iLoc)
        endif

        if (present(tLocO)) then
            tLocO = tLoc
        endif
    end subroutine getSGCellILoc


    subroutine getSGCellJLoc(shGr, pin, jLoc, pLocO)
        type(ShellGrid_T), intent(in) :: shGr
        real(rp), intent(in) :: pin
        integer, intent(out) :: jLoc
        real(rp), optional, intent(out) :: pLocO


        real(rp) :: p, deltap, dJ

        p = modulo(pin,2*PI)

        ! note, shellGrid only implements [0,2pi] grids
        ! but do this check here in case it's needed in the future
        if ( (p>shGr%maxPhi) .or. (p<shGr%minPhi) ) then
            ! Point not on this grid, get outta here
            write(*,*) "ERROR in getShellJLoc, phi outside of bounds"
            write(*,*) p, shGr%minPhi, shGr%maxPhi
            stop
        endif

        if (shGr%isPhiUniform) then
            ! note this is faster, thus preferred
            deltap = shGr%phc(2)-shGr%phc(1)
            dJ = p/deltap
            jLoc = floor(dJ) + 1
        else
            jLoc = shGr%jsg
            do while (p > shGr%ph(jLoc+1))
                jLoc = jLoc + 1
            enddo
            !jLoc = minloc( abs(shGr%phc-p),dim=1 ) ! Find closest lat cell center
        endif

        if (present(pLocO)) then
            pLocO = shGr%phc(jLoc)
        endif
    end subroutine getSGCellJLoc


    subroutine iLocCC2Corner(shGr, t, iLocInout, iLocCornerO, tLocO)
        !! Takes ShellGrid cell with index i, and returns i index of closest corner to theta t
        !! iLocInout is the i index for reference cell (i,j) 
        !! If iLocCornerO is present, we return corner i loc through this variable and leave iLocInout alone
        !! If iLocCornerO not present, we return corner i loc by modifying iLocInout
        type(ShellGrid_T), intent(in) :: shGr
        real(rp), intent(in) :: t
        integer , intent(inout) :: iLocInout
        integer , intent(out), optional :: iLocCornerO
        real(rp), intent(out), optional :: tLocO

        integer :: iLocCorner
        real(rp) :: tLoc, dTheta
        
        if ( (t>shGr%maxTheta) ) then
            iLocCorner = shGr%ieg+1 + floor( 0.5 + (t-shGr%maxGTheta)/(shGr%th(shGr%ieg+1)-shGr%th(shGr%ieg)) )
            tLoc = shGr%th(shGr%ieg+1)  ! Just return the last available theta value
            dTheta = shGr%th(shGr%ieg) - shGr%th(shGr%ieg-1)
            tLoc = shGr%th(shGr%ieg) + dTheta*(iLocCorner - shGr%ieg)
            !write(*,*)"theta going out of bounds",t,shGr%maxGTheta
        else if ( (t < shGr%minTheta)) then
            iLocCorner = shGr%isg   - floor( 0.5 + (shGr%minGTheta-t)/(shGr%th(shGr%isg+1)-shGr%th(shGr%isg)) )
            tLoc = shGr%th(shGr%isg)
            dTheta = shGr%th(shGr%isg+1) - shGr%th(shGr%isg)
            tLoc = shGr%th(shGr%isg) - dTheta*(shGr%isg - iLocCorner)
            !write(*,*)"theta going out of bounds",t,shGr%maxGTheta
        else
            ! If still here then the lat bounds are okay, find closest lat cell corner
            if ( (shGr%th(iLocInout+1) - t) < (t - shGr%th(iLocInout)) ) then
                iLocCorner = iLocInout + 1
            else
                iLocCorner = iLocInout
            endif
            tLoc = shGr%th(iLocCorner)
        endif

        if (present(tLocO)) then
            tLocO = tLoc
        endif
        if (present(iLocCornerO)) then
            iLocCornerO = iLocCorner
        else
            iLocInout = iLocCorner
        endif

    end subroutine iLocCC2Corner


    subroutine jLocCC2Corner(shGr, p, jLocInout, jLocCornerO, pLocO)
        !! Takes ShellGrid cell with index j, and returns j index of closest corner to phi p
        !! jLocInout is the j index for reference cell (i,j) 
        !! If jLocCornerO is present, we return corner j loc through this variable and leave jLocInout alone
        !! If jLocCornerO not present, we return corner j loc by modifying jLocInout
        type(ShellGrid_T), intent(in) :: shGr
        real(rp), intent(in) :: p
        integer , intent(inout) :: jLocInout
        integer , intent(out), optional :: jLocCornerO
        real(rp), intent(out), optional :: pLocO

        integer :: jLocCorner
        real(rp) :: pLoc

        if ( (shGr%ph(jLocInout+1) - p) < (p - shGr%ph(jLocInout)) ) then
            jLocCorner = jLocInout + 1
        else
            jLocCorner = jLocInout
        endif

        if (present(pLocO)) then
            pLocO = shGr%ph(jLocCorner)
        endif
        if (present(jLocCornerO)) then
            jLocCornerO = jLocCorner
        else
            jLocInout = jLocCorner
        endif

    end subroutine jLocCC2Corner

    subroutine printShellGridInfo(shGr)
        type(ShellGrid_T), intent(in) :: shGr

        write(*,*)'------SG-------'
        write(*,*)' Name: ',trim(shGr%name)
        write(*,*)' Ghosts(N,S,W,E): ',shGr%Ngn,shGr%Ngs,shGr%Ngw,shGr%Nge
        write(*,*)' Nt,Np: ',shGr%Nt, shGr%Np
        write(*,*)' isChild:', shGr%isChild, trim(shGr%parentName)
        write(*,*)' Theta corners:'
        write(*,*)shGr%th*rad2deg
        write(*,*)' Phi corners:'
        write(*,*)shGr%ph*rad2deg
        write(*,*)'------SG-------'
    end subroutine printShellGridInfo

end module shellUtils
