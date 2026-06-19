!Main field interpolation routines

module ebinterp
    use chmpdefs
    use ebtypes
    use gridloc
    use gridinterp
    use ebutils

    implicit none

    logical, parameter, private :: doCurldbdt = .false.
    logical, parameter, private :: doAxisFix  = .true. !Axis fix doesn't seem to fix much

    contains

    !Interpolate field type FLD at position/time (xyz,t)
    !isInO (Optional): Is this point in the domain
    !ijkO (Optional): Guess for location finder
    !wO (Optional): Weights for this interpolation
    function fldInterp(xyz,t,Model,ebState,FLD,isInO,ijkO,wO) result(V)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        integer, intent(in) :: FLD
        integer, intent(in) , optional :: ijkO(NDIM)
        logical, intent(out), optional :: isInO
        real(rp), intent(in), optional :: wO(Nw,Nw,Nw)

        real(rp) :: V(NDIM)
        real(rp) :: dt,w1,w2,V0(NDIM)
        logical :: isIn
        integer :: ijk(NDIM)
        integer :: i0,j0,k0,n
        integer :: i1,i2,i3
        real(rp), dimension(Nw,Nw,Nw) :: W
        real(rp), dimension(Nw,Nw,Nw,NDIM) :: v1b,v2b
        real(rp),dimension(NDIM) :: ezp
        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        V = 0.0

        ! if within inner boundary of grid, use pure dipole
        if (inGap(xyz,Model,ebGr)) then
            !Handle cases here then return to avoid rest
            select case(FLD)
            case(EFLD,DBFLD)
                V = 0.0
            case(BFLD)
                V = Model%B0(xyz)
            end select
            return
        endif

        !Start by doing localization
        if (present(ijkO)) then
            !Use supplied guess
            call locate(xyz,ijk,Model,ebGr,isIn,ijkO)
        else
            call locate(xyz,ijk,Model,ebGr,isIn)
        endif

        !Set isInO if necessary
        if (present(isInO)) then
            isInO = isIn
        endif

        if (.not. isIn) then
            !Not inside domain
            !TODO: Handle E/B differently
            V = 0.0
            return
        endif

        !If still here, do actual interpolation
        !Find time weighting
        if (ebState%doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = eb2%time-eb1%time
            w1 = (eb2%time-t)/dt
            w2 = (t-eb1%time)/dt
        endif

        !Calculate spatial weighting or use provided
        if (present(wO)) then
            W = wO
        else
            call GetWeights(xyz,ijk,W,Model,ebGr)
        endif

        !Get stencils, set initial offset (ie for B0)
        V0 = 0.0
        i0 = ijk(IDIR) ; j0 = ijk(JDIR) ; k0 = ijk(KDIR)

        v2b = 0.0 !By default for static case
        select case(FLD)
        !---------------
        case(BFLD,DBFLD)
            v1b = eb1%dB(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,XDIR:ZDIR)
            if (.not. ebState%doStatic) then
                v2b = eb2%dB(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,XDIR:ZDIR)
            endif
            if (FLD == BFLD) V0 = Model%B0(xyz)
        case(EFLD)
            v1b = eb1% E(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,XDIR:ZDIR)
            if (.not. ebState%doStatic) then
                v2b = eb2% E(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,XDIR:ZDIR)
            endif
        end select

        !Now accumulate weighted average
        do n=1,NDIM
            !V(n) = V0(n) + w1*sum(W*v1b(:,:,:,n)) + w2*sum(W*v2b(:,:,:,n))
            V(n) = V0(n)
            do i3=1,Nw
                do i2=1,Nw
                    do i1=1,Nw
                        V(n) = V(n) + w1*W(i1,i2,i3)*v1b(i1,i2,i3,n) + w2*W(i1,i2,i3)*v2b(i1,i2,i3,n)
                    enddo
                enddo
            enddo
        enddo
        
        end associate
    end function fldInterp

    !Optimized routine to get multiple fields all at once
    !Always return E,B
    !Optional: vExB, gcFields (fields/derivatives necessary for GC update)

    !ijkO is optional guess to pass to locate routine
    !     if present will return correct location
    recursive subroutine ebFields(xyz,t,Model,ebState,E,B,ijkO,vExB,gcFields)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), dimension(NDIM), intent(out) :: E,B
        real(rp), dimension(NDIM), intent(out), optional :: vExB
        type(gcFields_T), intent(out), optional :: gcFields
        integer, intent(inout) , optional :: ijkO(NDIM)

        integer, dimension(NDIM) :: ijk
        integer :: i,j,k,n,m, i0,j0,k0
        real(rp) :: dt,wT1,wT2
        real(rp), dimension(NDIM,NDIM) :: Tix
        real(rp), dimension(NDIM) :: ezp,B0,wE,wZ,wP,wEp,wZp,wPp
        real(rp), dimension(Nw,Nw,Nw,NDIM) :: E1,E2,dB1,dB2 !Interpolation stencils
        real(rp), dimension(Nw,Nw,Nw) :: W,eW,zW,pW !Interpolation weights
        logical :: isIn,doJacob,isAxis,isAxisS,isAxisE
        
        !Handling axis
        integer :: ip,jp,kp,ijkAx(NDIM)
        type(gcFields_T) :: gcFieldsAxP,gcFieldsAxM
        real(rp), dimension(NDIM) :: Xp,Xm,Xc,AxE,AxB
        real(rp) :: wAx
        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )

    !Initialize fields
        E = 0.0
        
        B0 = Model%B0(xyz) !Get background B
        B = B0

        if (present(vExB)) vExB = 0.0
        if (present(gcFields)) then
            gcFields%DotE = 0.0
            gcFields%DotB = 0.0
            gcFields%JacE = 0.0
            gcFields%JacB = Model%JacB0(xyz) !Background Jacobian
            doJacob = .true.
        else
            doJacob = .false.
        endif

        ! if within inner boundary of grid, use pure dipole
        if (inGap(xyz,Model,ebGr)) return


    !Start w/ localization
        if (present(ijkO)) then
            call locate(xyz,ijk,Model,ebGr,isIn,ijkO)
            !Return correct location
            ijkO = ijk
        else
            call locate(xyz,ijk,Model,ebGr,isIn)
        endif

        !Bail (returning 0 fields) if not in grid
        if (.not. isIn) return
        i0 = ijk(IDIR) ; j0 = ijk(JDIR) ; k0 = ijk(KDIR)

    !Trap for special cases (i.e. axis)
        isAxisS = .false.
        isAxisE = .false.
        if (doAxisFix .and. (ebGr%GrID == EGGGRID .or. ebGr%GrID == LFMGRID) )then
            !Check for +x
            if ( j0 <= ebGr%js+1 ) then !js/js+1
                isAxisS = .true.
            endif
            if ( j0 >= ebGr%je-1 ) then
                isAxisE = .true.
            endif
        endif

        isAxis = isAxisS .or. isAxisE

    !Get mapping and weights
        !Map to ezp
        ezp = Map2ezp(xyz,ijk,Model,ebGr)

        !Get time weights
        if (ebState%doStatic) then
            wT1 = 1.0
            wT2 = 0.0
        else
            dt = eb2%time-eb1%time
            wT1 = (eb2%time-t)/dt
            wT2 = (t-eb1%time)/dt
        endif

        !Get 1D spatial weights
        wE = Wgt1D(ezp(IDIR))
        wZ = Wgt1D(ezp(JDIR))
        wP = Wgt1D(ezp(KDIR))

        !Turn 1D weights into 3D weights
        do k=1,Nw
            do j=1,Nw
                do i=1,Nw
                    W(i,j,k) = wE(i)*wZ(j)*wP(k)
                enddo
            enddo
        enddo


    !Do E-B fields
        !Pull stencils
        dB1 = eb1%dB(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,:)
        E1  = eb1%E (i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,:)
        if (.not. ebState%doStatic) then
            dB2 = eb2%dB(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,:)
            E2  = eb2%E (i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,:)
        else
            dB2 = 0.0
            E2  = 0.0
        endif

        do n=1,NDIM
            B(n) = wT1*sum(W*dB1(:,:,:,n)) &
                 + wT2*sum(W*dB2(:,:,:,n)) &
                 + B0(n)
            E(n) = wT1*sum(W* E1(:,:,:,n)) &
                 + wT2*sum(W* E2(:,:,:,n))

        enddo

        !Clean E field (guarantee E.B = 0)
        if (Model%doEBFix) then
            call CleanE(E,B)
        endif
        
    !Do ExB field if necessary
        if (present(vExB)) then
            vExB = EBDrift(E,B)
        endif

    !Do star fields if necessary
        if (doJacob) then
        !Do Jacobians and time derivatives
            associate( JacB=>gcFields%JacB, JacE=>gcFields%JacE )

            if (isAxis) then
                !If on axis then do some trickery
                Xc = ebGr%xyzcc(i0,j0,k0,XDIR:ZDIR)
                if (isAxisS) then
                !Positive displacement
                    ip = i0;jp = ebGr%js+2;kp = k0
                    Xp = ebGr%xyzcc(ip,jp,kp,XDIR:ZDIR)
                    ijkAx = [ip,jp,kp]
                    call ebFields(Xp,t,Model,ebState,AxE,AxB,ijkAx,gcFields=gcFieldsAxP)
                !Negative displacement
                    call ijk2Active(Model,ebGr,i0,ebGr%js-3,k0,ip,jp,kp)
                    Xm = ebGr%xyzcc(ip,jp,kp,XDIR:ZDIR)
                    ijkAx = [ip,jp,kp]
                    call ebFields(Xm,t,Model,ebState,AxE,AxB,ijkAx,gcFields=gcFieldsAxM)
                !Weight for P (closer)
                    wAx = norm2(Xm-Xc)/norm2(Xp-Xm)

                else !Negative axis
                !Positive displacement, flipping positive (to closer point)   
                    ip = i0;jp = ebGr%je-2;kp = k0
                    Xp = ebGr%xyzcc(ip,jp,kp,XDIR:ZDIR)
                    ijkAx = [ip,jp,kp]
                    call ebFields(Xp,t,Model,ebState,AxE,AxB,ijkAx,gcFields=gcFieldsAxP)
                !Negative displacement
                    call ijk2Active(Model,ebGr,i0,ebGr%je+3,k0,ip,jp,kp)
                    Xm = ebGr%xyzcc(ip,jp,kp,XDIR:ZDIR)
                    ijkAx = [ip,jp,kp]
                    call ebFields(Xm,t,Model,ebState,AxE,AxB,ijkAx,gcFields=gcFieldsAxM)
                !Weight for P (closer)
                    wAx = norm2(Xm-Xc)/norm2(Xp-Xm)
                endif !isAxisS
                
                !Interpolate dJac across the axis
                !Add JacB0 at true point
                JacE = wAx*gcFieldsAxP%JacE + (1-wAx)*gcFieldsAxM%JacE
                JacB =     wAx *( gcFieldsAxP%JacB - Model%JacB0(Xp) ) + &
                        (1-wAx)*( gcFieldsAxM%JacB - Model%JacB0(Xm) ) + &
                                  Model%JacB0(xyz)
            else
                !Otherwise do standard thing

                !Jacobians
                !---------
                !Get 1D weights for Jacobians            
                wEp = Wgt1Dp(ezp(IDIR))
                wZp = Wgt1Dp(ezp(JDIR))
                wPp = Wgt1Dp(ezp(KDIR))

                !Turn 1D weights into 3D weights
                do k=1,Nw
                    do j=1,Nw
                        do i=1,Nw
                            !Partial derivatives of weights wrt eta,zeta,psi
                            eW(i,j,k) = wEp(i)*wZ (j)*wP (k)
                            zW(i,j,k) = wE (i)*wZp(j)*wP (k)
                            pW(i,j,k) = wE (i)*wZ (j)*wPp(k)
                        enddo
                    enddo
                enddo

                !Calculate Jacobians
                !JacA(i,j) = d B_Xi / dXj
                !Tix(i0,j0,k0,ezp,xyz) = ezp derivs wrt xyz

                !Pull metric terms
                Tix = ebGr%Tix(i0,j0,k0,:,:)

                !Do main calculation
                do m=1,NDIM !Derivative direction (x,y,z)
                    do n=1,NDIM !Vector component
                        JacB(n,m) = wT1*( Tix(IDIR,m)*sum(eW*dB1(:,:,:,n))   &
                                         +Tix(JDIR,m)*sum(zW*dB1(:,:,:,n))   &
                                         +Tix(KDIR,m)*sum(pW*dB1(:,:,:,n)) ) &
                                  + wT2*( Tix(IDIR,m)*sum(eW*dB2(:,:,:,n))   & 
                                         +Tix(JDIR,m)*sum(zW*dB2(:,:,:,n))   &
                                         +Tix(KDIR,m)*sum(pW*dB2(:,:,:,n)) )

                        JacE(n,m) = wT1*( Tix(IDIR,m)*sum(eW* E1(:,:,:,n))   &
                                         +Tix(JDIR,m)*sum(zW* E1(:,:,:,n))   &
                                         +Tix(KDIR,m)*sum(pW* E1(:,:,:,n)) ) &
                                  + wT2*( Tix(IDIR,m)*sum(eW* E2(:,:,:,n))   & 
                                         +Tix(JDIR,m)*sum(zW* E2(:,:,:,n))   &
                                         +Tix(KDIR,m)*sum(pW* E2(:,:,:,n)) )


                    enddo
                enddo
                !Add background to dJacB
                JacB = Div0Jac(JacB) + Div0Jac(Model%JacB0(xyz))
                
            endif !isAxis

            !Time derivatives
            !Either use Curl(E) or linear derivative for bdot
            !----------------
            if (ebState%doStatic) then
                gcFields%DotE = 0.0
                gcFields%DotB = 0.0
            else
                do n=1,NDIM
                    gcFields%DotE(n) = (1/dt)*( sum(W* E2(:,:,:,n)) &
                                              - sum(W* E1(:,:,:,n)) )
                    gcFields%DotB(n) = (1/dt)*( sum(W*dB2(:,:,:,n)) &
                                              - sum(W*dB1(:,:,:,n)) )
                enddo
                !Replace with CurlE if option says to
                if (doCurldbdt) then
                    gcFields%DotB = -Jac2Curl(JacE)
                endif
            endif
            
            end associate !Jacobians
        endif !doJacob
    
        end associate !Main associate
        
    end subroutine ebFields

    !Interpolate all bulk MHD variables @ (xyz,t) and return
    !Optionally accepts guess for ijk localization
    function mhdInterp(xyz,t,Model,ebState,ijkO) result(iQ)
        real(rp)         , intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        integer, intent(in) , optional :: ijkO(NDIM)
        
        real(rp) :: iQ(NVARMHD)
        logical :: isIn
        real(rp), dimension(Nw,Nw,Nw) :: Wijk
        real(rp), dimension(Nw,Nw,Nw,NVARMHD) :: Q1b,Q2b
        integer :: i0,j0,k0,n
        integer :: ijk(NDIM)
        real(rp) :: dt,w1,w2

        associate( ebGr=>ebState%ebGr,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        iQ = 0.0 !Interpolated values

        ! if within inner boundary set everything to 0 and exit
        if (inGap(xyz,Model,ebGr)) return
        if (.not. Model%doMHD)     return

        if (present(ijkO)) then
            !Use supplied guess
            call locate(xyz,ijk,Model,ebGr,isIn,ijkO)
        else
            call locate(xyz,ijk,Model,ebGr,isIn)
        endif

        if (.not. isIn) return

        !Get time weighting
        if (ebState%doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = eb2%time-eb1%time
            w1 = (eb2%time-t)/dt
            w2 = (t-eb1%time)/dt
        endif

        !Get weight coefficients
        call GetWeights(xyz,ijk,Wijk,Model,ebGr)

        !Get stencils
        i0 = ijk(IDIR) ; j0 = ijk(JDIR) ; k0 = ijk(KDIR)
        Q1b = eb1%W(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,1:NVARMHD,BLK)
        if (.not. ebState%doStatic) then
            Q2b = eb2%W(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,1:NVARMHD,BLK)
        else
            Q2b = 0.0
        endif

        !Do contraction of weights and stencils via accumulation
        do n=1,NVARMHD
            iQ(n) = w1*sum(Wijk*Q1b(:,:,:,n)) + w2*sum(Wijk*Q2b(:,:,:,n))
        enddo

        end associate

    end function mhdInterp


    !Interpolate moments (D,P) for a specific species
    subroutine mhdInterpDP(xyz,t,Model,ebState,D,P,sOpt)
        real(rp)         , intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        real(rp), intent(out) :: D,P
        integer , intent(in ), optional :: sOpt

        integer :: s0,i0,j0,k0,ijk(NDIM)
        logical :: isIn
        real(rp), dimension(Nw,Nw,Nw) :: Wijk,Db,Pb
        logical , dimension(Nw,Nw,Nw) :: isG
        real(rp) :: dt,w1,w2,wScl
        logical , parameter :: doRescale = .false.

        associate( ebGr=>ebState%ebGr,eb1=>ebState%eb1,eb2=>ebState%eb2 )
        D = 0.0
        P = 0.0
        if (present(sOpt)) then
            if (sOpt > Model%nSpc) then
                !Got nothing for ya
                return
            else
                s0 = sOpt
            endif
        else
            s0 = BLK
        endif

        if ( inGap(xyz,Model,ebGr) .or. (.not. Model%doMHD) ) return 

        call locate(xyz,ijk,Model,ebGr,isIn)
        if (.not. isIn) return

        !Get time weighting
        if (ebState%doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            write(*,*) "This won't work!"
            stop
            dt = eb2%time-eb1%time
            w1 = (eb2%time-t)/dt
            w2 = (t-eb1%time)/dt
        endif

        !Get weight coefficients
        call GetWeights(xyz,ijk,Wijk,Model,ebGr)

        !Get stencils
        i0 = ijk(IDIR) ; j0 = ijk(JDIR) ; k0 = ijk(KDIR)
        Db = eb1%W(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,DEN     ,s0)
        Pb = eb1%W(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,PRESSURE,s0)

        isG = (Db > TINY) !Check for good values
        !If no good values or cell point is in (center cell) is bad get outta here
        if ( (.not. any(isG)) .or. (.not. isG(2,2,2)) ) return

        if (doRescale) then
            !If still here, rescale weights to account for what we have
            wScl = sum(Wijk,mask=isG)
            if (wScl<TINY) return

            where (isG)
                Wijk = Wijk/wScl
            elsewhere
                Wijk = 0.0
                Db   = 0.0
                Pb   = 0.0
            endwhere
        endif
        
        D = sum(Wijk*Db)
        P = sum(Wijk*Pb)

        end associate
    end subroutine mhdInterpDP

    !Interpolate all species all MHD variables @ (xyz,t) and return
    !Optionally accepts guess for ijk localization
    function mhdInterpMF(xyz,t,Model,ebState,ijkO) result(iQ)
        real(rp)         , intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        integer, intent(in) , optional :: ijkO(NDIM)

        real(rp) :: iQ(NVARMHD,0:Model%nSpc)
        logical :: isIn
        real(rp), dimension(Nw,Nw,Nw) :: Wijk
        real(rp), dimension(Nw,Nw,Nw,NVARMHD) :: Q1b,Q2b
        integer :: i0,j0,k0,n,s
        integer :: ijk(NDIM)
        real(rp) :: dt,w1,w2

        associate( ebGr=>ebState%ebGr,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        iQ = 0.0 !Interpolated values
        ! if within inner boundary set everything to 0 and exit
        if (inGap(xyz,Model,ebGr)) return
        if (.not. Model%doMHD)     return

        if (present(ijkO)) then
            !Use supplied guess
            call locate(xyz,ijk,Model,ebGr,isIn,ijkO)
        else
            call locate(xyz,ijk,Model,ebGr,isIn)
        endif

        if (.not. isIn) return

        !Get time weighting
        if (ebState%doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = eb2%time-eb1%time
            w1 = (eb2%time-t)/dt
            w2 = (t-eb1%time)/dt
        endif

        !Get weight coefficients
        call GetWeights(xyz,ijk,Wijk,Model,ebGr)

        !Get stencils
        i0 = ijk(IDIR) ; j0 = ijk(JDIR) ; k0 = ijk(KDIR)

        !Loop over species
        do s=0,Model%nSpc
            Q1b = eb1%W(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,1:NVARMHD,s)
            if (.not. ebState%doStatic) then
                Q2b = eb2%W(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,1:NVARMHD,s)
            else
                Q2b = 0.0
            endif

            !Do contraction of weights and stencils via accumulation
            do n=1,NVARMHD
                iQ(n,s) = w1*sum(Wijk*Q1b(:,:,:,n)) + w2*sum(Wijk*Q2b(:,:,:,n))
            enddo !nvar

        enddo !spc

        end associate

    end function mhdInterpMF

    !Interpolate Jxyz vector @ (xyz,t) and return
    !Optionally accepts guess for ijk localization
    !TODO: Remove all the redundant code betweent his and mhdInterp above
    function jInterp(xyz,t,Model,ebState,ijkO) result(iJxyz)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        integer, intent(in) , optional :: ijkO(NDIM)
        
        real(rp) :: iJxyz(NDIM)
        logical :: isIn
        real(rp), dimension(Nw,Nw,Nw) :: Wijk
        real(rp), dimension(Nw,Nw,Nw,NDIM) :: Q1b,Q2b
        integer :: i0,j0,k0,n
        integer :: ijk(NDIM)
        real(rp) :: dt,w1,w2

        associate( ebGr=>ebState%ebGr,eb1=>ebState%eb1,eb2=>ebState%eb2 )

        iJxyz = 0.0 !Interpolated values

        ! if within inner boundary set everything to 0 and exit
        if (inGap(xyz,Model,ebGr)) return
        if (.not. Model%doJ)       return

        if (present(ijkO)) then
            !Use supplied guess
            call locate(xyz,ijk,Model,ebGr,isIn,ijkO)
        else
            call locate(xyz,ijk,Model,ebGr,isIn)
        endif

        if (.not. isIn) return

        !Get time weighting
        if (ebState%doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = eb2%time-eb1%time
            w1 = (eb2%time-t)/dt
            w2 = (t-eb1%time)/dt
        endif

        !Get weight coefficients
        call GetWeights(xyz,ijk,Wijk,Model,ebGr)

        !Get stencils
        i0 = ijk(IDIR) ; j0 = ijk(JDIR) ; k0 = ijk(KDIR)
        Q1b = eb1%Jxyz(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,1:NDIM)
        if (.not. ebState%doStatic) then
            Q2b = eb2%Jxyz(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,1:NDIM)
        else
            Q2b = 0.0
        endif

        !Do contraction of weights and stencils via accumulation
        do n=1,NDIM
            iJxyz(n) = w1*sum(Wijk*Q1b(:,:,:,n)) + w2*sum(Wijk*Q2b(:,:,:,n))
        enddo

        end associate
    end function jInterp


    !Get time weights for time t (assuming proper bracketing)
    subroutine GetTWgts(Model,ebState,t,w1,w2)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in)  :: t
        real(rp), intent(out) :: w1,w2

        real(rp) :: dt
        !Start by getting time weight
        if (ebState%doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = ebState%eb2%time-ebState%eb1%time
            w1 = (ebState%eb2%time-t)/dt
            w2 = (t-ebState%eb1%time)/dt
        endif

    end subroutine GetTWgts
    
end module ebinterp