!Module to hold some improved EUV formulations

module euvhelper
    use kdefs
    
    implicit none

!Things for LOMPE EUV
    !Min/Max zenith angle [DEG] for LOMPE's EUV
    real(rp), parameter, private :: xMin = 0.0
    real(rp), parameter, private :: xMax = 120.0
    !Resolution of data
    real(rp), parameter, private :: dX = 1.0
    integer , parameter, private :: Nx = 121

    !Store tabulated values for q' (add paper and eqn references when published)
    !Jeff isn't gonna like this, if he asks tell him it was Slava's idea

    real(rp), parameter, private, dimension(Nx) :: lompeQP = &
        [1.000000e+00,9.998500e-01,9.994000e-01,9.986500e-01,9.976000e-01, &
         9.962510e-01,9.946030e-01,9.926570e-01,9.904130e-01,9.878720e-01, &
         9.850340e-01,9.819010e-01,9.784740e-01,9.747540e-01,9.707420e-01, &
         9.664380e-01,9.618460e-01,9.569650e-01,9.517980e-01,9.463470e-01, &
         9.406120e-01,9.345960e-01,9.283010e-01,9.217290e-01,9.148810e-01, &
         9.077610e-01,9.003700e-01,8.927110e-01,8.847860e-01,8.765980e-01, &
         8.681500e-01,8.594440e-01,8.504830e-01,8.412700e-01,8.318080e-01, &
         8.221000e-01,8.121500e-01,8.019600e-01,7.915350e-01,7.808770e-01, &
         7.699900e-01,7.588780e-01,7.475450e-01,7.359940e-01,7.242290e-01, &
         7.122560e-01,7.000760e-01,6.876960e-01,6.751190e-01,6.623500e-01, &
         6.493930e-01,6.362540e-01,6.229360e-01,6.094460e-01,5.957870e-01, &
         5.819660e-01,5.679880e-01,5.538590e-01,5.395830e-01,5.251680e-01, &
         5.106190e-01,4.959420e-01,4.811460e-01,4.662350e-01,4.512190e-01, &
         4.361030e-01,4.208970e-01,4.056090e-01,3.902480e-01,3.748240e-01, &
         3.593460e-01,3.438260e-01,3.282750e-01,3.127080e-01,2.971370e-01, &
         2.815790e-01,2.660500e-01,2.505700e-01,2.351610e-01,2.198450e-01, &
         2.046500e-01,1.896070e-01,1.747500e-01,1.601180e-01,1.457550e-01, &
         1.317130e-01,1.180480e-01,1.048230e-01,9.211030e-02,7.998570e-02, &
         6.853160e-02,5.783310e-02,4.797420e-02,3.903300e-02,3.107490e-02, &
         2.414490e-02,1.826040e-02,1.340510e-02,9.526090e-03,6.535640e-03, &
         4.318240e-03,2.741420e-03,1.668790e-03,9.722760e-04,5.412880e-04, &
         2.875250e-04,1.455210e-04,7.008060e-05,3.207010e-05,1.392600e-05, &
         5.729740e-06,2.230220e-06,8.198480e-07,2.841240e-07,9.264680e-08, &
         2.836570e-08,8.136250e-09,2.181100e-09,5.450290e-10,1.266040e-10,2.725590e-11]

    contains

!These functions take solar zenith angle (IN RADIANS) and f107 and return SigP/SigH

    elemental function SigP_EUV_LOMPE(x,f107) result(SigP)
        real(rp), intent(in) :: x,f107
        real(rp) :: SigP

        real(rp) :: qpr
        
        !Get Qp
        qpr = InterpQP(x*rad2deg)
        SigP = (f107**0.49)*(0.34*qpr + 0.93*sqrt(qpr))

    end function SigP_EUV_LOMPE

    elemental function SigH_EUV_LOMPE(x,f107) result(SigH)
        real(rp), intent(in) :: x,f107
        real(rp) :: SigH
        
        real(rp) :: qpr

        !Get Qp
        qpr = InterpQP(x*rad2deg)
        SigH = (f107**0.53)*(0.81*qpr + 0.54*sqrt(qpr))

    end function SigH_EUV_LOMPE

!Return LOMPE's q'(x) for a solar zenith angle x (in DEG)
    elemental function InterpQP(x) result(qpr)
        real(rp), intent(in) :: x
        real(rp) :: qpr

        real(rp) :: w0,w1
        integer  :: i0,i1

        if (x <= xMin) then
            qpr = 1.0
        else if (x>=xMax) then
            qpr = 0.0
        else
            !Now we do work
            i0 = floor(x/dX) + 1!Get lower interval
            i1 = i0+1

            w1 = (x - (i0-1)*dX)/dX
            w0 = 1.0 - w1

            qpr = w0*lompeQP(i0) + w1*lompeQP(i1)
        endif
        qpr = max(qpr,0.0)

    end function InterpQP

end module euvhelper

