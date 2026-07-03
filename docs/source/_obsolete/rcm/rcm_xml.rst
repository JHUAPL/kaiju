XML input files for RCM
=======================

Introduction
------------

These XML elements and attributes were extracted from ``src/rcm/modules.F90``,
``src/rcm/rcm_mhd_io.F90``, and ``src/rcm/rcm_subs.F90``.

``<Kaiju>``

``<RCM>``

``<advect>``

``doSmoothDDV`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to
smooth ij derivatives of residual flux tube volume (FTV).

``epsPk`` (optional, float, default ``"1.0e-3"``): Cumulative pressure
fraction threshold used to limit number of energy channels evolved

``<chargex>``

``kill_fudge`` (optional, Boolean, default ``"F"``: Set to ``"T"`` means no
loss.

``L_dktime`` (optional, Boolean, default ``"T"``: Set to ``"T"`` to read the
loss from the dktime table.

``sunspot_number`` (optional, float, default ``"96.0"``: No longer in use.

``<clawpack>``

``doOMPClaw`` (optional, Boolean, default ``"T"``: Set to ``"T"`` to use
OpenMP threading on the clawpack solver.

``<eflux>``

``icorrect`` (optional, Boolean, default ``"T"``: Set to ``"T"`` make lat.
correction to EFLUX.

``ifloor`` (optional, Boolean, default ``"T"``: Set to ``"T"`` to install a
floor for EFLUX.

``<ellipse>``

``isDynamic`` (optional, Boolean, default ``"T"``):Set to ``"T"`` to use the
ellipse boundary on RCM.

``xSun`` (optional, float, default ``"12.5"``): The Sun-Earth distance.

``xTail`` (optional, float, default ``"-15.0"``): The Earth-tail distance.

``yDD`` (optional, float, default ``"15.0"``): Positive/negative y-axis
bounds.

``<experimental>``

``doNoBndFlow`` (optional, Boolean, default ``"F"``: Option to restrict
inward boundary flow. Recommended to be left as false

``NBFLayers`` (optional, integer, default ``"2"``: Number of cells from the
open-closed boudnary to restrict inward flow over. Only used if
doNoBndFlow="T"

``<grid>``

``doLatStretch`` (optional, Boolean, default ``"F"``): Set to ``"T"`` to use
non-uniform RCM grid

``HiLat`` (optional, float, default ``"75.0"``): High-latitude grid boundary
(in degrees)

``LowLat`` (optional, float, default ``"30.0"``): Low-latitude grid boundary
(in degrees)

``ibnd_type`` (optional, integer, default ``"4"``): Type of bndy (1-eq.p,
2-iono).

``ipcp_type`` (optional, integer, default ``"13"``): Type of bndy (1-eq.p,
2-iono).

`L_move_plasma_grid`` (optional, Boolean, default ``"T"``): No longer in use

``nsmthi`` (optional, integer, default ``"0"``): How much to smooth cond in
I

``nsmthj`` (optional, integer, default ``"0"``): How much to smooth cond in
J

``<loss>``

``doFLCLoss`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to use
FLC losses.

``doNewCX`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to use
newer CX loss estimate.

``doRelax`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to relax
energy distribution.

``doTDSLoss`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to use
TDS losses.

``eLossMethod`` (optional, string, default ``"FDG"``: Choose the electron
loss method within FDG,  WM, SS, C05,C19

``<output>``

``debug`` (optional, Boolean, default ``"F"``): Enable debug console output

``doDebug`` (optional, Boolean, default ``"F"``): Set to ``"T"`` to print
out diagnostic messages

``doDebugH5`` (optional, Boolean, default ``"F"``): Enable debug output to
hdf5 file. Greatly increases file size

``doFLOut`` (optional, Boolean, default ``"F"``): Set to ``"T"`` to output
field lines (slow).

``iDebug`` (optional, integer, default ``"1"``): Set to ``0`` to do disk
printout.

``nSkipFL`` (optional, integer, default ``"8"``): Stride for outputting
field lines (UNITS?)

``toMHD`` (optional, Boolean, default ``"F"``): No longer in use

``toRCM`` (optional, Boolean, default ``"F"``): No longer in use

``<params>``

``cmax`` (optional, float, default ``"3.0"``): In rcm_mod_balgn.

``eeta_cutoff`` (optional, float, default ``"3.0"``): As a fraction.

``ipot`` (optional, integer, default ``"-1"``): Which potential solver to
use.

``icond`` (optional, integer, default ``"3"``): 1 is active conductances, 2
is Hardy with kp, 3 is input.

``iwind`` (optional, integer, default ``"0"``): 0 is no neutral winds.

``<plasmasphere>``

``DenPP0`` (optional, float, default ``"0.0"``): Plasmasphere density
cutoff, [#/cc]

``doPPSmooth`` (optional, Boolean, default ``"T"``): Try to smooth
plasmapause

``doRefill`` (optional, Boolean, default ``"F"``): Set to ``"T"`` to refill
plasmasphere.

``initKp`` (optional, integer, default ``"1"``): Initial Kp condition for
the Gallagher model

``isDynamic`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to use
the dynamic plasmasphere

``staticR`` (optional, float, default ``"0.0"``): Set the static part (in
radius) of the plasmasphere

``tAvg`` (optional, float, default ``"0.0"``): Averaging timescale applied to
the electrostatic potential used to evolve the plasmasphere (seconds).

``<restart>``

``nRes`` (optional, integer, default ``"-1"``): Restart number.

``resID`` (optional, string, default ``"msphere"``): Run ID for restart.

``<sim>``

``nSubstep`` (optional, integer, default ``"4"``): Number of substeps in each
MHD-RCM coupling cycle.

``runid`` (optional, string, default ``"MAGE sim"``): ID string to use for
the run

``<tilt>``

``isTilt`` (optional, Boolean, default ``"F"``): Set to ``"T"`` to tilt the
dipole, must turn off corotation also.

``<tomhd>``

``doAvg2MHD`` (optional, Boolean, default ``"T"``): Determines whether the
moments given back to MHD are from instantaneous values or averaged over the]
coupling duration

``doQ0`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to include
implicit cold ions in tomhd moments.

``doRelax`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to relax
energy distribution.

``<torcm>``

``doKappa`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to do kappa
by default.

``doRescale`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to rescale
D, P => eta by default.

``doSmoothBNDLOC`` (optional, Boolean, default ``"T"``): Set to ``"T"`` to do
bndloc smoothing.
