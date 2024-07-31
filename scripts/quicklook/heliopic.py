#!/usr/bin/env python


"""Make a quick-look plot for a gamhelio run.

Make a quick-look plot for a gamhelio run.

Several different sets of plots are supported, and are distinguished by the
value of the "pic" argument.

pic1 (default): A 4-panel display showing pcolormesh plots in the z = 0
(equatorial) plane of the gamhelio frame used in the simulation. The plots
are (r0 is the inner radius of the grid, which should be 21.5 Rsun):

    Upper left: Solar wind speed (km/s)
    Upper right: Solar wind number density scaled by (r/r0)**2 (cm**-3)
    Lower left: Solar wind temperature scaled by r/r0 (MK)
    Lower right: Solar wind radial magnetic field scaled by (r/r0)**2 (nT)

pic2: A 4-panel display showing pcolormesh plots in the y = 0 (meridional,
containing Earth) plane of the gamhelio frame used in the simulation. The
plots are (r0 is the inner radius of the grid, which should be 21.5 Rsun):

    Upper left: Solar wind speed (km/s)
    Upper right: Solar wind number density scaled by (r/r0)**2 (cm**-3)
    Lower left: Solar wind temperature scaled by r/r0 (MK)
    Lower right: Solar wind radial magnetic field scaled by (r/r0)**2 (nT)

pic3: A 4-panel display showing pcolormesh plots in the r = 1 AU slice of the
gamhelio frame used in the simulation. The plots are:

    Upper left: Solar wind speed (km/s)
    Upper right: Solar wind number density (cm**-3)
    Lower left: Solar wind temperature (MK)
    Lower right: Solar wind radial magnetic field (nT)

pic4: A pcolormesh plot in the innermost radial slice (r = 22 Rsun) of the
gamhelio frame used in the simulation. The plot shows the radial magnetic
field in nT, in a coordinate frame rotating with the Sun.

pic5: A 3-panel display showing solar wind variables as a function of radius,
22 Rsun <= r <= 220 Rsun. The plots are:

    Upper left: Solar wind number density (cm**-3) 
    Upper right: Solar wind speed (km/s) 
    Lower left: Solar wind radial momentum flux (km**2/s**2/cm**3)

pic6: A 4-panel display showing components of the solar wind magnetic field
in the solar equatorial plane (z=0), for -200 Rsun <= X, Y <= +200 Rsun.

    Upper left: Radial component of magnetic field (nT)
    Upper right: x-component of magnetic field (nT)
    Lower left: y-component of magnetic field (nT)
    Lower right: z-component of magnetic field (nT)

pic7: A 4-panel display showing pcolormesh plots in a j-slice. A j-slice is
a slice through the gamhelio data cube at a fixed colatitude. j = 0 corresponds
to the YZ plane of the gamhelio frame used in the simulation. The j = Nj/2-1
slice corresponds to the equatorial plane. The plots are:

    Upper left: Solar wind speed (km/s)
    Upper right: Solar wind number density scaled by (r/r0)**2 (cm**-3)
    Lower left: Solar wind temperature scaled by r/r0 (MK)
    Lower right: Solar wind radial magnetic field scaled by r/r0 (nT)

Authors
-------
Elena Provornikova (elena.provornikova@jhuapl.edu)
Andrew McCubbin (andrew.mccubbin@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import os
import time

# Import supplemental modules.
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
import spacepy.datamodel as dm
from sunpy.coordinates import frames

# Import project-specific modules.
from kaipy import cdaweb_utils
import kaipy.gamhelio.helioViz as hviz
import kaipy.gamhelio.heliosphere as hsph
import kaipy.kaiH5 as kh5
import kaipy.kaiTools as ktools
import kaipy.kaiViz as kv
from kaipy.satcomp import scutils


# Program constants and defaults

# Program description.
DESCRIPTION = "Make a quicklook plot for a gamhelio run."

# Default identifier for results to read.
DEFAULT_RUNID = "wsa"

# List of steps
DEFAULT_STEPS = "1"

# Default slices
DEFAULT_SLICE = None

# Code for default picture type.
DEFAULT_PICTYPE = "pic1"

# Colors to use for spacecraft position symbols.
SPACECRAFT_COLORS = list(mpl.colors.TABLEAU_COLORS.keys())


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.

    Raises
    ------
    None
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--debug", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--directory", "-d", type=str, metavar="directory",
        default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "--hgsplot", action="store_true",
        help="Plot in the Heliographic Stonyhurst frame corresponding to the "
             "date of the plot (default: %(default)s)."
    )
    parser.add_argument(
        "-id", type=str, metavar="runid", default=DEFAULT_RUNID,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-jslice", type=int, metavar="jSlice", default=None,
        help="Index of j-slice for pic7 (default: Nj/2-1)"
    )
    parser.add_argument(
        "-lon", type=float, metavar="lon", default=0.0,
        help="Longitude of meridian slice (pic2) (default: %(default)s)"
    )
    parser.add_argument(
        "--nlist", type=lambda n: [int(item) for item in n.split(',')],
        metavar="list of steps", default=DEFAULT_STEPS,
        help="List of time slice(s) n1,n2,... to plot (default: %(default)s)"
    )
    parser.add_argument(
        "--nslice", type=lambda n: [int(item) for item in n.split(':')],
        metavar="step slice", default=DEFAULT_SLICE,
        help="Slice for range of time slice(s) n1:n2 to plot "
             "(default: %(default)s)"
    )
    parser.add_argument(
        "--nworkers", "-nw", type=int, metavar="nworkers", default=4,
        help="Number of parallel workers (default: %(default)s)"
    )
    parser.add_argument(
        "--parallel", "-p", action="store_true",
        help="Read from HDF5 in parallel (default: %(default)s)."
    )
    parser.add_argument(
        "-pic", type=str, metavar="pictype",
        default=DEFAULT_PICTYPE,
        help="Code for plot type (pic1-pic7) (default: %(default)s)"
    )
    parser.add_argument(
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Names of spacecraft to plot positions, separated by commas "
             "(default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        "--inner", action="store_true",
        help="Plot inner i-slice for pic3 (default: %(default)s)."
    )

    # Return the parser.
    return parser


def initFig(pic):
    """Create the matplotlib figure for a plot.

    Determine figure size (width, height) (inches) based on the pic type.

    Parameters
    ----------
    pic : str
        String representing picture type.

    Returns
    -------
    fig : mpl.Figure
        Matplotlib Figure to use for plots.

    Raises
    ------
    None
    """
    # Figure dimensions are in inches.
    if pic == "pic1" or pic == "pic2" or pic == "pic7":
        figSz = (10, 12.5)
    elif pic == "pic3":
        figSz = (10, 6.5)
    elif pic == "pic4":
        figSz = (10, 6)
    elif pic == "pic5":
        figSz = (12, 12)
    elif pic == "pic6":
        figSz = (10, 12.5)

    # Create the figure.
    fig = plt.figure(figsize=figSz, layout="constrained")
    return fig


def fOut(runid, pic, nStp, hgsplot):
    """Compute the name of the output file.

    Compute the name of the output file.

    Parameters
    ----------
    runid : str
        ID string for run
    pic : str
        String representing picture type.
    nStp : int
        Simulation step number used in plot.
    hgsplot : bool
        True if plot is in HGS frame at the date of the plot

    Returns
    -------
    s : str
        Name of file to receive the plot.

    Raises
    ------
    None
    """
    if hgsplot:
        s = f"qkpic_{runid}_{pic}_HGS_n{nStp}.png"
    else:
        s = f"qkpic_{runid}_{pic}_n{nStp}.png"
    return s


def GHtoHGS(mjd_gh, x_gh, y_gh, z_gh, mjd_hgs):
    """Convert Cartesin GH coordinates to HGS.

    Convert Cartesian coordinates in the gamhelio frame at time mjdc to
    the Heliographic Sonyhurst frame at time mjd.

    NOTE: The gamhelio frame at time t is related to the Heliographic
    Stonyhurst frame at time t by the reflection of the x- and y-axes:

    x_gh(t) = -x_hgs(t)
    y_gh(t) = -y_hgs(t)
    z_gh(t) = z_hgs(t)

    Since HGS is a time-dependent frame, a time must be provided for each set
    of coordinates.

    Parameters
    ----------
    mjd_gh : float
        MJD of source gamhelio frame
    x_gh, y_gh, z_gh : np.array of float (any shape) or scalar float
        Cartesian coordinates in GH(mjdc) frame. All three arrays x, y, z must
        have identical shapes.
    mjd_hgs : float
        MJD of target HGS frame

    Returns
    -------
    x_hgs, y_hgs, z_hgs : np.array of float (same shape as x_gh, y_gh, z_gh)
        Cartesian coordinates converted to HGS(mjd) frame.

    Raises
    ------
    None
    """
    # Load the source coordinates (originially in the GH(mjd_gh) frame) into
    #  the equivalent HGS(mjd_gh) frame.
    c_gh = SkyCoord(
        -x_gh*u.Rsun, -y_gh*u.Rsun, z_gh*u.Rsun,
        frame=frames.HeliographicStonyhurst,
        obstime=ktools.MJD2UT(mjd_gh),
        representation_type="cartesian"
    )

    # Create the target Heliographic Stonyhurst frame.
    hgs_frame = frames.HeliographicStonyhurst(
        obstime=ktools.MJD2UT(mjd_hgs)
    )

    # Convert the coordinates from GH(mjd_gh) to HGS(mjd_hgs).
    c_hgs = c_gh.transform_to(hgs_frame)

    # Extract and return the converted coordinates.
    x_hgs = dm.dmarray(c_hgs.cartesian.x)
    y_hgs = dm.dmarray(c_hgs.cartesian.y)
    z_hgs = dm.dmarray(c_hgs.cartesian.z)
    return x_hgs, y_hgs, z_hgs


def main():
    """Make a quick-look plot for a gamhelio run."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    debug = args.debug
    fdir = args.directory
    hgsplot = args.hgsplot
    ftag = args.id
    jslice = args.jslice
    pic2lon = args.lon
    steps = args.nlist
    slices = args.nslice
    nWorkers = args.nworkers
    doParallel = args.parallel
    pic = args.pic
    spacecraft = args.spacecraft
    verbose = args.verbose
    inner = args.inner

    if slices:
        print(f"Slice selected {slice(slices[0], slices[1], slices[2])}")

    # Fetch the plot domain based on the picture type.
    tic = time.perf_counter()
    xyBds = hviz.GetSizeBds(pic)
    toc = time.perf_counter()
    print(xyBds)
    print(f"Get bounds took {toc - tic} s")

    # Do work?
    doFast = False

    # Open a pipe to the results data.
    tic = time.perf_counter()
    gsph = hsph.GamsphPipe(fdir, ftag, doFast=doFast, doParallel=doParallel,
                           nWorkers=nWorkers)
    toc = time.perf_counter()
    print(f"Open pipe took {toc-tic} s")

    # Compute the range of time steps to use.
    if slices and steps[0] == 1:
        steps = range(gsph.s0,gsph.sFin + 1)[slice(slices[0], slices[1], slices[2])]
    print(f"steps = {steps}")

    # Get the MJDc value for use in computing the gamhelio frame.
    fname = gsph.f0
    MJDc = scutils.read_MJDc(fname)

    # Split the list into individual spacecraft names.
    if spacecraft:
        spacecraft = spacecraft.split(',')

    # Create figures in a memory buffer.
    mpl.use("Agg")

    # Make a plot for each time step in the list of time steps.
    for nStp in steps:
        if debug:
            print(f"nStp = {nStp}")

        tic = time.perf_counter()
        print(f"Generating {pic} for time step {nStp}.")
        fig = initFig(pic)

        # Extract the MJD for the step.
        if any(gsph.MJDs):
            mjd = gsph.MJDs[nStp-gsph.s0]
            time_stamp = ktools.MJD2UT(mjd)
        else:
            mjd = gsph.T[nStp-gsph.s0]/(3600./gsph.tScl)
            time_stamp = f"{mjd:0.2f} [hr]"
        if debug:
            print(f"mjd = {mjd}")

        # Lay out the subplots.
        if pic in ["pic1", "pic2", "pic3", "pic6", "pic7"]:
            gs = gridspec.GridSpec(4, 6, height_ratios=[20, 1, 20, 1], figure=fig)
            # Axes for plots.
            AxL0 = fig.add_subplot(gs[0, 0:3])
            AxR0 = fig.add_subplot(gs[0, 3:])
            AxL1 = fig.add_subplot(gs[2, 0:3])
            AxR1 = fig.add_subplot(gs[2, 3:])
            # Axes for colorbars.
            AxC1_0 = fig.add_subplot(gs[1, 0:3])
            AxC2_0 = fig.add_subplot(gs[1, 3:])
            AxC1_1 = fig.add_subplot(gs[3, 0:3])
            AxC2_1 = fig.add_subplot(gs[3, 3:])
        elif pic == "pic4":
            gs = gridspec.GridSpec(2, 1, height_ratios=[20, 1], figure=fig)
            Ax = fig.add_subplot(gs[0, 0])
            AxC = fig.add_subplot(gs[1, 0])
        elif pic == "pic5":
            gs = gridspec.GridSpec(2, 2, figure=fig)
            Ax = fig.add_subplot(gs[0, 0])
            AxC = fig.add_subplot(gs[0, 1])
            AxC1 = fig.add_subplot(gs[1, 0])
        else:
            raise TypeError(f"Invalid figure type: {pic}!")

        # If the step is -1, use the last step.
        if nStp < 0:
            nStp = gsph.sFin
            print(f"Using Step {nStp}")

        # Now create the actual plots.
        if pic == "pic1":
            # Equatorial plots in the XY plane of the modified HGS frame used
            # by gamhelio. If hgsplot is True, then the plot frame is the true
            # HGS frame at the time of the plot.
            hviz.PlotEqMagV(gsph, nStp, xyBds, AxL0, AxC1_0,
                            hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotEqD(gsph, nStp, xyBds, AxR0, AxC2_0,
                         hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotEqTemp(gsph, nStp, xyBds, AxL1, AxC1_1,
                            hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotEqBr(gsph, nStp, xyBds, AxR1, AxC2_1,
                          hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            if hgsplot:
                fig.suptitle("Heliographic Stonyhurst frame for "
                             f"{time_stamp}")
            else:
                fig.suptitle(f"GAMERA-Helio frame for {time_stamp}")
        elif pic == "pic2":
            # Meridional plots in the XZ plane of the modified HGS frame used
            # by gamhelio. If hgsplot is True, then the plot frame is the true
            # HGS frame at the time of the plot.
            hviz.PlotMerMagV(gsph, nStp, xyBds, AxL0, AxC1_0,
                             indx=(None, pic2lon),
                             hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotMerDNorm(gsph, nStp, xyBds, AxR0, AxC2_0,
                              indx=(None, pic2lon),
                              hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotMerTemp(gsph, nStp, xyBds, AxL1, AxC1_1,
                             indx=(None, pic2lon),
                             hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotMerBrNorm(gsph, nStp, xyBds, AxR1, AxC2_1,
                               indx=(None, pic2lon),
                               hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            if hgsplot:
                fig.suptitle("Heliographic Stonyhurst frame for "
                             f"{time_stamp}")
            else:
                fig.suptitle(f"GAMERA-Helio frame for {time_stamp}")
        elif pic == "pic3":
            # Lat/lon plot at 1 AU (the outer edge of the gamhelio grid), in
            # the modified HGS frame rotating with the Sun.
            AU_RSUN = 215.0
            radius = AU_RSUN
            I_RSUN = 21.5
            if inner : radius = I_RSUN
            hviz.PlotiSlMagV(gsph, nStp, xyBds, AxL0, AxC1_0, idx=radius,
                             idx_is_radius=True,
                             hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotiSlD(gsph, nStp, xyBds, AxR0, AxC2_0, idx=radius,
                          idx_is_radius=True,
                          hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd,
                          use_outer_range=(not inner) )
            hviz.PlotiSlTemp(gsph, nStp, xyBds, AxL1, AxC1_1, idx=radius,
                             idx_is_radius=True,
                             hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd,
                             use_outer_range=(not inner))
            hviz.PlotiSlBr(gsph, nStp, xyBds, AxR1, AxC2_1, idx=radius,
                           idx_is_radius=True,
                           hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd,
                           use_outer_range=(not inner))
            if hgsplot:
                fig.suptitle(f"Heliographic Stonyhurst frame at {radius} [RS] for {time_stamp}")
            else:
                fig.suptitle(f"GAMERA-Helio frame at {radius} [RS] for {time_stamp}")
        elif pic == "pic4":
            # Plot at 1 AU in frame rotating with Sun.
            hviz.PlotiSlBrRotatingFrame(gsph, nStp, xyBds, Ax, AxC)
        elif pic == "pic5":
            hviz.PlotDensityProf(gsph, nStp, xyBds, Ax)
            hviz.PlotSpeedProf(gsph, nStp, xyBds, AxC)
            hviz.PlotFluxProf(gsph, nStp, xyBds, AxC1)
        elif pic == "pic6":
            hviz.PlotEqBr(gsph, nStp, xyBds, AxL0, AxC1_0, hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotEqBx(gsph, nStp, xyBds, AxR0, AxC2_0, hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotEqBy(gsph, nStp, xyBds, AxL1, AxC1_1, hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotEqBz(gsph, nStp, xyBds, AxR1, AxC2_1, hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            fig.suptitle("GAMERA-Helio frame for "
                         f"{time_stamp}")
        elif pic == "pic7":
            if jslice is None:
                jidx = gsph.Nj//2 - 1
            else:
                jidx = jslice
            hviz.PlotjMagV(gsph, nStp, xyBds, AxL0, AxC1_0, jidx=jidx, hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotjD(gsph, nStp, xyBds, AxR0, AxC2_0, jidx=jidx, hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotjTemp(gsph, nStp, xyBds, AxL1, AxC1_1, jidx=jidx, hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            hviz.PlotjBr(gsph, nStp, xyBds, AxR1, AxC2_1, jidx=jidx, hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
            fig.suptitle("GAMERA-Helio frame for "
                         f"{time_stamp}")
        else:
            raise TypeError(f"Invalid figure type: {pic}!")

        # Add time in the upper left (if not in figure title).
        # if pic == "pic1" or pic == "pic2" or pic == "pic3" or pic == "pic6":
        if pic == "pic4" or pic == "pic5":
            gsph.AddTime(nStp, Ax, xy=[0.015, 0.92], fs="small")

        # Overlay the spacecraft positions.
        if spacecraft:

            # Fetch the MJD at start and end of the model results.
            MJD_start = kh5.tStep(fname, gsph.s0, aID="MJD")
            MJD_end = kh5.tStep(fname, gsph.sFin, aID="MJD")

            # Convert the start and stop MJD to a datetime object in UT.
            ut_start = ktools.MJD2UT(MJD_start)
            ut_end = ktools.MJD2UT(MJD_end)

            # Fetch the trajectory of each spacecraft from CDAWeb. Then
            # interpolate the position at the time of the plot, and plot the
            # spacecraft at the interpolated position.
            for (i_sc, sc_id) in enumerate(spacecraft):
                if verbose:
                    print(f"Fetching trajectory for {sc_id}.")

                # Fetch the spacecraft trajectory in whatever frame is
                # available from CDAWeb.
                sc_data = cdaweb_utils.fetch_helio_spacecraft_trajectory(
                    sc_id, ut_start, ut_end
                )
                if sc_data is None:
                    print(f"No trajectory found for {sc_id}.")
                    continue

                # Ingest the trajectory by converting it to the GH(MJDc) frame.
                if verbose:
                    print(f"Converting ephemeris for {sc_id} to the gamhelio "
                          f"frame at MJD {MJDc}.")
                t_strings = np.array([str(t) for t in sc_data["Epoch"]])
                t = astropy.time.Time(t_strings, scale='utc').mjd
                x, y, z = cdaweb_utils.ingest_helio_spacecraft_trajectory(
                    sc_id, sc_data, MJDc)

                # Interpolate the spacecraft position at the time for the plot.
                t_sc = mjd
                x_sc = np.interp(t_sc, t, x)
                y_sc = np.interp(t_sc, t, y)
                z_sc = np.interp(t_sc, t, z)

                # If needed, convert the position to HGS(mjd).
                if hgsplot:
                    x_sc, y_sc, z_sc = GHtoHGS(MJDc, x_sc, y_sc, z_sc, mjd)

                # If needed, compute heliocentric spherical coordinates
                # for the interpolated spacecraft position. Longitude is in
                # the -180 to +180 range. Convert to 0-360 if not using
                # hgsplot.
                if pic == "pic3" or pic == "pic4":
                    rxy = np.sqrt(x_sc**2 + y_sc**2)
                    theta = np.arctan2(rxy, z_sc)
                    phi = np.arctan2(y_sc, x_sc)
                    lat = np.degrees(np.pi/2 - theta)
                    lon = np.degrees(phi)
                    lat_sc = lat
                    lon_sc = lon
                    if not hgsplot:
                        if lon_sc < 0:
                            lon_sc += 360

                # Plot the position of the spacecraft at the plot time. Each
                # spacecraft is plotted as a colored dot with a black outline.
                color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
                x_nudge = 0.0
                y_nudge = 8.0
                if pic == "pic1":
                    for ax in (AxL0, AxR0, AxL1, AxR1):
                        ax.plot(x_sc, y_sc, 'o', c=color)
                        ax.plot(x_sc, y_sc, 'o', c="black", fillstyle="none")
                        ax.text(x_sc + x_nudge, y_sc + y_nudge, sc_id,
                                c="black", horizontalalignment="center")
                elif pic == "pic2":
                    for ax in (AxL0, AxR0, AxL1, AxR1):
                        ax.plot(x_sc, z_sc, 'o', c=color)
                        ax.plot(x_sc, z_sc, 'o', c="black", fillstyle="none")
                        ax.text(x_sc + x_nudge, z_sc + y_nudge, sc_id,
                                c="black", horizontalalignment="center")
                elif pic == "pic3":
                    for ax in (AxL0, AxR0, AxL1, AxR1):
                        ax.plot(lon_sc, lat_sc, 'o', c=color)
                        ax.plot(lon_sc, lat_sc, 'o', c="black",
                                fillstyle="none")
                        ax.text(lon_sc + x_nudge, lat_sc + y_nudge, sc_id,
                                c="black", horizontalalignment="center")
                elif pic == "pic4":
                    ax = Ax
                    ax.plot(lon_sc, lat_sc, 'o', c=color)
                    ax.plot(lon_sc, lat_sc, 'o', c="black", fillstyle="none")
                    ax.text(lon_sc + x_nudge, lat_sc + y_nudge, sc_id,
                            c="black", horizontalalignment="center")
                elif pic == "pic5":
                    pass
                elif pic == "pic6":
                    for ax in (AxL0, AxR0, AxL1, AxR1):
                        ax.plot(x_sc, y_sc, 'o', c=color)
                        ax.plot(x_sc, y_sc, 'o', c="black", fillstyle="none")
                        ax.text(x_sc + x_nudge, y_sc + y_nudge, sc_id,
                                c="black", horizontalalignment="center")
                elif pic == "pic7":
                    raise TypeError("Spacecraft not supported for pic7!")
                else:
                    raise TypeError(f"Invalid plot code: {pic}!")

        # Save the figure to a file.
        path = os.path.join(fdir, fOut(ftag, pic, nStp, hgsplot))
        kv.savePic(path, bLenX=40)
        plt.close()
        toc = time.perf_counter()
        print(f"Step {nStp} took {toc-tic} s")


if __name__ == "__main__":
    main()
