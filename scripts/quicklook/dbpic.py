#!/usr/bin/env python


"""Plot the ground magnetic field perturbations from a magnetosphere run.

Plot the ground magnetic field perturbations from a magnetosphere run.

Author
------
Kareem Sorathia (kareem.sorathia@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import datetime
import os

# Import 3rd-party modules.
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from psutil import cpu_count

# Import project-specific modules.
import kaipy.cdaweb_utils as cdaweb_utils
import kaipy.cmaps.kaimaps as kmaps
import kaipy.gamera.deltabViz as dbViz
import kaipy.gamera.gampp as gampp
import kaipy.kaiH5 as kh5
import kaipy.kaiTools as ktools
import kaipy.kaiViz as kv


# Program constants and defaults

# Program description.
description = "Plot the ground magnetic field perturbations for a MAGE magnetosphere run."

# Default identifier for results to read.
default_runid = "msphere"

# Plot the last step by default.
default_step = -1

# Default vertical layer to plot.
default_k0 = 0

# Do not show anomalous currents by default.
default_Jr = False

# Default output filename.
default_output_filename = "qkdbpic.png"

# Size of figure in inches (width x height).
figSz = (12, 6)

# Color to use for magnetic footprint positions.
FOOTPRINT_COLOR = 'red'


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
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-d", type=str, metavar="directory", default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "-n", type=int, metavar="step", default=default_step,
        help="Time slice to plot (default: %(default)s)"
    )
    parser.add_argument(
        "-id", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-Jr", action="store_true", default=default_Jr,
        help="Show radial component of anomalous current (default: %(default)s)."
    )
    parser.add_argument(
        '-k0', type=int, metavar="layer", default=default_k0,
        help="Vertical layer to plot (default: %(default)s)")
    parser.add_argument(
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Names of spacecraft to plot magnetic footprints, separated by commas (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        '-vid', action='store_true', default=False,
        help="Make a video and store in mixVid directory (default: %(default)s)"
    )
    parser.add_argument(
        '-overwrite', action='store_true', default=False,
        help="Overwrite existing vid files (default: %(default)s)"
    )
    parser.add_argument(
        '--ncpus', type=int, metavar="ncpus", default=1,
        help="Number of threads to use with --vid (default: %(default)s)"
    )
    return parser

def makePlot(i, fname, dbdata, k0, Re, spacecraft, nStp):

    plt.clf()

    if do_vid:
        filename = "{}.{:0>{n}d}.png".format("dbpic", i, n=n_pad)
        outPath = os.path.join(outDir, filename)
        if not do_overwrite and os.path.exists(outPath):
            return
    else:
        outPath = default_output_filename

    # Check the vertical level.
    Z0 = dbViz.CheckLevel(dbdata, k0, Re)
    if debug:
        print("Z0 = %s" % Z0)

    # If currents were requested, read them. Otherwise, read the ground
    # magnetic field perturbations.
    if (doJr):
        print("Reading Jr ...")
        Jr = dbdata.GetVar("dbJ", nStp, doVerb=False)[:, :, k0]
        Q = Jr
    else:
        dBn = dbdata.GetVar("dBn", nStp, doVerb=True)[:, :, k0]
        Q = dBn

    # Convert MJD to UT.
    MJD = kh5.tStep(fname, nStp, aID="MJD")
    if debug:
        print("MJD = %s" % MJD)
    utS = ktools.MJD2UT([MJD])
    if debug:
        print("utS = %s" % utS)
    utDT= utS[0]
    if debug:
        print("utDT = %s" % utDT)

    # Create the mapping grid.
    crs = ccrs.PlateCarree()
    if debug:
        print("ccrs = %s" % ccrs)
    LatI, LonI, LatC, LonC = dbViz.GenUniformLL(dbdata, k0)
    if debug:
        print("LatI = %s" % LatI)
        print("LonI = %s" % LonI)
        print("LatC = %s" % LatC)
        print("LonC = %s" % LonC)

    # Fetch the color map.
    cmap = kmaps.cmDiv
    if debug:
        print("cmap = %s" % cmap)

    # Determine color bar settings.
    if (doJr):
        vQ = kv.genNorm(jMag)
        cbStr = "Anomalous current"
    else:
        vQ = kv.genNorm(bMag, doSymLog=True, linP=bLin)
        cbStr = r"$\Delta B_N$ [nT]"
    if debug:
        print("vQ = %s" % vQ)
        print("cbStr = %s" % cbStr)

    # Specify the grid for the subplots.
    gs = gridspec.GridSpec(3, 1, height_ratios=[20, 1.0, 1.0], hspace=0.025)

    # Create the subplots.
    AxM = fig.add_subplot(gs[0, 0], projection=crs)
    AxCB = fig.add_subplot(gs[-1, 0])

    # Make the plot.
    AxM.pcolormesh(LonI, LatI, Q, norm=vQ, cmap=cmap)

    # If requested, overlay the spacecraft magnetic footprints.
    if spacecraft:
        print("Overplotting magnetic footprints of %s." % spacecraft)

        # Split the list into individual spacecraft names.
        spacecraft = spacecraft.split(',')

        # Fetch the position of each footprint pair from CDAWeb.
        for sc in spacecraft:

            # Fetch the northern footprint position.
            fp_nlat, fp_nlon = cdaweb_utils.fetch_satellite_magnetic_northern_footprint_position(
                sc, utS[0]
            )
            if debug:
                print("fp_nlat, fp_nlon = %s, %s" % (fp_nlat, fp_nlon))

            # Fetch the southern footprint position.
            fp_slat, fp_slon = cdaweb_utils.fetch_satellite_magnetic_southern_footprint_position(
                sc, utS[0]
            )
            if debug:
                print("fp_slat, fp_slon = %s, %s" % (fp_slat, fp_slon))

            # Plot a labelled dot at the location of each footprint.
            # Skip if no footprint position found.
            if fp_nlon is not None:
                AxM.plot(fp_nlon, fp_nlat, 'o', c=FOOTPRINT_COLOR)
                lon_nudge = 2.0
                lat_nudge = 2.0
                AxM.text(fp_nlon + lon_nudge, fp_nlat + lat_nudge, sc + ' (N)')
            else:
                print("No northern footprint found for spacecraft %s." % sc)
            if fp_slon is not None:
                AxM.plot(fp_slon, fp_slat, 'o', c=FOOTPRINT_COLOR)
                lon_nudge = 2.0
                lat_nudge = 2.0
                AxM.text(fp_slon + lon_nudge, fp_slat + lat_nudge, sc + ' (S)')
            else:
                print("No southern footprint found for spacecraft %s." % sc)

    # Make the colorbar.
    kv.genCB(AxCB, vQ, cbStr, cM=cmap)

    # Add labels and other decorations.
    tStr = dbViz.GenTStr(AxM, fname, nStp)
    if debug:
        print("tStr = %s" % tStr)
    dbViz.DecorateDBAxis(AxM, crs, utDT)

    # Save the figure.
    if debug:
        print("fOut = %s" % outPath)
    kv.savePic(outPath)

if __name__ == "__main__":
    """Plot the ground magnetic field perturbations."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    fdir = args.d
    runid = args.id
    nStp = args.n
    k0 = args.k0
    doJr = args.Jr
    spacecraft = args.spacecraft
    verbose = args.verbose
    do_vid = args.vid
    do_overwrite = args.overwrite
    ncpus = args.ncpus

    if debug:
        print("args = %s" % args)

    # Fetch constants.
    bLin = dbViz.dbLin
    bMag = dbViz.dbMag
    jMag = dbViz.jMag

    # Compute the name of the file containing the ground magnetic field perturbations.
    ftag = runid + ".deltab"
    if debug:
        print("ftag = %s" % ftag)

    # Read the ground magnetic field perturbations.
    fname = os.path.join(fdir, ftag + ".h5")
    if debug:
        print("fname = %s" % fname)
    dbdata = gampp.GameraPipe(fdir, ftag)
    if debug:
        print("dbdata = %s" % dbdata)
    print("---")

    # Get the ID of the coordinate system, and the Earth radius.
    CoordID, Re = dbViz.GetCoords(fname)
    print("Found %s coordinate data ..." % CoordID)
    if debug:
        print("CoordID = %s" % CoordID)
        print("Re = %s" % Re)

    # Create plot in memory.
    mpl.use("Agg")

    # Set global plot font options.
    mpl.rc('mathtext', fontset='stixsans', default='regular')
    mpl.rc('font', size=10)

    # Create the figure to hold the plot.
    fig = plt.figure(figsize=figSz)

    if not do_vid: # If we are making a single image, keep original functionality
        # If the last simulation step was requested, get the step number.
        if nStp < 0:
            nStp = dbdata.sFin
            print("Using Step %d" % nStp)
        makePlot(nStp, fname, dbdata, k0, Re, spacecraft, nStp)
    else:  # then we make a video, i.e. series of images saved to msphVid

        # Get video loop parameters
        s0 = max(dbdata.s0,1) # Skip Step#0
        sFin = dbdata.sFin
        nsteps = sFin - s0
        sIds = np.array(range(s0,sFin))
        outDir = 'dbVid'
        kh5.CheckDirOrMake(outDir)

        # How many 0's do we need for filenames?
        n_pad = int(np.log10(nsteps)) + 1

        if ncpus == 1:
            for i, nStp in enumerate(sIds):
                makePlot(i, spacecraft, nStp)
        else:
            # Make list of parallel arguments
            ag = ((i,fname, dbdata, k0, Re, spacecraft,nStp) for i, nStp in enumerate(sIds) )
            # Check we're not exceeding cpu_count on computer
            ncpus = min(int(ncpus),cpu_count(logical=False))
            print('Doing multithreading on ',ncpus,' threads')
            # Do parallel job
            with Pool(processes=ncpus) as pl:
                pl.starmap(makePlot,ag)
            print("Done making all the images. Go to mixVid folder")
