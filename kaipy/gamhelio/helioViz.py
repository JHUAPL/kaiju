#Various helper routines for heliosphere quick look plots

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
import numpy as np

import kaipy.kaiViz as kv
import kaipy.gamhelio.heliosphere as hsph
from kaipy.kdefs import *
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunpy.coordinates import frames
import kaipy.kaiTools as ktools
import spacepy.datamodel as dm


VMax = 800.
VMin = 300.
MagVCM = "inferno"

DMax = 1200.
DMin = 300.
DCM = "copper_r"

# limits for iSlice
# 21.5 R_S
D0Max = 2000.
D0Min = 300.
# 220 R_S
D0Max_outer = 1.
D0Min_outer = 20.
D0CM = "copper_r"

TMin = 0.2
TMax = 2.
TCM = "copper"

# 21.5 R_S
T0Min = 0.2
T0Max = 2.0
# 220 R_S
T0Min_outer = 0.0
T0Max_outer = 0.3

# 21.5 R_S
BMax = 150.
BMin = -150.
# 220 R_S
BCM = "coolwarm"

# 21.5 R_S
B0Min = -150.
B0Max = 150.
# 220 R_S
B0Min_outer = -5.0
B0Max_outer = 5.0

BZMin = -50
BZMax = 50


colorProf = "tab:orange"
#Function to Add different size options to argument
#not used for helio right now
def AddSizeArgs(parser):
    parser.add_argument('-small' , action='store_true', default=False,help="Use smaller domain bounds (default: %(default)s)")
    parser.add_argument('-big'   , action='store_true', default=False,help="Use larger domain bounds (default: %(default)s)")
    parser.add_argument('-bigger', action='store_true', default=False,help="Use larger-er domain bounds (default: %(default)s)")
    parser.add_argument('-huge'  , action='store_true', default=False,help="Use huge domain bounds (default: %(default)s)")

#Return domain size from parsed arguments; see msphViz for options
def GetSizeBds(pic):
    if (pic == "pic1" or pic == "pic2"):
                #for inner helio
        xyBds = [-220.,220.,-220.,220.]
                #for 1-10 au helio
                #xyBds = [-10.,10.,-10.,10.]
    elif (pic == "pic3"):
        xyBds = [0.,360.,-75.,75.]
    elif (pic == "pic4"):
        xyBds = [0.,360.,-90.,90.]
    elif (pic == "pic5"):
        xyBds = [20.,220.,1.,2000.]
    elif (pic == "pic6" or pic == "pic7"):
        xyBds = [-220.,220.,-220.,220.]
    else:        
        raise RuntimeError("No compatible pic type specified.")
    return xyBds


def PlotEqMagV(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        hgsplot=False, MJDc=None, MJD_plot=None
):
    """Plot solar wind speed in the solar equatorial plane.

    Plot solar wind speed in the solar equatorial plane. By default, the plot
    is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If hgsplot is True, MJDc and MJD_plot must be
    specified. In that case, the coordinates are mapped from the GH(MJDc) frame
    to the HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that, at any particular MJD:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc (the MJD of the central meridian in the WSA
    file used for initial conditions), meaning it is fixed in spatial
    orientation at that time. The HGS frame is defined at MJD_plot, also
    producing a (different) fixed spatial orientation. The conversion maps
    points in the GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a
    rotation about the z-axis, but also accounting for the eccentricity of
    Earth's orbit and other astronomical parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting
    doDeco : bool
        If True, add axis labels and other decorations to the plot
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot

    Returns
    -------
    MagV : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vMagV, cbT=None, cM=MagVCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    MagV = gsph.eqMagV(nStp)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originally in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame. All z values are
        # set to 0 since we are plotting in the equatorial (XY) plane.
        zzi = np.zeros_like(gsph.xxi)
        c = SkyCoord(
            -gsph.xxi*u.Rsun, -gsph.yyi*u.Rsun, zzi*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(x, y, MagV, cmap=MagVCM, norm=vMagV)

    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, MagV, cmap=MagVCM, norm=vMagV)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Speed [$km/s$]")
        Ax.set_xlabel(r"$X$ [$R_S$]")
        Ax.set_ylabel(r"$Y$ [$R_S$]")
        Ax.yaxis.tick_left()
        Ax.yaxis.set_label_position('left')
    # Return the data.
    return MagV


def PlotjMagV(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True, jidx=-1,
        hgsplot=False, MJDc=None, MJD_plot=None
):
    """Plot solar wind speed in a specific j plane.

    Plot solar wind speed in the a specific j plane. By default, the plot
    is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If MJD_plot is specified, MJDc must also be specified.
    In that case, the coordinates are mapped from the GH(MJDc) frame to the
    HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc, meaning it is fixed in spatial orientation
    at that time. The HGS frame is defined at MJD_plot, also producing a
    (different) fixed spatial orientation. The conversion maps points in the
    GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a rotation about
    the z-axis, but also accounting for the Earth's orbit and other
    astronommical and geodetic parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If true, clear the plot Axes before further plotting.
    doDeco : bool
        If true, add axis labels to the plot.
    jidx : int
        Index of j-plane to plot.
    hgsplot : bool
        If true, plot in HGS(MJD_plot) frame.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.

    Returns
    -------
    MagV : np.array of float
        Data for solar wind speed in selected j-plane, same shape as the
        j-plane in the gamhelio results.

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vMagV, cbT=None, cM=MagVCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    MagV = gsph.jMagV(nStp, jidx=jidx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:
        raise TypeError("HGS frame not supported for pic7!")
    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, MagV, cmap=MagVCM, norm=vMagV)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Speed [$km/s$]")
        Ax.set_xlabel(r"$X [R_S]$")
        Ax.set_ylabel(r"$Y [R_S]$")

    # Return the data.
    return MagV


def PlotMerMagV(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        indx=(None, None),
        hgsplot=False, MJDc=None, MJD_plot=None
):
    """Plot solar wind speed in a meridional plane.

    Plot solar wind speed in a meridional plane. By default, the plot
    is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If hgsplot is True, MJDc and MJD_plot must be
    specified. In that case, the coordinates are mapped from the GH(MJDc) frame
    to the HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that, at any particular MJD:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc (the MJD of the central meridian in the WSA
    file used for initial conditions), meaning it is fixed in spatial
    orientation at that time. The HGS frame is defined at MJD_plot, also
    producing a (different) fixed spatial orientation. The conversion maps
    points in the GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a
    rotation about the z-axis, but also accounting for the eccentricity of
    Earth's orbit and other astronomical parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting
    doDeco : bool
        If True, add axis labels and other decorations to the plot
    indx : tuple of 2 int or float
        Index or angle of meridional slice to plot
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot

    Returns
    -------
    Vr, Vl : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vMagV, cbT=None, cM=MagVCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Determine the angle of the meridional slice.
    phi = 0.0
    for idx in indx:
        if type(idx) is not None:
            if type(idx) is int:
                phi = idx/gsph.Nk*2*np.pi
            elif type(idx) is float:
                phi = idx
        else: 
            phi = ""

    # Fetch the data.
    # r(ight) and l(eft) data arrays
    Vr, Vl = gsph.MerMagV(nStp, indx=indx)

    # Fetch the coordinates of the grid cell corners in the specified
    # meridional plane, in the GH(MJDc) frame.
    xr, yr, zr, xl, yl, zl, r = gsph.MeridGridHalfs(*indx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the meridional grid cell vertex coordinates (originially in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame.
        cr = SkyCoord(
            -xr*u.Rsun, -yr*u.Rsun, zr*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )
        cl = SkyCoord(
            -xl*u.Rsun, -yl*u.Rsun, zl*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        cr = cr.transform_to(hgs_frame)
        cl = cl.transform_to(hgs_frame)

        # Extract the converted coordinates.
        xr = dm.dmarray(cr.cartesian.x)
        yr = dm.dmarray(cr.cartesian.y)
        zr = dm.dmarray(cr.cartesian.z)
        xl = dm.dmarray(cl.cartesian.x)
        yl = dm.dmarray(cl.cartesian.y)
        zl = dm.dmarray(cl.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(np.sqrt(xr**2 + yr**2), zr, Vl, cmap=MagVCM, norm=vMagV)
        Ax.pcolormesh(-np.sqrt(xl**2 + yl**2), zl, Vr, cmap=MagVCM, norm=vMagV)

    else:
        Ax.pcolormesh(np.sqrt(xr**2 + yr**2), zr, Vr, cmap=MagVCM, norm=vMagV)
        Ax.pcolormesh(-np.sqrt(xl**2 + yl**2), zl, Vl, cmap=MagVCM, norm=vMagV)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Speed [$km/s$]")
        Ax.set_xlabel(r"$R_{XY}$ [$R_S$] at $\phi=" +
                      f"{phi:{2}.{2}}$ [$rad$]")
        Ax.set_ylabel("$Z$ [$R_S$]")

    # Return the data.
    return Vr, Vl


def PlotMerDNorm(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        indx=(None, None),
        hgsplot=False, MJDc=None, MJD_plot=None
):
    """Plot normalized solar wind number density in a meridional plane.

    Plot normalized solar wind number density in a meridional plane. By
    default, the plot is produced in the GH(MJDc) frame (the gamhelio frame
    used for the simulation results). If hgsplot is True, MJDc and MJD_plot
    must be specified. In that case, the coordinates are mapped from the
    GH(MJDc) frame to the HGS(MJD_plot) frame.

    The density is normalized with the factor (r/r0)**2, where r0 is 21.5 Rsun
    (the inner edge of the gamhelio grid).

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that, at any particular MJD:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc (the MJD of the central meridian in the WSA
    file used for initial conditions), meaning it is fixed in spatial
    orientation at that time. The HGS frame is defined at MJD_plot, also
    producing a (different) fixed spatial orientation. The conversion maps
    points in the GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a
    rotation about the z-axis, but also accounting for the eccentricity of
    Earth's orbit and other astronomical parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting.
    doDeco : bool
        If True, add axis labels to the plot.
    indx : tuple of 2 int or float
        Index or angle of meridional slice to plot
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.

    Returns
    -------
    Dr, Dl : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vD = kv.genNorm(DMin, DMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vD, cbT=None, cM=DCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Determine the angle of the meridional slice.
    phi = 0.0
    for idx in indx:
        if type(idx) is not None:
            if type(idx) is int:
                phi = idx/gsph.Nk*2*np.pi
            elif type(idx) is float:
                phi = idx
        else:
            phi = ""

    # Fetch the data.
    # r(ight) and l(eft) data arrays
    Dr, Dl = gsph.MerDNrm(nStp, indx=indx)

    # Fetch the coordinates of the grid cell corners in the specified
    # meridional plane, in the GH(MJDc) frame.
    xr, yr, zr, xl, yl, zl, r = gsph.MeridGridHalfs(*indx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originially in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame.
        cr = SkyCoord(
            -xr*u.Rsun, -yr*u.Rsun, zr*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )
        cl = SkyCoord(
            -xl*u.Rsun, -yl*u.Rsun, zl*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        cr = cr.transform_to(hgs_frame)
        cl = cl.transform_to(hgs_frame)

        # Extract the converted coordinates.
        xr = dm.dmarray(cr.cartesian.x)
        yr = dm.dmarray(cr.cartesian.y)
        zr = dm.dmarray(cr.cartesian.z)
        xl = dm.dmarray(cl.cartesian.x)
        yl = dm.dmarray(cl.cartesian.y)
        zl = dm.dmarray(cl.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(np.sqrt(xr**2 + yr**2), zr, Dl, cmap=DCM, norm=vD,
                      shading='auto')
        Ax.pcolormesh(-np.sqrt(xl**2 + yl**2), zl, Dr, cmap=DCM, norm=vD,
                      shading='auto')

    else:
        Ax.pcolormesh(np.sqrt(xr**2 + yr**2), zr, Dr, cmap=DCM, norm=vD,
                      shading='auto')
        Ax.pcolormesh(-np.sqrt(xl**2 + yl**2), zl, Dl, cmap=DCM, norm=vD,
                      shading='auto')

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Number density $n$ [$(r/r_0)^2 cm^{-3}$]")
        Ax.set_xlabel(r"$R_{XY}$ [$R_S$] at $\phi=" +
                      f"{phi:{2}.{2}}$ [$rad$]")
        Ax.set_ylabel("Z [$R_S$]")
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the data.
    return Dr, Dl


def PlotMerBrNorm(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        indx=(None, None),
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot normalized solar wind radial magnetic field in a meridional plane.

    Plot normalized solar wind radial magnetic field in a meridional plane. By
    default, the plot is produced in the GH(MJDc) frame (the gamhelio frame
    used for the simulation results). If hgsplot is True, MJDc and MJD_plot
    must be specified. In that case, the coordinates are mapped from the
    GH(MJDc) frame to the HGS(MJD_plot) frame.

    The temperature is normalized with the factor (r/r0)**2, where r0 is
    21.5 Rsun (the inner edge of the gamhelio grid).

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that, at any particular MJD:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc (the MJD of the central meridian in the WSA
    file used for initial conditions), meaning it is fixed in spatial
    orientation at that time. The HGS frame is defined at MJD_plot, also
    producing a (different) fixed spatial orientation. The conversion maps
    points in the GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a
    rotation about the z-axis, but also accounting for the eccentricity of
    Earth's orbit and other astronomical parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting.
    doDeco : bool
        If True, add axis labels to the plot.
    indx : tuple of 2 int or float
        Index or angle of meridional slice to plot
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.

    Returns
    -------
    Br_r, Br_l : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, cbT=None, cM=BCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Determine the angle of the meridional slice.
    phi = 0.0
    for idx in indx:
        if type(idx) is not None:
            if type(idx) is int:
                phi = idx/gsph.Nk*2*np.pi
            elif type(idx) is float:
                phi = idx
        else:
            phi = ""

    # Fetch the data.
    # r(ight) and l(eft) data arrays
    Br_r, Br_l = gsph.MerBrNrm(nStp, indx=indx)

    # Fetch the coordinates of the grid cell corners in the specified
    # meridional plane, in the GH(MJDc) frame.
    xr, yr, zr, xl, yl, zl, r = gsph.MeridGridHalfs(*indx)

    # Compute the cell center coordinates (used to plot current sheet).
    xr_c = 0.25*(xr[:-1, :-1] + xr[:-1, 1:] + xr[1:, :-1] + xr[1:, 1:])
    yr_c = 0.25*(yr[:-1, :-1] + yr[:-1, 1:] + yr[1:, :-1] + yr[1:, 1:])
    zr_c = 0.25*(zr[:-1, :-1] + zr[:-1, 1:] + zr[1:, :-1] + zr[1:, 1:])
    xl_c = 0.25*(xl[:-1, :-1] + xl[:-1, 1:] + xl[1:, :-1] + xl[1:, 1:])
    yl_c = 0.25*(yl[:-1, :-1] + yl[:-1, 1:] + yl[1:, :-1] + yl[1:, 1:])
    zl_c = 0.25*(zl[:-1, :-1] + zl[:-1, 1:] + zl[1:, :-1] + zl[1:, 1:])

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex and center coordinates
        # (originially in the GH(MJDc) frame) in the equivalent HGS(MJDc)
        # frame. Set all z values to 0 since we are using the solar equatorial
        # plane.
        cr = SkyCoord(
            -xr*u.Rsun, -yr*u.Rsun, zr*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )
        cl = SkyCoord(
            -xl*u.Rsun, -yl*u.Rsun, zl*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )
        cr_c = SkyCoord(
            -xr_c*u.Rsun, -yr_c*u.Rsun, zr_c*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )
        cl_c = SkyCoord(
            -xl_c*u.Rsun, -yl_c*u.Rsun, zl_c*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        cr = cr.transform_to(hgs_frame)
        cl = cl.transform_to(hgs_frame)
        cr_c = cr_c.transform_to(hgs_frame)
        cl_c = cl_c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        xr = dm.dmarray(cr.cartesian.x)
        yr = dm.dmarray(cr.cartesian.y)
        zr = dm.dmarray(cr.cartesian.z)
        xl = dm.dmarray(cl.cartesian.x)
        yl = dm.dmarray(cl.cartesian.y)
        zl = dm.dmarray(cl.cartesian.z)
        xr_c = dm.dmarray(cr_c.cartesian.x)
        yr_c = dm.dmarray(cr_c.cartesian.y)
        zr_c = dm.dmarray(cr_c.cartesian.z)
        xl_c = dm.dmarray(cl_c.cartesian.x)
        yl_c = dm.dmarray(cl_c.cartesian.y)
        zl_c = dm.dmarray(cl_c.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(np.sqrt(xr**2 + yr**2), zr, Br_l, cmap=BCM, norm=vB,
                      shading="auto")
        Ax.pcolormesh(-np.sqrt(xl**2 + yl**2), zl, Br_r, cmap=BCM, norm=vB,
                      shading="auto")

        # Plot the heliospheric current sheet.
        Ax.contour(np.sqrt(xr_c**2 + yr_c**2), zr_c, Br_l, [0.],
                   colors='black')
        Ax.contour(-np.sqrt(xl_c**2 + yl_c**2), zl_c, Br_r, [0.],
                   colors='black')

    else:
        Ax.pcolormesh(np.sqrt(xr**2 + yr**2), zr, Br_r, cmap=BCM, norm=vB,
                      shading='auto')
        Ax.pcolormesh(-np.sqrt(xl**2 + yl**2), zl, Br_l, cmap=BCM, norm=vB,
                      shading='auto')
        Ax.contour(np.sqrt(xr_c**2 + yr_c**2), zr_c, Br_r, [0.],
                   colors='black')
        Ax.contour(-np.sqrt(xl_c**2 + yl_c**2), zl_c, Br_l, [0.],
                   colors='black')

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r'Radial magnetic field $B_r$ [$(r/r_0)^2 nT$]')
        Ax.set_xlabel(r"$R_{XY}$ [$R_S$] at $\phi=" +
                      f"{phi:{2}.{2}}$ [$rad$]")
        Ax.set_ylabel('Z [$R_S$]')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the data.
    return Br_r, Br_l


def PlotMerTemp(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        indx=(None, None),
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot normalized solar wind temperature in a meridional plane.

    Plot normalized solar wind temperature in a meridional plane. By default,
    the plot is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If hgsplot is True, MJDc and MJD_plot must be
    specified. In that case, the coordinates are mapped from the GH(MJDc) frame
    to the HGS(MJD_plot) frame.

    The temperature is normalized with the factor (r/r0), where r0 is 21.5 Rsun
    (the inner edge of the gamhelio grid).

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that, at any particular MJD:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc (the MJD of the central meridian in the WSA
    file used for initial conditions), meaning it is fixed in spatial
    orientation at that time. The HGS frame is defined at MJD_plot, also
    producing a (different) fixed spatial orientation. The conversion maps
    points in the GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a
    rotation about the z-axis, but also accounting for the eccentricity of
    Earth's orbit and other astronomical parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting.
    doDeco : bool
        If True, add axis labels to the plot.
    indx : tuple of 2 int or float
        Index or angle of meridional slice to plot
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.

    Returns
    -------
    Tempr, Templ : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vT, cbT=None, cM=TCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Determine the angle of the meridional slice.
    phi = 0.0
    for idx in indx:
        if type(idx) is not None:
            if type(idx) is int:
                phi = idx/gsph.Nk*2*np.pi
            elif type(idx) is float:
                phi = idx
        else:
            phi = ""

    # Fetch the data.
    # r(ight) and l(eft) data arrays
    Tempr, Templ = gsph.MerTemp(nStp, indx=indx)

    # Fetch the coordinates of the grid cell corners in the specified
    # meridional plane, in the GH(MJDc) frame.
    xr, yr, zr, xl, yl, zl, r = gsph.MeridGridHalfs(*indx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originially in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame.
        cr = SkyCoord(
            -xr*u.Rsun, -yr*u.Rsun, zr*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )
        cl = SkyCoord(
            -xl*u.Rsun, -yl*u.Rsun, zl*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        cr = cr.transform_to(hgs_frame)
        cl = cl.transform_to(hgs_frame)

        # Extract the converted coordinates.
        xr = dm.dmarray(cr.cartesian.x)
        yr = dm.dmarray(cr.cartesian.y)
        zr = dm.dmarray(cr.cartesian.z)
        xl = dm.dmarray(cl.cartesian.x)
        yl = dm.dmarray(cl.cartesian.y)
        zl = dm.dmarray(cl.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(np.sqrt(xr**2 + yr**2), zr, Templ, cmap=TCM, norm=vT)
        Ax.pcolormesh(-np.sqrt(xl**2 + yl**2), zl, Tempr, cmap=TCM, norm=vT)

    else:
        Ax.pcolormesh(np.sqrt(xr**2 + yr**2), zr, Tempr, cmap=TCM, norm=vT)
        Ax.pcolormesh(-np.sqrt(xl**2 + yl**2), zl, Templ, cmap=TCM, norm=vT)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r'Temperature $T$ [$(r/r_0) MK]$')
        Ax.set_xlabel(r"$R_{XY}$ [$R_S$] at $\phi=" +
                      f"{phi:{2}.{2}}$ [$rad$]")
        Ax.set_ylabel("Z [$R_S$]")

    # Return the data.
    return Tempr, Templ


def PlotEqD(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        hgsplot=False, MJDc=None, MJD_plot=None
):
    """Plot normalized solar wind number density in the solar equatorial plane.

    Plot normalized solar wind number density in the solar equatorial plane. By
    default, the plot is produced in the GH(MJDc) frame (the gamhelio frame
    used for the simulation results). If hgsplot is True, MJDc and MJD_plot
    must be specified. In that case, the coordinates are mapped from the
    GH(MJDc) frame to the HGS(MJD_plot) frame.

    The density is normalized with the factor (r/r0)**2, where r0 is radius of
    the inner edge of the gamhelio grid (should be 21.5 Rsun).

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that, at any particular MJD:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc (the MJD of the central meridian in the WSA
    file used for initial conditions), meaning it is fixed in spatial
    orientation at that time. The HGS frame is defined at MJD_plot, also
    producing a (different) fixed spatial orientation. The conversion maps
    points in the GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a
    rotation about the z-axis, but also accounting for the eccentricity of
    Earth's orbit and other astronomical parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting
    doDeco : bool
        If True, add axis labels and other decorations to the plot
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot

    Returns
    -------
    MagV : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vD = kv.genNorm(DMin, DMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vD, cbT=None, cM=DCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    NormD = gsph.eqNormD(nStp)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originally in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame. All z values are
        # set to 0 since we are plotting in the equatorial (XY) plane.
        zzi = np.zeros_like(gsph.xxi)
        c = SkyCoord(
            -gsph.xxi*u.Rsun, -gsph.yyi*u.Rsun, zzi*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(x, y, NormD, cmap=DCM, norm=vD)

    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, NormD, cmap=DCM, norm=vD)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Number density $n$ [$(r/r_0)^2 cm^{-3}$]")
        Ax.set_xlabel(r"$X$ [$R_S$]")
        Ax.set_ylabel(r"$Y$ [$R_S$]")
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position("right")

    # Return the data.
    return NormD


def PlotjD(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True, jidx=-1,
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot normalized density in a specific j plane.

    Plot normalized density in the a specific j plane. By default, the plot
    is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If MJD_plot is specified, MJDc must also be specified.
    In that case, the coordinates are mapped from the GH(MJDc) frame to the
    HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc, meaning it is fixed in spatial orientation
    at that time. The HGS frame is defined at MJD_plot, also producing a
    (different) fixed spatial orientation. The conversion maps points in the
    GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a rotation about
    the z-axis, but also accounting for the Earth's orbit and other
    astronommical and geodetic parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If true, clear the plot Axes before further plotting.
    doDeco : bool
        If true, add axis labels to the plot.
    jidx : int
        Index of j-plane to plot.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.
    hgsplot : bool
        If true, plot in HGS(MJD_plot) frame.

    Returns
    -------
    NormD : np.array of float
        Normalized number density in selected j-plane, same shape as the
        j-plane in the gamhelio results.

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vD = kv.genNorm(DMin, DMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vD, cbT=None, cM=DCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    NormD = gsph.jNormD(nStp, jidx=jidx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:
        raise TypeError("HGS frame not supported for pic7!")
    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, NormD, cmap=DCM, norm=vD)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Number density $n$ [$(r/r_0)^2 cm^{-3}$]")
        Ax.set_xlabel('$X [R_S]$')
        Ax.set_ylabel('$Y [R_S]$')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the data.
    return NormD

def PlotEqTemp(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot normalized solar wind temperature in the solar equatorial plane.

    Plot normalized solar wind temperature in the solar equatorial plane. By
    default, the plot is produced in the GH(MJDc) frame (the gamhelio frame
    used for the simulation results). If hgsplot is True, MJDc and MJD_plot
    must be specified. In that case, the coordinates are mapped from the
    GH(MJDc) frame to the HGS(MJD_plot) frame.

    The temperature is normalized with the factor (r/r0), where r0 is 21.5 Rsun
    (the inner edge of the gamhelio grid).

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that, at any particular MJD:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc (the MJD of the central meridian in the WSA
    file used for initial conditions), meaning it is fixed in spatial
    orientation at that time. The HGS frame is defined at MJD_plot, also
    producing a (different) fixed spatial orientation. The conversion maps
    points in the GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a
    rotation about the z-axis, but also accounting for the eccentricity of
    Earth's orbit and other astronomical parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting
    doDeco : bool
        If True, add axis labels and other decorations to the plot
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot

    Returns
    -------
    Temp : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vT, cbT=None, cM=TCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    Temp = gsph.eqTemp(nStp)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originially in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame. Set all z values
        # to 0 since we are using the solar equatorial plane.
        zzi = np.zeros_like(gsph.xxi)
        c = SkyCoord(
            -gsph.xxi*u.Rsun, -gsph.yyi*u.Rsun, zzi*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(x, y, Temp, cmap=TCM, norm=vT)

    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, Temp, cmap=TCM, norm=vT)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Temperature $T$ [$(r/r_0) MK$]")
        Ax.set_xlabel(r"$X$ [$R_S$]")
        Ax.set_ylabel(r"$Y$ [$R_S$]")
        Ax.yaxis.tick_left()
        Ax.yaxis.set_label_position('left')
    # Return the data.
    return Temp


def PlotjTemp(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True, jidx=-1,
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot normalized temperature in a specific j plane.

    Plot normalized temperature in the a specific j plane. By default, the plot
    is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If MJD_plot is specified, MJDc must also be specified.
    In that case, the coordinates are mapped from the GH(MJDc) frame to the
    HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc, meaning it is fixed in spatial orientation
    at that time. The HGS frame is defined at MJD_plot, also producing a
    (different) fixed spatial orientation. The conversion maps points in the
    GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a rotation about
    the z-axis, but also accounting for the Earth's orbit and other
    astronommical and geodetic parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If true, clear the plot Axes before further plotting.
    doDeco : bool
        If true, add axis labels to the plot.
    jidx : int
        Index of j-plane to plot.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.
    hgsplot : bool
        If true, plot in HGS(MJD_plot) frame.

    Returns
    -------
    Temp : np.array of float
        Normalized temperature in selected j-plane, same shape as the
        j-plane in the gamhelio results.

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vT, cbT=None, cM=TCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    Temp = gsph.jTemp(nStp, jidx=jidx)
 
    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:
        raise TypeError("HGS frame not supported for pic7!")
    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, Temp, cmap=TCM, norm=vT)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Temperature $T$ [$(r/r_0) MK$]")
        Ax.set_xlabel('$X [R_S]$')
        Ax.set_ylabel('$Y [R_S]$')

    # Return the data.
    return Temp


def PlotEqBr(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot normalized solar wind radial magnetic field in the solar equatorial plane.

    Plot normalized solar wind radial magnetic field in the solar equatorial
    plane. By default, the plot is produced in the GH(MJDc) frame (the gamhelio
    frame used for the simulation results). If hgsplot is True, MJDc and
    MJD_plot must be specified. In that case, the coordinates are mapped from
    the GH(MJDc) frame to the HGS(MJD_plot) frame.

    The radial magnetic field is normalized with the factor (r/r0)**2, where r0
    is 21.5 Rsun (the inner edge of the gamhelio grid).

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that, at any particular MJD:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc (the MJD of the central meridian in the WSA
    file used for initial conditions), meaning it is fixed in spatial
    orientation at that time. The HGS frame is defined at MJD_plot, also
    producing a (different) fixed spatial orientation. The conversion maps
    points in the GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a
    rotation about the z-axis, but also accounting for the eccentricity of
    Earth's orbit and other astronomical parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting
    doDeco : bool
        If True, add axis labels and other decorations to the plot
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot
    Returns
    -------
    Br : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, cbT=None, cM=BCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    Br = gsph.eqNormBr(nStp)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originially in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame. Set all z values
        # to 0 since we are using the solar equatorial plane.
        zzi = np.zeros_like(gsph.xxi)
        c = SkyCoord(
            -gsph.xxi*u.Rsun, -gsph.yyi*u.Rsun, zzi*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(x, y, Br, cmap=BCM, norm=vB)

    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, Br, cmap=BCM, norm=vB)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Radial magnetic field $B_r$ [$(r/r_0)^2 nT$]")
        Ax.set_xlabel(r"$X$ [$R_S$]")
        Ax.set_ylabel(r"$Y$ [$R_S$]")
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')
    # Return the data.
    return Br


def PlotjBr(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True, jidx=-1,
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot normalized radial magnetic field in a specific j plane.

    Plot normalized radial magnetic field in the a specific j plane. By default, the plot
    is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If MJD_plot is specified, MJDc must also be specified.
    In that case, the coordinates are mapped from the GH(MJDc) frame to the
    HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc, meaning it is fixed in spatial orientation
    at that time. The HGS frame is defined at MJD_plot, also producing a
    (different) fixed spatial orientation. The conversion maps points in the
    GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a rotation about
    the z-axis, but also accounting for the Earth's orbit and other
    astronommical and geodetic parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If true, clear the plot Axes before further plotting.
    doDeco : bool
        If true, add axis labels to the plot.
    jidx : int
        Index of j-plane to plot.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.
    hgsplot : bool
        If true, plot in HGS(MJD_plot) frame.

    Returns
    -------
    Br : np.array of float
        Normalized radial magnetic field in selected j-plane, same shape as the
        j-plane in the gamhelio results.

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, cbT=None, cM=BCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    Br = gsph.jNormBr(nStp, jidx=jidx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:
        raise TypeError("HGS frame not supported for pic7!")
    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, Br, cmap=BCM, norm=vB)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Radial magnetic field $B_r$ [$(r/r_0)^2 nT$]")
        Ax.set_xlabel('$X [R_S]$')
        Ax.set_ylabel('$Y [R_S]$')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the data.
    return Br


def PlotEqBx(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot solar wind x-magnetic field in the solar equatorial plane.

    Plot solar wind x-magnetic field in the solar equatorial plane. By default,
    the plot is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If MJD_plot is specified, MJDc must also be specified.
    In that case, the coordinates are mapped from the GH(MJDc) frame to the
    HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc, meaning it is fixed in spatial orientation
    at that time. The HGS frame is defined at MJD_plot, also producing a
    (different) fixed spatial orientation. The conversion maps points in the
    GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a rotation about
    the z-axis, but also accounting for the Earth's orbit and other
    astronommical and geodetic parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If true, clear the plot Axes before further plotting.
    doDeco : bool
        If true, add axis labels to the plot.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.
    hgsplot : bool
        If true, plot in HGS(MJD_plot) frame.

    Returns
    -------
    Bx : np.array of float
        Data for solar wind x-magnetic field in equatorial plane, same shape as
        the equatorial grid in the gamhelio results.

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, cbT=None, cM=BCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    Bx = gsph.eqBx(nStp)

    # Plot the data in the solar equatorial plane.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originially in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame. Set all z values
        # to 0 since we are using the solar equatorial plane.
        zzi = np.zeros_like(gsph.xxi)
        c = SkyCoord(
            -gsph.xxi*u.Rsun, -gsph.yyi*u.Rsun, zzi*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(x, y, Bx, cmap=BCM, norm=vB)

    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, Bx, cmap=BCM, norm=vB)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"x-magnetic field $B_x$ [$(r/r_0)^2 nT$]")
        Ax.set_xlabel('$X [R_S]$')
        Ax.set_ylabel('$Y [R_S]$')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the data.
    return Bx


def PlotEqBy(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot solar wind y-magnetic field in the solar equatorial plane.

    Plot solar wind y-magnetic field in the solar equatorial plane. By default,
    the plot is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If MJD_plot is specified, MJDc must also be specified.
    In that case, the coordinates are mapped from the GH(MJDc) frame to the
    HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc, meaning it is fixed in spatial orientation
    at that time. The HGS frame is defined at MJD_plot, also producing a
    (different) fixed spatial orientation. The conversion maps points in the
    GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a rotation about
    the z-axis, but also accounting for the Earth's orbit and other
    astronommical and geodetic parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If true, clear the plot Axes before further plotting.
    doDeco : bool
        If true, add axis labels to the plot.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.
    hgsplot : bool
        If true, plot in HGS(MJD_plot) frame.

    Returns
    -------
    By : np.array of float
        Data for solar wind y-magnetic field in equatorial plane, same shape as
        the equatorial grid in the gamhelio results.

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, cbT=None, cM=BCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    By = gsph.eqBy(nStp)

    # Plot the data in the solar equatorial plane.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originially in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame. Set all z values
        # to 0 since we are using the solar equatorial plane.
        zzi = np.zeros_like(gsph.xxi)
        c = SkyCoord(
            -gsph.xxi*u.Rsun, -gsph.yyi*u.Rsun, zzi*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(x, y, By, cmap=BCM, norm=vB)

    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, By, cmap=BCM, norm=vB)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"y-magnetic field $B_y$ [$(r/r_0)^2 nT$]")
        Ax.set_xlabel('$X [R_S]$')
        Ax.set_ylabel('$Y [R_S]$')

    # Return the data.
    return By


def PlotEqBz(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True,
        MJDc=None, MJD_plot=None, hgsplot=False
):
    """Plot solar wind z-magnetic field in the solar equatorial plane.

    Plot solar wind z-magnetic field in the solar equatorial plane. By default,
    the plot is produced in the GH(MJDc) frame (the gamhelio frame used for the
    simulation results). If MJD_plot is specified, MJDc must also be specified.
    In that case, the coordinates are mapped from the GH(MJDc) frame to the
    HGS(MJD_plot) frame.

    The gamhelio frame GH is based on the Heliographic Stonyhurst frame (HGS)
    frame. The difference is that:

    x (GH) = -x (HGS)
    y (GH) = -y (HGS)
    z (GH) = z (HGS)

    The GH frame is defined at MJDc, meaning it is fixed in spatial orientation
    at that time. The HGS frame is defined at MJD_plot, also producing a
    (different) fixed spatial orientation. The conversion maps points in the
    GH(MJDc) frame to the HGS(MJD_plot) frame, which is almost a rotation about
    the z-axis, but also accounting for the Earth's orbit and other
    astronommical and geodetic parameters.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If true, clear the plot Axes before further plotting.
    doDeco : bool
        If true, add axis labels to the plot.
    MJDc : float
        MJD used for the coordinate GH frame of the simulation.
    MJD_plot : float
        MJD to use for the HGS frame of the plot.
    hgsplot : bool
        If true, plot in HGS(MJD_plot) frame.

    Returns
    -------
    Bz : np.array of float
        Data for solar wind z-magnetic field in equatorial plane, same shape as
        the equatorial grid in the gamhelio results.

    Raises
    ------
    None
    """
    # Fetch the data.
    Bz = gsph.eqBz(nStp)

    # Create a normalizer object for the colorbar.
    maxBz = np.max(np.abs(Bz))
    print(f"maxBz = {maxBz}")
    # vB = kv.genNorm(-maxBz, maxBz, doLog=False, midP=None)
    vB = kv.genNorm(BZMin, BZMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, cbT=None, cM=BCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Plot the data in the solar equatorial plane.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Load the equatorial grid cell vertex coordinates (originially in the
        # GH(MJDc) frame) in the equivalent HGS(MJDc) frame. Set all z values
        # to 0 since we are using the solar equatorial plane.
        zzi = np.zeros_like(gsph.xxi)
        c = SkyCoord(
            -gsph.xxi*u.Rsun, -gsph.yyi*u.Rsun, zzi*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(x, y, Bz, cmap=BCM, norm=vB)

    else:
        Ax.pcolormesh(gsph.xxi, gsph.yyi, Bz, cmap=BCM, norm=vB)

    # Set the plot boundaries.
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"z-magnetic field $B_z$ [$(r/r_0)^2 nT$]")
        Ax.set_xlabel('$X [R_S]$')
        Ax.set_ylabel('$Y [R_S]$')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the data.
    return Bz


def find_radial_slice(gsph, radius):
    """Find the index of the radial slice containing a radius.

    Find the index of the radial slice containing a radius. This is the index
    of the radial grid cell edge just less than the given radius.

    This routine should work for LFM grids of any resolution.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    radius : float
        Radius in Rsun

    Returns
    -------
    idx : int
        Index of radial slice containing radius.

    Raises
    ------
    None
    """
    if( radius == 21.5):
        idx = 0
    else: 
        # Starting with the last radial layer, work inward until the first layer
        # is found with grid radius less than the specified radius.
        idx = -1
        r = np.sqrt(gsph.X[idx][0][0]**2 + gsph.Y[idx][0][0]**2 +
                    gsph.Z[idx][0][0]**2)
        while r > radius:
            idx -= 1
            r = np.sqrt(gsph.X[idx][0][0]**2 + gsph.Y[idx][0][0]**2 +
                        gsph.Z[idx][0][0]**2)

        # Convert the index-from-end to an index-from-start.
        idx += gsph.Ni + 1

    # Return the index of the slice containing the radius.
    return idx


def PlotiSlMagV(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True, idx=-1,
        idx_is_radius=False, hgsplot=False, MJDc=None, MJD_plot=None
):
    """Plot solar wind speed at a specified radial slice.

    Plot solar wind speed at a specified radial slice. The plot uses a
    latitude/longitude coordinate system in the GH(MJDc) frame. The resulting
    plot is, in effect, a view *toward* the origin (the Sun) along the X-axis
    in the GH(MJDc) frame. That is, the X-axis is positive *out* of the page,
    and longitude increases to the right. Think of it as the +X axis out of the
    figure at the middle of the left edge of the image, and the Sun is
    "unrolled" to the right.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting.
    doDeco : bool
        If True, add axis labels to the plot.
    idx : int OR float
        Index of radial slice to plot, OR radius in Rsun.
    idx_is_radius : bool
        If True, interpret idx as a radius in Rsun.
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot

    Returns
    -------
    V : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar.
    vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vMagV, cbT=None, cM=MagVCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    if idx_is_radius:
        radius = idx

        # Determine the index of the radial slice which is just less than the
        # specified radius, and compute the next higher, to bracket the
        # radius.
        i1 = find_radial_slice(gsph, radius)
        i2 = i1 + 1

        # Compute the bracketing radii.
        r1 = np.sqrt(gsph.X[i1][0][0]**2 + gsph.Y[i1][0][0]**2 +
                     gsph.Z[i1][0][0]**2)
        r2 = np.sqrt(gsph.X[i2][0][0]**2 + gsph.Y[i2][0][0]**2 +
                     gsph.Z[i2][0][0]**2)

        # Compute the interpolation slope for each grid point in the layer.
        m = (radius - r1)/(r2 - r1)

        # Fetch the values from the bracketing slices, and interpolate to the
        # specified radius.
        V1 = gsph.iSliceMagV(nStp, idx=i1)
        V2 = gsph.iSliceMagV(nStp, idx=i2)
        V = (1 - m)*V1 + m*V2

        # Fetch the latitude and longitude of the grid cells in the radial
        # slice. This is in the GH(MJDc) frame.
        lat, lon = gsph.iSliceGrid(idx=i1)

    else:

        # Fetch the values from the specified layer.
        V = gsph.iSliceMagV(nStp, idx=idx)

        # Fetch the latitude and longitude of the grid cells in the radial
        # slice. This is in the GH(MJDc) frame.
        lat, lon = gsph.iSliceGrid(idx=idx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Use the specified radius, or compute the radius of this i-slice from
        # the coordinates of the first point in the slice.
        if idx_is_radius:
            rg = idx
        else:
            rg = np.sqrt(gsph.X[idx, 0, 0]**2 + gsph.Y[idx, 0, 0]**2 +
                         gsph.Z[idx, 0, 0]**2)

        # Convert the lat/lon at the radial distance to Cartesian coordinates.
        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)
        xg = rg*np.cos(lat_rad)*np.cos(lon_rad)
        yg = rg*np.cos(lat_rad)*np.sin(lon_rad)
        zg = rg*np.sin(lat_rad)

        # Load the grid cell vertex coordinates (originally in the GH(MJDc)
        # frame) in the equivalent HGS(MJDc) frame.
        c = SkyCoord(
            -xg*u.Rsun, -yg*u.Rsun, zg*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Convert back to lat/lon coordinates.
        # -180 <= lon <= +180 deg
        rxy = np.sqrt(x**2 + y**2)
        lat_rad = np.arctan2(z, rxy)
        lon_rad = np.arctan2(y, x)
        lat = np.degrees(lat_rad)
        lon = np.degrees(lon_rad)

        # The arrays of lon and lat and data must now be rotated about axis 1
        # so that lon values increase monotonically to the right. Find the
        # index of the first *drop* in lon, and shift that index back to index
        # 0 for all 3 arrays.
        dlon = lon[0, 1:] - lon[0, :-1]
        i_shift = np.where(dlon < 0)[0][0]
        i_shift += 1
        lon = np.roll(lon, -i_shift, axis=1)
        lat = np.roll(lat, -i_shift, axis=1)
        V = np.roll(V, -i_shift, axis=1)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(lon, lat, V, cmap=MagVCM, norm=vMagV)

    else:
        Ax.pcolormesh(lon, lat, V, cmap=MagVCM, norm=vMagV)

    # Set the plot boundaries.
    if hgsplot:
        xyBds[0] = -180.0
        xyBds[1] = 180.0
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Speed [$km/s$]")
        Ax.set_xlabel(r"Longitude [$deg$]")
        Ax.set_ylabel(r"Latitude [$deg$]")

    # Return the data.
    return V


def PlotiSlD(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True, idx=-1,
        idx_is_radius=False, hgsplot=False, MJDc=None, MJD_plot=None,
        use_outer_range=False
):
    """Plot solar wind number density at a specified radial slice.

    Plot solar wind number density at a specified radial slice. The plot uses a
    latitude/longitude coordinate system in the GH(MJDc) frame. The resulting
    plot is, in effect, a view *toward* the origin (the Sun) along the X-axis
    in the GH(MJDc) frame. That is, the X-axis is positive *out* of the page,
    and longitude increases to the right. Think of it as the +X axis out of the
    figure at the middle of the left edge of the image, and the Sun is
    "unrolled" to the right.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting.
    doDeco : bool
        If True, add axis labels to the plot.
    idx : int OR float
        Index of radial slice to plot, OR radius in Rsun.
    idx_is_radius : bool
        If True, interpret idx as a radius in Rsun.
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot
    use_outer_range : bool
        If True, use a colorbar range optimized for use at large radius.

    Returns
    -------
    D : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar. If the radial slice is in
    # the outer regions of the grid, switch to a smaller data range.
    if use_outer_range:
        vD = kv.genNorm(D0Min_outer, D0Max_outer, doLog=False, midP=None)
    else:
        vD = kv.genNorm(D0Min, D0Max, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vD, cbT=None, cM=D0CM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    if idx_is_radius:
        radius = idx

        # Determine the index of the radial slice which is just less than the
        # specified radius, and compute the next higher, to bracket the
        # radius.
        i1 = find_radial_slice(gsph, radius)
        i2 = i1 + 1

        # Compute the bracketing radii.
        r1 = np.sqrt(gsph.X[i1][0][0]**2 + gsph.Y[i1][0][0]**2 +
                     gsph.Z[i1][0][0]**2)
        r2 = np.sqrt(gsph.X[i2][0][0]**2 + gsph.Y[i2][0][0]**2 +
                     gsph.Z[i2][0][0]**2)

        # Compute the interpolation slope for each grid point in the layer.
        m = (radius - r1)/(r2 - r1)

        # Fetch the values from the bracketing slices, and interpolate to the
        # specified radius.
        D1 = gsph.iSliceD(nStp, idx=i1)
        D2 = gsph.iSliceD(nStp, idx=i2)
        D = (1 - m)*D1 + m*D2

        # Fetch the latitude and longitude of the grid cells in the radial
        # slice. This is in the GH(MJDc) frame.
        lat, lon = gsph.iSliceGrid(idx=i1)

    else:

        # Fetch the values from the specified layer.
        D = gsph.iSliceD(nStp, idx=idx)

        # Fetch the latitude and longitude of the grid cells in the radial
        # slice. This is in the GH(MJDc) frame.
        lat, lon = gsph.iSliceGrid(idx=idx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Use the specified radius, or compute the radius of this i-slice from
        # the coordinates of the first point in the slice.
        if idx_is_radius:
            rg = idx
        else:
            rg = np.sqrt(gsph.X[idx, 0, 0]**2 + gsph.Y[idx, 0, 0]**2 +
                         gsph.Z[idx, 0, 0]**2)

        # Convert the lat/lon at the radial distance to Cartesian coordinates.
        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)
        xg = rg*np.cos(lat_rad)*np.cos(lon_rad)
        yg = rg*np.cos(lat_rad)*np.sin(lon_rad)
        zg = rg*np.sin(lat_rad)

        # Load the grid cell vertex coordinates (originally in the GH(MJDc)
        # frame) in the equivalent HGS(MJDc) frame.
        c = SkyCoord(
            -xg*u.Rsun, -yg*u.Rsun, zg*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Convert back to lat/lon coordinates.
        # -180 <= lon <= +180 deg
        rxy = np.sqrt(x**2 + y**2)
        lat_rad = np.arctan2(z, rxy)
        lon_rad = np.arctan2(y, x)
        lat = np.degrees(lat_rad)
        lon = np.degrees(lon_rad)

        # The arrays of lon and lat and data must now be rotated about axis 1
        # so that lon values increase monotonically to the right. Find the
        # index of the first *drop* in lon, and shift that index back to index
        # 0 for all 3 arrays.
        dlon = lon[0, 1:] - lon[0, :-1]
        i_shift = np.where(dlon < 0)[0][0]
        i_shift += 1
        lon = np.roll(lon, -i_shift, axis=1)
        lat = np.roll(lat, -i_shift, axis=1)
        D = np.roll(D, -i_shift, axis=1)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(lon, lat, D, cmap=D0CM, norm=vD)

    else:
        Ax.pcolormesh(lon, lat, D, cmap=D0CM, norm=vD)

    # Set the plot boundaries.
    if hgsplot:
        xyBds[0] = -180.0
        xyBds[1] = 180.0
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title("Number density $n$ [$cm^{-3}$]")
        Ax.set_xlabel(r"Longitude [$deg$]")
        Ax.set_ylabel(r"Latitude [$deg$]")
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the data.
    return D


def PlotiSlBr(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True, idx=-1,
        idx_is_radius=False, hgsplot=False, MJDc=None, MJD_plot=None,
        use_outer_range=False
):
    """Plot solar wind radial magnetic field and current sheet at a specified radial slice.

    Plot solar wind radial magnetic field and current sheet at a specified
    radial slice.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting.
    doDeco : bool
        If True, add axis labels to the plot.
    idx : int OR float
        Index of radial slice to plot, OR radius in Rsun.
    idx_is_radius : bool
        If True, interpret idx as a radius in Rsun.
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot
    use_outer_range : bool
        If True, use a colorbar range optimized for use at large radius.

    Returns
    -------
    Br : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar. If the radial slice is in
    # the outer regions of the grid, switch to a smaller data range.
    if use_outer_range:
        vB = kv.genNorm(B0Min_outer, B0Max_outer, doLog=False, midP=None)
    else:
        vB = kv.genNorm(B0Min, B0Max, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, cbT=None, cM=BCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    if idx_is_radius:
        radius = idx

        # Determine the index of the radial slice which is just less than the
        # specified radius, and compute the next higher, to bracket the
        # radius.
        i1 = find_radial_slice(gsph, radius)
        i2 = i1 + 1

        # Compute the bracketing radii.
        r1 = np.sqrt(gsph.X[i1][0][0]**2 + gsph.Y[i1][0][0]**2 +
                     gsph.Z[i1][0][0]**2)
        r2 = np.sqrt(gsph.X[i2][0][0]**2 + gsph.Y[i2][0][0]**2 +
                     gsph.Z[i2][0][0]**2)

        # Compute the interpolation slope for each grid point in the layer.
        m = (radius - r1)/(r2 - r1)

        # Fetch the values from the bracketing slices, and interpolate to the
        # specified radius.
        Br1 = gsph.iSliceBr(nStp, idx=i1)
        Br2 = gsph.iSliceBr(nStp, idx=i2)
        Br = (1 - m)*Br1 + m*Br2

        # Fetch the latitude and longitude of the grid cells in the radial
        # slice. This is in the GH(MJDc) frame.
        lat, lon = gsph.iSliceGrid(idx=i1)

    else:

        # Fetch the values from the specified layer.
        Br = gsph.iSliceBr(nStp, idx=idx)

        # Fetch the latitude and longitude of the grid cells in the radial
        # slice. This is in the GH(MJDc) frame.
        lat, lon = gsph.iSliceGrid(idx=idx)

    # Compute cell-centered lon lat coordinates for contour plot.
    lonc = 0.25*(lon[:-1, :-1] + lon[:-1, 1:] + lon[1:, :-1] + lon[1:, 1:])
    latc = 0.25*(lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Use the specified radius, or compute the radius of this i-slice from
        # the coordinates of the first point in the slice.
        if idx_is_radius:
            rg = idx
        else:
            rg = np.sqrt(gsph.X[idx, 0, 0]**2 + gsph.Y[idx, 0, 0]**2 +
                         gsph.Z[idx, 0, 0]**2)

        # Convert the lat/lon at the radial distance to Cartesian coordinates.
        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)
        xg = rg*np.cos(lat_rad)*np.cos(lon_rad)
        yg = rg*np.cos(lat_rad)*np.sin(lon_rad)
        zg = rg*np.sin(lat_rad)

        # Now convert the cell centers to Cartesian coordinates.
        latc_rad = np.radians(latc)
        lonc_rad = np.radians(lonc)
        xc = rg*np.cos(latc_rad)*np.cos(lonc_rad)
        yc = rg*np.cos(latc_rad)*np.sin(lonc_rad)
        zc = rg*np.sin(latc_rad)

        # Load the grid cell vertex coordinates (originally in the GH(MJDc)
        # frame) in the equivalent HGS(MJDc) frame.
        c = SkyCoord(
            -xg*u.Rsun, -yg*u.Rsun, zg*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Load the grid cell center coordinates (originally in the GH(MJDc)
        # frame) in the equivalent HGS(MJDc) frame.
        cc = SkyCoord(
            -xc*u.Rsun, -yc*u.Rsun, zc*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)
        cc = cc.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)
        xc = dm.dmarray(cc.cartesian.x)
        yc = dm.dmarray(cc.cartesian.y)
        zc = dm.dmarray(cc.cartesian.z)

        # Convert back to lat/lon coordinates.
        # -180 <= lon <= +180 deg
        rxy = np.sqrt(x**2 + y**2)
        lat_rad = np.arctan2(z, rxy)
        lon_rad = np.arctan2(y, x)
        lat = np.degrees(lat_rad)
        lon = np.degrees(lon_rad)

        # Do the same for cell centers.
        rxyc = np.sqrt(xc**2 + yc**2)
        latc_rad = np.arctan2(zc, rxyc)
        lonc_rad = np.arctan2(yc, xc)
        latc = np.degrees(latc_rad)
        lonc = np.degrees(lonc_rad)

        # The arrays of lon and lat and data must now be rotated about axis 1
        # so that lon values increase monotonically to the right. Find the
        # index of the first *drop* in lon, and shift that index back to index
        # 0 for all 3 arrays.
        dlon = lon[0, 1:] - lon[0, :-1]
        i_shift = np.where(dlon < 0)[0][0]
        i_shift += 1
        lon = np.roll(lon, -i_shift, axis=1)
        lat = np.roll(lat, -i_shift, axis=1)
        lonc = np.roll(lonc, -i_shift + 1, axis=1)
        latc = np.roll(latc, -i_shift + 1, axis=1)
        Br = np.roll(Br, -i_shift, axis=1)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(lon, lat, Br, cmap=BCM, norm=vB)

        # Draw the Br=0 contour line representing the current sheet.
        # <BUG>
        # This call creates a spurious horizontal line.
        #Ax.contour(lonc, latc, np.roll(Br, 1, axis=1), [0.], colors='black')
        # </BUG>

    else:
        Ax.pcolormesh(lon, lat, Br, cmap=BCM, norm=vB)

        # Draw the Br=0 contour line representing the current sheet.
        Ax.contour(lonc, latc, np.roll(Br, 1), [0.], colors='black')

    # Set the plot boundaries.
    if hgsplot:
        xyBds[0] = -180.0
        xyBds[1] = 180.0
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title(r"Radial magnetic field $B_r$ [$nT$]")
        Ax.set_xlabel(r"Longitude [$deg$]")
        Ax.set_ylabel(r"Latitude [$deg$]")
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the data.
    return Br

#Plot Br and current sheet (Br=0) at certain distance set in iSliceBr
def PlotiSlBrRotatingFrame(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True,idx=-1):
    BMin = -5.
    BMax = 5.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)
    if (AxCB is not None):
        AxCB.clear()
        kv.genCB(AxCB,vB,"Radial magnetic field [nT]",cM=BCM,Ntk=7)
    if (doClear):
        Ax.clear()

    #Br from the i=0
    Br = gsph.iSliceBrBound(nStp,idx=idx)
    lat, lon = gsph.iSliceGrid(idx=idx)
    
    #transform into rotating frame
    #Julian date of the initial map
    jd0 = gsph.MJDs.min()
    jd_c = gsph.MJDs[nStp]
    print (jd0, jd_c)
    #Julian date of the current solution
    time_days = (jd_c - jd0)
    print (time_days)
    omega=2*180./Tsolar

    #for contour cell-centered lon lat coordinates
    lon_c = 0.25*( lon[:-1,:-1]+lon[:-1,1:]+lon[1:,:-1]+lon[1:,1:] )
    lat_c = 0.25*( lat[:-1,:-1]+lat[:-1,1:]+lat[1:,:-1]+lat[1:,1:] )

    phi = lon_c[0,:] 
    phi_prime = (phi-omega*time_days)%(2*180.)

    if np.where(np.ediff1d(phi_prime)<0)[0].size!=0: #for the first map size =0, for other maps size=1
        ind0=np.where(np.ediff1d(phi_prime)<0)[0][0]+1
        #print 'ind = ', ind0
    else:
        ind0=0 # this is for the first map
    print('ind0 = ', ind0)

    Br = np.roll(Br, -ind0, axis = 1)

    Ax.pcolormesh(lon,lat,Br,cmap=BCM,norm=vB)
    Ax.contour(lon_c, lat_c,Br,[0.],colors='black')
    kv.SetAx(xyBds,Ax)

    if (doDeco):
        Ax.set_xlabel('Longitude')
        Ax.set_ylabel('Latitude')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')
        #for pic4
        Ax.set_aspect('equal')
    return Br


def PlotiSlTemp(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True, idx=-1,
        idx_is_radius=False, hgsplot=False, MJDc=None, MJD_plot=None,
        use_outer_range=False
):
    """Plot solar wind temperature at a specified radial slice.

    Plot solar wind temperature at a specified radial slice.

    Parameters
    ----------
    gsph : kaipy.gamhelio.heliosphere.GamsphPipe
        Pipe to simulation results
    nStp : int
        Index of simulation step to use in plot
    xyBds : list of 4 float
        Minimum and maximum values to plot for x- and y-axes
    Ax : matplotlib.axes.Axes
        Axes object to use for plot
    AxCB : matplotlib.axes.Axes
        Axes object to use for color bar
    doClear : bool
        If True, clear the plot Axes before further plotting.
    doDeco : bool
        If True, add axis labels to the plot.
    idx : int OR float
        Index of radial slice to plot, OR radius in Rsun.
    idx_is_radius : bool
        If True, interpret idx as a radius in Rsun.
    hgsplot : bool
        If True, plot in HGS(MJD_plot) frame
    MJDc : float
        MJD used for the GH frame of the simulation
    MJD_plot : float
        MJD to use for the HGS frame of the plot
    use_outer_range : bool
        If True, use a colorbar range optimized for use at large radius.

    Returns
    -------
    Temp : np.array of float
        Data values used in plot

    Raises
    ------
    None
    """
    # Create a normalizer object for the colorbar. If the radial slice is in
    # the outer regions of the grid, switch to a smaller data range.
    if use_outer_range:
        vT = kv.genNorm(T0Min_outer, T0Max_outer, doLog=False, midP=None)
    else:
        vT = kv.genNorm(T0Min, T0Max, doLog=False, midP=None)

    # Create the color bar.
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vT, cbT=None, cM=TCM, Ntk=7)

    # Clear the plot Axes.
    if doClear:
        Ax.clear()

    # Fetch the data.
    if idx_is_radius:
        radius = idx

        # Determine the index of the radial slice which is just less than the
        # specified radius, and compute the next higher, to bracket the
        # radius.
        i1 = find_radial_slice(gsph, radius)
        i2 = i1 + 1

        # Compute the bracketing radii.
        r1 = np.sqrt(gsph.X[i1][0][0]**2 + gsph.Y[i1][0][0]**2 +
                     gsph.Z[i1][0][0]**2)
        r2 = np.sqrt(gsph.X[i2][0][0]**2 + gsph.Y[i2][0][0]**2 +
                     gsph.Z[i2][0][0]**2)

        # Compute the interpolation slope for each grid point in the layer.
        m = (radius - r1)/(r2 - r1)

        # Fetch the values from the bracketing slices, and interpolate to the
        # specified radius.
        Temp1 = gsph.iSliceT(nStp, idx=i1)
        Temp2 = gsph.iSliceT(nStp, idx=i2)
        Temp = (1 - m)*Temp1 + m*Temp2

        # Fetch the latitude and longitude of the grid cells in the radial
        # slice. This is in the GH(MJDc) frame.
        lat, lon = gsph.iSliceGrid(idx=i1)

    else:

        # Fetch the values from the specified layer.
        Temp = gsph.iSliceT(nStp, idx=idx)

        # Fetch the latitude and longitude of the grid cells in the radial
        # slice. This is in the GH(MJDc) frame.
        lat, lon = gsph.iSliceGrid(idx=idx)

    # Plot the data.
    # If the HGS frame was requested, map the grid corner coordinates from the
    # GH(MJDc) frame to the HGS(MJD_plot) frame.
    if hgsplot:

        # Use the specified radius, or compute the radius of this i-slice from
        # the coordinates of the first point in the slice.
        if idx_is_radius:
            rg = idx
        else:
            rg = np.sqrt(gsph.X[idx, 0, 0]**2 + gsph.Y[idx, 0, 0]**2 +
                         gsph.Z[idx, 0, 0]**2)

        # Convert the lat/lon at the radial distance to Cartesian coordinates.
        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)
        xg = rg*np.cos(lat_rad)*np.cos(lon_rad)
        yg = rg*np.cos(lat_rad)*np.sin(lon_rad)
        zg = rg*np.sin(lat_rad)

        # Load the grid cell vertex coordinates (originally in the GH(MJDc)
        # frame) in the equivalent HGS(MJDc) frame.
        c = SkyCoord(
            -xg*u.Rsun, -yg*u.Rsun, zg*u.Rsun,
            frame=frames.HeliographicStonyhurst,
            obstime=ktools.MJD2UT(MJDc),
            representation_type="cartesian"
        )

        # Create a HGS frame for the plot time.
        hgs_frame = frames.HeliographicStonyhurst(
            obstime=ktools.MJD2UT(MJD_plot)
        )

        # Convert the coordinates from HGS(MJDc) to HGS(MJD_plot).
        c = c.transform_to(hgs_frame)

        # Extract the converted coordinates.
        x = dm.dmarray(c.cartesian.x)
        y = dm.dmarray(c.cartesian.y)
        z = dm.dmarray(c.cartesian.z)

        # Convert back to lat/lon coordinates.
        # -180 <= lon <= +180 deg
        rxy = np.sqrt(x**2 + y**2)
        lat_rad = np.arctan2(z, rxy)
        lon_rad = np.arctan2(y, x)
        lat = np.degrees(lat_rad)
        lon = np.degrees(lon_rad)

        # The arrays of lon and lat and data must now be rotated about axis 1
        # so that lon values increase monotonically to the right. Find the
        # index of the first *drop* in lon, and shift that index back to index
        # 0 for all 3 arrays.
        dlon = lon[0, 1:] - lon[0, :-1]
        i_shift = np.where(dlon < 0)[0][0]
        i_shift += 1
        lon = np.roll(lon, -i_shift, axis=1)
        lat = np.roll(lat, -i_shift, axis=1)
        Temp = np.roll(Temp, -i_shift, axis=1)

        # Plot the data in the HGS(MJD_plot) frame.
        Ax.pcolormesh(lon, lat, Temp, cmap=TCM, norm=vT)

    else:
        Ax.pcolormesh(lon, lat, Temp, cmap=TCM, norm=vT)

    # Set the plot boundaries.
    if hgsplot:
        xyBds[0] = -180.0
        xyBds[1] = 180.0
    kv.SetAx(xyBds, Ax)

    # Decorate the plots.
    if doDeco:
        Ax.set_title("Temperature $T$ [$MK$]")
        Ax.set_xlabel(r"Longitude [$deg$]")
        Ax.set_ylabel(r"Latitude [$deg$]")

    # Return the data.
    return Temp

#Plot Density as a function of distance
def PlotDensityProf(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
    if (doClear):
        Ax.clear()

    D = gsph.RadProfDen(nStp)
    rad  = gsph.RadialProfileGrid()

    Ax.plot(rad,D,colorProf)

    if (doDeco):
        Ax.set_xlabel('Radial distance [R_sun]')
        Ax.set_ylabel('Density [cm-3]')
        Ax.set_ylim(250.,450.)
        Ax.set_xlim(20.,220.)
                #Ax.yaxis.tick_right()
                #Ax.yaxis.set_label_position('right')
    return D

#Plot speed as a function of distance
def PlotSpeedProf(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
    if (doClear):
        Ax.clear()
    V = gsph.RadProfSpeed(nStp)
    rad  = gsph.RadialProfileGrid()
    Ax.plot(rad,V,colorProf)

    if (doDeco):
        Ax.set_xlabel('Radial distance [R_sun]')
        Ax.set_ylabel('Speed [km/s]')
        Ax.set_ylim(600.,750.)
        Ax.set_xlim(20.,220.)
    return V

def PlotFluxProf(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
    if (doClear):
        Ax.clear()
    F = gsph.RadProfFlux(nStp)
    rad  = gsph.RadialProfileGrid()
    Ax.plot(rad,F,colorProf)
    
    if (doDeco):
        Ax.set_xlabel('Radial distance [R_sun]')
        Ax.set_ylabel('RhoVr^2')
        Ax.set_ylim(180000.,280000.)
        Ax.set_xlim(20.,220.)
    return F

#Adds MPI contours
#this function is from magnetosphere Viz script. PlotMPI is not used for helio as of now 
def PlotMPI(gsph,Ax,ashd=0.5):
    gCol = mpiCol
    for i in range(gsph.Ri):
        i0 = i*gsph.dNi
        Ax.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],mpiCol,linewidth=cLW,alpha=ashd)

    if (gsph.Rj>1):
        for j in range(1,gsph.Rj):
            j0 = j*gsph.dNj
            Ax.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)
            Ax.plot(gsph.xxi[:,j0],-gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)
        #X-axis (+)
        Ax.plot(gsph.xxi[:,0], gsph.yyi[:,0],gCol,linewidth=cLW,alpha=ashd)
        #X-axis (-)
        j0 = (gsph.Rj)*gsph.dNj
        Ax.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)
            
