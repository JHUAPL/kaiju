#!/usr/bin/env python
"""
nudge_gen.py -- pure-Python port of the NCL Gen_Data_ne30 MERRA2 nudging-file generator.

Replaces  Gen_MERRA2_ne30.csh -> ncl makeIC_se_002.ncl (+ makeIC_se_procedures.ncl,
makeIC_extract_analyses_info.ncl, makeIC_Create_field_Master_List.ncl, ESMF_regridding.ncl,
MAKEIC.f90 WRAPIT'd to MAKEIC.so) for the one configuration GRX runs use:

    fields U,V,T,PS,Q from the CESM-prepared MERRA2 f09 archive (RDA d313003, 8 slots/day),
    horizontal ESMF first-order CONSERVE regrid to ne30pg3 (ncol=48600) with log(PS) interpolated,
    ADJUST_STATE_FROM_TOPO (psadj_se), Q RH clamp (q2rh_se/rh2q_se), vertical interpolation to the
    CAM L145 hybrid levels (vert_quad_opt1_se for T, vert_int_opt1_se for Q, vert_int_opt2_se for
    U/V), MASS_FIX (mass_fixer_se) twice, float32 NETCDF3_64BIT_OFFSET output named
    <CASE>.cam2.i.YYYY-MM-DD-SSSSS.nc with the same variables/attributes as the NCL.

Per slot the processing sequence is exactly the NCL's (see LOG/LogNCL.* of the original):

    PHIS_in = R(MERRA2 PHIS)                          R = ESMF conserve weights / frac_b
    PS_yy   = exp(R(float32 log(PS)))                 (unadjusted regridded PS)
    T_h = R(T);  Q_h = R(Q);  Q_h <- RH clamp on the 72 MERRA2 levels
    PS_adj  = psadj_se(T_h, p_in, p_in_half, PHIS_in, PHIS_topo, PS_yy)
    2x { T_out = vert_quad_opt1_se(...PS_adj);  Q_out = vert_int_opt1_se(...);
         RH clamp Q_out on L145;  PS_adj *= fixmas (mass_fixer_se) }
    T_out, Q_out again from the twice-fixed PS_adj + RH clamp (no further mass fix)
    U_out, V_out = vert_int_opt2_se(R(U)), vert_int_opt2_se(R(V))
    write float32 PHIS(topo) PS T Q U V

Numerics: everything is float64 as in MAKEIC.f90 (real*8) with every Fortran REAL literal
promoted from single precision (MAKEIC.so was built by gfortran 7.5 without -fdefault-real-8),
log(PS) is rounded to float32 before regridding (NCL log of a float is a float), and the only
other float32 stage is the final doubletofloat cast.  The sparse regrid is applied with
scipy.sparse from a pre-generated ESMF weight file (col/row/S/frac_b), so the module runs in an
environment without esmpy/xesmf/numba (numpy + scipy + netCDF4 only).  Weights are generated
only with --weights auto and only when the cached file is absent (ESMF_RegridWeightGen if
available, else esmpy), and are written in the ESMF weight-file format for reuse.

CLI:
    python nudge_gen.py --ref-date 20130301 --num-days 19 --outdir OUT \
        --inputdir /glade/campaign/collections/rda/data/d313003/0.9x1.25 \
        --grid-info <cam.i file with hyam/hybm/hyai/hybi/lat/lon/area on ncol> \
        --topo <ne30pg3 topo file with PHIS> --weights <WgtGrid.nc|auto> [--tmpdir DIR]
        [--times-per-day 8] [--case MERRA2_ne30_L145] [--workers N] [--force] [--slot YYYYMMDD-SSSSS ...]
        [--brief | --quiet]

Importable:  generate(Config(...)) or process_slot(static, day, sec, outpath).
"""

import argparse
import math
import os
import shutil
import socket
import subprocess
import sys
import time as _time
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from netCDF4 import Dataset

__all__ = ["Config", "Static", "generate", "process_slot", "compute_slot", "write_slot", "output_ok",
           "load_static", "load_weights",
           "regrid", "regrid_ps", "build_press", "psadj_se", "vert_quad_opt1_se", "vert_int_opt1_se",
           "vert_int_opt2_se", "q2rh_rh2q_se", "estblf", "mass_fixer_se", "date_to_time", "time_meta",
           "slot_filename", "write_scrip_rectilinear", "ensure_weights"]

TIMES_PER_DAY = 8
DEFAULT_CASE = "MERRA2_ne30_L145"
DEFAULT_INPUTDIR = "/glade/campaign/collections/rda/data/d313003/0.9x1.25"
RDA_TEMPLATE = "{YYYY}/MERRA2_0.9x1.25_{YYYYMMDD}.nc"
DEFAULT_DST_SCRIP = "/glade/campaign/cesm/cesmdata/cseg/inputdata/share/scripgrids/ne30pg3_scrip_170417.nc"
# 8.0.1 first: the NCL reference used it, and 8.6.0 gives slightly different weights
ESMF_VERSION_PINNED = "8.0.1"
KNOWN_ESMF_RWG = [
    "/glade/u/apps/derecho/23.09/spack/opt/spack/esmf/8.0.1/gcc/7.5.0/jr5q/bin/ESMF_RegridWeightGen",
    "/glade/u/apps/derecho/23.09/spack/opt/spack/esmf/8.6.0/gcc/13.2.0/2bya/bin/ESMF_RegridWeightGen",
]

# MAKEIC.so stores every REAL literal as a float32 promoted to real*8; reproduce that exactly.


def c32(v):
    """A single-precision Fortran literal promoted to double."""
    return float(np.float32(v))


GRAVIT_F = c32(9.80616)
DTDZ = c32(-0.0065)
LAPSE0 = -DTDZ
BOLTZ = c32(1.38065e-23)
AVOGAD = c32(6.02214e26)
MWDAIR = c32(28.966)
RGAS = AVOGAD * BOLTZ
RDAIR = RGAS / MWDAIR
RDG = RDAIR / GRAVIT_F               # constant-folded in the .so
LRDG = LAPSE0 * RDG
T_REF1, T_REF2, T_REF3 = 290.5, 255.0, 298.0
Z_REF_BOT, Z_REF_TOP, Z_MIN = 2000.0, 2500.0, 150.0
THRESHOLD = c32(0.001)
PT_LOG = math.log(c32(1.e-10))       # ghost top pressure (log) for the above-top quadratics
BIG36 = c32(1.e+36)
EPSILO = c32(0.622)
OMEPS = 1.0 - EPSILO
ES_TMIN, ES_TMAX, ES_TTRICE = c32(173.16), c32(375.16), c32(20.00)
ES_TFREEZ = c32(273.16)
# NCL-side constants are exact doubles
GRAVIT_NCL = 9.80616
PS0 = 100000.0
TMASS0 = 98222.0 / GRAVIT_NCL
P0 = 100000.0

# the folded values found in MAKEIC.so
assert RDAIR == 287.0423119769592 and RDG == 29.271632602676544 and LRDG == 0.19026561584302937
assert PT_LOG == -23.025850916589025 and OMEPS == 0.37800002098083496

# abort outside these bounds
FIELD_BOUNDS = {"U": (-8000.0, 8000.0), "V": (-8000.0, 8000.0), "T": (-1000.0, 3600.0),
                "Q": (0.0, 0.1), "PS": (30000.0, 109000.0), "PHIS": (-2.e3, 1.e5)}
FIELD_META = {
    "U": ("Zonal wind", "m/s", "conserve", "linear_logp_opt2"),
    "V": ("Meridional wind", "m/s", "conserve", "linear_logp_opt2"),
    "T": ("Temperature", "K", "conserve", "quad_logp_opt1"),
    "PS": ("Surface pressure", "Pa", "conserve", None),
    "Q": ("Specific humidity", "kg/kg", "conserve", "linear_opt1"),
    "PHIS": ("Surface geopotential", "m2/s2", "no_interp", None),
}


def log(msg):
    print(f"[nudge_gen] {msg}", flush=True)


def _vlevel(v):
    """Verbosity: 0 silent, 1 one line per slot, 2 full diagnostics; True/False map to 2/0."""
    if v is True:
        return 2
    if v is False or v is None:
        return 0
    return int(v)


def _gffgch(t, itype):
    """MAKEIC.f90 gffgch: Goff-Gratch es(T) in Pa."""
    tfreez, one = ES_TFREEZ, 1.0
    if itype < 0:
        tr = abs(float(itype))
        itype = 1
    else:
        tr = 0.0
    eswtr = None
    if not (t < (tfreez - tr) and itype == 1):
        ps = c32(1013.246)
        ts = c32(373.16)
        e1 = c32(11.344) * (1.0 - t / ts)
        e2 = c32(-3.49149) * (ts / t - 1.0)
        f1 = c32(-7.90298) * (ts / t - 1.0)
        f2 = c32(5.02808) * math.log10(ts / t)
        f3 = c32(-1.3816) * (10.0 ** e1 - 1.0) / 10000000.0
        f4 = c32(8.1328) * (10.0 ** e2 - 1.0) / 1000.0
        f5 = math.log10(ps)
        f = f1 + f2 + f3 + f4 + f5
        es = (10.0 ** f) * 100.0
        eswtr = es
        if t >= tfreez or itype == 0:
            return es
    t0 = tfreez
    term1 = c32(2.01889049) / (t0 / t)
    term2 = c32(3.56654) * math.log(t0 / t)
    term3 = c32(20.947031) * (t0 / t)
    es = c32(575.185606e10) * math.exp(-(term1 + term2 + term3))
    if t < (tfreez - tr):
        return es
    weight = min((tfreez - t) / tr, one)
    return weight * es + (1.0 - weight) * eswtr


def _build_estbl():
    """gestbl: the es table at 1 K steps from tmin (t accumulated as in the Fortran)."""
    lentbl = int(ES_TMAX - ES_TMIN + c32(2.000001))
    itype = int(-ES_TTRICE)
    estbl = np.full(250, -99999.0)
    t = ES_TMIN - 1.0
    for n in range(lentbl):
        t = t + 1.0
        estbl[n] = _gffgch(t, itype)
    return estbl, lentbl


ESTBL, ES_LENTBL = _build_estbl()


def estblf(td):
    """Vectorized statement function estblf: es [Pa] of temperature td."""
    td = np.asarray(td, dtype=np.float64)
    tl = np.maximum(np.minimum(td, ES_TMAX), ES_TMIN)
    i = (tl - ES_TMIN).astype(np.int64)
    w = (ES_TMIN + i) - tl
    return (w + 1.0) * ESTBL[i] - w * ESTBL[i + 1]


def _qsat(t, press):
    es = estblf(t)
    with np.errstate(divide="ignore", invalid="ignore"):
        qs = (EPSILO * es) / (press - OMEPS * es)
    qs = np.minimum(1.0, qs)
    return np.where(qs < 0.0, 1.0, qs)


def q2rh_rh2q_se(q, t, press, label="", verbose=True):
    """q2rh_se then rh2q_se in place: q <- clip(q/qs, 0, 1)*qs. Returns (nneg, ngt1, rhmax, rhmin)."""
    qs = _qsat(t, press)
    np.divide(q, qs, out=q)
    imax = int(np.argmax(q.ravel(order="C")))
    imin = int(np.argmin(q.ravel(order="C")))
    qmax, qmin = float(q.ravel()[imax]), float(q.ravel()[imin])
    neg = q < 0.0
    gt1 = q > 1.0
    countneg, countgt1 = int(neg.sum()), int(gt1.sum())
    q[neg] = 0.0
    q[gt1] = 1.0
    np.multiply(q, qs, out=q)
    if verbose:
        npts = q.size
        ncol = q.shape[1]
        log(f"{label}Maximum RH = {qmax * 100.:12.5f}% at n,k = {imax % ncol + 1:5d}{imax // ncol + 1:5d}")
        log(f"{label}Minimum RH = {qmin * 100.:12.5f}% at n,k = {imin % ncol + 1:5d}{imin // ncol + 1:5d}")
        if countneg:
            log(f"{label}WARNING in Q2RH: {countneg:10d} points (out of {npts:10d} "
                f"({countneg / npts * 100.:12.5f}%) ) were negative.  All were set to 0")
        if countgt1:
            log(f"{label}WARNING in Q2RH: {countgt1:10d} points (out of {npts:10d} "
                f"({countgt1 / npts * 100.:12.5f}%) ) were  .gt. 1.   All were set to 1")
    return countneg, countgt1, qmax, qmin


def build_press(hya, hyb, ps):
    """p[k,n] = hya[k]*1e5 + hyb[k]*ps[n], as NCL pres_hybrid_ccm."""
    hya = np.asarray(hya, dtype=np.float64)[:, None]
    hyb = np.asarray(hyb, dtype=np.float64)[:, None]
    return hya * P0 + hyb * np.asarray(ps, dtype=np.float64)[None, :]


def load_weights(path):
    """ESMF weight file -> (csr matrix (n_b, n_a), frac_b, src_grid_dims, n_a, n_b).

    The CSR is built directly so each row keeps the weight-file term order NCL sums in;
    scipy's COO path would sort by column and change the rounding."""
    with Dataset(path) as w:
        S = np.asarray(w.variables["S"][:], dtype=np.float64)
        row = np.asarray(w.variables["row"][:], dtype=np.int64) - 1
        col = np.asarray(w.variables["col"][:], dtype=np.int64) - 1
        n_a = len(w.dimensions["n_a"])
        n_b = len(w.dimensions["n_b"])
        src_dims = np.asarray(w.variables["src_grid_dims"][:]) if "src_grid_dims" in w.variables else None
        frac_b = np.asarray(w.variables["frac_b"][:], dtype=np.float64) if "frac_b" in w.variables else None
    if np.any(np.diff(row) < 0):
        order = np.argsort(row, kind="stable")
        row, col, S = row[order], col[order], S[order]
    if row.size and (row[0] < 0 or row[-1] >= n_b or col.min() < 0 or col.max() >= n_a):
        raise ValueError(f"{path}: row/col out of range for n_a={n_a}, n_b={n_b}")
    indptr = np.concatenate([[0], np.cumsum(np.bincount(row, minlength=n_b))]).astype(np.int64)
    M = sp.csr_matrix((S, col.astype(np.int64), indptr), shape=(n_b, n_a))
    if frac_b is None:
        # NCL divides by frac_b, which equals the row sums of S
        frac_b = np.asarray(M.sum(axis=1)).ravel()
    if np.any(frac_b == 0.0):
        raise ValueError(f"{path}: {int((frac_b == 0).sum())} destination cells receive no weight")
    return M, frac_b, src_dims, n_a, n_b


def regrid(M, frac_b, field):
    """Conservative regrid of a (lat,lon) or (lev,lat,lon) field to (ncol,) or (lev,ncol)."""
    x = np.asarray(field, dtype=np.float64)
    if x.ndim == 2:
        return (M @ x.ravel()) / frac_b
    nlev = x.shape[0]
    y = M @ np.ascontiguousarray(x.reshape(nlev, -1).T)      # (ncol, nlev)
    y = y / frac_b[:, None]
    return np.ascontiguousarray(y.T)


def regrid_ps(M, frac_b, ps32):
    """exp(regrid(log(PS))); a float32 PS has its log rounded to float32, as NCL does."""
    ps = np.asarray(ps32)
    lps = np.log(ps.astype(np.float64))
    if ps.dtype != np.float64:
        lps = lps.astype(np.float32).astype(np.float64)
    return np.exp(regrid(M, frac_b, lps))


def _take(a, idx):
    """a[idx[k, n], n] for every (k, n)."""
    ncol = a.shape[1]
    flat = a.reshape(-1) if a.flags["C_CONTIGUOUS"] else np.ascontiguousarray(a).reshape(-1)
    return flat[idx * ncol + np.arange(ncol)]


def _count_le(sorted_cols, queries):
    """For each (k,n), the count of sorted_cols[:,n] <= queries[k,n].

    Uses a two-pointer sweep when queries increase with k, else per-column searchsorted."""
    K, ncol = queries.shape
    L = sorted_cols.shape[0]
    if K > 1 and not np.all(np.diff(queries, axis=0) >= 0.0):
        out = np.empty((K, ncol), dtype=np.int64)
        sc = np.asfortranarray(sorted_cols)
        qc = np.asfortranarray(queries)
        for n in range(ncol):
            out[:, n] = np.searchsorted(sc[:, n], qc[:, n], side="right")
        return out
    sc = sorted_cols if sorted_cols.flags["C_CONTIGUOUS"] else np.ascontiguousarray(sorted_cols)
    flat = sc.reshape(-1)
    cols = np.arange(ncol)
    c = np.zeros(ncol, dtype=np.int64)
    out = np.empty((K, ncol), dtype=np.int64)
    for k in range(K):
        q = queries[k]
        while True:
            live = c < L
            if not live.any():
                break
            idx = np.where(live, c, 0)
            adv = live & (flat[idx * ncol + cols] <= q)
            if not adv.any():
                break
            c = c + adv
        out[k] = c
    return out


def _kbot_tsurf(t, pm, pi, ps_old, required=None):
    """(kbot, tbot, pbot, tsurf) shared by psadj_se and vert_quad_opt1_se, in the Fortran's sum order.

    Only columns marked required may fail the 150 m search; psadj_se skips the others."""
    L, ncol = t.shape
    hkk = (0.5 * (pi[1:] - pi[:-1])) / pm
    zi = (RDG * t) * hkk
    seq = np.empty((2 * L, ncol), dtype=np.float64)
    rev = zi[::-1]
    seq[0::2] = rev
    seq[1::2] = rev
    acc = np.add.accumulate(seq, axis=0)
    ztest = acc[0::2]                       # z after the first half-increment of each level
    above = ztest > Z_MIN
    ok = above.any(axis=0)
    if required is not None:
        ok = ok | ~np.asarray(required, dtype=bool)
    if not ok.all():
        bad = np.flatnonzero(~ok)
        raise RuntimeError(f"could not find model level above {Z_MIN} m in {bad.size} column(s), e.g. {bad[:5]}")
    j = np.argmax(above, axis=0)
    kbot = L - 1 - j
    ar = np.arange(ncol)
    tbot = t[kbot, ar]
    pbot = pm[kbot, ar]
    tmp = LRDG * (ps_old / pbot - 1.0)
    tsurf = tbot * (1.0 + tmp)
    return kbot, tbot, pbot, tsurf


def psadj_se(t, press_m, press_i, phis_old, phis_new, ps_old):
    """MAKEIC.f90 psadj_se: ECMWF-style PS adjustment for the topography change."""
    del_phis = phis_old - phis_new
    small = np.abs(del_phis) <= THRESHOLD
    kbot, tbot, pbot, tsurf = _kbot_tsurf(t, press_m, press_i, ps_old, required=~small)
    lapse = np.full_like(tsurf, LAPSE0)
    t0 = tsurf + (lapse * phis_old) / GRAVIT_F
    c1 = (t0 > T_REF1) & (tsurf <= T_REF1)
    c2 = (t0 > T_REF1) & (tsurf > T_REF1)
    with np.errstate(divide="ignore", invalid="ignore"):
        lapse_c1 = ((T_REF1 - tsurf) * GRAVIT_F) / phis_old
    lapse = np.where(c1, lapse_c1, np.where(c2, 0.0, lapse))
    tsurf = np.where(c2, (T_REF1 + tsurf) * 0.5, tsurf)
    c3 = tsurf < T_REF2
    lapse = np.where(c3, LAPSE0, lapse)
    tsurf = np.where(c3, (T_REF2 + tsurf) * 0.5, tsurf)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        x = (lapse * del_phis) / (GRAVIT_F * tsurf)
        tmp = (1.0 - x / 2.0) + (x * x) / 3.0
        tmp = (del_phis / (RDAIR * tsurf)) * tmp
        adjusted = ps_old * np.exp(tmp)
    return np.where(small, ps_old, adjusted)


def _quad_bracket(pi_m, po_m):
    """The 0-based bracket index kk of vert_quad_opt1_se for every (k,n), capped at L-3.

    Estimated from the level midpoints, then corrected with the exact Fortran test."""
    L, ncol = pi_m.shape
    K = po_m.shape[0]
    a = pi_m[1:L - 1]
    b = pi_m[2:L]
    mid = a + 0.5 * (b - a)
    m = _count_le(mid, po_m)

    def cond(mm):
        mm = np.clip(mm, 0, L - 3)
        aa = _take(pi_m, mm + 1)
        bb = _take(pi_m, mm + 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = (po_m - aa) / (bb - aa)
        return (aa <= po_m) & (beta >= 0.5)

    up = (m <= L - 3) & cond(m)
    m = np.where(up, m + 1, m)
    down = (m >= 1) & ~cond(m - 1)
    m = np.where(down, m - 1, m)
    if ((m <= L - 3) & cond(m)).any() or ((m >= 1) & ~cond(m - 1)).any():
        raise RuntimeError("vert_quad_opt1_se bracket search: midpoint estimate off by more than one")
    return np.minimum(m, L - 3)


def vert_quad_opt1_se(t_old, pressi_m, pressi_i, presso_m, phis_old, ps_old):
    """MAKEIC.f90 vert_quad_opt1_se, quadratic in log p (used for T)."""
    L, ncol = t_old.shape
    _, tbot, pbot_raw, tsurf = _kbot_tsurf(t_old, pressi_m, pressi_i, ps_old)
    kk0 = _quad_bracket(pressi_m, presso_m)
    x1 = _take(t_old, kk0)
    x2 = _take(t_old, kk0 + 1)
    x3 = _take(t_old, kk0 + 2)
    p1 = np.log(_take(pressi_m, kk0))
    p2 = np.log(_take(pressi_m, kk0 + 1))
    p3 = np.log(_take(pressi_m, kk0 + 2))
    pbot = np.log(pbot_raw)[None, :]
    psurf = np.log(ps_old)[None, :]
    px = np.log(presso_m)
    pt = PT_LOG
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # above the first level: quadratic through the ghost top
        tmpt = ((px - p1) * (px - p2)) / ((pt - p1) * (pt - p2))
        tmp1 = ((px - pt) * (px - p2)) / ((p1 - pt) * (p1 - p2))
        tmp2 = ((px - pt) * (px - p1)) / ((p2 - pt) * (p2 - p1))
        b1 = (x1 * tmpt + x1 * tmp1) + x2 * tmp2
        # between pbot and the surface: linear
        b2 = (tbot[None, :] * (psurf - px) + tsurf[None, :] * (px - pbot)) / (psurf - pbot)
        # below the surface: heuristic lapse rate
        t_ref3_bot = tsurf - Z_REF_BOT * DTDZ
        t_ref3_top = np.minimum(tsurf - Z_REF_TOP * DTDZ, T_REF3)
        z = phis_old / GRAVIT_F
        t0 = np.where(z >= Z_REF_TOP, np.minimum(tsurf - z * DTDZ, T_REF3),
                      (t_ref3_bot * (Z_REF_TOP - z) + t_ref3_top * (z - Z_REF_BOT)) / (Z_REF_TOP - Z_REF_BOT))
        lapse = np.where(z >= Z_REF_BOT, np.maximum((t0 - tsurf) / z, 0.0), LAPSE0)
        x = ((lapse * RDAIR) / GRAVIT_F)[None, :] * np.log(presso_m / ps_old[None, :])
        b3 = tsurf[None, :] * (((1.0 + x) + (x * x) / 2.0) + ((x * x) * x) / 6.0)
        # interior: 3-point Lagrange quadratic
        q1 = ((px - p2) * (px - p3)) / ((p1 - p2) * (p1 - p3))
        q2 = ((px - p1) * (px - p3)) / ((p2 - p1) * (p2 - p3))
        q3 = ((px - p1) * (px - p2)) / ((p3 - p1) * (p3 - p2))
        b5 = (x1 * q1 + x2 * q2) + x3 * q3
    c1 = px < p1
    c2 = (px >= pbot) & (px <= psurf)
    c3 = px > psurf
    c4 = px >= p3
    if (c4 & ~c1 & ~c2 & ~c3).any():
        raise RuntimeError("vert_quad_opt1_se: extrapolation below input levels not allowed")
    return np.select([c1, c2, c3], [b1, b2, b3], default=b5)


def vert_int_opt1_se(pressi, presso, xxi):
    """MAKEIC.f90 vert_int_opt1_se, linear in p (used for Q)."""
    L = pressi.shape[0]
    kk0 = np.clip(_count_le(pressi, presso), 1, L - 1) - 1
    x1 = _take(xxi, kk0)
    x2 = _take(xxi, kk0 + 1)
    p1 = _take(pressi, kk0)
    p2 = _take(pressi, kk0 + 1)
    px = presso
    lin = x1 * (p2 - px) + x2 * (px - p1)
    lin = lin / (p2 - p1)
    return np.select([px < p1, px >= p2], [x1, x2], default=lin)


def vert_int_opt2_se(pressi, presso, xxi):
    """MAKEIC.f90 vert_int_opt2_se, in log p (used for U and V)."""
    L = pressi.shape[0]
    kk0 = np.clip(_count_le(pressi, presso), 1, L - 1) - 1
    x1 = _take(xxi, kk0)
    x2 = _take(xxi, kk0 + 1)
    p1 = np.log(_take(pressi, kk0))
    p2 = np.log(_take(pressi, kk0 + 1))
    has3 = (kk0 + 2) <= (L - 1)
    kk2 = np.minimum(kk0 + 2, L - 1)
    x3 = np.where(has3, _take(xxi, kk2), BIG36)
    p3 = np.where(has3, np.log(_take(pressi, kk2)), math.log(BIG36))
    px = np.log(presso)
    pt = PT_LOG
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        xt = (x1 * (p2 - pt) - x2 * (p1 - pt)) / (p2 - p1)
        tmpt = ((px - p1) * (px - p2)) / ((pt - p1) * (pt - p2))
        tmp1 = ((px - pt) * (px - p2)) / ((p1 - pt) * (p1 - p2))
        tmp2 = ((px - pt) * (px - p1)) / ((p2 - pt) * (p2 - p1))
        b1 = (xt * tmpt + x1 * tmp1) + x2 * tmp2
        q1 = ((px - p2) * (px - p3)) / ((p1 - p2) * (p1 - p3))
        q2 = ((px - p1) * (px - p3)) / ((p2 - p1) * (p2 - p3))
        q3 = ((px - p1) * (px - p2)) / ((p3 - p1) * (p3 - p2))
        b3 = (x1 * q1 + x2 * q2) + x3 * q3
        b4 = (x1 * (p2 - px) + x2 * (px - p1)) / (p2 - p1)
    return np.select([px < p1, px >= p2, kk0 == 0], [b1, x2, b3], default=b4)


def mass_fixer_se(q, hyai, hybi, area, ps, verbose=True):
    """MAKEIC.f90 mass_fixer_se: scale ps in place so the dry mass equals TMASS0; returns fixmas.

    Sums are accumulated sequentially in the Fortran's order."""
    area = np.asarray(area, dtype=np.float64)
    hyai = np.asarray(hyai, dtype=np.float64)
    hybi = np.asarray(hybi, dtype=np.float64)
    anorm = np.add.accumulate(area)[-1]
    tmassf = np.add.accumulate(area * ps)[-1]
    qT = np.ascontiguousarray(q.T)                        # (ncol, lev)
    aq = area[:, None] * qT
    term1 = aq * ((hyai[1:] - hyai[:-1]) * PS0)[None, :]
    term2 = aq * ((hybi[1:] - hybi[:-1])[None, :] * ps[:, None])
    qmass1 = np.add.accumulate(term1.ravel())[-1]
    qmass2 = np.add.accumulate(term2.ravel())[-1]
    tmassf = tmassf / (GRAVIT_NCL * anorm)
    qmass1 = qmass1 / (GRAVIT_NCL * anorm)
    qmass2 = qmass2 / (GRAVIT_NCL * anorm)
    fixmas = (TMASS0 + qmass1) / (tmassf - qmass2)
    np.multiply(ps, fixmas, out=ps)
    if verbose:
        log(f"             FIXMAS = {fixmas:18.14f}")
    return fixmas


def date_to_time(d):
    """NCL date_to_time: days of YYYYMMDDSSSSS, keeping its year*365.25 rule so 'time' matches."""
    d = float(d)
    dayspm = [31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]
    if d > 2100000000000.:
        raise ValueError(f"date must be YYYYMMDDSSSSS: {d}")
    year = int(d / 1000000000.)
    tmp = float(year * 100)
    month = int(d / 10000000. - tmp)
    tmp = tmp * 100 + month * 100
    day = int(d / 100000. - tmp)
    tmp = tmp * 100000 + day * 100000
    fraction = (d - tmp) / 86400.
    if year == 4 * int(year / 4.):
        dayspm[1] = 29.
    if month < 1 or month > 12:
        raise ValueError(f"month out of range in {d}")
    if day < 1 or day > dayspm[month - 1]:
        raise ValueError(f"day out of range in {d}")
    t = year * 365.25
    for m in range(month - 1):
        t = t + dayspm[m]
    t = t + day
    t = t + fraction
    return t


def time_meta(day, sec, ref_date):
    """The NCL initialize_IC_dates time variables of one slot."""
    ymd = day.year * 10000 + day.month * 100 + day.day
    ref = int(ref_date)
    t = date_to_time(ymd * 100000. + sec) - date_to_time(ref * 100000.)
    return {"time": t, "ndbase": 0, "nsbase": 0, "nbdate": ref, "nbsec": 0, "ndcur": int(t),
            "nscur": int(sec), "date": ymd, "datesec": int(sec)}


def slot_filename(case, day, sec):
    return f"{case}.cam2.i.{day.year:04d}-{day.month:02d}-{day.day:02d}-{int(sec):05d}.nc"


def parse_ref_date(s):
    s = str(s)
    if len(s) != 8 or not s.isdigit():
        raise ValueError(f"--ref-date must be YYYYMMDD, got {s!r}")
    return date(int(s[:4]), int(s[4:6]), int(s[6:8]))


def write_scrip_rectilinear(path, lat, lon):
    """NCL rectilinear_to_SCRIP for a global lat/lon grid that includes the poles."""
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    nlat, nlon = lat.size, lon.size
    lat2d = np.broadcast_to(lat[:, None], (nlat, nlon))
    lon2d = np.broadcast_to(lon[None, :], (nlat, nlon))

    def corners(v):
        c = np.zeros((nlat + 1, nlon + 1))
        c[1:nlat, 1:nlon] = (v[:-1, :-1] + v[1:, :-1] + v[:-1, 1:] + v[1:, 1:]) / 4.0
        c[0, 1:nlon] = (v[0, 1:] + v[0, :-1]) / 2.0
        c[nlat, 1:nlon] = (v[nlat - 1, 1:] + v[nlat - 1, :-1]) / 2.0
        c[1:nlat, 0] = (v[1:, 0] + v[:-1, 0]) / 2.0
        c[1:nlat, nlon] = (v[1:, nlon - 1] + v[:-1, nlon - 1]) / 2.0
        c[0, 0], c[0, nlon] = v[0, 0], v[0, nlon - 1]
        c[nlat, nlon], c[nlat, 0] = v[nlat - 1, nlon - 1], v[nlat - 1, 0]
        return np.stack([c[:-1, :-1], c[:-1, 1:], c[1:, 1:], c[1:, :-1]], axis=-1).reshape(-1, 4)

    with Dataset(path, "w", format="NETCDF3_CLASSIC") as ds:
        ds.createDimension("grid_size", nlat * nlon)
        ds.createDimension("grid_corners", 4)
        ds.createDimension("grid_rank", 2)
        v = ds.createVariable("grid_dims", "i4", ("grid_rank",))
        v[:] = np.array([nlon, nlat], dtype=np.int32)
        v = ds.createVariable("grid_center_lat", "f8", ("grid_size",))
        v.units = "degrees"
        v[:] = lat2d.ravel()
        v = ds.createVariable("grid_center_lon", "f8", ("grid_size",))
        v.units = "degrees"
        v[:] = lon2d.ravel()
        v = ds.createVariable("grid_imask", "i4", ("grid_size",))
        v.units = "unitless"
        v[:] = np.ones(nlat * nlon, dtype=np.int32)
        v = ds.createVariable("grid_corner_lat", "f8", ("grid_size", "grid_corners"))
        v.units = "degrees"
        v[:] = corners(lat2d)
        v = ds.createVariable("grid_corner_lon", "f8", ("grid_size", "grid_corners"))
        v.units = "degrees"
        v[:] = corners(lon2d)
        ds.date_created = _time.strftime("%a %b %d %I:%M:%S %p %Z %Y")
        ds.Createdby = "nudge_gen.py (port of ESMF_regridding.ncl rectilinear_to_SCRIP)"
        ds.Conventions = "SCRIP"
        ds.title = f"rectilinear_to_SCRIP ({nlat},{nlon})"
    return path


def _find_esmf_rwg(explicit=None):
    """ESMF_RegridWeightGen: explicit or $ESMF_REGRIDWEIGHTGEN, else pinned 8.0.1, else PATH.

    A missing explicit binary raises rather than falling back, since another build changes the output."""
    for given, what in ((explicit, "--esmf-bin"), (os.environ.get("ESMF_REGRIDWEIGHTGEN"),
                                                   "$ESMF_REGRIDWEIGHTGEN")):
        if given:
            if os.path.isfile(given) and os.access(given, os.X_OK):
                return os.path.abspath(given)   # ESMF runs with cwd=tmpdir
            raise FileNotFoundError(f"{what}={given} is not an executable file; the pinned reference "
                                    f"build is {KNOWN_ESMF_RWG[0]}")
    for c in KNOWN_ESMF_RWG + [shutil.which("ESMF_RegridWeightGen")]:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    return None


def weights_version(path):
    try:
        with Dataset(path) as w:
            return str(getattr(w, "CVS_revision", getattr(w, "ESMF_version", "unknown")))
    except Exception:
        return "unknown"


def _check_weights_version(path, strict, what):
    """Warn, or raise when strict, if a weight file was not built by the pinned ESMF.

    Other builds change the output in the last float32 bit while every validator still passes."""
    ver = weights_version(path)
    if ver == ESMF_VERSION_PINNED:
        return ver
    msg = (f"{what} {path} was built with ESMF {ver}, not the pinned {ESMF_VERSION_PINNED}; "
           f"outputs would differ from the NCL reference set in the last float32 bit")
    if strict:
        raise RuntimeError(msg + " (--strict-esmf). Delete the cached weight file and/or pass "
                           f"--esmf-bin {KNOWN_ESMF_RWG[0]}")
    log("WARNING: " + msg)
    return ver


def _gen_weights_esmf_cli(rwg, src, dst, wgt, cwd):
    """Run ESMF_RegridWeightGen as NCL's ESMF_regrid_gen_weights does, in cwd.

    No --pole: the NCL's Pole attribute never reached its command line, so adding it would change
    the weights."""
    src, dst, wgt, cwd = (os.path.abspath(str(p)) for p in (src, dst, wgt, cwd))
    cmd = [rwg, "--source", src, "--destination", dst, "--weight", wgt, "--method", "conserve", "-i"]
    log("running: " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    out = r.stdout + r.stderr
    if r.returncode != 0 or "Completed weight generation successfully" not in out or not os.path.exists(wgt):
        raise RuntimeError(f"ESMF_RegridWeightGen failed (rc={r.returncode}):\n{out[-3000:]}")
    return wgt


def _gen_weights_esmpy(lat, lon, dst_scrip, wgt):
    """Fallback weight generation with esmpy, written in the ESMF weight-file layout.

    Not bitwise equal to ESMF_RegridWeightGen 8.0.1 (weights differ at ~1e-6)."""
    import esmpy
    esmpy.Manager(debug=False)
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    nlat, nlon = lat.size, lon.size
    latc = np.concatenate([[lat[0]], (lat[1:] + lat[:-1]) / 2.0, [lat[-1]]])
    lonc = np.concatenate([[lon[0]], (lon[1:] + lon[:-1]) / 2.0])
    g = esmpy.Grid(np.array([nlon, nlat], dtype=np.int32),
                   staggerloc=[esmpy.StaggerLoc.CENTER, esmpy.StaggerLoc.CORNER],
                   coord_sys=esmpy.CoordSys.SPH_DEG, num_peri_dims=1, periodic_dim=0, pole_dim=1)
    x = g.get_coords(0, staggerloc=esmpy.StaggerLoc.CENTER)
    y = g.get_coords(1, staggerloc=esmpy.StaggerLoc.CENTER)
    x[...] = lon[:, None]
    y[...] = lat[None, :]
    xc = g.get_coords(0, staggerloc=esmpy.StaggerLoc.CORNER)
    yc = g.get_coords(1, staggerloc=esmpy.StaggerLoc.CORNER)
    xc[...] = lonc[:xc.shape[0], None]
    yc[...] = latc[None, :yc.shape[1]]
    mesh = esmpy.Mesh(filename=str(dst_scrip), filetype=esmpy.FileFormat.SCRIP)
    sf = esmpy.Field(g, staggerloc=esmpy.StaggerLoc.CENTER)
    df = esmpy.Field(mesh, meshloc=esmpy.MeshLoc.ELEMENT)
    rh = esmpy.Regrid(sf, df, regrid_method=esmpy.RegridMethod.CONSERVE,
                      unmapped_action=esmpy.UnmappedAction.IGNORE, norm_type=esmpy.NormType.DSTAREA, factors=True)
    fl, fi = rh.get_factors(deep_copy=True)
    fi = np.asarray(fi)
    fl = np.asarray(fl, dtype=np.float64)
    col = fi[:, 0] if fi.ndim == 2 and fi.shape[1] == 2 else fi[0]
    row = fi[:, 1] if fi.ndim == 2 and fi.shape[1] == 2 else fi[1]
    n_a, n_b = nlat * nlon, df.data.size
    order = np.lexsort((col, row))
    col, row, fl = col[order], row[order], fl[order]
    M = sp.csr_matrix((fl, (row - 1, col - 1)), shape=(n_b, n_a))
    frac_b = np.asarray(M.sum(axis=1)).ravel()
    with Dataset(wgt, "w", format="NETCDF3_64BIT_OFFSET") as ds:
        ds.createDimension("n_a", n_a)
        ds.createDimension("n_b", n_b)
        ds.createDimension("n_s", fl.size)
        ds.createDimension("src_grid_rank", 2)
        ds.createDimension("dst_grid_rank", 1)
        ds.createVariable("src_grid_dims", "i4", ("src_grid_rank",))[:] = np.array([nlon, nlat], dtype=np.int32)
        ds.createVariable("dst_grid_dims", "i4", ("dst_grid_rank",))[:] = np.array([n_b], dtype=np.int32)
        ds.createVariable("frac_b", "f8", ("n_b",))[:] = frac_b
        ds.createVariable("col", "i4", ("n_s",))[:] = col.astype(np.int32)
        ds.createVariable("row", "i4", ("n_s",))[:] = row.astype(np.int32)
        ds.createVariable("S", "f8", ("n_s",))[:] = fl
        ds.title = "nudge_gen.py esmpy weight file"
        ds.normalization = "destarea"
        ds.map_method = "Conservative remapping"
        ds.ESMF_regrid_method = "First-order Conservative"
        ds.grid_file_dst = str(dst_scrip)
        ds.ESMF_version = getattr(esmpy, "__version__", "unknown")
    for o in (rh, sf, df, mesh, g):
        try:
            o.destroy()
        except Exception:
            pass
    return wgt


def ensure_weights(weights, tmpdir, lat, lon, ncol, dst_scrip=DEFAULT_DST_SCRIP, esmf_bin=None,
                   strict_esmf=False):
    """The weight file: an explicit path, or for 'auto' a cached file in tmpdir, generated if absent."""
    if weights != "auto":
        if not os.path.isfile(weights):
            raise FileNotFoundError(f"weight file not found: {weights} (use --weights auto to generate one)")
        _check_weights_version(weights, strict_esmf, "the weight file")
        return str(weights)
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    nlat, nlon = len(lat), len(lon)
    wgt = tmpdir / f"WgtGrid_{nlat}x{nlon}_{ncol:08d}_conserve_none.nc"
    if wgt.is_file():
        # the cache key is only the shape; Static checks the grids themselves
        _check_weights_version(wgt, strict_esmf, "the cached weight file")
        return str(wgt)
    src = tmpdir / f"SrcGrid_{nlat}x{nlon}.nc"
    if not src.is_file():
        log(f"creating source SCRIP grid {src}")
        stmp = f"{src}.tmp.{os.getpid()}"       # atomic: concurrent runs must not see a partial file
        write_scrip_rectilinear(stmp, lat, lon)
        os.replace(stmp, src)
    if not os.path.isfile(dst_scrip):
        raise FileNotFoundError(f"destination SCRIP grid not found: {dst_scrip}")
    rwg = _find_esmf_rwg(esmf_bin)
    tmp = str(wgt) + f".tmp.{os.getpid()}"
    if rwg:
        log(f"creating ESMF weights with {rwg}")
        _gen_weights_esmf_cli(rwg, src, dst_scrip, tmp, tmpdir)
    else:
        log("ESMF_RegridWeightGen not found; creating weights with esmpy (not bitwise vs the NCL's ESMF 8.0.1)")
        _gen_weights_esmpy(lat, lon, dst_scrip, tmp)
    try:
        _check_weights_version(tmp, strict_esmf, "the generated weights")
    except RuntimeError:
        os.remove(tmp)
        raise
    os.replace(tmp, wgt)
    return str(wgt)


class Config:
    def __init__(self, ref_date, num_days, outdir, grid_info, topo, weights="auto",
                 inputdir=DEFAULT_INPUTDIR, tmpdir=None, times_per_day=TIMES_PER_DAY, case=DEFAULT_CASE,
                 workers=1, force=False, slots=None, dst_scrip=DEFAULT_DST_SCRIP, esmf_bin=None,
                 verbose=True, strict_esmf=False):
        self.ref_date = parse_ref_date(ref_date) if not isinstance(ref_date, date) else ref_date
        self.num_days = int(num_days)
        self.outdir = str(outdir)
        self.grid_info = str(grid_info)
        self.topo = str(topo)
        self.weights = str(weights)
        self.inputdir = str(inputdir)
        self.tmpdir = str(tmpdir) if tmpdir else os.path.join(self.outdir, "TMP")
        self.times_per_day = int(times_per_day)
        self.case = case
        self.workers = max(1, int(workers))
        self.force = bool(force)
        self.slots = slots          # optional [(date, sec), ...] subset
        self.dst_scrip = dst_scrip
        self.esmf_bin = esmf_bin
        self.strict_esmf = bool(strict_esmf)
        self.verbose = _vlevel(verbose)
        if 86400 % self.times_per_day:
            raise ValueError(f"times-per-day {self.times_per_day} does not divide 86400")
        if self.num_days < 1:
            raise ValueError("num-days must be >= 1")
        last = self.ref_date + timedelta(days=self.num_days - 1)
        if last.year != self.ref_date.year and self.verbose:
            log(f"WARNING: the window {self.ref_date} .. {last} crosses a year boundary.  The NCL's "
                f"date_to_time uses t = year*365.25, so 'time'/'ndcur' are wrong (and not monotonic: "
                f"2012-12-31-75600 gives 0.875 while 2013-01-01-00000 gives 0.25) for slots on the "
                f"far side of Jan 1.  This port reproduces the artifact deliberately -- re-base "
                f"nbdate at Jan 1 (one --ref-date per calendar year; "
                f"nudge_solver.gen_chunks(split_years=True) does this automatically).")

    @property
    def ref_date_int(self):
        return self.ref_date.year * 10000 + self.ref_date.month * 100 + self.ref_date.day

    def input_path(self, day):
        return os.path.join(self.inputdir, RDA_TEMPLATE.format(YYYY=f"{day.year:04d}", YYYYMMDD=day.strftime("%Y%m%d")))

    def all_slots(self):
        step = 86400 // self.times_per_day
        out = []
        for d in range(self.num_days):
            day = self.ref_date + timedelta(days=d)
            for i in range(self.times_per_day):
                out.append((day, i * step))
        if self.slots:
            want = {(d, int(s)) for d, s in self.slots}
            out = [s for s in out if s in want]
            missing = want - set(out)
            if missing:
                raise ValueError(f"requested slot(s) outside the ref-date/num-days range: {sorted(missing)}")
        return out


class Static:
    """Everything shared by all slots: weights, target grid/levels, topography."""

    def __init__(self, cfg):
        with Dataset(cfg.grid_info) as g:
            self.ncol = len(g.dimensions["ncol"])
            self.plevo = len(g.dimensions["lev"])
            self.hyam = np.asarray(g.variables["hyam"][:], dtype=np.float64)
            self.hybm = np.asarray(g.variables["hybm"][:], dtype=np.float64)
            self.hyai = np.asarray(g.variables["hyai"][:], dtype=np.float64)
            self.hybi = np.asarray(g.variables["hybi"][:], dtype=np.float64)
            self.lat = np.asarray(g.variables["lat"][:], dtype=np.float64)
            self.lon = np.asarray(g.variables["lon"][:], dtype=np.float64)
            self.area = np.asarray(g.variables["area"][:], dtype=np.float64)
        for name in ("hyam", "hybm", "hyai", "hybi", "lat", "lon", "area"):
            n = getattr(self, name).size
            want = self.plevo if name in ("hyam", "hybm") else self.plevo + 1 if name in ("hyai", "hybi") else self.ncol
            if n != want:
                raise ValueError(f"{cfg.grid_info}: {name} has {n} entries, expected {want}")
        with Dataset(cfg.topo) as t:
            if len(t.dimensions["ncol"]) != self.ncol:
                raise ValueError(f"{cfg.topo}: ncol={len(t.dimensions['ncol'])} != grid_info ncol={self.ncol}")
            self.phis_out = np.asarray(t.variables["PHIS"][:], dtype=np.float64)
        if not np.isfinite(self.phis_out).all() or np.any(np.abs(self.phis_out) >= 1.e35):
            raise ValueError(f"{cfg.topo}: PHIS has missing values")
        # read the grid from the first requested day, not ref_date, so a --slot retry needs only its own file
        slots = cfg.all_slots()
        if not slots:
            raise ValueError("no slots selected")
        first = cfg.input_path(slots[0][0])
        if not os.path.isfile(first):
            raise FileNotFoundError(f"MERRA2 input missing: {first}")
        with Dataset(first) as m:
            self.mlat = np.asarray(m.variables["lat"][:], dtype=np.float64)
            self.mlon = np.asarray(m.variables["lon"][:], dtype=np.float64)
        if self.mlat[0] >= self.mlat[-1] or self.mlon[0] < -179.0:
            raise ValueError("MERRA2 grid must be S->N with lon starting at 0 (the NCL flips are not ported)")
        self.weights_path = ensure_weights(cfg.weights, cfg.tmpdir, self.mlat, self.mlon, self.ncol,
                                           cfg.dst_scrip, cfg.esmf_bin,
                                           getattr(cfg, "strict_esmf", False))
        self.M, self.frac_b, src_dims, n_a, n_b = load_weights(self.weights_path)
        if n_a != self.mlat.size * self.mlon.size or n_b != self.ncol:
            raise ValueError(f"{self.weights_path}: n_a={n_a}, n_b={n_b} do not match the MERRA2 grid "
                             f"({self.mlat.size}x{self.mlon.size}) / ncol={self.ncol}")
        if src_dims is not None and tuple(int(x) for x in src_dims) != (self.mlon.size, self.mlat.size):
            raise ValueError(f"{self.weights_path}: src_grid_dims {src_dims} != (nlon,nlat)=({self.mlon.size},{self.mlat.size})")
        self._check_weights_provenance(cfg)
        self.cfg = cfg

    def _check_weights_provenance(self, cfg):
        """Check the weight file's source grid and destination SCRIP match this run."""
        with Dataset(self.weights_path) as w:
            w.set_auto_mask(False)
            ver = str(getattr(w, "CVS_revision", getattr(w, "ESMF_version", "unknown")))
            dst = str(getattr(w, "grid_file_dst", getattr(w, "domain_b", "")))
            src = str(getattr(w, "grid_file_src", getattr(w, "domain_a", "")))
            xc = np.asarray(w.variables["xc_a"][:], dtype=np.float64) if "xc_a" in w.variables else None
            yc = np.asarray(w.variables["yc_a"][:], dtype=np.float64) if "yc_a" in w.variables else None
        if cfg.verbose:
            log(f"weights {self.weights_path}: ESMF {ver}, src {src or 'unknown'}, dst {dst or 'unknown'}")
        if ver != ESMF_VERSION_PINNED:
            _check_weights_version(self.weights_path, getattr(cfg, "strict_esmf", False),
                                   "these weights")
        if dst and os.path.basename(dst) != os.path.basename(cfg.dst_scrip):
            raise ValueError(f"{self.weights_path} was generated against {dst}, not the requested "
                             f"{cfg.dst_scrip} -- pass the matching --dst-scrip, or (for "
                             f"--weights auto) delete the cached weight file so it is rebuilt")
        if dst and os.path.realpath(dst) != os.path.realpath(cfg.dst_scrip) and cfg.verbose:
            log(f"NOTE: weight file records dst {dst}; this run asked for {cfg.dst_scrip} "
                f"(same basename, different path)")
        if xc is not None and yc is not None and xc.size == self.mlat.size * self.mlon.size:
            nlat, nlon = self.mlat.size, self.mlon.size
            wx = np.broadcast_to(self.mlon[None, :], (nlat, nlon)).ravel()
            wy = np.broadcast_to(self.mlat[:, None], (nlat, nlon)).ravel()
            if not (np.allclose(xc, wx, rtol=0.0, atol=1.e-6) and np.allclose(yc, wy, rtol=0.0, atol=1.e-6)):
                raise ValueError(f"{self.weights_path}: the source-grid centres in the weight file "
                                 f"are not this MERRA2 grid (max dlon {np.abs(xc - wx).max():.3g}, "
                                 f"max dlat {np.abs(yc - wy).max():.3g}) -- stale cached weights?")


def load_static(cfg):
    return Static(cfg)


def _field(ds, name, idx):
    """One MERRA2 field at idx, in the file's own precision as NCL reads it."""
    a = np.asarray(ds.variables[name][idx])
    return a if a.dtype in (np.float32, np.float64) else a.astype(np.float64)


def _read_slot(path, day, sec):
    fdate = (day.year * 10000 + day.month * 100 + day.day) * 100000 + int(sec)
    with Dataset(path) as m:
        m.set_auto_mask(False)
        dts = np.asarray(m.variables["date"][:], dtype=np.int64) * 100000 + np.asarray(m.variables["datesec"][:], dtype=np.int64)
        hits = np.flatnonzero(dts == fdate)
        if hits.size != 1:
            raise ValueError(f"{path}: desired date {fdate} (YYYYMMDDSSSSS) not found; file holds {dts.tolist()}")
        it = int(hits[0])
        lev = np.asarray(m.variables["lev"][:], dtype=np.float64)
        if lev[0] >= lev[-1]:
            raise ValueError(f"{path}: lev must increase top->bottom")
        d = {
            "it": it,
            "hyam": np.asarray(m.variables["hyam"][:], dtype=np.float64),
            "hybm": np.asarray(m.variables["hybm"][:], dtype=np.float64),
            "hyai": np.asarray(m.variables["hyai"][:], dtype=np.float64),
            "hybi": np.asarray(m.variables["hybi"][:], dtype=np.float64),
            "PS": _field(m, "PS", it),
            "PHIS": _field(m, "PHIS", 0),      # NCL reads PHIS_input at index 0
            "T": _field(m, "T", it),
            "Q": _field(m, "Q", it),
            "U": _field(m, "U", it),
            "V": _field(m, "V", it),
        }
    return d


def _check_field(name, x):
    lo, hi = FIELD_BOUNDS[name]
    if not np.isfinite(x).all():
        raise RuntimeError(f"output field {name} has non-finite values")
    if np.any(x == -9999.0) or np.any(x == -999.0):
        raise RuntimeError(f"output field {name} contains the -9999/-999 sentinel (output2CAM abort)")
    mn, mx = float(x.min()), float(x.max())
    if mn < lo or mx > hi:
        raise RuntimeError(f"output field {name} range [{mn}, {mx}] outside [{lo}, {hi}] (output2CAM abort)")


def compute_slot(st, day, sec, verbose=None):
    """(fields, diagnostics) of the full numerical pipeline for one slot."""
    cfg = st.cfg
    verbose = cfg.verbose if verbose is None else _vlevel(verbose)
    path = cfg.input_path(day)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"MERRA2 input missing for {day}: {path}")
    d = _read_slot(path, day, sec)
    diag = {"input": path, "time_index": d["it"], "rh": [], "fixmas": []}
    M, fb = st.M, st.frac_b

    phis_in = regrid(M, fb, d["PHIS"])
    ps_yy = regrid_ps(M, fb, d["PS"])
    T_h = regrid(M, fb, d["T"])
    Q_h = regrid(M, fb, d["Q"])
    for name, arr in (("PHIS", phis_in), ("PS", ps_yy), ("T", T_h), ("Q", Q_h)):
        if not np.isfinite(arr).all():
            raise RuntimeError(f"regridded {name} has missing values")

    press_in = build_press(d["hyam"], d["hybm"], ps_yy)
    press_in_half = build_press(d["hyai"], d["hybi"], ps_yy)
    if verbose >= 2:
        log(f"{day} {sec:05d}: input {os.path.basename(path)} time index {d['it']}; adjusting Q (RH clamp) on MERRA2 levels")
    diag["rh"].append(q2rh_rh2q_se(Q_h, T_h, press_in, verbose=verbose >= 2))

    if verbose >= 2:
        log("Adjusting PS")
    ps_adj = psadj_se(T_h, press_in, press_in_half, phis_in, st.phis_out, ps_yy)

    for it in (1, 2):
        if verbose >= 2:
            log(f"******* Atmospheric Mass Fixer, iteration {it} *******")
        press_out = build_press(st.hyam, st.hybm, ps_adj)
        T_out = vert_quad_opt1_se(T_h, press_in, press_in_half, press_out, phis_in, ps_yy)
        Q_out = vert_int_opt1_se(press_in, press_out, Q_h)
        diag["rh"].append(q2rh_rh2q_se(Q_out, T_out, press_out, verbose=verbose >= 2))
        diag["fixmas"].append(mass_fixer_se(Q_out, st.hyai, st.hybi, st.area, ps_adj, verbose=verbose >= 2))

    press_out = build_press(st.hyam, st.hybm, ps_adj)
    T_out = vert_quad_opt1_se(T_h, press_in, press_in_half, press_out, phis_in, ps_yy)
    Q_out = vert_int_opt1_se(press_in, press_out, Q_h)
    diag["rh"].append(q2rh_rh2q_se(Q_out, T_out, press_out, verbose=verbose >= 2))

    U_h = regrid(M, fb, d["U"])
    V_h = regrid(M, fb, d["V"])
    if not (np.isfinite(U_h).all() and np.isfinite(V_h).all()):
        raise RuntimeError("regridded U/V has missing values")
    U_out = vert_int_opt2_se(press_in, press_out, U_h)
    V_out = vert_int_opt2_se(press_in, press_out, V_h)

    out = {"PHIS": st.phis_out, "PS": ps_adj, "T": T_out, "Q": Q_out, "U": U_out, "V": V_out}
    for k, v in out.items():
        _check_field(k, v)
    return out, diag


def write_slot(st, day, sec, outpath, fields, input_path):
    """Write the slot file with the NCL output's variables, attributes and order."""
    cfg = st.cfg
    meta = time_meta(day, sec, cfg.ref_date_int)
    tmp = f"{outpath}.tmp.{os.getpid()}"
    with Dataset(tmp, "w", format="NETCDF3_64BIT_OFFSET") as ds:
        ds.createDimension("ncol", st.ncol)
        ds.createDimension("lev", st.plevo)
        ds.createDimension("ilev", st.plevo + 1)
        ds.createDimension("time", None)
        ds.createDimension("scalar", 1)
        for i in range(5, 9):          # unused dimensions the NCL output also has
            ds.createDimension(f"ncl{i}", 1)

        v = ds.createVariable("time", "f8", ("time",), fill_value=9.969209968386869e+36)
        v.long_name, v.units, v.calendar = "time", "days since base date (nbdate)", "noleap"
        for name, ln in (("ndbase", "base day"), ("nsbase", "seconds of base day"),
                         ("nbdate", "base date (YYYYMMDD)"), ("nbsec", "seconds of base date")):
            ds.createVariable(name, "i4", ("scalar",)).long_name = ln
        v = ds.createVariable("ndcur", "i4", ("time",), fill_value=-2147483647)
        v.long_name = "current day (from base day)"
        ds.createVariable("nscur", "i4", ("time",)).long_name = "current seconds of current day"
        ds.createVariable("date", "i4", ("time",)).long_name = "current date (YYYYMMDD)"
        ds.createVariable("datesec", "i4", ("time",)).long_name = "current seconds of current date"
        v = ds.createVariable("P0", "f8", ("scalar",))
        v.long_name, v.units = "reference pressure", "Pa"
        for name, dim, ln, ft in (("lev", "lev", "hybrid level at midpoints (1000*(A+B))", "a: hyam b: hybm p0: P0 ps: PS"),
                                  ("ilev", "ilev", "hybrid level at interfaces (1000*(A+B))", "a: hyai b: hybi p0: P0 ps: PS")):
            v = ds.createVariable(name, "f8", (dim,))
            v.long_name, v.units, v.positive = ln, "level", "down"
            v.standard_name = "atmosphere_hybrid_sigma_pressure_coordinate"
            v.formula_terms = ft
            v.source_file = cfg.grid_info
        for name, dim, ln in (("hyai", "ilev", "hybrid A coefficient at layer interfaces"),
                              ("hybi", "ilev", "hybrid B coefficient at layer interfaces"),
                              ("hyam", "lev", "hybrid A coefficient at layer midpoints"),
                              ("hybm", "lev", "hybrid B coefficient at layer midpoints")):
            v = ds.createVariable(name, "f8", (dim,))
            v.long_name, v.source_file = ln, cfg.grid_info
        for name, ln, un in (("lat", "latitude", "degrees_north"), ("lon", "longitude", "degrees_east"),
                             ("area", "gridpoint area", "normalized")):
            v = ds.createVariable(name, "f8", ("ncol",))
            v.long_name, v.units = ln, un
        for name, ln in (("ntrm", "spectral truncation parameter M"), ("ntrn", "spectral truncation parameter N"),
                         ("ntrk", "spectral truncation parameter K")):
            ds.createVariable(name, "i4", ("scalar",)).long_name = ln
        for name in ("U", "V", "T", "PS", "Q", "PHIS"):
            ln, un, hflag, vflag = FIELD_META[name]
            dims = ("time", "ncol") if name in ("PS", "PHIS") else ("time", "lev", "ncol")
            v = ds.createVariable(name, "f4", dims)
            v.long_name, v.units = ln, un
            v.source_file = cfg.topo if name == "PHIS" else input_path
            v.source_field_name = name
            v.horz_interp_flag = hflag
            if vflag:
                v.vert_interp_flag = vflag

        ds.Conventions = "NCAR-CSM"
        ds.Description = "Re-gridded CAM-SE Initial Conditions from the following source(s):"
        ds.source1 = input_path
        ds.source2 = cfg.topo
        ds.source3 = input_path
        ds.source_grid_info = cfg.grid_info
        ds.hybrid_sigma_pressure = "p(i,j,k) = A(k)*PO + B(k)*PS(i,j)"
        ds.logname = os.environ.get("USER", os.environ.get("LOGNAME", "unknown"))
        ds.host = socket.gethostname()
        ds.hosttype = "x86_64-linux"
        ds.Command = "nudge_gen.py (python port of makeIC_se_002.ncl)"
        ds.Date_Written = _time.strftime("%a %b %d %I:%M:%S %p %Z %Y")
        ds.case = cfg.case
        ds.ps_adjust = "True; Adjust Ps per ECMWF algorithm"
        ds.mass_fix = "True; Fix mass of atmosphere to be consistent with model mass fixer"
        ds.q_adjust = "True; Eliminate super-saturation"
        ds.interp_log_ps = "True; Interpolate ln(Ps) rather than Ps"

        ds.variables["time"][0] = meta["time"]
        for k in ("ndbase", "nsbase", "nbdate", "nbsec"):
            ds.variables[k][0] = meta[k]
        for k in ("ndcur", "nscur", "date", "datesec"):
            ds.variables[k][0] = meta[k]
        ds.variables["P0"][0] = P0
        ds.variables["lev"][:] = (st.hyam + st.hybm) * 1000.
        ds.variables["ilev"][:] = (st.hyai + st.hybi) * 1000.
        ds.variables["hyai"][:] = st.hyai
        ds.variables["hybi"][:] = st.hybi
        ds.variables["hyam"][:] = st.hyam
        ds.variables["hybm"][:] = st.hybm
        ds.variables["lat"][:] = st.lat
        ds.variables["lon"][:] = st.lon
        ds.variables["area"][:] = st.area
        for k in ("ntrm", "ntrn", "ntrk"):
            ds.variables[k][0] = 42
        for name in ("U", "V", "T", "PS", "Q", "PHIS"):
            ds.variables[name][0] = fields[name].astype(np.float32)
    os.replace(tmp, outpath)
    return outpath


def output_ok(path, ncol, nlev):
    """True when path looks like a complete slot file: size, variable shapes, a finite last value."""
    try:
        if os.path.getsize(path) < 4 * ncol * (4 * nlev + 2):
            return False
        with Dataset(path) as ds:
            ds.set_auto_mask(False)
            for name in ("U", "V", "T", "Q", "PS", "PHIS"):
                v = ds.variables[name]
                want = (1, ncol) if name in ("PS", "PHIS") else (1, nlev, ncol)
                if v.dtype != np.float32 or tuple(v.shape) != want:
                    return False
                last = v[0, -1] if len(want) == 2 else v[0, -1, -1]
                if not np.isfinite(np.asarray(last)).all():
                    return False
    except Exception:
        return False
    return True


def process_slot(st, day, sec, outpath=None, force=False):
    """Compute and write one slot; returns (outpath, 'written' or 'skipped', seconds)."""
    cfg = st.cfg
    outpath = outpath or os.path.join(cfg.outdir, slot_filename(cfg.case, day, sec))
    if os.path.exists(outpath) and not (force or cfg.force):
        if output_ok(outpath, st.ncol, st.plevo):
            return outpath, "skipped", 0.0
        if cfg.verbose:
            log(f"WARNING: existing {os.path.basename(outpath)} is truncated or unreadable -- regenerating")
    t0 = _time.perf_counter()
    fields, diag = compute_slot(st, day, sec)
    write_slot(st, day, sec, outpath, fields, diag["input"])
    dt = _time.perf_counter() - t0
    if cfg.verbose:
        log(f"wrote {outpath}  ({dt:.2f} s; FIXMAS {diag['fixmas'][0]:.14f} {diag['fixmas'][1]:.14f})")
    return outpath, "written", dt


_STATIC = None


def _worker(args):
    day, sec = args
    try:
        # generate() has already checked this slot's existing output
        return process_slot(_STATIC, day, sec, force=True)
    except Exception as e:
        raise RuntimeError(f"slot {day} {sec:05d}: {type(e).__name__}: {e}") from e


def generate(cfg):
    """Run every slot of cfg and return [(path, status, seconds)]; complete outputs are skipped."""
    global _STATIC
    os.makedirs(cfg.outdir, exist_ok=True)
    slots = cfg.all_slots()
    days = sorted({d for d, _ in slots})
    missing = [cfg.input_path(d) for d in days if not os.path.isfile(cfg.input_path(d))]
    if missing:
        raise FileNotFoundError("MERRA2 input file(s) missing:\n  " + "\n  ".join(missing))
    st = load_static(cfg)
    if cfg.verbose:
        log(f"{len(slots)} slot(s) {cfg.ref_date} + {cfg.num_days} day(s) -> {cfg.outdir}")
    todo, results = [], []
    for d, sec in slots:
        p = os.path.join(cfg.outdir, slot_filename(cfg.case, d, sec))
        if not cfg.force and os.path.exists(p):
            if output_ok(p, st.ncol, st.plevo):
                results.append((p, "skipped", 0.0))
                continue
            if cfg.verbose:
                log(f"WARNING: existing {os.path.basename(p)} is truncated or unreadable -- regenerating")
        todo.append((d, sec))
    if cfg.workers > 1 and len(todo) > 1:
        import multiprocessing as mp
        _STATIC = st
        with mp.get_context("fork").Pool(min(cfg.workers, len(todo))) as pool:
            results += pool.map(_worker, todo, chunksize=1)
    else:
        for d, s in todo:
            results.append(process_slot(st, d, s, force=True))
    done = [r for r in results if r[1] == "written"]
    if cfg.verbose:
        if done:
            log(f"{len(done)} file(s) written, {len(results) - len(done)} skipped; "
                f"mean {sum(r[2] for r in done) / len(done):.2f} s per slot")
        else:
            log(f"nothing to do ({len(results)} complete file(s) already exist; use --force to regenerate)")
    return results


def _parse_slot(s):
    try:
        d, sec = s.split("-")
        return parse_ref_date(d), int(sec)
    except Exception:
        raise argparse.ArgumentTypeError(f"--slot must be YYYYMMDD-SSSSS, got {s!r}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref-date", required=True, help="first day YYYYMMDD (= nbdate of every output)")
    p.add_argument("--num-days", required=True, type=int)
    p.add_argument("--outdir", required=True)
    p.add_argument("--inputdir", default=DEFAULT_INPUTDIR, help=f"root of {RDA_TEMPLATE} (default {DEFAULT_INPUTDIR})")
    p.add_argument("--grid-info", required=True, help="CAM file with hyam/hybm/hyai/hybi/lat/lon/area on the target grid")
    p.add_argument("--topo", required=True, help="target topography file with PHIS(ncol)")
    p.add_argument("--weights", default="auto", help="ESMF conserve weight file, or 'auto' (cache in --tmpdir)")
    p.add_argument("--tmpdir", default=None, help="cache for SrcGrid/WgtGrid (default <outdir>/TMP)")
    p.add_argument("--times-per-day", type=int, default=TIMES_PER_DAY)
    p.add_argument("--case", default=DEFAULT_CASE)
    p.add_argument("--workers", type=int, default=1, help="process this many slots in parallel (fork)")
    p.add_argument("--force", action="store_true", help="overwrite existing outputs")
    p.add_argument("--slot", action="append", type=_parse_slot, default=None, metavar="YYYYMMDD-SSSSS",
                   help="only generate this slot (repeatable; must lie inside ref-date/num-days)")
    p.add_argument("--dst-scrip", default=DEFAULT_DST_SCRIP, help="destination SCRIP grid for --weights auto")
    p.add_argument("--esmf-bin", default=None, help="ESMF_RegridWeightGen executable for --weights auto")
    p.add_argument("--strict-esmf", action="store_true",
                   help="refuse weights not built by the pinned ESMF %s (generated or cached); the "
                        "PBS wrapper nudge_solver.gen_job writes always passes this" % ESMF_VERSION_PINNED)
    p.add_argument("--quiet", action="store_true", help="print nothing but errors")
    p.add_argument("--brief", action="store_true",
                   help="one line per slot (for a PBS log) instead of the full NCL-style diagnostics")
    a = p.parse_args(argv)
    verbosity = 0 if a.quiet else 1 if a.brief else 2
    t0 = _time.perf_counter()
    try:
        cfg = Config(a.ref_date, a.num_days, a.outdir, a.grid_info, a.topo, a.weights, a.inputdir, a.tmpdir,
                     a.times_per_day, a.case, a.workers, a.force, a.slot, a.dst_scrip, a.esmf_bin,
                     verbosity, a.strict_esmf)
        generate(cfg)
    except (FileNotFoundError, ValueError, RuntimeError, OSError) as e:
        log(f"ERROR: {e}")
        return 1
    if cfg.verbose:
        log(f"total wall time {_time.perf_counter() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
