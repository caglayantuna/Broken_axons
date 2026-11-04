"""Neuron morphology for the 2WPT grafted-neuron images (bmi1 vs scr).

Shared by morphology_2WPT_czi.ipynb (regular CZI, max-Z projection) and
morphology_2WPT_airyscan.ipynb (Airyscan-processed, single plane).

Pipeline per image:
  1. neuron channel (AF594-T1) max-Z projection, resampled to TARGET_UM pixels
  2. background removal (large Gaussian) and noise normalisation -> SNR image
  3. foreground = FAST U-Net++ filament mask (thin extensions) united with SNR > low, kept only
     where the piece contains clearly bright AF594 signal (SNR > high)
  4. cell bodies = Cellpose-SAM at 0.12 µm/px, kept only if soma-sized, compact and AF594-positive
  5. extensions = foreground minus cell bodies -> skeleton -> branches (skan)
  6. extension diameter = FWHM of the intensity profile across each branch (threshold-free)
"""
import os, re, glob
import numpy as np
import pandas as pd
from pylibCZIrw import czi as pyczi
from scipy import ndimage
from skimage.transform import resize
from skimage.filters import gaussian, apply_hysteresis_threshold
from skimage.morphology import remove_small_objects, remove_small_holes, skeletonize, disk, binary_dilation
from skimage.measure import regionprops, label, find_contours

TARGET_UM = 0.0706          # common pixel size (the coarsest native one; Airyscan 35 nm is downsampled)
NEURON_CH = 'AF594-T1'
NUCLEUS_PREFIX = 'DAPI'     # AF405 is not a nuclear stain in this data set

# foreground
BG_SIGMA_UM   = 5.0         # background = Gaussian of this sigma
SNR_LOW, SNR_HIGH = 3.0, 6.0
MIN_OBJ_UM2   = 4.0
FAST_UM       = 0.0877      # FAST U-Net++ runs at the pixel size validated on the in-vitro axons
FAST_DIR      = '/Users/ctuna/Desktop/axon_codes/pretrained_models/fast_actin'
# cell bodies
CP_UM         = 0.12        # Cellpose-SAM pixel size; at coarser scales it misses most somata here
SOMA_MIN_UM2, SOMA_MAX_UM2 = 30.0, 900.0
SOMA_MIN_SNR  = 2.0         # median SNR inside the soma (somata are dimmer than extensions)
SOMA_MIN_SOLIDITY = 0.85
# extensions
SPUR_UM       = 2.0         # end branches shorter than this are pruned
SAMPLE_UM     = 0.5         # diameter sampled every SAMPLE_UM along a branch
PROFILE_UM    = 3.0         # profile half-length across the branch
MIN_SEG_UM    = 3.0         # branches shorter than this get no per-branch statistics
# swellings (Flamier et al. 2018, Fig 3F: % of axons with >= 1 swelling)
SWELL_FACTOR  = 2.0         # local width > SWELL_FACTOR x the branch median width ...
SWELL_MIN_UM  = 0.3         # ... and at least this much wider than the median (ignores noise on thin branches)
SWELL_END_UM  = 2.0         # samples this close to a branch end are not tested (junction blobs)
SWELL_BRIGHT  = 1.2         # a swelling is also brighter than the branch (peak > SWELL_BRIGHT x branch median peak)
# dense graft cores vs isolated neurons
DENSITY_UM    = 20.0        # window for the local AF594 density
DENSE_THR     = 0.20        # local mask fraction above this = dense cluster


# ── files ───────────────────────────────────────────────────────────────────

def czi_meta(path):
    with pyczi.open_czi(path) as cz:
        bb, md = cz.total_bounding_box, cz.raw_metadata
    if bb['X'][1] - bb['X'][0] <= 0:
        return None
    names = re.findall(r'<Channel Id="[^"]*" Name="([^"]+)"', md)[:bb['C'][1]]
    px = re.findall(r'<Distance Id="X">\s*<Value>([^<]+)</Value>', md)
    obj = re.findall(r'<Objective Id="[^"]*" Name="([^"]+)"', md)
    return dict(channels=names, px_um=float(px[0]) * 1e6, n_z=bb['Z'][1],
                size_px=bb['X'][1] - bb['X'][0], objective=obj[0] if obj else '')


def list_images(root, airyscan):
    """One row per readable CZI: mouse, condition, path. airyscan=True -> only 'Airyscan Processing' files."""
    rows = []
    for p in sorted(glob.glob(os.path.join(root, '**', '*.czi'), recursive=True)):
        name = os.path.basename(p)
        if ('Airyscan' in name) != airyscan:
            continue
        meta = czi_meta(p)
        if meta is None or NEURON_CH not in meta['channels']:
            continue
        rel = os.path.relpath(p, root).split(os.sep)
        rows.append(dict(mouse=rel[0], condition='bmi1' if 'bmi1' in rel[1].lower() else 'scr',
                         filename=name, path=p, n_z=meta['n_z'], px_nm=round(meta['px_um'] * 1000, 1),
                         objective=meta['objective'].split(' ')[0] + ' ' + meta['objective'].split(' ')[-1]
                         if meta['objective'] else '', channels=','.join(meta['channels'])))
    return pd.DataFrame(rows)


def _mip(cz, c, n_z):
    return np.stack([np.squeeze(cz.read(plane={'C': c, 'Z': z, 'T': 0})) for z in range(n_z)]).max(0)


def load_neuron(path):
    """AF594 max-Z projection (+ DAPI if present), both resampled to TARGET_UM. Returns (neuron, nucleus|None)."""
    meta = czi_meta(path)
    with pyczi.open_czi(path) as cz:
        neu = _mip(cz, meta['channels'].index(NEURON_CH), meta['n_z']).astype(np.float32)
        nuc_i = [i for i, n in enumerate(meta['channels']) if n.startswith(NUCLEUS_PREFIX) and '#' not in n]
        nuc = _mip(cz, nuc_i[0], meta['n_z']).astype(np.float32) if nuc_i else None
    s = meta['px_um'] / TARGET_UM
    shape = tuple(int(round(d * s)) for d in neu.shape)
    rs = lambda a: resize(a, shape, order=1, anti_aliasing=s < 1, preserve_range=True).astype(np.float32)
    return rs(neu), (rs(nuc) if nuc is not None else None)


# ── foreground ──────────────────────────────────────────────────────────────

def snr_image(img):
    """Background-subtracted image divided by a robust noise estimate."""
    sig = BG_SIGMA_UM / TARGET_UM
    small = resize(img, (img.shape[0] // 8, img.shape[1] // 8), anti_aliasing=True, preserve_range=True)
    bg = resize(gaussian(small, sig / 8, preserve_range=True), img.shape, order=1, preserve_range=True)
    sub = gaussian(img, 1, preserve_range=True) - bg
    noise = 1.4826 * np.median(np.abs(sub - np.median(sub))) + 1e-6
    return (sub / noise).astype(np.float32)


_fast = None
def fast_mask(img):
    """FAST U-Net++ filament mask (classes: actin & stress fibers + filopodia), run at FAST_UM."""
    global _fast
    import sys, torch
    if _fast is None:
        sys.path.insert(0, FAST_DIR)
        from model import build_model
        _fast = build_model()
        _fast.load_state_dict(torch.load(os.path.join(FAST_DIR, 'fast_model.pth'), map_location='cpu'))
        dev = 'mps' if torch.backends.mps.is_available() else ('cuda' if torch.cuda.is_available() else 'cpu')
        _fast = _fast.to(dev).eval()
    dev = next(_fast.parameters()).device
    s = TARGET_UM / FAST_UM
    sh = (int(img.shape[0] * s), int(img.shape[1] * s))
    x = resize(img, sh, anti_aliasing=True, preserve_range=True).astype(np.float32)
    x = np.pad(x, ((0, -sh[0] % 32), (0, -sh[1] % 32)), mode='reflect')
    with torch.no_grad():
        lg = _fast(torch.from_numpy(x)[None, None].to(dev)).float().cpu()
    p = torch.softmax(lg, 1)[0, [1, 4]].sum(0).numpy()[:sh[0], :sh[1]]
    return resize(p, img.shape, order=1) > 0.5


def foreground_mask(snr, img):
    """Pieces of (FAST | SNR > low) that contain at least one pixel with SNR > high."""
    low = fast_mask(img) | (snr > SNR_LOW)
    lbl = label(low, connectivity=2)
    seeds = np.unique(lbl[snr > SNR_HIGH])
    fg = np.isin(lbl, seeds[seeds > 0])
    fg = remove_small_objects(fg, int(MIN_OBJ_UM2 / TARGET_UM ** 2))
    return remove_small_holes(fg, int(2.0 / TARGET_UM ** 2))


# ── cell bodies (Cellpose-SAM) ──────────────────────────────────────────────

_cp = None
def _cellpose():
    global _cp
    if _cp is None:
        from cellpose import models
        _cp = models.CellposeModel(gpu=True)
    return _cp


def soma_mask(img, snr, nucleus=None):
    """Cellpose-SAM at CP_UM pixels, then keep soma-sized, compact, AF594-positive objects. Returns a label image."""
    f = TARGET_UM / CP_UM
    small_shape = (int(snr.shape[0] * f), int(snr.shape[1] * f))
    lo, hi = np.percentile(img, [1, 99.9])
    x = resize(np.clip((img - lo) / (hi - lo + 1e-6), 0, 1), small_shape, anti_aliasing=True, preserve_range=True)
    if False and nucleus is not None:        # DAPI is missing in many images; AF594 alone keeps it consistent
        n = resize(nucleus, small_shape, anti_aliasing=True, preserve_range=True)
        n = np.clip((n - np.percentile(n, 1)) / (np.percentile(n, 99.8) - np.percentile(n, 1) + 1e-6), 0, 1)
        x = np.stack([x, n, np.zeros_like(x)], -1)
    masks, _, _ = _cellpose().eval(x, normalize=True)
    masks = resize(masks, snr.shape, order=0, preserve_range=True, anti_aliasing=False).astype(np.int32)
    keep = np.zeros_like(masks)
    k = 0
    for p in regionprops(masks, intensity_image=snr):
        a = p.area * TARGET_UM ** 2
        if not (SOMA_MIN_UM2 <= a <= SOMA_MAX_UM2) or p.solidity < SOMA_MIN_SOLIDITY:
            continue
        if np.median(snr[p.coords[:, 0], p.coords[:, 1]]) < SOMA_MIN_SNR:
            continue
        k += 1
        keep[p.coords[:, 0], p.coords[:, 1]] = k
    return keep


def soma_table(somata, snr, fg):
    """Shape of each kept cell body + number of primary extensions (AF594 pieces crossing a 1-3 µm ring)."""
    rows = []
    r_in, r_out = int(round(1.0 / TARGET_UM)), int(round(3.0 / TARGET_UM))
    min_piece = int(0.5 / TARGET_UM ** 2)
    for p in regionprops(somata):
        r0, c0, r1, c1 = p.bbox                                   # work in a crop around the soma
        m = r_out + 2
        cr = (slice(max(r0 - m, 0), r1 + m), slice(max(c0 - m, 0), c1 + m))
        sub = (somata[cr] == p.label)
        ring = binary_dilation(sub, disk(r_out)) & ~binary_dilation(sub, disk(r_in)) & ~(somata[cr] > 0)
        pieces = regionprops(label(fg[cr] & ring, connectivity=2))
        n_prim = int(sum(q.area >= min_piece for q in pieces))
        per = p.perimeter * TARGET_UM
        area = p.area * TARGET_UM ** 2
        sub = (somata == p.label)
        rows.append(dict(soma_id=p.label, soma_area_um2=area, soma_perimeter_um=per,
                         soma_circularity=4 * np.pi * area / per ** 2 if per > 0 else np.nan,
                         soma_solidity=p.solidity,
                         soma_aspect_ratio=p.major_axis_length / max(p.minor_axis_length, 1e-6),
                         soma_feret_um=p.feret_diameter_max * TARGET_UM,
                         soma_snr=float(np.median(snr[sub])), n_primary_extensions=n_prim))
    return pd.DataFrame(rows)


# ── extensions ──────────────────────────────────────────────────────────────

def extension_skeleton(fg, somata):
    ext = fg & ~binary_dilation(somata > 0, disk(int(round(1.0 / TARGET_UM))))
    ext = remove_small_objects(ext, int(MIN_OBJ_UM2 / TARGET_UM ** 2))
    skel = skeletonize(ext)
    # prune short end branches once
    from skan import Skeleton, summarize
    if skel.sum() < 10:
        return ext, skel
    sk = Skeleton(skel)
    df = summarize(sk, separator='_')
    spur = (df['branch_type'] == 1) & (df['branch_distance'] * TARGET_UM < SPUR_UM)
    for i in df.index[spur]:
        coords = sk.path_coordinates(i)
        # keep the junction pixel (degree > 1) so the parent branch stays connected
        deg = sk.degrees[sk.path(i)]
        c = coords[deg <= 2]
        skel[c[:, 0].astype(int), c[:, 1].astype(int)] = False
    return ext, skeletonize(skel)


def _fwhm(profile):
    """Width at half maximum around the profile centre (in samples); NaN if it never drops to half on a side."""
    n = len(profile); c = n // 2
    w = max(1, n // 12)
    pk_i = c - w + int(np.argmax(profile[c - w:c + w + 1]))
    peak = profile[pk_i]
    base = np.percentile(profile, 10)
    if peak - base <= 0:
        return np.nan
    half = base + (peak - base) / 2
    left = np.where(profile[:pk_i] < half)[0]
    right = np.where(profile[pk_i:] < half)[0]
    if len(left) == 0 or len(right) == 0:
        return np.nan
    l, r = left[-1], pk_i + right[0]
    # linear interpolation of both crossings
    lf = l + (half - profile[l]) / (profile[l + 1] - profile[l] + 1e-9)
    rf = (r - 1) + (profile[r - 1] - half) / (profile[r - 1] - profile[r] + 1e-9)
    return rf - lf


def branch_table(skel, img_s, samples=None, density=None):
    """One row per skeleton branch: length, tortuosity, FWHM diameter (median, CV) along it."""
    from skan import Skeleton, summarize
    if skel.sum() < 10:
        return pd.DataFrame()
    sk = Skeleton(skel)
    df = summarize(sk, separator='_')
    step = max(1, int(round(SAMPLE_UM / TARGET_UM)))
    half = int(round(PROFILE_UM / TARGET_UM))
    t = np.arange(-half, half + 0.5, 0.5)
    rows = []
    for i in df.index:
        L = df.at[i, 'branch_distance'] * TARGET_UM
        e = df.at[i, 'euclidean_distance'] * TARGET_UM
        row = dict(branch_id=i, branch_type=int(df.at[i, 'branch_type']), length_um=L,
                   tortuosity=L / e if e > 0 else np.nan, diam_median_um=np.nan, diam_cv=np.nan)
        coords = sk.path_coordinates(i)
        if density is not None:
            row['local_density'] = float(np.median(density[coords[:, 0].astype(int), coords[:, 1].astype(int)]))
        if L >= MIN_SEG_UM and len(coords) > 2 * step:
            ws, ks, local, peaks = [], [], [], []
            for k in range(step, len(coords) - step, step):
                d = coords[min(k + step, len(coords) - 1)] - coords[max(k - step, 0)]
                nrm = np.hypot(*d)
                if nrm == 0:
                    continue
                nvec = np.array([-d[1], d[0]]) / nrm
                pts = coords[k][:, None] + nvec[:, None] * t[None, :]
                prof = ndimage.map_coordinates(img_s, pts, order=1, mode='nearest')
                w = _fwhm(prof)
                if np.isfinite(w):
                    ws.append(w * 0.5 * TARGET_UM); ks.append(k)
                    c_ = len(prof) // 2; w_ = max(1, len(prof) // 12)
                    peaks.append(float(prof[c_ - w_:c_ + w_ + 1].max()))
                    dv = density[int(coords[k][0]), int(coords[k][1])] if density is not None else np.nan
                    local.append([coords[k][0], coords[k][1], nvec[0], nvec[1], ws[-1], dv, 0.0])
            if len(ws) >= 3:
                ws = np.array(ws)
                med = float(np.median(ws))
                row['diam_median_um'] = med
                row['diam_cv'] = float(ws.std() / ws.mean())
                row['n_diam_samples'] = len(ws)
                # swellings: runs of consecutive samples much wider than the branch median
                pos = np.array(ks) * L / max(len(coords) - 1, 1)
                ok = (pos >= SWELL_END_UM) & (pos <= L - SWELL_END_UM)
                pk = np.array(peaks)
                nb = np.maximum(np.r_[ws[1:], -np.inf], np.r_[-np.inf, ws[:-1]])   # widest neighbour
                big = (ok & (ws > SWELL_FACTOR * med) & (ws - med >= SWELL_MIN_UM)
                       & (pk > SWELL_BRIGHT * np.median(pk)))
                # keep runs that peak locally: the run's widest sample beats the samples just outside it
                lbl_, n_runs = ndimage.label(big)
                for rid in range(1, n_runs + 1):
                    idx = np.where(lbl_ == rid)[0]
                    outside = [j for j in (idx[0] - 1, idx[-1] + 1) if 0 <= j < len(ws)]
                    if outside and ws[idx].max() <= max(ws[j] for j in outside):
                        big[idx] = False
                n_sw = int(np.sum(big[1:] & ~big[:-1]) + big[0]) if len(big) else 0
                row['n_swellings'] = n_sw
                row['swelling_tested_um'] = float(max(L - 2 * SWELL_END_UM, 0))
                for j in np.where(big)[0]:
                    local[j][6] = 1.0
            if samples is not None:
                samples.extend(tuple(v) for v in local)
        rows.append(row)
    return pd.DataFrame(rows)


# ── isolated neurons: dendrite length and area per neuron (Flamier et al. 2018, Fig 2B) ──

BORDER_UM = 2.0   # a neuron whose dendrites come this close to the image edge is flagged (length is a lower bound)

RING_FG_THR = 0.085  # AF594 mask fraction in the 1-20 µm ring at or above which a cell body counts as clustered
# calibrated by eye (2026-11-03): clearly isolated neurons 0.05-0.08, cells inside graft cores mostly 0.09-0.29

def soma_surroundings(somata, fg, r_in_um=1.0, r_out_um=20.0, r_near_um=10.0):
    """Per cell body: AF594 mask fraction in a 1-20 µm ring (own body excluded by construction) and
    whether another cell body lies within 10 µm."""
    r1, r2, rn = int(r_in_um / TARGET_UM), int(r_out_um / TARGET_UM), int(r_near_um / TARGET_UM)
    fr, near = [], []
    for q in regionprops(somata):
        y0, x0, y1, x1 = q.bbox
        cr = (slice(max(y0 - r2, 0), y1 + r2), slice(max(x0 - r2, 0), x1 + r2))
        d = ndimage.distance_transform_edt(somata[cr] != q.label)
        ring = (d > r1) & (d <= r2)
        fr.append(float(fg[cr][ring].mean()) if ring.any() else np.nan)
        others = somata[cr][(d > 0) & (d <= rn)]
        near.append(bool(((others > 0) & (others != q.label)).any()))
    return fr, near


SEED_GROW_UM = 1.5   # cell bodies are grown by this much before assigning extensions
GAP_UM       = 1.0   # mask gaps up to this size (dim stretches of one dendrite) are bridged for assignment

def territories(somata, fg):
    """Each AF594 mask pixel goes to the nearest cell body along the mask (watershed on a flat image).
    The search runs on the mask dilated by GAP_UM/2, so a dendrite interrupted by a short dim stretch stays
    with its neuron; the result is then cut back to the mask."""
    from skimage.segmentation import watershed, expand_labels
    seeds = expand_labels(somata, int(round(SEED_GROW_UM / TARGET_UM)))
    seeds = np.where(fg | (somata > 0), seeds, 0)
    bridge = binary_dilation(fg, disk(max(1, int(round(GAP_UM / 2 / TARGET_UM)))))
    terr = watershed(np.zeros(fg.shape, np.float32), markers=seeds, mask=bridge | (seeds > 0), connectivity=2)
    return np.where(fg | (somata > 0), terr, 0)


def neuron_table(somata, fg, skel, S):
    """Assign extensions to cell bodies and measure each neuron's dendrites.

    Territory: watershed from the cell bodies through the AF594 mask (flat image -> each mask pixel goes
    to the nearest cell body along the mask). Per neuron:
      dend_length_um        skeleton length in the territory (cell body excluded)
      dend_surface_area_um2 pi x mean FWHM diameter x length (what IMARIS calls dendrite area)
      dend_field_area_um2   convex hull of cell body + dendrites (area the arbor covers)
      n_branch_points, n_tips, branch_points_per_100um
      touches_border        dendrites reach the image edge -> length/area are lower bounds
    """
    from skimage.morphology import convex_hull_image
    if somata.max() == 0:
        return pd.DataFrame()
    terr = territories(somata, fg)
    nb = ndimage.convolve(skel.astype(np.uint8), np.ones((3, 3), np.uint8), mode='constant') - 1
    tips = skel & (nb == 1)
    junc = skel & (nb >= 3)
    b = int(round(BORDER_UM / TARGET_UM))
    edge = np.zeros(fg.shape, bool); edge[:b] = edge[-b:] = True; edge[:, :b] = edge[:, -b:] = True
    s_lab = terr[S[:, 0].astype(int), S[:, 1].astype(int)] if len(S) else np.array([])
    rows = []
    for p in regionprops(terr):
        r0, c0, r1, c1 = p.bbox
        cr = (slice(r0, r1), slice(c0, c1))
        own = terr[cr] == p.label
        sk = skel[cr] & own
        L = float(sk.sum() * TARGET_UM)
        w = S[s_lab == p.label, 4] if len(S) else np.array([])
        mean_d = float(w.mean()) if len(w) >= 3 else np.nan
        pts = sk | (somata[cr] == p.label)
        hull = convex_hull_image(pts) if pts.sum() >= 3 else pts
        n_j = int(label(junc[cr] & own, connectivity=2).max())
        rows.append(dict(soma_id=p.label, dend_length_um=L,
                         dend_mean_diam_um=mean_d,
                         dend_surface_area_um2=np.pi * mean_d * L if np.isfinite(mean_d) else np.nan,
                         dend_field_area_um2=float(hull.sum() * TARGET_UM ** 2),
                         n_branch_points=n_j, n_tips=int((tips[cr] & own).sum()),
                         branch_points_per_100um=100 * n_j / L if L >= 20 else np.nan,
                         touches_border=bool((sk & edge[cr]).any())))
    return pd.DataFrame(rows)


# ── one image ───────────────────────────────────────────────────────────────

def analyze_image(path):
    """Returns (image_row, branches_df, somata_df, arrays-for-plotting)."""
    neu, nuc = load_neuron(path)
    snr = snr_image(neu)
    fg = foreground_mask(snr, neu)
    somata = soma_mask(neu, snr, nuc)
    ext, skel = extension_skeleton(fg, somata)
    density = ndimage.uniform_filter(fg.astype(np.float32), size=int(DENSITY_UM / TARGET_UM))
    samples = []
    br = branch_table(skel, snr, samples, density)
    so = soma_table(somata, snr, fg)
    S = np.array(samples).reshape(-1, 7)
    if len(so) and len(S):
        # proximal (primary) extension diameter: width samples 1-5 µm from each cell body
        prox = []
        props = {q.label: q for q in regionprops(somata)}
        m = int(6 / TARGET_UM)
        for sid in so['soma_id']:
            r0, c0, r1, c1 = props[sid].bbox
            R0, C0 = max(r0 - m, 0), max(c0 - m, 0)
            crop = somata[R0:r1 + m, C0:c1 + m] != sid
            d = ndimage.distance_transform_edt(crop) * TARGET_UM
            rr, cc = S[:, 0].astype(int) - R0, S[:, 1].astype(int) - C0
            inside = (rr >= 0) & (cc >= 0) & (rr < crop.shape[0]) & (cc < crop.shape[1])
            dd = np.full(len(S), np.inf)
            dd[inside] = d[rr[inside], cc[inside]]
            w = S[(dd >= 1) & (dd <= 5), 4]
            prox.append(float(np.median(w)) if len(w) >= 3 else np.nan)
        so['primary_diam_um'] = prox
    if len(so):
        ring_fg, near = soma_surroundings(somata, fg)
        so['ring_fg_frac'] = ring_fg
        so['other_soma_within_10um'] = near
        so['in_dense_cluster'] = [bool(f >= RING_FG_THR or n) for f, n in zip(ring_fg, near)]
        nt = neuron_table(somata, fg, skel, S)
        if len(nt):
            so = so.merge(nt, on='soma_id', how='left')
    area_mm2 = neu.size * TARGET_UM ** 2 / 1e6
    total_len = float(skel.sum() * TARGET_UM)
    n_junc = 0
    if len(br):
        from skan import Skeleton
        n_junc = int((Skeleton(skel).degrees > 2).sum()) if skel.sum() >= 10 else 0
    good = br[br['length_um'] >= MIN_SEG_UM] if len(br) else br
    # isolated extensions only (outside dense graft cores)
    iso_px = skel & (density < DENSE_THR)
    iso_len = float(iso_px.sum() * TARGET_UM)
    iso_junc = 0
    if skel.sum() >= 10:
        from skan import Skeleton
        sk_ = Skeleton(skel)
        jc = sk_.coordinates[sk_.degrees > 2].astype(int)
        iso_junc = int((density[jc[:, 0], jc[:, 1]] < DENSE_THR).sum()) if len(jc) else 0
    Si = S[S[:, 5] < DENSE_THR] if len(S) else S
    gi = good[good['local_density'] < DENSE_THR] if len(good) else good
    def _swell(g):
        g = g.dropna(subset=['n_swellings']) if len(g) and 'n_swellings' in g else g.iloc[0:0]
        L = g['swelling_tested_um'].sum() if len(g) else 0
        return (100 * g['n_swellings'].sum() / L if L >= 50 else np.nan,
                100 * (g['n_swellings'] > 0).mean() if len(g) >= 5 else np.nan)
    sw_all, pct_all = _swell(good)
    sw_iso, pct_iso = _swell(gi)
    iso = dict(dense_frac=1 - iso_len / total_len if total_len > 0 else np.nan,
               swellings_per_100um=sw_all, pct_branches_swollen=pct_all,
               iso_swellings_per_100um=sw_iso, iso_pct_branches_swollen=pct_iso,
               iso_ext_length_um=iso_len,
               iso_diam_median_um=float(np.median(Si[:, 4])) if len(Si) >= 10 else np.nan,
               iso_diam_p90_um=float(np.percentile(Si[:, 4], 90)) if len(Si) >= 10 else np.nan,
               iso_diam_cv_median=float(gi['diam_cv'].median()) if len(gi) >= 3 else np.nan,
               iso_branch_points_per_100um=100 * iso_junc / iso_len if iso_len >= 50 else np.nan,
               iso_branch_length_median_um=float(gi['length_um'].median()) if len(gi) >= 3 else np.nan,
               iso_tortuosity_median=float(gi['tortuosity'].median()) if len(gi) >= 3 else np.nan)
    row = dict(field_area_mm2=area_mm2, fg_frac=float(fg.mean()),
               n_somata=int(somata.max()),
               ext_length_um=total_len,
               ext_density_um_per_mm2=total_len / area_mm2,
               branch_points_per_100um=100 * n_junc / total_len if total_len > 0 else np.nan,
               n_branches=len(good),
               branch_length_median_um=float(good['length_um'].median()) if len(good) else np.nan,
               tortuosity_median=float(good['tortuosity'].median()) if len(good) else np.nan,
               diam_median_um=float(np.median(S[:, 4])) if len(S) else np.nan,
               diam_p90_um=float(np.percentile(S[:, 4], 90)) if len(S) else np.nan,
               diam_cv_median=float(good['diam_cv'].median()) if len(good) else np.nan,
               ext_snr_median=float(np.median(snr[skel])) if skel.any() else np.nan, **iso)
    if len(so):
        for k in ['soma_area_um2', 'soma_circularity', 'soma_solidity', 'soma_aspect_ratio',
                  'soma_feret_um', 'n_primary_extensions', 'primary_diam_um']:
            row[k] = float(so[k].median()) if k in so else np.nan
    terr = territories(somata, fg) if somata.max() else np.zeros(fg.shape, np.int32)
    arrays = dict(neuron=neu, snr=snr, fg=fg, somata=somata, skel=skel, diam_samples=S, territory=terr,
                  dense=density >= DENSE_THR)
    return row, br, so, arrays


def analyze_all(files, keep_arrays=()):
    """files: output of list_images. Returns (images_df, branches_df, somata_df, kept arrays by filename key)."""
    img_rows, brs, sos, kept = [], [], [], {}
    for _, f in files.iterrows():
        row, br, so, arr = analyze_image(f['path'])
        key = f"{f['mouse']}/{f['condition']}/{f['filename']}"
        meta = dict(mouse=f['mouse'], condition=f['condition'], filename=f['filename'])
        img_rows.append({**meta, 'n_z': f['n_z'], 'px_nm': f['px_nm'], 'objective': f['objective'], **row})
        if len(br): brs.append(br.assign(**meta))
        if len(so): sos.append(so.assign(**meta))
        if key in keep_arrays:
            kept[key] = arr
    cat = lambda L: pd.concat(L, ignore_index=True) if L else pd.DataFrame()
    return pd.DataFrame(img_rows), cat(brs), cat(sos), kept


# ── statistics and plots ────────────────────────────────────────────────────

GROUP_ORDER = ['ID06902 scr', 'ID06902 bmi1', 'ID06914 scr']
GROUP_COLORS = {'ID06902 scr': '#4C72B0', 'ID06902 bmi1': '#DD8452', 'ID06914 scr': '#8FAADC'}

NEURON_METRICS = [
    ('dend_length_um',          'Dendrite length per neuron (µm)'),
    ('dend_surface_area_um2',   'Dendrite surface area per neuron (µm²)'),
    ('dend_field_area_um2',     'Dendritic field area per neuron (µm²)'),
    ('dend_mean_diam_um',       'Mean dendrite diameter per neuron (µm)'),
    ('n_branch_points',         'Branch points per neuron'),
    ('n_tips',                  'Dendrite tips per neuron'),
    ('branch_points_per_100um', 'Branch points per 100 µm dendrite'),
]
SWELL_METRICS = [
    ('iso_pct_branches_swollen', 'Isolated: % branches with >= 1 swelling'),
    ('iso_swellings_per_100um',  'Isolated: swellings per 100 µm'),
    ('pct_branches_swollen',     'All: % branches with >= 1 swelling'),
    ('swellings_per_100um',      'All: swellings per 100 µm'),
]
ISO_METRICS = [
    ('dense_frac',                  'Share of extensions in dense clusters'),
    ('iso_diam_median_um',          'Isolated: diameter, median (µm)'),
    ('iso_diam_p90_um',             'Isolated: diameter, 90th pct (µm)'),
    ('iso_diam_cv_median',          'Isolated: diameter CV along branch'),
    ('iso_branch_points_per_100um', 'Isolated: branch points per 100 µm'),
    ('iso_branch_length_median_um', 'Isolated: branch length, median (µm)'),
    ('iso_tortuosity_median',       'Isolated: tortuosity, median'),
]
IMAGE_METRICS = [
    ('diam_median_um',          'Extension diameter, median (µm)'),
    ('diam_p90_um',             'Extension diameter, 90th pct (µm)'),
    ('diam_cv_median',          'Diameter variation along branch (CV)'),
    ('ext_density_um_per_mm2',  'Extension length per field (µm / mm²)'),
    ('branch_points_per_100um', 'Branch points per 100 µm'),
    ('branch_length_median_um', 'Branch length, median (µm)'),
    ('tortuosity_median',       'Tortuosity, median'),
    ('ext_snr_median',          'AF594 signal on extensions (SNR)'),
]
SOMA_METRICS = [
    ('soma_area_um2',        'Soma area (µm²)'),
    ('soma_feret_um',        'Soma max length (µm)'),
    ('soma_circularity',     'Soma circularity'),
    ('soma_aspect_ratio',    'Soma aspect ratio'),
    ('n_primary_extensions', 'Primary extensions per soma'),
    ('primary_diam_um',      'Primary extension diameter (µm)'),
]


def add_group(df):
    return df.assign(group=df['mouse'] + ' ' + df['condition'])


def stats_table(df, metrics):
    """Median per group + Mann-Whitney p for bmi1 vs scr (same mouse ID06902) and bmi1 vs all scr."""
    from scipy.stats import mannwhitneyu
    df = add_group(df)
    rows = []
    for key, name in metrics:
        if key not in df:
            continue
        g = {k: df.loc[df['group'] == k, key].dropna() for k in GROUP_ORDER}
        allscr = df.loc[df['condition'] == 'scr', key].dropna()
        bmi = g['ID06902 bmi1']
        p_same = mannwhitneyu(bmi, g['ID06902 scr']).pvalue if len(bmi) > 1 and len(g['ID06902 scr']) > 1 else np.nan
        p_all = mannwhitneyu(bmi, allscr).pvalue if len(bmi) > 1 and len(allscr) > 1 else np.nan
        p_mouse = (mannwhitneyu(g['ID06902 scr'], g['ID06914 scr']).pvalue
                   if len(g['ID06902 scr']) > 1 and len(g['ID06914 scr']) > 1 else np.nan)
        rows.append({'metric': name, **{k: v.median() for k, v in g.items()},
                     **{f'n {k}': len(v) for k, v in g.items()},
                     'bmi1/scr (ID06902)': bmi.median() / g['ID06902 scr'].median() if len(bmi) and len(g['ID06902 scr']) else np.nan,
                     'p bmi1 vs scr, ID06902': p_same, 'p bmi1 vs all scr': p_all,
                     'p scr ID06902 vs ID06914': p_mouse})
    return pd.DataFrame(rows).set_index('metric')


def plot_strips(df, metrics, save=None, ncols=4, title=''):
    """One panel per metric: points per image (or per soma), median bar, p for bmi1 vs scr in ID06902."""
    import matplotlib.pyplot as plt
    from scipy.stats import mannwhitneyu
    df = add_group(df)
    metrics = [m for m in metrics if m[0] in df]
    nrows = int(np.ceil(len(metrics) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.8 * nrows), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, (key, name) in zip(axes.ravel(), metrics):
        for xi, gname in enumerate(GROUP_ORDER):
            v = df.loc[df['group'] == gname, key].dropna().to_numpy()
            ax.scatter(xi + rng.uniform(-0.15, 0.15, len(v)), v, s=28, alpha=0.8,
                       color=GROUP_COLORS[gname], zorder=3)
            if len(v):
                ax.hlines(np.median(v), xi - 0.28, xi + 0.28, color='k', lw=2, zorder=4)
        a = df.loc[df['group'] == 'ID06902 bmi1', key].dropna()
        b = df.loc[df['group'] == 'ID06902 scr', key].dropna()
        p = mannwhitneyu(a, b).pvalue if len(a) > 1 and len(b) > 1 else np.nan
        ax.set_title(f'{name}\nbmi1 vs scr (ID06902): p={p:.3f}', fontsize=9.5)
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels([g.replace(' ', '\n') for g in GROUP_ORDER], fontsize=9)
        ax.grid(axis='y', alpha=0.3, ls='--')
    for ax in axes.ravel()[len(metrics):]:
        ax.axis('off')
    if title:
        fig.suptitle(title, fontsize=12)
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')
    plt.show()


def plot_branch_distributions(branches, save=None, isolated_only=False):
    """Per-branch diameter, diameter CV, length and tortuosity as cumulative distributions per group."""
    import matplotlib.pyplot as plt
    br = add_group(branches)
    br = br[br['length_um'] >= MIN_SEG_UM]
    if isolated_only:
        br = br[br['local_density'] < DENSE_THR]
    cols = [('diam_median_um', 'Branch diameter (µm)', (0, 3)),
            ('diam_cv', 'Diameter CV along branch', (0, 1)),
            ('length_um', 'Branch length (µm)', (0, 40)),
            ('tortuosity', 'Tortuosity', (1, 2))]
    fig, axes = plt.subplots(1, 4, figsize=(19, 4))
    for ax, (key, name, lim) in zip(axes, cols):
        for gname in GROUP_ORDER:
            v = np.sort(br.loc[br['group'] == gname, key].dropna().to_numpy())
            if len(v):
                ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=GROUP_COLORS[gname], lw=2,
                        label=f'{gname} (n={len(v)})')
        ax.set_xlim(*lim); ax.set_xlabel(name); ax.set_ylabel('cumulative fraction')
        ax.grid(alpha=0.3, ls='--')
    axes[0].legend(fontsize=8)
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=150, bbox_inches='tight')
    plt.show()


def plot_examples(examples, save=None):
    """examples: list of (title, arrays). Raw AF594 | segmentation overlay with width samples coloured by diameter."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    n = len(examples)
    fig, axes = plt.subplots(n, 2, figsize=(12, 6 * n), squeeze=False)
    norm = Normalize(0.3, 2.0)
    for (title, a), (ax0, ax1) in zip(examples, axes):
        v = np.percentile(a['neuron'], 99.7)
        ax0.imshow(a['neuron'], cmap='gray', vmax=v); ax0.set_title(title, fontsize=10)
        ax1.imshow(a['neuron'], cmap='gray', vmax=v)
        ax1.imshow(np.ma.masked_where(~a['fg'], a['fg']), cmap='Greens', alpha=0.25, vmin=0, vmax=1.5)
        if a['somata'].any():
            ax1.contour(a['somata'] > 0, [0.5], colors='red', linewidths=1.5)
        if 'dense' in a and a['dense'].any():
            ax1.contour(a['dense'], [0.5], colors='cyan', linewidths=1.2, linestyles='--')
        S = a['diam_samples']
        if len(S):
            sc = ax1.scatter(S[:, 1], S[:, 0], c=S[:, 4], s=6, cmap='plasma', norm=norm)
            if S.shape[1] > 6 and (S[:, 6] > 0).any():
                sw = S[S[:, 6] > 0]
                ax1.scatter(sw[:, 1], sw[:, 0], s=70, facecolors='none', edgecolors='white', linewidths=1.2)
        ax1.set_title('green: AF594 mask   red: cell bodies   cyan dashed: dense cluster\ndots: extension diameter   white circles: swellings', fontsize=9)
        for ax in (ax0, ax1):
            ax.axis('off')
    if len(examples):
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='plasma'), ax=axes[:, 1].tolist(),
                          shrink=0.4, label='diameter (µm)')
    if save:
        fig.savefig(save, dpi=110, bbox_inches='tight')
    plt.show()


def plot_diameter_zoom(a, center, half=200, save=None, title=''):
    """Zoom with the FWHM width measured at each sample drawn as a line across the extension."""
    import matplotlib.pyplot as plt
    r0, c0 = center
    sy = slice(max(r0 - half, 0), r0 + half); sx = slice(max(c0 - half, 0), c0 + half)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(a['neuron'][sy, sx], cmap='gray', vmax=np.percentile(a['neuron'], 99.8))
    ax.contour((a['somata'] > 0)[sy, sx], [0.5], colors='red')
    for r, c, nr, nc, w, *rest in a["diam_samples"]:
        if sy.start <= r < sy.stop and sx.start <= c < sx.stop:
            h = w / TARGET_UM / 2
            x, y = c - sx.start, r - sy.start
            swollen = len(rest) > 1 and rest[1] > 0
            ax.plot([x - nc * h, x + nc * h], [y - nr * h, y + nr * h], color='red' if swollen else 'yellow',
                    lw=2.0 if swollen else 1.2)
            ax.text(x + nc * h + 2, y + nr * h, f'{w:.1f}', color='yellow', fontsize=7)
    ax.set_title(title or 'FWHM width samples (µm); red = swelling', fontsize=10); ax.axis('off')
    if save:
        fig.savefig(save, dpi=130, bbox_inches='tight')
    plt.show()


def plot_channels(path, save=None):
    """All channels of one CZI (max-Z), with their names."""
    import matplotlib.pyplot as plt
    meta = czi_meta(path)
    with pyczi.open_czi(path) as cz:
        ims = [_mip(cz, c, meta['n_z']) for c in range(len(meta['channels']))]
    fig, axes = plt.subplots(1, len(ims), figsize=(4 * len(ims), 4), squeeze=False)
    for ax, im, nm in zip(axes[0], ims, meta['channels']):
        ax.imshow(im, cmap='gray', vmax=max(np.percentile(im, 99.7), 1)); ax.set_title(nm, fontsize=10); ax.axis('off')
    fig.suptitle(os.path.relpath(path, os.path.dirname(os.path.dirname(os.path.dirname(path)))), fontsize=10)
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=100, bbox_inches='tight')
    plt.show()


def plot_neuron_gallery(items, save=None, pad_um=5):
    """items: list of (title, arrays, soma_id). Crop around each neuron's territory: raw | assigned dendrites."""
    import matplotlib.pyplot as plt
    n = len(items)
    if n == 0:
        return
    fig, axes = plt.subplots(n, 2, figsize=(9, 4.5 * n), squeeze=False)
    pad = int(pad_um / TARGET_UM)
    for (title, a, sid), (ax0, ax1) in zip(items, axes):
        own = a['territory'] == sid
        rr, cc = np.where(own)
        sy = slice(max(rr.min() - pad, 0), rr.max() + pad); sx = slice(max(cc.min() - pad, 0), cc.max() + pad)
        img = a['neuron'][sy, sx]; v = np.percentile(a['neuron'], 99.7)
        ax0.imshow(img, cmap='gray', vmax=v); ax0.set_title(title, fontsize=9)
        ax1.imshow(img, cmap='gray', vmax=v)
        sk = a['skel'][sy, sx] & own[sy, sx]
        other = a['skel'][sy, sx] & ~own[sy, sx]
        ax1.imshow(np.ma.masked_where(~other, other), cmap='Greys', alpha=0.6, interpolation='nearest')
        ax1.imshow(np.ma.masked_where(~sk, sk), cmap='autumn', interpolation='nearest')
        ax1.contour((a['somata'] == sid)[sy, sx], [0.5], colors='cyan', linewidths=1.5)
        ax1.set_title('orange: dendrites assigned to this neuron   cyan: cell body', fontsize=8)
        for ax in (ax0, ax1):
            ax.axis('off')
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=110, bbox_inches='tight')
    plt.show()
