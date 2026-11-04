"""Figures and summary statistics for the image-analysis protocol (paper_protocol/).

Figure 1  axon breaks (in vitro, FAST U-Net++ mask)      -> figures/Figure1_axon_breaks.{png,pdf}
Figure 2  tubulin puncta (in vitro, Spotiflow)           -> figures/Figure2_puncta.{png,pdf}
Figure 3  grafted-neuron morphology (in vivo, 2WPT, CZI) -> figures/Figure3_morphology.{png,pdf}
All numbers quoted in the protocol text           -> figures/stats_summary.csv

Run with the `axon` environment:  ~/axon/bin/python make_figures.py
"""
import os, re, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
from scipy.stats import mannwhitneyu, ttest_1samp
from scipy import ndimage
from skimage.measure import label
from pylibCZIrw import czi as pyczi

BASE = '/Users/ctuna/Desktop/axon_codes'
OUT = os.path.join(BASE, 'paper_protocol')
FIG = os.path.join(OUT, 'figures')
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 300})

COL = {'scr': '#4d4d4d', 'Ni': '#6f9fc8', 'bmi': '#c0392b'}
LAB = {'scr': 'SCR', 'Ni': 'NI', 'bmi': 'BMI1'}
EXPS = [  # (label, FAST break csv, spot csv, image directories by condition)
    ('Exp. 1', 'N1/csv/axon_break_metrics_FAST_N1.csv', 'N1/csv/spot_counts_N1.csv',
     {'bmi': 'N1/x63 bmi tubulin iii_data', 'Ni': 'N1/x63 Ni tubulin_data', 'scr': 'N1/x63 scr tubIII ID8221225_data'}),
    ('Exp. 2', 'N2/csv/axon_break_metrics_FAST_N2_without_airyscan.csv', 'N2/csv/spot_counts_N2_without_airyscan.csv',
     {'bmi': 'N2/BMI!/ID090426 j14 bmi1_data', 'Ni': 'N2/NI/id090426  j14 NI _data', 'scr': 'N2/SCR/SCR ID090426 J14 _data'}),
    ('Exp. 3', 'N3/csv/axon_break_metrics_FAST_N3.csv', 'N3/csv/spot_counts_N3.csv',
     {'bmi': 'N3/BMI!/id100426 bmi1 j14_data', 'Ni': 'N3/NI/ni n=3 _data', 'scr': 'N3/SCR/scr n=3_data'}),
    ('Exp. 4', 'N4/csv/axon_break_metrics_FAST_N4.csv', 'N4/csv/spot_counts_N4.csv',
     {'bmi': 'N4/bmi1 n=4/bmi1 n4_data', 'scr': 'N4/SCR n=4_data'}),
]
PX_FAST_UM = 0.08768704530423278          # break analysis runs at this pixel size
# Example images (Exp. 3) and zoom windows (top-left row, col in image pixels), chosen by eye from fields
# whose value is close to the group median: SCR shows long continuous axons, BMI1 fragmented axons.
FIG1_EXAMPLES = {'scr': ('New-12.czi', (1311, 912)), 'bmi': ('New-10.czi', (912, 513))}       # 87.7 nm px, 15 µm
FIG2_EXAMPLES = {'scr': ('New-18.czi', (0, 630)), 'bmi': ('New-10.czi', (980, 1050))}      # 70.6 nm px, 10 µm
STATS = []                                 # rows for stats_summary.csv


def stars(p):
    return 'ns' if not np.isfinite(p) or p >= 0.05 else ('*' if p >= 0.01 else ('**' if p >= 0.001 else '***'))


def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    return mannwhitneyu(a, b, alternative='two-sided').pvalue if len(a) > 1 and len(b) > 1 else np.nan


def panel_letter(fig, ax, letter, dx=-0.045, dy=0.01):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=10, fontweight='bold', va='bottom', ha='left')


def bracket(ax, x1, x2, y, text, h_frac=0.03):
    lo, hi = ax.get_ylim(); h = (hi - lo) * h_frac
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.6, color='k', clip_on=False)
    ax.text((x1 + x2) / 2, y + h, text, ha='center', va='bottom', fontsize=6.5 if text == 'ns' else 8, clip_on=False)


def strip(ax, data, order, colors, labels, ylabel=None, comparisons=(), seed=0, ymax=None):
    """data: dict key -> 1D values. Dots = images (or cells), bar = mean, error = SEM."""
    rng = np.random.default_rng(seed)
    top = 0
    for i, k in enumerate(order):
        v = np.asarray(data.get(k, []), float); v = v[np.isfinite(v)]
        if not len(v):
            continue
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=7, color=colors[k], alpha=0.55,
                   linewidths=0, zorder=2)
        m = v.mean(); se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0
        ax.hlines(m, i - 0.28, i + 0.28, color='k', lw=1.0, zorder=3)
        ax.errorbar(i, m, yerr=se, color='k', lw=0.8, capsize=2.5, capthick=0.8, zorder=3)
        top = max(top, v.max())
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([labels[k] for k in order])
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_ylim(0, ymax if ymax else top * (1.15 + 0.13 * len(comparisons)))
    lo, hi = ax.get_ylim()
    for j, (a, b, p) in enumerate(comparisons):
        bracket(ax, order.index(a), order.index(b), top + (hi - lo) * (0.05 + 0.11 * j), stars(p))


def ratio_panel(ax, ratios, ylabel, title_p):
    """ratios: dict cond -> list of (exp_label, ratio to SCR). One marker shape per experiment."""
    marks = {'Exp. 1': 'o', 'Exp. 2': 's', 'Exp. 3': '^', 'Exp. 4': 'D'}
    order = [k for k in ['scr', 'Ni', 'bmi'] if k in ratios]
    for i, k in enumerate(order):
        vals = np.array([r for _, r in ratios[k]])
        ax.bar(i, vals.mean(), width=0.6, color=COL[k], alpha=0.25 if k != 'scr' else 0.15,
               edgecolor=COL[k], linewidth=0.8)
        if k != 'scr':
            ax.errorbar(i, vals.mean(), yerr=vals.std(ddof=1) / np.sqrt(len(vals)), color='k', lw=0.8,
                        capsize=2.5, capthick=0.8)
        for e, r in ratios[k]:
            ax.scatter(i + {'Exp. 1': -0.15, 'Exp. 2': -0.05, 'Exp. 3': 0.05, 'Exp. 4': 0.15}[e], r,
                       marker=marks[e], s=12, facecolor='white', edgecolor=COL[k], linewidths=0.8, zorder=3)
    ax.axhline(1, color='k', lw=0.5, ls=':')
    ax.set_xticks(range(len(order))); ax.set_xticklabels([LAB[k] for k in order])
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylabel(ylabel)
    top = max(r for k in order for _, r in ratios[k])
    ax.set_ylim(0, top * 1.3)
    bracket(ax, order.index('scr'), order.index('bmi'), top * 1.06, stars(title_p))
    handles = [plt.Line2D([], [], marker=m, ls='', mfc='white', mec='k', mew=0.8, ms=3.5, label=e)
               for e, m in marks.items()]
    ax.legend(handles=handles, fontsize=5.5, loc='upper left', bbox_to_anchor=(0.0, 0.8), handletextpad=0.1,
              borderaxespad=0.1, labelspacing=0.2)


def scalebar(ax, um, px_um, img_shape, label=True, color='white', frac_y=0.07, lw=2.0):
    h, w = img_shape
    L = um / px_um
    x1 = w * 0.95; x0 = x1 - L; y = h * (1 - frac_y)
    ax.plot([x0, x1], [y, y], color=color, lw=lw, solid_capstyle='butt')
    if label:
        ax.text((x0 + x1) / 2, y - h * 0.025, f'{um:g} µm', color=color, ha='center', va='bottom', fontsize=6)


def show_gray(ax, img, lo=0.5, hi=99.7, ref=None):
    a, b = np.percentile(img if ref is None else ref, [lo, hi])
    ax.imshow(np.clip((img - a) / (b - a + 1e-9), 0, 1), cmap='gray', interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def image_row(fig, y, h_in, items, fig_h_in, fig_w_in=7.2, left=0.075, gap=0.012, zoom_frac=0.6):
    """Lay out [full, zoom, full, zoom] image axes at fixed positions. Returns list of (ax_full, ax_zoom)."""
    full_w = h_in / fig_w_in; zoom_w = full_w * zoom_frac
    h = h_in / fig_h_in
    out, x = [], left
    for i in range(items):
        af = fig.add_axes([x, y, full_w, h]); x += full_w + gap
        az = fig.add_axes([x, y + h * (1 - zoom_frac), zoom_w, h * zoom_frac]); x += zoom_w + 3.2 * gap
        out.append((af, az))
    return out


def example_rows(fig, FH, reps, overlay, zoom_bar_um, col_titles, full_marks=None, side=1.95):
    """Two image rows (SCR, BMI1): full field | raw zoom | zoom with the analysis drawn on it.
    Returns the field axes and the bottom of the image block (figure fraction)."""
    left_in, gap_in, top_in, row_gap_in = 0.42, 0.10, 0.22, 0.12
    axA = []
    for i, k in enumerate(['scr', 'bmi']):
        r = reps[k]; img = r['img']; px = r['px']
        y0, x0, w = r['box']
        r['clim'] = (0.5, 99.8)
        y_in = FH - top_in - (i + 1) * side - i * row_gap_in
        axes = [fig.add_axes([(left_in + j * (side + gap_in)) / 7.2, y_in / FH, side / 7.2, side / FH])
                for j in range(3)]
        show_gray(axes[0], img, *r['clim'])
        if full_marks:
            full_marks(axes[0], r)
        axes[0].add_patch(Rectangle((x0, y0), w, w, fill=False, ec='white', lw=0.8, ls='--'))
        sub = img[y0:y0 + w, x0:x0 + w]
        show_gray(axes[1], sub, *r['clim'], ref=img)
        overlay(axes[2], r, sub, y0, x0, w)
        scalebar(axes[0], 20, px, img.shape)
        for ax in axes[1:]:
            scalebar(ax, zoom_bar_um, px, sub.shape)
        if i == 0:
            for ax, t in zip(axes, col_titles):
                ax.set_title(t, fontsize=7, pad=3)
        fig.text((left_in - 0.12) / 7.2, (y_in + side / 2) / FH, LAB[k], rotation=90, ha='center', va='center',
                 fontsize=8, fontweight='bold', color=COL[k] if k == 'bmi' else 'k')
        axA.append(axes[0])
    return axA, (FH - top_in - 2 * side - row_gap_in) / FH


def load_mip(path, channel=0):
    with pyczi.open_czi(path) as cz:
        nz = cz.total_bounding_box['Z'][1]; md = cz.raw_metadata
        img = np.stack([np.squeeze(cz.read(plane={'C': channel, 'Z': z, 'T': 0})) for z in range(nz)]).max(0)
    px = float(re.findall(r'<Distance Id="X">\s*<Value>([^<]+)</Value>', md)[0]) * 1e6
    return img.astype(np.float32), px


def record(section, metric, comparison, **kw):
    STATS.append(dict(section=section, metric=metric, comparison=comparison, **kw))


# ════════════════════════════════════════════════════════════════════════════
# Figure 1: axon breaks
# ════════════════════════════════════════════════════════════════════════════

def figure1():
    to100 = 100 / (1000 * PX_FAST_UM)          # confirmed_rate is per 1,000 skeleton px at 87.7 nm
    per_exp = []
    for lab, bcsv, _, _ in EXPS:
        d = pd.read_csv(os.path.join(BASE, 'results_with_Fast-Unet', bcsv))
        d['val'] = d['confirmed_rate'] * to100
        per_exp.append((lab, d))
        g = {k: d.loc[d.condition == k, 'val'].to_numpy() for k in ['scr', 'Ni', 'bmi']}
        p = mwu(g['bmi'], g['scr'])
        record('Axon breaks', 'confirmed breaks per 100 µm', f'{lab} BMI1 vs SCR (per image, Mann-Whitney)',
               n_bmi=len(g['bmi']), n_scr=len(g['scr']), n_ni=len(g['Ni']),
               mean_bmi=g['bmi'].mean(), mean_scr=g['scr'].mean(),
               mean_ni=g['Ni'].mean() if len(g['Ni']) else np.nan, ratio=g['bmi'].mean() / g['scr'].mean(), p=p)

    # representative images: Exp. 3, image closest to the group median
    ns = {}
    nb = __import__('nbformat').read(os.path.join(BASE, 'results_with_Fast-Unet/N3/axon_break_detection_FAST_N3.ipynb'), as_version=4)
    for c in nb.cells[1:8]:
        if c.cell_type == 'code':
            exec(c.source, ns)
    d3 = per_exp[2][1]
    dirs = EXPS[2][3]
    reps = {}
    for k in ['scr', 'bmi']:
        fn, (y0, x0) = FIG1_EXAMPLES[k]
        img = ns['_load'](os.path.join(BASE, dirs[k], fn))
        skel, binary, _ = ns['_make_skeleton'](img)
        sc = ns['filter_small_components'](skel, min_size=15).astype(bool)
        br, _ = ns['find_broken_lines'](skel, binary_mask=binary)
        conf = ns['confirm_breaks'](img, binary, sc, br)
        reps[k] = dict(fn=fn, img=img, px=PX_FAST_UM, skel=sc, breaks=[b for b, c in zip(br, conf) if c],
                       box=(y0, x0, int(15 / PX_FAST_UM)),
                       val=float(d3.loc[(d3.condition == k) & (d3.filename == fn), 'val'].iloc[0]))

    rng = np.random.default_rng(3)
    pal = plt.cm.tab20(np.linspace(0, 1, 20))[:, :3][rng.permutation(20)]

    def overlay(ax, r, sub, y0, x0, w):
        """Each continuous axon segment in its own colour; triangles mark confirmed breaks."""
        a, b = np.percentile(r['img'], r['clim'])
        g = np.clip((sub - a) / (b - a + 1e-9), 0, 1)
        lab = label(r['skel'][y0:y0 + w, x0:x0 + w], connectivity=2)
        idx = ndimage.distance_transform_edt(lab == 0, return_distances=False, return_indices=True)
        near = lab[idx[0], idx[1]]
        thick = ndimage.binary_dilation(lab > 0, iterations=1)
        rgb = np.dstack([g * 0.4] * 3)
        rgb[thick] = pal[(near[thick] - 1) % 20]
        ax.imshow(rgb, interpolation='nearest'); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        for bk in r['breaks']:
            cy = (bk['pt1'][0] + bk['pt2'][0]) / 2 - y0; cx = (bk['pt1'][1] + bk['pt2'][1]) / 2 - x0
            if 0 <= cy < w and 0 <= cx < w:
                ax.plot(cx, cy, marker='o', ms=5.5, mfc='none', mec='white', mew=0.9)

    def full_marks(ax, r):
        for bk in r['breaks']:
            cy = (bk['pt1'][0] + bk['pt2'][0]) / 2; cx = (bk['pt1'][1] + bk['pt2'][1]) / 2
            ax.add_patch(plt.Circle((cx, cy), 22, fill=False, ec='#ffd400', lw=0.5))

    FH = 6.9
    fig = plt.figure(figsize=(7.2, FH))
    axA, img_bottom = example_rows(fig, FH, reps, overlay, 5, ['Field', 'Zoom', 'Continuous axon segments'],
                                   full_marks=full_marks)
    gs = gridspec.GridSpec(1, 1, left=0.075, right=0.985, top=img_bottom - 0.075, bottom=0.06)
    gs = [None, gs[0]]

    bot = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[1], wspace=0.55, width_ratios=[1, 1, 1, 0.75, 1.15])
    axes_b = []
    ymax = max(np.nanmax(d['val']) for _, d in per_exp) * 1.32
    for j, (lab, d) in enumerate(per_exp):
        ax = fig.add_subplot(bot[0, j])
        order = [k for k in ['scr', 'Ni', 'bmi'] if (d.condition == k).any()]
        g = {k: d.loc[d.condition == k, 'val'].to_numpy() for k in order}
        strip(ax, g, order, COL, LAB, ylabel='Confirmed breaks per 100 µm axon' if j == 0 else None,
              comparisons=[('scr', 'bmi', mwu(g['bmi'], g['scr']))], seed=j, ymax=ymax)
        ax.set_title(lab, fontsize=7)
        if j:
            ax.set_yticklabels([])
        axes_b.append(ax)
    # experiment-level summary
    ratios = {'scr': [], 'Ni': [], 'bmi': []}
    for lab, d in per_exp:
        ms = d.loc[d.condition == 'scr', 'val'].mean()
        for k in ['scr', 'Ni', 'bmi']:
            if (d.condition == k).any():
                ratios[k].append((lab, d.loc[d.condition == k, 'val'].mean() / ms))
    lr = np.log2([r for _, r in ratios['bmi']])
    p_exp = ttest_1samp(lr, 0).pvalue
    record('Axon breaks', 'confirmed breaks per 100 µm, BMI1/SCR ratio of experiment means',
           'BMI1 vs SCR across N=4 experiments (one-sample t-test on log2 ratio)',
           ratio=float(2 ** lr.mean()), ratios=';'.join(f'{e}:{r:.2f}' for e, r in ratios['bmi']), p=p_exp)
    axC = fig.add_subplot(bot[0, 4])
    ratio_panel(axC, ratios, 'Breaks relative to SCR', p_exp)
    axC.set_title('N = 4 experiments', fontsize=7)

    panel_letter(fig, axA[0], 'A', dx=-0.055, dy=0.012)
    panel_letter(fig, axes_b[0], 'B', dx=-0.06)
    panel_letter(fig, axC, 'C', dx=-0.06)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, f'Figure1_axon_breaks.{ext}'))
    plt.close(fig)
    return reps


# ════════════════════════════════════════════════════════════════════════════
# Figure 2: tubulin puncta (Spotiflow)
# ════════════════════════════════════════════════════════════════════════════

def field_area_mm2(path):
    with pyczi.open_czi(path) as cz:
        bb, md = cz.total_bounding_box, cz.raw_metadata
    px = float(re.findall(r'<Distance Id="X">\s*<Value>([^<]+)</Value>', md)[0]) * 1e6
    return (bb['X'][1] - bb['X'][0]) * (bb['Y'][1] - bb['Y'][0]) * px * px / 1e6


def figure2():
    per_exp = []
    for lab, _, scsv, dirs in EXPS:
        d = pd.read_csv(os.path.join(BASE, 'results', scsv))
        d['area_mm2'] = [field_area_mm2(os.path.join(BASE, dirs[c], f)) for c, f in zip(d.condition, d.filename)]
        d['val'] = d['n_spots'] / (d['area_mm2'] * 1e6) * 1000          # puncta per 1,000 µm²
        per_exp.append((lab, d))
        g = {k: d.loc[d.condition == k, 'val'].to_numpy() for k in ['scr', 'Ni', 'bmi']}
        record('Puncta', 'puncta per 1,000 µm²', f'{lab} BMI1 vs SCR (per image, Mann-Whitney)',
               n_bmi=len(g['bmi']), n_scr=len(g['scr']), n_ni=len(g['Ni']), mean_bmi=g['bmi'].mean(),
               mean_scr=g['scr'].mean(), mean_ni=g['Ni'].mean() if len(g['Ni']) else np.nan,
               ratio=g['bmi'].mean() / g['scr'].mean(), p=mwu(g['bmi'], g['scr']))

    from spotiflow.model import Spotiflow
    model = Spotiflow.from_pretrained('general')
    d3 = per_exp[2][1]; dirs = EXPS[2][3]
    reps = {}
    for k in ['scr', 'bmi']:
        fn, (y0, x0) = FIG2_EXAMPLES[k]
        img, px = load_mip(os.path.join(BASE, dirs[k], fn))
        pts, _ = model.predict(img, verbose=False)
        reps[k] = dict(fn=fn, img=img, px=px, pts=pts, box=(y0, x0, int(10 / px)),
                       val=float(d3.loc[(d3.condition == k) & (d3.filename == fn), 'val'].iloc[0]))

    def overlay(ax, r, sub, y0, x0, w):
        show_gray(ax, sub, *r['clim'], ref=r['img'])
        P = r['pts']
        q = (P[:, 0] >= y0) & (P[:, 0] < y0 + w) & (P[:, 1] >= x0) & (P[:, 1] < x0 + w)
        ax.scatter(P[q, 1] - x0, P[q, 0] - y0, s=30, facecolors='none', edgecolors='#ffd400', linewidths=0.7)

    FH = 6.9
    fig = plt.figure(figsize=(7.2, FH))
    axA, img_bottom = example_rows(fig, FH, reps, overlay, 2, ['Field', 'Zoom', 'Detected puncta'])
    gs = gridspec.GridSpec(1, 1, left=0.075, right=0.985, top=img_bottom - 0.075, bottom=0.06)
    gs = [None, gs[0]]

    bot = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[1], wspace=0.55, width_ratios=[1, 1, 1, 0.75, 1.15])
    axes_b = []
    ymax = max(np.nanmax(d['val']) for _, d in per_exp) * 1.32
    for j, (lab, d) in enumerate(per_exp):
        ax = fig.add_subplot(bot[0, j])
        order = [k for k in ['scr', 'Ni', 'bmi'] if (d.condition == k).any()]
        g = {k: d.loc[d.condition == k, 'val'].to_numpy() for k in order}
        strip(ax, g, order, COL, LAB, ylabel='Puncta per 1,000 µm²' if j == 0 else None,
              comparisons=[('scr', 'bmi', mwu(g['bmi'], g['scr']))], seed=j, ymax=ymax)
        ax.set_title(lab, fontsize=7)
        if j:
            ax.set_yticklabels([])
        axes_b.append(ax)
    ratios = {'scr': [], 'Ni': [], 'bmi': []}
    for lab, d in per_exp:
        ms = d.loc[d.condition == 'scr', 'val'].mean()
        for k in ['scr', 'Ni', 'bmi']:
            if (d.condition == k).any():
                ratios[k].append((lab, d.loc[d.condition == k, 'val'].mean() / ms))
    lr = np.log2([r for _, r in ratios['bmi']])
    p_exp = ttest_1samp(lr, 0).pvalue
    record('Puncta', 'puncta per 1,000 µm², BMI1/SCR ratio of experiment means',
           'BMI1 vs SCR across N=4 experiments (one-sample t-test on log2 ratio)',
           ratio=float(2 ** lr.mean()), ratios=';'.join(f'{e}:{r:.2f}' for e, r in ratios['bmi']), p=p_exp)
    axC = fig.add_subplot(bot[0, 4])
    ratio_panel(axC, ratios, 'Puncta relative to SCR', p_exp)
    axC.set_title('N = 4 experiments', fontsize=7)
    panel_letter(fig, axA[0], 'A', dx=-0.055, dy=0.012)
    panel_letter(fig, axes_b[0], 'B', dx=-0.06)
    panel_letter(fig, axC, 'C', dx=-0.06)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, f'Figure2_puncta.{ext}'))
    plt.close(fig)
    return reps


# ════════════════════════════════════════════════════════════════════════════
# Figure 3: grafted-neuron morphology (2WPT, regular CZI)
# ════════════════════════════════════════════════════════════════════════════

def figure3():
    sys.path.insert(0, os.path.join(BASE, 'results_2WPT'))
    import morphology_utils as mu
    im = pd.read_csv(os.path.join(BASE, 'results_2WPT/czi/csv/morphology_images_czi.csv'))
    im['grp'] = im['mouse'] + ' ' + im['condition']
    G = ['ID06902 scr', 'ID06902 bmi1', 'ID06914 scr']
    GL = {'ID06902 scr': 'SCR\nmouse 1', 'ID06902 bmi1': 'BMI1\nmouse 1', 'ID06914 scr': 'SCR\nmouse 2'}
    GC = {'ID06902 scr': COL['scr'], 'ID06902 bmi1': COL['bmi'], 'ID06914 scr': '#9a9a9a'}

    examples = [('SCR, isolated neuron', 'ID06902/SCR/scr_data/New-11.czi', 'crop'),
                ('BMI1, isolated neuron', 'ID06902/BMI1/bmi1_data/New-28.czi', 'crop'),
                ('BMI1, dense cluster', 'ID06902/BMI1/bmi1_data/New-03.czi', 'full')]
    arrs = []
    for t, f, mode in examples:
        row, br, so, a = mu.analyze_image(os.path.join(BASE, '2WPT', f))
        arrs.append((t, a, mode))

    fig = plt.figure(figsize=(7.2, 5.4))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.1, 1], hspace=0.42, left=0.075, right=0.985, top=0.95, bottom=0.11)
    top = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], wspace=0.06)
    axA = []
    for j, (t, a, mode) in enumerate(arrs):
        ax = fig.add_subplot(top[0, j])
        img = a['neuron']
        if mode == 'crop' and a['somata'].any():
            # crop to the cell body and the dendrites assigned to it
            sid = np.bincount(a['territory'][a['territory'] > 0]).argmax()
            rr, cc = np.where(a['territory'] == sid)
            pad = int(6 / mu.TARGET_UM)
            cy, cx = (rr.min() + rr.max()) // 2, (cc.min() + cc.max()) // 2
            half = max(rr.max() - rr.min(), cc.max() - cc.min()) // 2 + pad
            y0, x0 = max(cy - half, 0), max(cx - half, 0)
            sl = (slice(y0, cy + half), slice(x0, cx + half))
        else:
            sl = (slice(None), slice(None))
        show_gray(ax, img[sl], hi=99.8)
        if a['somata'].any():
            ax.contour((a['somata'] > 0)[sl], [0.5], colors='#00e5ff', linewidths=0.7)
        if mode == 'full' and a['dense'][sl].any():
            ax.contour(a['dense'][sl], [0.5], colors='#ffd400', linewidths=0.7, linestyles='--')
        scalebar(ax, 50 if mode == 'full' else 10, mu.TARGET_UM, img[sl].shape)
        ax.set_title(t, fontsize=7)
        axA.append(ax)

    bot = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[1], wspace=0.6)
    panels = [('dense_frac', 'Fraction of extension length\nin dense clusters', 1.0),
              ('iso_branch_points_per_100um', 'Branch points per 100 µm\n(isolated extensions)', 1.0),
              ('iso_diam_median_um', 'Extension diameter, µm\n(isolated extensions)', 1.0),
              ('iso_pct_branches_swollen', 'Branches with ≥ 1 swelling, %\n(isolated extensions)', 1.0)]
    axes_b = []
    for j, (key, ylab, _) in enumerate(panels):
        ax = fig.add_subplot(bot[0, j])
        g = {k: im.loc[im.grp == k, key].to_numpy() for k in G}
        p1 = mwu(g['ID06902 bmi1'], g['ID06902 scr']); p2 = mwu(g['ID06902 bmi1'], g['ID06914 scr'])
        strip(ax, g, G, GC, GL, ylabel=ylab,
              comparisons=[('ID06902 scr', 'ID06902 bmi1', p1), ('ID06902 bmi1', 'ID06914 scr', p2)], seed=j)
        ax.tick_params(axis='x', labelsize=6)
        axes_b.append(ax)
        n = {k: int(np.isfinite(v).sum()) for k, v in g.items()}
        med = {k: float(np.nanmedian(v)) for k, v in g.items()}
        mean = {k: float(np.nanmean(v)) for k, v in g.items()}
        record('Morphology (2WPT CZI)', key, 'BMI1 vs SCR same mouse / vs SCR second mouse (per image, Mann-Whitney)',
               n_bmi=n['ID06902 bmi1'], n_scr=n['ID06902 scr'], n_scr2=n['ID06914 scr'],
               mean_bmi=mean['ID06902 bmi1'], mean_scr=mean['ID06902 scr'], mean_scr2=mean['ID06914 scr'],
               median_bmi=med['ID06902 bmi1'], median_scr=med['ID06902 scr'], median_scr2=med['ID06914 scr'],
               p=p1, p_vs_scr2=p2)
    panel_letter(fig, axA[0], 'A', dx=-0.055, dy=0.012)
    for ax, L in zip(axes_b, 'BCDE'):
        panel_letter(fig, ax, L, dx=-0.065)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, f'Figure3_morphology.{ext}'))
    plt.close(fig)
    so = pd.read_csv(os.path.join(BASE, 'results_2WPT/czi/csv/morphology_somata_czi.csv'))
    so['grp'] = so['mouse'] + ' ' + so['condition']
    for k in G:
        s = so[so.grp == k]
        record('Morphology (2WPT CZI)', 'cell bodies', k, n_somata=len(s), n_in_cluster=int(s.in_dense_cluster.sum()))


if __name__ == '__main__':
    which = set(sys.argv[1:]) or {'1', '2', '3'}
    r1 = figure1() if '1' in which else {}
    r2 = figure2() if '2' in which else {}
    if '3' in which:
        figure3()
    st = pd.DataFrame(STATS)
    if which == {'1', '2', '3'}:
        st.to_csv(os.path.join(FIG, 'stats_summary.csv'), index=False)
    print(st.round(4).to_string())
    print('Figure 1 examples:', {k: (v['fn'], round(v['val'], 2)) for k, v in r1.items()})
    print('Figure 2 examples:', {k: (v['fn'], round(v['val'], 2)) for k, v in r2.items()})
