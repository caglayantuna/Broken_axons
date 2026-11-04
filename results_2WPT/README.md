# Grafted neuron morphology — 2WPT (bmi1 vs scr)

In vivo grafted neurons, AF594 label. Question: any morphological difference between bmi1 and scr,
in particular the diameter of the extensions.

- **ID06902:** scr and bmi1.
- **ID06914:** scr only (second-mouse reference).

| Notebook | Input | Output |
|---|---|---|
| `morphology_2WPT_czi.ipynb` | regular CZI files, max-Z projection (72 images) | `czi/csv`, `czi/figures` |
| `morphology_2WPT_airyscan.ipynb` | `*Airyscan Processing*.czi`, single planes (58 images) | `airyscan/csv`, `airyscan/figures` |

Both use `morphology_utils.py`. Kernel: `axon`; needs `cellpose` (Cellpose-SAM) and `skan`. The FAST
model files are read from `../pretrained_models/fast_actin/`.

## Method, in short
1. AF594 channel, picked by name (channel order changes between files), resampled to 70.6 nm.
2. Background subtraction and noise normalisation, then the AF594 mask: the FAST U-Net++ filament mask
   plus SNR > 3, kept where the piece contains SNR > 6.
3. Cell bodies: Cellpose-SAM at 0.12 µm/px, filtered on size (30–900 µm²), solidity and AF594 signal.
4. Extensions: skeleton (skan), branch length, tortuosity, branch points.
5. **Extension diameter: FWHM of the intensity profile across each branch, every 0.5 µm.** This is
   threshold-free.
6. **Swellings** (Flamier et al. 2018, Fig 3F). A swelling is a run of width samples that is:
   - wider than 2× the branch's median width, and at least 0.3 µm wider;
   - brighter than 1.2× the branch's median peak intensity;
   - wider than its neighbours;
   - at least 2 µm from either branch end.

   Reported as % of branches with ≥ 1 swelling and as swellings per 100 µm.
7. **Dense clusters vs isolated neurons:** local AF594 mask fraction in a 20 µm window, above 0.20 =
   dense graft core. Extension shape is compared on isolated extensions only.

## Files
- `csv/morphology_images_*.csv`: one row per image (all metrics, including `dense_frac` and `iso_*`).
- `csv/morphology_branches_*.csv`: one row per skeleton branch (length, tortuosity, diameter, CV,
  local density).
- `csv/morphology_somata_*.csv`: one row per cell body.
- `csv/stats_*.csv`: group medians and Mann-Whitney p values.
- `figures/`: channel checks, segmentation examples, FWHM diameter check, metric plots and branch
  distributions.

## Main result
- **Clustering (strongest):** bmi1 neurons sit in dense clusters (29–37% of extension length) while
  scr neurons are isolated (0–2%).
- **Branching:** isolated bmi1 extensions have more branch points (about 1.4x). This is significant
  in CZI and a trend in Airyscan.
- **Diameter and soma size:** bmi1 differs from scr in the same mouse, but the second scr mouse
  matches or exceeds bmi1, so these differences cannot be credited to bmi1.
- **Swellings:** rare in all groups (about 1 per 100 µm, 2–3% of branches), with no bmi1 vs scr
  difference. The detector is strict, and some detections sit at bends or branch points.
- **Main limit:** bmi1 comes from a single mouse.
