# Axon break detection — FAST U-Net++ masks

This folder repeats the **broken axon detection** from `results/`, but builds the axon mask with a
pretrained deep-learning model instead of Frangi + Otsu. The **spot detection** results
(`spot_detection_*.ipynb`, `csv/spot_counts_*.csv`, `figures/spot_count_comparison_*.png`) are copied
here unchanged from `results/`, so that each dataset folder holds both analyses. Spotiflow does not
use the axon mask. The copied notebooks write to this folder if re-run. The layout mirrors `results/`: one subfolder
per dataset (`N1/`, `N2/`, `N3/`, `N4/`), each with its notebooks, `csv/` and `figures/`.

## Why a different mask

`pretrained_axon_detection.ipynb` (top level) compared Frangi + Otsu with three pretrained models
(DeepD3, FAST U-Net++, vesselFM) and two ridge filters on N1. The Frangi + Otsu skeleton covers only
about 16% of the axon network that the other methods find. Most of its "breaks" are places where the
mask lost a dim but continuous axon. The **FAST U-Net++** filament model (BioImage.IO
`reliable-waterbuffalo`, doi 10.1242/jcs.264265) traced almost the whole network and gave the largest
share of real (dark-gap) breaks. It is the mask used here.

## Method, in short

1. Load the CZI image and max-project across all Z planes (channel 0 = tubulin).
2. **Resample to 87.7 nm pixels** (the N1 pixel size). N2 mixes 87.7, 73.1 and 35.3 nm pixels
   (Airyscan is 35.3 nm); N3/N4 are 70.6 nm. Resampling keeps the break parameters, which are in pixels, the same
   physical size everywhere. It also matters for the model: on native-resolution Airyscan images,
   FAST misses many axons.
3. Run **FAST U-Net++** on the raw projection; the model normalizes internally. No background
   subtraction is applied, for any dataset. Axon mask = summed softmax of the *actin & stress fibers*
   and *filopodia* classes > 0.5.
4. Skeletonize, drop fragments < 15 px, and detect break pairs with the **same** endpoint-direction +
   gap-clarity detector and parameters as `results/`.
5. **Confirm each break against the raw image:**
   `gap contrast = (gap intensity − background) / (axon intensity at the two ends − background)`.
   A break is *confirmed* when the contrast is below 0.35, i.e. the gap is genuinely dark rather than
   a dim stretch of axon.

## Figures per dataset

- `figures/axon_break_comparison_FAST_*.png`: the same 4-panel plot as `results/` (fragment rate,
  endpoint rate, gap-bridge rate, continuity), with Kruskal-Wallis tests.
- `figures/axon_break_counts_FAST_*.png`: break rate, **confirmed break rate** (the recommended
  readout), and % of breaks confirmed.
- `csv/axon_break_metrics_FAST_*.csv`: per-image table, including `pixel_size_nm`, `n_breaks`,
  `n_confirmed`, `confirmed_rate` and the skeleton metrics.

## Dataset notes

- **N1**: standard confocal, 17 images (bmi 9 / Ni 4 / scr 4).
- **N2**: regular acquisition (`*_without_airyscan`) and Airyscan-processed (`*_airyscan`) results
  are kept **separate**, as in `results/`.
- **N3**: 3-channel acquisition (AF647 tubulin, AF488, DAPI), 70.6 nm pixels, 54 images
  (bmi 18 / Ni 18 / scr 18). Channel 0 (AF647) is tubulin, as in N1/N2.
- **N4**: same acquisition as N3, but **only bmi and scr** (9 / 8 images, no Ni). The images are
  sparse, with thin isolated axons and large cell-body clusters. Two versions:
  - `axon_break_detection_FAST_N4.ipynb`: same pipeline as the other datasets.
  - `axon_break_detection_FAST_N4_soma_masked.ipynb` (outputs `*_soma_masked.*`): cell bodies are
    removed first. Candidates come from a tubulin opening (radius 12 px, ≥ 25 µm²), kept only where
    they touch a DAPI nucleus. Endpoints and break pairs at cell-body edges are dropped. This removes
    about 30% of the skeleton (filament texture inside cell bodies) and half of the breaks.
  - In both versions only 4–17% of N4 breaks are confirmed, against 30–45% in N1–N3. The main reason
    is not cell bodies: FAST breaks up thin, beaded, isolated axons into dotted segments, so most raw
    break pairs sit on an axon that is still visible (gap contrast around 1). The confirmation step
    rejects these, so use `confirmed_rate` for N4 and do not use the raw break count.

## Requirements

- Kernel: `axon`.
- Model files: `../pretrained_models/fast_actin/` (`model.py` + `fast_model.pth`).
