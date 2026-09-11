# Axon analysis — results

This directory holds the per-dataset outputs (comparison plots + CSV tables) for the
two analyses run on the CZI images in `N1/`, `N2/`, `N3/`. Each dataset has its own
subfolder (`N1/`, `N2/`, `N3/`) containing the result notebooks and the files they produce.

Every condition is one of three siRNA/treatment groups: **bmi**, **Ni**, **scr**.

## 1. Broken axon detection

**Goal:** count how many axon breaks (discontinuities) are visible per image, and compare
that across conditions.

**Method, in short:**
1. Load the CZI image, max-project across all Z planes (channel 0 = tubulin).
2. Enhance tube-like structures with a **Frangi vesselness filter**, then threshold
   (Otsu) to get a binary axon mask. No morphological closing is applied, so real gaps
   in the axon are preserved rather than filled in.
3. **Skeletonize** the mask to 1-pixel-wide lines and drop tiny spur fragments
   (< 15 px).
4. Find the skeleton **endpoints** (tips of each fragment). For every pair of endpoints
   that (a) point toward each other and (b) has a dark, axon-free gap between them, the
   pair is counted as one **detected break**.
5. Per image this gives a **break count**, plus supporting skeleton metrics
   (fragment count, endpoint rate, gap-bridge rate, mean fragment length) that describe
   how fragmented the axon network is overall.

The comparison plot shows, per condition, the fragment rate, endpoint rate, gap-bridge
rate, and mean fragment length (continuity), with a Kruskal-Wallis test between groups.

## 2. Spot detection

**Goal:** count puncta (e.g. synaptic/axonal spots) per image and compare counts across
conditions.

**Method, in short:**
1. Load the CZI image, max-project across all Z planes (channel 0).
2. Run the pretrained **Spotiflow** ("general") deep-learning spot detector on the
   projected image to get spot coordinates.
3. The **number of spots** detected is the per-image result.

The comparison plot shows a boxplot of spot counts per condition, with each image's
count shown as a labeled point.

## Dataset notes

- **N1** — standard confocal acquisition, 3 conditions (bmi/Ni/scr), single image set per
  condition.
- **N2** — has both a regular acquisition and an **Airyscan-processed** acquisition for
  each condition. Results for the two are kept **separate** (different notebooks/files),
  since Airyscan processing changes image resolution/noise characteristics enough that
  the two should not be pooled. Spot detection on N2 is run **without preprocessing**
  only (raw max-projection, no background subtraction) — broken axon detection on N2
  still uses the background-subtraction preprocessing step, matching what worked in
  the source notebooks.
- **N3** — not yet acquired. Placeholder notebooks are included with the expected
  folder layout (mirroring N1); update the `directories` paths once the data lands and
  re-run.

## Files per dataset folder

Each `N1/`, `N2/`, `N3/` folder contains:

- `*.ipynb` — the notebook that runs the analysis for that dataset/variant (no inline
  image display; text/table output only).
- `figures/*_comparison_*.png` — the saved comparison plot.
- `csv/*_metrics_*.csv` / `csv/*_counts_*.csv` — per-image, per-condition results table.
