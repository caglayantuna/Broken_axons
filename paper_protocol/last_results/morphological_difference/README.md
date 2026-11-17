# Morphological difference

Two image sets of neurons grafted into the mouse brain (AF594 label). Each has its own folder with a `morphology_per_image.csv`, an `isolated_vs_dense_comparison.csv` and an `original_csv` folder.

**Isolated and dense are measured separately.** A stretch of extension is *dense* when the stained area covers at least 20% of the 20 µm window around it (inside a cluster of grafted cells). Otherwise it is *isolated*. Branch points and diameter are measured for each kind on its own, and BMI1 is compared with SCR in each. In a dense cluster, neighbouring processes cross and overlap. This adds branch points and widens the extensions in **every** group, so a dense value must not be compared with an isolated one.

Mouse 1 (ID06902) has BMI1 and SCR neurons. Mouse 2 (ID06914) has SCR neurons only and is a second control.

## Columns of `morphology_per_image.csv` (one row per image)
| Column | Meaning |
|---|---|
| `acquisition` | (2WPT only) Regular confocal (used in the paper) or Airyscan |
| `mouse`, `condition`, `image` | Mouse, BMI1 or SCR, image file name |
| `z_stack_thickness_um`, `same_field_in_2WPT` | (new images only) Thickness of the z-stack, and whether the same field is also in the 2WPT set |
| `dense_cluster_fraction` | Share of extension length inside dense clusters (0 to 1) |
| `isolated_...` and `dense_...` | The same measure for isolated and for dense extensions: |
| `..._branch_points_per_100um` | Branch points per 100 µm of extension |
| `..._diameter_um` | Median width of the extensions, in µm |
| `..._branch_length_um`, `..._tortuosity` | Median branch length. Median tortuosity (1 means perfectly straight). |
| `..._branches_with_swelling_pct`, `..._swellings_per_100um` | Percent of branches with at least one swelling. Swellings per 100 µm. |
| `..._extension_length_um` | Length of extension of that kind in the image |
| `cell_bodies`, `field_area_um2` | Number of cell bodies found. Area of the image in µm². |

An empty cell means the image had too little extension of that kind to measure the value (for example, SCR images almost never have dense extension).

`isolated_vs_dense_comparison.csv` has one row per measure and regime: the median per group, the number of images that could be measured, the ratio BMI1 / SCR (mouse 1) and the P values.

## 1. `in_vivo_grafted_neurons_2WPT` (regular confocal, 72 images)
Median per image, BMI1 against SCR from the same mouse:

| Measure | Regime | SCR mouse 1 | BMI1 | SCR mouse 2 | P (BMI1 vs SCR mouse 1) |
|---|---|---|---|---|---|
| Branch points per 100 µm | Isolated | 6.4 | 9.5 | 6.9 | 0.013 |
| | Dense | 14.8 | 16.5 | 19.5 | 0.41 |
| Diameter (µm) | Isolated | 0.72 | 0.88 | 0.86 | 0.013 |
| | Dense | 1.02 | 1.21 | 1.28 | 0.43 |
| Swellings | Isolated | 3.1% of branches | 2.3% | 3.1% | 0.69 |

- **Dense clusters:** a median of 37% of BMI1 extension length lies in dense clusters, against 0% (SCR, mouse 1) and 2% (SCR, mouse 2). 35 of 41 BMI1 cell bodies sit in a cluster.
- **Isolated extensions:** BMI1 branches more (x1.5) and is thicker. SCR in mouse 2 is as thick as BMI1, so the diameter cannot be credited to BMI1.
- **Dense extensions:** no BMI1 vs SCR difference, but only 4 (mouse 1) and 5 (mouse 2) SCR images have enough dense extension to measure, against 19 BMI1 images. The test has little power.
- **The dense regime changes every group.** Inside BMI1 images, dense extensions have more branch points (16.3 vs 10.7) and are thicker (1.21 vs 0.90 µm) than isolated ones in the same image. SCR in mouse 2 is thicker in dense areas too (0.91 to 1.28 µm).
- **Airyscan images (58)** give the same pattern. They are in the same files, in the `acquisition` column.
- All BMI1 neurons come from one mouse.

## 2. `dendrite_imaris_new_images` (new Airyscan z-stacks, 29 images)
New z-stacks of the same two mice, prepared for dendrite analysis in Imaris: 9 BMI1, 13 SCR (mouse 1), 7 SCR (mouse 2). Four BMI1 fields (`New-28` to `New-31`) are the same fields as 2WPT images, here as full stacks.

| Measure | Regime | SCR mouse 1 | BMI1 | SCR mouse 2 | P (BMI1 vs SCR mouse 1) |
|---|---|---|---|---|---|
| Branch points per 100 µm | Isolated | 6.1 | 6.1 | 7.7 | 0.59 |
| | Dense (2 / 5 / 2 images) | 12.5 | 15.7 | 14.8 | 0.38 |
| Diameter (µm) | Isolated | 0.63 | 0.74 | 0.77 | 0.16 |
| | Dense (9 / 7 / 5 images) | 1.19 | 0.99 | 1.00 | 0.41 |

- With dense and isolated extensions pooled, BMI1 is thicker (0.83 vs 0.64 µm, P = 0.038) and has more branch points (9.0 vs 6.0, P = 0.053). For isolated extensions alone there is no difference. The pooled difference comes from dense areas (share of extension length in dense clusters: BMI1 0.31, SCR 0.06 and 0.09).
- The stacks are 2.6 to 21 µm thick. If only a 3 µm slab around the sharpest plane is used, isolated branch points become 7.3 vs 6.1 (P = 0.055) and diameter 0.91 vs 0.73 µm (P = 0.071). The isolated result therefore depends on how much of the stack is projected. Both versions are in `isolated_vs_dense_comparison.csv` (column `projection`).
- **Per-neuron dendrite length and area cannot be tested.** Only 8 isolated cell bodies were found in 29 images (2 BMI1, 6 SCR). They are in `original_csv/isolated_neurons_dendrite_imaris.csv`. The fields are 78 µm wide and many show processes without a complete cell body.
- The dense label also appears around isolated cell bodies, because the cell body itself fills the 20 µm window.
- Other folders in `Dendrite Imaris` (`BMI1/3 mice`, `NI`, `SCR/Raw data for Caglayan`) are older 2023 acquisitions from other mice and were not analysed.

`original_csv` holds the full per-image files, the slab run, the statistics tables and the per-neuron table. The per-branch tables are large and stay in `results_2WPT` and `results_dendrite_imaris`.
