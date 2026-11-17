# Spot detection (tubulin puncta)

Cultured neurons stained for βIII-tubulin. We count the bright puncta along the axons.

## How puncta are counted
Spotiflow, a pretrained deep learning spot detector, is run on the tubulin image with its default "general" model. Each detected spot is counted once. Fields differ in size between experiments, so counts are divided by the imaged area.

## File to open: `puncta_per_image.csv`
One row per image.

| Column | Meaning |
|---|---|
| `experiment` | Exp. 1 to 4. "Exp. 2 (Airyscan)" is a separate set of Airyscan images of Exp. 2. |
| `condition` | BMI1, SCR or NI |
| `image` | Image file name |
| `puncta` | Number of puncta in the image |
| `imaged_area_um2` | Area of the image in µm² |
| `puncta_per_1000um2` | **The number used in the paper.** Puncta per 1,000 µm². |

## Result (mean per condition, puncta per 1,000 µm²)
| Experiment | BMI1 | SCR | BMI1 / SCR | P |
|---|---|---|---|---|
| Exp. 1 (9 vs 4 images) | 44.9 | 18.1 | 2.5 | 0.034 |
| Exp. 2 (23 vs 19) | 33.2 | 14.6 | 2.3 | 0.003 |
| Exp. 3 (18 vs 18) | 61.9 | 23.6 | 2.6 | < 0.001 |
| Exp. 4 (9 vs 8) | 12.5 | 5.1 | 2.5 | 0.004 |
| All four | | | **2.5** | **< 0.001** |

The Airyscan images of Exp. 2 show no BMI1 increase (BMI1 91.6, NI 104.5, SCR 102.7). Airyscan images have a different resolution, so keep them apart from the regular images.

## `original_csv`
The raw count files, one per experiment (N1 to N4 are Exp. 1 to 4). They hold only `condition`, `filename` and `n_spots`. The area is added in the per-image file.
