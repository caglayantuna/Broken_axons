# Broken axons

We count breaks along axons. There are two data sets, with the same method:
- `in_vitro_tubulin_Exp1-4`: cultured neurons stained for βIII-tubulin (Exp. 1 to 4). This is the data set of the paper.
- `in_vivo_MAP2`: new images of mouse brain tissue, stained for MAP2.

## How a break is counted
1. A pretrained model (FAST) finds the axons in each image.
2. Each axon is thinned to a centre line. The free ends of the lines are found.
3. Two ends that face each other across a gap of 0.9 to 7 µm count as a break.
4. The break is **confirmed** only if the gap is dim: at most about one third as bright as the axon next to it. This removes axons that are only faint.

The number we use is **confirmed breaks per 100 µm of axon**.

## Example for the poster
`example_pretrained_model_axon_detection.png` (and `.pdf`) shows the pretrained model on one image (Exp. 3, BMI1, `New-10.czi`). Left: the βIII-tubulin staining. Right: the axons found by the model (cyan). The bottom row is a 20 µm zoom of the yellow square. Made by `../make_poster_example.py`.

## Columns of `broken_axons_per_image.csv` (both sub-folders)
One row per image.

| Column | Meaning |
|---|---|
| `experiment` | (in vitro) Exp. 1 to 4. "Exp. 2 (Airyscan)" is a separate set of Airyscan images of Exp. 2. |
| `condition` | BMI1, SCR or NI (in vitro). Condition 1 or Condition 2 (in vivo). |
| `acquisition` | (in vivo) Regular confocal or Airyscan |
| `image` | Image file name |
| `confirmed_breaks_per_100um` | **The number used in the paper.** Confirmed breaks divided by axon length. |
| `confirmed_breaks` | Number of confirmed breaks in the image |
| `candidate_breaks` | Number of gaps found before the brightness check |
| `axon_length_um` | Total length of axon centre lines in the image, in µm |

## 1. In vitro, Exp. 1 to 4 (mean per condition, confirmed breaks per 100 µm)
| Experiment | BMI1 | SCR | BMI1 / SCR | P |
|---|---|---|---|---|
| Exp. 1 (9 vs 4 images) | 0.75 | 0.43 | 1.7 | 0.076 |
| Exp. 2 (23 vs 19) | 0.75 | 0.28 | 2.7 | < 0.001 |
| Exp. 3 (18 vs 18) | 1.57 | 0.53 | 3.0 | < 0.001 |
| Exp. 4 (9 vs 8) | 0.43 | 0.10 | 4.5 | 0.002 |
| All four | | | **2.8** | **0.013** |

NI neurons are close to SCR. The Airyscan images of Exp. 2 point the same way (BMI1 0.40, NI 0.33, SCR 0.08).

`original_csv`: the raw output of the best method (FAST), one file per experiment (N1 to N4 are Exp. 1 to 4). Useful columns: `n_confirmed`, `n_breaks`, `confirmed_rate` (confirmed breaks per 1,000 axon pixels, converted to per 100 µm in the per-image file) and `skel_length` (axon length in pixels, 1 pixel = 0.0877 µm). The other columns are extra descriptors that are not used in the paper.

## 2. In vivo, MAP2 images
Mouse brain tissue stained for MAP2 (a dendrite marker, which also stains cell bodies). **The labels were not given**, so the two conditions are called Condition 1 and Condition 2. Four of the nine Condition 2 images are the same files as in the `BMI1` folder, so Condition 2 is probably BMI1. **Please confirm.**

| | Condition 1 | Condition 2 |
|---|---|---|
| Images | 13 | 9 |
| Microscope mode | Regular confocal, 8-bit | Airyscan, 16-bit |
| Pixel size | 70.6 nm | 35 to 55 nm |
| Axon length per image (mean) | 5,428 µm | 1,977 µm |

The two conditions were imaged differently, so a difference can come from the microscope and not from the biology. To check this, Condition 2 was blurred to look like confocal and the analysis was repeated.

Result (mean per condition, confirmed breaks per 100 µm):

| | Condition 1 | Condition 2 | Condition 2 / Condition 1 | P |
|---|---|---|---|---|
| Confirmed breaks | 3.36 | 2.48 | 0.74 | 0.033 |
| Same, Condition 2 blurred | 3.36 | 2.38 | 0.71 | 0.053 |
| All candidate breaks (before the brightness check) | 7.50 | 6.80 | 0.91 | 0.29 |

- Condition 1 has **more** confirmed breaks. The effect is modest. The 95% interval of the ratio is 0.57 to 0.94, and after correcting for the nine measures tested, q = 0.15.
- Three Condition 2 images have less than 400 µm of axon, so their rates are noisy. Without them the P value is 0.046 (images with at least 300 µm) and 0.087 (at least 500 µm).
- If Condition 2 is BMI1, this is the **opposite** direction to the cultured neurons. In tissue, MAP2 also stains dendrites and cell bodies, so a break count here is not the same measure as in the tubulin images of the cultures.

Files: `comparison_Condition1_vs_Condition2.csv` (the three rows above, plain names). `original_csv`: the raw per-image output, the blurred-Condition-2 run, the full statistics table and the robustness table.
