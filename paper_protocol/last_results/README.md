# Last results

Final numbers from the image analysis, for the biology team. There are three folders, one per question.

| Folder | Question | Main result |
|---|---|---|
| `broken_axons` | Do BMI1 neurons have more axon breaks? | In cultured neurons, yes: 2.8 times more than SCR across 4 experiments (P = 0.013). In the new in vivo MAP2 images, Condition 1 has more breaks than Condition 2 (labels to be confirmed). |
| `spot_detection` | Do BMI1 neurons have more bright tubulin puncta? | Yes. 2.5 times more than SCR across 4 experiments (P < 0.001). |
| `morphological_difference` | Does the shape of the neurons change? | In the mouse brain, BMI1 neurons sit in dense clusters. In isolated extensions they branch more and are thicker. Inside dense clusters no difference can be shown. |

## Words used in all folders
- **BMI1**: neurons with BMI1 knocked down. **SCR**: scrambled control. **NI**: non-infected neurons.
- **Exp. 1 to Exp. 4**: the four independent cultured-neuron experiments. In the file names they are called N1 to N4.
- **In vitro**: cultured neurons. **In vivo**: neurons in mouse brain tissue.
- **Image**: one microscope field. Every dot or row is one image. Images are not independent animals or cultures.
- **P**: Mann-Whitney test on images, BMI1 against SCR. Across experiments we compare the ratio BMI1 / SCR of the four experiments.
- **Isolated** and **dense** (morphology only): a stretch of extension is *dense* when the stained area covers at least 20% of the 20 µm window around it. Otherwise it is *isolated*.

## What is in each folder
- One sub-folder per data set. Inside it, a file ending in `_per_image.csv` (or `morphology_per_image.csv`). **Open this one first.** One row per image, plain column names, and the numbers used in the paper.
- A comparison file (`comparison_...csv` or `isolated_vs_dense_comparison.csv`) with the group medians and the P values.
- A folder `original_csv` with the raw output files of the analysis, unchanged. They have many extra columns for technical checks. You do not need them to read the results.
- A `README.md` that explains every column.

The figures (Figure 1 to 3) and the methods text are one level up, in `paper_protocol/figures` and `Image_analysis_protocol.docx`. The code and notebooks of each analysis are in the result folders of the project (`results_with_Fast-Unet`, `results`, `results_2WPT`, `results_MAP2`, `results_dendrite_imaris`).

## Opening the files
The files use a comma between columns and a dot for decimals. In Excel with French settings, use Data > From Text/CSV and choose "comma" and "dot", or the columns will not split correctly. An empty cell means there was not enough signal in that image to measure the value.

## Important limits
- Exp. 1 has only 4 SCR images, so its P value is weak.
- The mouse-brain results come from **one** BMI1 mouse. They need more animals.
- **MAP2 images: please help us.** The two conditions are called Condition 1 and Condition 2 because the labels were not given. Four of the nine Condition 2 images are the same files as in the `BMI1` folder, so Condition 2 is probably BMI1. Please confirm. The two conditions were also imaged with different microscope modes (see `broken_axons/README.md`).
- In the dense clusters, neighbouring neurons cross and overlap. This adds branch points and widens the extensions, so dense values must not be compared with isolated values.
