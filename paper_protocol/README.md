# paper_protocol

Image-analysis methods, figures and references for the BMI1 manuscript.

- `Image_analysis_protocol.docx`: methods text, short results, figure legends, references. The PDF is the same document.
  Red italic notes mark wet-lab details for the biologists to complete.
- `figures/Figure1_axon_breaks`, `Figure2_puncta`, `Figure3_morphology` (`.png` 300 dpi, `.pdf` vector).
- `figures/stats_summary.csv`: every number quoted in the text.
- `make_figures.py`: redraws the figures from the CSV files in `results/`, `results_with_Fast-Unet/` and
  `results_2WPT/`. `make_protocol_docx.py`: rebuilds the Word file from `stats_summary.csv`.

Run both with the `axon` env (`~/axon/bin/python`), in that order. `python-docx` was added to the env for the
Word file.
