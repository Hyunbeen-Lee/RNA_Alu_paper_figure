# Lee, Cao, Parada et al. 2025

This repository contains the code and analyses supporting our paper:

**Alu-mediated RNA duplexes are associated with widespread exon skipping across primate transcriptomes**

## Overview

Our work integrates computational prediction of intronic RNA stem-loop structures with experimental validation using SPLASH to reveal a prominent role for inverted Alu repeat elements in driving exon skipping.

## Repository Structure

- **README.md** – This file.
- **Data/** – Contains input files such as:
  - `alu_meta_ext.txt` – Annotation of intronic Alu pairs.
  - Plotting data files (e.g. `Chimp_plot_data.tsv`, `orangutan_plot_data.tsv`, etc.).
  - EvoPSI subfolder with compressed tables for pairwise comparisons.
- **Figures/** – Generated figures (e.g. `Fig_1A.PNG`, `Figure1.pdf`, `Figure2.pdf`, etc.).
- **Notebooks/** – Jupyter notebooks for data analysis and plotting:
  - [`Figure2.ipynb`](Notebooks/Figure2.ipynb) – Contains code for primate enrichment and pairwise comparisons.
  - [`Figure1.ipynb`](Notebooks/Figure1.ipynb) – Provides code for enrichment plots and additional comparisons.

## Dependencies

Our analyses require R (with libraries such as `ggplot2`, `data.table`, `cowplot`, etc.) and Python (with packages such as `pandas`, `numpy`, `scipy`, and others). Ensure that you have both environments properly set up and integrated via rpy2 (used in our notebooks).


## Citation

If you use our code or data in your work, please cite our paper.
