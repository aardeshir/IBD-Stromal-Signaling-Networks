# IBD Stromal Cell Analysis

This repository contains analysis scripts for characterizing stromal cell populations and BMP2 signaling in Inflammatory Bowel Disease (IBD) using single-cell RNA-sequencing data.

## Overview

The analysis focuses on understanding stromal cell heterogeneity and BMP2 signaling pathways in the context of IBD using the Kinchen et al. (2018) dataset of colonic mesenchymal cells from healthy controls and ulcerative colitis patients.

## Repository Structure

```
IBD_Stromal_Analysis/
├── data/                  # Data directory (not included in repo)
├── scripts/               # Analysis scripts
│   ├── 01-preprocessing.R # Data loading and quality control
│   ├── 02-stromal_subtyping.R # Identification of stromal subtypes
│   ├── 03-pathway_correlation_analysis.R # Analysis of pathway correlations
│   ├── 04-receptor_ligand_analysis.R # Receptor-ligand interaction analysis
│   └── 05-statistical_tests.R # Statistical testing for BMP2 analysis
├── results/               # Analysis results (generated from scripts)
└── figures/               # Generated figures (generated from scripts)
```

## Analysis Workflow

1. **Preprocessing**: Load and prepare the Kinchen dataset
2. **Stromal Subtyping**: Identify stromal cell subtypes using marker genes
3. **Pathway Correlation Analysis**: Examine correlations between BMP2 and various signaling pathways
4. **Receptor-Ligand Analysis**: Analyze BMP2 receptor expression and potential signaling interactions
5. **Statistical Tests**: Perform statistical tests on BMP2 expression and relationships

## Requirements

- R (version 4.0.0 or later)
- Seurat (version 5.0.0 or later)
- dplyr
- ggplot2
- patchwork
- ComplexHeatmap
- circlize
- RColorBrewer

## Usage

1. Place the Kinchen dataset in the `data/` directory
2. Run the scripts in numerical order:

```R
Rscript scripts/01-preprocessing.R
Rscript scripts/02-stromal_subtyping.R
Rscript scripts/03-pathway_correlation_analysis.R
Rscript scripts/04-receptor_ligand_analysis.R
Rscript scripts/05-statistical_tests.R
```

## Data Source

This analysis uses the dataset from:

Kinchen J, Chen HH, Parikh K, et al. Structural Remodeling of the Human Colonic Mesenchyme in Inflammatory Bowel Disease. Cell. 2018;175(2):372-386.e17. doi:10.1016/j.cell.2018.08.067

## Citation

If you use this code in your research, please cite the original data source (Kinchen et al. 2018) and our analysis as appropriate.

## License

MIT License 