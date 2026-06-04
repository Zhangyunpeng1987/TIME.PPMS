# TIME.PPMS

## Overview

This repository provides the R code used in the study:

"Pan-ancer single-ell dissection of coordinated cellular programs reveals immune remodeling and potential drug targets in tumorigenesis.

The project aims to systematically characterize tumor immune
microenvironment (TIME) cellular programs (NMF1-4) between paired tumor-normal samples, and to construct a pan-cancer prognostic model signature (PPMS).

------------------------------------------------------------------------

## Main Analyses Included

-   Single-cell RNA-seq data integration and quality control (Seurat + scib-metrics)
-   TIME cellular programs identification (NMF)
-   Differential expression and enrichment analysis (tidyverse + clusterProfiler)
-   Gene set scoring (GSVA)
-   Cancer cell state analysis (copykat + infercnv)
-   TF activity (decoupleR)
-   Cell-cell communication (CellChat)
-   KEGG metabolic activity analysis (AddModuleScore)
-   TCGA pan-cancer analysis (TCGAplot)
-   Prognostic model construction (Mime1)
-   Immune infiltration analysis 
-   Drug sensitivity prediction (oncoPredict)
-   Molecular docking analysis
-   Immunohistochemistry validation (Human Protein Atlas)
-   CRISPR screening data analysis (DepMap)
-   Spatial transcriptomic data analysis (SpaCET + SpaGene)


------------------------------------------------------------------------

## Code Organization

The scripts are organized by figure to facilitate reproducibility. Each folder contains the code used to generate the corresponding panels in the manuscript.

For example:
- `code/Result1/` contains scripts for main Figures and supplementary Figures related to Result 1.
- `code/Result2/` contains scripts for main Figures and supplementary Figures related to Result 2.
- Subsequent folders correspond to the remaining figures.

Users can follow these scripts to reproduce the results step by step.

------------------------------------------------------------------------

## Repository Structure

NSCLC-DMSPsig/
├── README.md
├── LICENSE
├── code/
│   ├── Remove_double_cells.R
│   ├── Result1/
│   │   ├── Result1.R
│   │
│   ├── Result2/
│   │   ├── Result2.R
│   │
│   ├── Result3/
│   │   ├── Result3.R
│   │
│   ├── Result4/
│   │   ├── Result4.R
│   │
│   ├── Result5/
│   │   ├── Result5.R
│   │
│   ├── Result6/
│   │   └── Result6.R
│   │
│   ├── Result7/
│   │   ├── Result7.R
│   │
│   ├── Result8/
│   │   ├── Result8.R
│
└── sessionInfo.txt

------------------------------------------------------------------------

## Data Availability

The datasets used in this study are publicly available:

- Single-cell RNA-seq and spatial transcriptomic datasets were obtained from the Gene Expression Omnibus (GEO) database (https://www.ncbi.nlm.nih.gov/geo/), with accession numbers provided in the manuscript.
- TCGA bulk RNA-seq data: https://portal.gdc.cancer.gov/
- Human Protein Atlas: https://www.proteinatlas.org/
- CRISPR screening data from the DepMap database: https://depmap.org/portal/

------------------------------------------------------------------------

## Requirements

R = 4.4.1

Major packages: Seurat, scib-metrics, NMF, tidyverse, GSVA, glmnet, CellChat, decoupleR, copykat, infercnv, oncoPredict, TCGAplot, Mime1, SpaCET, SpaGene, survival, survminer, etc.

------------------------------------------------------------------------

## Contact

Yunpeng Zhang: zhangyp@hrbmu.edu.cn\
Congxue Hu: hucx1996@hrbmu.edu.cn\
Xia Li: lixia@hrbmu.edu.cn

------------------------------------------------------------------------

## License

MIT License
