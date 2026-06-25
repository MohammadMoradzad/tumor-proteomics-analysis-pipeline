# tumor-proteomics-analysis-pipeline

This repository contains a R pipeline for the analysis of quantitative proteomics data generated in the context of IL-15 and IL-15Rα-dependent tumour immune surveillance. The workflow was developed to analyze proteomic changes in MC A -induced tumours and tumor-derived cell lines from WT, Il15−/− and Il15ra−/− mice with downstream pathway-level interpretation using KEGG and Hallmark gene-set enrichment analysis.

This analysis supports the proteomics component of a study to investigate if IL-15 and IL-15 trans-presentation via IL-15Rα have different roles in antitumor immunity, tumour immunosurveillance, tumour immunoediting and antigen-presentation-associated pathways.

## Project overview

IL-15 is an important cytokine for the development, maintenance, and activation of several immune cell populations, including NK cells and CD8+ T cells. IL-15 can be presented by IL-15Rα to cells expressing IL-15Rβ/γc, a process known as trans-presentation. In this project, quantitative proteomics was used to investigate how loss of IL-15 or IL-15Rα affects tumor proteomic programs, particularly pathways related to interferon signaling, antigen processing and presentation, immune regulation, and tumor immunoediting.

The pipeline includes:

1) Proteomics data import and quality control
2) Log2 transformation and matrix formatting
3) Differential protein abundance analysis
4) Volcano plot generation
5) KEGG gene-set enrichment analysis
6) Hallmark gene-set enrichment analysis
7) Targeted pathway heatmaps using GSEA-derived gene lists

## Repository structure

```text
tumor-proteomics-analysis-pipeline/
│
├── annotation/
│   └── h.all.v2026.1.Mm.symbols.gmt
│
├── data/
│   └── processed/
│       ├── Raw proteomics intensity_log2 transformed.xlsx
│       └── tumor_proteomics_log2_matrix.xlsx
│
├── docs/
│
├── outputs/
│   ├── Differential_expression/
│   ├── GSEA/
│   ├── GSEA_results/
│   │   ├── Hallmark/
│   │   └── KEGG/
│   ├── Heatmaps/
│   ├── QC/
│   └── Volcano/
│
├── scripts/
│   ├── 01_import_and_QC.R
│   ├── 02_differential_expression_limma_volcano.R
│   ├── 03_KEGG_GSEA_analysis.R
│   ├── 04_Hallmark_GSEA_analysis.R
│   └── 07_targeted_GSEA_pathway_heatmaps.R
│
└── README.md
```

## Input data

The main input files are protein abundance matrices generated from DIA mass spectrometry data. Processed log2-transformed protein abundance matrices are stored in:

```text
data/processed/

These files are used for quality control, differential protein abundance analysis, and pathway-level visualization.

Gene-set annotation files used for enrichment analysis are stored in:

```text
annotation/
```

For Hallmark GSEA, the mouse MSigDB Hallmark GMT file was used:

```text
annotation/h.all.v2026.1.Mm.symbols.gmt
```

## Scripts

### 01_import_and_QC.R

This script imports the processed proteomics matrix and performs quality control analysis. The QC workflow includes sample-level inspection, detection summaries, missing-value assessment, PCA, sample-to-sample correlation analysis, Euclidean distance analysis, and log2 intensity distribution plots.

Main outputs:

```text
outputs/QC/
```

### 02_differential_expression_limma_volcano.R

This script performs differential protein abundance analysis using limma-based statistical modeling. It compares experimental groups and generates volcano plots based on protein-level log2 fold changes and statistical significance.

Main outputs:

```text
outputs/Differential_expression/
outputs/Volcano/
```

### 03_KEGG_GSEA_analysis.R

This script performs KEGG pathway enrichment analysis using ranked protein/gene lists generated from differential protein abundance results. It identifies positively and negatively enriched KEGG pathways and generates summary plots.

Main outputs:

```text
outputs/GSEA/
outputs/GSEA_results/KEGG/
```

### 04_Hallmark_GSEA_analysis.R

This script performs Hallmark gene-set enrichment analysis using MSigDB Hallmark gene sets. It identifies enriched biological programs such as interferon response, inflammatory signaling, cell-cycle regulation, apoptosis, oxidative phosphorylation, and other cancer-relevant pathways.

Main outputs:

```text
outputs/GSEA/
outputs/GSEA_results/Hallmark/
```

### 07_targeted_GSEA_pathway_heatmaps.R

This script uses GSEA output files to extract pathway-associated core-enrichment genes and generate targeted heatmaps from the proteomics expression matrix. This was used to visualize selected KEGG and Hallmark pathways, including interferon-associated pathways and other immune or tumor-biology-related signatures.

Main outputs:

```text
outputs/GSEA_results/Hallmark/
outputs/GSEA_results/KEGG/
outputs/Heatmaps/
```

## Output files

The repository contains both tabular results and figures generated from the proteomics workflow.

### Quality control outputs

```text
outputs/QC/
```

Includes PCA, detected protein summaries, missing-value plots, correlation heatmaps, Euclidean distance heatmaps, log2 intensity boxplots, and combined QC panels.

### Differential protein abundance outputs

```text
outputs/Differential_expression/
```

Includes differential protein abundance tables, statistical summaries, and comparison-specific result files.

### Volcano plots

```text
outputs/Volcano/
```

Includes volcano plots generated from differential protein abundance analysis.

### GSEA outputs

```text
outputs/GSEA/
outputs/GSEA_results/
```

Includes KEGG and Hallmark enrichment results, top enriched pathway plots, and Excel files containing pathway-level statistics and core-enrichment gene information.

### Targeted pathway heatmaps

```text
outputs/GSEA_results/Hallmark/
outputs/GSEA_results/KEGG/
outputs/Heatmaps/
```

Includes pathway-specific heatmaps generated using genes extracted from KEGG and Hallmark GSEA results.

## Required R packages

The analysis uses the following R packages:

```r
readxl
openxlsx
dplyr
tidyr
tibble
stringr
ggplot2
ggrepel
pheatmap
limma
clusterProfiler
fgsea
msigdbr
enrichplot
patchwork
ggplotify
```

Some packages may need to be installed from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "limma",
  "clusterProfiler",
  "fgsea",
  "enrichplot"
))
```

## How to run the pipeline

Run the scripts in numerical order:

```r
source("scripts/01_import_and_QC.R")
source("scripts/02_differential_expression_limma_volcano.R")
source("scripts/03_KEGG_GSEA_analysis.R")
source("scripts/04_Hallmark_GSEA_analysis.R")
source("scripts/07_targeted_GSEA_pathway_heatmaps.R")
```

Each script saves its corresponding figures and tables into the `outputs/` directory.

## Notes on reproducibility

This repository is organized so that input matrices, annotation files, analysis scripts, and generated outputs are separated clearly. The goal is to make the proteomics workflow transparent, reproducible, and easy to inspect.

Large raw mass spectrometry files are not stored directly in this repository. Mass spectrometry raw data will be deposited in a public proteomics repository such as PRIDE/ProteomeXchange when available. This GitHub repository is intended to document the downstream R-based analysis workflow.

## Citation

If using this workflow, please cite the associated manuscript:

Rexhepi et al. Trans-presentation of IL-15 by IL-15Rα attenuates tumor immune surveillance and is dispensable for IL-15-dependent tumor growth control.

## Author

Mohammad Moradzad
Department of Immunology and Cell Biology
Université de Sherbrooke
