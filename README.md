# Microbiome Analysis Pipeline

This repository contains a comprehensive and modular R-based pipeline designed for the analysis of microbiome data (16S or Metagenomics). It leverages the `phyloseq` framework to perform a complete ecological analysis, ranging from data quality control to advanced multivariate statistics, co-occurrence networks, and indicator species analysis.

## Features

This pipeline automates the entire workflow, generating publication-ready figures (TIFF/600dpi) and statistical tables (TSV) for the following analyses:

* **Quality Control:** Rarefaction curves and Good’s Coverage estimation.
* **Alpha Diversity:** Richness, Shannon, Simpson, and Pielou's Evenness indices (with Kruskal-Wallis statistical testing).
* **Beta Diversity:** PCoA based on Bray-Curtis dissimilarity (with PERMANOVA testing).
* **Taxonomic Composition:** Abundance bar plots at multiple taxonomic levels.
* **Core Microbiome:** Heatmap visualization of core taxa.
* **Differential Abundance:** DESeq2 analysis for finding statistically significant biomarkers.
* **Environmental Drivers:** Redundancy Analysis (RDA) with forward selection of significant environmental variables.
* **Co-occurrence Networks:** Spearman-based networks with topological analysis (Degree, Betweenness, Closeness) and visual clustering.
* **Indicator Species (IndVal):** Identification of bioindicator taxa using the Indicator Species Analysis (indicspecies package).

## Getting Started

### Prerequisites
You need R and RStudio installed. If you are on Windows, ensure you have Rtools installed to build the required packages.

### Installation
Run the installation block provided in the script once to set up all necessary dependencies, including Bioconductor packages (phyloseq, DESeq2).

## Configuration (Control Panel)
The script is designed for ease of use. You do not need to modify the processing logic. Simply update the "Section 3: Control Panel" at the top of the script with your parameters:

* **File Paths:** Set `ARQUIVO_ABUNDANCIA` and `ARQUIVO_METADADOS`.
* **Grouping:** Set `VAR_AGRUPAMENTO` to the column name in your metadata that defines your experimental groups.
* **Taxonomic Rank:** Customize `NETWORK_RANK`, `CORE_MICROBIOME_RANK`, and `INDVAL_RANK` to perform analyses at the level that best fits your study.

## Input Data Structure
To use this script, your data files should follow these structures:

1. **Abundance Table (abundance.tsv):**
   - Tab-separated file.
   - Must contain a column named "tax" with semicolon-separated lineages (Domain;Kingdom;Phylum;Class;Order;Family;Genus;Specie).
2. **Metadata Table (metadados.tsv):**
   - Tab-separated file.
   - The first column must contain the Sample IDs (which must match the column headers of your abundance table).

## Outputs
The script automatically generates high-quality files in your working directory:
* **TIFF Figures:** All plots are saved in TIFF format (600 dpi, LZW compression), ready for manuscript submission.
* **TSV Tables:** Statistical logs, including full DESeq2 results, RDA scores, and Indicator Species tables, are exported for transparency and supplementary material.

## Author
Ângelo Felipe Barbosa de Oliveira
PhD Candidate, Genetics and Molecular Biology (PPGBM)
Federal University of Pará (UFPA)
Biological Engineering Laboratory (ENGBIO)
