# 🧬 16S-NANOPORE: High-Resolution Microbial Ecology Pipeline

This repository contains a modular **R pipeline** developed for the statistical analysis and high-quality visualization of 16S rRNA metagenomic data. This script is specifically optimized for **Oxford Nanopore Technologies (ONT)** long-read data, leveraging nearly full-length 16S sequences (~1580 bp) to achieve superior taxonomic resolution at the genus and species levels.

## 🌟 Key Features

* **Standardized Normalization**: Implements **rarefaction** to the minimum library depth to ensure unbiased diversity comparisons across samples.
* **Statistical Validation**:
    * **Alpha Diversity**: Non-parametric comparisons using the **Kruskal-Wallis test** to evaluate richness and evenness.
    * **Beta Diversity**: Community structure validation via **PERMANOVA (999 permutations)** based on the **Bray-Curtis dissimilarity matrix**.
    * **Biomarker Discovery**: Utilizes **DESeq2** with high-stringency filters (Adjusted P < 0.05 and Log2 Fold Change > 5) to identify robust environmental indicators.
* **Ecological Networks**: Constructs co-occurrence networks based on **Spearman’s rank correlation** ($\rho > 0.6, p < 0.05$) to reveal biotic interactions between key taxa.
* **Core Microbiome Analysis**: Identifies persistent microbial members based on custom prevalence (e.g., 90%) and abundance thresholds.

## 🛠️ Getting Started

### 1. Prerequisites
The script automatically manages dependencies. It utilizes **CRAN** and **Bioconductor** to install and load the necessary packages:
* `phyloseq`, `vegan`, `DESeq2` (Core microbiome and statistical analysis).
* `ggplot2`, `patchwork`, `ggrepel` (Advanced data visualization).
* `igraph`, `ggraph`, `Hmisc` (Network topology and correlation).

### 2. Input Files
The pipeline requires two main tab-separated files in the working directory:
* `ABUND.tsv`: Taxonomic abundance table (must include a `tax` column with the full taxonomic string).
* `DATA.tsv`: Sample metadata (the first column must contain Sample IDs matching the abundance table).

### 3. Configuration (Control Panel)
Global parameters can be adjusted in **Section 3** of the script to ensure reproducibility without modifying the core logic:

```r
VAR_AGRUPAMENTO  <- "Beach_type"  # Experimental variable (e.g., Urban vs. Island)
DESEQ_LFC        <- 5             # Minimum Log2 Fold Change (Stringent: 32x difference)
CORE_PREVALENCIA <- 0.90          # Taxa present in 90% of samples
SEED_GERAL       <- 123           # Global seed for reproducible results
```
## 📊 Analytical Workflow

The pipeline follows a rigorous statistical framework to transform raw taxonomic counts into ecological insights, specifically addressing the high resolution provided by Nanopore long-reads (~1580 bp):

1.  **Quality Assessment**: Validates sequencing depth and sampling sufficiency using **Good's Coverage** ($>95\%$) and **Rarefaction Curves** to ensure the microbial diversity of the sampled urban and island beaches is adequately represented.
2.  **Diversity Profiling**: Calculates a comprehensive suite of alpha diversity indices (Observed OTUs, Chao1, ACE, Shannon, Simpson, and Pielou’s Evenness). Statistical significance between groups is determined using **Kruskal-Wallis non-parametric tests**.
3.  **Community Structure**: Visualizes group separation using **Principal Coordinates Analysis (PCoA)** based on the **Bray-Curtis dissimilarity matrix**. The statistical significance of the environmental clustering is validated via **PERMANOVA** with 999 permutations.
4.  **Differential Abundance (Biomarkers)**: Identifies significant **genomic enrichment** (reflecting genomic potential, not gene expression) of specific genera using **DESeq2**. Highly stringent filters ($|LFC| > 5$ and $p_{adj} < 0.05$) highlight robust biomarkers such as *Pseudomonas* and *Vibrio* in urban areas versus *Synechococcus* in island sites.
5.  **Biotic Interactions**: Maps co-occurrence patterns using **Spearman’s rank correlations** ($\rho > 0.6, p < 0.05$) to identify "Hub" taxa and understand the connectivity of the coastal microbiome under different anthropogenic pressures.

---

## 📝 Citation

If you utilize this pipeline or the associated findings in your research, please cite the main manuscript:

> **under review**

---

## ⚖️ License

Distributed under the **MIT License**. This allows for free use and modification of the scripts for academic and commercial purposes, provided original authorship is credited.

---

*This pipeline was developed by Angelo Felipe Barbosa de Oliveira as part of PhD research in Bioinformatics at the Federal University of Pará (UFPA), Belém, Brazil.*
