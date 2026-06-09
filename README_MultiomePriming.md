# Multiome Priming: Single-Cell Regulatory Priming Analysis

## Overview

**Multiome Priming** is a user-friendly **Shiny application** and computational framework for analyzing **single-cell multiome RNA + ATAC datasets**. The method is designed to identify **regulatory priming states** by measuring the degree of discordance between transcriptional activity and chromatin accessibility at the single-cell level.

The tool introduces the **Regulatory Difference Score (RDS)**, a cell-level score derived from RNA–ATAC decoupling:

```r
Regulatory Difference Score = 1 - RNA_ATAC_Correlation
```

Cells with high RNA–ATAC concordance are interpreted as transcriptionally and epigenetically coordinated, while cells with high regulatory difference scores are interpreted as potentially **primed**, transitional, or regulatory-decoupled states.

Unlike RNA-only approaches, Multiome Priming directly incorporates chromatin accessibility information, allowing detection of cells where chromatin may be open before full transcriptional activation.

---

## Key Features

- **Single-cell RNA + ATAC integration** using matched multiome data
- **Cell-wise RNA–ATAC correlation analysis**
- **Regulatory Difference Score calculation**
- **High-, Middle-, and Low-Priming cell classification**
- **Gene-level regulatory decoupling analysis**
- **Top primed genes and cells identification**
- **Hallmark pathway enrichment analysis**
- **Pathway robustness validation across multiple gene thresholds**
- **Cell-type enrichment analysis of priming states**
- **Chi-square and odds-ratio statistical enrichment testing**
- **External validation on independent multiome datasets**
- **Publication-ready plots exported as PNG at 1000 dpi**
- **CSV output files for downstream analysis**

---

## Method Summary

Multiome Priming quantifies regulatory priming by comparing RNA expression and ATAC-derived gene activity for the same cells and genes.

For each cell:

1. Shared genes between RNA and ATAC gene activity matrices are identified.
2. RNA expression and ATAC accessibility vectors are extracted.
3. RNA–ATAC correlation is calculated.
4. Regulatory Difference Score is computed as:

```r
RDS = 1 - RNA_ATAC_Correlation
```

Interpretation:

| Score Pattern | Biological Interpretation |
|--------------|----------------------------|
| Low RDS | RNA and ATAC are coordinated |
| Intermediate RDS | Transitional regulatory state |
| High RDS | RNA–ATAC decoupling / potential regulatory priming |

Cells are then grouped into:

- **High_Priming**
- **Middle**
- **Low_Priming**

based on score distribution cutoffs.

---

## Input Requirements

The Shiny application expects a prepared `.rds` object containing matched RNA and ATAC gene-level matrices.

Required object structure:

```r
list(
  RNA = RNA_matrix,        # genes × cells
  ATAC = ATAC_matrix,      # genes × cells; gene activity matrix
  genes = shared_genes,    # character vector
  cells = shared_cells     # character vector
)
```

Important requirements:

- RNA and ATAC matrices must contain matched cell names.
- RNA and ATAC matrices must contain matched gene names.
- ATAC should be represented as a gene-level accessibility or gene activity matrix, not raw peak-level counts.
- Sparse matrices are recommended for large datasets.

---

## Preparing Input Data

### Example: Creating an App-Ready Multiome Object

```r
rna_mat <- your_rna_matrix
atac_mat <- your_gene_activity_matrix

common_genes <- intersect(rownames(rna_mat), rownames(atac_mat))
common_cells <- intersect(colnames(rna_mat), colnames(atac_mat))

rna_mat <- rna_mat[common_genes, common_cells]
atac_mat <- atac_mat[common_genes, common_cells]

multiome_object <- list(
  RNA = rna_mat,
  ATAC = atac_mat,
  genes = common_genes,
  cells = common_cells
)

saveRDS(multiome_object, "MultiomePriming_input.rds")
```

---

## Application Modules

### 1. Multiome Priming Analysis

This module calculates:

- Cell-wise RNA–ATAC correlation
- Regulatory Difference Score
- High-, Middle-, and Low-Priming groups
- Gene-level RNA–ATAC correlation
- Gene-level Regulatory Difference Score

Main outputs:

- `Multiome_Cell_Priming_Scores.csv`
- `Multiome_Gene_Priming_Scores.csv`

---

### 2. Multiome Priming Visualization

This module generates publication-ready figures, including:

- Cell Priming UMAP
- Cell Priming Score Distribution
- Top 20 Primed Cells
- Top 20 Primed Genes
- High vs Low Priming Boxplot
- Priming Group Counts
- Primed Pathway Enrichment Plot

All figures can be downloaded as high-resolution PNG files at **1000 dpi**.

---

### 3. Pathway Enrichment Analysis

Gene-level regulatory difference scores are used to identify biological programs associated with high regulatory decoupling.

The pathway analysis supports:

- MSigDB Hallmark pathway analysis
- Top primed gene enrichment
- Pathway robustness testing across gene-selection thresholds

In the project analysis, recurrent pathways included:

- TNFα signaling via NF-κB
- P53 pathway
- Oxidative phosphorylation
- Heme metabolism
- Apical junction
- DNA repair
- Hypoxia

---

### 4. Robustness Validation

To test whether biological conclusions were stable across different gene-selection thresholds, pathway enrichment was performed using multiple top-gene cutoffs.

Thresholds tested:

- Top 10% primed genes
- Top 20% primed genes
- Top 30% primed genes

Recurring pathway programs across thresholds supported the stability of the Multiome Priming framework.

---

### 5. Cell-Type Enrichment Analysis

Cell-type annotations can be integrated with cell-level priming scores to determine whether High-Priming cells are enriched in specific biological populations.

Statistical testing includes:

- Pearson's Chi-square test
- Fisher's exact test
- Odds-ratio enrichment analysis

In the discovery analysis, High-Priming cells were enriched in immune-responsive populations such as:

- CD8 T cells
- NK cells
- Macrophages
- Plasma-like IRF4 cells

Low-Priming states were enriched in more stable or less regulatory-decoupled populations such as:

- Naive CD4 T cells
- Endothelial cells
- Monocytes
- Plasma cells

---

## Validation Strategy

Multiome Priming was evaluated through multiple validation layers:

1. **Discovery cohort analysis** using a lung cancer multiome dataset
2. **Pathway enrichment** of highly primed genes
3. **Pathway robustness analysis** across multiple gene-selection thresholds
4. **Cell-type enrichment analysis** using annotated cellular populations
5. **Statistical enrichment testing** using Chi-square and odds-ratio analyses
6. **External validation** in an independent lung adenocarcinoma dataset
7. **Cross-cancer validation** in a human ATRT multiome dataset

These analyses support the robustness, reproducibility, and portability of the Multiome Priming framework.

---

## External Validation Datasets

### Independent Lung Adenocarcinoma Validation

An independent human lung adenocarcinoma multiome dataset was used for external validation.

Processing steps included:

- Reconstruction of RNA and ATAC matrices from Cell Ranger outputs
- Generation of ATAC gene activity scores using Signac and GENCODE gene annotations
- Feature harmonization between RNA and ATAC matrices
- Construction of an app-ready Multiome Priming input object
- Analysis using the same Shiny application workflow

This validation confirmed that Multiome Priming can be applied to an independent lung cancer cohort.

---

### Cross-Cancer ATRT Validation

A human atypical teratoid/rhabdoid tumor (ATRT) multiome dataset was used for cross-cancer validation.

Processing steps included:

- Loading matched RNA and ATAC peak-count matrices
- Converting ATAC peaks to gene-level accessibility by nearest-gene assignment
- Aggregating peaks assigned to the same gene
- Matching RNA and ATAC gene-level matrices
- Creating a reduced app-ready object using informative genes while preserving all cells
- Running the full Multiome Priming workflow in the Shiny app

This analysis demonstrated that the framework can generalize beyond lung cancer to a biologically distinct pediatric brain tumor context.

---

## Comparison with Existing Tools

| Method | Main Purpose | Detects Priming States | Requires Trajectory Inference | Key Output |
|--------|--------------|------------------------|-------------------------------|------------|
| Multiome Priming | Quantify regulatory discordance between RNA and ATAC | Yes | No | Regulatory Difference Score and priming-state classification |
| TREASMO | Identify RNA–ATAC correlations | No | No | Gene–peak regulatory correlations |
| ArchR | Analyze chromatin accessibility and regulatory networks | No | Optional | Gene activity scores and peak-to-gene links |
| Signac | Single-cell ATAC-seq analysis and integration | No | Optional | Chromatin accessibility and gene activity profiles |
| Chromatin Potential | Predict future cellular states | Indirectly | Yes | Predicted cell-state trajectories |

Existing multiomic frameworks primarily focus on regulatory correlation analysis, chromatin accessibility profiling, or trajectory inference. In contrast, Multiome Priming directly quantifies regulatory discordance between transcriptional output and chromatin accessibility at the single-cell level. This enables the identification of primed cellular states without requiring pseudotime analysis, lineage reconstruction, or prior biological knowledge.

---

## Outputs

### CSV Outputs

- `Multiome_Cell_Priming_Scores.csv`
- `Multiome_Gene_Priming_Scores.csv`
- Pathway enrichment tables
- Robustness analysis tables
- Cell-type enrichment tables
- Odds-ratio statistics

### Figure Outputs

- Cell Priming UMAP
- Cell Priming Score Distribution
- Top Primed Cells
- Top Primed Genes
- High vs Low Priming Boxplot
- Priming Group Counts
- Pathway Enrichment Plot
- Pathway Robustness Plot
- Cell-Type Enrichment Plot
- Odds-Ratio Bubble Plot

Figures are exported in high-resolution PNG format at **1000 dpi**.

---

## Quick Start

### Run the Shiny App

```r
shiny::runApp("path_to_MultiomePriming_app_folder")
```

### Upload Data

Upload a prepared `.rds` file with the required list structure:

```r
list(
  RNA = RNA_matrix,
  ATAC = ATAC_matrix,
  genes = shared_genes,
  cells = shared_cells
)
```

### Run Analysis

1. Upload the multiome object.
2. Click **Run Multiome Priming Analysis**.
3. Download cell and gene score tables.
4. Generate visualizations.
5. Download figures as 1000 dpi PNG files.

---

## Recommended Citation Statement

If you use Multiome Priming in your work, please cite the associated manuscript or repository once available.

Suggested description:

> Multiome Priming quantifies cell-wise regulatory discordance between RNA expression and chromatin accessibility to identify primed cellular states in single-cell multiome datasets.

---

## License

MIT License — free for academic and commercial use with attribution.

---

## Contact

- **Lead developers:** Ali Haider (alihaider32@gcuf.edu.pk)
- For questions, collaborations, or feature requests, please open a GitHub Issue in this repository.
