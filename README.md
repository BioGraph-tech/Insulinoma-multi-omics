# Data and Code Availability — Insulinoma single-cell multi-omics study

This folder contains the processed data and analysis code needed to
reproduce the results of the manuscript:

> Single-cell multi-omics reveals selective neuroendocrine preservation and
> ribosome-associated program attenuation in human insulinoma


## Raw and public data

- Raw scRNA-seq / scATAC-seq / WES data of the five insulinoma specimens:
- Public normal islet scRNA-seq: HPAP via the Gaulton Lab islet expression
  portal (www.gaultonlab.org/pages/Islet_expression_HPAP/).
- Public normal islet scATAC-seq: GEO **GSE169453**.

These public datasets are not redistributed here; download them from the
sources above.

## Data/

Processed data objects (hg38).

### Data/scRNA/
| File | Description | Used for |
|---|---|---|
| `insulinoma_scRNA_annotated.h5ad` | QC-passed, annotated scRNA-seq object of 5 insulinoma samples (59,595 cells) | Fig. 1b–e |
| `insulinoma_scRNA_inferCNV.h5ad` | scRNA-seq object with transcriptome-inferred CNV profiles and CNV scores | Fig. 1d, Fig. S2a,b |
| `beta_normal_insulinoma_integrated_scored.h5ad` | Integrated HPAP normal beta cells + insulinoma endocrine cells (Harmony), with AUCell program scores (secretory-granule, neuroendocrine/synapse, translation/ribosome, DNA repair) | Fig. 2, Fig. 3, Fig. S4 |

### Data/scATAC/
| File | Description | Used for |
|---|---|---|
| `insulinoma_scATAC_annotated.h5ad` | QC-passed, annotated scATAC-seq object (14,361 nuclei, P4–P5), 500-bp peak matrix | Fig. 1f–h |
| `beta_normal_insulinoma_peak_integrated.h5ad` | Integrated normal beta-cell / insulinoma endocrine-cell cell-by-peak matrix (consensus peaks, MACS3) | Fig. 4 |
| `beta_chromVAR_motif_deviation.h5ad` | chromVAR motif deviation scores per cell | Fig. 5 |
| `beta_motif_match.h5ad` | Binary peak-by-motif match matrix (FigR human motifs) | GRN inference |
| `fragments/260119T_fragments.tsv.gz(.tbi)` | Cell Ranger ATAC fragment file, sample P4 | scATAC reprocessing |
| `fragments/260129T_fragments.tsv.gz(.tbi)` | Cell Ranger ATAC fragment file, sample P5 | scATAC reprocessing |

### Data/GRN/
| File | Description |
|---|---|
| `GRN_final.Rds` | Inferred TF–target edge table (TF, targetgene, log2FoldChange, region support) used for Fig. 5g |
| `tf_go_top10_table.Rds` | Top-10 GO terms per TF regulon |

## Script/

Analysis code, in pipeline order. Python notebooks use Scanpy / SnapATAC2 /
omicverse / infercnvpy; R scripts use Seurat / DESeq2 / clusterProfiler /
igraph.

| File | Step |
|---|---|
| `01_scRNA_QC.ipynb` | scRNA-seq QC (Scanpy/omicverse thresholds, Scrublet) |
| `02_scRNA_dimred_annotation.ipynb` | Normalization, Harmony, UMAP, Leiden, cell-type annotation |
| `03_scRNA_inferCNV.ipynb` | Transcriptome-inferred CNV (infercnvpy) |
| `04_scRNA_integrate_HPAP_normal.ipynb` | Integration with HPAP normal islet scRNA-seq; beta-cell subset |
| `05_scRNA_gene_program_scoring_AUCell.ipynb` | AUCell gene-program scoring |
| `05b_scRNA_cNMF_optional.ipynb` | cNMF program analysis (auxiliary) |
| `06_scATAC_get_normal_islet_data.ipynb` | Download/prepare public normal islet scATAC-seq (GSE169453) |
| `07_scATAC_QC_dimred_annotation.ipynb` | scATAC-seq QC (SnapATAC2), spectral embedding, annotation |
| `08_scATAC_integrate_normal_peak_calling.ipynb` | ATAC integration, MACS3 peak calling, consensus peak matrix |
| `09_scATAC_chromVAR_motif.ipynb` | chromVAR motif activity, differential motif tests, GRN inputs |
| `10_Figure1_plot.R` | Figure 1 panels |
| `11_Figure2_DESeq2_GO.R` | Pseudobulk DESeq2, gene clustering, GO enrichment (Fig. 2, S3) |
| `11b_Figure2_main_heatmap.R` | Figure 2 heatmap assembly |
| `11c_Figure2_gene_scores.R` | Figure 3 program-score plots |
| `12_Figure5_GRN_network.R` | GRN network figure (Fig. 5g); uses Data/GRN/ |
| `assemble_publication_supplementary_figures.R`, `assemble_supplementary_figure_3.R` | Supplementary figure/table assembly |
| `utils/h5ad2rds.R` | AnnData → Seurat conversion (sceasy) for R-based figure scripts |
| `utils/utils.R`, `utils/utils.py` | Shared helper functions |

Note: paths inside the scripts refer to the original working directory
layout (`h5ad/`, `integrate_h5ad/`, `fragment/`, `Rds/`, `table/`); adjust
paths to the files in `Data/` when rerunning.

## Not included (by design)

- Raw FASTQ (scRNA, scATAC, WES): deposited at GSA PRJCA072502.
- Public HPAP / GSE169453 data: available from the original sources.

