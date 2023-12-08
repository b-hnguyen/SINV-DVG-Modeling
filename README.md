# sinv-pbmc-scseq

This repository contains the processed single-cell RNA sequencing data for the following publication:

In this study, we sequenced PBMCs from uninfected and Sindbis virus-infected mice that were inoculated intracerebrally. Cells from the infected mice were collected at 7 days post-infection.

====

(1) control_filtered_feature_bc_matrix: cellranger output for uninfected sample
* can be fed directly into Seurat for analysis

(2) infected_filtered_feature_bc_matrix: cellranger output for SINV-infected 7dpi sample
* can be fed directly into Seurat for analysis

====
Processing pipeline:
* FASTQ files
* Run cellranger count for alignment, filtering, barcode counting etc --> filtered_feature_bc_matrix
* Run seurat pipeline in R
