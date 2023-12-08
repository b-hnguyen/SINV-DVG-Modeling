# sinv-pbmc-scseq
====
(1) control_filtered_feature_bc_matrix: cellranger output for uninfected sample
* can be fed directly into Seurat for analysis

(2) infected_filtered_feature_bc_matrix: cellranger output for SINV-infected 7dpi sample
* can be fed directly into Seurat for analysis

====
Processing pipeline:
(1) FASTQ files 
(2) Run cellranger count for alignment, filtering, barcode counting etc --> filtered_feature_bc_matrix
(3) Run seurat pipeline in R
