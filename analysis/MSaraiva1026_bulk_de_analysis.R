## ============================================================
## BULK RNA-SEQ - MSaraiva1026 - IL1R KO vs B6
## ============================================================

## Load packages
  library(DESeq2)
  library(SummarizedExperiment)
  library(apeglm)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(fgsea)
  library(msigdbr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(matrixStats)
  library(readxl)
  library(openxlsx)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(dplyr)

## 1 - Paths
if (!file.exists(file.path("functions", "utils_funct.R"))) {
  if (file.exists(file.path("..", "functions", "utils_funct.R"))) {
    setwd("..")
  } else {
    stop("Set the working directory to the repository root ")
  }
}

functions_dir <- "functions"

source(file.path(functions_dir, "utils_funct.R"))
source(file.path(functions_dir, "dge_bulk_funct.R"))
source(file.path(functions_dir, "pca_funct.R"))
source(file.path(functions_dir, "kmeans_funct.R"))
source(file.path(functions_dir, "gsea_funct.R"))
source(file.path(functions_dir, "heatmap_funct.R"))


## 2 - Experiment set-up ----

experiment <- "MSaraiva1026"

## Metadata and raw gene counts
input_dir <- "data"
metadata_file <- "Metadata.xlsx"
counts_file <- "Gene_counts.tsv"
output_dir <- "results"

## 2.1 - Samples and groups

## Strain (metadata) -> group code used in DESeq2
group_map <- c("B6" = "B6", "IL1R KO" = "IL1R_KO")
group_order <- c("B6", "IL1R_KO")
group_display <- c("B6" = "B6", "IL1R_KO" = "IL1R KO")

## 2.2 - Colours

## PCA and gene-list heatmaps
group_colors <- c("B6" = "#61a5c2", "IL1R KO" = "#014f86")
pca_ellipse_colors <- c("B6" = "#61a5c2", "IL1R KO" = "#01497c")

## K-means heatmap annotation
kmeans_group_colors <- c("B6" = "#5184B2", "IL1R KO" = "#D8621D")

## Heatmap z-score palette
hm_colors <- c("#304C84", "#5184B2", "#73BCE1", "#FFFFFF", "#D8621D", "#CF382C", "#942924")

## 2.3 - Comparison and DEG criteria

comparison_name <- paste0(experiment, "_IL1RKO_vs_B6")
comparison_display <- "IL1R KO vs B6"

numerator <- "IL1R_KO"
denominator <- "B6"

padj_cut <- 0.05
fold_change_cut <- 1.5

## 2.4 - Analyses to run

run_pca <- TRUE
run_volcano <- TRUE
run_kmeans <- TRUE
run_gsea <- TRUE
run_heatmaps <- TRUE

## K-means cluster annotation: GO-BP, Hallmark and Reactome.
## FALSE runs only clustering, heatmaps and gene lists per cluster.
annotate_kmeans <- TRUE

## 2.5 - Plot and analysis parameters

## PCA (NA = all genes)
pca_top_n <- c(500, 1000, NA)
pca_x_limits <- c(-50, 50)
pca_y_limits <- c(-30, 30)

## K-means
k_values <- 5:15
kmeans_seed <- 42L
kmeans_nstart <- 50L
kmeans_iter_max <- 100L

## GSEA
gsea_fdr <- 0.05
gsea_top_n <- 20L

## Gene lists for heatmaps
am_identity <- c(
  "Flt1", "Fabp4", "Car4", "Krt79", "Pparg",
  "Siglecf", "Mafb", "Maf", "Bhlhe40", "Egr2",
  "Csf2ra", "Csf2rb", "Spi1", "Marco", "Mertk")

self_renewal <- c("Smc1a", "Stag2", "Rad21", "Wapl", "Smc3")


## 3 - Differential gene expression with DESeq2 ----

data <- read_rnaseq_data(
  input_dir = input_dir,
  metadata_file = metadata_file,
  counts_file = counts_file,
  group_map = group_map,
  group_order = group_order)

meta <- data$metadata
meta$group <- relevel(meta$group, ref = denominator)

count_mat <- data$counts
target_tbl <- data$targets

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = meta,
  design = ~ group)

dds <- DESeq2::DESeq(dds)


## 4 - DEG lists and volcano plot ----

## Excel: summary, all_genes, DE_genes, DE_up, DE_down and GeneRanking.
deg <- run_degs(
  dds = dds,
  target_tbl = target_tbl,
  comp_name = comparison_name,
  numerator = numerator,
  denominator = denominator,
  output_dir = file.path(output_dir, "DEGs"),
  padj_cut = padj_cut,
  fold_change_cut = fold_change_cut)

if (run_volcano) {
  volcano <- plot_volcano(
    deg = deg,
    comp_name = comparison_name,
    comparison_display = comparison_display,
    output_dir = file.path(output_dir, "DEGs"),
    padj_cut = padj_cut,
    fold_change_cut = fold_change_cut)
}


## 5 - PCA (rlog, blind = TRUE) ----

if (run_pca) {
  rlog_pca <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))

  pca <- plot_rlog_pca(
    rlog_mat = rlog_pca,
    metadata = meta,
    experiment = experiment,
    group_display = group_display,
    group_cols = group_colors,
    ellipse_cols = pca_ellipse_colors,
    output_dir = file.path(output_dir, "PCA"),
    top_n = pca_top_n,
    x_limits = pca_x_limits,
    y_limits = pca_y_limits)
}


## rlog blind = FALSE for k-means and heatmaps
if (run_kmeans || run_heatmaps) {
  rlog_mat <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = FALSE))
}


## 6 - K-means clustering of DEGs ----

if (run_kmeans) {
  clusters <- run_deg_kmeans(
    deg = deg,
    rlog_mat = rlog_mat,
    metadata = meta,
    comparison_name = comparison_name,
    comparison_display = comparison_display,
    numerator = numerator,
    denominator = denominator,
    group_order = group_order,
    group_display = group_display,
    output_dir = file.path(output_dir, "Kmeans"),
    hm_cols = hm_colors,
    group_cols = kmeans_group_colors,
    k_values_default = k_values,
    kmeans_nstart = kmeans_nstart,
    kmeans_iter_max = kmeans_iter_max,
    seed = kmeans_seed,
    deg_padj_cutoff = padj_cut,
    deg_lfc_cutoff = log2(fold_change_cut),
    run_enrichment = annotate_kmeans)
}


## 7 - GSEA with GO-BP ----

if (run_gsea) {
  gsea <- run_gsea_go_bp(
    deg = deg,
    comparison_name = comparison_name,
    comparison_display = comparison_display,
    output_dir = file.path(output_dir, "GSEA_GO_BP"),
    padj_cutoff = gsea_fdr,
    top_n = gsea_top_n)
}


## 8 - Heatmaps of gene lists ----

if (run_heatmaps) {

  ## 8.1 - AM identity
  heatmap_am <- plot_gene_heatmap(
    rlog_mat = rlog_mat,
    count_mat = count_mat,
    metadata = meta,
    genes = am_identity,
    signature_label = "AM identity",
    file_prefix = "AM_identity",
    group_display = group_display,
    group_cols = group_colors,
    hm_cols = hm_colors,
    output_dir = file.path(output_dir, "Heatmaps", "AM_identity"),
    row_title_size = 9)

  ## 8.2 - Self-renewal
  heatmap_self_renewal <- plot_gene_heatmap(
    rlog_mat = rlog_mat,
    count_mat = count_mat,
    metadata = meta,
    genes = self_renewal,
    signature_label = "Self-renewal",
    file_prefix = "Self_renewal",
    group_display = group_display,
    group_cols = group_colors,
    hm_cols = hm_colors,
    output_dir = file.path(output_dir, "Heatmaps", "Self_renewal"),
    row_title_size = 8)
}


## Session info for reproducibility
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("Analysis complete. Results in: ", normalizePath(output_dir))
