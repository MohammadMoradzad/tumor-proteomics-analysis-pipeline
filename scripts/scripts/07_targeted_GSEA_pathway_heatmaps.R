# 07_targeted_GSEA_pathway_heatmaps.R
# Targeted heatmaps from GSEA pathway genes
# Project: IL-15 proteomics analysis
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(openxlsx)

# 1. Set paths

expr_file <- "data/processed/tumor_proteomics_imputed_matrix.xlsx"
gsea_file <- "outputs/tables/pathway_analysis/Hallmark_results.xlsx"

out_dir <- "outputs/figures/heatmaps/targeted_GSEA_pathways"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

table_out_dir <- "outputs/tables/pathway_analysis/targeted_heatmap_gene_tables"
dir.create(table_out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Read expression matrix

expr <- read_excel(expr_file)

# Change this if your gene/protein column has another name
gene_col <- "Gene"

if (!gene_col %in% colnames(expr)) {
  stop("Gene column not found. Please check the gene_col name.")
}

expr <- expr %>%
  filter(!is.na(.data[[gene_col]])) %>%
  distinct(.data[[gene_col]], .keep_all = TRUE)

expr_mat <- expr %>%
  column_to_rownames(gene_col) %>%
  as.data.frame()

expr_mat[] <- lapply(expr_mat, as.numeric)
expr_mat <- as.matrix(expr_mat)

# 3. Sample annotation
sample_annotation <- data.frame(
  Sample = colnames(expr_mat),
  Condition = case_when(
    str_detect(colnames(expr_mat), regex("WT|wildtype|wild_type", ignore_case = TRUE)) ~ "WT",
    str_detect(colnames(expr_mat), regex("IL15", ignore_case = TRUE)) ~ "IL15",
    str_detect(colnames(expr_mat), regex("KO|knockout|IL15RaKO|IL15R", ignore_case = TRUE)) ~ "KO",
    TRUE ~ "Unknown"
  )
)

rownames(sample_annotation) <- sample_annotation$Sample
sample_annotation$Sample <- NULL

# 4. Read GSEA results

gsea <- read_excel(gsea_file)

# Try to automatically detect pathway and leading-edge/core-enrichment columns
pathway_col <- intersect(
  c("pathway", "Pathway", "NAME", "Description", "ID"),
  colnames(gsea)
)[1]

gene_list_col <- intersect(
  c("core_enrichment", "leadingEdge", "leading_edge", "genes", "Genes"),
  colnames(gsea)
)[1]

if (is.na(pathway_col)) {
  stop("Pathway/name column not found in GSEA file.")
}

if (is.na(gene_list_col)) {
  stop("Gene-list column not found. Expected core_enrichment, leadingEdge, or genes.")
}

# 5. Select targeted interferon pathways

target_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

gsea_target <- gsea %>%
  filter(.data[[pathway_col]] %in% target_pathways)

if (nrow(gsea_target) == 0) {
  stop("None of the target pathways were found in the GSEA file.")
}
# 6. Function to extract genes from GSEA gene column

extract_genes <- function(x) {
  x <- as.character(x)
  genes <- unlist(str_split(x, "/|;|,|\\s+"))
  genes <- genes[genes != ""]
  unique(genes)
}

# 7. Generate heatmap for each pathway

for (i in seq_len(nrow(gsea_target))) {
  
  pathway_name <- gsea_target[[pathway_col]][i]
  pathway_genes <- extract_genes(gsea_target[[gene_list_col]][i])
  
  matched_genes <- intersect(pathway_genes, rownames(expr_mat))
  
  if (length(matched_genes) < 2) {
    warning(paste("Too few matched genes for:", pathway_name))
    next
  }
  
  heatmap_mat <- expr_mat[matched_genes, , drop = FALSE]
  
  # Z-score by protein/gene
  heatmap_z <- t(scale(t(heatmap_mat)))
  heatmap_z[is.na(heatmap_z)] <- 0
  
  safe_name <- gsub("[^A-Za-z0-9_]", "_", pathway_name)
 
  # Save matched gene table
  
  gene_table <- data.frame(
    Pathway = pathway_name,
    Gene = matched_genes
  )
  
  write.xlsx(
    gene_table,
    file = file.path(table_out_dir, paste0(safe_name, "_matched_genes.xlsx")),
    overwrite = TRUE
  )
  
  # Save heatmap PDF
  
  pdf(
    file = file.path(out_dir, paste0(safe_name, "_heatmap.pdf")),
    width = 8,
    height = max(5, length(matched_genes) * 0.18)
  )
  
  pheatmap(
    heatmap_z,
    annotation_col = sample_annotation,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 7,
    fontsize_col = 9,
    main = pathway_name,
    border_color = NA
  )
  
  dev.off()
}

# End of script
