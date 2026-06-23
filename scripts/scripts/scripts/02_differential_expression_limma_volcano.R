# Differential protein expression analysis using limma
# Input: log2-transformed tumor proteomics matrix
# Groups: WT_MCA, RAKO_MCA, 15KO_MCA
# Output: DEP tables, volcano plots, and GSEA ranked files

packages <- c(
  "readxl", "dplyr", "tidyr", "stringr", "ggplot2",
  "ggrepel", "openxlsx", "svglite", "scales"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}

library(limma)

input_file <- "data/processed/tumor_proteomics_log2_matrix.xlsx"

output_dep_dir <- "outputs/Differential_expression"
output_volcano_dir <- "outputs/Volcano"
output_gsea_dir <- "outputs/GSEA"

dir.create(output_dep_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_gsea_dir, recursive = TRUE, showWarnings = FALSE)

sheets <- readxl::excel_sheets(input_file)
sheet_name <- sheets[1]

raw <- readxl::read_excel(input_file, sheet = sheet_name)

sample_cols <- c(
  "WT_MCA_1", "WT_MCA_2", "WT_MCA_3", "WT_MCA_4",
  "RAKO_MCA_1", "RAKO_MCA_2", "RAKO_MCA_3", "RAKO_MCA_4",
  "15KO_MCA_1", "15KO_MCA_2", "15KO_MCA_3", "15KO_MCA_4"
)

required_cols <- c("Protein.Names", "Genes", sample_cols)
missing_cols <- setdiff(required_cols, colnames(raw))

if (length(missing_cols) > 0) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}

dat <- raw %>%
  dplyr::select(
    Protein.Names,
    Genes,
    dplyr::all_of(sample_cols)
  )

dat[sample_cols] <- lapply(dat[sample_cols], as.numeric)

sample_info <- data.frame(
  Sample = sample_cols,
  Condition = c(
    rep("WT_MCA", 4),
    rep("RAKO_MCA", 4),
    rep("KO15_MCA", 4)
  ),
  Display_Condition = c(
    rep("WT_MCA", 4),
    rep("RAKO_MCA", 4),
    rep("15KO_MCA", 4)
  ),
  Replicate = c(1:4, 1:4, 1:4),
  stringsAsFactors = FALSE
)

sample_info$Condition <- factor(
  sample_info$Condition,
  levels = c("WT_MCA", "RAKO_MCA", "KO15_MCA")
)

log2_mat <- as.matrix(dat[, sample_cols])
rownames(log2_mat) <- make.unique(as.character(dat$Genes))

design <- model.matrix(~ 0 + Condition, data = sample_info)
colnames(design) <- levels(sample_info$Condition)

contrast_matrix <- makeContrasts(
  KO15_MCA_vs_WT_MCA = KO15_MCA - WT_MCA,
  RAKO_MCA_vs_WT_MCA = RAKO_MCA - WT_MCA,
  RAKO_MCA_vs_KO15_MCA = RAKO_MCA - KO15_MCA,
  levels = design
)

clean_gene <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA
  x <- sapply(strsplit(x, ";|,"), `[`, 1)
  trimws(x)
}

save_volcano <- function(plot, filename, width = 8, height = 7, dpi = 600) {
  ggplot2::ggsave(
    file.path(output_volcano_dir, paste0(filename, ".jpg")),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    file.path(output_volcano_dir, paste0(filename, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )

  ggplot2::ggsave(
    file.path(output_volcano_dir, paste0(filename, ".svg")),
    plot,
    width = width,
    height = height,
    device = svglite::svglite,
    bg = "white"
  )
}

run_dep <- function(contrast_name, display_name, group_a, group_b) {

  group_a_samples <- sample_info$Sample[sample_info$Condition == group_a]
  group_b_samples <- sample_info$Sample[sample_info$Condition == group_b]

  keep <- rowSums(!is.na(log2_mat[, group_a_samples, drop = FALSE])) >= 2 &
    rowSums(!is.na(log2_mat[, group_b_samples, drop = FALSE])) >= 2

  mat_use <- log2_mat[keep, , drop = FALSE]

  fit <- lmFit(mat_use, design)
  fit2 <- contrasts.fit(fit, contrast_matrix[, contrast_name])
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

  res <- topTable(
    fit2,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )

  meta <- dat[keep, c("Protein.Names", "Genes")]
  meta <- meta[match(rownames(res), rownames(mat_use)), ]

  res_final <- cbind(meta, res) %>%
    mutate(
      Contrast = display_name,
      Direction_FDR_0.05 = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "Up",
        adj.P.Val < 0.05 & logFC < 0 ~ "Down",
        TRUE ~ "Not significant"
      ),
      Gene_Label = clean_gene(Genes),
      GSEA_gene = ifelse(is.na(Gene_Label), rownames(res), Gene_Label)
    ) %>%
    relocate(Contrast, .before = Protein.Names)

  gsea_rank <- res_final %>%
    dplyr::select(GSEA_gene, t) %>%
    filter(!is.na(GSEA_gene), !is.na(t)) %>%
    group_by(GSEA_gene) %>%
    slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(t))

  write.table(
    gsea_rank,
    file = file.path(output_gsea_dir, paste0("GSEA_ranked_", display_name, ".rnk")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  write.csv(
    gsea_rank,
    file = file.path(output_gsea_dir, paste0("GSEA_ranked_", display_name, ".csv")),
    row.names = FALSE
  )

  volcano_df <- res_final %>%
    mutate(
      neg_log10_FDR = -log10(adj.P.Val),
      Volcano_Group = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "Up",
        adj.P.Val < 0.05 & logFC < 0 ~ "Down",
        TRUE ~ "Not significant"
      )
    )

  label_df <- bind_rows(
    volcano_df %>%
      filter(adj.P.Val < 0.05, logFC > 0) %>%
      arrange(adj.P.Val, desc(logFC)) %>%
      slice_head(n = 10),

    volcano_df %>%
      filter(adj.P.Val < 0.05, logFC < 0) %>%
      arrange(adj.P.Val, logFC) %>%
      slice_head(n = 10)
  ) %>%
    mutate(Label = ifelse(is.na(Gene_Label), Protein.Names, Gene_Label))

  p <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_FDR)) +
    geom_point(aes(color = Volcano_Group), size = 1.7, alpha = 0.8) +
    scale_color_manual(
      values = c(
        "Up" = "green3",
        "Down" = "#4B2E1A",
        "Not significant" = "grey75"
      )
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 0.5) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Label),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.linewidth = 0.25,
      show.legend = FALSE
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(
      title = paste0("Volcano Plot: ", display_name),
      subtitle = "Green = upregulated; dark brown = downregulated; labels = top FDR < 0.05 proteins",
      x = "log2 fold change",
      y = expression(-log[10]("FDR"))
    )

  save_volcano(p, paste0("Volcano_", display_name))

  return(res_final)
}

res_KO15_vs_WT <- run_dep(
  contrast_name = "KO15_MCA_vs_WT_MCA",
  display_name = "15KO_MCA_vs_WT_MCA",
  group_a = "KO15_MCA",
  group_b = "WT_MCA"
)

res_RAKO_vs_WT <- run_dep(
  contrast_name = "RAKO_MCA_vs_WT_MCA",
  display_name = "RAKO_MCA_vs_WT_MCA",
  group_a = "RAKO_MCA",
  group_b = "WT_MCA"
)

res_RAKO_vs_KO15 <- run_dep(
  contrast_name = "RAKO_MCA_vs_KO15_MCA",
  display_name = "RAKO_MCA_vs_15KO_MCA",
  group_a = "RAKO_MCA",
  group_b = "KO15_MCA"
)

all_results <- list(
  "15KO_MCA_vs_WT_MCA" = res_KO15_vs_WT,
  "RAKO_MCA_vs_WT_MCA" = res_RAKO_vs_WT,
  "RAKO_MCA_vs_15KO_MCA" = res_RAKO_vs_KO15
)

combined_all <- bind_rows(all_results)

dep_count_summary <- bind_rows(lapply(names(all_results), function(nm) {
  x <- all_results[[nm]]

  data.frame(
    Contrast = nm,
    Total_tested = nrow(x),
    FDR_0.05_total = sum(x$adj.P.Val < 0.05, na.rm = TRUE),
    FDR_0.05_up = sum(x$adj.P.Val < 0.05 & x$logFC > 0, na.rm = TRUE),
    FDR_0.05_down = sum(x$adj.P.Val < 0.05 & x$logFC < 0, na.rm = TRUE),
    FDR_0.10_total = sum(x$adj.P.Val < 0.10, na.rm = TRUE),
    FDR_0.10_up = sum(x$adj.P.Val < 0.10 & x$logFC > 0, na.rm = TRUE),
    FDR_0.10_down = sum(x$adj.P.Val < 0.10 & x$logFC < 0, na.rm = TRUE),
    FDR_0.25_total = sum(x$adj.P.Val < 0.25, na.rm = TRUE),
    FDR_0.25_up = sum(x$adj.P.Val < 0.25 & x$logFC > 0, na.rm = TRUE),
    FDR_0.25_down = sum(x$adj.P.Val < 0.25 & x$logFC < 0, na.rm = TRUE)
  )
}))

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Sample_Design")
openxlsx::writeData(wb, "Sample_Design", sample_info)

openxlsx::addWorksheet(wb, "DEP_Count_Summary")
openxlsx::writeData(wb, "DEP_Count_Summary", dep_count_summary)

openxlsx::addWorksheet(wb, "Analysis_Description")
analysis_description <- data.frame(
  Step = c(
    "Input file",
    "Input transformation",
    "Condition naming",
    "Missing-value handling",
    "Filtering",
    "Statistical model",
    "Contrasts",
    "Multiple testing",
    "DEP thresholds",
    "Volcano plots",
    "GSEA ranked files"
  ),
  Description = c(
    input_file,
    "The input matrix was already log2-transformed. No additional log2 transformation was applied.",
    "The condition 15KO_MCA was internally renamed to KO15_MCA only for R model compatibility because model variable names should not start with a number. Output labels remain 15KO_MCA.",
    "Missing values were not globally imputed for differential expression.",
    "For each pairwise comparison, proteins were retained only if detected in at least 2 replicates in each compared condition.",
    "Differential protein expression was tested using limma with empirical Bayes moderation, trend=TRUE and robust=TRUE.",
    "Pairwise contrasts were 15KO_MCA_vs_WT_MCA, RAKO_MCA_vs_WT_MCA, and RAKO_MCA_vs_15KO_MCA.",
    "P values were adjusted using Benjamini-Hochberg FDR.",
    "DEP tables were generated for FDR < 0.05, FDR < 0.10, and FDR < 0.25.",
    "Volcano plots show upregulated proteins in green and downregulated proteins in dark brown. Top FDR < 0.05 proteins were labelled.",
    "Ranked files for GSEA were generated using the moderated t-statistic."
  )
)
openxlsx::writeData(wb, "Analysis_Description", analysis_description)

for (contrast_name in names(all_results)) {
  res <- all_results[[contrast_name]]

  openxlsx::addWorksheet(wb, paste0(contrast_name, "_All"))
  openxlsx::writeData(wb, paste0(contrast_name, "_All"), res)

  openxlsx::addWorksheet(wb, paste0(contrast_name, "_FDR005"))
  openxlsx::writeData(
    wb,
    paste0(contrast_name, "_FDR005"),
    res %>% filter(adj.P.Val < 0.05) %>% arrange(adj.P.Val)
  )

  openxlsx::addWorksheet(wb, paste0(contrast_name, "_FDR010"))
  openxlsx::writeData(
    wb,
    paste0(contrast_name, "_FDR010"),
    res %>% filter(adj.P.Val < 0.10) %>% arrange(adj.P.Val)
  )

  openxlsx::addWorksheet(wb, paste0(contrast_name, "_FDR025"))
  openxlsx::writeData(
    wb,
    paste0(contrast_name, "_FDR025"),
    res %>% filter(adj.P.Val < 0.25) %>% arrange(adj.P.Val)
  )
}

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAD3",
  border = "Bottom"
)

for (s in names(wb)) {
  openxlsx::addStyle(
    wb,
    s,
    header_style,
    rows = 1,
    cols = 1:100,
    gridExpand = TRUE
  )
  openxlsx::freezePane(wb, s, firstRow = TRUE)
  openxlsx::setColWidths(wb, s, cols = 1:100, widths = "auto")
}

excel_out <- file.path(output_dep_dir, "DEP_Analysis_Report_tumor_proteomics.xlsx")
openxlsx::saveWorkbook(wb, excel_out, overwrite = TRUE)

write.csv(
  combined_all,
  file.path(output_dep_dir, "All_DEP_results_all_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  combined_all %>% filter(adj.P.Val < 0.05),
  file.path(output_dep_dir, "All_DEPs_FDR_less_than_0.05.csv"),
  row.names = FALSE
)

write.csv(
  combined_all %>% filter(adj.P.Val < 0.10),
  file.path(output_dep_dir, "All_DEPs_FDR_less_than_0.10.csv"),
  row.names = FALSE
)

write.csv(
  combined_all %>% filter(adj.P.Val < 0.25),
  file.path(output_dep_dir, "All_DEPs_FDR_less_than_0.25.csv"),
  row.names = FALSE
)

message("Differential protein expression analysis completed successfully.")
message("Input matrix was treated as already log2-transformed.")
message("DEP tables saved in: ", output_dep_dir)
message("Volcano plots saved in: ", output_volcano_dir)
message("GSEA ranked files saved in: ", output_gsea_dir)
