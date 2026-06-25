# Proteomics QC pipeline
# Input: log2-transformed tumor proteomics matrix
# Groups: WT_MCA, RAKO_MCA, 15KO_MCA
# Output: QC plots, PCA plots, distance/correlation heatmaps, and QC report

packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "openxlsx", "svglite", "scales", "grid"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

input_file <- "data/processed/tumor_proteomics_log2_matrix.xlsx"

output_qc_dir <- "outputs/QC"
output_pca_dir <- "outputs/PCA"

dir.create(output_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_pca_dir, recursive = TRUE, showWarnings = FALSE)

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
    rep("15KO_MCA", 4)
  ),
  Replicate = c(1:4, 1:4, 1:4),
  stringsAsFactors = FALSE
)

sample_info$Condition <- factor(
  sample_info$Condition,
  levels = c("WT_MCA", "RAKO_MCA", "15KO_MCA")
)

sample_info$Sample <- factor(sample_info$Sample, levels = sample_cols)

log2_mat <- as.matrix(dat[, sample_cols])
rownames(log2_mat) <- make.unique(as.character(dat$Genes))

detected_summary <- data.frame(
  Sample = sample_cols,
  Condition = sample_info$Condition,
  Replicate = sample_info$Replicate,
  Detected_Proteins = colSums(!is.na(log2_mat)),
  Missing_Values = colSums(is.na(log2_mat)),
  Total_Proteins = nrow(log2_mat),
  Missing_Percentage = colMeans(is.na(log2_mat)) * 100
)

condition_colors <- c(
  "WT_MCA" = "#1B9E77",
  "RAKO_MCA" = "#7570B3",
  "15KO_MCA" = "#D95F02"
)

nature_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

save_all_formats <- function(plot, output_dir, filename, width = 8, height = 6, dpi = 600) {
  ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".jpg")),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )

  ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".svg")),
    plot,
    width = width,
    height = height,
    device = svglite::svglite,
    bg = "white"
  )
}

box_df <- dat %>%
  dplyr::select(dplyr::all_of(sample_cols)) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Log2_Intensity"
  ) %>%
  dplyr::left_join(sample_info, by = "Sample")

box_df$Sample <- factor(box_df$Sample, levels = sample_cols)

p_box <- ggplot(box_df, aes(x = Sample, y = Log2_Intensity, fill = Condition)) +
  geom_boxplot(outlier.size = 0.35, linewidth = 0.35, na.rm = TRUE) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Log2 Protein Intensity Distribution",
    subtitle = "Tumor proteomics matrix; log2-transformed input data",
    x = NULL,
    y = "Log2 intensity"
  ) +
  nature_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_all_formats(p_box, output_qc_dir, "01_log2_intensity_boxplot", 10, 6)

density_df <- box_df %>%
  filter(!is.na(Log2_Intensity))

p_density <- ggplot(density_df, aes(x = Log2_Intensity, color = Condition, group = Sample)) +
  geom_density(linewidth = 0.8, alpha = 0.85, na.rm = TRUE) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "Log2 Protein Intensity Density",
    subtitle = "Distribution of log2 protein abundance values across samples",
    x = "Log2 intensity",
    y = "Density"
  ) +
  nature_theme

save_all_formats(p_density, output_qc_dir, "02_log2_intensity_density", 9, 6)

p_detected <- ggplot(
  detected_summary,
  aes(x = Sample, y = Detected_Proteins, fill = Condition)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  geom_text(aes(label = Detected_Proteins), vjust = -0.4, size = 3.4) +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Number of Detected Proteins per Sample",
    subtitle = "Detected proteins after log2-transformed matrix import",
    x = NULL,
    y = "Detected proteins"
  ) +
  nature_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_all_formats(p_detected, output_qc_dir, "03_detected_proteins_per_sample", 10, 6)

p_missing <- ggplot(
  detected_summary,
  aes(x = Sample, y = Missing_Percentage, fill = Condition)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = paste0(round(Missing_Percentage, 1), "%")),
    vjust = -0.4,
    size = 3.4
  ) +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Missing Values per Sample",
    subtitle = "Percentage of missing protein abundance values",
    x = NULL,
    y = "Missing values"
  ) +
  nature_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_all_formats(p_missing, output_qc_dir, "04_missing_values_per_sample", 10, 6)

impute_row_median <- function(x) {
  med <- median(x, na.rm = TRUE)
  x[is.na(x)] <- med
  x
}

min_detected_samples <- ceiling(ncol(log2_mat) * 0.70)

pca_mat <- log2_mat[
  rowSums(!is.na(log2_mat)) >= min_detected_samples,
  ,
  drop = FALSE
]

pca_mat_imp <- t(apply(pca_mat, 1, impute_row_median))
pca_mat_imp <- pca_mat_imp[
  apply(pca_mat_imp, 1, sd, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

pca_res <- prcomp(t(pca_mat_imp), center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(sample_info, by = "Sample")

pc_var <- round((pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  stat_ellipse(
    aes(group = Condition),
    type = "norm",
    linetype = "dotted",
    linewidth = 0.8,
    level = 0.68,
    show.legend = FALSE
  ) +
  geom_point(size = 4.2, alpha = 0.95) +
  ggrepel::geom_text_repel(
    aes(label = Sample),
    size = 4,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA of Tumor Proteomics Samples",
    subtitle = "Proteins detected in at least 70% of samples; missing values imputed by row median",
    x = paste0("PC1: ", pc_var[1], "% variance"),
    y = paste0("PC2: ", pc_var[2], "% variance")
  ) +
  nature_theme

save_all_formats(p_pca, output_pca_dir, "01_PCA_tumor_proteomics", 8, 7)

cor_mat <- cor(log2_mat, use = "pairwise.complete.obs", method = "pearson")
cor_mat <- cor_mat[sample_cols, sample_cols]

annotation_col <- data.frame(
  Condition = sample_info$Condition
)
rownames(annotation_col) <- sample_info$Sample

ann_colors <- list(Condition = condition_colors)

save_pheatmap_all <- function(ph, output_dir, filename, width = 9, height = 8) {
  pdf(
    file.path(output_dir, paste0(filename, ".pdf")),
    width = width,
    height = height,
    useDingbats = FALSE
  )
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  jpeg(
    file.path(output_dir, paste0(filename, ".jpg")),
    width = width,
    height = height,
    units = "in",
    res = 600,
    quality = 100
  )
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  svglite::svglite(
    file.path(output_dir, paste0(filename, ".svg")),
    width = width,
    height = height
  )
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
}

ph_cor <- pheatmap::pheatmap(
  cor_mat,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  main = "Sample-to-Sample Pearson Correlation",
  silent = TRUE,
  border_color = "grey70"
)

save_pheatmap_all(ph_cor, output_qc_dir, "05_sample_correlation_heatmap", 9, 8)

dist_mat <- as.matrix(dist(t(pca_mat_imp), method = "euclidean"))
dist_mat <- dist_mat[sample_cols, sample_cols]

ph_dist <- pheatmap::pheatmap(
  dist_mat,
  color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.1f",
  fontsize_number = 8,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = ann_colors,
  main = "Sample-to-Sample Euclidean Distance",
  silent = TRUE,
  border_color = "grey70"
)

save_pheatmap_all(ph_dist, output_qc_dir, "06_euclidean_distance_heatmap", 9, 8)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Sample_Info")
openxlsx::writeData(wb, "Sample_Info", sample_info)

openxlsx::addWorksheet(wb, "QC_Summary")
openxlsx::writeData(wb, "QC_Summary", detected_summary)

openxlsx::addWorksheet(wb, "PCA_Coordinates")
openxlsx::writeData(wb, "PCA_Coordinates", pca_df)

openxlsx::addWorksheet(wb, "Correlation_Matrix")
openxlsx::writeData(
  wb,
  "Correlation_Matrix",
  data.frame(Sample = rownames(cor_mat), cor_mat, check.names = FALSE)
)

openxlsx::addWorksheet(wb, "Euclidean_Distance")
openxlsx::writeData(
  wb,
  "Euclidean_Distance",
  data.frame(Sample = rownames(dist_mat), dist_mat, check.names = FALSE)
)

openxlsx::addWorksheet(wb, "Analysis_Notes")
notes <- data.frame(
  Step = c(
    "Input file",
    "Input transformation",
    "Sample design",
    "Missing values",
    "PCA filtering",
    "PCA imputation",
    "Correlation heatmap",
    "Euclidean distance heatmap",
    "Output formats"
  ),
  Description = c(
    input_file,
    "The input matrix was already log2-transformed. No additional log2 transformation was applied in this script.",
    "The final design contains four WT_MCA, four RAKO_MCA, and four 15KO_MCA tumor samples.",
    "Missing values were not imputed for boxplots, detected protein counts, missing-value percentages, or Pearson correlation.",
    "For PCA and Euclidean distance analysis, proteins detected in at least 70% of samples were retained.",
    "Remaining missing values for PCA and Euclidean distance analysis were imputed using protein-wise median values.",
    "Pearson correlation was calculated using pairwise complete observations.",
    "Euclidean distance was calculated from the filtered, log2-transformed, median-imputed matrix.",
    "Figures were exported as high-resolution JPG, PDF, and SVG files."
  )
)

openxlsx::writeData(wb, "Analysis_Notes", notes)

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
    cols = 1:80,
    gridExpand = TRUE
  )
  openxlsx::freezePane(wb, s, firstRow = TRUE)
  openxlsx::setColWidths(wb, s, cols = 1:80, widths = "auto")
}

excel_out <- file.path(output_qc_dir, "QC_Report_tumor_proteomics.xlsx")
openxlsx::saveWorkbook(wb, excel_out, overwrite = TRUE)

message("Proteomics QC analysis completed successfully.")
message("Input matrix was treated as already log2-transformed.")
message("QC outputs saved in: ", output_qc_dir)
message("PCA outputs saved in: ", output_pca_dir)
