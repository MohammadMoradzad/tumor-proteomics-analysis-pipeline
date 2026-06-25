############################################################
# Compare Cell-line vs Tumor GSEA Results
# Databases: Hallmark and KEGG
# Outputs: Excel + scatter plots + dumbbell plots + panels
############################################################

# =========================
# 1) Install/load packages
# =========================
packages <- c(
  "readr", "readxl", "dplyr", "tidyr", "stringr",
  "ggplot2", "ggrepel", "openxlsx", "svglite",
  "patchwork", "forcats", "scales"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  library(p, character.only = TRUE)
}

# =========================
# 2) Manually choose files
# =========================
message("Choose CELL LINE Hallmark GSEA CSV file")
cell_hallmark_file <- file.choose()

message("Choose CELL LINE KEGG GSEA CSV file")
cell_kegg_file <- file.choose()

message("Choose TUMOR Hallmark GSEA Excel file")
tumor_hallmark_file <- file.choose()

message("Choose TUMOR KEGG GSEA Excel file")
tumor_kegg_file <- file.choose()

# =========================
# 3) Output folders
# =========================
out_dir <- "D:/Fjolla_cell line/GSEA_CellLine_vs_Tumor_Comparison"
plot_dir <- file.path(out_dir, "Figures")
excel_dir <- file.path(out_dir, "Excel_results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 4) Helper functions
# =========================
clean_pathway_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "^HALLMARK_", "")
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_squish(x)
  x <- stringr::str_to_title(x)
  x
}

normalize_comparison <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "X15KO", "15KO")
  x <- stringr::str_replace_all(x, "RAKO", "RaKO")
  x <- stringr::str_replace_all(x, "Rako", "RaKO")
  x
}

flip_direction <- function(x) {
  dplyr::case_when(
    x == "Up" ~ "Down",
    x == "Down" ~ "Up",
    TRUE ~ x
  )
}

make_direction_from_nes <- function(nes) {
  dplyr::case_when(
    nes > 0 ~ "Up",
    nes < 0 ~ "Down",
    TRUE ~ "Neutral"
  )
}

extract_fdr_col <- function(df) {
  if ("p.adjust" %in% colnames(df)) return("p.adjust")
  if ("FDR" %in% colnames(df)) return("FDR")
  if ("qvalue" %in% colnames(df)) return("qvalue")
  stop("No FDR column found. Expected p.adjust, FDR, or qvalue.")
}

extract_p_col <- function(df) {
  if ("pvalue" %in% colnames(df)) return("pvalue")
  if ("P_value" %in% colnames(df)) return("P_value")
  if ("p.value" %in% colnames(df)) return("p.value")
  return(NA_character_)
}

extract_gene_col <- function(df) {
  candidates <- c(
    "Core_enrichment_Symbols",
    "Core_enrichment_Genes",
    "core_enrichment_symbols",
    "Gene_Symbols",
    "core_enrichment",
    "Core_enrichment"
  )
  found <- candidates[candidates %in% colnames(df)]
  if (length(found) == 0) return(NA_character_)
  found[1]
}

standardize_gsea_df <- function(df, source, database) {
  
  fdr_col <- extract_fdr_col(df)
  p_col <- extract_p_col(df)
  gene_col <- extract_gene_col(df)
  
  if (!"Description" %in% colnames(df)) {
    if ("Pathway" %in% colnames(df)) df$Description <- df$Pathway
    else stop("No Description or Pathway column found.")
  }
  
  if (!"Comparison" %in% colnames(df)) {
    stop("No Comparison column found.")
  }
  
  out <- df %>%
    dplyr::mutate(
      Source = source,
      Database = database,
      Comparison_raw = as.character(Comparison),
      Comparison = normalize_comparison(Comparison_raw),
      Pathway = clean_pathway_name(Description),
      NES = as.numeric(NES),
      FDR = as.numeric(.data[[fdr_col]]),
      P_value = if (!is.na(p_col)) as.numeric(.data[[p_col]]) else NA_real_,
      Gene_List = if (!is.na(gene_col)) as.character(.data[[gene_col]]) else NA_character_,
      Direction = make_direction_from_nes(NES)
    ) %>%
    dplyr::select(
      Source, Database, Comparison, Comparison_raw,
      Pathway, NES, P_value, FDR, Direction,
      dplyr::everything()
    )
  
  out
}

read_tumor_excel_sheets <- function(file, database) {
  
  sheets <- readxl::excel_sheets(file)
  
  if (database == "Hallmark") {
    use_sheets <- sheets[stringr::str_detect(sheets, "Hallmark$|_Hallmark$")]
  } else {
    use_sheets <- sheets[stringr::str_detect(sheets, "allKEGG$|_allKEGG$")]
  }
  
  if (length(use_sheets) == 0) {
    stop("No usable sheets found for ", database, " in tumor Excel file.")
  }
  
  all <- lapply(use_sheets, function(s) {
    tmp <- readxl::read_excel(file, sheet = s)
    if (!"Comparison" %in% colnames(tmp)) {
      tmp$Comparison <- s
    }
    tmp
  })
  
  dplyr::bind_rows(all)
}

harmonize_contrasts <- function(df) {
  
  df <- df %>%
    dplyr::mutate(
      Target_Comparison = dplyr::case_when(
        Comparison == "15KO_vs_WT" ~ "15KO_vs_WT",
        Comparison == "RaKO_vs_WT" ~ "RaKO_vs_WT",
        Comparison == "15KO_vs_RaKO" ~ "15KO_vs_RaKO",
        Comparison == "RaKO_vs_15KO" ~ "15KO_vs_RaKO",
        TRUE ~ Comparison
      ),
      Needs_flip = Comparison == "RaKO_vs_15KO",
      NES = ifelse(Needs_flip, -NES, NES),
      Direction = make_direction_from_nes(NES),
      Comparison = Target_Comparison
    ) %>%
    dplyr::select(-Target_Comparison, -Needs_flip)
  
  df
}

save_plot_all <- function(plot, filename, width = 9, height = 7, dpi = 600) {
  ggsave(file.path(plot_dir, paste0(filename, ".jpg")),
         plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")),
         plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(plot_dir, paste0(filename, ".svg")),
         plot, width = width, height = height, device = svglite::svglite, bg = "white")
}

# =========================
# 5) Read and standardize files
# =========================
cell_hallmark <- readr::read_csv(cell_hallmark_file, show_col_types = FALSE) %>%
  standardize_gsea_df(source = "Cell line", database = "Hallmark") %>%
  harmonize_contrasts()

cell_kegg <- readr::read_csv(cell_kegg_file, show_col_types = FALSE) %>%
  standardize_gsea_df(source = "Cell line", database = "KEGG") %>%
  harmonize_contrasts()

tumor_hallmark_raw <- read_tumor_excel_sheets(tumor_hallmark_file, database = "Hallmark")
tumor_kegg_raw <- read_tumor_excel_sheets(tumor_kegg_file, database = "KEGG")

tumor_hallmark <- tumor_hallmark_raw %>%
  standardize_gsea_df(source = "Tumor", database = "Hallmark") %>%
  harmonize_contrasts()

tumor_kegg <- tumor_kegg_raw %>%
  standardize_gsea_df(source = "Tumor", database = "KEGG") %>%
  harmonize_contrasts()

# =========================
# 6) Compare common pathways
# =========================
compare_cell_tumor <- function(cell_df, tumor_df, database_name) {
  
  cell_small <- cell_df %>%
    dplyr::select(
      Database, Comparison, Pathway,
      Cell_NES = NES,
      Cell_FDR = FDR,
      Cell_P_value = P_value,
      Cell_Direction = Direction,
      Cell_Gene_List = Gene_List
    )
  
  tumor_small <- tumor_df %>%
    dplyr::select(
      Database, Comparison, Pathway,
      Tumor_NES = NES,
      Tumor_FDR = FDR,
      Tumor_P_value = P_value,
      Tumor_Direction = Direction,
      Tumor_Gene_List = Gene_List
    )
  
  common <- dplyr::inner_join(
    cell_small,
    tumor_small,
    by = c("Database", "Comparison", "Pathway")
  ) %>%
    dplyr::mutate(
      Direction_Relationship = dplyr::case_when(
        Cell_NES > 0 & Tumor_NES > 0 ~ "Same direction",
        Cell_NES < 0 & Tumor_NES < 0 ~ "Same direction",
        Cell_NES * Tumor_NES < 0 ~ "Opposite direction",
        TRUE ~ "Neutral"
      ),
      Concordance_Class = dplyr::case_when(
        Direction_Relationship == "Same direction" & Cell_NES > 0 ~ "Common Up",
        Direction_Relationship == "Same direction" & Cell_NES < 0 ~ "Common Down",
        Direction_Relationship == "Opposite direction" ~ "Opposite",
        TRUE ~ "Neutral"
      ),
      Min_FDR = pmin(Cell_FDR, Tumor_FDR, na.rm = TRUE),
      Max_FDR = pmax(Cell_FDR, Tumor_FDR, na.rm = TRUE),
      Combined_strength = abs(Cell_NES) * abs(Tumor_NES),
      FDR_Label = dplyr::case_when(
        Cell_FDR < 0.05 & Tumor_FDR < 0.05 ~ "**",
        Cell_FDR < 0.10 & Tumor_FDR < 0.10 ~ "#",
        TRUE ~ ""
      )
    ) %>%
    dplyr::arrange(Comparison, Direction_Relationship, dplyr::desc(Combined_strength))
  
  common
}

common_hallmark <- compare_cell_tumor(cell_hallmark, tumor_hallmark, "Hallmark")
common_kegg <- compare_cell_tumor(cell_kegg, tumor_kegg, "KEGG")

common_all <- dplyr::bind_rows(common_hallmark, common_kegg)

# =========================
# 7) Plot functions
# =========================
color_values <- c(
  "Same direction" = "green3",
  "Opposite direction" = "#4B2E1A",
  "Neutral" = "grey70"
)

plot_theme <- theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold")
  )

make_scatter_plot <- function(df, database_name, comparison_name, top_label_n = 12) {
  
  plot_df <- df %>%
    dplyr::filter(Database == database_name, Comparison == comparison_name)
  
  if (nrow(plot_df) == 0) return(NULL)
  
  label_df <- plot_df %>%
    dplyr::filter(Cell_FDR < 0.10 | Tumor_FDR < 0.10) %>%
    dplyr::arrange(dplyr::desc(Combined_strength)) %>%
    dplyr::slice_head(n = top_label_n)
  
  p <- ggplot(plot_df, aes(x = Cell_NES, y = Tumor_NES)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_point(
      aes(
        fill = Direction_Relationship,
        size = Combined_strength
      ),
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.90
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Pathway),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.25,
      segment.linewidth = 0.25,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = color_values) +
    scale_size_continuous(range = c(3, 10)) +
    labs(
      title = paste0(database_name, " GSEA Concordance: ", comparison_name),
      subtitle = "Green = same NES direction; dark brown = opposite NES direction",
      x = "Cell line NES",
      y = "Tumor NES",
      fill = "Direction",
      size = "|Cell NES| × |Tumor NES|"
    ) +
    plot_theme
  
  p
}

make_dumbbell_plot <- function(df, database_name, comparison_name, top_n = 30) {
  
  plot_df <- df %>%
    dplyr::filter(Database == database_name, Comparison == comparison_name) %>%
    dplyr::arrange(dplyr::desc(Combined_strength)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Pathway = forcats::fct_reorder(Pathway, Combined_strength)
    )
  
  if (nrow(plot_df) == 0) return(NULL)
  
  p <- ggplot(plot_df, aes(y = Pathway)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_segment(
      aes(
        x = Cell_NES,
        xend = Tumor_NES,
        yend = Pathway,
        color = Direction_Relationship
      ),
      linewidth = 0.8,
      alpha = 0.85
    ) +
    geom_point(aes(x = Cell_NES), size = 3.2, shape = 21, fill = "white", color = "black") +
    geom_point(aes(x = Tumor_NES, fill = Direction_Relationship), size = 3.8, shape = 21, color = "black") +
    scale_color_manual(values = color_values) +
    scale_fill_manual(values = color_values) +
    labs(
      title = paste0(database_name, " Common Pathways: ", comparison_name),
      subtitle = "Open circle = cell line NES; filled circle = tumor NES",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      color = "Direction",
      fill = "Direction"
    ) +
    plot_theme +
    theme(axis.text.y = element_text(size = 8.5))
  
  p
}

make_heatmap_like_plot <- function(df, database_name, comparison_name, top_n = 30) {
  
  plot_df <- df %>%
    dplyr::filter(Database == database_name, Comparison == comparison_name) %>%
    dplyr::arrange(dplyr::desc(Combined_strength)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::select(Pathway, Cell_NES, Tumor_NES, Direction_Relationship) %>%
    tidyr::pivot_longer(
      cols = c(Cell_NES, Tumor_NES),
      names_to = "Sample_Type",
      values_to = "NES"
    ) %>%
    dplyr::mutate(
      Sample_Type = dplyr::recode(
        Sample_Type,
        "Cell_NES" = "Cell line",
        "Tumor_NES" = "Tumor"
      ),
      Pathway = forcats::fct_reorder(Pathway, NES, .fun = mean)
    )
  
  if (nrow(plot_df) == 0) return(NULL)
  
  p <- ggplot(plot_df, aes(x = Sample_Type, y = Pathway, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(
      title = paste0(database_name, " NES Heatmap: ", comparison_name),
      x = NULL,
      y = NULL,
      fill = "NES"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 8.5),
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  p
}

# =========================
# 8) Generate individual figures
# =========================
comparisons <- c("15KO_vs_WT", "RaKO_vs_WT", "15KO_vs_RaKO")
databases <- c("Hallmark", "KEGG")

scatter_plots <- list()
dumbbell_plots <- list()
heatmap_plots <- list()

for (db in databases) {
  for (comp in comparisons) {
    
    p_scatter <- make_scatter_plot(common_all, db, comp)
    p_dumbbell <- make_dumbbell_plot(common_all, db, comp, top_n = 30)
    p_heat <- make_heatmap_like_plot(common_all, db, comp, top_n = 30)
    
    if (!is.null(p_scatter)) {
      name <- paste0(db, "_Scatter_", comp)
      scatter_plots[[name]] <- p_scatter
      save_plot_all(p_scatter, name, width = 8.5, height = 7)
    }
    
    if (!is.null(p_dumbbell)) {
      name <- paste0(db, "_Dumbbell_", comp)
      dumbbell_plots[[name]] <- p_dumbbell
      save_plot_all(p_dumbbell, name, width = 10, height = 9)
    }
    
    if (!is.null(p_heat)) {
      name <- paste0(db, "_NES_Heatmap_", comp)
      heatmap_plots[[name]] <- p_heat
      save_plot_all(p_heat, name, width = 7, height = 9)
    }
  }
}

# =========================
# 9) Generate panel figures
# =========================
make_scatter_panel <- function(df, database_name) {
  
  plot_df <- df %>%
    dplyr::filter(Database == database_name, Comparison %in% comparisons)
  
  ggplot(plot_df, aes(x = Cell_NES, y = Tumor_NES)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_point(
      aes(fill = Direction_Relationship, size = Combined_strength),
      shape = 21,
      color = "black",
      stroke = 0.30,
      alpha = 0.90
    ) +
    scale_fill_manual(values = color_values) +
    scale_size_continuous(range = c(2.5, 9)) +
    facet_wrap(~ Comparison, nrow = 1) +
    labs(
      title = paste0(database_name, " GSEA Cell Line vs Tumor Concordance"),
      subtitle = "Green = same NES direction; dark brown = opposite NES direction",
      x = "Cell line NES",
      y = "Tumor NES",
      fill = "Direction",
      size = "Concordance strength"
    ) +
    plot_theme
}

p_hallmark_panel <- make_scatter_panel(common_all, "Hallmark")
p_kegg_panel <- make_scatter_panel(common_all, "KEGG")

save_plot_all(p_hallmark_panel, "PANEL_Hallmark_scatter_all_comparisons", width = 15, height = 5.5)
save_plot_all(p_kegg_panel, "PANEL_KEGG_scatter_all_comparisons", width = 15, height = 5.5)

p_combined_panel <- p_hallmark_panel / p_kegg_panel +
  patchwork::plot_annotation(
    title = "Common GSEA Pathway Concordance Between Cell Line and Tumor",
    subtitle = "NES concordance across Hallmark and KEGG pathway databases"
  )

save_plot_all(p_combined_panel, "PANEL_Hallmark_KEGG_combined_scatter", width = 15, height = 11)

# =========================
# 10) Summary tables
# =========================
summary_counts <- common_all %>%
  dplyr::group_by(Database, Comparison, Direction_Relationship) %>%
  dplyr::summarise(
    N_common_pathways = dplyr::n(),
    N_both_FDR_less_0.05 = sum(Cell_FDR < 0.05 & Tumor_FDR < 0.05, na.rm = TRUE),
    N_both_FDR_less_0.10 = sum(Cell_FDR < 0.10 & Tumor_FDR < 0.10, na.rm = TRUE),
    .groups = "drop"
  )

top_conserved <- common_all %>%
  dplyr::filter(Direction_Relationship == "Same direction") %>%
  dplyr::arrange(Database, Comparison, dplyr::desc(Combined_strength))

top_opposite <- common_all %>%
  dplyr::filter(Direction_Relationship == "Opposite direction") %>%
  dplyr::arrange(Database, Comparison, dplyr::desc(Combined_strength))

# =========================
# 11) Excel output
# =========================
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Author_Info")
author_info <- data.frame(
  Field = c(
    "Author",
    "Position",
    "Host lab",
    "Department",
    "Faculty",
    "University",
    "Email",
    "Project",
    "Analysis",
    "Output folder"
  ),
  Information = c(
    "Mohammad Moradzad",
    "PhD student in Immunology",
    "Subburaj Ilangumaran",
    "Department of Immunology and Cell Biology",
    "Faculty of Medicine and Health Sciences",
    "University of Sherbrooke",
    "Mohammad.Moradzad@USherbrooke.ca",
    "Cell line data_proteomics",
    "Comparison of common Hallmark and KEGG GSEA pathways between cell line and tumor samples",
    out_dir
  )
)

openxlsx::writeData(wb, "Author_Info", author_info)

openxlsx::addWorksheet(wb, "Analysis_Description")
description <- data.frame(
  Step = c(
    "Purpose",
    "Databases",
    "Comparisons",
    "Contrast harmonization",
    "Common pathway definition",
    "Same direction",
    "Opposite direction",
    "Concordance strength",
    "Figures"
  ),
  Description = c(
    "Identify common GSEA pathways between cell line and tumor samples.",
    "Hallmark and KEGG were analyzed separately.",
    "The comparisons were 15KO_vs_WT, RaKO_vs_WT, and 15KO_vs_RaKO.",
    "Cell line RaKO_vs_15KO was reversed to 15KO_vs_RaKO by multiplying NES by -1, so it could be compared directly with tumor 15KO_vs_RaKO.",
    "A pathway was considered common when it was present in both cell line and tumor GSEA result tables for the same database and contrast.",
    "Same direction means cell-line NES and tumor NES have the same sign.",
    "Opposite direction means cell-line NES and tumor NES have opposite signs.",
    "Concordance strength was calculated as abs(Cell_NES) multiplied by abs(Tumor_NES).",
    "Scatter, dumbbell, NES heatmap-style plots, and multi-panel figures were generated as JPG, PDF, and SVG."
  )
)

openxlsx::writeData(wb, "Analysis_Description", description)

openxlsx::addWorksheet(wb, "Summary_Counts")
openxlsx::writeData(wb, "Summary_Counts", summary_counts)

openxlsx::addWorksheet(wb, "Common_All")
openxlsx::writeData(wb, "Common_All", common_all)

openxlsx::addWorksheet(wb, "Common_Hallmark")
openxlsx::writeData(wb, "Common_Hallmark", common_hallmark)

openxlsx::addWorksheet(wb, "Common_KEGG")
openxlsx::writeData(wb, "Common_KEGG", common_kegg)

openxlsx::addWorksheet(wb, "Top_Conserved")
openxlsx::writeData(wb, "Top_Conserved", top_conserved)

openxlsx::addWorksheet(wb, "Top_Opposite")
openxlsx::writeData(wb, "Top_Opposite", top_opposite)

for (db in databases) {
  for (comp in comparisons) {
    tmp <- common_all %>%
      dplyr::filter(Database == db, Comparison == comp) %>%
      dplyr::arrange(Direction_Relationship, dplyr::desc(Combined_strength))
    
    sheet_name <- paste0(substr(db, 1, 4), "_", comp)
    sheet_name <- stringr::str_sub(sheet_name, 1, 31)
    
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, tmp)
  }
}

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAD3",
  border = "Bottom"
)

for (s in names(wb)) {
  openxlsx::addStyle(wb, s, header_style, rows = 1, cols = 1:100, gridExpand = TRUE)
  openxlsx::freezePane(wb, s, firstRow = TRUE)
  openxlsx::setColWidths(wb, s, cols = 1:100, widths = "auto")
}

excel_out <- file.path(excel_dir, "CellLine_vs_Tumor_GSEA_Common_Pathway_Comparison.xlsx")
openxlsx::saveWorkbook(wb, excel_out, overwrite = TRUE)

# =========================
# 12) CSV outputs
# =========================
write.csv(common_all, file.path(out_dir, "Common_GSEA_pathways_all.csv"), row.names = FALSE)
write.csv(summary_counts, file.path(out_dir, "Summary_common_pathway_counts.csv"), row.names = FALSE)
write.csv(top_conserved, file.path(out_dir, "Top_conserved_same_direction_pathways.csv"), row.names = FALSE)
write.csv(top_opposite, file.path(out_dir, "Top_opposite_direction_pathways.csv"), row.names = FALSE)

message("Analysis completed successfully.")
message("Excel report saved to: ", excel_out)
message("Figures saved to: ", plot_dir)
