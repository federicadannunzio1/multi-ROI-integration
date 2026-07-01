# =============================================================================
# 07_visualization.R
# Visualizzazioni finali publication-ready ed export matrici
#
# Input:  output/data/04_seurat_clustered.rds
# Output: output/plots/Final_*.pdf/png
#         output/data/07_integrated_matrix_arcsinh.tsv
#
# NOTA UMAP:
#   DimPlot di Seurat auto-rasterizza con >100k cellule ma la qualita'
#   del raster dipende dal backend grafico. Si usa ggplot manuale su
#   subsample per garantire plot corretti indipendentemente dall'ambiente.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  library(pheatmap)
  library(scales)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== STEP 7: Visualizzazioni finali ed export ===")

seurat_obj  <- readRDS(file.path(OUT_DATA, "04_seurat_clustered.rds"))
SHORT_NAMES <- rownames(seurat_obj)
n_clusters  <- length(unique(seurat_obj$cluster))
COFACTOR    <- seurat_obj@misc$cofactor

# Palette dinamiche
all_patients <- sort(unique(seurat_obj$patient_id))
g3_patients  <- sort(unique(seurat_obj$patient_id[seurat_obj$group == "G3"]))
shh_patients <- sort(unique(seurat_obj$patient_id[seurat_obj$group == "SHH"]))
palette_patient <- c(
  setNames(PALETTE_PATIENTS_G3[seq_along(g3_patients)],   g3_patients),
  setNames(PALETTE_PATIENTS_SHH[seq_along(shh_patients)], shh_patients)
)
palette_cluster <- colorRampPalette(brewer.pal(9, "Set3"))(n_clusters)

message(sprintf("Cellule: %d | Cluster: %d | Cofactor: %d",
                ncol(seurat_obj), n_clusters, COFACTOR))

mat_all <- GetAssayData(seurat_obj, layer = "data", assay = "MICS")

# --------------------------------------------------------------------------
# STEP 7a — UMAP publication-ready (ggplot manuale su subsample)
# --------------------------------------------------------------------------

set.seed(SEED)
N_PLOT   <- min(100000, ncol(seurat_obj))
idx_plot <- sample(seq_len(ncol(seurat_obj)), N_PLOT)

umap_df            <- as.data.frame(Embeddings(seurat_obj, "umap")[idx_plot, ])
colnames(umap_df)  <- c("UMAP1", "UMAP2")
umap_df$cluster    <- factor(seurat_obj$cluster[idx_plot])
umap_df$group      <- seurat_obj$group[idx_plot]
umap_df$patient_id <- seurat_obj$patient_id[idx_plot]

cluster_centers <- umap_df %>%
  group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

p_cluster <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
  geom_point(size = 0.1, alpha = 0.2) +
  geom_label(data = cluster_centers, aes(label = cluster),
             size = 4, fontface = "bold", colour = "black",
             fill = "white", alpha = 0.8, label.size = 0.3) +
  scale_colour_manual(values = palette_cluster) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Cluster") +
  theme_classic(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

p_group <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = group)) +
  geom_point(size = 0.1, alpha = 0.2) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Group") +
  theme_classic(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

p_patient <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = patient_id)) +
  geom_point(size = 0.1, alpha = 0.2) +
  scale_colour_manual(values = palette_patient) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Patient") +
  theme_classic(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggsave(file.path(OUT_PLOTS, "Final_01_UMAP_cluster_group.pdf"),
       p_cluster + p_group, width = 14, height = 6)
ggsave(file.path(OUT_PLOTS, "Final_01_UMAP_cluster_group.png"),
       p_cluster + p_group, width = 14, height = 6, dpi = 300)
message("  Salvato: Final_01_UMAP_cluster_group.pdf/png")

ggsave(file.path(OUT_PLOTS, "Final_01b_UMAP_patient.pdf"),
       p_cluster + p_patient, width = 14, height = 6)
message("  Salvato: Final_01b_UMAP_patient.pdf")

# UMAP split per gruppo
p_split_list <- lapply(sort(unique(seurat_obj$group)), function(grp) {
  sub <- umap_df[umap_df$group == grp, ]
  ggplot(sub, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 0.08, alpha = 0.2) +
    scale_colour_manual(values = palette_cluster) +
    labs(title = grp, x = "UMAP 1", y = "UMAP 2") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
})
ggsave(file.path(OUT_PLOTS, "Final_02_UMAP_split_group.pdf"),
       wrap_plots(p_split_list, nrow = 1),
       width = 6 * length(p_split_list), height = 5)
message("  Salvato: Final_02_UMAP_split_group.pdf")

# --------------------------------------------------------------------------
# STEP 7b — Dot plot
# --------------------------------------------------------------------------

p_dot <- DotPlot(seurat_obj, features = SHORT_NAMES, group.by = "cluster",
  assay = "MICS", scale = TRUE, col.min = -2.5, col.max = 2.5, dot.scale = 6) +
  scale_colour_gradient2(low = "#1A5276", mid = "white", high = "#C0392B",
                         midpoint = 0) +
  labs(title = "Dot plot by cluster", x = "Protein", y = "Cluster") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_PLOTS, "Final_03_dotplot_clusters.pdf"),
       p_dot, width = 9, height = max(5, n_clusters * 0.6 + 2))
message("  Salvato: Final_03_dotplot_clusters.pdf")

# --------------------------------------------------------------------------
# STEP 7c — Heatmap single-cell (subsample)
# --------------------------------------------------------------------------

set.seed(SEED)
n_per_cluster <- 100
cells_heatmap <- lapply(levels(seurat_obj$cluster), function(cl) {
  all_cl <- WhichCells(seurat_obj, idents = cl)
  sample(all_cl, min(n_per_cluster, length(all_cl)))
}) |> unlist()

mat_sample <- mat_all[, cells_heatmap]
ann_col <- data.frame(
  Cluster = seurat_obj$cluster[cells_heatmap],
  Group   = seurat_obj$group[cells_heatmap],
  row.names = cells_heatmap
)
ann_colors <- list(
  Cluster = setNames(palette_cluster, levels(seurat_obj$cluster)),
  Group   = PALETTE_GROUP
)

pheatmap(mat_sample,
  annotation_col  = ann_col,
  annotation_colors = ann_colors,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_colnames   = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  scale = "row", fontsize_row = 10,
  main = paste0("Single-cell heatmap (", n_per_cluster, " cells per cluster, Z-score per protein)"),
  filename = file.path(OUT_PLOTS, "Final_04_heatmap_singlecell.pdf"),
  width = 12, height = 5
)
message("  Salvato: Final_04_heatmap_singlecell.pdf")

# --------------------------------------------------------------------------
# STEP 7d — Composizione cluster per paziente: barplot stacked + tabella
# --------------------------------------------------------------------------

message("Composizione cluster per paziente...")

meta_df <- as.data.frame(seurat_obj@meta.data[, c("patient_id", "group", "cluster")])

# Totali marginali (paziente e cluster) — necessari per le percentuali
tot_patient <- meta_df %>%
  group_by(patient_id) %>%
  summarise(total_cells_patient = n(), .groups = "drop")

tot_cluster <- meta_df %>%
  group_by(cluster) %>%
  summarise(total_cells_cluster = n(), .groups = "drop")

# ── Tabella principale: una riga per ogni combinazione paziente x cluster ──

comp_table <- meta_df %>%
  group_by(cluster, patient_id, group) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  left_join(tot_patient, by = "patient_id") %>%
  left_join(tot_cluster, by = "cluster") %>%
  mutate(
    pct_of_cluster = round(100 * n_cells / total_cells_cluster, 2),
    pct_of_patient = round(100 * n_cells / total_cells_patient, 2)
  ) %>%
  arrange(cluster, group, patient_id)

# ── Subtotali per gruppo all'interno di ogni cluster ──

group_subtotals <- comp_table %>%
  group_by(cluster, group) %>%
  summarise(
    patient_id           = paste0("[subtotal_", first(group), "]"),
    n_cells              = sum(n_cells),
    total_cells_patient  = NA_real_,
    total_cells_cluster  = first(total_cells_cluster),
    pct_of_cluster       = round(sum(pct_of_cluster), 2),
    pct_of_patient       = NA_real_,
    .groups = "drop"
  )

# ── Total per cluster ──

cluster_totals <- comp_table %>%
  group_by(cluster) %>%
  summarise(
    group                = "[total]",
    patient_id           = "[total]",
    n_cells              = sum(n_cells),
    total_cells_patient  = NA_real_,
    total_cells_cluster  = first(total_cells_cluster),
    pct_of_cluster       = 100,
    pct_of_patient       = NA_real_,
    .groups = "drop"
  )

# ── TABLE 1 (long): one row per patient x cluster, with subtotals ──

comp_full <- bind_rows(comp_table, group_subtotals, cluster_totals) %>%
  arrange(cluster,
          !grepl("^\\[", patient_id),   # patients before subtotals/totals
          group, patient_id)

write.csv(comp_full,
          file.path(OUT_DATA, "07_cluster_composition_long.csv"),
          row.names = FALSE, quote = FALSE)
message("  Salvato: 07_cluster_composition_long.csv")

# ── TABLE 2 (wide): rows = cluster, cols = patients + group subtotals ──
# Easiest to read at a glance.

cells_wide <- pivot_wider(
  comp_table[, c("cluster", "patient_id", "n_cells")],
  names_from = "patient_id", values_from = "n_cells", values_fill = 0L
)
# Aggiungi colonne G3_total, SHH_total e grand total
cells_wide$G3_total  <- rowSums(cells_wide[, g3_patients,  drop = FALSE])
cells_wide$SHH_total <- rowSums(cells_wide[, shh_patients, drop = FALSE])
cells_wide$grand_total <- cells_wide$G3_total + cells_wide$SHH_total
cells_wide$pct_G3    <- round(100 * cells_wide$G3_total  / cells_wide$grand_total, 1)
cells_wide$pct_SHH   <- round(100 * cells_wide$SHH_total / cells_wide$grand_total, 1)

write.csv(cells_wide,
          file.path(OUT_DATA, "07_cluster_composition_wide.csv"),
          row.names = FALSE, quote = FALSE)
message("  Salvato: 07_cluster_composition_wide.csv")

# ── TABLE 3 (patient summary): rows = patient, cols = clusters ──
# How each patient distributes across clusters.

patient_wide <- pivot_wider(
  comp_table[, c("patient_id", "group", "cluster", "n_cells", "pct_of_patient")],
  names_from  = "cluster",
  values_from = c("n_cells", "pct_of_patient"),
  values_fill = 0
)
# Add patient total
patient_totals <- comp_table %>%
  group_by(patient_id, group) %>%
  summarise(total_cells = sum(n_cells), .groups = "drop")
patient_wide <- left_join(patient_wide, patient_totals, by = c("patient_id", "group")) %>%
  arrange(group, patient_id)

write.csv(patient_wide,
          file.path(OUT_DATA, "07_patient_cluster_distribution.csv"),
          row.names = FALSE, quote = FALSE)
message("  Salvato: 07_patient_cluster_distribution.csv")

# Print log
message("\nCells per cluster and patient (wide):")
print(as.data.frame(cells_wide[, c("cluster", g3_patients, shh_patients,
                                   "G3_total", "SHH_total", "grand_total",
                                   "pct_G3", "pct_SHH")]))

message("\n% of cluster occupied per patient:")
pct_wide <- pivot_wider(
  comp_table[, c("cluster", "patient_id", "pct_of_cluster")],
  names_from = "patient_id", values_from = "pct_of_cluster", values_fill = 0
)
print(as.data.frame(pct_wide))

# ── Ordine numerico dei cluster per i plot ──

cluster_levels <- sort(unique(as.integer(as.character(comp_table$cluster))))
comp_table$cluster <- factor(comp_table$cluster, levels = cluster_levels)

# ── Plot A: barplot stacked assoluto (n cellule) ──

p_stack_abs <- ggplot(comp_table,
    aes(x = cluster, y = n_cells, fill = patient_id)) +
  geom_col(position = "stack", width = 0.8, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = palette_patient) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Cluster composition by patient — cell counts",
    subtitle = sprintf("Total: %s cells | %d clusters | %d patients (%d G3, %d SHH)",
                       format(ncol(seurat_obj), big.mark = ","),
                       n_clusters,
                       length(all_patients),
                       length(g3_patients),
                       length(shh_patients)),
    x = "Cluster", y = "Number of cells", fill = "Patient"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.major.x = element_blank())

# ── Plot B: barplot stacked proporzionale (% per cluster) ──

p_stack_pct <- ggplot(comp_table,
    aes(x = cluster, y = pct_of_cluster, fill = patient_id)) +
  geom_col(position = "stack", width = 0.8, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = palette_patient) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  labs(
    title    = "Cluster composition by patient — proportions",
    subtitle = "Each bar sums to 100%: allows comparison across clusters of different sizes",
    x = "Cluster", y = "% of cells in cluster", fill = "Patient"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.major.x = element_blank())

ggsave(file.path(OUT_PLOTS, "Final_05_stacked_barplot_cluster_composition.pdf"),
       p_stack_abs / p_stack_pct,
       width = max(10, n_clusters * 0.9 + 3), height = 14)
message("  Salvato: Final_05_stacked_barplot_cluster_composition.pdf")

# --------------------------------------------------------------------------
# STEP 7e — Export matrice integrata completa
# --------------------------------------------------------------------------

message("Export matrice completa...")

mat_export  <- t(as.matrix(mat_all))
umap_coords <- Embeddings(seurat_obj, "umap")

export_df <- data.frame(
  cell_id    = rownames(mat_export),
  patient_id = seurat_obj$patient_id,
  group      = seurat_obj$group,
  roi_id     = seurat_obj$roi_id,
  cluster    = seurat_obj$cluster,
  UMAP_1     = round(umap_coords[, 1], 4),
  UMAP_2     = round(umap_coords[, 2], 4),
  round(mat_export, 4),
  check.names = FALSE
)

out_tsv <- file.path(OUT_DATA, "07_integrated_matrix_arcsinh.tsv")
write.table(export_df, file = out_tsv, sep = "\t",
            row.names = FALSE, quote = FALSE)
message(sprintf("  Salvato: 07_integrated_matrix_arcsinh.tsv"))
message(sprintf("  Dimensioni: %d cellule x %d colonne", nrow(export_df), ncol(export_df)))

message("=== STEP 7 completato ===\n")
