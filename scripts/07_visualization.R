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
  setNames(colorRampPalette(c("#E41A1C", "#FC8D59"))(length(g3_patients)),  g3_patients),
  setNames(colorRampPalette(c("#377EB8", "#91BFDB"))(length(shh_patients)), shh_patients)
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
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Gruppo") +
  theme_classic(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

p_patient <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = patient_id)) +
  geom_point(size = 0.1, alpha = 0.2) +
  scale_colour_manual(values = palette_patient) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Paziente") +
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
  scale_colour_gradient2(low = "steelblue", mid = "white", high = "darkred",
                         midpoint = 0) +
  labs(title = "Dot Plot per cluster", x = "Proteina", y = "Cluster") +
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
  Gruppo  = seurat_obj$group[cells_heatmap],
  row.names = cells_heatmap
)
ann_colors <- list(
  Cluster = setNames(palette_cluster, levels(seurat_obj$cluster)),
  Gruppo  = PALETTE_GROUP
)

pheatmap(mat_sample,
  annotation_col  = ann_col,
  annotation_colors = ann_colors,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_colnames   = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  scale = "row", fontsize_row = 10,
  main = paste0("Heatmap single-cell (", n_per_cluster, " per cluster)"),
  filename = file.path(OUT_PLOTS, "Final_04_heatmap_singlecell.pdf"),
  width = 12, height = 5
)
message("  Salvato: Final_04_heatmap_singlecell.pdf")

# --------------------------------------------------------------------------
# STEP 7d — Export matrice integrata completa
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
