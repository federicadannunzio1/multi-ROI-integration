# =============================================================================
# 04_clustering.R
# Clustering Louvain e UMAP su tutte le cellule
#
# Input:  output/data/03_seurat_integrated.rds
# Output: output/data/04_seurat_clustered.rds
#         output/data/04_cell_cluster_assignments.csv
#         output/plots/Clustering_*.pdf
#
# NOTA CLUSTER vs LOCALE:
#   Su cluster con RAM sufficiente (>=128GB) si usa TUTTE le cellule
#   per FindNeighbors e FindClusters, senza subsampling.
#   Su Mac 8GB era necessario subsampliare — qui non piu'.
#
# SCELTA RESOLUTION:
#   Si eseguono piu' resolutions e si usa clustree per scegliere.
#   La resolution finale e' definita in config.R (FINAL_RESOLUTION).
#   Si consiglia di ispezionare Clustering_01_clustree.pdf prima di
#   procedere con gli step successivi e aggiornare FINAL_RESOLUTION.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(clustree)
  library(pheatmap)
  library(dplyr)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

# Con ~5M cellule, Seurat usa future per parallelizzare FindNeighbors/FindClusters.
# Il limite default (500MB) e' insufficiente — lo portiamo a 20GB.
options(future.globals.maxSize = 20 * 1024^3)

message("=== STEP 4: Clustering e UMAP ===")

seurat_obj <- readRDS(file.path(OUT_DATA, "03_seurat_integrated.rds"))
message("Oggetto caricato: ", ncol(seurat_obj), " cellule")

n_pcs_harmony <- ncol(Embeddings(seurat_obj, "harmony"))
message(sprintf("Dimensioni Harmony: %d", n_pcs_harmony))

# Palette paziente dinamica
all_patients <- sort(unique(seurat_obj$patient_id))
g3_patients  <- sort(unique(seurat_obj$patient_id[seurat_obj$group == "G3"]))
shh_patients <- sort(unique(seurat_obj$patient_id[seurat_obj$group == "SHH"]))
palette_patient <- c(
  setNames(colorRampPalette(c("#E41A1C", "#FC8D59"))(length(g3_patients)),  g3_patients),
  setNames(colorRampPalette(c("#377EB8", "#91BFDB"))(length(shh_patients)), shh_patients)
)

# --------------------------------------------------------------------------
# STEP 4a — FindNeighbors su tutte le cellule
# --------------------------------------------------------------------------

message(sprintf("FindNeighbors (k=%d, dims=1:%d) su %d cellule...",
                KNN_K, n_pcs_harmony, ncol(seurat_obj)))

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims      = seq_len(n_pcs_harmony),
  k.param   = KNN_K,
  verbose   = FALSE
)

# --------------------------------------------------------------------------
# STEP 4b — FindClusters: scansione multiple resolutions
# --------------------------------------------------------------------------

message(sprintf("FindClusters (resolutions: %s)...",
                paste(CLUSTERING_RESOLUTIONS, collapse = ", ")))

for (res in CLUSTERING_RESOLUTIONS) {
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution  = res,
    algorithm   = 1,  # Louvain
    random.seed = SEED,
    verbose     = FALSE
  )
  n_cl <- length(unique(seurat_obj$seurat_clusters))
  message(sprintf("  resolution=%.2f -> %d cluster", res, n_cl))
}

# --------------------------------------------------------------------------
# STEP 4c — Clustree: visualizzazione stabilita' dei cluster
# --------------------------------------------------------------------------

assay_prefix <- paste0(DefaultAssay(seurat_obj), "_snn_res.")

p_clustree <- clustree(seurat_obj, prefix = assay_prefix) +
  labs(title = "Clustree: stabilita' dei cluster per resolution",
       subtitle = "Scegliere la resolution dove i cluster si stabilizzano") +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_PLOTS, "Clustering_01_clustree.pdf"),
       p_clustree, width = 10, height = 12)
message("  Salvato: Clustering_01_clustree.pdf")
message(sprintf("  >> Ispeziona il clustree e verifica che FINAL_RESOLUTION=%.1f sia appropriata",
                FINAL_RESOLUTION))

# --------------------------------------------------------------------------
# STEP 4d — Assegna cluster dalla resolution finale
# --------------------------------------------------------------------------

final_cluster_col <- paste0(assay_prefix, FINAL_RESOLUTION)

if (!final_cluster_col %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Colonna '%s' non trovata nei metadati. ",
               "Verifica che FINAL_RESOLUTION (%.1f) sia in CLUSTERING_RESOLUTIONS.",
               final_cluster_col, FINAL_RESOLUTION))
}

seurat_obj$cluster <- seurat_obj[[final_cluster_col]][, 1]
Idents(seurat_obj) <- "cluster"
n_clusters <- length(unique(seurat_obj$cluster))
message(sprintf("Resolution finale: %.1f -> %d cluster (%d cellule totali)",
                FINAL_RESOLUTION, n_clusters, ncol(seurat_obj)))

palette_cluster <- colorRampPalette(brewer.pal(9, "Set3"))(n_clusters)

# --------------------------------------------------------------------------
# STEP 4e — UMAP su tutte le cellule
# --------------------------------------------------------------------------

message("UMAP in corso (su tutte le cellule — puo' richiedere tempo)...")
set.seed(SEED)
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction      = "harmony",
  dims           = seq_len(n_pcs_harmony),
  n.neighbors    = 30,
  min.dist       = 0.3,
  spread         = 1,
  reduction.name = "umap",
  verbose        = FALSE
)
message("UMAP completato.")

# --------------------------------------------------------------------------
# STEP 4f — Plot UMAP (subsample per visualizzazione — i punti non si vedono tutti)
# --------------------------------------------------------------------------

set.seed(SEED)
N_PLOT   <- min(100000, ncol(seurat_obj))
idx_plot <- sample(seq_len(ncol(seurat_obj)), N_PLOT)

umap_df            <- as.data.frame(Embeddings(seurat_obj, "umap")[idx_plot, ])
colnames(umap_df)  <- c("UMAP1", "UMAP2")
umap_df$cluster    <- factor(seurat_obj$cluster[idx_plot])
umap_df$group      <- seurat_obj$group[idx_plot]
umap_df$patient_id <- seurat_obj$patient_id[idx_plot]

# Centroidi cluster per label
cluster_centers <- umap_df %>%
  group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

p_umap_cluster <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
  geom_point(size = 0.05, alpha = 0.2) +
  geom_label(data = cluster_centers, aes(label = cluster),
             size = 3.5, fontface = "bold", colour = "black",
             fill = "white", alpha = 0.8, label.size = 0.2) +
  scale_colour_manual(values = palette_cluster) +
  labs(title = sprintf("UMAP — Cluster (res=%.1f)", FINAL_RESOLUTION),
       colour = "Cluster") +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

p_umap_group <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = group)) +
  geom_point(size = 0.05, alpha = 0.2) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(title = "UMAP — Gruppo (G3 vs SHH)", colour = "Gruppo") +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

p_umap_patient <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = patient_id)) +
  geom_point(size = 0.05, alpha = 0.2) +
  scale_colour_manual(values = palette_patient) +
  labs(title = "UMAP — Paziente", colour = "Paziente") +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(OUT_PLOTS, "Clustering_02_UMAP.pdf"),
       p_umap_cluster + p_umap_group + p_umap_patient,
       width = 21, height = 7)
message("  Salvato: Clustering_02_UMAP.pdf")

# Feature plots per proteina
seurat_sub_plot <- seurat_obj[, idx_plot]
p_features <- FeaturePlot(
  seurat_sub_plot,
  features   = rownames(seurat_obj),
  reduction  = "umap",
  pt.size    = 0.01,
  order      = FALSE,
  raster     = FALSE,
  ncol       = 3,
  cols       = c("lightgrey", "darkred")
) & theme_bw(base_size = 9) &
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(OUT_PLOTS, "Clustering_03_UMAP_features.pdf"),
       p_features, width = 12, height = 8)
message("  Salvato: Clustering_03_UMAP_features.pdf")

# --------------------------------------------------------------------------
# STEP 4g — Profili proteici per cluster (heatmap mediane)
# --------------------------------------------------------------------------

medians <- lapply(levels(seurat_obj$cluster), function(cl) {
  cells <- WhichCells(seurat_obj, idents = cl)
  data  <- GetAssayData(seurat_obj, layer = "data", assay = "MICS")[, cells]
  apply(data, 1, median)
})
median_mat <- do.call(cbind, medians)
colnames(median_mat) <- paste0("Cl_", levels(seurat_obj$cluster))

pheatmap(
  median_mat,
  scale                    = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize = 10,
  main = "Profili proteici mediani per cluster (Z-score per proteina)",
  filename = file.path(OUT_PLOTS, "Clustering_04_heatmap_cluster_profiles.pdf"),
  width = 8, height = 5
)
message("  Salvato: Clustering_04_heatmap_cluster_profiles.pdf")

# --------------------------------------------------------------------------
# Salvataggio
# --------------------------------------------------------------------------

saveRDS(seurat_obj, file.path(OUT_DATA, "04_seurat_clustered.rds"))
message("Salvato: 04_seurat_clustered.rds")

cluster_table <- data.frame(
  cell_id    = colnames(seurat_obj),
  patient_id = seurat_obj$patient_id,
  group      = seurat_obj$group,
  roi_id     = seurat_obj$roi_id,
  cluster    = seurat_obj$cluster
)
write.csv(cluster_table, file.path(OUT_DATA, "04_cell_cluster_assignments.csv"),
          row.names = FALSE, quote = FALSE)
message("Salvato: 04_cell_cluster_assignments.csv")
message("=== STEP 4 completato ===\n")
