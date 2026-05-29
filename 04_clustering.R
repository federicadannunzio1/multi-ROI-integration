# =============================================================================
# 04_clustering.R
# Clustering cellulare e UMAP
#
# Input:  03_seurat_integrated.rds
# Output: oggetto Seurat con cluster e UMAP
#         salvato come RDS in output/data/
#         plot in output/plots/
#
# Cosa fa questo script:
#   1. FindNeighbors: costruisce il grafo k-NN sull'embedding Harmony
#   2. FindClusters: algoritmo di Louvain per community detection
#      con scansione di piu' resolution per scegliere la granularita'
#   3. RunUMAP: proiezione 2D per visualizzazione
#   4. clustree: visualizzazione della stabilita' dei cluster
#      al variare della resolution
#   5. Heatmap dei profili proteici per cluster
#
# SCELTA DELL'ALGORITMO DI CLUSTERING:
#   Louvain (default Seurat) vs Leiden:
#   - Leiden (Traag et al., Scientific Reports 2019) produce cluster piu'
#     stabili e risolve i problemi di disconnessione del Louvain.
#     E' raccomandato nella letteratura recente (Luecken & Theis 2019).
#   - In Seurat, Leiden richiede il pacchetto 'leiden' e igraph>=1.2.
#   - Per semplicita', usiamo Louvain (algorithm=1) come default e
#     offriamo Leiden (algorithm=4) come alternativa commentata.
#
# RESOLUTION:
#   Con solo 6 proteine, pochi cluster biologicamente distinti sono attesi
#   (probabilmente 3-8 popolazioni cellulari basate su questi marker).
#   Facciamo una scansione da 0.1 a 1.5 e usiamo clustree per scegliere.
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(clustree)
library(pheatmap)
library(dplyr)

message("=== STEP 4: Clustering e UMAP ===")

BASE_OUT <- "/Users/federicadannunzio/Documents/projects/collaborations/IgnazioCaruana/MACSima_pipeline/output"
PLOT_DIR <- file.path(BASE_OUT, "plots")
DATA_DIR <- file.path(BASE_OUT, "data")

# --------------------------------------------------------------------------
# Caricamento
# --------------------------------------------------------------------------

seurat_obj <- readRDS(file.path(DATA_DIR, "03_seurat_integrated.rds"))
message("Oggetto Seurat caricato: ", ncol(seurat_obj), " cellule")

SHORT_NAMES <- rownames(seurat_obj)

# --------------------------------------------------------------------------
# STEP 4a — FindNeighbors
# --------------------------------------------------------------------------
# Costruisce il grafo k-NN (k-nearest neighbors) nello spazio Harmony.
# k=30 e' un buon compromesso per dataset grandi (>50K cellule):
#   - k piccolo -> cluster piccoli e molto granulari
#   - k grande  -> cluster piu' robusti ma potenzialmente sovrasmoothati
# Con 461K cellule e 6 proteine, k=30 e' appropriato.

n_pcs_harmony <- ncol(Embeddings(seurat_obj, "harmony"))
message(sprintf("FindNeighbors (k=30, dims=1:%d da Harmony)...", n_pcs_harmony))
seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims      = seq_len(n_pcs_harmony),
  k.param   = 30,
  verbose   = FALSE
)

# --------------------------------------------------------------------------
# STEP 4b — Scansione resolution con FindClusters
# --------------------------------------------------------------------------
# Eseguiamo FindClusters a piu' resolution e salviamo ogni risultato.
# Poi usiamo clustree per scegliere la resolution piu' stabile.
# La resolution giusta e' quella dove l'aggiunta di una nuova resolution
# non sposta molte cellule tra cluster esistenti.

resolutions <- c(0.1, 0.2, 0.3, 0.5)

message("Scansione resolution: ", paste(resolutions, collapse = ", "))
for (res in resolutions) {
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution     = res,
    algorithm      = 1,   # Louvain; usa 4 per Leiden (richiede pkg 'leiden')
    random.seed    = 42,
    verbose        = FALSE
  )
  message(sprintf("  resolution=%.1f -> %d cluster",
                  res, length(unique(seurat_obj$seurat_clusters))))
}

# --------------------------------------------------------------------------
# STEP 4c — Clustree: visualizzazione stabilita' dei cluster
# --------------------------------------------------------------------------
# clustree mostra come le cellule si spostano tra cluster al variare
# della resolution. Scegli la resolution dove le frecce (migrazioni
# di cellule) si stabilizzano.

p_clustree <- clustree(seurat_obj, prefix = "MICS_snn_res.")
ggsave(file.path(PLOT_DIR, "Clustering_01_clustree.pdf"),
       p_clustree, width = 12, height = 10)
message("  Salvato: Clustering_01_clustree.pdf")
message("AZIONE RICHIESTA: Guarda Clustering_01_clustree.pdf e scegli la resolution")
message("  (suggerito: resolution con minime migrazioni inter-cluster)")

# --------------------------------------------------------------------------
# STEP 4d — Scelta resolution finale e impostazione cluster
# --------------------------------------------------------------------------
# MODIFICA QUESTA VARIABILE in base all'ispezione del clustree plot.
# Default: 0.5 (buon punto di partenza per pannelli piccoli)

FINAL_RESOLUTION <- 0.1


# Imposta i cluster finali come colonna 'leiden_clusters' (nome generico
# anche se usiamo Louvain, per compatibilita' futura con Leiden)
final_cluster_col <- paste0("MICS_snn_res.", FINAL_RESOLUTION)
seurat_obj$cluster <- seurat_obj[[final_cluster_col]][, 1]
Idents(seurat_obj) <- "cluster"

n_clusters <- length(unique(seurat_obj$cluster))
message(sprintf("Resolution finale: %.1f -> %d cluster", FINAL_RESOLUTION, n_clusters))

# --------------------------------------------------------------------------
# STEP 4e — UMAP
# --------------------------------------------------------------------------
# L'UMAP viene calcolato sull'embedding Harmony (non sulla PCA raw).
# min.dist e spread controllano la "compattezza" della proiezione:
#   - min.dist=0.3, spread=1 (default): cluster compatti
#   - min.dist=0.5, spread=1.5: piu' spazio tra cluster, utile con molte cellule

message("UMAP in corso...")
set.seed(42)
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
# STEP 4f — Plot UMAP base
# --------------------------------------------------------------------------

palette_roi     <- brewer.pal(4, "Set1")
palette_cluster <- colorRampPalette(brewer.pal(9, "Set3"))(n_clusters)

# UMAP colorato per cluster
p_umap_cluster <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "cluster",
  pt.size = 0.005, alpha = 0.3, label = TRUE, label.size = 4
) +
  scale_color_manual(values = palette_cluster) +
  labs(title = sprintf("UMAP — Cluster (resolution=%.1f)", FINAL_RESOLUTION)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

# UMAP colorato per ROI (verifica che i ROI siano ben miscelati)
p_umap_roi <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "roi_id",
  pt.size = 0.005, alpha = 0.3
) +
  scale_color_manual(values = palette_roi) +
  labs(title = "UMAP — ROI (verifica mixing dopo Harmony)") +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

p_umap_combined <- p_umap_cluster + p_umap_roi
ggsave(file.path(PLOT_DIR, "Clustering_02_UMAP_cluster_roi.pdf"),
       p_umap_combined, width = 14, height = 6)
message("  Salvato: Clustering_02_UMAP_cluster_roi.pdf")

# UMAP espresso con intensita' per ogni proteina (FeaturePlot)
# Usa i valori arcsinh (layer "data")
p_features <- FeaturePlot(
  seurat_obj,
  features  = SHORT_NAMES,
  reduction = "umap",
  pt.size   = 0.005,
  order     = FALSE,
  raster    = FALSE,
  ncol      = 3,
  cols      = c("lightgrey", "darkred")
) &
  theme_bw(base_size = 9) &
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(PLOT_DIR, "Clustering_03_UMAP_feature_plots.pdf"),
       p_features, width = 12, height = 8)
message("  Salvato: Clustering_03_UMAP_feature_plots.pdf")

# --------------------------------------------------------------------------
# STEP 4g — Profili proteici per cluster (heatmap + violin)
# --------------------------------------------------------------------------

# Calcola mediana arcsinh per proteina per cluster
medians <- lapply(levels(seurat_obj$cluster), function(cl) {
  cells <- WhichCells(seurat_obj, idents = cl)
  data  <- GetAssayData(seurat_obj, layer = "data", assay = "MICS")[, cells]
  apply(data, 1, median)
})
median_mat <- do.call(cbind, medians)
colnames(median_mat) <- paste0("Cl_", levels(seurat_obj$cluster))

# Heatmap gerarchica dei profili cluster
p_heatmap <- pheatmap(
  median_mat,
  scale        = "row",     # Z-score per riga (proteina) per comparare profili
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  color        = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize     = 10,
  main         = "Profili proteici mediani per cluster (Z-score per proteina)",
  filename     = file.path(PLOT_DIR, "Clustering_04_heatmap_cluster_profiles.pdf"),
  width        = 8, height = 5
)
message("  Salvato: Clustering_04_heatmap_cluster_profiles.pdf")

# Violin plot per ogni proteina e cluster
p_violin <- VlnPlot(
  seurat_obj,
  features  = SHORT_NAMES,
  group.by  = "cluster",
  pt.size   = 0,         # Niente punti (troppe cellule)
  ncol      = 3,
  assay     = "MICS",
  layer     = "data"
) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(PLOT_DIR, "Clustering_05_violin_per_cluster.pdf"),
       p_violin, width = 12, height = 8)
message("  Salvato: Clustering_05_violin_per_cluster.pdf")

# --------------------------------------------------------------------------
# STEP 4h — Composizione per ROI (proporzioni cluster per ROI)
# --------------------------------------------------------------------------
# Questo plot mostra quante cellule di ogni cluster ci sono per ROI.
# Utile per valutare se la composizione cellulare e' diversa tra i pezzi
# di tessuto (eterogeneita' intra-paziente).

comp_df <- seurat_obj@meta.data %>%
  group_by(roi_id, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(roi_id) %>%
  mutate(freq = n / sum(n))

p_composition <- ggplot(comp_df, aes(x = roi_id, y = freq, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = palette_cluster) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Composizione cluster per ROI",
       subtitle = "Eterogeneita' spaziale intra-paziente",
       x = "ROI", y = "Proporzione cellule", fill = "Cluster") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PLOT_DIR, "Clustering_06_composition_per_ROI.pdf"),
       p_composition, width = 8, height = 6)
message("  Salvato: Clustering_06_composition_per_ROI.pdf")

# --------------------------------------------------------------------------
# Salvataggio
# --------------------------------------------------------------------------

saveRDS(seurat_obj, file.path(DATA_DIR, "04_seurat_clustered.rds"))
message("Salvato: 04_seurat_clustered.rds")

# Esporta anche tabella cell ID -> cluster per uso esterno
cluster_table <- data.frame(
  cell_id  = colnames(seurat_obj),
  roi_id   = seurat_obj$roi_id,
  cluster  = seurat_obj$cluster
)
write.csv(cluster_table,
          file.path(DATA_DIR, "04_cell_cluster_assignments.csv"),
          row.names = FALSE, quote = FALSE)
message("Salvato: 04_cell_cluster_assignments.csv")

message("=== STEP 4 completato ===\n")
