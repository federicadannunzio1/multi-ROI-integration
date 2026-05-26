# =============================================================================
# 05_visualization.R
# Visualizzazioni finali e interpretazione dei cluster
#
# Input:  04_seurat_clustered.rds
# Output: plot diagnostici e di presentazione in output/plots/
#         tabella marker per cluster in output/data/
#
# Cosa fa questo script:
#   1. Dot plot: espressione media + % cellule positive per proteina/cluster
#   2. Heatmap a singola cellula (subsample) per validazione visiva
#   3. UMAP con label grandi (version publication-ready)
#   4. Tabella marker: mediana, media, % positive per ogni cluster
#   5. Annotazione automatica suggerita in base ai profili (da validare)
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(scales)

message("=== STEP 5: Visualizzazioni finali ===")

BASE_OUT <- "/Users/federicadannunzio/Documents/projects/collaborations/IgnazioCaruana/MACSima_pipeline/output"
PLOT_DIR <- file.path(BASE_OUT, "plots")
DATA_DIR <- file.path(BASE_OUT, "data")

# --------------------------------------------------------------------------
# Caricamento
# --------------------------------------------------------------------------

seurat_obj  <- readRDS(file.path(DATA_DIR, "04_seurat_clustered.rds"))
SHORT_NAMES <- rownames(seurat_obj)
n_clusters  <- length(unique(seurat_obj$cluster))

message("Oggetto caricato: ", ncol(seurat_obj), " cellule, ",
        n_clusters, " cluster")

palette_cluster <- colorRampPalette(brewer.pal(9, "Set3"))(n_clusters)
palette_roi     <- brewer.pal(4, "Set1")

# --------------------------------------------------------------------------
# STEP 5a — Dot Plot
# --------------------------------------------------------------------------
# Combina due informazioni in un unico grafico:
#   - Dimensione del punto: % di cellule nel cluster con espressione rilevabile
#   - Colore del punto:     espressione media (arcsinh)
# E' il plot piu' informativo per pannelli piccoli.

p_dot <- DotPlot(
  seurat_obj,
  features      = SHORT_NAMES,
  group.by      = "cluster",
  assay         = "MICS",
  scale         = TRUE,
  col.min       = -2.5,
  col.max       = 2.5,
  dot.scale     = 6
) +
  scale_colour_gradient2(
    low  = "steelblue",
    mid  = "white",
    high = "darkred",
    midpoint = 0
  ) +
  labs(title = "Dot Plot: profili proteici per cluster",
       subtitle = "Dimensione = % cellule positive | Colore = espressione media scalata",
       x = "Proteina", y = "Cluster") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PLOT_DIR, "Final_01_dotplot.pdf"),
       p_dot, width = 9, height = max(5, n_clusters * 0.6 + 2))
message("  Salvato: Final_01_dotplot.pdf")

# --------------------------------------------------------------------------
# STEP 5b — Heatmap a singola cellula (subsample)
# --------------------------------------------------------------------------
# Mostra i valori reali delle singole cellule (non le mediane).
# Con ~461K cellule non e' possibile visualizzarle tutte:
# facciamo subsample di 200 cellule per cluster.

set.seed(42)
n_per_cluster <- 200
cells_sample  <- lapply(levels(seurat_obj$cluster), function(cl) {
  all_cells <- WhichCells(seurat_obj, idents = cl)
  sample(all_cells, min(n_per_cluster, length(all_cells)))
}) |> unlist()

mat_sample <- GetAssayData(seurat_obj, layer = "data", assay = "MICS")[, cells_sample]

# Annotazione colonne (cluster e ROI)
ann_col <- data.frame(
  Cluster = seurat_obj$cluster[cells_sample],
  ROI     = seurat_obj$roi_id[cells_sample],
  row.names = cells_sample
)

ann_colors <- list(
  Cluster = setNames(palette_cluster, levels(seurat_obj$cluster)),
  ROI     = setNames(palette_roi, c("ROI1", "ROI2", "ROI3", "ROI4"))
)

pheatmap(
  mat_sample,
  annotation_col       = ann_col,
  annotation_colors    = ann_colors,
  cluster_rows         = TRUE,
  cluster_cols         = TRUE,
  show_colnames        = FALSE,
  color                = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  scale                = "row",
  fontsize_row         = 10,
  main                 = paste0("Heatmap single-cell (subsample ", n_per_cluster, " per cluster)"),
  filename             = file.path(PLOT_DIR, "Final_02_heatmap_singlecell.pdf"),
  width = 12, height = 5
)
message("  Salvato: Final_02_heatmap_singlecell.pdf")

# --------------------------------------------------------------------------
# STEP 5c — UMAP publication-ready
# --------------------------------------------------------------------------

p_pub <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "cluster",
  pt.size = 0.003, alpha = 0.2, label = TRUE, label.size = 5,
  repel = TRUE
) +
  scale_color_manual(values = palette_cluster) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2", colour = "Cluster") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.line = element_line(linewidth = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggsave(file.path(PLOT_DIR, "Final_03_UMAP_publication.pdf"),
       p_pub, width = 7, height = 6)
ggsave(file.path(PLOT_DIR, "Final_03_UMAP_publication.png"),
       p_pub, width = 7, height = 6, dpi = 300)
message("  Salvato: Final_03_UMAP_publication.pdf/png")

# UMAP split per ROI (vedi eterogeneita' spaziale)
p_split <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "cluster",
  split.by = "roi_id",
  pt.size  = 0.003, alpha = 0.2, label = FALSE
) +
  scale_color_manual(values = palette_cluster) +
  labs(title = "UMAP split per ROI") +
  theme_bw(base_size = 10)

ggsave(file.path(PLOT_DIR, "Final_04_UMAP_split_roi.pdf"),
       p_split, width = 14, height = 4)
message("  Salvato: Final_04_UMAP_split_roi.pdf")

# --------------------------------------------------------------------------
# STEP 5d — Tabella marker per cluster
# --------------------------------------------------------------------------
# Per ogni cluster: media, mediana, % cellule positive (arcsinh > 0.5)
# per ogni proteina. Questa tabella e' utile per l'annotazione manuale.

POSITIVITY_THRESHOLD <- 0.5  # arcsinh > 0.5 = cellula "positiva"

mat_all <- GetAssayData(seurat_obj, layer = "data", assay = "MICS")

marker_stats <- lapply(levels(seurat_obj$cluster), function(cl) {
  cells   <- WhichCells(seurat_obj, idents = cl)
  submat  <- mat_all[, cells, drop = FALSE]
  n_cells <- ncol(submat)

  stats <- lapply(SHORT_NAMES, function(prot) {
    vals    <- submat[prot, ]
    data.frame(
      cluster    = cl,
      n_cells    = n_cells,
      protein    = prot,
      mean       = round(mean(vals), 3),
      median     = round(median(vals), 3),
      pct_pos    = round(100 * mean(vals > POSITIVITY_THRESHOLD), 1)
    )
  })
  do.call(rbind, stats)
})
marker_stats <- do.call(rbind, marker_stats)

write.csv(marker_stats,
          file.path(DATA_DIR, "05_cluster_marker_stats.csv"),
          row.names = FALSE, quote = FALSE)
message("  Salvato: 05_cluster_marker_stats.csv")

# Stampa a schermo riepilogo
message("\nProfilo cluster (mediana arcsinh):")
pivot_table <- tidyr::pivot_wider(
  marker_stats[, c("cluster", "protein", "median")],
  names_from  = "protein",
  values_from = "median"
)
print(pivot_table)

# --------------------------------------------------------------------------
# STEP 5e — Suggerimento annotazione biologica
# --------------------------------------------------------------------------
# Con questo pannello di 6 marker, i tipi cellulari attesi sono:
#   Ki67+     : cellule proliferanti (tumorali o stromali)
#   LC3B+ P62+: cellule con autofagia attiva (o blocco autofagico se P62 alto)
#   LC3B+ P62-: autofagia funzionale (flux attivo)
#   Nectin2+  : cellule che esprimono CD112 (ligando checkpoint TIGIT)
#   PVR+      : cellule che esprimono CD155 (ligando checkpoint TIGIT/CD226)
#   pMTOR+    : mTOR pathway attivo (inibisce autofagia)
#
# L'interpretazione finale richiede validazione con il biologo/clinico.
# Usa la tabella 05_cluster_marker_stats.csv come punto di partenza.

message("\n--- NOTE PER ANNOTAZIONE BIOLOGICA ---")
message("Marker nel pannello:")
message("  Ki67    -> proliferazione")
message("  LC3B    -> autofagosome (marker autofagia)")
message("  P62     -> substrato autofagia (alto = blocco autofagico)")
message("  pMTOR   -> mTOR attivo (soppressore autofagia)")
message("  Nectin2 -> CD112 (ligando TIGIT, immune checkpoint)")
message("  PVR     -> CD155 (ligando TIGIT/CD226, immune checkpoint)")
message("  LC3B+/P62- = autofagia attiva (flux)")
message("  LC3B+/P62+ = blocco autofagico")
message("  pMTOR+/LC3B- = inibizione autofagia via mTOR")
message("---------------------------------------")

# --------------------------------------------------------------------------
# STEP 5f — Export matrice integrata completa
# --------------------------------------------------------------------------
# Esporta un'unica tabella TSV con:
#   - cell_id       : identificatore univoco della cellula (ROI1_cell_0, ...)
#   - roi_id        : ROI di provenienza
#   - cluster       : cluster assegnato
#   - UMAP_1, UMAP_2: coordinate UMAP
#   - Ki67 ... PVR  : valori arcsinh-trasformati (layer "data")
#
# Formato TSV (tab-separated):
#   - Leggibile da R (read.table), Python (pd.read_csv sep='\t'), Excel
#   - Preferito a xlsx per dataset grandi: nessun limite di righe,
#     file piu' leggeri, nessuna dipendenza da pacchetti extra
#
# NOTA: i valori esportati sono arcsinh(MFI / cofactor=1), NON i raw MFI.
# Se vuoi anche i raw, cambia layer = "counts" e salva un secondo file.

message("Export matrice integrata...")

# Matrice arcsinh (proteine x cellule) -> trasponi a cellule x proteine
mat_export <- t(as.matrix(
  GetAssayData(seurat_obj, layer = "data", assay = "MICS")
))

# Coordinate UMAP
umap_coords <- Embeddings(seurat_obj, "umap")

# Assembla dataframe finale
export_df <- data.frame(
  cell_id = rownames(mat_export),
  roi_id  = seurat_obj$roi_id,
  cluster = seurat_obj$cluster,
  UMAP_1  = round(umap_coords[, 1], 4),
  UMAP_2  = round(umap_coords[, 2], 4),
  round(mat_export, 4),
  check.names = FALSE
)

out_tsv <- file.path(DATA_DIR, "06_integrated_matrix_arcsinh.tsv")
write.table(export_df,
            file      = out_tsv,
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)

message(sprintf("  Salvato: 06_integrated_matrix_arcsinh.tsv"))
message(sprintf("  Dimensioni: %d cellule x %d colonne", nrow(export_df), ncol(export_df)))
message(sprintf("  Colonne: cell_id, roi_id, cluster, UMAP_1, UMAP_2, %s",
                paste(SHORT_NAMES, collapse = ", ")))

# Matrice minima: solo cellule x proteine (arcsinh), senza metadati
out_mat <- file.path(DATA_DIR, "06_protein_matrix_arcsinh.tsv")
write.table(round(mat_export, 4),
            file      = out_mat,
            sep       = "\t",
            row.names = TRUE,   # cell_id come rownames
            col.names = TRUE,
            quote     = FALSE)
message("  Salvato: 06_protein_matrix_arcsinh.tsv (cellule x proteine, solo valori)")

message("=== STEP 5 completato ===\n")
message("Pipeline completata. Output in:")
message("  Plot:  ", PLOT_DIR)
message("  Dati:  ", DATA_DIR)
