# =============================================================================
# run_pipeline.R
# Script master: esegue l'intera pipeline in sequenza
#
# UTILIZZO:
#   Da terminale:
#     Rscript run_pipeline.R
#   Da RStudio:
#     source("run_pipeline.R")
#
# STRUTTURA OUTPUT:
#   output/
#   ├── data/
#   │   ├── protein_metadata.rds         <- mappa nomi proteine
#   │   ├── 01_combined_raw.rds          <- dati unificati pre-QC
#   │   ├── 02_seurat_normalized.rds     <- Seurat object post-arcsinh
#   │   ├── 03_seurat_integrated.rds     <- Seurat object post-Harmony
#   │   ├── 04_seurat_clustered.rds      <- Seurat object finale con cluster
#   │   ├── 04_cell_cluster_assignments.csv  <- cell_id -> cluster
#   │   └── 05_cluster_marker_stats.csv  <- profili proteici per cluster
#   └── plots/
#       ├── QC_01_raw_density.pdf
#       ├── QC_02_raw_boxplot.pdf
#       ├── QC_03_arcsinh_density.pdf
#       ├── QC_04_arcsinh_ridgeplot.pdf
#       ├── Integration_01_PCA_pre_harmony.pdf
#       ├── Integration_02_Harmony_post_correction.pdf
#       ├── Integration_03_comparison_pre_post.pdf
#       ├── Clustering_01_clustree.pdf       <- GUARDA QUI per scegliere resolution
#       ├── Clustering_02_UMAP_cluster_roi.pdf
#       ├── Clustering_03_UMAP_feature_plots.pdf
#       ├── Clustering_04_heatmap_cluster_profiles.pdf
#       ├── Clustering_05_violin_per_cluster.pdf
#       ├── Clustering_06_composition_per_ROI.pdf
#       ├── Final_01_dotplot.pdf
#       ├── Final_02_heatmap_singlecell.pdf
#       ├── Final_03_UMAP_publication.pdf/png
#       └── Final_04_UMAP_split_roi.pdf
#
# PARAMETRI MODIFICABILI:
#   In 02_qc_normalization.R:
#     COFACTOR (default=1)   -> cofactor arcsinh
#   In 04_clustering.R:
#     FINAL_RESOLUTION       -> resolution Leiden dopo aver visto clustree
#     resolutions            -> range di resolution da esplorare
# =============================================================================

# Imposta working directory alla posizione dello script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Se eseguito da terminale (Rscript), usa invece:
# setwd(getSrcDirectory(function(){}))

PIPELINE_DIR <- file.path(
  "/Users/federicadannunzio/Documents/projects/collaborations/IgnazioCaruana",
  "MACSima_pipeline"
)

SCRIPT_DIR <- file.path(PIPELINE_DIR, "scripts")

# --------------------------------------------------------------------------
# Crea struttura output se non esiste
# --------------------------------------------------------------------------

dirs_to_create <- c(
  file.path(PIPELINE_DIR, "output"),
  file.path(PIPELINE_DIR, "output", "data"),
  file.path(PIPELINE_DIR, "output", "plots")
)
for (d in dirs_to_create) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --------------------------------------------------------------------------
# Esecuzione pipeline
# --------------------------------------------------------------------------

steps <- list(
  list(script = "01_load_data.R",       label = "Step 1: Caricamento dati"),
  list(script = "02_qc_normalization.R", label = "Step 2: QC e normalizzazione"),
  list(script = "03_integration.R",      label = "Step 3: Integrazione Harmony"),
  list(script = "04_clustering.R",       label = "Step 4: Clustering e UMAP"),
  list(script = "05_visualization.R",    label = "Step 5: Visualizzazioni finali")
)

# Opzionale: esegui solo dal passo N (utile se vuoi ripartire da un punto)
# Imposta START_FROM = 1 per eseguire tutto dall'inizio
START_FROM <- 1

t_total_start <- proc.time()

for (i in seq_along(steps)) {
  if (i < START_FROM) {
    message("--- Saltato: ", steps[[i]]$label, " ---")
    next
  }

  message("\n", strrep("=", 60))
  message(steps[[i]]$label)
  message(strrep("=", 60))

  t_step_start <- proc.time()
  source(file.path(SCRIPT_DIR, steps[[i]]$script), local = FALSE)

  elapsed <- (proc.time() - t_step_start)["elapsed"]
  message(sprintf("[%s completato in %.1f minuti]",
                  steps[[i]]$label, elapsed / 60))
}

t_total <- (proc.time() - t_total_start)["elapsed"]
message("\n", strrep("=", 60))
message(sprintf("PIPELINE COMPLETATA in %.1f minuti totali", t_total / 60))
message(strrep("=", 60))
message("\nOutput salvato in:")
message("  ", file.path(PIPELINE_DIR, "output"))
