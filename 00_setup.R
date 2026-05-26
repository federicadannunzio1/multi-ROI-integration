# =============================================================================
# 00_setup.R
# Verifica e installazione dei pacchetti necessari alla pipeline
#
# Eseguire SOLO la prima volta oppure dopo un aggiornamento di R.
# Tutti gli altri script assumono che i pacchetti siano gia' installati.
# =============================================================================

message("=== Verifica dipendenze ===")

# --------------------------------------------------------------------------
# Pacchetti CRAN
# --------------------------------------------------------------------------
cran_packages <- c(
  "Seurat",       # Framework principale per single-cell analysis
  "harmony",      # Batch correction (Korsunsky et al., Nature Methods 2019)
  "ggplot2",      # Visualizzazione base
  "dplyr",        # Manipolazione dataframe
  "tidyr",        # Reshape dataframe (pivot_longer/wider)
  "data.table",   # Lettura rapida CSV grandi (fread >> read.csv per file >50MB)
  "patchwork",    # Composizione multipanel di ggplot
  "RColorBrewer", # Palette colori
  "viridis",      # Scale colori per heatmap
  "scales",       # Utility per assi e scale
  "ggridges",     # Density ridge plots (QC distribuzioni)
  "pheatmap",     # Heatmap gerarchica per profili cluster
  "clustree"      # Visualizzazione stabilita' clustering a diverse resolution
)

# Installa solo i pacchetti mancanti
missing_cran <- cran_packages[!cran_packages %in% installed.packages()[, "Package"]]

if (length(missing_cran) > 0) {
  message("Installazione pacchetti CRAN mancanti: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
} else {
  message("Tutti i pacchetti CRAN sono gia' installati.")
}

# --------------------------------------------------------------------------
# Verifica versioni critiche
# --------------------------------------------------------------------------
seurat_ver <- packageVersion("Seurat")
message("Seurat version: ", seurat_ver)

if (as.numeric(seurat_ver$major) < 5) {
  warning(
    "Rilevato Seurat v", seurat_ver, ". ",
    "Questa pipeline e' ottimizzata per Seurat v5. ",
    "Alcune chiamate (SetAssayData layer=) potrebbero richiedere la sostituzione ",
    "di 'layer' con 'slot'. Aggiornamento consigliato: install.packages('Seurat')"
  )
}

message("Setup completato.")
