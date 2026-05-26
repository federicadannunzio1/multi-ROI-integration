# =============================================================================
# 01_load_data.R
# Caricamento e pre-processing iniziale dei dati MACSima (WO Normalization)
#
# Input:  4 file CSV raw (WO Normalization), uno per ROI
# Output: un singolo data.frame unificato con cell ID e roi_id
#         salvato come RDS in output/data/
#
# Cosa fa questo script:
#   1. Carica i file raw (NON quelli con Z-score — vedi motivazione sotto)
#   2. Standardizza i nomi e l'ordine delle colonne (ROI2 e ROI3 hanno ordine
#      diverso nei file raw — bug silenzioso se non corretto)
#   3. Assegna ID cellulari unici del tipo ROI1_cell_0, ROI1_cell_1, ...
#   4. Aggiunge la variabile roi_id come metadato
#   5. Concatena tutto in un unico oggetto
#
# PERCHE' I FILE RAW E NON GLI Z-SCORE:
#   I file Z-score di ROI1 e ROI2 hanno un numero diverso di righe rispetto
#   ai file raw (ROI1: 93489 vs 174528). Questo indica che il Z-score e' stato
#   esportato da una sottopopolazione filtrata in MACS iQ View (probabilmente
#   un singolo cluster). Usare i file Z-score escluderebbe quindi una parte
#   delle cellule. Inoltre, lo Z-score e' calcolato per-ROI separatamente:
#   Z=1 in ROI1 non e' comparabile con Z=1 in ROI2. Per l'integrazione
#   multi-campione, partiamo sempre dai valori raw e normalizziamo noi
#   in modo coerente su tutti i ROI insieme.
# =============================================================================

library(data.table)
library(dplyr)

message("=== STEP 1: Caricamento dati ===")

# --------------------------------------------------------------------------
# Percorsi
# --------------------------------------------------------------------------
BASE_RAW  <- "/Users/federicadannunzio/Documents/projects/collaborations/IgnazioCaruana/Gr3_24-0268"
BASE_OUT  <- "/Users/federicadannunzio/Documents/projects/collaborations/IgnazioCaruana/MACSima_pipeline/output/data"

raw_files <- list(
  ROI1 = file.path(BASE_RAW, "ROI1", "2 Heatmap_18_data_WO Normalization.csv"),
  ROI2 = file.path(BASE_RAW, "ROI2", "2 Heat Map in Segmentation 0 WO Normalization.csv"),
  ROI3 = file.path(BASE_RAW, "ROI3", "2 Heat Map in Segmentation 0 WO Normalization.csv"),
  ROI4 = file.path(BASE_RAW, "ROI4", "1 Heat Map in Segmentation 0 WO Normalization.csv")
)

# Verifica esistenza file
for (roi in names(raw_files)) {
  if (!file.exists(raw_files[[roi]])) {
    stop("File non trovato per ", roi, ": ", raw_files[[roi]])
  }
}
message("Tutti i file trovati.")

# --------------------------------------------------------------------------
# Nomi canonici delle proteine
# Definiti qui in modo centralizzato: questo e' l'unico posto
# dove mappare i nomi raw -> nomi corti per le visualizzazioni.
# --------------------------------------------------------------------------

# Nomi originali cosi' come compaiono nei file raw
RAW_PROTEIN_NAMES <- c(
  "Ki_67 REA1123 Biomarker Exp",
  "LC3B Cyto Exp",
  "Nectin 2 Biomarker Exp",
  "P62 Cyto Exp",
  "Phospo mTor 49F9 Cyto Exp",
  "Poliovirus Receptor Biomarker Exp"
)

# Nomi corti per visualizzazioni (stesso ordine di RAW_PROTEIN_NAMES)
SHORT_NAMES <- c("Ki67", "LC3B", "Nectin2", "P62", "pMTOR", "PVR")

# Mappa: nome raw -> nome corto
protein_map <- setNames(SHORT_NAMES, RAW_PROTEIN_NAMES)

# --------------------------------------------------------------------------
# Caricamento e standardizzazione
# --------------------------------------------------------------------------

df_list <- lapply(names(raw_files), function(roi_id) {

  message("  Caricamento ", roi_id, " ...")
  path <- raw_files[[roi_id]]

  # fread e' molto piu' veloce di read.csv per file grandi (>50MB)
  dt <- fread(path, header = TRUE, data.table = FALSE)

  n_cells   <- nrow(dt)
  n_cols    <- ncol(dt)
  message("    -> ", n_cells, " cellule, ", n_cols, " proteine")

  # Verifica che tutte le proteine attese siano presenti
  missing_proteins <- setdiff(RAW_PROTEIN_NAMES, colnames(dt))
  if (length(missing_proteins) > 0) {
    stop("Proteine mancanti in ", roi_id, ": ", paste(missing_proteins, collapse = ", "))
  }

  # Riordina le colonne in ordine canonico
  # IMPORTANTE: ROI2 e ROI3 hanno un ordine diverso nel file raw.
  # Se non si corregge, la concatenazione assegna valori sbagliati alle proteine.
  dt <- dt[, RAW_PROTEIN_NAMES, drop = FALSE]

  # Rinomina con nomi corti
  colnames(dt) <- SHORT_NAMES

  # Assegna cell ID univoci: ROIx_cell_0, ROIx_cell_1, ...
  # Il formato zero-padded e' utile per ordinamento lessicografico coerente
  n_digits  <- nchar(as.character(n_cells))
  cell_ids  <- paste0(roi_id, "_cell_", formatC(seq(0, n_cells - 1),
                                                 width = n_digits,
                                                 flag  = "0"))
  rownames(dt) <- cell_ids

  # Aggiungi metadato ROI
  dt$roi_id <- roi_id

  dt
})

names(df_list) <- names(raw_files)

# --------------------------------------------------------------------------
# Riepilogo pre-concatenazione
# --------------------------------------------------------------------------

message("\nRiepilogo per ROI:")
for (roi_id in names(df_list)) {
  cat(sprintf("  %-6s: %7d cellule\n", roi_id, nrow(df_list[[roi_id]])))
}

# --------------------------------------------------------------------------
# Concatenazione
# --------------------------------------------------------------------------

df_combined <- do.call(rbind, df_list)

message(sprintf("\nTotale cellule combinate: %d", nrow(df_combined)))
message(sprintf("Proteine: %s", paste(SHORT_NAMES, collapse = ", ")))

# Verifica unicita' cell ID
if (any(duplicated(rownames(df_combined)))) {
  stop("ERRORE: cell ID duplicati dopo concatenazione. Verificare la logica di assegnazione.")
}
message("Verifica cell ID: tutti univoci.")

# --------------------------------------------------------------------------
# Salvataggio
# --------------------------------------------------------------------------

out_path <- file.path(BASE_OUT, "01_combined_raw.rds")
saveRDS(df_combined, out_path)
message("Salvato: ", out_path)

# Salva anche la mappa proteine per uso negli script successivi
saveRDS(list(
  raw_names   = RAW_PROTEIN_NAMES,
  short_names = SHORT_NAMES,
  protein_map = protein_map
), file.path(BASE_OUT, "protein_metadata.rds"))

message("=== STEP 1 completato ===\n")
