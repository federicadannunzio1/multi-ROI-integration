# =============================================================================
# 01_load_data.R
# Caricamento dati MACSima multi-campione
#
# Input:  cartelle paziente in data/ (definite in config.R)
# Output: output/data/01_combined_raw.rds   — data.frame cellule + metadati
#         output/data/protein_metadata.rds  — mappa nomi proteine
#
# Logica selezione file raw:
#   Per ogni ROI vengono cercati file CSV che:
#     1. Contengono colonne "Exp" (MFI non normalizzato)
#     2. NON contengono "Z-Score" (file normalizzato)
#     3. Contengono TUTTE le proteine attese (RAW_PROTEIN_NAMES)
#   Se piu' file soddisfano i criteri, viene preso il primo.
#   Se nessun file contiene tutte le proteine, il ROI viene saltato
#   con un warning (non interrompe la pipeline).
#
# Validazione:
#   - Cell ID univoci garantiti per costruzione
#   - Righe con valori non numerici segnalate e rimosse
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== STEP 1: Caricamento dati multi-campione ===")
message(sprintf("Campioni attesi: %d", length(SAMPLES)))

# --------------------------------------------------------------------------
# Funzione: trova il file raw in una cartella ROI
# Restituisce il path del file oppure NULL se non trovato (non fa stop).
# --------------------------------------------------------------------------

find_raw_file <- function(roi_dir, required_proteins) {
  if (!dir.exists(roi_dir)) stop("Cartella ROI non trovata: ", roi_dir)
  csvs <- list.files(roi_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) stop("Nessun CSV in: ", roi_dir)

  candidates <- character(0)
  for (f in csvs) {
    header     <- colnames(fread(f, nrows = 0))
    has_exp    <- any(grepl("Exp$", header))
    has_zscore <- any(grepl("Z-Score|Z Score|Z score", header, ignore.case = TRUE))
    has_all    <- all(required_proteins %in% header)
    if (has_exp && !has_zscore && has_all) candidates <- c(candidates, f)
  }

  if (length(candidates) == 0) return(NULL)
  if (length(candidates) > 1)
    message(sprintf("     [INFO] Trovati %d file candidati, uso il primo: %s",
                    length(candidates), basename(candidates[1])))
  return(candidates[1])
}

# --------------------------------------------------------------------------
# Caricamento
# --------------------------------------------------------------------------

df_list  <- list()
skipped  <- character(0)

for (samp in SAMPLES) {
  pid  <- samp$patient_id
  grp  <- samp$group
  rois <- samp$rois

  message(sprintf("\nCampione: %s (%s)", pid, grp))

  for (roi_name in names(rois)) {
    roi_dir  <- rois[[roi_name]]
    raw_file <- tryCatch(
      find_raw_file(roi_dir, RAW_PROTEIN_NAMES),
      error = function(e) stop(sprintf("[%s %s] %s", pid, roi_name, e$message))
    )

    if (is.null(raw_file)) {
      msg <- sprintf("[%s %s] Nessun file con tutte le proteine — ROI saltato", pid, roi_name)
      warning(msg)
      skipped <- c(skipped, paste0(pid, "_", roi_name))
      next
    }

    message(sprintf("  %s -> %s", roi_name, basename(raw_file)))
    dt <- fread(raw_file, header = TRUE, data.table = FALSE)
    message(sprintf("     %d cellule lette", nrow(dt)))

    # Seleziona e rinomina colonne proteiche
    dt <- dt[, RAW_PROTEIN_NAMES, drop = FALSE]
    colnames(dt) <- SHORT_NAMES

    # Converti a numerico (fread puo' leggere come character se ci sono valori anomali)
    dt <- as.data.frame(lapply(dt, as.numeric))

    # Rimuovi righe con NA (valori non numerici nei CSV originali)
    na_mask <- rowSums(is.na(dt)) > 0
    if (any(na_mask)) {
      message(sprintf("     [WARN] %d righe con NA rimosse", sum(na_mask)))
      dt <- dt[!na_mask, , drop = FALSE]
    }

    n_cells <- nrow(dt)
    if (n_cells == 0) {
      warning(sprintf("[%s %s] Nessuna cellula valida dopo QC — ROI saltato", pid, roi_name))
      next
    }

    # Cell ID univoci: patient_roi_cell_N
    cell_ids <- paste0(pid, "_", roi_name, "_cell_",
                       formatC(seq_len(n_cells), width = nchar(n_cells), flag = "0"))
    rownames(dt) <- cell_ids

    # Metadati
    dt$patient_id <- pid
    dt$group      <- grp
    dt$roi_id     <- paste0(pid, "_", roi_name)

    df_list[[paste0(pid, "_", roi_name)]] <- dt
    message(sprintf("     -> OK (%d cellule valide)", n_cells))
  }
}

# --------------------------------------------------------------------------
# Concatenazione
# --------------------------------------------------------------------------

if (length(skipped) > 0) {
  message(sprintf("\n[WARN] ROI saltati per proteine mancanti (%d): %s",
                  length(skipped), paste(skipped, collapse = ", ")))
}

if (length(df_list) == 0) stop("Nessun ROI valido caricato. Controlla i dati.")

df_combined <- do.call(rbind, df_list)
message(sprintf("\nTotale cellule combinate: %d", nrow(df_combined)))

# Riepilogo
message("\nRiepilogo per paziente:")
print(as.data.frame(df_combined %>%
  group_by(group, patient_id) %>%
  summarise(n_cells = n(), n_roi = n_distinct(roi_id), .groups = "drop")))

# Verifica unicita' cell ID
stopifnot("Cell ID duplicati" = !any(duplicated(rownames(df_combined))))
message("Cell ID: tutti univoci.")

# --------------------------------------------------------------------------
# Salvataggio
# --------------------------------------------------------------------------

saveRDS(df_combined, file.path(OUT_DATA, "01_combined_raw.rds"))
saveRDS(list(
  raw_names   = RAW_PROTEIN_NAMES,
  short_names = SHORT_NAMES,
  protein_map = PROTEIN_MAP
), file.path(OUT_DATA, "protein_metadata.rds"))

message("Salvato: 01_combined_raw.rds")
message("=== STEP 1 completato ===\n")
