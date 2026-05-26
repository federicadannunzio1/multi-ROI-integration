# =============================================================================
# 02_qc_normalization.R
# Quality Control e normalizzazione arcsinh dei dati MACSima
#
# Input:  01_combined_raw.rds
# Output: oggetto Seurat con dati raw e normalizzati
#         salvato come RDS in output/data/
#         plot QC in output/plots/
#
# Cosa fa questo script:
#   1. QC: identifica e rimuove oggetti non-cellulari (artefatti segmentazione)
#   2. Visualizza la distribuzione delle intensita' prima/dopo normalizzazione
#   3. Applica la trasformazione arcsinh (cofactor=1, standard per MICS)
#   4. Costruisce l'oggetto Seurat con i layer corretti:
#        - counts: valori raw MFI originali
#        - data:   valori arcsinh-trasformati (usati per PCA/UMAP/clustering)
#
# PERCHE' ARCSINH CON COFACTOR=1:
#   I valori raw sono Mean Fluorescence Intensity (MFI) da immagini a 16-bit,
#   con range tipico 0-17000+.
#   La trasformazione e': x_norm = arcsinh(x / cofactor)
#   - E' asintoticamente equivalente al log10 per valori grandi
#   - A differenza del log, gestisce correttamente i valori = 0 e negativi
#     (possibili dopo background subtraction)
#   - Cofactor=1 per immunofluorescenza: separa bene i picchi neg/pos
#     per valori nel range 0-20000 (Bendall et al., Science 2011;
#     Chevrier et al., Cell Systems 2018)
#   - Per confronto: cofactor=5 per CyTOF, cofactor=150 per flow citometria
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggridges)
library(viridis)

message("=== STEP 2: QC e normalizzazione ===")

BASE_OUT   <- "/Users/federicadannunzio/Documents/projects/collaborations/IgnazioCaruana/MACSima_pipeline/output"
PLOT_DIR   <- file.path(BASE_OUT, "plots")
DATA_DIR   <- file.path(BASE_OUT, "data")

# --------------------------------------------------------------------------
# Caricamento dati
# --------------------------------------------------------------------------

df_combined  <- readRDS(file.path(DATA_DIR, "01_combined_raw.rds"))
prot_meta    <- readRDS(file.path(DATA_DIR, "protein_metadata.rds"))
SHORT_NAMES  <- prot_meta$short_names
COFACTOR     <- 1   # Modifica qui se vuoi testare altri cofactor

message("Dati caricati: ", nrow(df_combined), " cellule x ",
        length(SHORT_NAMES), " proteine")

# --------------------------------------------------------------------------
# STEP 2a — QC: rimozione artefatti di segmentazione
# --------------------------------------------------------------------------
# In MACSima, la segmentazione puo' includere oggetti non-cellulari:
# background, detriti, bordi del tessuto. Questi tendono ad avere
# MFI=0 su tutti i canali (nessun segnale proteico rilevabile).

# Estrai solo la matrice proteica
mat_raw <- df_combined[, SHORT_NAMES]

# Flag celle con TUTTI i valori = 0
all_zero_mask <- rowSums(mat_raw == 0) == ncol(mat_raw)
n_allzero     <- sum(all_zero_mask)
message(sprintf("QC: cellule con tutti i valori = 0 (artefatti): %d (%.2f%%)",
                n_allzero, 100 * n_allzero / nrow(mat_raw)))

# Rimozione
df_clean <- df_combined[!all_zero_mask, ]
mat_raw  <- df_clean[, SHORT_NAMES]
message("Cellule rimaste dopo QC: ", nrow(df_clean))

# Riepilogo per ROI dopo QC
message("Cellule per ROI dopo QC:")
for (roi in unique(df_clean$roi_id)) {
  cat(sprintf("  %-6s: %7d\n", roi, sum(df_clean$roi_id == roi)))
}

# --------------------------------------------------------------------------
# STEP 2b — Plot distribuzione PRIMA della normalizzazione
# --------------------------------------------------------------------------
# Importante per valutare la qualita' del pannello e scegliere il cofactor.
# Ci aspettiamo distribuzioni asimmetriche con un picco a valori bassi
# (cellule negative) e una coda a valori alti (cellule positive).

message("Generazione plot QC pre-normalizzazione...")

# Converti in formato long per ggplot
df_long_raw <- tidyr::pivot_longer(
  cbind(cell_id = rownames(mat_raw), mat_raw),
  cols      = all_of(SHORT_NAMES),
  names_to  = "protein",
  values_to = "intensity"
)
df_long_raw$roi_id <- df_clean[df_long_raw$cell_id, "roi_id"]

# Density plot per-proteina con sovrapposizione per ROI
p_raw_density <- ggplot(df_long_raw,
                         aes(x = intensity, fill = roi_id, colour = roi_id)) +
  geom_density(alpha = 0.3, linewidth = 0.4) +
  facet_wrap(~protein, scales = "free") +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Distribuzione MFI raw per proteina e ROI",
       x = "MFI raw",
       y = "Densita'",
       fill = "ROI", colour = "ROI") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90"))

ggsave(file.path(PLOT_DIR, "QC_01_raw_density.pdf"),
       p_raw_density, width = 12, height = 8)
message("  Salvato: QC_01_raw_density.pdf")

# Boxplot per rilevare outlier globali
p_raw_box <- ggplot(df_long_raw, aes(x = roi_id, y = intensity, fill = roi_id)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.3) +
  facet_wrap(~protein, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Boxplot MFI raw per ROI",
       x = "ROI", y = "MFI raw", fill = "ROI") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PLOT_DIR, "QC_02_raw_boxplot.pdf"),
       p_raw_box, width = 12, height = 8)
message("  Salvato: QC_02_raw_boxplot.pdf")

# --------------------------------------------------------------------------
# STEP 2c — Trasformazione arcsinh
# --------------------------------------------------------------------------
# Formula: x_norm = arcsinh(x / cofactor)
# Proprieta':
#   - arcsinh(0) = 0 (le celle negative rimangono a 0)
#   - Per x >> cofactor: arcsinh(x/c) ≈ log(2x/c) (comportamento log)
#   - Simmetrica attorno allo 0 (gestisce valori negativi)
# Con cofactor=1 e valori ~0-17000: output range ~0-9.8

message(sprintf("Applicazione arcsinh (cofactor = %g)...", COFACTOR))
mat_arcsinh <- asinh(as.matrix(mat_raw) / COFACTOR)

message(sprintf("  Range post-trasformazione: [%.3f, %.3f]",
                min(mat_arcsinh), max(mat_arcsinh)))

# --------------------------------------------------------------------------
# STEP 2d — Plot distribuzione DOPO normalizzazione
# --------------------------------------------------------------------------
# Le distribuzioni devono diventare piu' simmetriche e comparabili tra ROI.
# Se le distribuzioni per lo stesso marker sono molto diverse tra ROI,
# potrebbe essere necessario il batch correction (Step 3).

df_long_norm <- as.data.frame(mat_arcsinh)
df_long_norm$cell_id <- rownames(mat_arcsinh)
df_long_norm$roi_id  <- df_clean[rownames(mat_arcsinh), "roi_id"]

df_long_norm <- tidyr::pivot_longer(df_long_norm,
  cols     = all_of(SHORT_NAMES),
  names_to = "protein", values_to = "arcsinh")

p_norm_density <- ggplot(df_long_norm,
                          aes(x = arcsinh, fill = roi_id, colour = roi_id)) +
  geom_density(alpha = 0.3, linewidth = 0.4) +
  facet_wrap(~protein, scales = "free") +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  labs(title = sprintf("Distribuzione arcsinh (cofactor=%g) per proteina e ROI", COFACTOR),
       x = "arcsinh(MFI / cofactor)",
       y = "Densita'", fill = "ROI", colour = "ROI") +
  theme_bw(base_size = 11)

ggsave(file.path(PLOT_DIR, "QC_03_arcsinh_density.pdf"),
       p_norm_density, width = 12, height = 8)
message("  Salvato: QC_03_arcsinh_density.pdf")

# Ridge plot: piu' compatto, buono per presentazioni
p_ridge <- ggplot(df_long_norm, aes(x = arcsinh, y = roi_id, fill = roi_id)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  facet_wrap(~protein, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Ridge plot arcsinh per ROI",
       x = "arcsinh(MFI / cofactor)", y = "ROI") +
  theme_ridges(font_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "QC_04_arcsinh_ridgeplot.pdf"),
       p_ridge, width = 12, height = 8)
message("  Salvato: QC_04_arcsinh_ridgeplot.pdf")

# --------------------------------------------------------------------------
# STEP 2e — Costruzione oggetto Seurat
# --------------------------------------------------------------------------
# Struttura dell'assay "MICS" nell'oggetto Seurat:
#   - counts: matrice MFI raw (features x cells)
#   - data:   matrice arcsinh-trasformata (features x cells)
#   - scale.data: viene riempito in Step 3 da ScaleData()
#
# NOTA: Seurat si aspetta features x cells (proteine x cellule = trasposto)

message("Costruzione oggetto Seurat...")

# Seurat vuole la matrice come features x cells
mat_raw_t     <- t(as.matrix(mat_raw))       # 6 x N_cells
mat_arcsinh_t <- t(mat_arcsinh)              # 6 x N_cells

# Crea l'oggetto con i raw counts
# min.cells=0 e min.features=0: no filtering automatico (gia' fatto in QC)
seurat_obj <- CreateSeuratObject(
  counts      = mat_raw_t,
  assay       = "MICS",
  min.cells   = 0,
  min.features = 0,
  meta.data   = data.frame(
    row.names = colnames(mat_raw_t),
    roi_id    = df_clean[colnames(mat_raw_t), "roi_id"]
  )
)

# Imposta i dati arcsinh come layer "data" (usato da PCA, UMAP, etc.)
seurat_obj <- SetAssayData(
  seurat_obj,
  assay    = "MICS",
  layer    = "data",      # In Seurat v4: slot="data"
  new.data = mat_arcsinh_t
)

message("Oggetto Seurat creato:")
print(seurat_obj)

# Verifica
n_per_roi <- table(seurat_obj$roi_id)
message("Cellule per ROI nell'oggetto Seurat:")
print(n_per_roi)

# --------------------------------------------------------------------------
# Salvataggio
# --------------------------------------------------------------------------

saveRDS(seurat_obj,
        file.path(DATA_DIR, "02_seurat_normalized.rds"))
message("Salvato: 02_seurat_normalized.rds")

message("=== STEP 2 completato ===\n")
