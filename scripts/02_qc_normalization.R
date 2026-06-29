# =============================================================================
# 02_qc_normalization.R
# QC, selezione empirica del cofactor e normalizzazione arcsinh
#
# Input:  output/data/01_combined_raw.rds
# Output: output/data/02_seurat_normalized.rds
#         output/data/02_cofactor_quantiles.xlsx
#         output/plots/QC_*.pdf
#
# COFACTOR:
#   Scelto empiricamente come mediana dei p50 (mediana) di tutti i valori
#   MFI raw, per ogni combinazione proteina x paziente.
#   Questo garantisce che la transizione lineare->logaritmica avvenga
#   nell'intervallo biologicamente rilevante per l'intero dataset,
#   non ottimizzato per un singolo campione.
#
# TRASFORMAZIONE ARCSINH:
#   arcsinh(x / cofactor)
#   - Lineare per x << cofactor (range basso, spesso rumore)
#   - Logaritmica per x >> cofactor (range alto, segnale biologico)
#   - Gestisce valori = 0 e negativi (al contrario di log)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggridges)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== STEP 2: QC, selezione cofactor e normalizzazione ===")

df_combined <- readRDS(file.path(OUT_DATA, "01_combined_raw.rds"))
message(sprintf("Dati caricati: %d cellule x %d proteine",
                nrow(df_combined), length(SHORT_NAMES)))

# --------------------------------------------------------------------------
# STEP 2a — QC: rimozione artefatti
# --------------------------------------------------------------------------

mat_raw       <- df_combined[, SHORT_NAMES]
all_zero_mask <- rowSums(mat_raw == 0, na.rm = TRUE) == ncol(mat_raw)
message(sprintf("QC: cellule con tutti i valori = 0: %d (%.2f%%)",
                sum(all_zero_mask), 100 * mean(all_zero_mask)))

df_clean <- df_combined[!all_zero_mask, ]
mat_raw  <- df_clean[, SHORT_NAMES]

message(sprintf("Cellule rimaste dopo QC: %d", nrow(df_clean)))

message("\nCellule per paziente dopo QC:")
print(as.data.frame(df_clean %>%
  group_by(group, patient_id) %>%
  summarise(n_cells = n(), n_roi = n_distinct(roi_id), .groups = "drop")))

# --------------------------------------------------------------------------
# STEP 2b — Selezione cofactor empirico
# --------------------------------------------------------------------------

message("\nAnalisi distribuzione MFI raw per selezione cofactor...")

probs      <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
prob_names <- c("p1","p5","p10","p25","p50","p75","p90","p95","p99")

quantile_rows <- list()
for (pid in unique(df_clean$patient_id)) {
  grp  <- unique(df_clean$group[df_clean$patient_id == pid])
  mask <- df_clean$patient_id == pid
  for (prot in SHORT_NAMES) {
    vals <- as.numeric(mat_raw[mask, prot])
    vals <- vals[!is.na(vals)]
    q    <- quantile(vals, probs = probs, na.rm = TRUE)
    pcts <- sapply(c(50, 100, 200, 500, 1000),
                   function(t) round(100 * mean(vals < t, na.rm = TRUE), 1))
    quantile_rows[[length(quantile_rows) + 1]] <- c(
      list(patient_id = pid, group = grp, protein = prot,
           min = round(min(vals)), max = round(max(vals))),
      setNames(as.list(round(q)), prob_names),
      list(pct_below_50   = pcts[1], pct_below_100  = pcts[2],
           pct_below_200  = pcts[3], pct_below_500  = pcts[4],
           pct_below_1000 = pcts[5])
    )
  }
}
quantile_df <- do.call(rbind, lapply(quantile_rows, as.data.frame))

write.csv(quantile_df, file.path(OUT_DATA, "02_cofactor_quantiles.csv"),
          row.names = FALSE, quote = FALSE)
message("  Salvato: 02_cofactor_quantiles.csv")

COFACTOR <- round(median(quantile_df$p50))
message(sprintf("\nCofactor scelto: %d", COFACTOR))
message(sprintf("  (mediana dei p50 su %d combinazioni proteina x paziente)",
                nrow(quantile_df)))

p50_table <- pivot_wider(
  quantile_df[, c("patient_id","group","protein","p50")],
  names_from = "protein", values_from = "p50"
)
message("  p50 per proteina x paziente:")
print(as.data.frame(p50_table))

# Plot cofactor selection
p_cofactor <- ggplot(quantile_df,
                     aes(x = protein, y = p50, colour = group, shape = patient_id)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = COFACTOR, linetype = "dashed") +
  annotate("text", x = 1, y = COFACTOR * 1.05,
           label = paste0("cofactor = ", COFACTOR), hjust = 0, size = 3.5) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(title = "Distribuzione p50 MFI raw per proteina e paziente",
       subtitle = "Linea tratteggiata = cofactor scelto (mediana dei p50)",
       x = "Proteina", y = "p50 MFI raw",
       colour = "Gruppo", shape = "Paziente") +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_PLOTS, "QC_00_cofactor_selection.pdf"),
       p_cofactor, width = 10, height = 6)
message("  Salvato: QC_00_cofactor_selection.pdf")

# --------------------------------------------------------------------------
# STEP 2c — Plot distribuzione raw per paziente
# --------------------------------------------------------------------------

df_long_raw <- pivot_longer(
  cbind(df_clean[, c("patient_id", "group")], mat_raw),
  cols = all_of(SHORT_NAMES), names_to = "protein", values_to = "intensity"
)

p_raw <- ggplot(df_long_raw,
                aes(x = intensity, fill = patient_id, colour = patient_id)) +
  geom_density(alpha = 0.3, linewidth = 0.4) +
  facet_grid(protein ~ group, scales = "free") +
  labs(title = "Distribuzione MFI raw per proteina e gruppo",
       x = "MFI raw", y = "Densita'",
       fill = "Paziente", colour = "Paziente") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey90"))

ggsave(file.path(OUT_PLOTS, "QC_01_raw_density_all_samples.pdf"),
       p_raw, width = 14, height = 16)
message("  Salvato: QC_01_raw_density_all_samples.pdf")

# --------------------------------------------------------------------------
# STEP 2d — Trasformazione arcsinh
# --------------------------------------------------------------------------

message(sprintf("\nApplicazione arcsinh (cofactor = %d)...", COFACTOR))
mat_arcsinh <- asinh(as.matrix(mat_raw) / COFACTOR)
message(sprintf("  Range post-trasformazione: [%.3f, %.3f]",
                min(mat_arcsinh), max(mat_arcsinh)))

# --------------------------------------------------------------------------
# STEP 2e — Plot distribuzione arcsinh
# --------------------------------------------------------------------------

df_long_norm <- cbind(
  df_clean[, c("patient_id", "group")],
  as.data.frame(mat_arcsinh)
)
df_long_norm <- pivot_longer(df_long_norm,
  cols = all_of(SHORT_NAMES), names_to = "protein", values_to = "arcsinh"
)

p_norm <- ggplot(df_long_norm,
                 aes(x = arcsinh, fill = patient_id, colour = patient_id)) +
  geom_density(alpha = 0.3, linewidth = 0.4) +
  facet_grid(protein ~ group, scales = "free_y") +
  labs(title = sprintf("Distribuzione arcsinh (cofactor=%d)", COFACTOR),
       x = sprintf("arcsinh(MFI / %d)", COFACTOR), y = "Densita'",
       fill = "Paziente", colour = "Paziente") +
  theme_bw(base_size = 10)

ggsave(file.path(OUT_PLOTS, "QC_02_arcsinh_density_all_samples.pdf"),
       p_norm, width = 14, height = 16)
message("  Salvato: QC_02_arcsinh_density_all_samples.pdf")

p_ridge <- ggplot(df_long_norm,
                  aes(x = arcsinh, y = patient_id, fill = group)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  facet_wrap(~protein, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(title = "Ridge plot arcsinh per paziente e gruppo",
       x = sprintf("arcsinh(MFI / %d)", COFACTOR), y = "Paziente",
       fill = "Gruppo") +
  theme_ridges(font_size = 9) +
  theme(legend.position = "top")

ggsave(file.path(OUT_PLOTS, "QC_03_arcsinh_ridgeplot.pdf"),
       p_ridge, width = 12, height = 12)
message("  Salvato: QC_03_arcsinh_ridgeplot.pdf")

# --------------------------------------------------------------------------
# STEP 2f — Costruzione oggetto Seurat
# --------------------------------------------------------------------------

message("\nCostruzione oggetto Seurat...")

mat_raw_t     <- t(as.matrix(mat_raw))
mat_arcsinh_t <- t(mat_arcsinh)

seurat_obj <- CreateSeuratObject(
  counts       = mat_raw_t,
  assay        = "MICS",
  min.cells    = 0,
  min.features = 0,
  meta.data    = data.frame(
    row.names  = colnames(mat_raw_t),
    patient_id = df_clean[colnames(mat_raw_t), "patient_id"],
    group      = df_clean[colnames(mat_raw_t), "group"],
    roi_id     = df_clean[colnames(mat_raw_t), "roi_id"]
  )
)

seurat_obj <- SetAssayData(seurat_obj, assay = "MICS", layer = "data",
                           new.data = mat_arcsinh_t)

seurat_obj@misc$cofactor <- COFACTOR

message("Oggetto Seurat creato:")
print(seurat_obj)

message("\nCellule per paziente:")
print(table(seurat_obj$group, seurat_obj$patient_id))

# --------------------------------------------------------------------------
# Salvataggio
# --------------------------------------------------------------------------

saveRDS(seurat_obj, file.path(OUT_DATA, "02_seurat_normalized.rds"))
message("\nSalvato: 02_seurat_normalized.rds")
message("=== STEP 2 completato ===\n")
