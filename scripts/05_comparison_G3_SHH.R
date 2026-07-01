# =============================================================================
# 05_comparison_G3_SHH.R
# Confronto descrittivo G3 vs SHH a livello proteico
#
# Input:  output/data/04_seurat_clustered.rds
# Output: output/data/05_patient_medians_arcsinh.csv
#         output/plots/Comparison_*.pdf
#
# NOTA STATISTICA:
#   Con n paziente per gruppo variabile (non garantito n>=3), nessun test
#   statistico formale e' appropriato per confrontare i gruppi.
#   L'analisi e' DESCRITTIVA: mostriamo le distribuzioni single-cell
#   e le mediane per paziente (unita' biologica corretta).
#   Test su singola cellula (Mann-Whitney su milioni di cellule) sarebbero
#   pseudoreplicazione — p-value artificialmente bassi per campionamento.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== STEP 5: Confronto G3 vs SHH ===")

seurat_obj  <- readRDS(file.path(OUT_DATA, "04_seurat_clustered.rds"))
SHORT_NAMES <- rownames(seurat_obj)
COFACTOR    <- seurat_obj@misc$cofactor

# Palette paziente dinamica
all_patients <- sort(unique(seurat_obj$patient_id))
g3_patients  <- sort(unique(seurat_obj$patient_id[seurat_obj$group == "G3"]))
shh_patients <- sort(unique(seurat_obj$patient_id[seurat_obj$group == "SHH"]))
palette_patient <- c(
  setNames(PALETTE_PATIENTS_G3[seq_along(g3_patients)],   g3_patients),
  setNames(PALETTE_PATIENTS_SHH[seq_along(shh_patients)], shh_patients)
)

message(sprintf("Cellule: %d | Pazienti G3: %d | Pazienti SHH: %d",
                ncol(seurat_obj), length(g3_patients), length(shh_patients)))
print(table(seurat_obj$group, seurat_obj$patient_id))

mat_all <- GetAssayData(seurat_obj, layer = "data", assay = "MICS")

# --------------------------------------------------------------------------
# STEP 5a — Mediane per paziente
# --------------------------------------------------------------------------

message("Calcolo mediane per paziente...")

patient_medians <- lapply(unique(seurat_obj$patient_id), function(pid) {
  cells <- colnames(seurat_obj)[seurat_obj$patient_id == pid]
  grp   <- unique(seurat_obj$group[seurat_obj$patient_id == pid])
  meds  <- apply(mat_all[, cells], 1, median)
  data.frame(patient_id = pid, group = grp,
             protein = names(meds), median_arcsinh = round(meds, 4))
})
patient_medians_df <- do.call(rbind, patient_medians)

write.csv(patient_medians_df,
          file.path(OUT_DATA, "05_patient_medians_arcsinh.csv"),
          row.names = FALSE, quote = FALSE)
message("  Salvato: 05_patient_medians_arcsinh.csv")

patient_med_wide <- pivot_wider(
  patient_medians_df[, c("patient_id","group","protein","median_arcsinh")],
  names_from  = "protein",
  values_from = "median_arcsinh"
)
message("\nMediane per paziente (arcsinh):")
print(as.data.frame(patient_med_wide))

# --------------------------------------------------------------------------
# STEP 5b — Heatmap mediane per paziente
# --------------------------------------------------------------------------

mat_heatmap <- as.matrix(patient_med_wide[, SHORT_NAMES])
rownames(mat_heatmap) <- paste0(patient_med_wide$patient_id,
                                " (", patient_med_wide$group, ")")

ann_row <- data.frame(Group = patient_med_wide$group,
                      row.names = rownames(mat_heatmap))

pheatmap(
  mat_heatmap,
  scale                    = "column",
  annotation_row           = ann_row,
  annotation_colors        = list(Group = PALETTE_GROUP),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize = 11,
  main = sprintf("Patient arcsinh medians (Z-score per protein, cofactor = %d)", COFACTOR),
  filename = file.path(OUT_PLOTS, "Comparison_01_heatmap_patient_medians.pdf"),
  width = 9, height = 7
)
message("  Salvato: Comparison_01_heatmap_patient_medians.pdf")

# --------------------------------------------------------------------------
# STEP 5c — Violin plot per gruppo (subsample per visualizzazione)
# --------------------------------------------------------------------------

set.seed(SEED)
cells_sub <- unlist(lapply(unique(seurat_obj$group), function(grp) {
  cells <- colnames(seurat_obj)[seurat_obj$group == grp]
  sample(cells, min(100000, length(cells)))
}))

df_violin <- cbind(
  data.frame(group = seurat_obj$group[cells_sub],
             patient_id = seurat_obj$patient_id[cells_sub]),
  t(as.matrix(mat_all[, cells_sub]))
)

df_violin_long <- pivot_longer(df_violin,
  cols = all_of(SHORT_NAMES), names_to = "protein", values_to = "arcsinh"
)

p_violin <- ggplot(df_violin_long,
                   aes(x = group, y = arcsinh, fill = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0, fill = "white", alpha = 0.8) +
  geom_point(
    data = patient_medians_df,
    aes(x = group, y = median_arcsinh, shape = patient_id),
    colour = "black", size = 2.5, inherit.aes = FALSE
  ) +
  facet_wrap(~protein, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title    = "Protein expression: G3 vs SHH",
    subtitle = sprintf("Violin = single-cell distribution | Dots = patient median | arcsinh(MFI / %d)", COFACTOR),
    x = "Group", y = sprintf("arcsinh(MFI / %d)", COFACTOR),
    fill = "Group", shape = "Patient"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_PLOTS, "Comparison_02_violin_G3_SHH.pdf"),
       p_violin, width = 12, height = 10)
message("  Salvato: Comparison_02_violin_G3_SHH.pdf")

# --------------------------------------------------------------------------
# STEP 5d — Barplot mediane per paziente
# --------------------------------------------------------------------------
# palette_patient e' gia' definita sopra come gradiente per gruppo.
# Usiamo direttamente i colori per paziente nel fill.

p_bar <- ggplot(patient_medians_df,
                aes(x = protein, y = median_arcsinh,
                    fill = patient_id)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7 / length(all_patients)) +
  scale_fill_manual(values = palette_patient) +
  facet_wrap(~group) +
  labs(
    title    = "Patient median expression: G3 vs SHH",
    subtitle = sprintf("arcsinh(MFI / %d)", COFACTOR),
    x = "Protein", y = "Median arcsinh", fill = "Patient"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_PLOTS, "Comparison_03_barplot_patient_medians.pdf"),
       p_bar, width = 11, height = 6)
message("  Salvato: Comparison_03_barplot_patient_medians.pdf")

message("=== STEP 5 completato ===\n")
