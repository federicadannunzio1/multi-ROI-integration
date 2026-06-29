# =============================================================================
# 03_integration.R
# Integrazione multi-campione con Harmony (batch correction per paziente)
#
# Input:  output/data/02_seurat_normalized.rds
# Output: output/data/03_seurat_integrated.rds
#         output/plots/Integration_*.pdf
#
# STRATEGIA:
#   Harmony corregge la variabilita' tecnica inter-paziente (batch),
#   preservando la variabilita' biologica G3 vs SHH.
#   group.by.vars = "patient_id": ogni paziente e' trattato come un batch.
#   NON correggiamo per "group" (G3/SHH): e' la variabile di interesse.
#
# CON 6 PROTEINE:
#   La PCA produce al massimo 5 componenti significative (p-1).
#   Harmony opera su questo spazio 5D.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== STEP 3: Integrazione con Harmony ===")

seurat_obj <- readRDS(file.path(OUT_DATA, "02_seurat_normalized.rds"))
message("Oggetto caricato: ", ncol(seurat_obj), " cellule")
message("Pazienti: ", paste(sort(unique(seurat_obj$patient_id)), collapse = ", "))

# Palette paziente dinamica (G3 = toni rossi, SHH = toni blu)
all_patients  <- sort(unique(seurat_obj$patient_id))
g3_patients   <- sort(all_patients[all_patients %in%
                    colnames(seurat_obj)[seurat_obj$group == "G3"] |
                    all_patients %in% unique(seurat_obj$patient_id[seurat_obj$group == "G3"])])
shh_patients  <- sort(setdiff(all_patients, g3_patients))
palette_patient <- c(
  setNames(colorRampPalette(c("#E41A1C", "#FC8D59"))(length(g3_patients)),  g3_patients),
  setNames(colorRampPalette(c("#377EB8", "#91BFDB"))(length(shh_patients)), shh_patients)
)

# --------------------------------------------------------------------------
# STEP 3a — ScaleData e PCA
# --------------------------------------------------------------------------

message("ScaleData...")
seurat_obj <- ScaleData(seurat_obj, assay = "MICS",
                        features = rownames(seurat_obj), verbose = FALSE)

n_pcs <- length(rownames(seurat_obj)) - 1  # max PCA dims = n_features - 1
message(sprintf("PCA (%d componenti)...", n_pcs))
seurat_obj <- RunPCA(seurat_obj, assay = "MICS",
                     features = rownames(seurat_obj),
                     npcs = n_pcs, reduction.name = "pca", verbose = FALSE)

var_explained <- (seurat_obj[["pca"]]@stdev)^2 /
  sum((seurat_obj[["pca"]]@stdev)^2) * 100
message("Varianza spiegata per PC:")
for (i in seq_along(var_explained))
  cat(sprintf("  PC%d: %.1f%%\n", i, var_explained[i]))

# --------------------------------------------------------------------------
# STEP 3b — Plot PCA pre-Harmony (subsample per visualizzazione)
# --------------------------------------------------------------------------

set.seed(SEED)
N_PLOT   <- min(50000, ncol(seurat_obj))
idx_plot <- sample(seq_len(ncol(seurat_obj)), N_PLOT)

pca_df           <- as.data.frame(Embeddings(seurat_obj, "pca")[idx_plot, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
pca_df$patient_id <- seurat_obj$patient_id[idx_plot]
pca_df$group      <- seurat_obj$group[idx_plot]

p_pre_pat <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = patient_id)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_colour_manual(values = palette_patient) +
  labs(title = "PCA pre-Harmony", subtitle = "Colorato per paziente",
       colour = "Paziente") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw()

p_pre_grp <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(title = "PCA pre-Harmony", subtitle = "Colorato per gruppo",
       colour = "Gruppo") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw()

ggsave(file.path(OUT_PLOTS, "Integration_01_PCA_pre_harmony.pdf"),
       p_pre_pat + p_pre_grp, width = 14, height = 6)
message("  Salvato: Integration_01_PCA_pre_harmony.pdf")

# --------------------------------------------------------------------------
# STEP 3c — Harmony
# --------------------------------------------------------------------------

message(sprintf("Harmony (batch = patient_id, theta=%.1f, max_iter=%d)...",
                HARMONY_THETA, HARMONY_MAX_ITER))

set.seed(SEED)
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars    = "patient_id",
  reduction.use    = "pca",
  reduction.save   = "harmony",
  dims.use         = seq_len(n_pcs),
  theta            = HARMONY_THETA,
  max_iter         = HARMONY_MAX_ITER,
  plot_convergence = FALSE,
  verbose          = FALSE
)
message("Harmony completato.")

# --------------------------------------------------------------------------
# STEP 3d — Plot post-Harmony
# --------------------------------------------------------------------------

harm_df           <- as.data.frame(Embeddings(seurat_obj, "harmony")[idx_plot, 1:2])
colnames(harm_df) <- c("H1", "H2")
harm_df$patient_id <- seurat_obj$patient_id[idx_plot]
harm_df$group      <- seurat_obj$group[idx_plot]

p_post_pat <- ggplot(harm_df, aes(x = H1, y = H2, colour = patient_id)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_colour_manual(values = palette_patient) +
  labs(title = "Post-Harmony", subtitle = "I pazienti devono sovrapporsi",
       x = "Harmony 1", y = "Harmony 2", colour = "Paziente") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw()

p_post_grp <- ggplot(harm_df, aes(x = H1, y = H2, colour = group)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(title = "Post-Harmony", subtitle = "G3 e SHH devono restare separati",
       x = "Harmony 1", y = "Harmony 2", colour = "Gruppo") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw()

ggsave(file.path(OUT_PLOTS, "Integration_02_Harmony_post.pdf"),
       p_post_pat + p_post_grp, width = 14, height = 6)
message("  Salvato: Integration_02_Harmony_post.pdf")

p_comparison <- (p_pre_pat + p_pre_grp) / (p_post_pat + p_post_grp) +
  plot_annotation(title = "Batch correction: PCA vs Harmony",
                  subtitle = "Sopra: pre-Harmony | Sotto: post-Harmony")

ggsave(file.path(OUT_PLOTS, "Integration_03_comparison_pre_post.pdf"),
       p_comparison, width = 14, height = 12)
message("  Salvato: Integration_03_comparison_pre_post.pdf")

# --------------------------------------------------------------------------
# Salvataggio
# --------------------------------------------------------------------------

saveRDS(seurat_obj, file.path(OUT_DATA, "03_seurat_integrated.rds"))
message("Salvato: 03_seurat_integrated.rds")
message("=== STEP 3 completato ===\n")
