# =============================================================================
# 03_integration.R
# Integrazione multi-ROI con Harmony (batch correction)
#
# Input:  02_seurat_normalized.rds
# Output: oggetto Seurat con embedding PCA corretto per batch
#         salvato come RDS in output/data/
#         plot diagnostici in output/plots/
#
# Cosa fa questo script:
#   1. ScaleData: scala i valori arcsinh per marcatore (media=0, var=1)
#      necessario per evitare che proteine con range diversi dominino la PCA
#   2. PCA su tutti e 6 i marcatori
#   3. Harmony: corregge l'embedding PCA per batch effect inter-ROI,
#      producendo X_pca_harmony usato per clustering e UMAP
#   4. Plot diagnostici: confronto PCA prima/dopo Harmony
#
# PERCHE' HARMONY E NON ALTRI METODI:
#   - Harmony (Korsunsky et al., Nature Methods 2019) e' il metodo con
#     le migliori performance in benchmark per integrazione single-cell
#     (Luecken et al., Nature Methods 2022)
#   - Lavora sull'embedding PCA (non sui valori raw), quindi e' rapido
#     anche con milioni di cellule
#   - Conserva la variabilita' biologica separando il segnale tecnico
#     di batch dal segnale biologico
#   - Produce direttamente l'embedding corretto da usare per clustering/UMAP
#
# NOTA SU 6 PROTEINE:
#   Con solo 6 features, la PCA produce al massimo 6 componenti.
#   In questo caso usiamo tutte e 6 (nessuna selezione di PC significative
#   con elbow plot come si farebbe con dati ad alta dimensione).
#   Il segnale biologico sara' comunque catturato in questi 6 assi.
#
# BATCH EFFECT IN QUESTO STUDIO:
#   Anche se i 4 ROI vengono dallo stesso paziente e (presumibilmente)
#   dallo stesso run strumentale, la variabilita' tecnica inter-ROI
#   puo' esistere per:
#   - Differenze nello spessore della sezione
#   - Variabilita' locale nell'efficienza di staining/bleaching
#   - Posizione nel campo visivo durante l'acquisizione
#   Harmony corregge queste differenze senza oversmoothare la biologia.
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

message("=== STEP 3: Integrazione con Harmony ===")

BASE_OUT <- "/Users/federicadannunzio/Documents/projects/collaborations/IgnazioCaruana/MACSima_pipeline/output"
PLOT_DIR <- file.path(BASE_OUT, "plots")
DATA_DIR <- file.path(BASE_OUT, "data")

# --------------------------------------------------------------------------
# Caricamento
# --------------------------------------------------------------------------

seurat_obj <- readRDS(file.path(DATA_DIR, "02_seurat_normalized.rds"))
message("Oggetto Seurat caricato: ", ncol(seurat_obj), " cellule")

# --------------------------------------------------------------------------
# STEP 3a — ScaleData
# --------------------------------------------------------------------------
# Scala ogni proteina a media=0 e deviazione standard=1.
# Questo e' importante prima della PCA perche':
#   - Elimina l'effetto di scale diverse tra proteine
#     (es. Ki67 con range 0-8 vs pMTOR con range 0-10 dopo arcsinh)
#   - Assicura che ogni proteina contribuisca equamente alla PCA
# I valori scalati vengono salvati nel layer "scale.data" dell'assay.

message("ScaleData in corso...")
seurat_obj <- ScaleData(seurat_obj,
                         assay    = "MICS",
                         features = rownames(seurat_obj))

# --------------------------------------------------------------------------
# STEP 3b — PCA
# --------------------------------------------------------------------------
# Con 6 proteine, n_comps massimo = 6.
# Usiamo tutte e 6 per non perdere informazione.

message("PCA in corso (6 componenti)...")
seurat_obj <- RunPCA(
  seurat_obj,
  assay        = "MICS",
  features     = rownames(seurat_obj),
  npcs         = 6,
  reduction.name = "pca",
  verbose      = FALSE
)

# Mostra varianza spiegata
var_explained <- (seurat_obj[["pca"]]@stdev)^2 /
  sum((seurat_obj[["pca"]]@stdev)^2) * 100
message("Varianza spiegata per PC:")
for (i in seq_along(var_explained)) {
  cat(sprintf("  PC%d: %.1f%%\n", i, var_explained[i]))
}

# --------------------------------------------------------------------------
# STEP 3c — Plot PCA pre-Harmony
# --------------------------------------------------------------------------
# Plot diagnostici con ggplot2 diretto (no DimPlot).
# DimPlot puo' produrre plot vuoti con assay custom in Seurat v5.
# Estraiamo le coordinate con Embeddings() e plottiamo manualmente.
# Subsample di 50K cellule per leggibilita' del PDF.

set.seed(42)
N_PLOT     <- min(50000, ncol(seurat_obj))
idx_plot   <- sample(seq_len(ncol(seurat_obj)), N_PLOT)
message(sprintf("  Subsample per plot diagnostici: %d cellule su %d",
                N_PLOT, ncol(seurat_obj)))

# Dataframe per i plot: coordinate PCA + roi_id
pca_coords <- as.data.frame(Embeddings(seurat_obj, "pca")[idx_plot, 1:2])
colnames(pca_coords) <- c("PC1", "PC2")
pca_coords$roi_id <- seurat_obj$roi_id[idx_plot]

p_pca_roi <- ggplot(pca_coords, aes(x = PC1, y = PC2, colour = roi_id)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_colour_brewer(palette = "Set1") +
  labs(title    = "PCA pre-Harmony (colorato per ROI)",
       subtitle = sprintf("Subsample %dk celle — se i ROI si separano -> batch effect rilevante",
                          N_PLOT / 1000),
       colour = "ROI") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw()

ggsave(file.path(PLOT_DIR, "Integration_01_PCA_pre_harmony.pdf"),
       p_pca_roi, width = 8, height = 6)
message("  Salvato: Integration_01_PCA_pre_harmony.pdf")

# --------------------------------------------------------------------------
# STEP 3d — Harmony batch correction
# --------------------------------------------------------------------------
# Harmony prende l'embedding PCA e restituisce un embedding corretto
# (X_pca_harmony) in cui la variabilita' dovuta al batch (roi_id) e'
# minimizzata, preservando la struttura biologica.
#
# Parametri chiave:
#   - vars_use:       variabile di batch (roi_id)
#   - theta:          forza della correzione (default=2, aumentare se i batch
#                     sono ancora visibili, diminuire se si perde troppa biologia)
#   - max_iter_harmony: iterazioni massime (20 e' sufficiente per convergenza)
#   - plot_convergence: True per vedere la convergenza dell'algoritmo

message("Harmony in corso...")

# Rileva dinamicamente quante PC sono state effettivamente calcolate.
# Con pannelli piccoli (6 proteine) RunPCA puo' produrre meno di npcs
# componenti se alcune features hanno varianza quasi zero dopo ScaleData.
n_pcs_computed <- ncol(Embeddings(seurat_obj, "pca"))
message(sprintf("  PC calcolate da RunPCA: %d", n_pcs_computed))

set.seed(42)  # Riproducibilita'
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars    = "roi_id",
  reduction.use    = "pca",
  reduction.save   = "harmony",
  dims.use         = seq_len(n_pcs_computed),
  theta            = 2,
  max_iter         = 20,
  plot_convergence = FALSE,
  verbose          = FALSE
)

message("Harmony completato.")
message("Embedding corretto disponibile in: seurat_obj[['harmony']]")

# --------------------------------------------------------------------------
# STEP 3e — Plot PCA post-Harmony
# --------------------------------------------------------------------------

harmony_coords <- as.data.frame(Embeddings(seurat_obj, "harmony")[idx_plot, 1:2])
colnames(harmony_coords) <- c("H1", "H2")
harmony_coords$roi_id <- seurat_obj$roi_id[idx_plot]

p_harmony_roi <- ggplot(harmony_coords, aes(x = H1, y = H2, colour = roi_id)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_colour_brewer(palette = "Set1") +
  labs(title    = "Embedding Harmony (colorato per ROI)",
       subtitle = "I ROI devono sovrapporsi se il batch effect e' stato corretto",
       x = "Harmony 1", y = "Harmony 2",
       colour = "ROI") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw()

ggsave(file.path(PLOT_DIR, "Integration_02_Harmony_post_correction.pdf"),
       p_harmony_roi, width = 8, height = 6)
message("  Salvato: Integration_02_Harmony_post_correction.pdf")

# Pannello comparativo pre vs post
p_comparison <- p_pca_roi + p_harmony_roi +
  plot_annotation(title = "Batch correction: PCA vs Harmony",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(PLOT_DIR, "Integration_03_comparison_pre_post.pdf"),
       p_comparison, width = 14, height = 6)
message("  Salvato: Integration_03_comparison_pre_post.pdf")

# --------------------------------------------------------------------------
# Salvataggio
# --------------------------------------------------------------------------

saveRDS(seurat_obj, file.path(DATA_DIR, "03_seurat_integrated.rds"))
message("Salvato: 03_seurat_integrated.rds")

message("=== STEP 3 completato ===\n")
