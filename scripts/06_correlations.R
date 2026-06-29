# =============================================================================
# 06_correlations.R
# Correlazioni single-cell: marcatori autofagia vs ligandi DNAM-1
#
# Input:  output/data/04_seurat_clustered.rds
# Output: output/data/06_spearman_correlations.csv
#         output/plots/Correlations_*.pdf
#
# DOMANDA BIOLOGICA (Reviewer #3, Major Comment 2):
#   Le cellule con alta attivita' autofagica hanno bassa espressione dei
#   ligandi DNAM-1 (PVR, Nectin2) nel tessuto primario?
#   Questa relazione e' presente in G3? E' diversa da SHH?
#
# INTERPRETAZIONE DEI MARCATORI:
#   LC3B  — marcatore di autofagosomi (alto = molti autofagosomi)
#   P62   — substrato autofagico (alto = blocco del flusso o alta produzione)
#   pMTOR — mTOR attivo = autofagia SOPPRESSA (attesa correlazione positiva
#            con ligandi: mTOR alto -> autofagia off -> ligandi non degradati)
#
# NOTA STATISTICA:
#   Con centinaia di migliaia di cellule, qualsiasi cor.test() producera'
#   p-value << 0.05 per ragioni puramente matematiche (potenza statistica
#   infinita). L'effetto biologico rilevante e' la dimensione di rho,
#   non la sua significativita'. Riportare sempre rho e n_cells.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== STEP 6: Correlazioni autofagia vs ligandi DNAM-1 ===")

seurat_obj  <- readRDS(file.path(OUT_DATA, "04_seurat_clustered.rds"))
SHORT_NAMES <- rownames(seurat_obj)
COFACTOR    <- seurat_obj@misc$cofactor

AUTOPHAGY_MARKERS <- c("LC3B", "P62", "pMTOR")
DNAM1_LIGANDS     <- c("PVR", "Nectin2")
PAIRS <- expand.grid(autophagy = AUTOPHAGY_MARKERS,
                     ligand    = DNAM1_LIGANDS,
                     stringsAsFactors = FALSE)

mat_all <- GetAssayData(seurat_obj, layer = "data", assay = "MICS")

# --------------------------------------------------------------------------
# STEP 6a — Correlazioni Spearman per paziente
# --------------------------------------------------------------------------

message("Calcolo correlazioni Spearman per paziente...")

corr_results <- list()

for (pid in unique(seurat_obj$patient_id)) {
  grp   <- unique(seurat_obj$group[seurat_obj$patient_id == pid])
  cells <- colnames(seurat_obj)[seurat_obj$patient_id == pid]

  for (i in seq_len(nrow(PAIRS))) {
    x_name <- PAIRS$autophagy[i]
    y_name <- PAIRS$ligand[i]
    x_vals <- as.numeric(mat_all[x_name, cells])
    y_vals <- as.numeric(mat_all[y_name, cells])

    # exact=FALSE obbligatorio per n > 5000 (evita calcolo esatto lento)
    ct <- cor.test(x_vals, y_vals, method = "spearman", exact = FALSE)

    corr_results[[length(corr_results) + 1]] <- data.frame(
      patient_id = pid,
      group      = grp,
      autophagy  = x_name,
      ligand     = y_name,
      pair       = paste0(x_name, " ~ ", y_name),
      rho        = round(ct$estimate, 4),
      p_value    = ct$p.value,
      n_cells    = length(cells),
      stringsAsFactors = FALSE
    )
  }
}
corr_df <- do.call(rbind, corr_results)
rownames(corr_df) <- NULL

write.csv(corr_df, file.path(OUT_DATA, "06_spearman_correlations.csv"),
          row.names = FALSE, quote = FALSE)
message("  Salvato: 06_spearman_correlations.csv")

message("\nCorrelazioni Spearman (rho) per paziente:")
corr_wide <- pivot_wider(corr_df[, c("patient_id","group","pair","rho")],
  names_from = "pair", values_from = "rho")
print(as.data.frame(corr_wide))

# --------------------------------------------------------------------------
# STEP 6b — Scatter plots per coppia (subsample per visualizzazione)
# --------------------------------------------------------------------------

message("\nGenerazione scatter plots...")

set.seed(SEED)
cells_sample <- unlist(lapply(unique(seurat_obj$patient_id), function(pid) {
  cells <- colnames(seurat_obj)[seurat_obj$patient_id == pid]
  sample(cells, min(5000, length(cells)))
}))

df_scatter <- cbind(
  data.frame(
    patient_id = seurat_obj$patient_id[cells_sample],
    group      = seurat_obj$group[cells_sample]
  ),
  t(as.matrix(mat_all[c(AUTOPHAGY_MARKERS, DNAM1_LIGANDS), cells_sample]))
)

# rho medio per gruppo (per annotazioni)
rho_by_group <- corr_df %>%
  group_by(group, autophagy, ligand) %>%
  summarise(rho_mean = round(mean(rho), 3), .groups = "drop") %>%
  mutate(pair = paste0(autophagy, " ~ ", ligand))

plot_list <- list()

for (i in seq_len(nrow(PAIRS))) {
  x_name <- PAIRS$autophagy[i]
  y_name <- PAIRS$ligand[i]

  rho_ann <- rho_by_group %>%
    filter(autophagy == x_name, ligand == y_name) %>%
    mutate(label = paste0("rho=", rho_mean))

  # FIX: usa .data[[]] invece di aes_string (deprecato da ggplot2 >= 3.0)
  p <- ggplot(df_scatter,
              aes(x = .data[[x_name]], y = .data[[y_name]], colour = group)) +
    geom_point(size = 0.3, alpha = 0.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, formula = y ~ x) +
    facet_wrap(~group) +
    scale_colour_manual(values = PALETTE_GROUP) +
    geom_text(data = rho_ann,
              aes(label = label, colour = group),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
              size = 3.5, inherit.aes = FALSE) +
    labs(
      title    = paste0(x_name, " vs ", y_name),
      subtitle = "Spearman rho (media per gruppo) su subsample 5k/paziente",
      x = sprintf("%s [arcsinh(MFI/%d)]", x_name, COFACTOR),
      y = sprintf("%s [arcsinh(MFI/%d)]", y_name, COFACTOR)
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")

  plot_list[[paste0(x_name, "_", y_name)]] <- p
}

p_pvr    <- plot_list[["LC3B_PVR"]] / plot_list[["P62_PVR"]] / plot_list[["pMTOR_PVR"]]
p_nectin <- plot_list[["LC3B_Nectin2"]] / plot_list[["P62_Nectin2"]] / plot_list[["pMTOR_Nectin2"]]

ggsave(file.path(OUT_PLOTS, "Correlations_01_autophagy_vs_PVR.pdf"),
       p_pvr, width = 10, height = 14)
ggsave(file.path(OUT_PLOTS, "Correlations_02_autophagy_vs_Nectin2.pdf"),
       p_nectin, width = 10, height = 14)
message("  Salvato: Correlations_01_autophagy_vs_PVR.pdf")
message("  Salvato: Correlations_02_autophagy_vs_Nectin2.pdf")

# --------------------------------------------------------------------------
# STEP 6c — Heatmap rho media per gruppo
# --------------------------------------------------------------------------

corr_mean <- corr_df %>%
  group_by(group, autophagy, ligand) %>%
  summarise(rho_mean = round(mean(rho), 3), .groups = "drop")

p_rho <- ggplot(corr_mean,
                aes(x = ligand, y = autophagy, fill = rho_mean)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = rho_mean), size = 4) +
  facet_wrap(~group) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho") +
  labs(title = "Correlazioni Spearman: autofagia vs ligandi DNAM-1",
       subtitle = "Media dei rho tra pazienti per gruppo",
       x = "Ligando DNAM-1", y = "Marcatore autofagia") +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_PLOTS, "Correlations_03_rho_heatmap_by_group.pdf"),
       p_rho, width = 8, height = 5)
message("  Salvato: Correlations_03_rho_heatmap_by_group.pdf")

# Heatmap rho per singolo paziente
corr_pat_wide <- pivot_wider(
  corr_df[, c("patient_id","group","pair","rho")],
  names_from = "pair", values_from = "rho"
)
mat_rho_pat <- as.matrix(corr_pat_wide[, -(1:2)])
rownames(mat_rho_pat) <- paste0(corr_pat_wide$patient_id, " (", corr_pat_wide$group, ")")

p_rho_pat <- ggplot(corr_df,
                    aes(x = pair, y = patient_id, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = rho), size = 3.2) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Spearman rho per paziente",
       x = "Coppia", y = "Paziente") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_PLOTS, "Correlations_04_rho_per_patient.pdf"),
       p_rho_pat, width = 10, height = max(5, length(unique(seurat_obj$patient_id)) * 0.8 + 2))
message("  Salvato: Correlations_04_rho_per_patient.pdf")

message("=== STEP 6 completato ===\n")
