# =============================================================================
# skewness_analysis.R
# Analisi distribuzione skewness LC3B e P62: G3 vs SHH
#
# Analisi standalone, separata dalla pipeline principale.
# Carica i file *Skewness*.csv da TUTTE le ROI disponibili per ogni paziente
# (incluse le ROI escluse dalla pipeline MFI per proteine mancanti,
#  poiche' il motivo di esclusione non si applica ai dati di skewness).
#
# Output: output/data/skewness_cells.csv
#         output/data/skewness_patient_medians.csv
#         output/plots/Skewness_01_ridge_per_paziente.pdf
#         output/plots/Skewness_02_violin_G3vsSHH.pdf
#         output/plots/Skewness_03_scatter_LC3B_vs_P62.pdf
#         output/plots/Skewness_04_patient_medians.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== SKEWNESS ANALYSIS: LC3B e P62 ===")

# Colonne target nel file skewness (nomi esatti come da MACS iQ View)
SKEW_COLS  <- c("LC3B Cytoplasm Intensity Skewness",
                "P62 Cytoplasm Intensity Skewness")
SHORT_SKEW <- c("LC3B_skewness", "P62_skewness")

# --------------------------------------------------------------------------
# Funzione: trova file skewness con entrambe le colonne in una ROI dir
# Gestisce i diversi naming usati da MACS iQ View
# --------------------------------------------------------------------------

find_skewness_file <- function(roi_dir) {
  if (!dir.exists(roi_dir)) return(NULL)
  files <- list.files(roi_dir, pattern = "Skewness", full.names = TRUE,
                      ignore.case = TRUE)
  if (length(files) == 0) return(NULL)
  for (f in files) {
    header <- tryCatch(colnames(fread(f, nrows = 0)),
                       error = function(e) character(0))
    if (all(SKEW_COLS %in% header)) return(f)
  }
  return(NULL)
}

# --------------------------------------------------------------------------
# Caricamento: scansiona TUTTE le ROI di ogni paziente in DATA_DIR
# (non solo quelle in SAMPLES, per includere ROI escluse dalla pipeline MFI)
# --------------------------------------------------------------------------

# Recupera patient_id e gruppo da SAMPLES
patient_info <- unique(do.call(rbind, lapply(SAMPLES, function(s) {
  data.frame(patient_id = s$patient_id, group = s$group,
             stringsAsFactors = FALSE)
})))

df_list <- list()
skipped <- character(0)

for (i in seq_len(nrow(patient_info))) {
  pid <- patient_info$patient_id[i]
  grp <- patient_info$group[i]

  pdir <- file.path(DATA_DIR, pid)
  if (!dir.exists(pdir)) {
    message(sprintf("  [SKIP] Cartella non trovata: %s", pdir))
    next
  }

  roi_dirs <- list.dirs(pdir, recursive = FALSE, full.names = TRUE)
  message(sprintf("\n%s (%s):", pid, grp))

  for (rdir in roi_dirs) {
    roi_name  <- basename(rdir)
    skew_file <- find_skewness_file(rdir)

    if (is.null(skew_file)) {
      message(sprintf("  [SKIP] %s: nessun file skewness con LC3B+P62", roi_name))
      skipped <- c(skipped, paste0(pid, "_", roi_name))
      next
    }

    message(sprintf("  %s -> %s", roi_name, basename(skew_file)))

    dt <- tryCatch(
      fread(skew_file, select = SKEW_COLS, header = TRUE, data.table = FALSE),
      error = function(e) {
        message(sprintf("     [ERR] lettura fallita: %s", e$message))
        NULL
      }
    )
    if (is.null(dt) || nrow(dt) == 0) next

    colnames(dt) <- SHORT_SKEW
    dt <- as.data.frame(lapply(dt, as.numeric))

    na_mask <- rowSums(is.na(dt)) > 0
    if (any(na_mask)) {
      message(sprintf("     [WARN] %d righe con NA rimosse", sum(na_mask)))
      dt <- dt[!na_mask, , drop = FALSE]
    }
    if (nrow(dt) == 0) next

    n_cells  <- nrow(dt)
    cell_ids <- paste0(pid, "_", roi_name, "_cell_",
                       formatC(seq_len(n_cells), width = nchar(n_cells), flag = "0"))
    rownames(dt) <- cell_ids
    dt$patient_id <- pid
    dt$group      <- grp
    dt$roi_id     <- paste0(pid, "_", roi_name)

    df_list[[paste0(pid, "_", roi_name)]] <- dt
    message(sprintf("     -> OK (%d cellule)", n_cells))
  }
}

if (length(df_list) == 0) stop("Nessun dato skewness trovato.")

df_all <- do.call(rbind, df_list)
df_all$group <- factor(df_all$group, levels = c("G3", "SHH"))

message(sprintf("\nTotale cellule: %d da %d ROI", nrow(df_all), length(df_list)))
if (length(skipped) > 0)
  message(sprintf("ROI senza skewness: %d (%s)",
                  length(skipped), paste(skipped, collapse = ", ")))

message("\nCellule per paziente:")
print(as.data.frame(df_all %>%
  group_by(group, patient_id) %>%
  summarise(n_cells = n(), n_roi = n_distinct(roi_id), .groups = "drop")))

# --------------------------------------------------------------------------
# Salva dati cell-level
# --------------------------------------------------------------------------

write.csv(df_all, file.path(OUT_DATA, "skewness_cells.csv"),
          row.names = FALSE, quote = FALSE)
message("Salvato: skewness_cells.csv")

# --------------------------------------------------------------------------
# Palette
# --------------------------------------------------------------------------

g3_pats  <- sort(unique(df_all$patient_id[df_all$group == "G3"]))
shh_pats <- sort(unique(df_all$patient_id[df_all$group == "SHH"]))
palette_patient <- c(
  setNames(PALETTE_PATIENTS_G3[seq_along(g3_pats)],   g3_pats),
  setNames(PALETTE_PATIENTS_SHH[seq_along(shh_pats)], shh_pats)
)

# Ordine pazienti per plot: G3 in alto, SHH in basso
patient_order <- c(rev(g3_pats), rev(shh_pats))
df_all$patient_id <- factor(df_all$patient_id, levels = patient_order)

# --------------------------------------------------------------------------
# PLOT 1 — Ridge plot: distribuzione skewness per paziente
# --------------------------------------------------------------------------

df_long <- pivot_longer(df_all,
  cols = SHORT_SKEW, names_to = "protein", values_to = "skewness")
df_long$protein <- factor(df_long$protein,
  levels = SHORT_SKEW,
  labels = c("LC3B Skewness", "P62 Skewness"))

p_ridge <- ggplot(df_long,
    aes(x = skewness, y = patient_id, fill = group)) +
  geom_density_ridges(alpha = 0.75, scale = 1.2, rel_min_height = 0.005) +
  facet_wrap(~protein, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = PALETTE_GROUP) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  labs(
    title    = "LC3B and P62 skewness distribution by patient",
    subtitle = "Skewness > 0 = punctate distribution; skewness \u2248 0 = diffuse distribution",
    x = "Skewness", y = "Patient", fill = "Group"
  ) +
  theme_ridges(font_size = 10) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey90"))

ggsave(file.path(OUT_PLOTS, "Skewness_01_ridge_per_paziente.pdf"),
       p_ridge, width = 14, height = 8)
message("Salvato: Skewness_01_ridge_per_paziente.pdf")

# --------------------------------------------------------------------------
# PLOT 2 — Violin: G3 vs SHH (subsample per velocita')
# --------------------------------------------------------------------------

set.seed(SEED)
MAX_VIOLIN <- 200000
df_sub <- if (nrow(df_all) > MAX_VIOLIN) {
  df_all[sample(seq_len(nrow(df_all)), MAX_VIOLIN), ]
} else df_all

df_long_sub <- pivot_longer(df_sub,
  cols = SHORT_SKEW, names_to = "protein", values_to = "skewness")
df_long_sub$protein <- factor(df_long_sub$protein,
  levels = SHORT_SKEW,
  labels = c("LC3B Skewness", "P62 Skewness"))

p_violin <- ggplot(df_long_sub, aes(x = group, y = skewness, fill = group)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  facet_wrap(~protein, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = PALETTE_GROUP) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  labs(
    title = "LC3B and P62 skewness: G3 vs SHH (single-cell level)",
    x = "Group", y = "Skewness", fill = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"))

ggsave(file.path(OUT_PLOTS, "Skewness_02_violin_G3vsSHH.pdf"),
       p_violin, width = 10, height = 6)
message("Salvato: Skewness_02_violin_G3vsSHH.pdf")

# --------------------------------------------------------------------------
# PLOT 3 — Scatter 2D: LC3B skewness vs P62 skewness con densita'
# --------------------------------------------------------------------------

set.seed(SEED)
MAX_SCATTER <- 100000
df_sc <- if (nrow(df_all) > MAX_SCATTER) {
  df_all[sample(seq_len(nrow(df_all)), MAX_SCATTER), ]
} else df_all

p_scatter <- ggplot(df_sc, aes(x = LC3B_skewness, y = P62_skewness,
                                colour = group)) +
  geom_point(size = 0.08, alpha = 0.15) +
  geom_density_2d(linewidth = 0.5, alpha = 0.9) +
  scale_colour_manual(values = PALETTE_GROUP) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.3) +
  labs(
    title    = "LC3B skewness vs P62 skewness",
    subtitle = sprintf("Subsample: %s cells; contours = density per group",
                       format(nrow(df_sc), big.mark = ",")),
    x = "LC3B cytoplasm skewness",
    y = "P62 cytoplasm skewness",
    colour = "Group"
  ) +
  theme_bw(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(OUT_PLOTS, "Skewness_03_scatter_LC3B_vs_P62.pdf"),
       p_scatter, width = 8, height = 7)
message("Salvato: Skewness_03_scatter_LC3B_vs_P62.pdf")

# --------------------------------------------------------------------------
# PLOT 4 — Mediane per paziente (unita' biologica)
# --------------------------------------------------------------------------

patient_medians <- df_all %>%
  group_by(group, patient_id) %>%
  summarise(
    LC3B_skewness_median = median(LC3B_skewness, na.rm = TRUE),
    P62_skewness_median  = median(P62_skewness,  na.rm = TRUE),
    n_cells              = n(),
    .groups = "drop"
  )

message("\nMediane skewness per paziente:")
print(as.data.frame(patient_medians))

df_med_long <- pivot_longer(patient_medians,
  cols      = c("LC3B_skewness_median", "P62_skewness_median"),
  names_to  = "protein",
  values_to = "median_skewness")
df_med_long$protein <- factor(df_med_long$protein,
  levels = c("LC3B_skewness_median", "P62_skewness_median"),
  labels = c("LC3B Skewness", "P62 Skewness"))

p_medians <- ggplot(df_med_long,
    aes(x = group, y = median_skewness, colour = group)) +
  geom_jitter(aes(shape = patient_id), size = 3.5, width = 0.1) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.35, colour = "black", linewidth = 0.6) +
  facet_wrap(~protein, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_shape_manual(values = c(15, 16, 17, 18, 4, 8, 9, 10, 11, 12)[
    seq_along(levels(df_all$patient_id))]) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  labs(
    title    = "Median skewness per patient: G3 vs SHH",
    subtitle = "Each point = one patient; bar = group mean",
    x = "Group", y = "Median skewness",
    colour = "Group", shape = "Patient"
  ) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"))

ggsave(file.path(OUT_PLOTS, "Skewness_04_patient_medians.pdf"),
       p_medians, width = 10, height = 6)
message("Salvato: Skewness_04_patient_medians.pdf")

write.csv(patient_medians, file.path(OUT_DATA, "skewness_patient_medians.csv"),
          row.names = FALSE, quote = FALSE)
message("Salvato: skewness_patient_medians.csv")

message("\n=== SKEWNESS ANALYSIS completata ===")
message(sprintf("Plot in: %s", OUT_PLOTS))
