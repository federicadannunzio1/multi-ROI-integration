# =============================================================================
# config.R
# Configurazione centralizzata — tutti gli script fanno source() di questo file
#
# Modifica SOLO questo file per adattare la pipeline a un nuovo ambiente.
# =============================================================================

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

BASE_DIR  <- "/lustre/home/gfiscon/projects/MACSima_pipeline"
DATA_DIR  <- file.path(BASE_DIR, "data")
OUT_DATA  <- file.path(BASE_DIR, "output", "data")
OUT_PLOTS <- file.path(BASE_DIR, "output", "plots")
OUT_LOGS  <- file.path(BASE_DIR, "logs")

# Crea output directories se non esistono
for (d in c(OUT_DATA, OUT_PLOTS, OUT_LOGS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --------------------------------------------------------------------------
# Proteine
# --------------------------------------------------------------------------

RAW_PROTEIN_NAMES <- c(
  "Ki_67 REA1123 Biomarker Exp",
  "LC3B Cyto Exp",
  "Nectin 2 Biomarker Exp",
  "P62 Cyto Exp",
  "Phospo mTor 49F9 Cyto Exp",
  "Poliovirus Receptor Biomarker Exp"
)
SHORT_NAMES <- c("Ki67", "LC3B", "Nectin2", "P62", "pMTOR", "PVR")
PROTEIN_MAP <- setNames(SHORT_NAMES, RAW_PROTEIN_NAMES)

# --------------------------------------------------------------------------
# Campioni
# Struttura: patient_id, group (G3/SHH), lista di ROI (nome -> cartella)
#
# NOTE SUI ROI INCOMPLETI:
#   Alcuni ROI mancano di LC3B e P62 nell'export MACS iQ View.
#   Poiche' questi marcatori sono centrali all'analisi (autofagia),
#   i ROI incompleti sono ESCLUSI. Solo i ROI con tutte e 6 le proteine
#   sono inclusi qui sotto.
#
#   Gr3_25-1580  ROI1, ROI2 -> esclusi (mancano LC3B, P62); ROI3 -> incluso
#   SHH_24-2143  ROI2       -> escluso (mancano LC3B, P62); ROI1 -> incluso
#
# NOTE GENERALI:
#   - SHH_24-8477 ROI4: incluso (dati ri-esportati da MACS iQ View)
# --------------------------------------------------------------------------

SAMPLES <- list(

  # ── G3 (5 pazienti) ───────────────────────────────────────────────────────
  list(
    patient_id = "Gr3_23-3017",
    group      = "G3",
    rois       = list(
      ROI1 = file.path(DATA_DIR, "Gr3_23-3017", "ROI1"),
      ROI2 = file.path(DATA_DIR, "Gr3_23-3017", "ROI2")
    )
  ),
  list(
    patient_id = "Gr3_23-3106",
    group      = "G3",
    rois       = list(
      ROI1 = file.path(DATA_DIR, "Gr3_23-3106", "ROI1"),
      ROI2 = file.path(DATA_DIR, "Gr3_23-3106", "ROI2")
    )
  ),
  list(
    patient_id = "Gr3_24-0268",
    group      = "G3",
    rois       = list(
      ROI1 = file.path(DATA_DIR, "Gr3_24-0268", "ROI1"),
      ROI2 = file.path(DATA_DIR, "Gr3_24-0268", "ROI2"),
      ROI3 = file.path(DATA_DIR, "Gr3_24-0268", "ROI3"),
      ROI4 = file.path(DATA_DIR, "Gr3_24-0268", "ROI4")
    )
  ),
  list(
    patient_id = "Gr3_25-1580",
    group      = "G3",
    rois       = list(
      # ROI1, ROI2 esclusi: mancano LC3B e P62
      ROI3 = file.path(DATA_DIR, "Gr3_25-1580", "ROI3")
    )
  ),
  list(
    patient_id = "Gr3_25-7278",
    group      = "G3",
    rois       = list(
      ROI1 = file.path(DATA_DIR, "Gr3_25-7278", "ROI1"),
      ROI2 = file.path(DATA_DIR, "Gr3_25-7278", "ROI2"),
      ROI3 = file.path(DATA_DIR, "Gr3_25-7278", "ROI3"),
      ROI4 = file.path(DATA_DIR, "Gr3_25-7278", "ROI4")
    )
  ),

  # ── SHH (5 pazienti) ──────────────────────────────────────────────────────
  list(
    patient_id = "SHH_22-6172",
    group      = "SHH",
    rois       = list(
      ROI1 = file.path(DATA_DIR, "SHH_22-6172", "ROI1"),
      ROI2 = file.path(DATA_DIR, "SHH_22-6172", "ROI2"),
      ROI3 = file.path(DATA_DIR, "SHH_22-6172", "ROI3")
    )
  ),
  list(
    patient_id = "SHH_23-9574",
    group      = "SHH",
    rois       = list(
      ROI1 = file.path(DATA_DIR, "SHH_23-9574", "ROI1"),
      ROI2 = file.path(DATA_DIR, "SHH_23-9574", "ROI2")
    )
  ),
  list(
    patient_id = "SHH_24-2143",
    group      = "SHH",
    rois       = list(
      ROI1 = file.path(DATA_DIR, "SHH_24-2143", "ROI1")
      # ROI2 escluso: mancano LC3B e P62
    )
  ),
  list(
    patient_id = "SHH_24-8477",
    group      = "SHH",
    rois       = list(
      ROI1 = file.path(DATA_DIR, "SHH_24-8477", "ROI1"),
      ROI2 = file.path(DATA_DIR, "SHH_24-8477", "ROI2"),
      ROI3 = file.path(DATA_DIR, "SHH_24-8477", "ROI3"),
      ROI4 = file.path(DATA_DIR, "SHH_24-8477", "ROI4")
    )
  ),
  # SHH_25-6667 ESCLUSO — artefatto tecnico confermato:
  #   - P62 MFI raw = 0 su tutte le cellule (arcsinh = 0), biologicamente impossibile
  #     (P62 e' un recettore autofagico costitutivamente espresso)
  #   - LC3B MFI = 0.063 (6-44x inferiore a tutti gli altri pazienti, range: 2.8-9.3)
  #   - LC3B skewness = 3.29 (tutti gli altri pazienti < 0.30)
  #   Probabile failure dell'export MACS iQ View per questi canali.
  #   Confermato con collaboratori; escluso da tutte le analisi.
)

# --------------------------------------------------------------------------
# Parametri analisi
# --------------------------------------------------------------------------

# Clustering
CLUSTERING_RESOLUTIONS <- c(0.1)  # resolution validata con clustree nella run esplorativa
FINAL_RESOLUTION       <- 0.1    # confermata: struttura stabile, biologicamente interpretabile
KNN_K                  <- 30                        # k per FindNeighbors (standard Seurat)

# Harmony
HARMONY_THETA    <- 2     # forza della batch correction
HARMONY_MAX_ITER <- 20    # iterazioni massime

# Colori fissi per gruppo
PALETTE_GROUP <- c(G3 = "#C0392B", SHH = "#1A5276")

# Palette per paziente — colori qualitativi distinguibili, separati per gruppo.
# Toni caldi per G3 (rosso/arancione/viola), toni freddi per SHH (blu/verde).
# Supporta fino a 5 pazienti per gruppo.
PALETTE_PATIENTS_G3  <- c("#C0392B", "#E67E22", "#F4D03F", "#884EA0", "#D35400")
PALETTE_PATIENTS_SHH <- c("#1A5276", "#2E86C1", "#148F77", "#28B463")

# Seed riproducibilita'
SEED <- 42

message("config.R caricato — BASE_DIR: ", BASE_DIR)
