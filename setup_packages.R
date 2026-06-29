# =============================================================================
# setup_packages.R
# Verifica i pacchetti richiesti dalla pipeline MACSima e installa i mancanti.
#
# UTILIZZO SUL CLUSTER (TeraStat2, Sapienza):
#   Seurat e' installato nell'ambiente conda 'seurat_env'.
#   Attivare l'ambiente PRIMA di eseguire questo script:
#
#   source /lustre/software/anaconda/2022.10_all/etc/profile.d/conda.sh
#   conda activate seurat_env
#   Rscript setup_packages.R
#
# I pacchetti mancanti verranno installati nella libreria utente:
#   /lustre/home/gfiscon/R/x86_64-pc-linux-gnu-library/4.4/
# =============================================================================

# --------------------------------------------------------------------------
# 0. Configura percorso libreria utente
# --------------------------------------------------------------------------

USER_LIB <- "/lustre/home/gfiscon/R/x86_64-pc-linux-gnu-library/4.4"
if (!dir.exists(USER_LIB)) {
  dir.create(USER_LIB, recursive = TRUE)
  message("Creata libreria utente: ", USER_LIB)
}

# Aggiungi in testa a .libPaths() se non c'e' gia'
if (!USER_LIB %in% .libPaths()) {
  .libPaths(c(USER_LIB, .libPaths()))
}
message("Librerie attive:")
for (lp in .libPaths()) message("  ", lp)

# --------------------------------------------------------------------------
# 1. Pacchetti richiesti dalla pipeline
# --------------------------------------------------------------------------

# Pacchetti CRAN
CRAN_PKGS <- c(
  "data.table",    # 01_load_data.R
  "dplyr",         # tutti gli step
  "tidyr",         # step 2, 5, 6, 7
  "ggplot2",       # step 2-7
  "patchwork",     # step 3, 4, 5, 6, 7
  "RColorBrewer",  # step 4, 7
  "clustree",      # step 4 - richiede ggraph (installato come dipendenza)
  "pheatmap",      # step 4, 5, 7
  "ggridges",      # step 2 - ridge plot distribuzione arcsinh
  "writexl",       # step 2 - export tabella cofactor quantili
  "harmony",       # step 3
  "Seurat",        # step 2-7  (installa automaticamente SeuratObject)
  "igraph"         # dipendenza di clustree e Seurat
)

# Pacchetti Bioconductor (dipendenze indirette di Seurat v5)
BIOC_PKGS <- c(
  "BiocGenerics",
  "GenomicRanges",
  "IRanges",
  "S4Vectors"
)

# --------------------------------------------------------------------------
# 2. Verifica disponibilita'
# --------------------------------------------------------------------------

check_packages <- function(pkgs) {
  status <- sapply(pkgs, function(p) {
    installed <- requireNamespace(p, quietly = TRUE)
    list(installed = installed,
         version   = if (installed) as.character(packageVersion(p)) else NA_character_)
  }, simplify = FALSE)
  status
}

cat("\n--- Verifica pacchetti CRAN ---\n")
cran_status <- check_packages(CRAN_PKGS)
cran_missing <- c()
for (pkg in names(cran_status)) {
  s <- cran_status[[pkg]]
  if (s$installed) {
    cat(sprintf("  [OK]      %-20s  v%s\n", pkg, s$version))
  } else {
    cat(sprintf("  [MISSING] %-20s\n", pkg))
    cran_missing <- c(cran_missing, pkg)
  }
}

cat("\n--- Verifica pacchetti Bioconductor ---\n")
bioc_status <- check_packages(BIOC_PKGS)
bioc_missing <- c()
for (pkg in names(bioc_status)) {
  s <- bioc_status[[pkg]]
  if (s$installed) {
    cat(sprintf("  [OK]      %-20s  v%s\n", pkg, s$version))
  } else {
    cat(sprintf("  [MISSING] %-20s\n", pkg))
    bioc_missing <- c(bioc_missing, pkg)
  }
}

# --------------------------------------------------------------------------
# 3. Installazione pacchetti mancanti
# --------------------------------------------------------------------------

CRAN_MIRROR <- "https://cloud.r-project.org"

if (length(cran_missing) > 0) {
  cat(sprintf("\n--- Installazione %d pacchetti CRAN mancanti ---\n",
              length(cran_missing)))
  cat("  Pacchetti: ", paste(cran_missing, collapse = ", "), "\n\n")

  install.packages(
    cran_missing,
    lib      = USER_LIB,
    repos    = CRAN_MIRROR,
    Ncpus    = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))),
    dependencies = TRUE
  )
} else {
  cat("\nTutti i pacchetti CRAN sono gia' installati.\n")
}

if (length(bioc_missing) > 0) {
  cat(sprintf("\n--- Installazione %d pacchetti Bioconductor mancanti ---\n",
              length(bioc_missing)))

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = USER_LIB, repos = CRAN_MIRROR)
  }
  BiocManager::install(bioc_missing, lib = USER_LIB, update = FALSE, ask = FALSE)
} else {
  cat("Tutti i pacchetti Bioconductor sono gia' installati.\n")
}

# --------------------------------------------------------------------------
# 4. Verifica finale: carica tutti i pacchetti per confermare che funzionino
# --------------------------------------------------------------------------

cat("\n--- Verifica finale (library()) ---\n")

all_pkgs   <- c(CRAN_PKGS, BIOC_PKGS)
load_ok    <- c()
load_fail  <- c()

for (pkg in all_pkgs) {
  ok <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    TRUE
  }, error = function(e) FALSE)

  if (ok) {
    load_ok   <- c(load_ok, pkg)
    cat(sprintf("  [LOAD OK] %s\n", pkg))
  } else {
    load_fail <- c(load_fail, pkg)
    cat(sprintf("  [LOAD FAIL] %s  <-- controlla dipendenze di sistema\n", pkg))
  }
}

# --------------------------------------------------------------------------
# 5. Riepilogo
# --------------------------------------------------------------------------

cat("\n=============================================\n")
cat(sprintf("Caricati correttamente: %d / %d\n", length(load_ok), length(all_pkgs)))
if (length(load_fail) > 0) {
  cat("FALLITI:\n")
  for (p in load_fail) cat("  -", p, "\n")
  cat("\nPossibili cause:\n")
  cat("  - Dipendenze di sistema mancanti (libcurl, openssl, libxml2, zlib)\n")
  cat("    Prova:  module load curl openssl  (o equivalente sul tuo cluster)\n")
  cat("  - Versione R incompatibile (Seurat v5 richiede R >= 4.1)\n")
  cat("    Verifica con:  R.version\n")
  cat("  - Pacchetto non compilato per l'architettura del nodo\n")
  cat("    Riinstalla con:  install.packages('", paste(load_fail, collapse="','"),
      "', type='source')\n", sep="")
} else {
  cat("Tutti i pacchetti sono pronti. Puoi lanciare la pipeline.\n")
}
cat("=============================================\n")
