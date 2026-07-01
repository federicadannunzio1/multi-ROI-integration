# =============================================================================
# setup_packages.R
# Check required packages for the MACSima pipeline and install any missing ones.
#
# USAGE:
#   Run this script once before the first pipeline execution.
#   If working in a conda environment with Seurat pre-installed, activate it first:
#
#   conda activate <your_env_name>
#   Rscript setup_packages.R
#
# By default, missing packages are installed into the R user library detected
# automatically from the R_LIBS_USER environment variable (or ~/.R/library).
# Override by setting R_LIBS_USER before running:
#   export R_LIBS_USER=/your/custom/R/library
#   Rscript setup_packages.R
# =============================================================================

# --------------------------------------------------------------------------
# 0. Configure user library path
# --------------------------------------------------------------------------

USER_LIB <- Sys.getenv(
  "R_LIBS_USER",
  unset = file.path(Sys.getenv("HOME"), "R", "library")
)

if (!dir.exists(USER_LIB)) {
  dir.create(USER_LIB, recursive = TRUE)
  message("Created user library: ", USER_LIB)
}

# Prepend to .libPaths() if not already present
if (!USER_LIB %in% .libPaths()) {
  .libPaths(c(USER_LIB, .libPaths()))
}
message("Active library paths:")
for (lp in .libPaths()) message("  ", lp)

# --------------------------------------------------------------------------
# 1. Required packages
# --------------------------------------------------------------------------

# CRAN packages
CRAN_PKGS <- c(
  "data.table",    # 01_load_data.R
  "dplyr",         # all steps
  "tidyr",         # steps 2, 5, 6, 7
  "ggplot2",       # steps 2–7
  "scales",        # step 7 - color scales
  "patchwork",     # steps 3, 4, 5, 6, 7
  "RColorBrewer",  # steps 4, 7
  "clustree",      # step 4 — installs ggraph as dependency
  "pheatmap",      # steps 4, 5, 7
  "ggridges",      # step 2 — arcsinh distribution ridge plot
  "harmony",       # step 3
  "Seurat",        # steps 2–7 (installs SeuratObject automatically)
  "igraph"         # dependency of clustree and Seurat
)

# Bioconductor packages (indirect dependencies of Seurat v5)
BIOC_PKGS <- c(
  "BiocGenerics",
  "GenomicRanges",
  "IRanges",
  "S4Vectors"
)

# --------------------------------------------------------------------------
# 2. Check availability
# --------------------------------------------------------------------------

check_packages <- function(pkgs) {
  sapply(pkgs, function(p) {
    installed <- requireNamespace(p, quietly = TRUE)
    list(installed = installed,
         version   = if (installed) as.character(packageVersion(p)) else NA_character_)
  }, simplify = FALSE)
}

cat("\n--- Checking CRAN packages ---\n")
cran_status  <- check_packages(CRAN_PKGS)
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

cat("\n--- Checking Bioconductor packages ---\n")
bioc_status  <- check_packages(BIOC_PKGS)
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
# 3. Install missing packages
# --------------------------------------------------------------------------

CRAN_MIRROR <- "https://cloud.r-project.org"

if (length(cran_missing) > 0) {
  cat(sprintf("\n--- Installing %d missing CRAN packages ---\n", length(cran_missing)))
  cat("  Packages: ", paste(cran_missing, collapse = ", "), "\n\n")

  install.packages(
    cran_missing,
    lib          = USER_LIB,
    repos        = CRAN_MIRROR,
    Ncpus        = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))),
    dependencies = TRUE
  )
} else {
  cat("\nAll CRAN packages are already installed.\n")
}

if (length(bioc_missing) > 0) {
  cat(sprintf("\n--- Installing %d missing Bioconductor packages ---\n",
              length(bioc_missing)))

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = USER_LIB, repos = CRAN_MIRROR)
  }
  BiocManager::install(bioc_missing, lib = USER_LIB, update = FALSE, ask = FALSE)
} else {
  cat("All Bioconductor packages are already installed.\n")
}

# --------------------------------------------------------------------------
# 4. Final check: load all packages to confirm they work
# --------------------------------------------------------------------------

cat("\n--- Final load check (library()) ---\n")

all_pkgs  <- c(CRAN_PKGS, BIOC_PKGS)
load_ok   <- c()
load_fail <- c()

for (pkg in all_pkgs) {
  ok <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    TRUE
  }, error = function(e) FALSE)

  if (ok) {
    load_ok   <- c(load_ok, pkg)
    cat(sprintf("  [LOAD OK]   %s\n", pkg))
  } else {
    load_fail <- c(load_fail, pkg)
    cat(sprintf("  [LOAD FAIL] %s\n", pkg))
  }
}

# --------------------------------------------------------------------------
# 5. Summary
# --------------------------------------------------------------------------

cat("\n=============================================\n")
cat(sprintf("Successfully loaded: %d / %d\n", length(load_ok), length(all_pkgs)))
if (length(load_fail) > 0) {
  cat("Failed packages:\n")
  for (p in load_fail) cat("  -", p, "\n")
  cat("\nPossible causes:\n")
  cat("  - Missing system libraries (libcurl, openssl, libxml2, zlib)\n")
  cat("    Try loading the relevant modules on your cluster, e.g.:\n")
  cat("      module load curl openssl\n")
  cat("  - Incompatible R version (Seurat v5 requires R >= 4.1)\n")
  cat("    Check with:  R.version\n")
  cat("  - Package not compiled for the current architecture\n")
  cat("    Reinstall from source:\n")
  cat("      install.packages('", paste(load_fail, collapse = "','"),
      "', type = 'source')\n", sep = "")
} else {
  cat("All packages are ready. You can now run the pipeline.\n")
}
cat("=============================================\n")
