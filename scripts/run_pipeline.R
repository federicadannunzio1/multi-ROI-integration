# =============================================================================
# run_pipeline.R
# Script master: esegue tutti gli step della pipeline in sequenza
#
# Uso:
#   Rscript scripts_cluster/run_pipeline.R
#   Rscript scripts_cluster/run_pipeline.R --steps 1,2,3   # solo step selezionati
#   Rscript scripts_cluster/run_pipeline.R --from 3         # riparte da step 3
#
# Ogni step viene eseguito in un ambiente separato (sys.source) per evitare
# conflitti tra variabili globali. Il log di ogni step e' scritto sia su
# stdout (catturato da SLURM) sia su file in OUT_LOGS/.
# =============================================================================

# --------------------------------------------------------------------------
# Parsing argomenti da riga di comando
# --------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

parse_steps <- function(args) {
  steps_all <- 1:7
  steps_flag  <- grep("^--steps$",  args)
  from_flag   <- grep("^--from$",   args)

  if (length(steps_flag) > 0) {
    val <- args[steps_flag + 1]
    return(as.integer(unlist(strsplit(val, ","))))
  }
  if (length(from_flag) > 0) {
    val <- as.integer(args[from_flag + 1])
    return(val:7)
  }
  return(steps_all)
}

STEPS_TO_RUN <- parse_steps(args)

STEP_SCRIPTS <- c(
  "1" = "01_load_data.R",
  "2" = "02_qc_normalization.R",
  "3" = "03_integration.R",
  "4" = "04_clustering.R",
  "5" = "05_comparison_G3_SHH.R",
  "6" = "06_correlations.R",
  "7" = "07_visualization.R"
)

# dirname(sys.frame(1)$ofile) funziona solo con source(), non con Rscript.
# commandArgs() e' il metodo corretto per script lanciati direttamente.
args_full   <- commandArgs(trailingOnly = FALSE)
file_flag   <- grep("--file=", args_full, value = TRUE)
SCRIPTS_DIR <- if (length(file_flag) > 0) {
  dirname(normalizePath(sub("--file=", "", file_flag[1])))
} else {
  getwd()
}

# --------------------------------------------------------------------------
# Setup logging
# --------------------------------------------------------------------------

# Carica config solo per OUT_LOGS
source(file.path(SCRIPTS_DIR, "config.R"))

log_file <- file.path(OUT_LOGS,
  sprintf("pipeline_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))

log_message <- function(msg) {
  ts  <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  line <- paste(ts, msg)
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}

# --------------------------------------------------------------------------
# Esecuzione step
# --------------------------------------------------------------------------

timings <- list()

log_message("======================================================")
log_message("MACSima pipeline — START")
log_message(sprintf("  Steps da eseguire: %s", paste(STEPS_TO_RUN, collapse = ", ")))
log_message(sprintf("  R version: %s", R.version.string))
log_message(sprintf("  Log: %s", log_file))
log_message("======================================================")

# Rende SCRIPTS_DIR disponibile agli step script come variabile d'ambiente
# (fallback quando sys.frame(1)$ofile non e' disponibile via sys.source)
Sys.setenv(MACSIMA_SCRIPTS_DIR = SCRIPTS_DIR)

pipeline_start <- proc.time()

for (step in STEPS_TO_RUN) {
  script_name <- STEP_SCRIPTS[as.character(step)]

  if (is.na(script_name)) {
    log_message(sprintf("WARNING: step %d non riconosciuto, saltato.", step))
    next
  }

  script_path <- file.path(SCRIPTS_DIR, script_name)

  if (!file.exists(script_path)) {
    stop(sprintf("Script non trovato: %s", script_path))
  }

  log_message(sprintf("--- STEP %d: %s ---", step, script_name))
  t_start <- proc.time()

  tryCatch({
    # Ogni step gira in un environment pulito per evitare contaminazione
    # tra variabili dei diversi step
    step_env <- new.env(parent = globalenv())
    sys.source(script_path, envir = step_env, keep.source = TRUE)
  }, error = function(e) {
    log_message(sprintf("ERRORE in step %d (%s):", step, script_name))
    log_message(sprintf("  %s", conditionMessage(e)))
    stop(sprintf("Pipeline interrotta allo step %d. Controlla i log.", step))
  })

  elapsed <- (proc.time() - t_start)[["elapsed"]]
  timings[[script_name]] <- elapsed
  log_message(sprintf("  Completato in %.1f min", elapsed / 60))
}

# --------------------------------------------------------------------------
# Riepilogo finale
# --------------------------------------------------------------------------

total_elapsed <- (proc.time() - pipeline_start)[["elapsed"]]

log_message("======================================================")
log_message("MACSima pipeline — COMPLETATA")
log_message(sprintf("  Tempo totale: %.1f min", total_elapsed / 60))
log_message("  Tempi per step:")
for (nm in names(timings)) {
  log_message(sprintf("    %-35s %.1f min", nm, timings[[nm]] / 60))
}
log_message("======================================================")
