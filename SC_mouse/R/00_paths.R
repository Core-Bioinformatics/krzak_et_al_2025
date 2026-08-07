# Shared path configuration for this analysis kit.
#
# Defaults assume this file lives in <kit>/R/ and inputs live in <kit>/data/.
# Override with environment variables when needed:
#   SCRNA_METHODS_ROOT  — kit root (directory containing R/, data/, output/)
#   SCRNA_DATA_DIR      — Seurat / ClustAssess RDS and small auxiliaries
#   SCRNA_OUT_DIR       — analysis outputs
#   SCRNA_SEURAT_RDS    — Seurat RDS filename inside DATA_DIR (default: seurat_object.rds)
#   SCRNA_CLUSTASSESS_RDS — ClustAssess RDS filename (default: clustassess_object.rds)

.methods_publish_this_file <- (function() {
    ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(ofile) && nzchar(ofile)) {
        return(normalizePath(ofile, winslash = "/", mustWork = FALSE))
    }
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd_args, value = TRUE)
    if (length(file_arg) == 1) {
        return(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = FALSE))
    }
    NA_character_
})()

if (!is.na(.methods_publish_this_file) && nzchar(.methods_publish_this_file)) {
    METHODS_R_DIR <- dirname(.methods_publish_this_file)
    METHODS_ROOT_DEFAULT <- dirname(METHODS_R_DIR)
} else {
    METHODS_R_DIR <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    METHODS_ROOT_DEFAULT <- METHODS_R_DIR
}

METHODS_ROOT <- Sys.getenv("SCRNA_METHODS_ROOT", unset = METHODS_ROOT_DEFAULT)
METHODS_ROOT <- normalizePath(METHODS_ROOT, winslash = "/", mustWork = FALSE)

DATA_DIR <- Sys.getenv("SCRNA_DATA_DIR", unset = file.path(METHODS_ROOT, "data"))
OUT_DIR <- Sys.getenv("SCRNA_OUT_DIR", unset = file.path(METHODS_ROOT, "output"))
DATA_DIR <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)
OUT_DIR <- normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE)

PROJECT_DIR <- METHODS_ROOT
CACHE_DIR <- DATA_DIR
R_DIR <- file.path(METHODS_ROOT, "R")
if (!dir.exists(R_DIR)) {
    R_DIR <- METHODS_R_DIR
}

SEURAT_RDS <- Sys.getenv("SCRNA_SEURAT_RDS", unset = "seurat_object.rds")
CLUSTASSESS_RDS <- Sys.getenv("SCRNA_CLUSTASSESS_RDS", unset = "clustassess_object.rds")

if (!dir.exists(OUT_DIR)) {
    dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

source_methods <- function(script_name) {
    path <- file.path(R_DIR, script_name)
    if (!file.exists(path)) {
        stop("Missing Methods script: ", path)
    }
    source(path, local = FALSE)
    invisible(path)
}

ensure_out_dir <- function(analysis_name) {
    path <- file.path(OUT_DIR, analysis_name)
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    path
}

message("METHODS_ROOT=", METHODS_ROOT)
message("DATA_DIR=", DATA_DIR)
message("OUT_DIR=", OUT_DIR)
