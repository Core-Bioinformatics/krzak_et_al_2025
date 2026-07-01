Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  R_DEFAULT_NUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(bulkAnalyseR)
})

counts <- read.delim(
  "bulkAnalyseR_app/NatImm_bulkAnalyseR_input/expression_counts_raw.tsv",
  check.names = FALSE
)

expr <- as.matrix(counts[, -1, drop = FALSE])
rownames(expr) <- counts$Geneid
storage.mode(expr) <- "numeric"

expr <- expr[rowSums(expr) > 0, , drop = FALSE]

cat("Input dim after zero-sum removal:", dim(expr), "\n")

expr.proc <- preprocessExpressionMatrix(expr, output.plot = FALSE)

cat("Processed dim:", dim(expr.proc), "\n")
cat("NA:", sum(is.na(expr.proc)), "\n")
cat("min:", min(expr.proc), "\n")
cat("max:", max(expr.proc), "\n")

save(expr.proc, file = "bulkAnalyseR_app/test_preprocess_standard_expr_proc.rda")
cat("STANDARD_PREPROCESS_OK\n")
