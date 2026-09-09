#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
source(file.path(repo_root, "R", "source.R"))

# Break caught: mclust::Mclust() evaluates an unqualified mclustBIC call in
# its caller, so a namespace-only invocation fails unless the exported BIC
# function is deliberately bound in that caller frame.
fake_mclust_bic <- function(data, G = NULL, verbose = FALSE, ...) {
  structure(matrix(c(-10, -8), nrow = 2L), G = G, modelNames = "V")
}
fake_mclust <- function(data, G = NULL, verbose = FALSE, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- as.name("mclustBIC")
  mc[[2L]] <- data
  bic <- eval(mc, parent.frame())
  list(G = 2L, modelName = "V", BIC = bic)
}
core_fit <- run_mclust_with_binding(
  data = seq(-1, 1, length.out = 40L),
  G = 1:2,
  mclust_fun = fake_mclust,
  mclust_bic_fun = fake_mclust_bic
)
stopifnot(core_fit$G == 2L, core_fit$modelName == "V")

diagnostic <- run_mclust_diagnostic(seq(-1, 1, length.out = 40L), G = 1:3)
if (requireNamespace("mclust", quietly = TRUE)) {
  stopifnot(
    diagnostic$status == "PASS",
    diagnostic$n_finite == 40L,
    diagnostic$selected_G %in% 1:3,
    is.data.frame(diagnostic$bic_table)
  )
} else {
  stopifnot(diagnostic$status == "SKIPPED_PACKAGE_UNAVAILABLE")
}

small <- run_mclust_diagnostic(1:10, min_n = 20L)
if (requireNamespace("mclust", quietly = TRUE)) {
  stopifnot(small$status == "SKIPPED_INSUFFICIENT_DATA")
} else {
  stopifnot(small$status == "SKIPPED_PACKAGE_UNAVAILABLE")
}

cat("Extended Eosinophil helper tests passed.\n")
