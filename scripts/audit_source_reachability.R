#!/usr/bin/env Rscript

# Audit which top-level functions in R/source.R are reachable from committed
# notebooks, current tests, and current execution tooling. The audit is
# read-only and prints ACTIVE/UNUSED names for review before archival.

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
source_path <- file.path(repo_root, "R", "source.R")

environment <- new.env(parent = globalenv())
sys.source(source_path, envir = environment)
function_names <- sort(ls(environment, all.names = TRUE)[vapply(ls(environment, all.names = TRUE), function(name) {
  is.function(get(name, envir = environment, inherits = FALSE))
}, logical(1))])

notebook_paths <- list.files(file.path(repo_root, "notebooks"), pattern = "[.]ipynb$", full.names = TRUE)
notebook_text <- unlist(lapply(notebook_paths, function(path) {
  notebook <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  vapply(notebook$cells, function(cell) {
    if (identical(cell$cell_type, "code")) paste(unlist(cell$source), collapse = "") else ""
  }, character(1))
}), use.names = FALSE)

consumer_paths <- c(
  file.path(repo_root, "tests", c(
    "test_fixed_cell_qc.R", "test_scwat_initial_qc_contract.R",
    "test_source_contract.R", "test_source.R", "test_notebook_contracts.py"
  )),
  file.path(repo_root, "scripts", c(
    "execute_local_subset.ps1", "execute_r_notebook.R",
    "build_scwat_initial_qc_notebooks.py", "render_notebooks.py"
  ))
)
consumer_text <- c(notebook_text, unlist(lapply(consumer_paths[file.exists(consumer_paths)], readLines, warn = FALSE)))
consumer_blob <- paste(consumer_text, collapse = "\n")
roots <- function_names[vapply(function_names, function(name) {
  # Count executable calls only. Function names mentioned as character values
  # in source-contract allow/archive lists must not keep obsolete code active.
  grepl(sprintf("(?<![A-Za-z0-9._])%s[[:space:]]*\\(", name), consumer_blob, perl = TRUE)
}, logical(1))]

dependencies <- lapply(function_names, function(name) {
  globals <- codetools::findGlobals(get(name, envir = environment), merge = FALSE)$functions
  intersect(globals, function_names)
})
names(dependencies) <- function_names
active <- unique(roots)
repeat {
  expanded <- unique(c(active, unlist(dependencies[active], use.names = FALSE)))
  if (setequal(expanded, active)) break
  active <- expanded
}
active <- sort(active)
unused <- setdiff(function_names, active)
cat(sprintf("ACTIVE %d\n", length(active)))
cat(paste(active, collapse = "\n"), "\n")
cat(sprintf("UNUSED %d\n", length(unused)))
cat(paste(unused, collapse = "\n"), "\n")
