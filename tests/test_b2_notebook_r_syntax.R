#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
notebook_dir <- file.path(repo_root, "notebooks")
files <- list.files(
  notebook_dir,
  pattern = "^B2_Region[1-4]_(all_QCpass|adipose_only|lymph_node_only)_479[.]ipynb$",
  full.names = TRUE
)
stopifnot(length(files) == 12L)

for (file in files) {
  notebook <- jsonlite::fromJSON(file, simplifyVector = FALSE)
  for (index in seq_along(notebook$cells)) {
    cell <- notebook$cells[[index]]
    if (identical(cell$cell_type, "code")) {
      source <- paste(unlist(cell$source), collapse = "")
      tryCatch(
        parse(text = source),
        error = function(error) stop(
          basename(file), " code cell ", index, ": ", conditionMessage(error),
          call. = FALSE
        )
      )
    }
  }
}

cat("All split B2 notebook R cells parse successfully.\n")
