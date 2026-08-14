#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: execute_r_notebook.R input.ipynb output.executed.ipynb", call. = FALSE)
input_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_path <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.", call. = FALSE)

notebook <- jsonlite::read_json(input_path, simplifyVector = FALSE)
execution_environment <- new.env(parent = globalenv())
execution_count <- 0L
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

for (index in seq_along(notebook$cells)) {
  cell <- notebook$cells[[index]]
  if (!identical(cell$cell_type, "code")) next
  execution_count <- execution_count + 1L
  source_text <- paste(unlist(cell$source), collapse = "")
  output_text <- character()
  error <- NULL
  output_text <- tryCatch(
    capture.output({
      value <- withVisible(eval(parse(text = source_text), envir = execution_environment))
      if (isTRUE(value$visible)) print(value$value)
    }),
    error = function(condition) {
      error <<- condition
      character()
    }
  )
  notebook$cells[[index]]$execution_count <- execution_count
  if (is.null(error)) {
    notebook$cells[[index]]$outputs <- if (length(output_text)) list(list(name = "stdout", output_type = "stream", text = paste0(output_text, "\n"))) else list()
  } else {
    notebook$cells[[index]]$outputs <- list(list(
      output_type = "error", ename = class(error)[[1]], evalue = conditionMessage(error),
      traceback = conditionCall(error) |> deparse()
    ))
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(notebook, output_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    stop(sprintf("Notebook cell %d failed: %s", execution_count, conditionMessage(error)), call. = FALSE)
  }
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(notebook, output_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("Executed %d R code cells: %s\n", execution_count, output_path))
