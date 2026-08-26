#!/usr/bin/env Rscript

# Move functions with no executable consumer into R/source_bk.R. Reachability
# starts from function calls in committed notebooks, current tests, and current
# execution tooling, then follows the complete source.R call graph.

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
source_path <- file.path(repo_root, "R", "source.R")
backup_path <- file.path(repo_root, "R", "source_bk.R")

source_environment <- new.env(parent = globalenv())
sys.source(source_path, envir = source_environment)
function_names <- sort(ls(source_environment)[vapply(ls(source_environment), function(name) {
  is.function(get(name, envir = source_environment, inherits = FALSE))
}, logical(1))])

notebook_paths <- list.files(file.path(repo_root, "notebooks"), pattern = "[.]ipynb$", full.names = TRUE)
notebook_text <- unlist(lapply(notebook_paths, function(path) {
  notebook <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  vapply(notebook$cells, function(cell) {
    if (identical(cell$cell_type, "code")) paste(unlist(cell$source), collapse = "") else ""
  }, character(1))
}), use.names = FALSE)
consumer_paths <- c(
  list.files(file.path(repo_root, "tests"), full.names = TRUE),
  file.path(repo_root, "scripts", c(
    "execute_local_subset.ps1", "execute_r_notebook.R",
    "build_scwat_initial_qc_notebooks.py", "render_notebooks.py"
  ))
)
consumer_text <- c(notebook_text, unlist(lapply(consumer_paths[file.exists(consumer_paths)], readLines, warn = FALSE)))
consumer_blob <- paste(consumer_text, collapse = "\n")
roots <- function_names[vapply(function_names, function(name) {
  grepl(sprintf("(?<![A-Za-z0-9._])%s[[:space:]]*\\(", name), consumer_blob, perl = TRUE)
}, logical(1))]
dependencies <- lapply(function_names, function(name) {
  intersect(codetools::findGlobals(get(name, envir = source_environment), merge = FALSE)$functions, function_names)
})
names(dependencies) <- function_names
active <- unique(roots)
repeat {
  expanded <- unique(c(active, unlist(dependencies[active], use.names = FALSE)))
  if (setequal(expanded, active)) break
  active <- expanded
}
archive_names <- sort(setdiff(function_names, active))
if (!length(archive_names)) {
  cat("No newly unused source functions were found.\n")
  quit(save = "no", status = 0L)
}

source_lines <- readLines(source_path, warn = FALSE)
definition_pattern <- "^([A-Za-z%][A-Za-z0-9._%]*)[[:space:]]*<-[[:space:]]*function[[:space:]]*\\(.*$"
definition_lines <- grep(definition_pattern, source_lines)
definition_names <- sub(definition_pattern, "\\1", source_lines[definition_lines])
definition_ends <- c(definition_lines[-1L] - 1L, length(source_lines))

blocks <- lapply(archive_names, function(name) {
  position <- match(name, definition_names)
  if (is.na(position)) stop(sprintf("Cannot locate function definition: %s", name), call. = FALSE)
  start <- definition_lines[[position]]
  while (start > 1L && (grepl("^[[:space:]]*#", source_lines[[start - 1L]]) || !nzchar(trimws(source_lines[[start - 1L]])))) {
    start <- start - 1L
  }
  list(name = name, start = start, end = definition_ends[[position]])
})

remove <- rep(FALSE, length(source_lines))
for (block in blocks) remove[block$start:block$end] <- TRUE
active_lines <- source_lines[!remove]
while (length(active_lines) && !nzchar(tail(active_lines, 1L))) active_lines <- head(active_lines, -1L)

backup_lines <- readLines(backup_path, warn = FALSE)
existing_definitions <- sub(definition_pattern, "\\1", backup_lines[grep(definition_pattern, backup_lines)])
duplicates <- intersect(archive_names, existing_definitions)
if (length(duplicates)) stop(sprintf("Archive already defines: %s", paste(duplicates, collapse = ", ")), call. = FALSE)
archive_lines <- c(
  backup_lines,
  "",
  "# -----------------------------------------------------------------------------",
  "# Archived 2026-08-27 after Colon-parity initial-QC notebook replacement",
  "# -----------------------------------------------------------------------------",
  "",
  unlist(lapply(blocks, function(block) {
    c(sprintf("# Archived function: %s", block$name), source_lines[block$start:block$end], "")
  }), use.names = FALSE)
)
while (length(archive_lines) && !nzchar(tail(archive_lines, 1L))) archive_lines <- head(archive_lines, -1L)

writeLines(active_lines, source_path, useBytes = TRUE)
writeLines(archive_lines, backup_path, useBytes = TRUE)
invisible(parse(file = source_path))
invisible(parse(file = backup_path))
cat(sprintf("Kept %d active functions; moved %d newly unused functions to source_bk.R.\n", length(active), length(archive_names)))
