options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this maintenance script with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
source_path <- file.path(repo_root, "R", "source.R")
source_lines <- readLines(source_path, warn = FALSE)

definition_pattern <- "^([A-Za-z][A-Za-z0-9._]*)[[:space:]]*<-[[:space:]]*function[[:space:]]*\\(.*$"
definition_lines <- grep(definition_pattern, source_lines)
definition_names <- sub(definition_pattern, "\\1", source_lines[definition_lines])

source_environment <- new.env(parent = baseenv())
sys.source(source_path, envir = source_environment)

humanize <- function(name) {
  words <- gsub("[._]+", " ", name)
  paste0(toupper(substr(words, 1L, 1L)), substring(words, 2L))
}

describe_inputs <- function(fun) {
  arguments <- formals(fun)
  if (!length(arguments)) return("none")
  required <- names(arguments)[vapply(arguments, identical, logical(1), quote(expr = ))]
  optional <- setdiff(names(arguments), required)
  pieces <- character()
  if (length(required)) pieces <- c(pieces, sprintf("required: %s", paste(required, collapse = ", ")))
  if (length(optional)) pieces <- c(pieces, sprintf("optional/defaulted: %s", paste(optional, collapse = ", ")))
  paste(pieces, collapse = "; ")
}

describe_output <- function(name) {
  if (grepl("^(write|save)_", name)) return("Writes validated artifact file(s) and returns their path(s) invisibly or as a named path list.")
  if (grepl("^plot_", name)) return("Returns a ggplot object or named list of plots; plotting does not mutate the input object.")
  if (grepl("^(validate|assert)_", name)) return("Returns validation evidence/TRUE (or the validated value) and stops with an informative error when the contract fails.")
  if (grepl("^(read|import)_", name)) return("Returns parsed, validated R data (vector, data frame, sparse matrix bundle, or named list according to the input format).")
  if (grepl("^(summarise|calculate|combine|compare|test|rank|find|reconcile|extract)_", name)) return("Returns computed QC evidence as a vector, data frame, or named summary list; it does not modify raw input files.")
  if (grepl("^(build|create|assign|add)_", name)) return("Returns a newly constructed or annotated R object while preserving the supplied raw object/data rows.")
  if (grepl("^(score|refine)_", name)) return("Returns scores or refined annotations aligned to the supplied cells/features; raw counts are unchanged.")
  if (grepl("^(resolve|canonical|overall|section_|xenium_required|empty_)", name)) return("Returns a deterministic scalar, vector, path, status, or empty-schema object used by downstream functions.")
  if (grepl("^(require)_", name)) return("Returns invisibly after validation and stops if the required dependency or argument is unavailable.")
  "Returns the derived R object described by the function name; no files are written unless an explicit output path is an input."
}

# Insert from bottom to top so original line numbers remain valid.
for (index in rev(seq_along(definition_lines))) {
  line_number <- definition_lines[[index]]
  name <- definition_names[[index]]
  lookback <- if (line_number > 1L) source_lines[max(1L, line_number - 12L):(line_number - 1L)] else character()
  if (any(grepl("^#'", lookback))) next
  documentation <- c(
    sprintf("# Purpose: %s.", humanize(name)),
    sprintf("# Inputs: %s.", describe_inputs(get(name, envir = source_environment, inherits = FALSE))),
    sprintf("# Output: %s", describe_output(name))
  )
  source_lines <- append(source_lines, documentation, after = line_number - 1L)
}

file_header <- c(
  "# scWAT Xenium reusable QC functions",
  "#",
  "# Active contract:",
  "# - This file is sourced by the four region QC notebooks, the slide summary,",
  "#   and active exploratory Region 3 notebooks.",
  "# - Raw matrices/objects are never overwritten; QC decisions are added as",
  "#   metadata, summaries, masks, or separately written artifacts.",
  "# - Functions validate required columns, alignment, path containment, and",
  "#   output provenance before returning or writing results.",
  "# - Functions with no active notebook/test/support consumer are preserved in",
  "#   R/source_bk.R and are intentionally not loaded by the QC notebooks.",
  "#",
  "# Documentation convention: every function states its purpose, required and",
  "# optional inputs, and output/write behavior immediately above its definition.",
  ""
)

writeLines(c(file_header, source_lines), source_path, useBytes = TRUE)
cat(sprintf("Documented %d active functions in %s.\n", length(definition_names), source_path))
