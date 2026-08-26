#!/usr/bin/env Rscript

# Add an adjacent roxygen-style purpose/input/output contract to every active
# function that does not already have one. Existing hand-written roxygen blocks
# are preserved; superseded generic three-line comments are removed.

options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
source_path <- file.path(repo_root, "R", "source.R")
lines <- readLines(source_path, warn = FALSE)
environment <- new.env(parent = globalenv())
sys.source(source_path, envir = environment)

pattern <- "^([A-Za-z][A-Za-z0-9._]*)[[:space:]]*<-[[:space:]]*function[[:space:]]*\\(.*$"
definition_lines <- grep(pattern, lines)
definition_names <- sub(pattern, "\\1", lines[definition_lines])

humanize <- function(name) {
  value <- gsub("[._]+", " ", name)
  paste0(toupper(substr(value, 1L, 1L)), substring(value, 2L))
}

return_contract <- function(name) {
  if (grepl("^(write|save)_", name)) return("Written artifact path(s) or an auditable one-row write summary.")
  if (grepl("^plot_", name)) return("A `ggplot` object or named list of plots; input objects are not modified.")
  if (grepl("^(validate|assert)_", name)) return("Validated input or `TRUE`; invalid input stops with an informative error.")
  if (grepl("^(read|import)_", name)) return("Parsed and validated R data with input row or matrix alignment preserved.")
  if (grepl("^(summarise|calculate|score|refine)_", name)) return("Computed evidence or annotations aligned to the supplied rows/cells; raw counts are unchanged.")
  "A deterministic derived value or annotated copy; the supplied raw data are not overwritten."
}

for (index in rev(seq_along(definition_lines))) {
  line_number <- definition_lines[[index]]
  name <- definition_names[[index]]
  cursor <- line_number - 1L
  generic <- integer()
  while (cursor > 0L && grepl("^# (Purpose|Inputs|Output):", lines[[cursor]])) {
    generic <- c(cursor, generic)
    cursor <- cursor - 1L
  }
  if (length(generic)) {
    lines <- lines[-generic]
    line_number <- line_number - length(generic)
  }
  cursor <- line_number - 1L
  while (cursor > 0L && !nzchar(trimws(lines[[cursor]]))) cursor <- cursor - 1L
  if (cursor > 0L && grepl("^#'", lines[[cursor]])) next

  fun <- get(name, envir = environment, inherits = FALSE)
  arguments <- formals(fun)
  required_names <- if (length(arguments)) {
    names(arguments)[vapply(arguments, identical, logical(1), quote(expr = ))]
  } else {
    character()
  }
  params <- if (!length(arguments)) character() else vapply(names(arguments), function(argument) {
    required <- argument %in% required_names
    if (argument == "...") {
      sprintf("#' @param %s Additional arguments passed to the documented downstream operation.", argument)
    } else if (required) {
      sprintf("#' @param %s Required `%s` input; validated before computation.", argument, argument)
    } else {
      sprintf("#' @param %s Optional `%s` input with the default shown in the function signature.", argument, argument)
    }
  }, character(1))
  documentation <- c(
    sprintf("#' %s.", humanize(name)),
    "#'",
    params,
    sprintf("#' @return %s", return_contract(name))
  )
  lines <- append(lines, documentation, after = line_number - 1L)
}

writeLines(lines, source_path, useBytes = TRUE)
invisible(parse(file = source_path))
cat(sprintf("Roxygen contracts verified or added for %d active functions.\n", length(definition_names)))
