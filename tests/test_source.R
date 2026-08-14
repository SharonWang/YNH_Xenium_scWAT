options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this test with Rscript.", call. = FALSE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))

panel_path <- file.path(repo_root, "config", "eos_gene_sets.tsv")
stopifnot(file.exists(panel_path))
panel <- utils::read.delim(panel_path, check.names = FALSE)
stopifnot(identical(names(panel), c("gene", "gene_set")))
stopifnot(nrow(panel) == 100L, length(unique(panel$gene)) == 100L)
stopifnot(!anyNA(panel), !any(trimws(as.matrix(panel)) == ""))
group_sizes <- table(panel$gene_set)
stopifnot(
  identical(as.integer(group_sizes[c("common", "short_lived", "long_lived")]), c(7L, 47L, 46L))
)

source_path <- file.path(repo_root, "R", "source.R")
stopifnot(file.exists(source_path))

cat("Repository contracts passed.\n")
