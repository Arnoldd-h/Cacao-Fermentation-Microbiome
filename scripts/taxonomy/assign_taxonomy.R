#!/usr/bin/env Rscript
# Input: pinned taxonomy config plus independently validated one-study ASV tables.
# Output: original calls, bootstrap support, masked calls, flags and coverage.
# Fails on invalid input/mixed studies, incompatible ranks, changed reference or R errors.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
# Some Rscript builds encode spaces in their internal --file argument as ~+~.
# Preserve a literal path when it exists; decode only an otherwise absent path.
if (!file.exists(script_file)) script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
script_path <- normalizePath(script_file, mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."))
source(file.path(root, "scripts/dada2/helpers.R"))
source(file.path(dirname(script_path), "helpers.R"))

args <- commandArgs(TRUE)
if ("--help" %in% args || "-h" %in% args) {
  cat("Usage: Rscript scripts/taxonomy/assign_taxonomy.R --config PATH --input-dir PATH --output-dir PATH --threads INTEGER\n",
      "Inputs: taxonomy JSON/YAML, pinned ignored FASTA and one-study DADA2 ASV sequences/counts.\n",
      "Outputs: original/support/masked taxonomy, sensitivity, screening, coverage, versions and figures.\n",
      "Fails on missing or inconsistent input, mixed studies, invalid ranks or classifier failure.\n",
      "Run via run_taxonomy.py to validate inputs and record hashes/provenance/success.\n")
  quit(status = 0)
}
if (length(args) != 8L || anyDuplicated(args[seq(1L, 8L, 2L)])) fail("Invalid CLI; use --help")
options <- as.list(args[seq(2L, 8L, 2L)])
names(options) <- sub("^--", "", args[seq(1L, 8L, 2L)])
require_fields(options, c("config", "input-dir", "output-dir", "threads"), "CLI")
config <- yaml::read_yaml(options$config)
settings <- config$classification
threads <- as.integer(options$threads)
check_number(threads, "threads", minimum = 1, maximum = settings$threads, integer = TRUE)
sequences <- utils::read.delim(file.path(options[["input-dir"]], "asv_sequences.tsv"), colClasses = "character", check.names = FALSE)
counts_table <- utils::read.delim(file.path(options[["input-dir"]], "asv_counts.tsv"), check.names = FALSE)
require_fields(sequences, c("study_id", "asv_id", "sequence", "total_reads"), "sequences")
if (!nrow(sequences) || anyDuplicated(sequences$asv_id) || anyDuplicated(sequences$sequence) ||
    any(!grepl("^[ACGT]+$", sequences$sequence)) || length(unique(sequences$study_id)) != 1L ||
    any(counts_table$study_id != sequences$study_id[1])) fail("Invalid or mixed-study ASVs")
if (!setequal(colnames(counts_table)[-(1:2)], sequences$asv_id)) fail("Counts and sequences disagree")
counts <- as.matrix(counts_table[, sequences$asv_id, drop = FALSE])
rownames(counts) <- counts_table$sample_id
if (anyNA(counts) || any(counts < 0 | counts != floor(counts)) ||
    any(colSums(counts) != as.numeric(sequences$total_reads))) fail("Invalid ASV counts")
reference <- file.path(root, config$reference$path)
if (file.info(reference)$size != config$reference$expected_bytes ||
    unname(tools::md5sum(reference)) != config$reference$expected_md5) fail("Reference checksum mismatch")

dir.create(options[["output-dir"]], recursive = TRUE, showWarnings = FALSE)
warnings <- character()
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(settings$random_seed)
log_step("Classifying ", nrow(sequences), " ASVs against ", config$reference$name, " ", config$reference$version)
assignment <- withCallingHandlers(dada2::assignTaxonomy(sequences$sequence, reference, minBoot = 0,
  tryRC = settings$try_reverse_complement, outputBootstraps = TRUE,
  taxLevels = settings$tax_levels, multithread = threads, verbose = TRUE),
  warning = function(w) { warnings <<- c(warnings, conditionMessage(w)) })
if (!identical(colnames(assignment$tax), settings$tax_levels) ||
    !identical(rownames(assignment$tax), sequences$sequence)) fail("Classifier order or ranks disagree")
identifiers <- sequences[, c("study_id", "asv_id"), drop = FALSE]
output <- options[["output-dir"]]
export <- function(value, name) write_tsv(value, file.path(output, name))
export(cbind(identifiers, assignment$tax), "taxonomy_unfiltered.tsv")
export(cbind(identifiers, assignment$boot), "taxonomy_bootstraps.tsv")
thresholds <- c(settings$primary_min_boot, settings$sensitivity_min_boot)
masked <- lapply(thresholds, function(threshold) mask_taxonomy(assignment$tax, assignment$boot, threshold))
export(cbind(identifiers, masked[[1]]), "taxonomy.tsv")
export(do.call(rbind, lapply(seq_along(thresholds), function(i)
  cbind(identifiers, min_boot = thresholds[i], masked[[i]]))), "taxonomy_sensitivity.tsv")
screen <- screen_taxonomy(masked[[1]], config$screening)
export(cbind(identifiers, total_reads = sequences$total_reads, screen), "taxonomy_screening.tsv")
coverage <- do.call(rbind, lapply(seq_along(thresholds), function(i)
  taxonomy_coverage(masked[[i]], counts, sequences$study_id[1], thresholds[i])))
export(coverage, "assignment_coverage.tsv")
export(do.call(rbind, lapply(seq_along(thresholds), function(i)
  sample_coverage(masked[[i]], counts, sequences$study_id[1], thresholds[i]))), "sample_coverage.tsv")
packages <- c("dada2", "yaml", "Biostrings", "RcppParallel")
export(data.frame(software = c("R", packages), version = c(as.character(getRversion()),
  vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)))), "software_versions.tsv")
writeLines(capture.output(utils::sessionInfo()), file.path(output, "session_info.txt"))
writeLines(warnings, file.path(output, "warnings.txt"))

for (extension in c("pdf", "svg", "png")) {
  path <- file.path(output, paste0("assignment_coverage.", extension))
  if (extension == "pdf") grDevices::pdf(path, width = 10, height = 5)
  if (extension == "svg") grDevices::svg(path, width = 10, height = 5)
  if (extension == "png") grDevices::png(path, width = 10, height = 5, units = "in", res = 300, type = "cairo")
  graphics::par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
  for (metric in c("asvs", "reads")) {
    values <- t(vapply(seq_along(thresholds), function(i) {
      rows <- coverage[coverage$min_boot == thresholds[i], ]
      if (metric == "asvs") 100 * rows$assigned_asvs / rows$total_asvs else 100 * rows$assigned_reads / rows$total_reads
    }, numeric(length(settings$tax_levels))))
    graphics::barplot(values, beside = TRUE, names.arg = settings$tax_levels, las = 2,
      ylim = c(0, 110), ylab = "% assigned", main = paste("Pilot coverage:", metric),
      col = c("#27647B", "#BC7D38"), legend.text = paste("Bootstrap >=", thresholds),
      args.legend = list(bty = "n", cex = 0.8))
  }
  grDevices::dev.off()
}
log_step("Taxonomy tables and coverage exported; no ASVs excluded")
