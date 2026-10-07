#!/usr/bin/env Rscript
# Paired-end DADA2 execution; parameters must be committed before real analyses.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
source(file.path(dirname(script_path), "helpers.R"))

usage <- function() {
  cat(paste0(
    "Usage: Rscript scripts/dada2/run_pilot.R --config PATH --manifest PATH ",
    "--output-dir PATH --intermediate-dir PATH [--input-dir PATH] [--validate-only]\n",
    "Optional: --raw-quality PATH --trimmed-quality PATH to validate/count preceding QC.\n",
    "Optional: --threads INTEGER limits runtime threads to 1..configured threads.\n",
    "Inputs: configuration and one-study paired Illumina manifest; primer-trimmed FASTQ.\n",
    "Paths resolve from the repository root. Default input: data/interim/{study_id}/pilot.\n",
    "Outputs: read_tracking.tsv, ASV counts/sequences/FASTA, length and merge diagnostics,\n",
    "error-model PDF/SVG/PNG (300 dpi), versions, provenance, warnings and success marker.\n",
    "Heavy filtered FASTQ and RDS remain under --intermediate-dir (results/intermediate).\n",
    "Fails on invalid metadata/config, mixed studies, missing/mismatched reads,\n",
    "empty samples after critical stages, invalid count conservation or unusable errors.\n",
    "--validate-only checks inputs/settings without filtering or writing outputs.\n"))
}

parse_cli <- function(args) {
  if ("--help" %in% args || "-h" %in% args) {
    usage()
    quit(status = 0)
  }
  result <- list(validate_only = FALSE)
  allowed <- c("--config", "--manifest", "--output-dir", "--intermediate-dir", "--input-dir",
               "--raw-quality", "--trimmed-quality", "--threads")
  index <- 1L
  while (index <= length(args)) {
    key <- args[[index]]
    if (key == "--validate-only") {
      result$validate_only <- TRUE
      index <- index + 1L
    } else {
      if (!key %in% allowed || index == length(args) || startsWith(args[[index + 1L]], "--")) {
        fail("Unknown argument or missing value: ", key)
      }
      name <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
      if (!is.null(result[[name]])) fail("Duplicate CLI argument: ", key)
      result[[name]] <- args[[index + 1L]]
      index <- index + 2L
    }
  }
  require_fields(result, c("config", "manifest", "output_dir", "intermediate_dir"), "CLI")
  if (xor(is.null(result$raw_quality), is.null(result$trimmed_quality))) {
    fail("Provide both --raw-quality and --trimmed-quality, or neither")
  }
  result
}

resolve_path <- function(path) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(root, path)
}

git_output <- function(args) {
  output <- system2("git", c("-C", shQuote(root), args), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(output, "status"))) fail("Unable to read Git provenance: ", paste(output, collapse = "\n"))
  output
}

export_error_plot <- function(error_model, direction, output_dir) {
  plot <- dada2::plotErrors(error_model, nominalQ = TRUE)
  for (extension in c("pdf", "svg", "png")) {
    path <- file.path(output_dir, paste0("error_model_", direction, ".", extension))
    if (extension == "pdf") grDevices::pdf(path, width = 10, height = 8)
    if (extension == "svg") grDevices::svg(path, width = 10, height = 8)
    if (extension == "png") grDevices::png(path, width = 10, height = 8, units = "in", res = 300, type = "cairo")
    tryCatch(print(plot), finally = grDevices::dev.off())
  }
}

check_error_model <- function(model, direction) {
  errors <- dada2::getErrors(model)
  if (!is.matrix(errors) || nrow(errors) != 16 || anyNA(errors) ||
      any(!is.finite(errors) | errors <= 0 | errors > 1)) fail("Unusable learned error rates: ", direction)
}

run_pipeline <- function(args, config, manifest, settings, paths, warning_log) {
  output_dir <- resolve_path(args$output_dir)
  intermediate_dir <- resolve_path(args$intermediate_dir)
  if (!startsWith(normalizePath(intermediate_dir, mustWork = FALSE, winslash = "/"),
                  paste0(normalizePath(root, winslash = "/"), "/results/intermediate/"))) {
    fail("Heavy outputs must be under results/intermediate/")
  }
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  revision <- git_output(c("rev-parse", "HEAD"))
  status <- git_output(c("status", "--porcelain"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  success_path <- file.path(output_dir, "SUCCESS")
  if (file.exists(success_path)) unlink(success_path)
  validation_path <- file.path(output_dir, "validation.json")
  if (file.exists(validation_path)) unlink(validation_path)
  input_files <- c(resolve_path(args$config), resolve_path(args$manifest), script_path,
                   file.path(dirname(script_path), "helpers.R"), unlist(paths, use.names = FALSE))
  if (!is.null(args$raw_quality)) {
    input_files <- c(input_files, resolve_path(args$raw_quality), resolve_path(args$trimmed_quality))
  }
  hashes <- data.frame(path = vapply(input_files, relative_path, character(1), root = root),
                       md5 = unname(tools::md5sum(input_files)), stringsAsFactors = FALSE)
  write_tsv(hashes, file.path(output_dir, "input_checksums.tsv"))
  writeLines(yaml::as.yaml(config), file.path(output_dir, "config_snapshot.yaml"))
  write_tsv(manifest, file.path(output_dir, "sample_metadata.tsv"))
  packages <- c("dada2", "yaml", "ShortRead", "Biostrings", "RcppParallel", "ggplot2")
  versions <- data.frame(software = c("R", packages),
                          version = c(as.character(getRversion()), vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1))))
  write_tsv(versions, file.path(output_dir, "software_versions.tsv"))
  provenance <- list(status = "running", started_at_utc = started, git_commit = revision,
                     git_dirty = length(status) > 0L, git_status = status,
                     command = c("Rscript", relative_path(script_path, root), commandArgs(trailingOnly = TRUE)),
                     study_id = unique(manifest$study_id), bioproject = unique(manifest$bioproject),
                     manifest_sample_count = nrow(manifest), random_seed = settings$random_seed,
                     settings = settings, checksum_algorithm = "MD5", checksum_table = "input_checksums.tsv",
                     software_versions = "software_versions.tsv", host_platform = R.version$platform,
                     effective_threads = args$effective_threads)
  on.exit({
    provenance$finished_at_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    provenance$warnings <- warning_log$messages
    if (!is.null(warning_log$error)) provenance$failure_message <- warning_log$error
    yaml::write_yaml(provenance, file.path(output_dir, "run_provenance.yaml"))
    writeLines(warning_log$messages, file.path(output_dir, "warnings.txt"))
    writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
    if (identical(provenance$status, "complete")) writeLines(revision, success_path)
  }, add = TRUE)
  provenance$status <- "failed"
  set.seed(settings$random_seed)
  filtered_dir <- file.path(intermediate_dir, "filtered")
  dir.create(filtered_dir, showWarnings = FALSE)
  filtered <- lapply(c("F", "R"), function(direction) {
    stats::setNames(file.path(filtered_dir, paste0(manifest$run_accession, "_", direction, ".fastq.gz")), manifest$sample_id)
  })
  filt <- settings$filter
  threads <- if (args$effective_threads == 1) FALSE else as.integer(args$effective_threads)
  log_step("Filtering ", nrow(manifest), " paired samples")
  filtered_counts <- dada2::filterAndTrim(paths[[1]], filtered[[1]], paths[[2]], filtered[[2]],
    truncLen = filt$trunc_len, maxN = filt$max_n, maxEE = filt$max_ee, truncQ = filt$trunc_q,
    minLen = filt$min_len, trimLeft = filt$trim_left, rm.phix = filt$rm_phix,
    matchIDs = filt$match_ids, id.field = filt$id_field, id.sep = filt$id_sep,
    compress = TRUE, multithread = threads, verbose = TRUE)
  tracking <- manifest[, c("study_id", "bioproject", "sample_id", "run_accession", "fermentation_batch",
                           "fermentation_hours", "relative_time", "fermentation_stage"), drop = FALSE]
  tracking$input <- filtered_counts[, "reads.in"]
  tracking$filtered <- filtered_counts[, "reads.out"]
  if (!is.null(args$raw_quality)) {
    raw_counts <- quality_counts(resolve_path(args$raw_quality), manifest)
    trimmed_counts <- quality_counts(resolve_path(args$trimmed_quality), manifest)
    if (any(trimmed_counts != tracking$input) || any(raw_counts < trimmed_counts)) {
      fail("FastQC and DADA2 input read counts disagree")
    }
    tracking$raw_input <- raw_counts
  }
  write_tsv(tracking, file.path(output_dir, "read_tracking.tsv"))
  if (any(tracking$filtered == 0)) fail("Sample(s) with zero filtered pairs: ", paste(tracking$sample_id[tracking$filtered == 0], collapse = ", "))
  error_settings <- settings$error_learning
  error_diagnostics <- list()
  models <- lapply(seq_along(filtered), function(index) {
    log_step("Learning error rates for direction ", index)
    value <- dada2::learnErrors(filtered[[index]], nbases = error_settings$nbases,
      randomize = error_settings$randomize, MAX_CONSIST = error_settings$max_consist,
      OMEGA_C = error_settings$omega_c, multithread = threads, verbose = TRUE)
    check_error_model(value, c("F", "R")[[index]])
    diagnostic <- error_convergence(value, c("F", "R")[[index]], error_settings$max_consist)
    error_diagnostics[[index]] <<- diagnostic
    write_tsv(do.call(rbind, error_diagnostics), file.path(output_dir, "error_learning.tsv"))
    export_error_plot(value, c("F", "R")[[index]], output_dir)
    if (!diagnostic$converged_before_limit) {
      warning("Error-learning convergence requires review for direction ", c("F", "R")[[index]],
               ": rounds=", diagnostic$self_consistency_rounds,
               ", repeated matrix=", diagnostic$repeated_error_matrix, call. = FALSE)
    }
    value
  })
  provenance$error_models_converged_before_limit <- all(vapply(error_diagnostics, function(value) value$converged_before_limit, logical(1)))
  saveRDS(models, file.path(intermediate_dir, "error_models.rds"))
  log_step("Inferring samples independently")
  dereplicated <- lapply(filtered, dada2::derepFastq, verbose = TRUE)
  # derepFastq returns a single object for a one-file input; retain list shape.
  dereplicated <- lapply(dereplicated, function(value) if (inherits(value, "derep")) list(value) else value)
  dereplicated <- lapply(dereplicated, function(value) stats::setNames(value, manifest$sample_id))
  inferred <- lapply(seq_along(dereplicated), function(index) {
    dada2::dada(dereplicated[[index]], err = models[[index]], pool = settings$inference$pool,
                OMEGA_A = settings$inference$omega_a, OMEGA_C = settings$inference$omega_c,
                multithread = threads, verbose = TRUE)
  })
  inferred <- lapply(inferred, function(value) if (inherits(value, "dada")) list(value) else value)
  inferred <- lapply(inferred, function(value) stats::setNames(value, manifest$sample_id))
  count_reads <- function(value) sum(dada2::getUniques(value))
  tracking$denoised_f <- vapply(inferred[[1]], count_reads, numeric(1))
  tracking$denoised_r <- vapply(inferred[[2]], count_reads, numeric(1))
  log_step("Merging denoised read pairs")
  merge <- settings$merge
  mergers_all <- dada2::mergePairs(inferred[[1]], dereplicated[[1]], inferred[[2]], dereplicated[[2]],
    minOverlap = merge$min_overlap, maxMismatch = merge$max_mismatch,
    trimOverhang = merge$trim_overhang, justConcatenate = merge$just_concatenate,
    returnRejects = TRUE, verbose = TRUE)
  if (is.data.frame(mergers_all)) mergers_all <- stats::setNames(list(mergers_all), manifest$sample_id)
  mergers <- lapply(mergers_all, function(value) value[value$accept, , drop = FALSE])
  merge_diagnostics <- data.frame(study_id = manifest$study_id, sample_id = manifest$sample_id,
    accepted_pairs = vapply(mergers_all, function(value) sum(value$abundance[value$accept]), numeric(1)),
    rejected_pairs = vapply(mergers_all, function(value) sum(value$abundance[!value$accept]), numeric(1)))
  write_tsv(merge_diagnostics, file.path(output_dir, "merge_diagnostics.tsv"))
  tracking$merged <- merge_diagnostics$accepted_pairs
  write_tsv(tracking, file.path(output_dir, "read_tracking.tsv"))
  if (any(tracking$merged == 0)) fail("Sample(s) with zero merged pairs: ", paste(tracking$sample_id[tracking$merged == 0], collapse = ", "))
  sequence_table <- dada2::makeSequenceTable(mergers)
  log_step("Removing chimeras by within-study consensus")
  chim <- settings$chimera
  nonchim <- dada2::removeBimeraDenovo(sequence_table, method = chim$method,
    minFoldParentOverAbundance = chim$min_fold_parent_over_abundance,
    minParentAbundance = chim$min_parent_abundance, minSampleFraction = chim$min_sample_fraction,
    ignoreNNegatives = chim$ignore_n_negatives, allowOneOff = chim$allow_one_off,
    multithread = threads, verbose = TRUE)
  tracking$nonchim <- rowSums(nonchim)[tracking$sample_id]
  validate_tracking(tracking)
  tracking$filtered_fraction <- tracking$filtered / tracking$input
  tracking$merged_fraction_of_filtered <- tracking$merged / tracking$filtered
  tracking$nonchim_fraction_of_merged <- tracking$nonchim / tracking$merged
  tracking$nonchim_fraction_of_input <- tracking$nonchim / tracking$input
  tracking$asv_count <- rowSums(nonchim > 0)[tracking$sample_id]
  write_tsv(tracking, file.path(output_dir, "read_tracking.tsv"))
  if (any(tracking$nonchim == 0)) fail("Sample(s) with zero non-chimeric pairs: ", paste(tracking$sample_id[tracking$nonchim == 0], collapse = ", "))
  study_id <- unique(manifest$study_id)
  exports <- sequence_exports(nonchim, study_id)
  write_tsv(exports$sequences, file.path(output_dir, "asv_sequences.tsv"))
  write_tsv(exports$counts, file.path(output_dir, "asv_counts.tsv"))
  writeLines(as.vector(rbind(paste0(">", exports$sequences$asv_id), exports$sequences$sequence)),
             file.path(output_dir, "asv_sequences.fasta"))
  write_tsv(rbind(length_summary(sequence_table, "merged", study_id), length_summary(nonchim, "nonchim", study_id)),
             file.path(output_dir, "sequence_length_distribution.tsv"))
  summary <- data.frame(study_id = study_id, bioproject = unique(manifest$bioproject), samples = nrow(manifest),
    input_pairs = sum(tracking$input), filtered_pairs = sum(tracking$filtered), merged_pairs = sum(tracking$merged),
    nonchim_pairs = sum(tracking$nonchim), merged_asvs = ncol(sequence_table), nonchim_asvs = ncol(nonchim),
    nonchim_fraction_of_input = sum(tracking$nonchim) / sum(tracking$input))
  write_tsv(summary, file.path(output_dir, "summary.tsv"))
  saveRDS(list(sequence_table = sequence_table, nonchim = nonchim, mergers = mergers_all,
               inferred = inferred), file.path(intermediate_dir, "dada2_objects.rds"))
  if (!identical(unname(tools::md5sum(input_files)), hashes$md5)) {
    fail("An input, configuration, or script changed during this analysis; outputs are not certified")
  }
  provenance$input_checksums_verified_at_completion <- TRUE
  output_files <- list.files(output_dir, full.names = TRUE)
  output_files <- output_files[!file.info(output_files)$isdir]
  output_files <- output_files[!basename(output_files) %in%
    c("run_provenance.yaml", "output_checksums.tsv", "warnings.txt", "session_info.txt", "SUCCESS", "validation.json")]
  write_tsv(data.frame(path = basename(output_files), md5 = unname(tools::md5sum(output_files))),
             file.path(output_dir, "output_checksums.tsv"))
  provenance$output_checksums <- "output_checksums.tsv"
  provenance$status <- "complete"
  log_step("Complete: ", sum(tracking$nonchim), " non-chimeric pairs; ", ncol(nonchim), " ASVs")
}

main <- function() {
  options(warn = 1)
  args <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("yaml", quietly = TRUE)) fail("The declared r-yaml dependency is unavailable")
  config <- yaml::read_yaml(resolve_path(args$config))
  manifest <- utils::read.delim(resolve_path(args$manifest), colClasses = "character", check.names = FALSE,
                               quote = "", comment.char = "", na.strings = character())
  validate_manifest(manifest)
  study_config <- config$amplicon_processing[[unique(manifest$bioproject)]]
  if (is.null(study_config) || !identical(study_config$study_id, unique(manifest$study_id))) fail("Manifest study is not configured")
  if (!identical(study_config$marker, "16S rRNA")) fail("Only configured 16S rRNA studies are supported")
  settings <- normalize_settings(study_config$dada2)
  validate_settings(settings, config$project$default_random_seed)
  if (!is.null(args$threads) && !grepl("^[0-9]+$", args$threads)) fail("--threads must be an integer")
  requested_threads <- if (is.null(args$threads)) settings$threads else as.numeric(args$threads)
  check_number(requested_threads, "runtime threads", minimum = 1, maximum = settings$threads, integer = TRUE)
  args$effective_threads <- if (.Platform$OS.type == "windows") 1L else as.integer(requested_threads)
  if (is.null(args$input_dir)) args$input_dir <- file.path("data", "interim", unique(manifest$study_id), "pilot")
  paths <- input_paths(manifest, resolve_path(args$input_dir))
  if (!is.null(args$raw_quality)) {
    raw_counts <- quality_counts(resolve_path(args$raw_quality), manifest)
    trimmed_counts <- quality_counts(resolve_path(args$trimmed_quality), manifest)
    if (any(raw_counts < trimmed_counts)) fail("Trimmed QC counts exceed raw QC counts")
  }
  if (args$validate_only) {
    log_step("Validated ", nrow(manifest), " manifest rows, settings and paired input paths")
    return(invisible(NULL))
  }
  if (!requireNamespace("dada2", quietly = TRUE)) fail("The declared DADA2 dependency is unavailable")
  warning_log <- new.env(parent = emptyenv())
  warning_log$messages <- character()
  warning_log$error <- NULL
  withCallingHandlers(run_pipeline(args, config, manifest, settings, paths, warning_log), warning = function(value) {
    warning_log$messages <- c(warning_log$messages, conditionMessage(value))
    # Do not muffle: scientific warnings remain visible and are archived as well.
  }, error = function(value) warning_log$error <- conditionMessage(value))
}

tryCatch(main(), error = function(value) {
  message("ERROR: ", conditionMessage(value))
  quit(status = 1)
})
