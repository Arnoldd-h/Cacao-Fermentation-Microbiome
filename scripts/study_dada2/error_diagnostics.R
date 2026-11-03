#!/usr/bin/env Rscript
# Read-only diagnostics of saved DADA2 error models and execution settings.

diagnostic_fail <- function(...) stop(..., call. = FALSE)

diagnostic_usage <- function() {
  cat(paste0(
    "Usage: Rscript scripts/study_dada2/error_diagnostics.R --config PATH --manifest PATH ",
    "--run-dir PATH --intermediate-dir PATH --output PATH --expected-r VERSION ",
    "--expected-dada2 VERSION --expected-threads INTEGER\n",
    "Inputs: committed processing profile, manifest, completed DADA2 run and error_models.rds.\n",
    "Output: one TSV with F/R matrix-repeat diagnostics; no model or sequence edits.\n",
    "Relative paths resolve from the repository root. The .tsv output must be outside\n",
    "the run/intermediate directories and must not overwrite an input.\n",
    "Fails on inconsistent settings, versions, seed, sample count, threads or malformed models.\n",
    "Repeated indices are 1-based err_in entries; repetition need not mean a fixed point.\n",
    "DADA2 learnErrors does not store training sample/read/base totals in its model RDS;\n",
    "those columns remain empty rather than being inferred from transition counts.\n"))
}

parse_diagnostic_cli <- function(args) {
  required <- c("config", "manifest", "run-dir", "intermediate-dir", "output",
                "expected-r", "expected-dada2", "expected-threads")
  if (!length(args)) diagnostic_fail("All documented arguments are required")
  if (length(args) %% 2L) diagnostic_fail("Each argument requires a value")
  result <- list()
  for (index in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[[index]])
    if (!startsWith(args[[index]], "--") || !key %in% required || !is.null(result[[key]]) ||
        !nzchar(args[[index + 1L]]) || startsWith(args[[index + 1L]], "--")) {
      diagnostic_fail("Unknown, duplicate or empty argument: ", args[[index]])
    }
    result[[key]] <- args[[index + 1L]]
  }
  if (!setequal(names(result), required)) diagnostic_fail("All documented arguments are required")
  if (!grepl("^[1-9][0-9]*$", result[["expected-threads"]])) diagnostic_fail("Invalid expected threads")
  result[["expected-threads"]] <- as.numeric(result[["expected-threads"]])
  if (!is.finite(result[["expected-threads"]])) diagnostic_fail("Invalid expected threads")
  if (!grepl("[.]tsv$", result$output)) diagnostic_fail("Output must end in .tsv")
  result
}

check_diagnostic_matrix <- function(value, label, counts = FALSE) {
  if (!is.matrix(value) || !is.numeric(value) || nrow(value) != 16L || ncol(value) < 1L ||
      anyNA(value) || any(!is.finite(value)) || any(value < 0) ||
      (!counts && any(value > 1)) || (counts && any(value != floor(value)))) {
    diagnostic_fail("Invalid ", label, " matrix")
  }
}

summarize_error_model <- function(model, direction) {
  if (!is.list(model) || !is.list(model$err_in) || !length(model$err_in)) {
    diagnostic_fail("Missing error iteration history: ", direction)
  }
  check_diagnostic_matrix(model$err_out, "output error")
  if (any(model$err_out <= 0)) diagnostic_fail("Output error rates must be positive")
  for (value in model$err_in) check_diagnostic_matrix(value, "input error")
  check_diagnostic_matrix(model$trans, "transition", counts = TRUE)
  matches <- which(vapply(model$err_in, identical, logical(1), model$err_out))
  data.frame(read_direction = direction, self_consistency_rounds = length(model$err_in),
    repeated_any = length(matches) > 0L,
    last_input_equals_output = identical(model$err_in[[length(model$err_in)]], model$err_out),
    repeated_iteration_indices = paste(matches, collapse = ";"),
    transition_base_observations = sum(model$trans),
    training_samples = "", training_reads = "", training_bases = "",
    training_metadata_source = "not_stored_by_learnErrors", stringsAsFactors = FALSE)
}

summarize_error_models <- function(models) {
  if (!is.list(models) || length(models) != 2L) diagnostic_fail("Expected exactly two F/R error models")
  if (!is.null(names(models))) {
    if (anyDuplicated(names(models)) || !setequal(names(models), c("F", "R"))) {
      diagnostic_fail("Ambiguous error model direction names")
    }
    models <- models[c("F", "R")]
  }
  # run_pilot.R saves an unnamed list in forward/reverse order.
  do.call(rbind, Map(summarize_error_model, models, c("F", "R")))
}

check_dada2_execution <- function(config_path, manifest_path, run_dir, expected_r,
                                  expected_dada2, expected_threads, helpers_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) diagnostic_fail("The registered r-yaml dependency is unavailable")
  helpers <- new.env(parent = baseenv())
  sys.source(helpers_path, envir = helpers)
  config <- yaml::read_yaml(config_path)
  snapshot <- yaml::read_yaml(file.path(run_dir, "config_snapshot.yaml"))
  if (!identical(config, snapshot)) diagnostic_fail("Configuration snapshot differs from the registered profile")
  manifest <- utils::read.delim(manifest_path, colClasses = "character", check.names = FALSE,
                              quote = "", comment.char = "", na.strings = character())
  helpers$validate_manifest(manifest)
  profile <- config$amplicon_processing[[unique(manifest$bioproject)]]
  if (is.null(profile) || !identical(profile$study_id, unique(manifest$study_id))) {
    diagnostic_fail("Manifest study is not configured")
  }
  settings <- helpers$normalize_settings(profile$dada2)
  helpers$validate_settings(settings, config$project$default_random_seed)
  provenance <- yaml::read_yaml(file.path(run_dir, "run_provenance.yaml"))
  if (!identical(provenance$status, "complete") || !isTRUE(provenance$input_checksums_verified_at_completion)) {
    diagnostic_fail("A completed run with verified inputs is required")
  }
  if (!identical(provenance$settings, settings)) diagnostic_fail("Executed settings differ from the registered profile")
  for (field in c("study_id", "bioproject")) {
    if (!identical(provenance[[field]], unique(manifest[[field]]))) diagnostic_fail("Provenance ", field, " mismatch")
  }
  if (!identical(provenance$manifest_sample_count, nrow(manifest))) diagnostic_fail("Provenance sample count mismatch")
  if (!identical(as.numeric(provenance$random_seed), as.numeric(settings$random_seed))) diagnostic_fail("Provenance seed mismatch")
  helpers$check_number(expected_threads, "expected threads", minimum = 1, maximum = settings$threads, integer = TRUE)
  if (!identical(as.numeric(provenance$effective_threads), as.numeric(expected_threads))) diagnostic_fail("Effective threads mismatch")
  versions <- utils::read.delim(file.path(run_dir, "software_versions.tsv"), colClasses = "character")
  if (!all(c("software", "version") %in% names(versions)) || anyDuplicated(versions$software)) diagnostic_fail("Invalid software version table")
  for (software in c("R", "dada2")) {
    expected <- if (software == "R") expected_r else expected_dada2
    observed <- versions$version[versions$software == software]
    if (!identical(observed, expected)) diagnostic_fail("Recorded ", software, " version mismatch")
  }
  list(settings = settings, provenance = provenance)
}

run_error_diagnostics <- function(args, root) {
  resolve <- function(path) if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(root, path)
  for (key in c("config", "manifest", "run-dir", "intermediate-dir", "output")) args[[key]] <- resolve(args[[key]])
  execution <- check_dada2_execution(args$config, args$manifest, args[["run-dir"]], args[["expected-r"]],
    args[["expected-dada2"]], args[["expected-threads"]], file.path(root, "scripts/dada2/helpers.R"))
  command <- execution$provenance$command
  for (flag in c("--config", "--manifest", "--output-dir", "--intermediate-dir")) {
    index <- which(command == flag)
    if (length(index) != 1L || index >= length(command)) diagnostic_fail("Missing or repeated run command flag: ", flag)
    key <- if (flag == "--output-dir") "run-dir" else substring(flag, 3L)
    observed <- normalizePath(resolve(command[[index + 1L]]), mustWork = TRUE)
    expected <- normalizePath(args[[key]], mustWork = TRUE)
    if (!identical(observed, expected)) diagnostic_fail("Run command path mismatch: ", flag)
  }
  output_path <- normalizePath(args$output, winslash = "/", mustWork = FALSE)
  if (output_path %in% vapply(args[c("config", "manifest")], normalizePath, character(1),
                             winslash = "/", mustWork = TRUE)) diagnostic_fail("Diagnostic output must not overwrite an input")
  for (key in c("run-dir", "intermediate-dir")) {
    directory <- normalizePath(args[[key]], winslash = "/", mustWork = TRUE)
    if (startsWith(output_path, paste0(directory, "/"))) diagnostic_fail("Diagnostic TSV must remain outside the run/intermediate directories")
  }
  rows <- summarize_error_models(readRDS(file.path(args[["intermediate-dir"]], "error_models.rds")))
  limit <- execution$settings$error_learning$max_consist
  if (any(rows$self_consistency_rounds > limit)) diagnostic_fail("Error history exceeds the registered iteration limit")
  rows$maximum_rounds <- limit
  rows$reached_iteration_limit <- rows$self_consistency_rounds >= limit
  rows$converged_before_limit <- rows$repeated_any & !rows$reached_iteration_limit
  if (!identical(all(rows$converged_before_limit), execution$provenance$error_models_converged_before_limit)) {
    diagnostic_fail("Saved models disagree with the recorded convergence indicator")
  }
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = "error_diagnostics_", tmpdir = dirname(args$output))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.table(rows, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (!file.rename(temporary, args$output)) diagnostic_fail("Could not atomically write diagnostic TSV")
  message("Validated execution and wrote F/R diagnostics: ", args$output)
  invisible(rows)
}

if (sys.nframe() == 0L) {
  tryCatch({
    arguments <- commandArgs(trailingOnly = TRUE)
    if (any(arguments %in% c("--help", "-h"))) {
      diagnostic_usage()
    } else {
      script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
      root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
      run_error_diagnostics(parse_diagnostic_cli(arguments), root)
    }
  }, error = function(value) {
    message("ERROR: ", conditionMessage(value))
    quit(status = 1L)
  })
}
