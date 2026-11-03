#!/usr/bin/env Rscript
# Synthetic models test diagnostics only; they are not biological results.
test_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]))
root <- normalizePath(file.path(dirname(test_path), "..", ".."))
source(file.path(root, "scripts/study_dada2/error_diagnostics.R"))
checks <- 0L
expect_true <- function(value) {
  stopifnot(isTRUE(value))
  checks <<- checks + 1L
}
expect_error <- function(expression, pattern) {
  value <- tryCatch({ force(expression); NULL }, error = identity)
  if (is.null(value) || !grepl(pattern, conditionMessage(value))) stop("Expected error: ", pattern)
  checks <<- checks + 1L
}

first <- matrix(.1, 16L, 4L)
second <- first * 2
transitions <- matrix(2L, 16L, 4L)
fixed <- list(err_in = list(first, second), err_out = second, trans = transitions)
cycle <- list(err_in = list(first, second), err_out = first, trans = transitions)
unrepeated <- list(err_in = list(first, second), err_out = first * 3, trans = transitions)
rows <- summarize_error_models(list(fixed, cycle))
expect_true(identical(rows$read_direction, c("F", "R")))
expect_true(all(rows$repeated_any))
expect_true(identical(rows$last_input_equals_output, c(TRUE, FALSE)))
expect_true(identical(rows$repeated_iteration_indices, c("2", "1")))
expect_true(all(rows$transition_base_observations == 128))
expect_true(all(rows$training_samples == "" & rows$training_reads == "" & rows$training_bases == ""))
expect_true(identical(summarize_error_models(list(R = cycle, F = fixed))$last_input_equals_output, c(TRUE, FALSE)))
expect_true(!summarize_error_model(unrepeated, "F")$repeated_any)
expect_true(summarize_error_model(list(err_in = list(first, second, first), err_out = first, trans = transitions), "F")$repeated_iteration_indices == "1;3")
expect_error(summarize_error_models(list(fixed)), "exactly two")
expect_error(summarize_error_models(list(A = fixed, B = cycle)), "direction names")
broken <- fixed; broken$err_out[1, 1] <- NA_real_
expect_error(summarize_error_model(broken, "F"), "Invalid output")
broken <- fixed; broken$trans[1, 1] <- .5
expect_error(summarize_error_model(broken, "F"), "Invalid transition")
expect_error(parse_diagnostic_cli(character()), "All documented")
expect_error(parse_diagnostic_cli(c("--config", "a", "--config", "b")), "duplicate")

# Use the registered settings as a schema, but all manifests/models below are fixtures.
config <- yaml::read_yaml(file.path(root, "config/prjna627078_processing.yaml"))
helpers <- new.env(parent = baseenv())
sys.source(file.path(root, "scripts/dada2/helpers.R"), envir = helpers)
study <- config$amplicon_processing[[1L]]$study_id
project <- names(config$amplicon_processing)[[1L]]
settings <- helpers$normalize_settings(config$amplicon_processing[[1L]]$dada2)
manifest <- data.frame(study_id = study, bioproject = project, sample_id = "fixture_sample",
  run_accession = "FIXTURE1", fermentation_stage = "early", fermentation_hours = "0",
  relative_time = "0", fermentation_batch = "fixture_batch", library_layout = "PAIRED",
  sequencing_platform = "ILLUMINA", sampling_stratum = "", stringsAsFactors = FALSE)
temporary <- tempfile("error_diagnostics_test_")
dir.create(temporary)
run_dir <- file.path(temporary, "run")
intermediate <- file.path(temporary, "intermediate")
dir.create(run_dir); dir.create(intermediate)
config_path <- file.path(temporary, "profile.yaml")
manifest_path <- file.path(temporary, "manifest.tsv")
yaml::write_yaml(config, config_path)
yaml::write_yaml(config, file.path(run_dir, "config_snapshot.yaml"))
utils::write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
expected_r <- as.character(getRversion())
versions <- data.frame(software = c("R", "dada2"), version = c(expected_r, "1.34.0"))
utils::write.table(versions, file.path(run_dir, "software_versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
provenance <- list(status = "complete", input_checksums_verified_at_completion = TRUE,
  study_id = study, bioproject = project, manifest_sample_count = 1L,
  random_seed = settings$random_seed, effective_threads = 2L, settings = settings,
  error_models_converged_before_limit = TRUE,
  command = c("Rscript", "scripts/dada2/run_pilot.R", "--config", config_path,
              "--manifest", manifest_path, "--output-dir", run_dir, "--intermediate-dir", intermediate))
provenance_path <- file.path(run_dir, "run_provenance.yaml")
yaml::write_yaml(provenance, provenance_path)
saveRDS(list(fixed, cycle), file.path(intermediate, "error_models.rds"))
args <- list(config = config_path, manifest = manifest_path, "run-dir" = run_dir,
  "intermediate-dir" = intermediate, output = file.path(temporary, "audit", "error_diagnostics.tsv"),
  "expected-r" = expected_r, "expected-dada2" = "1.34.0", "expected-threads" = 2L)
tryCatch({
  rows <- run_error_diagnostics(args, root)
  expect_true(nrow(rows) == 2L && all(rows$converged_before_limit) && !rows$last_input_equals_output[[2L]])
  written <- utils::read.delim(args$output, colClasses = "character", na.strings = character())
  expect_true(identical(written$training_bases, c("", "")))
  changed <- args; changed[["expected-dada2"]] <- "0.0.0"
  expect_error(run_error_diagnostics(changed, root), "version mismatch")
  changed <- args; changed[["expected-threads"]] <- 1L
  expect_error(run_error_diagnostics(changed, root), "threads mismatch")
  changed <- args; changed$output <- file.path(run_dir, "audit.tsv")
  expect_error(run_error_diagnostics(changed, root), "outside")
  changed <- args; changed$output <- manifest_path
  expect_error(run_error_diagnostics(changed, root), "overwrite an input")
  altered <- provenance; altered$settings$filter$min_len <- 149L
  yaml::write_yaml(altered, provenance_path)
  expect_error(run_error_diagnostics(args, root), "Executed settings")
  altered <- provenance; altered$random_seed <- 7L
  yaml::write_yaml(altered, provenance_path)
  expect_error(run_error_diagnostics(args, root), "seed mismatch")
  altered <- provenance; altered$manifest_sample_count <- 2L
  yaml::write_yaml(altered, provenance_path)
  expect_error(run_error_diagnostics(args, root), "sample count")
  altered <- provenance; altered$error_models_converged_before_limit <- FALSE
  yaml::write_yaml(altered, provenance_path)
  expect_error(run_error_diagnostics(args, root), "convergence indicator")
  altered <- provenance; altered$status <- "failed"
  yaml::write_yaml(altered, provenance_path)
  expect_error(run_error_diagnostics(args, root), "completed run")
  yaml::write_yaml(provenance, provenance_path)
  changed <- config; changed$project$default_random_seed <- 7L
  yaml::write_yaml(changed, file.path(run_dir, "config_snapshot.yaml"))
  expect_error(run_error_diagnostics(args, root), "snapshot differs")
}, finally = unlink(temporary, recursive = TRUE))
cat("Passed", checks, "synthetic error-diagnostic checks\n")
