#!/usr/bin/env Rscript
# Synthetic fixtures validate contracts, not evidence about biological samples.
args <- commandArgs(trailingOnly = FALSE)
test_path <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1]]))
root <- normalizePath(file.path(dirname(test_path), "..", ".."))
source(file.path(root, "scripts", "dada2", "helpers.R"))
checks <- 0L
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  if (is.null(error) || !grepl(pattern, conditionMessage(error))) stop("Expected error: ", pattern)
  checks <<- checks + 1L
}
expect_true <- function(value) {
  stopifnot(isTRUE(value))
  checks <<- checks + 1L
}

manifest <- data.frame(study_id = c("fixture", "fixture"), bioproject = "PRJNA0",
  sample_id = c("s2", "s1"), run_accession = c("SRR2", "SRR1"), fermentation_stage = c("early", "late"),
  fermentation_hours = c("0", "24"), relative_time = c("0", "1"), fermentation_batch = "fixture::batch",
  library_layout = "PAIRED", sequencing_platform = "ILLUMINA")
expect_true(identical(validate_manifest(manifest), manifest))
changed <- manifest; changed$study_id[[2]] <- "another"
expect_error(validate_manifest(changed), "exactly one study")
changed <- manifest; changed$sample_id[[2]] <- changed$sample_id[[1]]
expect_error(validate_manifest(changed), "Duplicate")
changed <- manifest; changed$run_accession[[1]] <- "../escape"
expect_error(validate_manifest(changed), "Unsafe")
changed <- manifest; changed$relative_time[[1]] <- "NaN"
expect_error(validate_manifest(changed), "Invalid manifest time")
changed <- manifest; changed$library_layout[[1]] <- "SINGLE"
expect_error(validate_manifest(changed), "paired Illumina")

tracking <- data.frame(sample_id = c("s2", "s1"), input = c(100, 200), filtered = c(90, 180),
  denoised_f = c(80, 170), denoised_r = c(85, 160), merged = c(75, 155), nonchim = c(70, 150))
expect_true(identical(validate_tracking(tracking), tracking))
changed <- tracking; changed$merged[[2]] <- 165
expect_error(validate_tracking(changed), "conservation")
changed <- tracking; changed$nonchim[[1]] <- 0.5
expect_error(validate_tracking(changed), "Invalid read counts")
changed <- tracking; changed$nonchim[[2]] <- 156
expect_error(validate_tracking(changed), "conservation")

counts <- matrix(c(8L, 2L, 1L, 5L), nrow = 2,
                  dimnames = list(c("s2", "s1"), c("TTTT", "ACGT")))
exports <- sequence_exports(counts, "fixture")
expect_true(identical(exports$sequences$sequence, c("ACGT", "TTTT")))
expect_true(identical(exports$counts$sample_id, c("s2", "s1")))
expect_true(sum(exports$counts[, -(1:2)]) == sum(counts))
reordered <- sequence_exports(counts[, 2:1], "fixture")
expect_true(identical(exports, reordered))
changed <- counts; colnames(changed)[[1]] <- "NAAA"
expect_error(sequence_exports(changed, "fixture"), "Invalid ASV")
changed <- counts; changed[[1]] <- -1
expect_error(sequence_exports(changed, "fixture"), "Invalid ASV")
lengths <- length_summary(counts, "nonchim", "fixture")
expect_true(lengths$asv_count == 2L && lengths$read_count == 16L)

quality <- manifest[rep(seq_len(nrow(manifest)), each = 2), c("study_id", "sample_id", "run_accession")]
quality$read_direction <- rep(c("R1", "R2"), nrow(manifest))
quality$total_sequences <- c(100, 100, 200, 200)
temp <- tempfile(fileext = ".tsv")
write_tsv(quality[4:1, ], temp)
expect_true(identical(quality_counts(temp, manifest), c(100, 200)))
quality$total_sequences[[2]] <- 99
write_tsv(quality, temp)
expect_error(quality_counts(temp, manifest), "paired counts disagree")
quality$total_sequences[[2]] <- 100
quality$run_accession[[2]] <- "wrong_run"
write_tsv(quality, temp)
expect_error(quality_counts(temp, manifest), "metadata disagree")
unlink(temp)

settings <- list(filter = list(trunc_len = c(0, 0), max_n = 0, max_ee = c(2, 2), trunc_q = 2,
  min_len = 150, trim_left = c(0, 0), rm_phix = TRUE, match_ids = TRUE, id_field = 1, id_sep = "\\s"),
  error_learning = list(nbases = 1e8, randomize = FALSE, max_consist = 10, omega_c = 0),
  inference = list(pool = FALSE, omega_a = 1e-40, omega_c = 1e-40),
  merge = list(min_overlap = 12, max_mismatch = 0, trim_overhang = FALSE, just_concatenate = FALSE),
  chimera = list(method = "consensus", min_fold_parent_over_abundance = 1.5, min_parent_abundance = 2,
    min_sample_fraction = 0.9, ignore_n_negatives = 1, allow_one_off = FALSE),
  length_filter = list(enabled = FALSE), random_seed = 17, threads = 1)
expect_true(identical(validate_settings(settings, 17), settings))
changed <- settings; changed$inference$omega_a <- "1e-40"
expect_true(identical(normalize_settings(changed), settings))
changed$inference$omega_a <- "not a number"
expect_error(validate_settings(normalize_settings(changed), 17), "Invalid omega_a")
changed <- settings; changed$filter$max_n <- 1
expect_error(validate_settings(changed, 17), "max_n")
changed <- settings; changed$filter$match_ids <- FALSE
expect_error(validate_settings(changed, 17), "matching must remain")
changed <- settings; changed$inference$pool <- "pseudo"
expect_error(validate_settings(changed, 17), "independent sample")
changed <- settings; changed$merge$just_concatenate <- TRUE
expect_error(validate_settings(changed, 17), "non-overlapping")
expect_error(validate_settings(settings, 42), "random seeds disagree")
changed <- settings; changed$filter$trunc_len <- c(100, 100)
expect_error(validate_settings(changed, 17), "below min_len")
error_matrix <- matrix(0.1, nrow = 16, ncol = 10)
error_model <- list(err_in = list(error_matrix * 2, error_matrix), err_out = error_matrix,
                    trans = matrix(1L, nrow = 16, ncol = 10))
diagnostic <- error_convergence(error_model, "F", 10)
expect_true(diagnostic$converged_before_limit && diagnostic$self_consistency_rounds == 2)
diagnostic <- error_convergence(error_model, "F", 2)
expect_true(!diagnostic$converged_before_limit && diagnostic$reached_iteration_limit)
error_model$err_out <- error_matrix * 3
expect_true(!error_convergence(error_model, "F", 10)$repeated_error_matrix)
cat("DADA2 helper checks passed:", checks, "\n")
