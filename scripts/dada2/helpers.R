# Validation and tabulation helpers for independently processed DADA2 studies.

fail <- function(...) stop(..., call. = FALSE)

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), " ", ...)
}

require_fields <- function(value, fields, label) {
  missing <- setdiff(fields, names(value))
  if (length(missing)) fail(label, " missing fields: ", paste(missing, collapse = ", "))
}

check_number <- function(value, label, minimum = 0, maximum = Inf,
                         size = 1L, integer = FALSE) {
  if (!is.numeric(value) || length(value) != size || anyNA(value) ||
      any(!is.finite(value)) || any(value < minimum | value > maximum) ||
      (integer && any(value != floor(value)))) fail("Invalid ", label)
}

check_boolean <- function(value, label) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) fail("Invalid ", label)
}

normalize_settings <- function(settings) {
  # JSON-compatible YAML uses 1e-40; YAML 1.1 resolves that spelling as text.
  # Convert only the three probability fields, rejecting arbitrary strings later.
  for (section in c("error_learning", "inference")) {
    for (key in intersect(c("omega_a", "omega_c"), names(settings[[section]]))) {
      value <- settings[[section]][[key]]
      if (is.character(value) && length(value) == 1L && !is.na(value) &&
          grepl("^[0-9]+([.][0-9]+)?[eE][-+]?[0-9]+$", value)) {
        settings[[section]][[key]] <- as.numeric(value)
      }
    }
  }
  settings
}

validate_settings <- function(settings, project_seed) {
  require_fields(settings, c("filter", "error_learning", "inference", "merge",
                            "chimera", "length_filter", "random_seed", "threads"), "dada2")
  filt <- settings$filter
  require_fields(filt, c("trunc_len", "max_n", "max_ee", "trunc_q", "min_len",
                        "trim_left", "rm_phix", "match_ids", "id_field", "id_sep"), "filter")
  check_number(filt$trunc_len, "trunc_len", size = 2L, integer = TRUE)
  check_number(filt$trim_left, "trim_left", size = 2L, integer = TRUE)
  check_number(filt$max_ee, "max_ee", minimum = .Machine$double.eps, size = 2L)
  check_number(filt$max_n, "max_n", maximum = 0, integer = TRUE)
  check_number(filt$trunc_q, "trunc_q", maximum = 93, integer = TRUE)
  check_number(filt$min_len, "min_len", minimum = 1, integer = TRUE)
  check_boolean(filt$rm_phix, "rm_phix")
  check_boolean(filt$match_ids, "match_ids")
  if (!filt$match_ids) fail("Pair identifier matching must remain enabled")
  check_number(filt$id_field, "id_field", minimum = 1, integer = TRUE)
  if (!is.character(filt$id_sep) || length(filt$id_sep) != 1L || !nzchar(filt$id_sep)) fail("Invalid id_sep")
  if (any(filt$trunc_len > 0 & filt$trunc_len - filt$trim_left < filt$min_len)) {
    fail("trunc_len minus trim_left cannot be below min_len")
  }
  err <- settings$error_learning
  require_fields(err, c("nbases", "randomize", "max_consist", "omega_c"), "error_learning")
  check_number(err$nbases, "nbases", minimum = 1, integer = TRUE)
  check_number(err$max_consist, "max_consist", minimum = 2, integer = TRUE)
  check_number(err$omega_c, "error_learning.omega_c", maximum = 1)
  check_boolean(err$randomize, "randomize")
  require_fields(settings$inference, c("pool", "omega_a", "omega_c"), "inference")
  if (!identical(settings$inference$pool, FALSE)) fail("Only independent sample inference is supported")
  for (key in c("omega_a", "omega_c")) {
    check_number(settings$inference[[key]], key, minimum = .Machine$double.xmin, maximum = 1)
  }
  merge <- settings$merge
  require_fields(merge, c("min_overlap", "max_mismatch", "trim_overhang", "just_concatenate"), "merge")
  check_number(merge$min_overlap, "min_overlap", minimum = 1, integer = TRUE)
  check_number(merge$max_mismatch, "max_mismatch", integer = TRUE)
  check_boolean(merge$trim_overhang, "trim_overhang")
  if (!identical(merge$just_concatenate, FALSE)) fail("Concatenating non-overlapping reads is unsupported")
  chim <- settings$chimera
  require_fields(chim, c("method", "min_fold_parent_over_abundance", "min_parent_abundance",
                        "min_sample_fraction", "ignore_n_negatives", "allow_one_off"), "chimera")
  if (!identical(chim$method, "consensus")) fail("Only consensus chimera detection is supported")
  check_number(chim$min_fold_parent_over_abundance, "min_fold_parent_over_abundance", minimum = 1)
  check_number(chim$min_parent_abundance, "min_parent_abundance", minimum = 1, integer = TRUE)
  check_number(chim$min_sample_fraction, "min_sample_fraction", minimum = 0, maximum = 1)
  check_number(chim$ignore_n_negatives, "ignore_n_negatives", integer = TRUE)
  check_boolean(chim$allow_one_off, "allow_one_off")
  if (!identical(settings$length_filter$enabled, FALSE)) fail("Post-merge length filtering is not implemented")
  check_number(settings$random_seed, "random_seed", maximum = .Machine$integer.max, integer = TRUE)
  check_number(settings$threads, "threads", minimum = 1, integer = TRUE)
  if (!identical(as.numeric(settings$random_seed), as.numeric(project_seed))) {
    fail("DADA2 and project random seeds disagree")
  }
  invisible(settings)
}

validate_manifest <- function(manifest) {
  required <- c("study_id", "bioproject", "sample_id", "run_accession", "fermentation_stage",
                "fermentation_hours", "relative_time", "fermentation_batch", "library_layout",
                "sequencing_platform")
  require_fields(manifest, required, "manifest")
  if (!nrow(manifest)) fail("Empty manifest")
  for (field in required) {
    if (anyNA(manifest[[field]]) || any(!nzchar(trimws(manifest[[field]])))) fail("Empty manifest field: ", field)
  }
  if (length(unique(manifest$study_id)) != 1 || length(unique(manifest$bioproject)) != 1) {
    fail("A DADA2 execution must contain exactly one study and BioProject")
  }
  if (anyDuplicated(manifest$sample_id) || anyDuplicated(manifest$run_accession)) {
    fail("Duplicate sample_id or run_accession; combine technical runs only via a documented upstream rule")
  }
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", manifest$run_accession)) ||
      any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", manifest$study_id))) fail("Unsafe manifest identifier")
  if (any(manifest$library_layout != "PAIRED") || any(manifest$sequencing_platform != "ILLUMINA")) {
    fail("This workflow accepts paired Illumina amplicons only")
  }
  hours <- suppressWarnings(as.numeric(manifest$fermentation_hours))
  relative_time <- suppressWarnings(as.numeric(manifest$relative_time))
  if (anyNA(hours) || any(!is.finite(hours) | hours < 0) || anyNA(relative_time) ||
      any(!is.finite(relative_time) | relative_time < 0 | relative_time > 1)) fail("Invalid manifest time")
  if (any(!manifest$fermentation_stage %in% c("early", "mid", "late"))) fail("Invalid fermentation stage")
  invisible(manifest)
}

write_tsv <- function(value, path) {
  utils::write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

validate_tracking <- function(tracking) {
  stages <- c("input", "filtered", "denoised_f", "denoised_r", "merged", "nonchim")
  require_fields(tracking, c("sample_id", stages), "tracking")
  counts <- as.matrix(tracking[, stages, drop = FALSE])
  if (anyNA(counts) || any(!is.finite(counts) | counts < 0 | counts != floor(counts))) fail("Invalid read counts")
  if (any(tracking$filtered > tracking$input) ||
      any(tracking$denoised_f > tracking$filtered | tracking$denoised_r > tracking$filtered) ||
      any(tracking$merged > pmin(tracking$denoised_f, tracking$denoised_r)) ||
      any(tracking$nonchim > tracking$merged)) fail("Read conservation failed")
  invisible(tracking)
}

sequence_exports <- function(sequence_table, study_id) {
  if (!is.matrix(sequence_table) || !ncol(sequence_table) || !nrow(sequence_table) ||
      is.null(colnames(sequence_table)) || is.null(rownames(sequence_table)) ||
      anyNA(sequence_table) || any(!is.finite(sequence_table) | sequence_table < 0 |
                                   sequence_table != floor(sequence_table)) ||
      any(!grepl("^[ACGT]+$", colnames(sequence_table))) ||
      anyDuplicated(colnames(sequence_table)) || anyDuplicated(rownames(sequence_table))) {
    fail("Invalid ASV count matrix")
  }
  # Lexical sequence order is deterministic and independent of sample abundance.
  sequence_table <- sequence_table[, order(colnames(sequence_table), method = "radix"), drop = FALSE]
  ids <- paste0(study_id, "__ASV", sprintf("%06d", seq_len(ncol(sequence_table))))
  sequences <- data.frame(study_id = study_id, asv_id = ids, sequence = colnames(sequence_table),
                          length = nchar(colnames(sequence_table)), total_reads = colSums(sequence_table))
  colnames(sequence_table) <- ids
  counts <- data.frame(study_id = study_id, sample_id = rownames(sequence_table),
                       sequence_table, check.names = FALSE, row.names = NULL)
  list(sequences = sequences, counts = counts)
}

length_summary <- function(sequence_table, stage, study_id) {
  lengths <- nchar(colnames(sequence_table))
  data.frame(study_id = study_id, stage = stage, length = sort(unique(lengths)),
             asv_count = vapply(sort(unique(lengths)), function(n) sum(lengths == n), integer(1)),
             read_count = vapply(sort(unique(lengths)), function(n) sum(sequence_table[, lengths == n, drop = FALSE]), numeric(1)))
}

relative_path <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  prefix <- paste0(normalizePath(root, winslash = "/", mustWork = TRUE), "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}

input_paths <- function(manifest, input_dir) {
  paths <- lapply(c("1", "2"), function(direction) {
    file.path(input_dir, manifest$run_accession,
              paste0(manifest$run_accession, "_", direction, ".fastq.gz"))
  })
  missing <- unlist(paths)[!file.exists(unlist(paths))]
  if (length(missing)) fail("Missing trimmed FASTQ: ", paste(missing, collapse = ", "))
  lapply(paths, function(value) stats::setNames(value, manifest$sample_id))
}

quality_counts <- function(path, manifest) {
  quality <- utils::read.delim(path, colClasses = "character", check.names = FALSE)
  require_fields(quality, c("study_id", "sample_id", "run_accession", "read_direction", "total_sequences"), "QC table")
  keys <- paste(quality$sample_id, quality$read_direction)
  expected_f <- paste(manifest$sample_id, "R1")
  expected_r <- paste(manifest$sample_id, "R2")
  if (anyDuplicated(keys) || !setequal(keys, c(expected_f, expected_r))) fail("QC sample/direction keys disagree with manifest")
  counts <- suppressWarnings(as.numeric(quality$total_sequences))
  if (anyNA(counts) || any(!is.finite(counts) | counts < 0 | counts != floor(counts))) fail("Invalid QC read counts")
  matched <- match(quality$sample_id, manifest$sample_id)
  if (any(quality$study_id != manifest$study_id[matched]) ||
      any(quality$run_accession != manifest$run_accession[matched])) fail("QC metadata disagree with manifest")
  forward <- counts[match(expected_f, keys)]
  if (any(forward != counts[match(expected_r, keys)])) fail("QC paired counts disagree")
  forward
}

error_convergence <- function(model, direction, max_consist) {
  if (!is.list(model$err_in) || !length(model$err_in) || !is.matrix(model$err_out)) {
    fail("Missing error iteration history: ", direction)
  }
  # DADA2 1.34 stops on an exactly repeated error matrix (including a cycle).
  repeated <- any(vapply(model$err_in, identical, logical(1), model$err_out))
  rounds <- length(model$err_in)
  data.frame(read_direction = direction, self_consistency_rounds = rounds,
    maximum_rounds = max_consist, repeated_error_matrix = repeated,
    reached_iteration_limit = rounds >= max_consist,
    converged_before_limit = repeated && rounds < max_consist,
    transition_base_observations = sum(model$trans))
}
