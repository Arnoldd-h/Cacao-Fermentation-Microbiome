#!/usr/bin/env Rscript
# Independently reconstruct CLR centers and check Python OLS/sign/BY results.
# Inputs: registered config, original bacterial counts/metadata, derived tables.
# Outputs: independent_validation.tsv, software_versions.tsv, session_info.txt.
# Fails on membership/coefficients/statistic discrepancies; no source is edited.

main <- function() {
  args <- commandArgs(TRUE)
  if ("--help" %in% args) {
    cat("Usage: validate_model.R --config PATH --counts TSV --metadata TSV --output-dir DIR\n",
        "Reconstructs models with lm.fit; checks signs/BY against exported tables.\n",
        "Writes validation/software/session summaries; fails on any disagreement.\n")
    return(invisible(NULL))
  }
  if (length(args) != 8L || !setequal(args[c(1, 3, 5, 7)], c("--config", "--counts", "--metadata", "--output-dir"))) stop("Invalid model validation CLI")
  cli <- stats::setNames(as.list(args[c(2, 4, 6, 8)]), args[c(1, 3, 5, 7)])
  config <- yaml::read_yaml(cli[["--config"]])
  output <- cli[["--output-dir"]]
  read <- function(path) utils::read.delim(path, colClasses = "character", check.names = FALSE, na.strings = character())
  source_counts <- read(cli[["--counts"]])
  meta <- read(cli[["--metadata"]])
  features <- setdiff(colnames(source_counts), c("study_id", "sample_id"))
  counts <- matrix(as.numeric(as.matrix(source_counts[, features, drop = FALSE])), nrow = nrow(source_counts), dimnames = list(source_counts$sample_id, features))
  stopifnot(!anyDuplicated(meta$sample_id), setequal(meta$sample_id, rownames(counts)), all(is.finite(counts)), all(counts >= 0))
  meta <- meta[match(rownames(counts), meta$sample_id), , drop = FALSE]
  meta$hours <- as.numeric(meta$fermentation_hours)
  meta$relative <- as.numeric(meta$relative_time)
  stopifnot(all(is.finite(meta$hours)), all(is.finite(meta$relative)))
  observations <- read(file.path(output, "observations.tsv"))
  membership <- read(file.path(output, "sample_membership.tsv"))
  models <- read(file.path(output, "batch_models.tsv"))
  slopes <- read(file.path(output, "slopes.tsv"))
  intercepts <- read(file.path(output, "intercepts.tsv"))
  globals <- read(file.path(output, "global_tests.tsv"))
  signs_exported <- read(file.path(output, "sign_statistics.tsv"))
  plot <- read(file.path(output, "trajectory_plot_data.tsv"))
  tolerance <- as.numeric(config$method$validation_tolerance)
  sign_tolerance <- as.numeric(config$method$comparison_tolerance)
  stopifnot(is.finite(tolerance), tolerance > 0, is.finite(sign_tolerance), sign_tolerance > 0)
  equal <- function(observed, expected, label) {
    observed <- as.numeric(observed)
    if (length(observed) != length(expected) || any(!is.finite(observed)) ||
        any(abs(observed - expected) > tolerance * pmax(1, abs(expected)))) stop(paste("Independent R mismatch:", label))
  }
  checks <- list()
  raw_p <- numeric(length(config$scenarios))
  for (scenario_index in seq_along(config$scenarios)) {
    scenario <- config$scenarios[[scenario_index]]
    id <- scenario$id
    raw_clr <- log(counts + scenario$pseudocount)
    raw_clr <- sweep(raw_clr, 1, rowMeans(raw_clr), "-")
    keep <- rep(TRUE, nrow(meta))
    for (batch in unique(meta$fermentation_batch)) for (hour in unique(meta$hours[meta$fermentation_batch == batch])) {
      index <- which(meta$fermentation_batch == batch & meta$hours == hour)
      if (scenario$strata_policy == "complete" && !setequal(meta$sampling_stratum[index], unlist(config$cohort$expected_strata))) keep[index] <- FALSE
    }
    if (scenario$strata_policy == "single") keep <- meta$sampling_stratum == scenario$stratum
    exported_members <- membership[membership$scenario_id == id, , drop = FALSE]
    stopifnot(nrow(exported_members) == nrow(meta), !anyDuplicated(exported_members$sample_id),
              setequal(exported_members$sample_id, meta$sample_id))
    exported_members <- exported_members[match(meta$sample_id, exported_members$sample_id), , drop = FALSE]
    stopifnot(identical(exported_members$included == "TRUE", keep))
    data <- observations[observations$scenario_id == id, , drop = FALSE]
    batches <- sort(unique(meta$fermentation_batch))
    stopifnot(setequal(data$fermentation_batch, batches))
    centers <- matrix(NA_real_, nrow = nrow(data), ncol = ncol(counts))
    x <- numeric(nrow(data))
    for (i in seq_len(nrow(data))) {
      index <- which(keep & meta$fermentation_batch == data$fermentation_batch[i] & meta$hours == as.numeric(data$fermentation_hours[i]))
      index <- index[order(meta$sampling_stratum[index], meta$sample_id[index])]
      stopifnot(length(index) > 0, identical(paste(meta$sample_id[index], collapse = ";"), data$source_sample_ids[i]))
      centers[i, ] <- colMeans(raw_clr[index, , drop = FALSE])
      x[i] <- (if (scenario$time_field == "relative_time") meta$relative[index[1]] else meta$hours[index[1]]) / scenario$time_scale
      equal(data$time_value[i], x[i], "time scale")
      equal(data$n_strata[i], length(index), "equal stratum count")
      equal(data$summed_source_reads[i], sum(counts[index, , drop = FALSE]), "source read audit")
    }
    betas <- matrix(NA_real_, nrow = length(batches), ncol = ncol(counts))
    for (b in seq_along(batches)) {
      index <- which(data$fermentation_batch == batches[b])
      response <- centers[index, , drop = FALSE]
      fit <- stats::lm.fit(cbind(1, x[index]), response)
      stopifnot(fit$rank == 2)
      betas[b, ] <- fit$coefficients[2, ]
      selector <- slopes$scenario_id == id & slopes$fermentation_batch == batches[b]
      stopifnot(sum(selector) == 1)
      equal(unlist(slopes[selector, features, drop = FALSE]), fit$coefficients[2, ], "QR slopes")
      selector <- intercepts$scenario_id == id & intercepts$fermentation_batch == batches[b]
      stopifnot(sum(selector) == 1)
      equal(unlist(intercepts[selector, features, drop = FALSE]), fit$coefficients[1, ], "QR intercepts")
      row <- models[models$scenario_id == id & models$fermentation_batch == batches[b], , drop = FALSE]
      stopifnot(nrow(row) == 1)
      rss <- sum(fit$residuals^2)
      total <- sum(sweep(response, 2, colMeans(response), "-")^2)
      equal(row$residual_sum_squares, rss, "QR residuals")
      equal(row$total_sum_squares, total, "total variation")
      if (row$r_squared != "") equal(row$r_squared, 1 - rss / total, "multivariate R squared")
      if (id == config$primary_scenario) {
        baseline <- which(as.numeric(data$fermentation_hours[index]) == 0)
        stopifnot(length(baseline) == 1)
        distances <- sqrt(rowSums(sweep(response, 2, response[baseline, ], "-")^2))
        for (j in seq_along(index)) {
          point <- plot[plot$fermentation_batch == batches[b] & as.numeric(plot$fermentation_hours) == as.numeric(data$fermentation_hours[index[j]]), , drop = FALSE]
          stopifnot(nrow(point) == 1)
          equal(point$distance_from_baseline, distances[j], "baseline Aitchison distance")
        }
      }
    }
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(batches))))
    null_statistics <- sqrt(rowSums((signs %*% betas / length(batches))^2))
    observed <- sqrt(sum(colMeans(betas)^2))
    p <- mean(null_statistics >= observed - sign_tolerance * max(1, observed))
    raw_p[scenario_index] <- p
    row <- globals[globals$scenario_id == id, , drop = FALSE]
    stopifnot(nrow(row) == 1)
    equal(row$norm_mean_slope, observed, "global norm")
    equal(row$exact_p_value, p, "exhaustive sign probability")
    equal(row$between_batch_rms, sqrt(sum(sweep(betas, 2, colMeans(betas), "-")^2) / nrow(betas)), "slope heterogeneity")
    exported <- signs_exported[signs_exported$scenario_id == id, , drop = FALSE]
    stopifnot(nrow(exported) == nrow(signs))
    equal(sort(as.numeric(exported$statistic)), sort(null_statistics), "complete sign distribution")
    checks[[scenario_index]] <- data.frame(scenario_id = id, status = "passed", asv_coefficients_checked = 2 * length(betas),
      observations_checked = nrow(data), batches_checked = length(batches))
  }
  adjusted <- stats::p.adjust(raw_p, method = config$method$fdr_method)
  for (i in seq_along(config$scenarios)) equal(globals$fdr_by_q_value[globals$scenario_id == config$scenarios[[i]]$id], adjusted[i], "BY FDR")
  utils::write.table(do.call(rbind, checks), file.path(output, "independent_validation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(software = c("R", "stats", "yaml"), version = c(as.character(getRversion()),
    as.character(utils::packageVersion("stats")), as.character(utils::packageVersion("yaml")))), file.path(output, "software_versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output, "session_info.txt"))
  message("Independent R QR, sign enumeration, BY and source-membership checks passed")
}

if (sys.nframe() == 0L) main()
