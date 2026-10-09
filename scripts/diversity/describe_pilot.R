#!/usr/bin/env Rscript
# Inputs: registered diversity config and validated bacterial ASV counts/metadata.
# Outputs: unrarefied alpha, CLR, distances, PCA, versions, warnings and figures.
# Fails on mixed studies, empty samples/features, mismatched metadata or R errors.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
if (!file.exists(script_file)) script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
root <- normalizePath(file.path(dirname(normalizePath(script_file, mustWork = TRUE)), "..", ".."))
source(file.path(root, "scripts/dada2/helpers.R"))
args <- commandArgs(TRUE)
if ("--help" %in% args || "-h" %in% args) {
  cat("Usage: Rscript scripts/diversity/describe_pilot.R --config PATH\n",
      "Inputs: registered JSON/YAML and one-study validated candidate ASVs/metadata.\n",
      "Outputs: descriptive alpha, CLR, Aitchison/Bray-Curtis, PCA and PDF/SVG/300-dpi PNG.\n",
      "Fails on inconsistent identifiers, empty libraries or calculation/device errors.\n",
      "Use run_diversity.py to check upstream provenance and register hashes/SUCCESS.\n")
  quit(status = 0)
}
if (length(args) != 2L || args[1] != "--config") fail("Invalid CLI; use --help")
config <- yaml::read_yaml(args[2])
settings <- config$analysis
input <- file.path(root, config$paths$bacterial_dir)
output <- file.path(root, config$paths$diversity_dir)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
rows <- utils::read.delim(file.path(input, "asv_counts.tsv"), check.names = FALSE)
metadata <- utils::read.delim(file.path(input, "sample_metadata.tsv"), colClasses = "character", check.names = FALSE)
require_fields(rows, c("study_id", "sample_id"), "counts")
require_fields(metadata, c("study_id", "sample_id", "fermentation_batch", "fermentation_hours", "relative_time", "fermentation_stage", "sampling_stratum", "run_accession"), "metadata")
if (length(unique(rows$study_id)) != 1L || anyDuplicated(rows$sample_id) || anyDuplicated(metadata$sample_id) ||
    !setequal(rows$sample_id, metadata$sample_id)) fail("Mixed studies or mismatched sample IDs")
metadata <- metadata[match(rows$sample_id, metadata$sample_id), , drop = FALSE]
if (any(metadata$study_id != rows$study_id)) fail("Metadata study mismatch")
counts <- as.matrix(rows[, -(1:2), drop = FALSE])
rownames(counts) <- rows$sample_id
if (!is.numeric(counts) || anyNA(counts) || any(counts < 0 | counts != floor(counts)) ||
    any(rowSums(counts) == 0) || any(colSums(counts) == 0) || nrow(counts) < 3 || ncol(counts) < 2) fail("Invalid counts or insufficient dimensions (at least three samples/two ASVs)")
set.seed(settings$random_seed)
warnings <- character()
export <- function(value, name) write_tsv(value, file.path(output, name))
withCallingHandlers({
  log_step("Describing ", nrow(counts), " samples and ", ncol(counts), " retained ASVs; inference disabled")
  sample_fields <- c("study_id", "sample_id", "run_accession", "fermentation_batch", "fermentation_hours", "relative_time", "fermentation_stage", "sampling_stratum")
  alpha <- metadata[, sample_fields, drop = FALSE]
  alpha$library_reads <- rowSums(counts)
  alpha$observed_asvs <- rowSums(counts > 0)
  alpha$shannon <- vegan::diversity(counts, index = "shannon", base = exp(1))
  alpha$gini_simpson <- vegan::diversity(counts, index = "simpson")
  alpha$inverse_simpson <- vegan::diversity(counts, index = "invsimpson")
  alpha$shannon_effective <- exp(alpha$shannon)
  export(alpha, "alpha_diversity.tsv")
  proportions <- counts / rowSums(counts)
  bray <- as.matrix(vegan::vegdist(proportions, method = "bray"))
  pairs <- expand.grid(first = seq_len(nrow(counts)), second = seq_len(nrow(counts)))
  distance_table <- function(matrix) data.frame(study_id = rows$study_id[1],
    sample_id_1 = rows$sample_id[pairs$first], sample_id_2 = rows$sample_id[pairs$second],
    distance = matrix[cbind(pairs$first, pairs$second)])
  export(distance_table(bray), "bray_curtis_distances.tsv")
  constants <- c(settings$primary_pseudocount, settings$sensitivity_pseudocounts)
  clr_tables <- distance_tables <- score_tables <- variance_tables <- fits <- list()
  for (i in seq_along(constants)) {
    logged <- log(counts + constants[i])
    clr <- logged - rowMeans(logged)
    cell <- expand.grid(sample = seq_len(nrow(counts)), feature = seq_len(ncol(counts)))
    clr_tables[[i]] <- data.frame(study_id = rows$study_id[1], pseudocount = constants[i],
      sample_id = rows$sample_id[cell$sample], asv_id = colnames(counts)[cell$feature],
      clr = clr[cbind(cell$sample, cell$feature)])
    distance_tables[[i]] <- cbind(pseudocount = constants[i], distance_table(as.matrix(stats::dist(clr))))
    fits[[i]] <- stats::prcomp(clr, center = settings$pca_center, scale. = settings$pca_scale)
    axes <- seq_len(min(nrow(counts) - 1L, ncol(counts)))
    # Resolve SVD sign ambiguity using the largest absolute loading per axis.
    for (axis in axes) {
      anchor <- which.max(abs(fits[[i]]$rotation[, axis]))
      if (fits[[i]]$rotation[anchor, axis] < 0) fits[[i]]$x[, axis] <- -fits[[i]]$x[, axis]
    }
    if (sum(fits[[i]]$sdev^2) <= 0) fail("PCA has no between-sample variation")
    axis_rows <- expand.grid(sample = seq_len(nrow(counts)), axis = axes)
    score_tables[[i]] <- data.frame(study_id = rows$study_id[1], pseudocount = constants[i],
      sample_id = rows$sample_id[axis_rows$sample], axis = axis_rows$axis,
      score = fits[[i]]$x[cbind(axis_rows$sample, axis_rows$axis)])
    variance_tables[[i]] <- data.frame(study_id = rows$study_id[1], pseudocount = constants[i], axis = axes,
      variance = fits[[i]]$sdev[axes]^2, explained_fraction = fits[[i]]$sdev[axes]^2 / sum(fits[[i]]$sdev^2))
  }
  export(do.call(rbind, clr_tables), "clr_coordinates.tsv")
  export(do.call(rbind, distance_tables), "aitchison_distances.tsv")
  export(do.call(rbind, score_tables), "pca_scores.tsv")
  export(do.call(rbind, variance_tables), "pca_variance.tsv")
  colors <- config$figures$stage_colors[metadata$fermentation_stage]
  colors <- unlist(lapply(colors, function(color) if (is.null(color)) "#777777" else color), use.names = FALSE)
  open_figure <- function(name, extension) {
    path <- file.path(output, paste0(name, ".", extension))
    width <- config$figures$width_inches; height <- config$figures$height_inches
    if (extension == "pdf") grDevices::pdf(path, width = width, height = height)
    if (extension == "svg") grDevices::svg(path, width = width, height = height)
    if (extension == "png") grDevices::png(path, width = width, height = height, units = "in", res = config$figures$png_dpi, type = "cairo")
  }
  for (extension in c("pdf", "svg", "png")) {
    open_figure("alpha_diversity", extension)
    graphics::par(mfrow = c(2, 2), mar = c(6, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (metric in c("library_reads", "observed_asvs", "shannon_effective", "inverse_simpson")) {
      graphics::plot(seq_len(nrow(alpha)), alpha[[metric]], xaxt = "n", xlab = "", ylab = metric,
        pch = 19, col = colors, main = gsub("_", " ", metric), las = 1)
      graphics::axis(1, at = seq_len(nrow(alpha)), labels = alpha$sample_id, las = 2, cex.axis = 0.7)
      graphics::legend("topright", legend = names(config$figures$stage_colors),
        col = unlist(config$figures$stage_colors), pch = 19, bty = "n", cex = 0.6)
    }
    graphics::mtext("Pilot descriptive diversity - unrarefied; samples remain separate", outer = TRUE)
    grDevices::dev.off()
    open_figure("aitchison_pca", extension)
    graphics::par(mfrow = c(1, length(constants)), mar = c(4, 4, 4, 2))
    for (i in seq_along(constants)) {
      xy <- fits[[i]]$x[, 1:2, drop = FALSE]
      fractions <- 100 * variance_tables[[i]]$explained_fraction[1:2]
      padded <- function(value) { span <- max(diff(range(value)), 1); range(value) + c(-1, 1) * 0.25 * span }
      graphics::plot(xy, xlim = padded(xy[, 1]), ylim = padded(xy[, 2]), pch = 19, col = colors,
        xlab = sprintf("PC1 (%.1f%%)", fractions[1]), ylab = sprintf("PC2 (%.1f%%)", fractions[2]),
        main = paste("CLR PCA - count pseudocount", constants[i]), las = 1)
      graphics::text(xy, labels = rows$sample_id, pos = 3, cex = 0.6)
      graphics::legend("bottomleft", legend = names(config$figures$stage_colors),
        col = unlist(config$figures$stage_colors), pch = 19, bty = "n", cex = 0.7)
    }
    grDevices::dev.off()
  }
}, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)) })
export(data.frame(software = c("R", "vegan", "yaml"), version = c(as.character(getRversion()),
  as.character(utils::packageVersion("vegan")), as.character(utils::packageVersion("yaml")))), "software_versions.tsv")
writeLines(capture.output(utils::sessionInfo()), file.path(output, "session_info.txt"))
writeLines(warnings, file.path(output, "warnings.txt"))
log_step("Descriptive diversity exported; no hypothesis tests performed")
