#!/usr/bin/env Rscript
# Plot observed values only; no fitting, aggregation or missing-point imputation.

adjacent_pairs <- function(times, batch_times) {
  if (anyDuplicated(times) || is.unsorted(times, strictly = TRUE) || any(!times %in% batch_times)) {
    stop("Invalid temporal series order or membership")
  }
  if (length(times) < 2L) return(matrix(integer(), ncol = 2))
  pairs <- cbind(seq_len(length(times) - 1L), seq.int(2L, length(times)))
  pairs[diff(match(times, batch_times)) == 1L, , drop = FALSE]
}

render_temporal <- function(data, config, profile, destination, extension) {
  figures <- config$figures
  if (extension == "pdf") grDevices::pdf(destination, width = figures$width_inches, height = figures$height_inches)
  if (extension == "svg") grDevices::svg(destination, width = figures$width_inches, height = figures$height_inches)
  if (extension == "png") grDevices::png(destination, width = figures$width_inches, height = figures$height_inches,
                                        units = "in", res = figures$png_dpi, type = "cairo")
  on.exit(grDevices::dev.off())
  batches <- sort(unique(data$fermentation_batch))
  strata <- sort(unique(data$sampling_stratum))
  colors <- stats::setNames(unlist(figures$stratum_colors)[seq_along(strata)], strata)
  metrics <- vapply(config$metrics, `[[`, character(1), "field")
  graphics::par(mfrow = c(length(batches), length(metrics)), mar = c(3, 4, 3, 1), oma = c(5, 1, 5, 1))
  graphics::par(cex = 1, ps = 10, mgp = c(2, 0.6, 0))
  study_prefix <- paste0(unique(data$study_id), "::")
  for (batch in batches) {
    subset <- data[data$fermentation_batch == batch, , drop = FALSE]
    batch_times <- sort(unique(subset$fermentation_hours))
    title <- if (startsWith(batch, study_prefix)) substring(batch, nchar(study_prefix) + 1L) else batch
    for (metric_index in seq_along(metrics)) {
      metric <- metrics[[metric_index]]
      limits <- range(c(0, data[[metric]]))
      single_time <- length(batch_times) == 1L
      time_limits <- if (single_time) c(max(0, batch_times - 1), batch_times + 1) else range(batch_times)
      graphics::plot(NA, xlim = time_limits, ylim = limits, xlab = "", ylab = "",
                     main = paste(title, config$metrics[[metric_index]]$label, sep = "\n"), cex.main = 0.85,
                     cex.lab = 0.8, cex.axis = 0.75, xaxt = if (single_time) "n" else "s")
      if (single_time) graphics::axis(1, at = batch_times, cex.axis = 0.75)
      for (stratum in strata) {
        series <- subset[subset$sampling_stratum == stratum, , drop = FALSE]
        if (!nrow(series)) next
        series <- series[order(series$fermentation_hours), , drop = FALSE]
        if (isTRUE(profile$connect_points)) {
          pairs <- adjacent_pairs(series$fermentation_hours, batch_times)
          if (nrow(pairs)) graphics::segments(series$fermentation_hours[pairs[, 1]], series[[metric]][pairs[, 1]],
                                             series$fermentation_hours[pairs[, 2]], series[[metric]][pairs[, 2]],
                                             col = colors[[stratum]], lwd = 1.1)
        }
        graphics::points(series$fermentation_hours, series[[metric]], col = colors[[stratum]], pch = 16, cex = 0.7)
      }
    }
  }
  graphics::mtext("Diversidad observada por lote y estrato", outer = TRUE, side = 3, line = 3, font = 2, cex = 1.1)
  graphics::mtext(profile$scope_label, outer = TRUE, side = 3, line = 1.6, cex = 0.85)
  note <- if (isTRUE(profile$connect_points)) "Líneas entre observaciones consecutivas del mismo estrato; huecos sin imputar. Sin pruebas estadísticas." else "Puntos aislados: el piloto no reconstruye trayectorias. Sin pruebas estadísticas."
  graphics::mtext(note, outer = TRUE, side = 3, line = 0.4, cex = 0.7)
  graphics::mtext("Horas de fermentación", outer = TRUE, side = 1, line = 1, cex = 0.85)
  graphics::par(fig = c(0, 1, 0, 0.035), new = TRUE, mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend("center", legend = paste("Estrato", strata), col = colors, pch = 16,
                   lty = if (isTRUE(profile$connect_points)) 1 else 0, horiz = TRUE, bty = "n", cex = 0.9)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat("Usage: plot_temporal_diversity.R --config PATH --analysis KEY --input TSV --output-dir DIR\n",
        "Input: validated, ordered plot_data.tsv and policy; outputs: PDF/SVG/PNG/session info.\n",
        "Fails on missing arguments, invalid data/config or device errors.\n")
    return(invisible(NULL))
  }
  if (length(args) != 8L || !setequal(args[c(1, 3, 5, 7)], c("--config", "--analysis", "--input", "--output-dir"))) stop("Invalid temporal plotting arguments")
  cli_options <- stats::setNames(as.list(args[c(2, 4, 6, 8)]), args[c(1, 3, 5, 7)])
  config <- yaml::read_yaml(cli_options[["--config"]])
  profile <- config$analyses[[cli_options[["--analysis"]]]]
  if (is.null(profile) || !isFALSE(config$inference) || !isTRUE(config$break_at_missing_batch_times) ||
      config$sample_aggregation != "none" || config$x_field != "fermentation_hours" ||
      !identical(unlist(config$group_fields), c("fermentation_batch", "sampling_stratum"))) stop("Unsupported temporal policy")
  data <- utils::read.delim(cli_options[["--input"]], colClasses = "character", check.names = FALSE, stringsAsFactors = FALSE)
  fields <- c("fermentation_hours", vapply(config$metrics, `[[`, character(1), "field"))
  for (field in fields) {
    data[[field]] <- as.numeric(data[[field]])
    if (any(!is.finite(data[[field]]))) stop("Nonfinite temporal data")
  }
  if (nrow(data) != profile$expected_samples || length(unique(data$study_id)) != 1L || length(unique(data$fermentation_batch)) != profile$expected_batches) stop("Temporal scope mismatch")
  for (extension in c("pdf", "svg", "png")) render_temporal(data, config, profile,
    file.path(cli_options[["--output-dir"]], paste0("temporal_diversity.", extension)), extension)
  writeLines(capture.output(utils::sessionInfo()), file.path(cli_options[["--output-dir"]], "session_info.txt"))
}

if (sys.nframe() == 0L) main()
