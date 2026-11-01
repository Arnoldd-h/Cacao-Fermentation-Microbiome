#!/usr/bin/env Rscript
# Inputs: validated primary plot data and config. Outputs: PDF/SVG/300-dpi PNG.
# Fails on missing baselines, invalid distances, duplicate times or device errors.

main <- function() {
  args <- commandArgs(TRUE)
  if ("--help" %in% args) {
    cat("Usage: plot_trajectories.R --config PATH --input TSV --output-dir DIR\n",
        "Plots observed Aitchison distance from baseline by batch; exports PDF/SVG/PNG.\n",
        "Fails on invalid data or device errors; no curves fitted or gaps imputed.\n")
    return(invisible(NULL))
  }
  if (length(args) != 6L || !setequal(args[c(1, 3, 5)], c("--config", "--input", "--output-dir"))) stop("Invalid plotting CLI")
  cli <- stats::setNames(as.list(args[c(2, 4, 6)]), args[c(1, 3, 5)])
  config <- yaml::read_yaml(cli[["--config"]])
  data <- utils::read.delim(cli[["--input"]], check.names = FALSE)
  batches <- sort(unique(data$fermentation_batch))
  stopifnot(length(batches) == config$cohort$expected_batches, all(is.finite(data$distance_from_baseline)), all(data$distance_from_baseline >= 0))
  figures <- config$figures
  primary <- config$scenarios[[match(config$primary_scenario, vapply(config$scenarios, `[[`, character(1), "id"))]]
  strata_label <- switch(primary$strata_policy, complete = "ambos estratos", available = "estratos disponibles", single = paste("estrato", primary$stratum))
  for (extension in c("pdf", "svg", "png")) {
    destination <- file.path(cli[["--output-dir"]], paste0("community_trajectories.", extension))
    if (extension == "pdf") grDevices::pdf(destination, width = figures$width_inches, height = figures$height_inches)
    if (extension == "svg") grDevices::svg(destination, width = figures$width_inches, height = figures$height_inches)
    if (extension == "png") grDevices::png(destination, width = figures$width_inches, height = figures$height_inches, units = "in", res = figures$png_dpi, type = "cairo")
    columns <- ceiling(sqrt(length(batches)))
    graphics::par(mfrow = c(ceiling(length(batches) / columns), columns), mar = c(4, 4, 3, 1), oma = c(4, 1, 5, 1), mgp = c(2.2, 0.7, 0))
    for (batch in batches) {
      rows <- data[data$fermentation_batch == batch, , drop = FALSE]
      rows <- rows[order(rows$relative_time), , drop = FALSE]
      stopifnot(!anyDuplicated(rows$relative_time), sum(rows$fermentation_hours == 0) == 1)
      title <- sub(paste0("^", config$cohort$study_id, "::"), "", batch)
      graphics::plot(rows$relative_time, rows$distance_from_baseline, pch = 16, col = figures$color,
        xlim = c(0, 1), ylim = c(0, max(1, data$distance_from_baseline) * 1.05),
        xlab = "Tiempo relativo observado", ylab = "Distancia Aitchison desde 0 h", main = title)
      if (nrow(rows) > 1) for (i in seq_len(nrow(rows) - 1L)) {
        if (rows$batch_time_index[i + 1] - rows$batch_time_index[i] == 1) graphics::segments(
          rows$relative_time[i], rows$distance_from_baseline[i], rows$relative_time[i + 1], rows$distance_from_baseline[i + 1], col = figures$color)
      }
    }
    graphics::mtext("Cambio comunitario dentro del primer estudio", outer = TRUE, side = 3, line = 3, font = 2, cex = 1.15)
    graphics::mtext(paste(length(batches), "lotes independientes; centro CLR con igual peso por estrato"), outer = TRUE, side = 3, line = 1.6, cex = 0.85)
    graphics::mtext(paste("Principal:", strata_label, "- pseudoconteo", primary$pseudocount, "- huecos sin imputar"), outer = TRUE, side = 3, line = 0.4, cex = 0.8)
    graphics::mtext("Distancia desde el origen no demuestra una dirección compartida ni conservación entre estudios.", outer = TRUE, side = 1, line = 1, cex = 0.75)
    grDevices::dev.off()
  }
}

if (sys.nframe() == 0L) main()
