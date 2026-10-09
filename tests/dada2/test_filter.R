#!/usr/bin/env Rscript
# Synthetic SRA-style paired FASTQ regression. Never used as biological evidence.
if (!requireNamespace("dada2", quietly = TRUE)) stop("DADA2 is required for this integration check")
workspace <- tempfile("dada2_synthetic_")
dir.create(workspace)
sequence <- paste(rep("ACGT", 50), collapse = "")
write_fixture <- function(direction) {
  sequences <- rep(sequence, 3)
  if (direction == 2) substr(sequences[[2]], 20, 20) <- "N"
  records <- unlist(lapply(seq_along(sequences), function(index) {
    c(paste0("@SRR0.", index, " ", index, "/", direction), sequences[[index]], "+",
      paste(rep("I", nchar(sequences[[index]])), collapse = ""))
  }))
  path <- file.path(workspace, paste0("R", direction, ".fastq"))
  writeLines(records, path)
  path
}
forward <- write_fixture(1)
reverse <- write_fixture(2)
result <- dada2::filterAndTrim(forward, file.path(workspace, "F_filtered.fastq.gz"),
  reverse, file.path(workspace, "R_filtered.fastq.gz"), truncLen = c(0, 0),
  maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 150, rm.phix = FALSE,
  matchIDs = TRUE, id.field = 1, id.sep = "\\s", multithread = FALSE, compress = TRUE)
stopifnot(result[1, "reads.in"] == 3, result[1, "reads.out"] == 2)
filtered_f <- dada2::derepFastq(file.path(workspace, "F_filtered.fastq.gz"))
filtered_r <- dada2::derepFastq(file.path(workspace, "R_filtered.fastq.gz"))
stopifnot(sum(dada2::getUniques(filtered_f)) == 2, sum(dada2::getUniques(filtered_r)) == 2)
unlink(workspace, recursive = TRUE)
cat("DADA2 paired SRA-ID regression passed (synthetic fixture)\n")
