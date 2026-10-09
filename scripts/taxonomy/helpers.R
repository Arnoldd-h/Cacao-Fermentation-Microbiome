# Taxonomy masks retain lineage order and distinguish labels from named genera.

mask_taxonomy <- function(taxa, boots, threshold) {
  if (!is.matrix(taxa) || !identical(dim(taxa), dim(boots)) ||
      anyNA(boots) || any(!is.finite(boots) | boots < 0 | boots > 100)) {
    fail("Invalid taxonomy/bootstrap matrices")
  }
  result <- taxa
  for (i in seq_len(nrow(taxa))) {
    supported <- cumprod(!is.na(taxa[i, ]) & nzchar(taxa[i, ]) & boots[i, ] >= threshold) == 1
    result[i, !supported] <- NA_character_
  }
  result
}

screen_taxonomy <- function(taxa, settings) {
  kingdom <- taxa[, "Kingdom"]
  known <- !is.na(kingdom) & nzchar(kingdom)
  status <- rep("unclassified", nrow(taxa))
  status[known] <- "non_bacterial"
  status[known & kingdom == settings$bacterial_kingdom] <- "bacterial"
  organelle <- apply(taxa, 1, function(row) {
    matches <- settings$organelle_labels[tolower(settings$organelle_labels) %in% tolower(row)]
    paste(matches, collapse = ";")
  })
  genus <- taxa[, "Genus"]
  named <- !is.na(genus) & nzchar(genus) &
    !grepl(settings$unresolved_genus_pattern, genus, ignore.case = TRUE)
  data.frame(kingdom_status = status, organelle_flag = organelle,
             named_genus = named, excluded = FALSE, stringsAsFactors = FALSE)
}

taxonomy_coverage <- function(taxa, counts, study_id, threshold) {
  totals <- colSums(counts)
  data.frame(study_id = study_id, min_boot = threshold, rank = colnames(taxa),
    assigned_asvs = colSums(!is.na(taxa)), total_asvs = nrow(taxa),
    assigned_reads = vapply(seq_len(ncol(taxa)), function(j) sum(totals[!is.na(taxa[, j])]), numeric(1)),
    total_reads = sum(counts), row.names = NULL)
}

sample_coverage <- function(taxa, counts, study_id, threshold) {
  do.call(rbind, lapply(seq_len(nrow(counts)), function(i) {
    data.frame(study_id = study_id, sample_id = rownames(counts)[i], min_boot = threshold,
      rank = colnames(taxa),
      assigned_reads = vapply(seq_len(ncol(taxa)), function(j) sum(counts[i, !is.na(taxa[, j])]), numeric(1)),
      total_reads = sum(counts[i, ]), row.names = NULL)
  }))
}
