#!/usr/bin/env Rscript
# Synthetic, deliberately tied references. No biological observations.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
source("scripts/taxonomy/seeded_classifier.R")
settings <- yaml::read_yaml("config/taxonomy.yaml")$classification
original_entry <- get("C_assign_taxonomy2", envir = asNamespace("dada2"))
classifier <- load_seeded_classifier(normalizePath("."), settings)
stopifnot(identical(original_entry, get("C_assign_taxonomy2", envir = asNamespace("dada2"))))
reference <- tempfile(fileext = ".fa")
# Two equal-sized genera with identical training sequences force exact score ties.
training <- paste(rep("ACGT", 30), collapse = "")
writeLines(unlist(lapply(seq_len(12), function(i)
  c(paste0(">Bacteria;P;C;O;F;G", 1L + (i %% 2L), ";"), training))), reference)
set.seed(17)
queries <- replicate(300, paste(sample(c("A", "C", "G", "T"), 120, TRUE), collapse = ""))
do.call(RNGkind, as.list(settings$rng_kind))
set.seed(settings$random_seed)
ambiguous <- classifier(queries, reference, minBoot = 0, tryRC = TRUE,
  outputBootstraps = TRUE, taxLevels = settings$tax_levels, multithread = as.integer(args[2]))
stopifnot(length(unique(ambiguous$tax[, "Genus"])) == 2L,
          all(ambiguous$boot[, "Family"] == 100L),
          all(ambiguous$boot[, "Genus"] > 10L & ambiguous$boot[, "Genus"] < 90L))
set.seed(settings$random_seed + 1L)
other_seed <- classifier(queries, reference, minBoot = 0, tryRC = TRUE,
  outputBootstraps = TRUE, taxLevels = settings$tax_levels, multithread = as.integer(args[2]))
stopifnot(!identical(ambiguous, other_seed))
# No tied lineages: exact agreement with the installed unmodified entry point.
writeLines(unlist(lapply(seq_len(12), function(i)
  c(">Bacteria;P;C;O;F;G;", training))), reference)
set.seed(settings$random_seed)
controlled <- classifier(queries[1:3], reference, minBoot = 0, tryRC = TRUE,
  outputBootstraps = TRUE, taxLevels = settings$tax_levels, multithread = as.integer(args[2]))
set.seed(settings$random_seed)
upstream <- dada2::assignTaxonomy(queries[1:3], reference, minBoot = 0, tryRC = TRUE,
  outputBootstraps = TRUE, taxLevels = settings$tax_levels, multithread = as.integer(args[2]))
stopifnot(identical(controlled, upstream))
saveRDS(ambiguous, args[1], compress = FALSE)
unlink(reference)
