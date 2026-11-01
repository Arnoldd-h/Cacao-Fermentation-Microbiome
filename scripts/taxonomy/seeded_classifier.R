# Isolated native entry point; never modifies the installed dada2 namespace.
load_seeded_classifier <- function(root, settings) {
  policy <- settings$native_rng
  if (is.null(policy)) return(dada2::assignTaxonomy) # Validate historical runs.
  if (!identical(policy$adapter, "r_seeded_query_mt19937_v1") ||
      as.character(utils::packageVersion("dada2")) != policy$dada2_version) {
    stop("Unsupported seeded taxonomy adapter or DADA2 version")
  }
  native <- new.env(parent = globalenv())
  cache <- file.path(root, "results/cache/taxonomy_seeded_v1")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(root, "scripts/taxonomy/native/taxonomy_seeded.cpp"),
                  env = native, cacheDir = cache, rebuild = TRUE, showOutput = FALSE)
  classifier <- dada2::assignTaxonomy
  adapter_environment <- new.env(parent = environment(classifier))
  adapter_environment$C_assign_taxonomy2 <- native$C_assign_taxonomy2
  environment(classifier) <- adapter_environment
  classifier
}
