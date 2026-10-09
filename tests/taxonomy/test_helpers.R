# Synthetic masks and flags only; no biological inference.
source("scripts/dada2/helpers.R")
source("scripts/taxonomy/helpers.R")
config <- yaml::read_yaml("config/taxonomy.yaml")
taxa <- matrix(c("Bacteria", "P", "C", "O", "F", "G"), nrow = 1,
               dimnames = list("fixture", config$classification$tax_levels))
boots <- matrix(c(100, 100, 100, 100, 100, 70), nrow = 1, dimnames = dimnames(taxa))
stopifnot(is.na(mask_taxonomy(taxa, boots, 80)[1, "Genus"]))
stopifnot(mask_taxonomy(taxa, boots, 50)[1, "Genus"] == "G")
boots[1, "Family"] <- 40
boots[1, "Genus"] <- 100
stopifnot(is.na(mask_taxonomy(taxa, boots, 50)[1, "Genus"]))
taxa[1, "Order"] <- "Chloroplast"
taxa[1, "Genus"] <- "uncultured bacterium"
flags <- screen_taxonomy(taxa, config$screening)
stopifnot(flags$organelle_flag == "Chloroplast", !flags$named_genus, !flags$excluded)
cat("6 taxonomy helper checks passed\n")
