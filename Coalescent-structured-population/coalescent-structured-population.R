library(ggplot2)
library(cowplot)
library(coala)
library(ggtree)
library(ape)

rm(list = ls())

pop_config <- c(5, 5)
migration_rate <- 2.0
mutation_rate <- 2.0
loci_number <- 1
locus_length <- 1000

get_char <- function(snps, haplos) {
  ret_val <- "NOT FOUND"
  for (i in 1:nrow(haplos)) {
    match <- TRUE
    for (j in 1:length(snps)) {
      ## j + 1 in haplos since letter is in the first column
      ##
      match <- match && (snps[j] == as.numeric(haplos[i, j + 1]))
    }
    if (match) {
      ret_val <- haplos[i, 1]
    }
  }
  return(ret_val)
}

model <- coal_model(sample_size = pop_config,
                    loci_number = loci_number,
                    loci_length = locus_length,
                    ploidy = 1) +
  feat_migration(rate = migration_rate,
                 symmetric = TRUE) +
  feat_mutation(rate = mutation_rate,
                model = "IFS") +
  sumstat_seg_sites() +
  sumstat_jsfs() +
  sumstat_trees()
results <- simulate(model)

## translate the tree into APE format
##
tree <- read.tree(text = results$trees[[1]])

snps <- results$seg_sites[[1]]$snps
haplos <- unique(snps)
## add first column as character to distinguish haplotypes
## Note: will fail if there are more than 26 haplotypes
##
haplos <- cbind(letters[seq(from = 1, to = nrow(haplos), by = 1)], haplos)
## associate haplotypes with tips
##
haplos_char <- character(nrow(snps))
for (i in 1:nrow(snps)) {
  haplos_char[i] <- get_char(snps[i,], haplos)
}
if (0) {
  ## collect information for the tree
  ##
  group_info <- split(tree$tip.label, haplos_char)
  snp_groups <- groupOTU(tree, group_info)
  ## plot the tree with different colors for different alleles
  ##
  p_snp <- ggtree(snp_groups, aes(color = group)) +
    scale_color_brewer(palette = "Dark2") +
    geom_tiplab(size = 1.5) +
    theme_tree2() +
    ggtitle("Colored by alleles")
}

## associate haplotypes with populations
##
pops <- character(0)
for (i in 1:length(pop_config)) {
  pops <- c(pops, rep(paste0("Pop", LETTERS[i]), pop_config[i]))
}
if (0) {
  ## collect information for the tree
  ##
  group_info <- split(tree$tip.label, pops)
  pop_groups <- groupOTU(tree, group_info)
  ## plot the tree with different colors for different populations
  ##
  p_pop <- ggtree(pop_groups, aes(color = group)) +
    scale_color_brewer(palette = "Dark2") +
    geom_tiplab(size = 1.5) +
    theme_tree2() +
    ggtitle("Colored by populations")
}

if (0) {
  print(plot_grid(p_snp, p_pop))
}

## construct base tree
##
p <- ggtree(tree) +
  theme_tree2()
## before adding allele and population information to it
##
for_plot <- data.frame(taxa = tree$tip.label,
                       Allele = haplos_char,
                       Population = pops)
p <- p %<+% for_plot +
  geom_tippoint(aes(color = Allele, shape = Population)) +
  theme(legend.position = "right")
print(p)
