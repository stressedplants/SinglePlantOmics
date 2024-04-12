library(ape)

set.seed(123)

stree = nj(dist.gene(t(only_geno_numeric)))

# decide on correct tip colours for bolting
tip_cols <- sapply(conds_filtered$bolting, function(bolt_status) {
  if (bolt_status == "N") {
    return(no_col)
  } 
  
  if (bolt_status == "Y") {
    return(yes_col)
  }
})

# save this and load if already calculated
phylogeny_bootstrap_file <- "outputs/variant_calling/boot_phylo.Rdata"

if (!file.exists(phylogeny_bootstrap_file)) {
  boot_output <- boot.phylo(
    stree,
    only_geno_numeric,
    FUN = function(bootstrappedData) {
      nj(dist.gene(t(bootstrappedData)))
    },
    B = 1000
  )
  
  save(boot_output, file = phylogeny_bootstrap_file)
} else {
  load(phylogeny_bootstrap_file)
}

gc()

svg(filename = "plots/variant_calling/bootstrap_phylogeny.svg",
    width = 10, height = 8, pointsize = 7.5)
plot(stree, tip.color = tip_cols)
nodelabels(boot_output, frame = "none")
add.scale.bar(0, 0, length = 20)
dev.off()

set.seed(Sys.time())
