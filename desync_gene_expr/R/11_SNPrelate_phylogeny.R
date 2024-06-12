library(gdsfmt)
library(SNPRelate)

vcf_filename <- "data/variant_calling/final_filtered.recode.vcf"
gds_filename <- "outputs/variant_calling/final_filtered.gds"

snpgdsVCF2GDS(vcf_filename, gds_filename, 
              method = "biallelic.only")
# Summary
snpgdsSummary(gds_filename)
genofile <- snpgdsOpen(gds_filename)
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread = 2))

# Determine groups of individuals automatically
rv <- snpgdsCutTree(ibs.hc, label.H = TRUE)

svg("plots/variant_calling/SNPRelate_tree.svg",
    width = 11, height = 8)
plot(rv$dendrogram, leaflab = "perpendicular")
dev.off()
