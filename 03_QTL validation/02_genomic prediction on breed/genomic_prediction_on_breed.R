library(ggplot2)
library(dplyr)
wkdir = "./Genomic_prediction/code/example"
setwd(wkdir)

#============================================
# transfer file format from .vcf to 012 using vcftools
system("vcftools --gzvcf leadSNP.vcf.gz --012 --out leadSNP")

# get genotype SNP info
geno = read.table("leadSNP.012")
iid = geno$V1
geno = geno[, -1]
geno = as.matrix(geno)
geno_indiID = read.table("leadSNP.012.indv")
geno_snpID = read.table("leadSNP.012.pos")
geno_snpID = paste0(geno_snpID$V1, "_", geno_snpID$V2)
#============================================
## get GWAS summary 
GWAS_file = "leadSNP_gwas_summary.txt"
GWASSummary = read.table(GWAS_file)
## match SNP to confirm the order between genotype and effect size
GWAS_snpID = paste0(gsub("chr", "", GWASSummary$V3), "_", GWASSummary$V4)
GWASSummary = GWASSummary[match(geno_snpID, GWAS_snpID), ]

beta = GWASSummary$V10
beta = as.matrix(beta, length(beta), 1)

predict_value = geno %*% beta

results = data.frame(id = iid, predict_value = predict_value[, 1])
