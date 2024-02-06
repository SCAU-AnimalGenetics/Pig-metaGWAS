
# 1. Make input file for metasoft
library('data.table')
library('stringr')
LeadSNP <- fread('./LeadSNP_info_5e_8.txt')
LeadSNP <- unique(LeadSNP$SNP) # all unique lead SNPs

All_Trait_info <- read.csv('./metaTable_232traits.csv')
colnames(All_Trait_info) <- c('TraitName') # 232 trait name info

for(j in 1:nrow(All_Trait_info)){
  # j=1
  if(j == 1){
    gwas <- fread(paste0('./GWAS_Summary/',All_Trait_info[j,]$TraitName,'.txt.gz'))
    gwas <- gwas[gwas$variant_id %in% LeadSNP, ]
    gwas <- gwas[,c('variant_id', 'effect_size', 'standard_error')]
    inputFile <- gwas
    colnames(inputFile) <- c('SNP', paste0('BETA_', j), paste0('SE_', j))
  }
  if(j != 1){
    gwas <- fread(paste0('./GWAS_Summary/',All_Trait_info[j,]$TraitName,'.txt.gz'))
    gwas <- gwas[gwas$variant_id %in% LeadSNP, ]
    gwas <- gwas[,c('variant_id', 'effect_size', 'standard_error')]
    colnames(gwas) <- c('SNP', paste0('BETA_', j), paste0('SE_', j))
    inputFile <- merge(inputFile, gwas, all = T)
  }
  
}

fwrite(inputFile, paste0('./inputFiles/232Traits.txt'),quote = F, row.names = F, col.names = F,eol = "\n", sep = '\t', na = 'NA')



# 2. perform metasoft
system(paste0('java -jar ./Metasoft/Metasoft.jar -input ./inputFiles/232Traits.txt -output ./outputFiles/232Traits.txt -pvalue_table ./Metasoft/HanEskinPvalueTable.txt -log ./LogFiles/232Traits.log -mvalue -mvalue_method mcmc -mvalue_p_thres 1 -seed 2023'))

