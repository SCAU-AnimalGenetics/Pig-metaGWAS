options(stringsAsFactors = F)
setwd('./SigSNPFiles')
library('data.table')
library(stringr)

# 01. SNPs Genomic Variants Annotation
SNPeff <- 'Annofile.txt.gz' # SNPeff annotation file
gwas <- 'D_ADG_5e_8.txt' # GWAS summary of significant SNPs (P<5e-8) from D_ADG trait
annoResult <- 'D_ADG_Genomic_Variants_Annotation.txt'
#-- read annotation file
geneAnnoFile <- fread(SNPeff)
colnames(geneAnnoFile) <- c('CHR', 'POS', 'GenomeType')
geneAnnoFile$SNP <- paste0(geneAnnoFile$CHR, '_', geneAnnoFile$POS)
#-- read GWAS file
gwasSummary <- fread(gwas)[,c('chromosome', 'position')]
colnames(gwasSummary) <- c('CHR', 'POS')
gwasSummary$CHR <- gsub('chr', '', gwasSummary$CHR)
#-- SNPs annotation
AnnoResult <- geneAnnoFile[geneAnnoFile$SNP %in% paste0(gwasSummary$CHR, '_', gwasSummary$POS), ]
write.table(AnnoResult, paste0('./GeneAnnoResult/', annoResult ),quote = F, row.names = F, col.names = T,eol = "\n", sep = '\t')



# 02. SNPs Genomic Variants Enrichment
library(fmsb)
allGWASSNPs <- 'metaGWAS_all_Unique_SNPs_annotation_result.txt'
enrichResult <- 'D_ADG_Genomic_Variants_Enrichment.txt'
#-- read annotation result of all unique SNPs in meta-GWAS
metaGWASanno <- fread(allGWASSNPs)

AnnoResult <- unique(AnnoResult)
SNPtype <- as.data.frame(table(AnnoResult$GenomeType))
colnames(SNPtype) <- c('snpType', 'Freq')
SNPtype$pro <- SNPtype$Freq/nrow(AnnoResult) # 计算基因组注释比例
SNPtype$Trait <- gsub('_5e_8.txt', '', gwas)

#-- enrichment analysis using the 'fmsb' package
for(i in 1:nrow(SNPtype)){
  #i=6
  OA <- oddsratio(as.numeric(SNPtype[i, 'Freq']), sum(SNPtype$Freq),as.numeric(metaGWASanno[metaGWASanno$snpType == SNPtype[i,'snpType'], 'Freq']) ,sum(metaGWASanno$Freq))
  
  SNPtype[i,'EnrichMultip'] <- OA$estimate
  SNPtype[i,'Pvalue'] <- OA$p.value
  SNPtype[i,'95 percent confidence interval 1'] <- OA$conf.int[1]
  SNPtype[i,'95 percent confidence interval 2'] <- OA$conf.int[2]
}

fwrite(SNPtype, paste0('./GeneAnnoResult/', enrichResult),quote = F, row.names = F, col.names = T,eol = "\n", sep = '\t')

