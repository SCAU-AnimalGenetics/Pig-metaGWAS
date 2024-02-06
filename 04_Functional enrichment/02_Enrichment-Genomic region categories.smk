# 1. Functional enrichment with bedtools: 
rule run_enrichment:
    input:
        category="CDS.bed" # (i). Genomic categories: CDS, promoter, 5’UTR + 2kb upstream, 3’UTR + 2kb downstream, protein-coding genes, non-protein-coding genes, intron regions, (ii). Conserved elements, (iii). Each of 14 active chromatin states detected from 14 major pig tissues, (iv). Tissue-specific functional regions for each of 34 tissues in FarmGTEx
        univese_bed="metaGWAS_all_Unique_SNP.bed" # All unique SNPs in the meta-GWAS analyses
        test="D_ADG_5e_8.bed" #Significant/lead SNPs (P<5e-8) for each meta-GWAS summary
    output:
        outfile= "CDS_D_ADG_5e_8.bedtools.txt"
    shell:
        univese_num=`cat ${univese_bed}|wc -l`
        univese_int=`bedtools intersect -a ${category} -b ${univese_bed} -wb |wc -l`
		    query_num=`cat ${test}|wc -l`
		    query_int=`bedtools intersect -a ${category} -b ${test} -wb |wc -l`
            A=`echo "scale=4; ${query_int}/${query_num}" | bc`
            B=`echo "scale=4; ${univese_int}/${univese_num}" | bc`
            enrichment=`echo "scale=4; ${A}/${B}" | bc`
		    echo "LoFs" ${category} ${univese_num} ${univese_int} ${query_num} ${query_int} ${enrichment} >> ${outfile}


# 2. Permutation test for functional enrichment
rule run_permutation:
    R:
#input:
  category="CDS.bed" # (i). Genomic categories: CDS, promoter, 5’UTR + 2kb upstream, 3’UTR + 2kb downstream, protein-coding genes, non-protein-coding genes, intron regions, (ii). Conserved elements, (iii). Each of 14 active chromatin states detected from 14 major pig tissues, (iv). Tissue-specific functional regions for each of 34 tissues in FarmGTEx
  univese_bed="metaGWAS_all_Unique_SNP.bed" # All unique SNPs in the meta-GWAS analyses
  test="D_ADG_5e_8.bed" #Significant/lead SNPs (P<5e-8) for each meta-GWAS summary
#output:
  outfile="CDS_D_ADG_5e_8.permutation.txt"

library(regioneR)
library(data.table)
set.seed(1234)

A1 <- read.table(test,header=F)
A<-toGRanges(A1)
as_universe1 <- fread(univese_bed,header = F,sep = "\t")
as_universe<-toGRanges(as_universe1)
regionr_res<-data.frame(Class="",pval="",zscore="",numOverlaps="")
QTL_region <- read.table(category,header = F,sep = "\t")
B1 <- QTL_region[,c(1:3)]
B<-toGRanges(B1)
Class <- gsub('.bed', '', category)
pt.5000.reg <-permTest(A=A, ntimes=10000, randomize.function=resampleRegions, universe=as_universe, 
                           evaluate.function=numOverlaps, B=B, verbose=TRUE,force.parallel=TRUE)
res<-c(Class,pt.5000.reg$numOverlaps$pval,pt.5000.reg$numOverlaps$zscore,pt.5000.reg$numOverlaps$observed)
regionr_res<-rbind(regionr_res,res)
write.table(regionr_res,outfile,sep = "\t",row.names = F,quote = F)


# 3. Functional enrichment with LOLA R package
rule run_enrichment:
    R:
#input:
   category_dir="./Genomic_categories" #(i). Genomic categories: CDS, promoter, 5’UTR + 2kb upstream, 3’UTR + 2kb downstream, protein-coding genes, non-protein-coding genes, intron regions, (ii). Conserved elements, (iii). Each of 14 active chromatin states detected from 14 major pig tissues, (iv). Tissue-specific functional regions for each of 34 tissues in FarmGTEx
   univese_bed="metaGWAS_all_Unique_SNP.bed",
   test="D_ADG_5e_8.bed"
#output:
   outfile="D_ADG_5e_8.lola.csv"
    
library(LOLA)
library(GenomicRanges)
library(data.table)
univese_file<-fread(univese_bed,header=F)
univese_bed<-GRanges(seqnames =univese_file$V1,ranges = IRanges(univese_file$V2, univese_file$V3),strand = "*")
regionDB<-loadRegionDB(category_dir)
pop_bed<-read.table(test,sep="\t",header=F)
query_bed<-GRanges(seqnames =pop_bed$V1,ranges = IRanges(pop_bed$V2, pop_bed$V3),strand = "*")
locResults <-runLOLA(query_bed, univese_bed, regionDB, cores=1)
write.csv( locResults,outfile,row.names = F,quote = F)


