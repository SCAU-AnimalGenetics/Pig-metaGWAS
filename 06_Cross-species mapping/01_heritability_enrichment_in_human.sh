# Input files and parameters
bedfile="" # including chr, start and end, see LDSC in detail
chr=1
sumstat="" # summary statistics, see LDSC in detail

#------------------------------------------------------------
# step 1: calculate the ld score using annotation
# Create annotfile -----
python ./make_annot.py \
--bed-file ${bedfile} \
--bimfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
--annot-file ./1000G.EUR.QC.${chr}.annot.gz

# Add base Column 
gunzip ${annotDir}/1000G.EUR.QC.${chr}.annot.gz
sed "s/^/1\t&/g" ${annotDir}/1000G.EUR.QC.${chr}.annot > ${annotDir}/1000G.EUR.QC.${chr}.annot.s1   # add base line
awk '{if(1==NR){gsub(/.*/,"base", $1)}print}' ${annotDir}/1000G.EUR.QC.${chr}.annot.s1 > ${annotDir}/1000G.EUR.QC.${chr}.annot   # modify the colname
gzip ${annotDir}/1000G.EUR.QC.${chr}.annot

# Compute LDscore -----
python ${ldscPath}/ldsc.py \
--l2 \
--bfile ./1000G.EUR.QC.${chr} \
--ld-wind-cm 1 \
--annot ${annotDir}/1000G.EUR.QC.${chr}.annot.gz \
--thin-annot \
--out ./1000G.EUR.QC.${chr} \
--print-snps ./hm.${chr}.snp

#-----------------------------------------------------------------
# step 2: heritability enrichment
${ldscPath}/ldsc.py \
--h2 ${sumstat} \
--ref-ld-chr ./1000G.EUR.QC. \
--w-ld-chr ./weights. \
--overlap-annot \
--frqfile-chr ./1000G.EUR.QC. \
--out  ./`basename ${sumstat} .sumstats.gz`

