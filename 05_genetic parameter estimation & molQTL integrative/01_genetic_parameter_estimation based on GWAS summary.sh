
# download the LD Score Regression (LDSC) Version 1.0.1
# wget https://github.com/bulik/ldsc.git

bfile=PGRP_genotype
ldscore_OutFile=PGRP_l2
GWASSumFile=M_ADG.txt # including snpid chr bp a1 a2 beta pval
GWAS_N=10000
Sumstat_OutFile=M_ADG.sumstats
h2_OutFile=M_ADG.heritability
# calculate genetic correlation
Trait1=M_ADG.sumstats
Trait2=M_BFT.sumstats
rg_OutFile=M_ADG.M_BFT.rg

# step 1: calculate the ld score for each snp
python2 ./ldsc.py \
--bfile ${bfile} \
--l2 \
--ld-wind-kb 1000 \
--out ${outfile}

# step 2: pre-processing the GWAS summary statistics
python2 ./munge_sumstats.py \
--sumstats ${GWASSumFile} \ 
--N ${GWAS_N} \
--out ${Sumstat_OutFile}

# step 3: estimate the heritability
python2 ./ldsc.py \
--h2 ${Sumstat_OutFile} \
--ref-ld ${ldscore_OutFile} \
--w-ld ${ldscore_OutFile} \
--out ${h2_OutFile}

# step 4: estimate the genetic correlation
python2 ./ldsc.py \
--rg ${Trait1} ${Trait2} \
--ref-ld ${ldscore_OutFile} \
--w-ld ${ldscore_OutFile} \
--out ${rg_OutFile}
