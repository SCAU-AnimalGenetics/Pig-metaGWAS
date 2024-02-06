# ldak v5.0
bfile=PGRP_genotype
tagging_OutFile=BLD-LDAK_model
annot_Prefix=Annot # SNP list with multiple type of annotation：① significant cis-molQTLs from five molecular phenotypes in 34 tissues; ② the significant variants on the tissue-sharing/specific egenes for PCG
Summary_InputFile=D_ADG..summarise # including Predictor A1 A2 Direction Stat n
SumHer_OutFile=D_ADG

# step 1: calculate the tagging file
./ldak5.linux.fast --calc-tagging ${tagging_OutFile} \
--bfile ${bfile} \
--ignore-weights YES \
--max-threads 23 \
--power 0 \
--window-kb 200 \
--annotation-number 3 \
--annotation-prefix ${annot_Prefix}

# step 2: implement the heritability enrichment
./ldak5.linux.fast --sum-hers ${SumHer_OutFile} \
--summary ${Summary_InputFile} \
--tagfile ${tagging_OutFile}.tagging \
--max-threads 23 \
--check-sums NO
