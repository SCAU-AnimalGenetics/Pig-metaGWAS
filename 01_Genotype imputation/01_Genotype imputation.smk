# 1. Genotype imputation with beagle
  input:
    Chr=1  # Chromosome 1~18
  params:
    threads: 22 # Maximum number of threads for beagle program
    MemoryUSE=120 # Number of memory for Beagle program
  output:
    out="Target_file.Chr${Chr}.CF"
 # shell:
  ## Allele correction of imputed data to match the reference panel
  java -Xmx${MemoryUSE}g -XX:ParallelGCThreads=${ThreadUSE} -jar conform-gt.24May16.cee.jar  \
  ref=PGRP.vcf.gz  \
  gt=Target_file.Chr${Chr}.vcf.gz  \
  chrom=${Chr}  \
  out=${out}  \
  match=POS    
  
  ## Imputation
  java -Xmx${MemoryUSE}g -jar beagle.17Apr21.dc6.jar  \
  ref=PGRP.vcf.gz  \
  gt=${out}.vcf.gz  \
  out=${out}.Impute  \
  impute=true  \
  ne=1000  \
  seed=9823  \
  nthreads=${ThreadUSE}
  

# 2. Genotype quality control 
  input:
    QCStandard="DR2>=0.80&MAF>=0.01" # Keep snps with DR2>=0.80 and minor allele frequency>=0.01
    Chr=1  # Chromosome 1~18
  output:  
    out="Target_file.Chr${Chr}.CF.Impute.QC"
  # shell: 
  ## Indexing
  bcftools index -t Target_file.Chr${Chr}.CF.Impute.vcf.gz
  
  ## DR2 and MAF quality control
  bcftools filter -i "${QCStandard}" -Oz Target_file.Chr${Chr}.CF.Impute.vcf.gz > ${out}.vcf.gz
  