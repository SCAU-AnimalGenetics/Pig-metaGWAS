# 1. GREML for GCTA in quality tarits
rule GREML_GCTA:
    input:
        GctaPath="./gcta64", # GCTA software executed file path
        GrmFilePath="./", # genomic-based relationship matrix prefix in GCTA
        PheFilePath="./phenotype.txt", # Phenotye record file
        covarfilePath="./covar.txt", # Covariate file
        QCovarFilePath="./qcovar.txt", # Quantitative covariate file
        ThreadsNum="20", # Specify the number of threads
    output:
        OutPrefix="./greml"
    shell:
        '''
        ${GctaPath} \
        --reml \
        --grm ${GrmFilePath} \
        --pheno ${PheFilePath} \
        --covar ${covarfilePath} \
        --qcovar ${QCovarFilePath} \
        --out ${OutPrefix} \
        --threads ${ThreadsNum}
        '''

# 2. GREML for GCTA in binary traits
rule GREML_GCTA:
    input:
        GctaPath="./gcta64", # GCTA software executed file path
        GrmFilePath="./", # genomic-based relationship matrix prefix in GCTA
        PheFilePath="./phenotype.txt", # Phenotye record file
        covarfilePath="./covar.txt", # Covariate file
        QCovarFilePath="./qcovar.txt", # Quantitative covariate file
        ThreadsNum="20", # Specify the number of threads
    output:
        OutPrefix="./greml"
    shell:
        '''
        ${GctaPath} \
        --reml \
        --grm ${GrmFilePath} \
        --pheno ${PheFilePath} \
        --covar ${covarfileORweightfilePath} \
        --qcovar ${QCovarFilePath} \
        --prevalence 0.01 \
        --out ${OutPrefix} --threads ${ThreadsNum}
        '''