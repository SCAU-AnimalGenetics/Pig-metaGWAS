# 1. GWAS for quantitative tarits in GCTA
rule GRM_GCTA
    input:
        GctaPath="./gcta64", # GCTA software executed file path
        GenoFilePath="./", # Prefixes for genotype data in PLINK binary format
    output:
        OutPrefix="./Grm"
    shell：
        '''
        ${GctaPath}\
        --bfile ./${GenoFilePath} \
        --make-grm-alg 1 \
        --make-grm \
        --out ./${OutPrefix}

        ${Gcta_path} \
        --grm ./${OutPrefix} \
        --make-bK-sparse 0.05 \
        --out ./${OutPrefix}_fastGWAS
        '''

rule GWAS_GCTA:
    input:
        GctaPath="./gcta64", # GCTA software executed file path
        GrmFilePath="./", # genomic-based relationship matrix prefix in GCTA
        PheFilePath="./phenotype.txt", # Phenotye record file
        covarfilePath="./covar.txt", # Covariate file
        QCovarFilePath="./qcovar.txt", # Quantitative covariate file
        ThreadsNum="20", # Specify the number of threads
    output:
        OutPrefix="./gwas"
    shell:
        '''
        ${GctaPath} \
            --grm-sparse ${GrmFilePath} \
            --bfile ${GenFileprefix} \
            --pheno ${PheFilePath} \
            --mpheno 1 \
            --fastGWA-mlm \
            --covar ${covarfilePath} \
            --qcovar ${QCovarFilePath} \
            --out ${OutPrefix} \
            --threads ${ThreadsNum} \
            --seed 555
        '''

# 2. GWAS for binary traits in GCTA
rule GRM_GCTA
    input:
        GctaPath="./gcta64", # GCTA software executed file path
        GenoFilePath="./", # Prefixes for genotype data in PLINK binary format
    output:
        OutPrefix="./Grm"
    shell：
        '''
        ${GctaPath}\
        --bfile ./${GenoFilePath} \
        --make-grm-alg 1 \
        --make-grm \
        --out ./${OutPrefix}

        ${Gcta_path} \
        --grm ./${OutPrefix} \
        --make-bK-sparse 0.05 \
        --out ./${OutPrefix}_fastGWAS
        '''

rule GWAS_GCTA:
    input:
        GctaPath="./gcta64", # GCTA software executed file path
        GrmFilePath="./", # genomic-based relationship matrix prefix in GCTA
        PheFilePath="./phenotype.txt", # Phenotye record file
        covarfileORweightfilePath="./covar.txt", # Covariate file
        QCovarFilePath="./qcovar.txt", # Quantitative covariate file
        ThreadsNum="20", # Specify the number of threads
    output:
        OutPrefix="./gwas"
    shell:
        '''
        ${GctaPath} \
            --grm-sparse ${GrmFilePath} \
            --bfile ${GenFileprefix} \
            --pheno ${PheFilePath} \
            --mpheno 1 \
            --fastGWA-mlm-binary \
            --covar ${covarfileORweightfilePath} \
            --qcovar ${QCovarFilePath} \
            --out ${OutPrefix} \
            --threads ${ThreadsNum} \
            --seed 555
        '''




# 3. GWAS for MT traits in MMAP 
rule GRM_MMAP
    input:
        MmapPath="./mmap.2021_08_19_22_30.intel", # MMAP software executed file path
        GenoFilePath="./", # Prefixes for genotype data in PLINK binary format
        ThreadsNum="20", # Specify the number of threads
        ${ChrNum}="1", # Specify the number of Chromosome
    output:
        OutPrefix="./Grm"
    shell：
        '''
        ${Mmap_path} \
        --plink_bfile2mmap \
        --plink_bfile ./${GenoFilePath} \
        --binary_output_prefix ./${GenoFilePath}.MS
        
        ${Mmap_path} \
        --transpose_binary_genotype_file \
        --binary_input_filename ./${GenoFilePath}.MS.bin \
        --binary_output_filename ./${GenoFilePath}.SM.bin

        ${Mmap_path} \
        --write_binary_gmatrix_file \
        --write_matrix_counts \
        --binary_genotype_filename ./${GenoFilePath}.SM.bin \
        --group_size 100 \
        --single_pedigree \
        --num_mkl_threads 1 \
        --chromosome ${ChrNum}  \
        --binary_output_filename ./${GenoFilePath}_chr${ChrNum}.SM.grm.bin 
        
        ${Mmap_path} \
        --combine_binary_matrix_files $(ls ./${GenoFilePath}_chr*.SM.grm.bin) \
        --binary_output_filename ./${GenoFilePath}.SM.grm.bin \
        --num_mkl_threads ${ThreadsNum}
        '''

rule GWAS_MMAP:
    input:
        MmapPath="./mmap.2021_08_19_22_30.intel", # MMAP software executed file path
        GenFileprefix="./", # MMAP genotype MS and SM format file path
        PheFilePath="./phenotype.txt", # Phenotye record file
        GrmFilePath="./covar.txt", # genomic-based relationship matrix prefix in GCTA
        covarfileORweightfilePath="./qcovar.txt" # Covariate file
        chr="", # Chromosome selected to be analyzed
    output:
        OutPrefix="gwas"
    shell:
        '''
        ${MmapPath} \
            --ped ${GenFileprefix}.MS.ped.csv \
            --phenotype_filename ${PheFilePath} \
            --read_binary_covariance_file ${GrmFilePath} \
            --error_variance_diagonal_filename ${covarfileORweightfilePath} \
            --file_suffix ${OutPrefix}.mmap \
            --binary_genotype_filename ${GenFileprefix}.MS.bin \
            --trait MT \
            --all_output \
            --num_mkl_threads 1 \
            --chromosome ${chr}
        '''