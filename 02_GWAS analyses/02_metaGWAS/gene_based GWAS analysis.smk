# Gene-based analysis of GWAS 
rule GREML_GCTA:
    input:
        MagmaPath="./magma", # MAGMA software executed file path
        ReferencePlinkbfilePath="./", # reference panel raw gentype 
        GeneAnnotationFilePath="./", # Gene annotation
        GWASFilePath="./", # GWAS summary file
        variant_id_colname="variant_id"
        pvalue_colname="pvalue"
        sample_size_colname="sample_size"
    output:
        OutFile="./"
    shell:
        '''
        ${MagmaPath} \
        --bfile ${ReferencePlinkbfilePath} \
        --gene-annot ${GeneAnnotationFilePath} \
        --pval ${GWASFilePath} \
        use=${variant_id_colname},${pvalue_colname} \
        ncol=${sample_size_colname} \
        --out ${OutFile}
        '''