# 1. parameter/file built for METAL
rule GWAS_GCTA:
    input:
        gwasFilelist=("./a.gwas" "./b.gwas") # A list of GWAS summary path

    output:
        output="meta.txt"

    shell:
        '''
        echo \
        "MARKER SNP
        ALLELE A1 A2
        FREQ AF1
        PVAL P
        EFFECT BETA
        STDERR SE

        SCHEME STDERR
        WEIGHT N
        AVERAGEFREQ ON
        MINMAXFREQ ON
        VERBOSE OFF
        GENOMICCONTROL ON
        " > ${output}

        for i in ${gwasFilelist[*]}
        do
            echo PROCESS ${i} >> ${output}
        done

        echo >> ${output}
        echo ANALYZE HETEROGENEITY >> ${output}
        '''

# 2. RUN METAL
    input:
        MetalPath="./metal" #METAL software executed file path
        input="meta.txt" # parameter file path

    output:
        output="meta.log" # output log file path

    shell:
        '''
        ${MetalPath} ${input} > ${output} 2>&1
        '''