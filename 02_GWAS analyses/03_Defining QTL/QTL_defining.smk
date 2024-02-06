# 1. QTL defining
rule QTL defining:
    input:
        gwasfilepath="./gcta.gwas", # GWAS summary path

    output:
        output="./gcta_qtl.txt"

    Rscript:
        '''
        library(data.table)
        library(stringr)

        Distance = 1000000
        HaonanFinemap = function(data,Distance,P_threshold){
        
        lead_snp = data.frame()
        data_sig = data[which(data$P<P_threshold),]
        for (chr in unique(data_sig$CHR)) {
            
            data_sig_chr = data_sig[data_sig$CHR==chr,]
            
            while(dim(data_sig_chr)[1]>0){
            
            top_snp = data_sig_chr[which(data_sig_chr$P == min(data_sig_chr$P)),][1,]
            top_snp$QTLR = top_snp$QTLL = top_snp$BP
            top_snp_QTL = data_sig_chr[which(data_sig_chr$BP >= top_snp$QTLL-550000 & data_sig_chr$BP <= top_snp$QTLR+550000 & data_sig_chr$P >= top_snp$P & data_sig_chr$P <= top_snp$P+1e-4),]
            if(top_snp$QTLL - min(top_snp_QTL$BP) >= 500000){top_snp$QTLL = top_snp$QTLL-500000}else{top_snp$QTLL = min(top_snp_QTL$BP)}
            if(max(top_snp_QTL$BP) - top_snp$QTLR >= 500000){top_snp$QTLR = top_snp$QTLR+500000}else{top_snp$QTLR = max(top_snp_QTL$BP)}
            
            top_snp$BPR = top_snp$BPL = top_snp$BP
            top_snp$BPL[top_snp$BPL<0] = 0
            
            top_snp_left = data_sig_chr[which(data_sig_chr$BP >= top_snp$BPL-Distance & data_sig_chr$BP <= top_snp$BPL),]
            while(dim(top_snp_left)[1]>0) {
                top_snp$BPL = min(top_snp_left$BP)
                data_sig_chr = data_sig_chr[-which(data_sig_chr$SNP%in%top_snp_left$SNP),]
                top_snp_left = data_sig_chr[which(data_sig_chr$BP >= top_snp$BPL-Distance & data_sig_chr$BP <= top_snp$BPL),]
            }
            
            top_snp_right = data_sig_chr[which(data_sig_chr$BP >= top_snp$BPR & data_sig_chr$BP <= top_snp$BPR+Distance),]
            while(dim(top_snp_right)[1]>0) {
                top_snp$BPR = max(top_snp_right$BP)
                data_sig_chr = data_sig_chr[-which(data_sig_chr$SNP%in%top_snp_right$SNP),]
                top_snp_right = data_sig_chr[which(data_sig_chr$BP >= top_snp$BPR & data_sig_chr$BP <= top_snp$BPR+Distance),]
            }
            
            lead_snp = rbind(lead_snp,top_snp)
            
            }
            
        }
        return(lead_snp)
        
        }

        FilePath="${gwasfilepath}"
        data = fread(FilePath,data.table = F,nThread=1)[,c("SNP","CHR","POS","P")]
        colnames(data) = c("SNP","CHR","BP","P")
        data[,2:4] = apply(data[,2:4],2,as.numeric)

        lead_snp = HaonanFinemap(data,Distance,5e-8)

        if(dim(lead_snp)[1]>0){
        OutPath="${output}"
        write.table(lead_snp,OutPath,sep="\t",row.names=F,quote=F)
        }
        '''