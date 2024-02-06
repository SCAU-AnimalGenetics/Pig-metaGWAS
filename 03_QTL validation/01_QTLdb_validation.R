# load package
library(parallel)
library(ggplot2)
library(tidyverse)

# set work path
wkdir = "./QTLdb_annotaion"
setwd(wkdir)

QTLdb_FullName = "Average daily gain"
RegionType = "SelfDefine"
ExpandWidth = 1000000 # define the window size
type = paste0(RegionType, "_", ExpandWidth)

filename = "D_ADG_LeadSNP.txt" # read QTL region
QTLdb_file = "Animal_QTLdb_release46_pigSS11.bed"

##=========== preprocess the QTL information =============
QTLdb = as.matrix(read.table(QTLdb_file, skip = 20, fill = T, sep = "\t"))
# recode chr information
QTLdb[, 1] = gsub("Chr.", "", QTLdb[, 1])
# rm the unknown start and end
QTLdb = QTLdb[QTLdb[, 1] != "X", ]
QTLdb <- QTLdb[!is.na(QTLdb[, 2]), ]
QTLdb <- QTLdb[!is.na(QTLdb[, 3]), ]
# remove the QTL region with length more than 1Mb
QTLdb = QTLdb[as.numeric(QTLdb[, 3]) - as.numeric(QTLdb[, 2]) >= 0 
              & as.numeric(QTLdb[, 3]) - as.numeric(QTLdb[, 2]) <= 1000000, ]
# get the trait name
QTLdb_TraitName = gsub(" QTL .*", "", QTLdb[, 4])
QTLdb_TraitName = gsub("[{}()]", "_", QTLdb_TraitName)

#============== read QTL region ======================
leadinfo = read.table(filename, header = T)
trait_name = gsub("_LeadSNP.txt", "", filename)

leadinfo = data.frame(leadinfo, NA)
colnames(leadinfo) = c(colnames(leadinfo)[1:(ncol(leadinfo) - 1)], type)

#------------------ QTL region define -------------------------
data = data.frame(SNPID = leadinfo$SNP, CHR = leadinfo$CHR, BP = leadinfo$BP, 
                  START = (leadinfo$BP - ExpandWidth), END = (leadinfo$BP + ExpandWidth))

##------------------ define the core region -----------------------------
qtl = data.frame(as.numeric(data$CHR), as.numeric(data$START), as.numeric(data$END))
colnames(qtl) = c("Chr", "Upstream", "Downstream")

# define the QTL region overlapped with the QTLdb or not.
QTLmatch = function(qtldb, qtl){
  # if(qtl[1] < qtldb[1] & qtl[2] < qtldb[1] | qtl[1] > qtldb[2] & qtl[2] > qtldb[2]){
  if(qtl[2] < qtldb[1] | qtl[1] > qtldb[2]){
    return(F)
  } else {
    return(T)
  }
}

# --------- QTL region annotation
# have trait name and have more than 1 QTL region in QTLdb
if(length(grep(QTLdb_FullName, QTLdb_TraitName)) > 1){
  qtldb_trait = data.frame(QTLdb[grep(QTLdb_FullName, QTLdb_TraitName), ])
  # qtl and qtldb match
  for(j in 1:nrow(qtl)){
    # j = 1 # nrow(qtl)
    qtl_j = as.numeric(qtl[j, ])
    qtldb_trait_j = qtldb_trait[qtldb_trait$V1 == qtl_j[1], ]
    if(nrow(qtldb_trait_j) > 0){
      overlap = c()
      for(k in 1:nrow(qtldb_trait_j)){
        # k = 1 # nrow(qtldb_trait_j)
        qtldb_trait_j_k = as.numeric(qtldb_trait_j[k, 2:3])
        overlap = c(overlap, QTLmatch(qtldb_trait_j_k, qtl_j[2:3]))
        if(sum(overlap) > 0){
          leadinfo[j, ncol(leadinfo)] = "TRUE"
        } else {
          leadinfo[j, ncol(leadinfo)] = "FALSE"
        }
      }
    } else {
      leadinfo[j, ncol(leadinfo)] = "FALSE"
    }
  }
}
# have trait name but only have 1 QTL region in QTLdb
if(length(grep(QTLdb_FullName, QTLdb_TraitName)) == 1){
  qtldb_trait = data.frame(matrix(QTLdb[grep(QTLdb_FullName, QTLdb_TraitName), ], ncol = 12))
  for(j in 1:nrow(qtl)){
    # j = 19 # nrow(qtl)
    qtl_j = as.numeric(qtl[j, ])
    qtldb_trait_j = qtldb_trait[qtldb_trait$X1 == qtl_j[1], ]
    if(nrow(qtldb_trait_j) > 0){
      # overlap = c()
      qtldb_trait_j_k = as.numeric(qtldb_trait_j[1, 2:3])
      
      if(sum(QTLmatch(qtldb_trait_j_k, qtl_j[2:3])) > 0){
        leadinfo[j, ncol(leadinfo)] = "TRUE"
      } else {
        leadinfo[j, ncol(leadinfo)] = "FALSE"
      }
    } else {
      leadinfo[j, ncol(leadinfo)] = "FALSE"
    }
  }
}
# have trait name but do not have any QTL region in QTLdb
if(length(grep(QTLdb_FullName, QTLdb_TraitName)) == 0){
  qtldb_trait = data.frame(matrix(QTLdb[grep(QTLdb_FullName, QTLdb_TraitName), ], ncol = 12))
  for(j in 1:nrow(qtl)){
    # j = 19 # nrow(qtl)
    leadinfo[j, ncol(leadinfo)] = "FALSE"
  }
}
# save the result
write.table(leadinfo, filename, row.names = F, col.names = T, quote = F)
