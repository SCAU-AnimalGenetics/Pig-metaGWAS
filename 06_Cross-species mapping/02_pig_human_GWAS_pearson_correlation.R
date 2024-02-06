library(data.table)

# set work path
wkdir = ""
setwd(wkdir)
# match the homologous variants before calculating correlation
filename = "D_ADG&AD_Kunkle_2019.txt"

# read data
data = fread(filename)

# calculate the correlation of Zscore
z_cor = cor.test(data$pig_Zscore, data$human_Zscore) 
# calculate the correlation of absolute Zscore
abs_z_cor = cor.test(abs(data$pig_Zscore), abs(data$human_Zscore))
# calculate the correlation of P value
p_cor = cor.test(data$pig_pvalue, data$human_pvalue)
res = data.frame(trait_pair = gsub(".txt", "", filename), 
                 z_cor = z_cor$estimate, z_cor_p = z_cor$p.value, 
                 abs_z_cor = abs_z_cor$estimate, abs_z_cor_p = abs_z_cor$p.value, 
                 p_cor = p_cor$estimate, p_cor_p = p_cor$p.value)
