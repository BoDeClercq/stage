library(readr)
gene_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s <- read_delim("C:/Users/bodec/Downloads/plasma_LR/plasma_LR/gene_counts-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s.tsv", 
                                                                                     delim = "\t", escape_double = FALSE, 
                                                                                     col_types = cols(begin = col_number(), 
                                                                                                      end = col_number(), `counts-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number(), 
                                                                                                      `nicounts-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number()), 
                                                                                     trim_ws = TRUE)


#dataset opmaken 
gene_d11345 <- mutate(`gene_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s`,Sample = 'd11345') 
View(gene_d11345)

colnames(gene_d11345)
names(gene_d11345)[names(gene_d11345) == 'counts-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s'] <- 'counts' 
names(gene_d11345)[names(gene_d11345) == 'nicounts-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s'] <- 'nicounts' 
colnames(gene_d11345)




#20 meest voorkomende genen 
gene_d11345_sorted <- gene_d11345[order(gene_d11345$nicounts,decreasing = TRUE),]
gene_d113445_20 <- head(gene_d11345_sorted, 20)
View(gene_d113445_20)
