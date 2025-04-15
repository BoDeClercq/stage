library(readr)
gene_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s <- read_delim("C:/Users/bodec/Downloads/plasma_LR/plasma_LR/gene_counts-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s.tsv", 
                                                                                     delim = "\t", escape_double = FALSE, 
                                                                                     col_types = cols(begin = col_number(), 
                                                                                                      end = col_number(), `counts-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number(), 
                                                                                                      `nicounts-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number()), 
                                                                                     trim_ws = TRUE)


#dataset opmaken 
gene_d11358 <- mutate(`gene_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s`,Sample = 'd11358') 
View(gene_d11358)

colnames(gene_d11358)
names(gene_d11358)[names(gene_d11358) == 'counts-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s'] <- 'counts' 
names(gene_d11358)[names(gene_d11358) == 'nicounts-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s'] <- 'nicounts' 
colnames(gene_d11358)



#20 meest voorkomende genen 
gene_d11358_sorted <- gene_d11358[order(gene_d11358$nicounts,decreasing = TRUE),]
gene_d113458_20 <- head(gene_d11358_sorted, 20)
View(gene_d113458_20)
