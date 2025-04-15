library(readr)
gene_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s <- read_delim("C:/Users/bodec/Downloads/plasma_LR/plasma_LR/gene_counts-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s.tsv", 
                                                                                     delim = "\t", escape_double = FALSE, 
                                                                                     col_types = cols(begin = col_number(), 
                                                                                                      end = col_number(), `counts-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number(),                                                                                         `nicounts-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number()), 
                                                                                     trim_ws = TRUE)



#packages + libraries 
install.packages('tidyverse')
library(tidyverse)
library(dplyr)
library(ggplot2)

if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('biomaRt')

library(biomaRt)
ensemblgrch38 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes_ens <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol", "gene_biotype"), mart=ensemblgrch38)


#dataset opmaken 
gene_d11344 <- mutate(gene_d11344,Sample = 'd11344') 
View(gene_d11344)

colnames(gene_d11344)
names(gene_d11344)[names(gene_d11344) == 'counts-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s'] <- 'counts' 
names(gene_d11344)[names(gene_d11344) == 'nicounts-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s'] <- 'nicounts' 
colnames(gene_d11344)


# 1 dataset van 3 samples 
gene_combined <- bind_rows(gene_d11344, gene_d11345, gene_d11358)
View(gene_combined)


#cumulatieve plot genes + bereking nicounts bij 50%
gene_combined_filtered <- subset(gene_combined, nicounts !=0)
View(gene_combined_filtered)


ggplot(gene_combined_filtered, aes(x = nicounts, color = Sample)) +
  stat_ecdf() + 
  labs (x = 'Non-intronic counts', y = 'Cumulative fraction') +
  scale_x_log10() +
  scale_color_manual(values = c('#00C1A7', '#C39BD3', '#A3C1DA')) +
  theme_minimal()

gene_combined_filtered <- gene_combined_filtered %>%
  group_by(Sample) %>%
  mutate(ecdf_value = ecdf(nicounts)(nicounts)) %>%
  ungroup()
threshold_50 <- aggregate(nicounts ~ Sample, gene_combined_filtered[gene_combined_filtered$ecdf_value >= 0.5, ], function(x) min(x))
print(threshold_50)


#20 meest voorkomende genen per sample 
gene_d11344_sorted <- gene_d11344[order(gene_d11344$nicounts,decreasing = TRUE),]
gene_d11344_20 <- head(gene_d11344_sorted, 20)
View(gene_d11344_20)

gene_20_combined <- bind_rows(gene_d113445_20, gene_d11344_20, gene_d113458_20)
View(gene_20_combined)

ggplot(gene_20_combined, aes(x = gene, y = nicounts, fill = Sample)) +
  geom_bar(stat = 'identity', position = 'stack' ) + 
  labs (x = 'Gene', y = 'Non-intronic counts') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c('#00C1A7', '#A3C1DA', '#C39BD3'))
 

#biotypes 

    #biotype kolom toevoegen aan database
gene_combined <- mutate(gene_combined, geneid = gsub("\\..*", "", geneid))
gene_combined <- mutate(gene_combined, gene_biotype = genes_ens$gene_biotype[match(gene_combined$geneid, genes_ens$ensembl_gene_id)]) 
View(gene_combined)
gene_combined_nicounts_excluded5 <- subset(gene_combined, nicounts > 5)
View(gene_combined_nicounts_excluded5)

#plot van biotype verdeling (kort)
gene_combined_biotype <- gene_combined_nicounts_excluded5 %>%
    group_by(Sample, gene_biotype) %>%
    summarise(Biotype_count = sum(nicounts), .groups = "drop") %>%  
    group_by(Sample) %>%
    mutate(Biotype_percent = Biotype_count / sum(Biotype_count) * 100) %>%
    mutate(gene_biotype = ifelse(Biotype_percent < 1, "less than 1%", gene_biotype))
View(gene_combined_biotype) 


ggplot(gene_combined_biotype, aes(x = Sample, y = Biotype_percent, fill = gene_biotype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  labs(x = NULL, y = "Fraction of reads", fill = 'Biotype') +
  theme_minimal() + 
  scale_fill_manual(values = c('#E9967A','#00C1A7','#C39BD3', 'grey', '#A3C1DA'))
  
#plot van biotype verdeling (lang, gemiddelde)
gene_combined_biotype_all <- gene_combined_acountss_excluded5 %>%
  group_by(Sample, gene_biotype) %>%
  summarise(Biotype_count = sum(nicounts), .groups = "drop") %>%  # Som van counts per biotype
  group_by(Sample) %>%
  mutate(Biotype_percent = Biotype_count / sum(Biotype_count) * 100) 

gene_combined_biotype_avg <- gene_combined_biotype_all %>%
  group_by(gene_biotype) %>%
  summarise(Average_Biotype_count = mean(Biotype_count))
View(gene_combined_biotype_avg)

ggplot(gene_combined_biotype_avg, aes(x = gene_biotype, y = Average_Biotype_count)) +
  geom_bar(stat = "identity", position = "dodge", fill = 'grey', color= 'black') +  
  labs(x = 'Biotype', y = "Average number of non-intronic counts") +  
  theme_minimal() +
  scale_y_log10() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    panel.spacing = unit(1.5, "lines"),
    legend.text = element_text(size = 10)
  )

