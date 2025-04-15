library(readr)
isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s <- read_delim("C:/Users/bodec/Downloads/plasma_LR/plasma_LR/isoform_counts-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s.tsv", 
                                                                                        delim = "\t", escape_double = FALSE, 
                                                                                        col_types = cols(begin = col_number(), 
                                                                                                         end = col_number(), exonStarts = col_character(), 
                                                                                                         exonEnds = col_character(), cdsStart = col_number(), 
                                                                                                         cdsEnd = col_number(), size = col_number(), 
                                                                                                         `counts_iqall-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number(), 
                                                                                                         `counts_weighed-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number(), 
                                                                                                         `counts_unique-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number(), 
                                                                                                         `counts_strict-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number(), 
                                                                                                         `counts_aweighed-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number(), 
                                                                                                         `counts_aunique-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number(), 
                                                                                                         `counts_astrict-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s` = col_number()), 
                                                                                        trim_ws = TRUE, skip = 27)


#packages + libraries 
install.packages('tidyverse')
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

install.packages("ggrepel")
library(ggrepel)

if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('biomaRt')

library(biomaRt)
ensemblgrch38 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes_ens <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol", "gene_biotype"), mart=ensemblgrch38)



#dataset ordenen + extra nuttige kolommen toevoegen
colnames(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s)
names(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s) == 'counts_aweighed-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s'] <- 'acountsw' 
names(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s) == 'counts_aunique-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s'] <- 'acountsu'
names(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s) == 'counts_astrict-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s'] <- 'acountss'
colnames(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s)

isoform_d11344_sorted <- isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s[order(isoform_counts_isoquant_joint_sminimap2_splice_d11344_v7_6_8_210332_hg38s$acountsw,decreasing = TRUE),]
isoform_d11344_sorted <- mutate(isoform_d11344_sorted,Sample = 'd11344') 
View(isoform_d11344_sorted)



#poly(A) data 
d11344_som_weighted <- sum(isoform_d11344_sorted$`counts_weighed-isoquant_joint-sminimap2_splice-d11344_v7.6.8_210332_hg38s`, na.rm = TRUE)
head(d11344_som_weighted)
d11344_som_aweigthed <- sum(isoform_d11344_sorted$acountsw, na.rm = TRUE)
head(d11344_som_aweigthed)
d11344_fraction <- d11344_som_aweigthed/d11344_som_weighted
head(d11344_fraction)



#data fractie van verschillende counts
d11344_som_aweigthed <- sum(isoform_d11344_sorted$acountsw, na.rm = TRUE)
head(d11344_som_aweigthed)
d11344_som_aunique <- sum(isoform_d11344_sorted$acountsu, na.rm = TRUE)
head(d11344_som_aunique)
d11344_som_astrict <- sum(isoform_d11344_sorted$acountss, na.rm = TRUE)
head(d11344_som_astrict)

d11344_counts <- data.frame(
  category = c("Unique", "Intact", "Weighted"),
  value = c(d11344_som_aunique, d11344_som_astrict, d11344_som_aweigthed),
  ID = c('d11344')
)
         
        #combinatie van de 3 stalen 
combined_counts <- bind_rows(d11344_counts, d11345_counts, d11358_counts)
View(combined_counts)


ggplot() + 
  geom_bar(data = d11344_counts[d11344_counts$category == "Weighted", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.7, alpha = 0.9) +  
  geom_bar(data = d11344_counts[d11344_counts$category == "Unique", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.6, alpha = 0.9) +  
  geom_bar(data = d11344_counts[d11344_counts$category == "Intact", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.5, alpha = 0.9) +
  geom_bar(data = d11345_counts[d11345_counts$category == "Weighted", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.7, alpha = 0.9) +  
  geom_bar(data = d11345_counts[d11345_counts$category == "Unique", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.6, alpha = 0.9) + 
  geom_bar(data = d11345_counts[d11345_counts$category == "Intact", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.5, alpha = 0.9) +  
  geom_bar(data = d11358_counts[d11358_counts$category == "Weighted", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.7, alpha = 0.9) +  
  geom_bar(data = d11358_counts[d11358_counts$category == "Unique", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.6, alpha = 0.9) +  
  geom_bar(data = d11358_counts[d11358_counts$category == "Intact", ], 
           aes(x = "Total", y = value, fill = category), 
           stat = "identity", position = "identity", width = 0.5, alpha = 0.9) +  
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs (y = 'Poly(A) counts') +
  scale_fill_manual(values = c('#A3C1DA','#C39BD3','#00C1A7'))+
  labs(fill = "Category") +
  facet_wrap(~ID, strip.position = "bottom") +
  guides(fill = guide_legend(title = NULL)) 


#transcript data 
isoform_d11344_transcripts <- summarise(group_by(isoform_d11344_sorted, gene), 'aantal transcripten' = n())
isoform_d11344_transcripts <- mutate(isoform_d11344_transcripts, Sample = 'd111344')
isoform_d11344_transcripts <- isoform_d11344_transcripts[order(isoform_d11344_transcripts$'aantal transcripten',decreasing = TRUE),]
View(isoform_d11344_transcripts)

isoform_combined_transcripts <- bind_rows(isoform_d11344_transcripts, isoform_d11345_transcripts, isoform_d11358_transcripts)
View(isoform_combined_transcripts)

        #cumulatieve plot van transcripten + waarden
ggplot(isoform_combined_transcripts, aes(x = `aantal transcripten`, color = Sample)) +
  stat_ecdf() + 
  scale_x_log10()+
  scale_color_manual(values = c('#A3C1DA','#C39BD3','#00C1A7')) +
  labs(y = 'Cumulative fraction', x = 'Number of transcripts') +
  theme_minimal()
  
isoform_combined_transcripts <- isoform_combined_transcripts %>%
  group_by(Sample) %>%
  mutate(ecdf_value = ecdf(`aantal transcripten`)(`aantal transcripten`)) %>%
  ungroup()
threshold_50 <- aggregate(`aantal transcripten` ~ Sample, isoform_combined_transcripts[isoform_combined_transcripts$ecdf_value >= 0.5, ], function(x) min(x))
print(threshold_50)


#cumulatieve plot van counts gescheiden per categorie 
isoform_d11344_acountss_excluded0 <- subset(isoform_d11344_sorted, acountss !=0)

isoform_combined_acountss_excluded0 <- bind_rows(isoform_d11344_acountss_excluded0, isoform_d11345_acountss_excluded0, isoform_d11358_acountss_excluded0)
View(isoform_combined_acountss_excluded0)


ggplot(isoform_combined_acountss_excluded0_avg, aes(x = Average_category_count, color = category)) +
  stat_ecdf() + 
  scale_x_log10() +
  labs (x = 'Intact poly(A) counts', y = 'Cumulative fraction', color = 'Category') +
  scale_color_manual(values = c('darkgrey','#C39BD3', '#00C1A7', '#A3C1DA'))
  theme_minimal()

isoform_combined_acountss_excluded0_avg <- isoform_combined_acountss_excluded0_avg %>%
  group_by(category) %>%
  mutate(ecdf_value = ecdf(`Average_category_count`)(`Average_category_count`)) %>%
  ungroup()
threshold_90 <- aggregate(`Average_category_count` ~ category, isoform_combined_acountss_excluded0_avg[isoform_combined_acountss_excluded0_avg$ecdf_value >= 0.9, ], function(x) min(x))
print(threshold_90)


#novel data
isoform_d11344_filtered_known <- subset(isoform_d11344_sorted, subset = category == 'known')
isoform_d11344_filtered_known_excluded0 <- subset(isoform_d11344_filtered_known, acountss !=0)
isoform_d11344_filtered_novelincatalog <- subset(isoform_d11344_sorted, subset = category == 'novel_in_catalog')
isoform_d11344_filtered_novelnotincatalog <- subset(isoform_d11344_sorted, subset = category == 'novel_not_in_catalog')

       #t
t <- rbind(isoform_d11344_filtered_novelincatalog %>%
             dplyr::select(gene, exonStarts, exonEnds, acountss, Sample), 
           isoform_d11358_filtered_novelincatalog %>%
             dplyr::select(gene, exonStarts, exonEnds, acountss, Sample), 
           isoform_d11345_filtered_novelincatalog %>%
             dplyr::select(gene, exonStarts, exonEnds, acountss, Sample))


t <- t %>% 
  arrange(desc(acountss)) %>%   
  group_by(gene, exonEnds, exonStarts) %>%  
  mutate(n = n()) %>%   
  ungroup() %>%
  select(-Sample) %>%    
  distinct(gene, exonStarts, exonEnds, .keep_all = TRUE) %>%     
  arrange(desc(n))  


t <- mutate(t, category = 'novel in catalog')
View(t)

         #x
x <- rbind(isoform_d11344_filtered_novelnotincatalog %>%
             dplyr::select(gene, exonStarts, exonEnds, acountss, Sample), 
           isoform_d11358_filtered_novelnotincatalog %>%
             dplyr::select(gene, exonStarts, exonEnds, acountss, Sample), 
           isoform_d11345_filtered_novelnotincatalog %>%
             dplyr::select(gene, exonStarts, exonEnds, acountss, Sample))

x <- x %>% 
  arrange(desc(acountss)) %>%   
  group_by(gene, exonEnds, exonStarts) %>%   
  mutate(n = n()) %>%   
  ungroup() %>%
  select(-Sample) %>%   
  distinct(gene, exonStarts, exonEnds, .keep_all = TRUE) %>%    
  arrange(desc(n))  

x <- mutate(x, category = 'novel not in catalog')
View(x)

        #samenzetten van x en t grafiek 
xt <- bind_rows(x, t)
View(xt)

ggplot(xt, aes(x = n)) +  
  geom_bar(fill = "grey", color = "black") +
  scale_x_discrete(limits = c("1", "2")) + 
  facet_wrap(~category) +
  theme_minimal() +
  scale_y_log10() +
  labs(x = "Number of samples", y = "Number of transcripts")


t_subset <- subset(t, n==2)
x_subset <- subset(x, n==2)


            #novel in catalog 
isoform_combined_filtered_known_excluded0 <- bind_rows(isoform_d11344_filtered_known_excluded0, isoform_d11345_filtered_known_excluded0, isoform_d11358_filtered_known_excluded0)
quantile(isoform_combined_filtered_known_excluded0$acountss, probs = 0.8)

t_subset <- subset(t_subset, acountss > 26.2)
sum(t_subset$acountss > 26.2)

t_subset <- t_subset %>%
  mutate(exonStarts = strsplit(as.character(exonStarts), ","), 
         exonEnds = strsplit(as.character(exonEnds), ",")) %>%
  unnest(cols = c(exonStarts, exonEnds))

t_subset$exonStarts <- as.numeric(t_subset$exonStarts)
t_subset$exonEnds <- as.numeric(t_subset$exonEnds)
View(t_subset)

t_subset_SRSF3 <- subset(t_subset, gene == 'SRSF3')
t_subset_SRSF3 <- t_subset_SRSF3 %>%
  mutate(exon_number = row_number())
View(t_subset_SRSF3)

t_subset_METAP2 <- subset(t_subset, gene == 'METAP2')
t_subset_METAP2 <- t_subset_METAP2 %>%
  mutate(exon_number = row_number())
View(t_subset_METAP2)

t_subset_CNFN <- subset(t_subset, gene == 'CNFN')
t_subset_CNFN <- t_subset_CNFN %>%
  mutate(exon_number = row_number())
View(t_subset_CNFN)

t_subset_combined <- bind_rows(
  mutate(t_subset_SRSF3, gene = "SRSF3"),
  mutate(t_subset_METAP2, gene = "METAP2"),
  mutate(t_subset_CNFN, gene = "CNFN")
)

ggplot(t_subset_combined, aes(x = exonStarts, xend = exonEnds, y = 1, yend = 1, color = as.factor(exon_number))) +
  geom_segment(size = 1) +
  geom_text(aes(x = exonStarts, label = exonStarts), 
            size = 3, vjust = -2, hjust = -0, angle = 45) +  
  geom_text(aes(x = exonEnds, label = exonEnds), 
            size = 3, vjust = 3, hjust = 1, angle = 45) +  
  labs(x = "Genomic Position",y = NULL) +
  theme_minimal() +
  facet_wrap(~ gene, scales = "free_x", strip.position = "top", ncol = 1) +  
  theme(legend.position = "none", 
        strip.background = element_blank(),  
        strip.text = element_text(size = 12, face = "bold"),  
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(), 
        plot.margin = margin(1, 20, 1, 20)) 



            # novel not in catalog
isoform_combined_filtered_known_excluded0 <- bind_rows(isoform_d11344_filtered_known_excluded0, isoform_d11345_filtered_known_excluded0, isoform_d11358_filtered_known_excluded0)
quantile(isoform_combined_filtered_known_excluded0$acountss, probs = 0.8)

x_subset <- subset(x_subset, acountss > 26.2)
sum(x_subset$acountss > 26.2)

x_subset <- x_subset %>%
  mutate(exonStarts = strsplit(as.character(exonStarts), ","), 
         exonEnds = strsplit(as.character(exonEnds), ",")) %>%
  unnest(cols = c(exonStarts, exonEnds))

x_subset$exonStarts <- as.numeric(x_subset$exonStarts)
x_subset$exonEnds <- as.numeric(x_subset$exonEnds)
View(x_subset)

x_subset_TRBC2 <- subset(x_subset, gene == 'TRBC2')
x_subset_TRBC2 <- x_subset_TRBC2 %>%
  mutate(exon_number = row_number())
View(x_subset_TRBC2)

x_subset_IGLL5 <- subset(x_subset, gene == 'IGLL5')
x_subset_IGLL5 <- x_subset_IGLL5 %>%
  mutate(exon_number = row_number())
View(x_subset_IGLL5)

x_subset_RNASEK <- subset(x_subset, gene == 'RNASEK')
x_subset_RNASEK <- x_subset_RNASEK %>%
  mutate(exon_number = row_number())
View(x_subset_RNASEK)

x_subset_combined <- bind_rows(
  mutate(x_subset_TRBC2, gene = "TRBC2"),
  mutate(x_subset_IGLL5, gene = "IGLL5"),
  mutate(x_subset_RNASEK, gene = "RNASEK")
)

ggplot(x_subset_combined, aes(x = exonStarts, xend = exonEnds, y = 1, yend = 1, color = as.factor(exon_number))) +
  geom_segment(size = 1) +
  # Voeg een tekstlaag voor exonStarts toe (schuin boven de lijn aan het begin)
  geom_text(aes(x = exonStarts, label = exonStarts), 
            size = 3, vjust = -2, hjust = -0, angle = 45) +  
  # Voeg een tekstlaag voor exonEnds toe (schuin boven de lijn aan het einde)
  geom_text(aes(x = exonEnds, label = exonEnds), 
            size = 3, vjust = 3, hjust = 1, angle = 45) +  
  labs(x = "Genomic Position",
       y = NULL) +
  theme_minimal() +
  facet_wrap(~ gene, scales = "free_x", strip.position = "top", ncol = 1) +  # Facet per gen, met aparte x-as per facet
  theme(legend.position = "none", 
        strip.background = element_blank(),  # Verberg de achtergrond van de strip
        strip.text = element_text(size = 12, face = "bold"),  # Maak de gen-namen duidelijker
        axis.text.y = element_blank(),  # Verwijder de y-as tekst (de gene namen)
        axis.ticks.y = element_blank(), 
        plot.margin = margin(1, 20, 1, 1)) 


#data fragmentatie met gemiddelde 
isoform_d11344_excluded0 <- subset(isoform_d11344_sorted, acountsu !=0)
View(isoform_d11344_sorted)

isoform_0excluded_combined <- bind_rows(isoform_d11344_excluded0, isoform_d11358_excluded0, isoform_d11345_excluded0)

isoform_0excluded_combined <- mutate(isoform_0excluded_combined, geneid = gsub("\\..*", "", geneid))
isoform_0excluded_combined <- mutate(isoform_0excluded_combined, gene_biotype = genes_ens$gene_biotype[match(isoform_0excluded_combined$geneid, genes_ens$ensembl_gene_id)])
View(isoform_0excluded_combined)


isoform_0excluded_combined_avg <- isoform_0excluded_combined %>%
  group_by(transcript) %>%
  summarise(
    Average_acountsu = mean(acountsu),
    Average_acountss = mean(acountss),
    biotype = first(gene_biotype),
    gene = first(gene),
    fragmentation = (Average_acountss/Average_acountsu)
    )
View(isoform_0excluded_combined_avg)

quantile(isoform_0excluded_combined_avg$Average_acountsu, probs = 0.999)
isoform_0excluded_combined_avg <-  mutate(isoform_0excluded_combined_avg, fragmentation_status = ifelse(fragmentation > 0.8 & Average_acountsu > 4277.067  , "Highlighted", "Normal"))
View(isoform_0excluded_combined_avg)

ggplot(isoform_0excluded_combined_avg, aes(x = Average_acountsu, y = fragmentation, color = fragmentation_status)) +
  geom_point() +
  scale_x_log10() +
  scale_color_manual(values = c("Normal" = grey(0.5), "Highlighted" = "#00C1A7")) +
  labs (x = 'Counts_a_unique (log)', y = 'Fragmentation', color = '') +
  theme_minimal()

table(isoform_0excluded_combined_avg$fragmentation_status)

ggplot(isoform_0excluded_combined_avg, aes(x = Average_acountsu, y = fragmentation, color = fragmentation_status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_x_log10() +
  scale_color_manual(values = c("Normal" = grey(0.5), "Highlighted" = "#00C1A7")) +
  labs(x = 'Unique poly(A) counts', y = 'Fragmentation', color = '') +
  theme_minimal() +
  geom_label_repel(
    data = isoform_0excluded_combined_avg %>% filter(fragmentation_status == "Highlighted"), 
    aes(label = gene),
    size = 5,                  
    box.padding = 0.8,         
    point.padding = 0.7,       
    force = 5,                 
    max.overlaps = Inf,        
    nudge_x = 0.3,             
    nudge_y = 0.3,
    fill = "white",            
    color = "black",          
    segment.size = 0.5,       
    segment.color = "black"   
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

isoform_0excluded_combined_avg <- summarise(group_by(isoform_0excluded_combined_avg, fragmentation_status, biotype), 'status' = n())
View(isoform_0excluded_combined_avg)


#data biotypes 
isoform_combined <- bind_rows(isoform_d11344_sorted, isoform_d11345_sorted, isoform_d11358_sorted)
View(isoform_combined)
isoform_combined <- mutate(isoform_combined, geneid = gsub("\\..*", "", geneid))
isoform_combined <- mutate(isoform_combined, gene_biotype = genes_ens$gene_biotype[match(isoform_combined$geneid, genes_ens$ensembl_gene_id)]) 
isoform_combined_acountss_excluded5 <- subset(isoform_combined, acountss > 5)

      #plot van biotype verdeling (kort)
isoform_combined_biotype <- isoform_combined_acountss_excluded5 %>%
  group_by(Sample, gene_biotype) %>%
  summarise(Biotype_count = sum(acountss), .groups = "drop") %>%  
  group_by(Sample) %>%
  mutate(Biotype_percent = Biotype_count / sum(Biotype_count) * 100) %>%
  mutate(gene_biotype = ifelse(Biotype_percent < 1, "less than 1%", gene_biotype))
View(isoform_combined_biotype) 

     
ggplot(isoform_combined_biotype, aes(x = Sample, y = Biotype_percent, fill = gene_biotype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  labs(x = NULL, y = "Fraction of reads", fill = '') +
  theme_minimal()+
  scale_fill_manual(values = c('darkgrey','#C39BD3', '#00C1A7', '#A3C1DA'))
                   

ggplot(isoform_combined_biotype, aes(x = gene_biotype, y = Biotype_count, fill = gene_biotype)) +
  geom_bar(stat = "identity", position = "dodge") +  
  scale_y_log10() +
  facet_wrap(~Sample) +
  labs(x = 'Gene biotype', y = "Counts per Biotype (log)", fill = NULL) +  
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    panel.spacing = unit(1.5, "lines"),
    legend.text = element_text(size = 10)
  )

       #plot van biotype verdeling (lang, gemiddelde)
isoform_combined_biotype_all <- isoform_combined_acountss_excluded5 %>%
  group_by(Sample, gene_biotype) %>%
  summarise(Biotype_count = sum(acountss), .groups = "drop") %>%  
  group_by(Sample) %>%
  mutate(Biotype_percent = Biotype_count / sum(Biotype_count) * 100) 
View(isoform_combined_biotype_all) 


isoform_combined_biotype_avg <- isoform_combined_biotype_all %>%
  group_by(gene_biotype) %>%
  summarise(Average_Biotype_count = mean(Biotype_count))
View(isoform_combined_biotype_avg)

ggplot(isoform_combined_biotype_avg, aes(x = gene_biotype, y = Average_Biotype_count)) +
  geom_bar(stat = "identity", position = "dodge", fill = 'grey', color = 'black') +  
  scale_y_log10() +
  labs(x = 'Biotype', y = "Average number of intact poly(A) counts ", fill = NULL) +  
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    panel.spacing = unit(1.5, "lines"),
    legend.text = element_text(size = 10)
  )
