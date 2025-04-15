library(readr)
isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s <- read_delim("C:/Users/bodec/Downloads/plasma_LR/plasma_LR/isoform_counts-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s.tsv", 
                                                                                        delim = "\t", escape_double = FALSE, 
                                                                                        col_types = cols(begin = col_number(), 
                                                                                                         end = col_number(), exonStarts = col_character(), 
                                                                                                         exonEnds = col_character(), cdsStart = col_number(), 
                                                                                                         cdsEnd = col_number(), size = col_number(), 
                                                                                                         `counts_iqall-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number(), 
                                                                                                         `counts_weighed-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number(), 
                                                                                                         `counts_unique-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number(), 
                                                                                                         `counts_strict-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number(), 
                                                                                                         `counts_aweighed-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number(), 
                                                                                                         `counts_aunique-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number(), 
                                                                                                         `counts_astrict-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s` = col_number()), 
                                                                                        trim_ws = TRUE, skip = 27)





#data ordenen 
colnames(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s)
names(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s) == 'counts_aweighed-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s'] <- 'acountsw' 
names(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s) == 'counts_aunique-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s'] <- 'acountsu'
names(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s) == 'counts_astrict-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s'] <- 'acountss'
colnames(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s)

isoform_d11358_sorted <- isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s[order(isoform_counts_isoquant_joint_sminimap2_splice_d11358_v7_6_8_210334_hg38s$acountsw,decreasing = TRUE),]
isoform_d11358_sorted <- mutate(isoform_d11358_sorted,Sample = 'd11358') 
View(isoform_d11358_sorted)


#poly(A) data 
d11358_som_weighted <- sum(isoform_d11358_sorted$`counts_weighed-isoquant_joint-sminimap2_splice-d11358_v7.6.8_210334_hg38s`, na.rm = TRUE)
head(d11358_som_weighted)
d11358_som_aweigthed <- sum(isoform_d11358_sorted$acountsw, na.rm = TRUE)
head(d11358_som_aweigthed)
d11358_fraction <- d11358_som_aweigthed/d11358_som_weighted
head(d11358_fraction)


#data fractie van verschillende counts
d11358_som_aweigthed <- sum(isoform_d11358_sorted$acountsw, na.rm = TRUE)
head(d11358_som_aweigthed)
d11358_som_aunique <- sum(isoform_d11358_sorted$acountsu, na.rm = TRUE)
head(d11358_som_aunique)
d11358_som_astrict <- sum(isoform_d11358_sorted$acountss, na.rm = TRUE)
head(d11358_som_astrict)

d11358_counts <- data.frame(
  category = c("Unique", "Intact", "Weighted"),
  value = c(d11358_som_aunique, d11358_som_astrict, d11358_som_aweigthed),
  ID = c('d11358')
)


#transcript data
isoform_d11358_transcripts <- summarise(group_by(isoform_d11358_sorted, gene), 'aantal transcripten' = n())
isoform_d11358_transcripts <- mutate(isoform_d11358_transcripts, Sample = 'd111358')
isoform_d11358_transcripts <- isoform_d11358_transcripts[order(isoform_d11358_transcripts$'aantal transcripten',decreasing = TRUE),]
View(isoform_d11358_transcripts)



#data verwerking voor cumulatieve plot van counts gescheiden per categrie 
isoform_d11358_acountss_excluded0 <- subset(isoform_d11358_sorted, acountss !=0)
View(isoform_d11358_acountss_excluded0)



#novel data
isoform_d11358_filtered_known <- subset(isoform_d11358_sorted, subset = category == 'known')
isoform_d11358_filtered_known_excluded0 <- subset(isoform_d11358_filtered_known, acountss !=0)
isoform_d11358_filtered_novelincatalog <- subset(isoform_d11358_sorted, subset = category == 'novel_in_catalog')
isoform_d11358_filtered_novelnotincatalog <- subset(isoform_d11358_sorted, subset = category == 'novel_not_in_catalog')

#data fragmentatie
isoform_d11358_excluded0 <- subset(isoform_d11358_sorted, acountsu !=0)



#datasets kakijken!!

isoform_d11358_20 <- head(isoform_d11358_filtered, 20)
View(isoform_d11358_20)



isoform_d11358_transcripts_20 <- head(isoform_d11358_transcripts, 20)
View(isoform_d11358_transcripts_20)


isoform_d11358_transcripts_20 <- mutate(isoform_d11358_transcripts_20,Sample = 'd11358') 
View(isoform_d11358_transcripts_20 )



