library(readr)
isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s <- read_delim("C:/Users/bodec/Downloads/plasma_LR/plasma_LR/isoform_counts-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s.tsv", 
                                                                                        delim = "\t", escape_double = FALSE, 
                                                                                        col_types = cols(begin = col_number(), 
                                                                                                         end = col_number(), exonStarts = col_character(), 
                                                                                                         exonEnds = col_character(), cdsStart = col_number(), 
                                                                                                         cdsEnd = col_number(), size = col_number(), 
                                                                                                         `counts_iqall-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number(), 
                                                                                                         `counts_weighed-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number(), 
                                                                                                         `counts_unique-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number(), 
                                                                                                         `counts_strict-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number(), 
                                                                                                         `counts_aweighed-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number(), 
                                                                                                         `counts_aunique-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number(), 
                                                                                                         `counts_astrict-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s` = col_number()), 
                                                                                        trim_ws = TRUE, skip = 27)



#data ordenen 
colnames(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s)
names(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s) == 'counts_aweighed-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s'] <- 'acountsw' 
names(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s) == 'counts_aunique-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s'] <- 'acountsu'
names(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s)[names(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s) == 'counts_astrict-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s'] <- 'acountss'
colnames(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s)

isoform_d11345_sorted <- isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s[order(isoform_counts_isoquant_joint_sminimap2_splice_d11345_v7_6_8_210333_hg38s$acountsw,decreasing = TRUE),]
isoform_d11345_sorted <- mutate(isoform_d11345_sorted,Sample = 'd11345') 
View(isoform_d11345_sorted)


#poly(A) data 
d11345_som_weighted <- sum(isoform_d11345_sorted$`counts_weighed-isoquant_joint-sminimap2_splice-d11345_v7.6.8_210333_hg38s`, na.rm = TRUE)
head(d11345_som_weighted)
d11345_som_aweigthed <- sum(isoform_d11345_sorted$acountsw, na.rm = TRUE)
head(d11345_som_aweigthed)
d11345_fraction <- d11345_som_aweigthed/d11345_som_weighted
head(d11345_fraction)

#data fractie van verschillende counts
d11345_som_aweigthed <- sum(isoform_d11345_sorted$acountsw, na.rm = TRUE)
head(d11345_som_aweigthed)
d11345_som_aunique <- sum(isoform_d11345_sorted$acountsu, na.rm = TRUE)
head(d11345_som_aunique)
d11345_som_astrict <- sum(isoform_d11345_sorted$acountss, na.rm = TRUE)
head(d11345_som_astrict)

d11345_counts <- data.frame(
  category = c("Unique", "Intact", "Weighted"),
  value = c(d11345_som_aunique, d11345_som_astrict, d11345_som_aweigthed),
  ID = c('d11345')
)


#transcript data 
isoform_d11345_transcripts <- summarise(group_by(isoform_d11345_sorted, gene), 'aantal transcripten' = n())
isoform_d11345_transcripts <- mutate(isoform_d11345_transcripts, Sample = 'd111345')
isoform_d11345_transcripts <- isoform_d11345_transcripts[order(isoform_d11345_transcripts$'aantal transcripten',decreasing = TRUE),]
View(isoform_d11345_transcripts)




#data verwerking voor cumulatieve plot van counts gescheiden per categrie 
isoform_d11345_acountss_excluded0 <- subset(isoform_d11345_sorted, acountss !=0)
View(isoform_d11345_acountss_excluded0)


#novel data
isoform_d11345_filtered_known <- subset(isoform_d11345_sorted, subset = category == 'known')
isoform_d11345_filtered_known_excluded0 <- subset(isoform_d11345_filtered_known, acountss !=0)
isoform_d11345_filtered_novelincatalog <- subset(isoform_d11345_sorted, subset = category == 'novel_in_catalog')
isoform_d11345_filtered_novelnotincatalog <- subset(isoform_d11345_sorted, subset = category == 'novel_not_in_catalog')


#data fragmentatie
isoform_d11345_excluded0 <- subset(isoform_d11345_sorted, acountsu !=0)





#vanaf hier ook zeker database namen nakijken want kloppen niet helemaal meer

isoform_d11345_20 <- head(isoform_d11345_filtered, 20)
View(isoform_d11345_20)




isoform_d11345_transcripts_20 <- head(isoform_d11345_transcripts, 20)
View(isoform_d11345_transcripts_20)


isoform_d11345_transcripts_20 <- mutate(isoform_d11345_transcripts_20,Sample = 'd11345') 
View(isoform_d11345_transcripts_20 )


