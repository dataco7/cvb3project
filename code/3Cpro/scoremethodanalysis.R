library(tidyverse)
library(readr)

ref_val <- read.delim("3C/provalidated.txt", header = F)
ref_val <- distinct (ref_val)

#write.table(both8Ruh, file = "8R_unique_both.txt", sep=" ", quote= F)


#Alignments analysis
met_BIC_a<- read_tsv("3C/pssmsearch method analysis/alignments/PSSMSearch_3C_results_BLASTIC_8R.txt")
met_B_a <-read_tsv("3C/pssmsearch method analysis/alignments/PSSMSearch_3C_results_BLAST_8R.txt")
met_M_a <- read_tsv("3C/pssmsearch method analysis/alignments/PSSMSearch_3C_results_MOTIPS_8R.txt")
met_bin_a <- read_tsv("3C/pssmsearch method analysis/alignments/PSSMSearch_3C_results_binomial_8R.txt")
met_ratio_a <- read_tsv("3C/pssmsearch method analysis/alignments/PSSMSearch_3C_results_ratio_8R.txt")
met_odds_a <- read_tsv("3C/pssmsearch method analysis/alignments/PSSMSearch_3C_results_odds_8R.txt")


nt_bic_a <- nrow(met_BIC_a)
nt_b_a <- nrow(met_B_a)
nt_m_a <- nrow(met_M_a)
nt_bin_a <- nrow(met_bin_a)
nt_ratio_a <- nrow(met_ratio_a)
nt_odds_a <- nrow(met_odds_a)


npr_bic_a<- sum(met_BIC_a$GeneName %in% ref_val$V1)
npr_b_a<- sum(met_B_a$GeneName %in% ref_val$V1)
npr_m_a<- sum(met_M_a$GeneName %in% ref_val$V1)
npr_bin_a<- sum(met_bin_a$GeneName %in% ref_val$V1)
npr_ratio_a<- sum(met_ratio_a$GeneName %in% ref_val$V1)
npr_odds_a<- sum(met_odds_a$GeneName %in% ref_val$V1)

#ratios
r_bic_a<- npr_bic_a/nt_bic_a
r_b_a<- npr_b_a/nt_b_a
r_m_a<- npr_m_a/nt_m_a
r_bin_a<- npr_bin_a/nt_bin_a
r_ratio_a<- npr_ratio_a/nt_ratio_a
r_odds_a<- npr_odds_a/nt_odds_a


#plotting

scoreplot<- data.frame(Sequences=rep(c("alignments"), each=6), meth=c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio', 'Log Odds'),
                     rat=c(r_bic_a, r_b_a, r_m_a, r_bin_a, r_ratio_a, r_odds_a), cant=c(npr_bic_a, npr_b_a, npr_m_a, npr_bin_a, npr_ratio_a, npr_odds_a))
ggplot(scoreplot, aes(meth, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background with 8R. Score method comparison")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio', 'Log Odds')) +
  labs(y="Ratio", x= "Score Method") +
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()


#DMS analysis
met_BIC_d<- read_tsv("3C/pssmsearch method analysis/dms/PSSMSearch_results_3c_BLASTIC_DMS.txt")
met_B_d <-read_tsv("3C/pssmsearch method analysis/dms/PSSMSearch_results_3c_BLAST_DMS.txt")
met_M_d <- read_tsv("3C/pssmsearch method analysis/dms/PSSMSearch_results_3c_MOTIPS_DMS.txt")
met_bin_d <- read_tsv("3C/pssmsearch method analysis/dms/PSSMSearch_results_3c_Binomial log10_DMS.txt")
met_ratio_d <- read_tsv("3C/pssmsearch method analysis/dms/PSSMSearch_results_3c_ratio_DMS.txt")
met_odds_d <- read_tsv("3C/pssmsearch method analysis/alignments/PSSMSearch_3C_results_odds_8R.txt")


nt_bic_d <- nrow(met_BIC_d)
nt_b_d <- nrow(met_B_d)
nt_m_d <- nrow(met_M_d)
nt_bin_d <- nrow(met_bin_d)
nt_ratio_d <- nrow(met_ratio_d)
nt_odds_d <- nrow(met_odds_d)


npr_bic_d<- sum(met_BIC_d$GeneName %in% ref_val$V1)
npr_b_d<- sum(met_B_d$GeneName %in% ref_val$V1)
npr_m_d<- sum(met_M_d$GeneName %in% ref_val$V1)
npr_bin_d<- sum(met_bin_d$GeneName %in% ref_val$V1)
npr_ratio_d<- sum(met_ratio_d$GeneName %in% ref_val$V1)
npr_odds_d<- sum(met_odds_d$GeneName %in% ref_val$V1)

#ratios
r_bic_d<- npr_bic_d/nt_bic_d
r_b_d<- npr_b_d/nt_b_d
r_m_d<- npr_m_d/nt_m_d
r_bin_d<- npr_bin_d/nt_bin_d
r_ratio_d<- npr_ratio_d/nt_ratio_d
r_odds_d<- npr_odds_d/nt_odds_d


#plotting

scoreplot<- data.frame(Sequences=rep(c("alignments", "DMS"), each=5), meth=c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio'),
                       rat=c(r_bic_a, r_b_a, r_m_a, r_bin_a, r_ratio_a, r_bic_d, r_b_d, r_m_d, r_bin_d, r_ratio_d), cant=c(npr_bic_a, npr_b_a, npr_m_a, npr_bin_a, npr_ratio_a,  npr_bic_d, npr_b_d, npr_m_d, npr_bin_d, npr_ratio_d))
ggplot(scoreplot, aes(meth, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background with 8R. Score method comparison")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio')) +
  labs(y="Ratio", x= "Score Method") +
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()

ug_bic_a <-distinct(met_BIC_a, GeneName, .keep_all = T)
ug_b_a<-distinct(met_B_a, GeneName, .keep_all = T)
ug_m_a<-distinct(met_M_a, GeneName, .keep_all = T)
ug_bin_a<-distinct(met_bin_a, GeneName, .keep_all = T)
ug_ratio_a<-distinct(met_ratio_a, GeneName, .keep_all = T)

ug_bic_d<-distinct(met_BIC_d, GeneName, .keep_all = T)
ug_b_d<-distinct(met_B_d, GeneName, .keep_all = T)
ug_m_d<-distinct(met_M_d, GeneName, .keep_all = T)
ug_bin_d<-distinct(met_bin_d, GeneName, .keep_all = T)
ug_ratio_d<-distinct(met_ratio_d, GeneName, .keep_all = T)


uh_bic_a<-distinct(met_BIC_a, GeneName, Hit, .keep_all = T)
uh_b_a<-distinct(met_B_a, GeneName,Hit, .keep_all = T)
uh_m_a<-distinct(met_M_a, GeneName,Hit, .keep_all = T)
uh_bin_a<-distinct(met_bin_a, GeneName, Hit, .keep_all = T)
uh_ratio_a<-distinct(met_ratio_a, GeneName, .keep_all = T)

uh_bic_d<-distinct(met_BIC_d, GeneName, Hit, .keep_all = T)
uh_b_d<-distinct(met_B_d, GeneName, Hit, .keep_all = T)
uh_m_d<-distinct(met_M_d, GeneName, Hit, .keep_all = T)
uh_bin_d<-distinct(met_bin_d, GeneName,Hit, .keep_all = T)
uh_ratio_d<-distinct(met_ratio_d, GeneName, Hit, .keep_all = T)

#Unique Genes

g_bic_a<- sum(ug_bic_a$GeneName %in% ref_val$V1)
g_b_a<- sum(ug_b_a$GeneName %in% ref_val$V1)
g_m_a<- sum(ug_m_a$GeneName %in% ref_val$V1)
g_bin_a<- sum(ug_bin_a$GeneName %in% ref_val$V1)
g_ratio_a<- sum(ug_ratio_a$GeneName %in% ref_val$V1)


g_bic_d<- sum(ug_bic_d$GeneName %in% ref_val$V1)
g_b_d<- sum(ug_b_d$GeneName %in% ref_val$V1)
g_m_d<- sum(ug_m_d$GeneName %in% ref_val$V1)
g_bin_d<- sum(ug_bin_d$GeneName %in% ref_val$V1)
g_ratio_d<- sum(met_ratio_d$GeneName %in% ref_val$V1)

#Unique Hits

h_bic_a<- sum(uh_bic_a$GeneName %in% ref_val$V1)
h_b_a<- sum(uh_b_a$GeneName %in% ref_val$V1)
h_m_a<- sum(uh_m_a$GeneName %in% ref_val$V1)
h_bin_a<- sum(uh_bin_a$GeneName %in% ref_val$V1)
h_ratio_a<- sum(uh_ratio_a$GeneName %in% ref_val$V1)


h_bic_d<- sum(uh_bic_d$GeneName %in% ref_val$V1)
h_b_d<- sum(uh_b_d$GeneName %in% ref_val$V1)
h_m_d<- sum(uh_m_d$GeneName %in% ref_val$V1)
h_bin_d<- sum(uh_bin_d$GeneName %in% ref_val$V1)
h_ratio_d<- sum(uh_ratio_d$GeneName %in% ref_val$V1)

#ratios
r_bic_a_g<-g_bic_a/nt_bic_a
r_b_a_g<-g_b_a/nt_b_a
r_m_a_g<-g_m_a/nt_m_a
r_bin_a_g<-g_bin_a/nt_bin_a
r_ratio_a_g<-g_ratio_a/nt_ratio_a
r_bic_d_g<-g_bic_d/nt_bic_d
r_b_d_g<-g_b_d/nt_b_d
r_m_d_g<-g_m_d/nt_m_d
r_bin_d_g<-g_bin_d/nt_bin_d
r_ratio_d_g<-g_ratio_d/nt_ratio_d

r_bic_a_h<-h_bic_a/nt_bic_a
r_b_a_h<-h_b_a/nt_b_a
r_m_a_h<-h_m_a/nt_m_a
r_bin_a_h<-h_bin_a/nt_bin_a
r_ratio_a_h<-h_ratio_a/nt_ratio_a
r_bic_d_h<-h_bic_d/nt_bic_d
r_b_d_h<-h_b_d/nt_b_d
r_m_d_h<-h_m_d/nt_m_d
r_bin_d_h<-h_bin_d/nt_bin_d
r_ratio_d_h<-h_ratio_d/nt_ratio_d


scoreplot_g<- data.frame(Sequences=rep(c("alignments", "DMS"), each=5), meth=c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio'),
                       rat=c(r_bic_a_g, r_b_a_g, r_m_a_g, r_bin_a_g, r_ratio_a_g, r_bic_d_g, r_b_d_g, r_m_d_g, r_bin_d_g, r_ratio_d_g), cant=c(g_bic_a, g_b_a, g_m_a, g_bin_a, g_ratio_a,  g_bic_d, g_b_d, g_m_d, g_bin_d, g_ratio_d))
ggplot(scoreplot_g, aes(meth, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background with 8R. Score method comparison (Unique Gene Names)")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio')) +
  labs(y="Ratio", x= "Score Method") +
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal(base_size = 18)

scoreplot_h<- data.frame(Sequences=rep(c("alignments", "DMS"), each=5), meth=c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio'),
                       rat=c(r_bic_a_h, r_b_a_h, r_m_a_h, r_bin_a_h, r_ratio_a_h, r_bic_d_h, r_b_d_h, r_m_d_h, r_bin_d_h, r_ratio_d_h), cant=c(h_bic_a, h_b_a, h_m_a, h_bin_a, h_ratio_a,  h_bic_d, h_b_d, h_m_d, h_bin_d, h_ratio_d))
ggplot(scoreplot_h, aes(meth, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background with 8R. Score method comparison (Unique Hits)")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio')) +
  labs(y="Ratio", x= "Score Method") +
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()

sum((uh_ratio_d$GeneName %in% ref_val$V1) %in% (uh_bic_a$GeneName %in% ref_val$V1))

test_ratio_d <- uh_ratio_d[uh_ratio_d$GeneName %in% ref_val$V1,]
test_bic_a <- uh_bic_a[uh_bic_a$GeneName %in% ref_val$V1,]
test_bic_d<- uh_bic_d[uh_bic_d$GeneName %in% ref_val$V1,]

sum(test_ratio_d$GeneName %in% test_bic_a$GeneName)
