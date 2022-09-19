library(tidyverse)
library(readr)

ref_val <- read.delim("2A/provalidated_2a.txt", header = F)
ref_val <- distinct (ref_val)

#write.table(both8Ruh, file = "8R_unique_both.txt", sep=" ", quote= F)


#Alignments analysis
met_BIC_a<- read_tsv("2A/score method analysis/alignments/PSSMsearch_2apro_9wR_BLASTIC_alignments.txt")
met_B_a <-read_tsv("2A/score method analysis/alignments/PSSMsearch_2apro_9wR_BLASTIC_alignments.txt")
met_M_a <- read_tsv("2A/score method analysis/alignments/PSSMSearch_2A_9R_MOTIPS_alignments.txt")
met_bin_a <- read_tsv("2A/score method analysis/alignments/PSSMSearch_2A_9R_binomial_alignments.txt")
met_ratio_a <- read_tsv("2A/score method analysis/alignments/PSSMSearch_2A_9R_ratio_alignments.txt")

nt_bic_a <- nrow(met_BIC_a)
nt_b_a <- nrow(met_B_a)
nt_m_a <- nrow(met_M_a)
nt_bin_a <- nrow(met_bin_a)
nt_ratio_a <- nrow(met_ratio_a)


npr_bic_a<- sum(met_BIC_a$GeneName %in% ref_val$V1)
npr_b_a<- sum(met_B_a$GeneName %in% ref_val$V1)
npr_m_a<- sum(met_M_a$GeneName %in% ref_val$V1)
npr_bin_a<- sum(met_bin_a$GeneName %in% ref_val$V1)
npr_ratio_a<- sum(met_ratio_a$GeneName %in% ref_val$V1)

#ratios
r_bic_a<- npr_bic_a/nt_bic_a
r_b_a<- npr_b_a/nt_b_a
r_m_a<- npr_m_a/nt_m_a
r_bin_a<- npr_bin_a/nt_bin_a
r_ratio_a<- npr_ratio_a/nt_ratio_a


#plotting

scoreplot<- data.frame(Sequences=rep(c("alignments"), each=5), meth=c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio'),
                     rat=c(r_bic_a, r_b_a, r_m_a, r_bin_a, r_ratio_a), cant=c(npr_bic_a, npr_b_a, npr_m_a, npr_bin_a, npr_ratio_a))
ggplot(scoreplot, aes(meth, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background with 8R. Score method comparison")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'BLAST-IC', 'MOTIPS', 'Binomial log10', 'Ratio')) +
  labs(y="Ratio", x= "Score Method") +
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()


#DMS analysis
met_BIC_d<- read_tsv("2A/score method analysis/DMS/PSSMsearch_2apro_9wR_BLASTIC_DMS.txt")
met_B_d <-read_tsv("2A/score method analysis/DMS/PSSMsearch_2apro_9wR_BLAST_DMS.txt")
met_M_d <- read_tsv("2A/score method analysis/DMS/PSSMsearch_2apro_9wR_MOTIPS_DMS.txt")
met_bin_d <- read_tsv("2A/score method analysis/DMS/PSSMsearch_2apro_9wR_binomialog10_DMS.txt")
met_ratio_d <- read_tsv("2A/score method analysis/DMS/PSSMsearch_2apro_9wR_ratio_DMS.txt")



nt_bic_d <- nrow(met_BIC_d)
nt_b_d <- nrow(met_B_d)
nt_m_d <- nrow(met_M_d)
nt_bin_d <- nrow(met_bin_d)
nt_ratio_d <- nrow(met_ratio_d)


npr_bic_d<- sum(met_BIC_d$GeneName %in% ref_val$V1)
npr_b_d<- sum(met_B_d$GeneName %in% ref_val$V1)
npr_m_d<- sum(met_M_d$GeneName %in% ref_val$V1)
npr_bin_d<- sum(met_bin_d$GeneName %in% ref_val$V1)
npr_ratio_d<- sum(met_ratio_d$GeneName %in% ref_val$V1)

#ratios
r_bic_d<- npr_bic_d/nt_bic_d
r_b_d<- npr_b_d/nt_b_d
r_m_d<- npr_m_d/nt_m_d
r_bin_d<- npr_bin_d/nt_bin_d
r_ratio_d<- npr_ratio_d/nt_ratio_d


#plotting

scoreplot<- data.frame(Sequences=rep(c("alignments", "DMS"), each=5), meth=c('PSI-BLAST-IC', 'PSI-BLAST', 'MOTIPS', 'Binomial log10', 'Ratio'),
                       rat=c(r_bic_a, r_b_a, r_m_a, r_bin_a, r_ratio_a, r_bic_d, r_b_d, r_m_d, r_bin_d, r_ratio_d), cant=c(npr_bic_a, npr_b_a, npr_m_a, npr_bin_a, npr_ratio_a,  npr_bic_d, npr_b_d, npr_m_d, npr_bin_d, npr_ratio_d))
ggplot(scoreplot, aes(meth, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background with 9R. Score method comparison (Raw data)")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'PSI-BLAST', 'MOTIPS', 'Binomial log10', 'Ratio')) +
  labs(y="Ratio", x= "Score Method") +
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .5,vjust=4, size = 4)+ 
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


scoreplot_g<- data.frame(Sequences=rep(c("alignments", "DMS"), each=5), meth=c('PSI-BLAST-IC', 'PSI-BLAST', 'MOTIPS', 'Binomial log10', 'Ratio'),
                       rat=c(r_bic_a_g, r_b_a_g, r_m_a_g, r_bin_a_g, r_ratio_a_g, r_bic_d_g, r_b_d_g, r_m_d_g, r_bin_d_g, r_ratio_d_g), cant=c(g_bic_a, g_b_a, g_m_a, g_bin_a, g_ratio_a,  g_bic_d, g_b_d, g_m_d, g_bin_d, g_ratio_d))
ggplot(scoreplot_g, aes(meth, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background with 9R. Score method comparison (Unique Gene Names)")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'PSI-BLAST', 'MOTIPS', 'Binomial log10', 'Ratio')) +
  labs(y="Ratio", x= "Score Method") +
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()

scoreplot_h<- data.frame(Sequences=rep(c("alignments", "DMS"), each=5), meth=c('PSI-BLAST-IC', 'PSI-BLAST', 'MOTIPS', 'Binomial log10', 'Ratio'),
                       rat=c(r_bic_a_h, r_b_a_h, r_m_a_h, r_bin_a_h, r_ratio_a_h, r_bic_d_h, r_b_d_h, r_m_d_h, r_bin_d_h, r_ratio_d_h), cant=c(h_bic_a, h_b_a, h_m_a, h_bin_a, h_ratio_a,  h_bic_d, h_b_d, h_m_d, h_bin_d, h_ratio_d))
ggplot(scoreplot_h, aes(meth, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background with 9R. Score method comparison (Unique Cleavage Sites)")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'PSI-BLAST', 'MOTIPS', 'Binomial log10', 'Ratio')) +
  labs(y="Ratio", x= "Score Method") +
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()

sum((uh_ratio_d$GeneName %in% ref_val$V1) %in% (uh_bic_a$GeneName %in% ref_val$V1))

test_ratio_d <- uh_ratio_d[uh_ratio_d$GeneName %in% ref_val$V1,]
test_bic_a <- uh_bic_a[uh_bic_a$GeneName %in% ref_val$V1,]
test_bic_d<- uh_bic_d[uh_bic_d$GeneName %in% ref_val$V1,]

sum(test_ratio_d$GeneName %in% test_bic_a$GeneName)


#combining A and D

both_bic<-rbind(met_BIC_a,met_BIC_d,by="ProteinAcc")
both_b<-rbind(met_B_a,met_B_d,by="ProteinAcc")
both_m<-rbind(met_M_a,met_M_d,by="ProteinAcc")
both_bin<-rbind(met_bin_a,met_bin_d,by="ProteinAcc")
both_ratio<-rbind(met_ratio_a,met_ratio_d,by="ProteinAcc")

both_bic_g<-distinct(both_bic,GeneName,.keep_all = T)
both_b_g<-distinct(both_b,GeneName,.keep_all = T)
both_m_g<-distinct(both_m,GeneName,.keep_all = T)
both_bin_g<-distinct(both_bin,GeneName,.keep_all = T)
both_ratio_g<-distinct(both_ratio,GeneName,.keep_all = T)

both_bic<- both_bic[grepl("TG|YG|FG|LG|MG", both_bic$Hit),]
both_b<-both_b[grepl("TG|YG|FG|LG|MG", both_b$Hit),]
both_m<-both_m[grepl("TG|YG|FG|LG|MG", both_m$Hit),]
both_bin<-both_bin[grepl("TG|YG|FG|LG|MG", both_bin$Hit),]
both_ratio<-both_ratio[grepl("TG|YG|FG|LG|MG", both_ratio$Hit),]

both_bic_h<-distinct(both_bic,GeneName, Hit,.keep_all = T)
both_b_h<-distinct(both_b,GeneName,Hit,.keep_all = T)
both_m_h<-distinct(both_m,GeneName,Hit,.keep_all = T)
both_bin_h<-distinct(both_bin,GeneName,Hit,.keep_all = T)
both_ratio_h<-distinct(both_ratio,GeneName,Hit,.keep_all = T)

total_bbic<-nrow(both_bic)
total_bb<-nrow(both_b)
total_bm<-nrow(both_m)
total_bbin<-nrow(both_bin)
total_bratio<-nrow(both_ratio)


nboth_bic<- sum(both_bic$GeneName %in% ref_val$V1)
nboth_b<- sum(both_b$GeneName %in% ref_val$V1)
nboth_m<- sum(both_m$GeneName %in% ref_val$V1)
nboth_bin<- sum(both_bin$GeneName %in% ref_val$V1)
nboth_ratio<- sum(both_ratio$GeneName %in% ref_val$V1)

rat_bth_ic<-nboth_bic/total_bbic
rat_bth_b<-nboth_b/total_bb
rat_bth_m<-nboth_m/total_bm
rat_bth_bin<-nboth_bin/total_bbin
rat_bth_ratio<-nboth_ratio/total_bratio


nboth_bic_g<- sum(both_bic_g$GeneName %in% ref_val$V1)
nboth_b_g<- sum(both_b_g$GeneName %in% ref_val$V1)
nboth_m_g<- sum(both_m_g$GeneName %in% ref_val$V1)
nboth_bin_g<- sum(both_bin_g$GeneName %in% ref_val$V1)
nboth_ratio_g<- sum(both_ratio_g$GeneName %in% ref_val$V1)

rat_bth_icg<-nboth_bic_g/total_bbic
rat_bth_bg<-nboth_b_g/total_bb
rat_bth_mg<-nboth_m_g/total_bm
rat_bth_bing<-nboth_bin_g/total_bbin
rat_bth_ratg<-nboth_ratio_g/total_bratio


nboth_bic_h<- sum(both_bic_h$GeneName %in% ref_val$V1)
nboth_b_h<- sum(both_b_h$GeneName %in% ref_val$V1)
nboth_m_h<- sum(both_m_h$GeneName %in% ref_val$V1)
nboth_bin_h<- sum(both_bin_h$GeneName %in% ref_val$V1)
nboth_ratio_h<- sum(both_ratio_h$GeneName %in% ref_val$V1)

rat_bth_ich<-nboth_bic_h/total_bbic
rat_bth_bh<-nboth_b_h/total_bb
rat_bth_mh<-nboth_m_h/total_bm
rat_bth_binh<-nboth_bin_h/total_bbin
rat_bth_rath<-nboth_ratio_h/total_bratio


ratplotb2<- data.frame(Source=rep(c("raw", "unique gene names", "unique cleavage sites"), each=5),meth=c('PSI-BLAST-IC', 'PSI-BLAST', 'MOTIPS', 'Binomial log10', 'Ratio'),
                      rat=c(rat_bth_ic,rat_bth_b, rat_bth_m, rat_bth_bin,rat_bth_ratio,rat_bth_icg, rat_bth_bg, rat_bth_mg, rat_bth_bing, rat_bth_ratg, rat_bth_ich, rat_bth_bh, rat_bth_mh, rat_bth_binh, rat_bth_rath), 
                      cant=c(nboth_bic, nboth_b, nboth_m, nboth_bin, nboth_ratio,nboth_bic_g, nboth_b_g, nboth_m_g, nboth_bin_g, nboth_ratio_g,nboth_bic_h, nboth_b_h, nboth_m_h, nboth_bin_h, nboth_ratio_h))

ggplot(ratplotb2, aes(meth, rat, cant, label=Source)) + ggtitle("Ratio hits-background (Alignments-DMS)")+
  geom_bar(aes(fill = Source), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('PSI-BLAST-IC', 'PSI-BLAST', 'MOTIPS', 'Binomial log10', 'Ratio')) +
  labs(y="Ratio", x= "Score Method")+
  geom_text(aes(label= cant, group = Source),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal(base_size = 19)


