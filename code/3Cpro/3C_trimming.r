#3CPRO

library(tidyverse)
library(readr)

ref_val <- read.delim("3C/provalidated.txt", header = F)
ref_val <- distinct (ref_val)

pep3align <- read_tsv("3C/EntB_3Cpro/3Cpro_EntB.txt", col_names = F)
pep3dms <- read_tsv("3C/DMS_3Cpro/3Cpro.txt", col_names = F)

### Assimetrical trimming 

#Alignments
#Right trimming
ali9R<-as.data.frame(gsub(".$","",pep3align$X1))
ali8R<-as.data.frame(gsub("..$","",pep3align$X1))
ali7R<-as.data.frame(gsub("...$","",pep3align$X1))
ali6R<-as.data.frame(gsub("....$","",pep3align$X1))
#Left trimming
ali9L<-as.data.frame(gsub("^.","",pep3align$X1))
ali8L<-as.data.frame(gsub("^..","",pep3align$X1))
ali7L<-as.data.frame(gsub("^...","",pep3align$X1))
ali6L<-as.data.frame(gsub("^....","",pep3align$X1))

write.table(ali9R, "3C/EntB_3Cpro/assimetrical trimming/ali3cpro_9R.txt", quote= F, col.names = F, row.names=F)
write.table(ali8R, "3C/EntB_3Cpro/assimetrical trimming/ali3cpro_8R.txt", quote= F, col.names = F, row.names=F)
write.table(ali7R, "3C/EntB_3Cpro/assimetrical trimming/ali3cpro_7R.txt", quote= F, col.names = F, row.names=F)
write.table(ali6R, "3C/EntB_3Cpro/assimetrical trimming/ali3cpro_6R.txt", quote= F, col.names = F, row.names=F)
write.table(ali9L, "3C/EntB_3Cpro/assimetrical trimming/ali3cpro_9L.txt", quote= F, col.names = F, row.names=F)
write.table(ali8L, "3C/EntB_3Cpro/assimetrical trimming/ali3cpro_8L.txt", quote= F, col.names = F, row.names=F)
write.table(ali7L, "3C/EntB_3Cpro/assimetrical trimming/ali3cpro_7L.txt", quote= F, col.names = F, row.names=F)
write.table(ali6L, "3C/EntB_3Cpro/assimetrical trimming/ali3cpro_6L.txt", quote= F, col.names = F, row.names=F)


#DMS
#Right trimming
dms9R<-as.data.frame(gsub(".$","",pep3dms$X1))
dms8R<-as.data.frame(gsub("..$","",pep3dms$X1))
dms7R<-as.data.frame(gsub("...$","",pep3dms$X1))
dms6R<-as.data.frame(gsub("....$","",pep3dms$X1))
#Left trimming
dms9L<-as.data.frame(gsub("^.","",pep3dms$X1))
dms8L<-as.data.frame(gsub("^..","",pep3dms$X1))
dms7L<-as.data.frame(gsub("^...","",pep3dms$X1))
dms6L<-as.data.frame(gsub("^....","",pep3dms$X1))

write.table(dms9R, "3C/DMS_3Cpro/assimetrical trimming/dms3cpro_9R.txt", quote= F, col.names = F, row.names=F)
write.table(dms8R, "3C/DMS_3Cpro/assimetrical trimming/dms3cpro_8R.txt", quote= F, col.names = F, row.names=F)
write.table(dms7R, "3C/DMS_3Cpro/assimetrical trimming/dms3cpro_7R.txt", quote= F, col.names = F, row.names=F)
write.table(dms6R, "3C/DMS_3Cpro/assimetrical trimming/dms3cpro_6R.txt", quote= F, col.names = F, row.names=F)
write.table(dms9L, "3C/DMS_3Cpro/assimetrical trimming/dms3cpro_9L.txt", quote= F, col.names = F, row.names=F)
write.table(dms8L, "3C/DMS_3Cpro/assimetrical trimming/dms3cpro_8L.txt", quote= F, col.names = F, row.names=F)
write.table(dms7L, "3C/DMS_3Cpro/assimetrical trimming/dms3cpro_7L.txt", quote= F, col.names = F, row.names=F)
write.table(dms6L, "3C/DMS_3Cpro/assimetrical trimming/dms3cpro_6L.txt", quote= F, col.names = F, row.names=F)





#COUNTING
##PSSMSearch outputs from trimmed cleavage sites sequences variants

entbic10 <- read_tsv("3C/EntB_3Cpro/3C_PSSM_results_EntB.txt")
entbic9R<- read_tsv("3C/EntB_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_EntB_w9R.txt")
entbic9L<-read_tsv("3C/EntB_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_EntB_w9L.txt")
entbic8R<- read_tsv("3C/EntB_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_EntB_w8R.txt")
entbic8L<- read_tsv("3C/EntB_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_EntB_w8L.txt")
entbic7R<- read_tsv("3C/EntB_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_EntB_w7R.txt")
entbic7L<- read_tsv("3C/EntB_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_EntB_w7L.txt")
entbic6R<- read_tsv("3C/EntB_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_EntB_w6R.txt")
entbic6L<- read_tsv("3C/EntB_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_EntB_w6L.txt")


dmsic10 <- read_tsv("3C/DMS_3Cpro/3Cpro_DMS_PSSMsearch_results.txt")
dmsic9R<- read_tsv("3C/DMS_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_w9R.txt")
dmsic9L<-read_tsv("3C/DMS_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_w9L.txt")
dmsic8R<- read_tsv("3C/DMS_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_w8R.txt")
dmsic8L<- read_tsv("3C/DMS_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_w8L.txt")
dmsic7R<- read_tsv("3C/DMS_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_w7R.txt")
dmsic7L<- read_tsv("3C/DMS_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_w7L.txt")
dmsic6R<- read_tsv("3C/DMS_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_w6R.txt")
dmsic6L<- read_tsv("3C/DMS_3Cpro/assimetrical trimming/3Cpro_PSSMsearch_results_w6L.txt")


#raw hit number from PSSMSearch output

aw10<-sum(entbic10$GeneName %in% ref_val$V1)
aw9R<-sum(entbic9R$GeneName %in% ref_val$V1)
aw9L<-sum(entbic9L$GeneName %in% ref_val$V1)
aw8R<-sum(entbic8R$GeneName %in% ref_val$V1)
aw8L<-sum(entbic8L$GeneName %in% ref_val$V1)
aw7R<-sum(entbic7R$GeneName %in% ref_val$V1)
aw7L<-sum(entbic7L$GeneName %in% ref_val$V1)
aw6R<-sum(entbic6R$GeneName %in% ref_val$V1)
aw6L<-sum(entbic6L$GeneName %in% ref_val$V1)

dw10 <-sum(dmsic10$GeneName %in% ref_val$V1)
dw9R<-sum(dmsic9R$GeneName %in% ref_val$V1)
dw9L<-sum(dmsic9L$GeneName %in% ref_val$V1)
dw8R<-sum(dmsic8R$GeneName %in% ref_val$V1)
dw8L<-sum(dmsic8L$GeneName %in% ref_val$V1)
dw7R<-sum(dmsic7R$GeneName %in% ref_val$V1)
dw7L<-sum(dmsic7L$GeneName %in% ref_val$V1)
dw6R<-sum(dmsic6R$GeneName %in% ref_val$V1)
dw6L<-sum(dmsic6L$GeneName %in% ref_val$V1)


#totals from PSSMSearch output
total10_d<- nrow(dmsic10)
total9R_d<- nrow(dmsic9R)
total9L_d<- nrow(dmsic9L)
total8R_d<- nrow(dmsic8R)
total8L_d<- nrow(dmsic8L)
total7R_d<- nrow(dmsic7R)
total7L_d<- nrow(dmsic7L)
total6R_d<- nrow(dmsic6R)
total6L_d<- nrow(dmsic6L)


total10_a<- nrow(entbic10)
total9R_a<- nrow(entbic9R)
total9L_a<- nrow(entbic9L)
total8R_a<- nrow(entbic8R)
total8L_a<- nrow(entbic8L)
total7R_a<- nrow(entbic7R)
total7L_a<- nrow(entbic7L)
total6R_a<- nrow(entbic6R)
total6L_a<- nrow(entbic6L)



#unique hits by name
entbic10u<-distinct(entbic10,GeneName, .keep_all = T)
entbic9Ru<-distinct(entbic9R,GeneName, .keep_all = T)
entbic9Lu<-distinct(entbic9L,GeneName, .keep_all = T)
entbic8Ru<-distinct(entbic8R,GeneName, .keep_all = T)
entbic8Lu<-distinct(entbic8L,GeneName, .keep_all = T)
entbic7Ru<-distinct(entbic7R,GeneName, .keep_all = T)
entbic7Lu<-distinct(entbic7L,GeneName, .keep_all = T)
entbic6Ru<-distinct(entbic6R,GeneName, .keep_all = T)
entbic6Lu<-distinct(entbic6L,GeneName, .keep_all = T)

dmsic10u<- distinct(dmsic10,GeneName,.keep_all = T)
dmsic9Ru<- distinct(dmsic9R,GeneName,.keep_all = T)
dmsic9Lu<- distinct(dmsic9L,GeneName,.keep_all = T)
dmsic8Ru<- distinct(dmsic8R,GeneName,.keep_all = T)
dmsic8Lu<- distinct(dmsic8L,GeneName,.keep_all = T)
dmsic7Ru<- distinct(dmsic7R,GeneName,.keep_all = T)
dmsic7Lu<- distinct(dmsic7L,GeneName,.keep_all = T)
dmsic6Ru<- distinct(dmsic6R,GeneName,.keep_all = T)
dmsic6Lu<- distinct(dmsic6L,GeneName,.keep_all = T)

  
dw10u <-sum(dmsic10u$GeneName %in% ref_val$V1)
dw9Ru<- sum(dmsic9Ru$GeneName %in% ref_val$V1)
dw9Lu<- sum(dmsic9Lu$GeneName %in% ref_val$V1)
dw8Ru<- sum(dmsic8Ru$GeneName %in% ref_val$V1)
dw8Lu<- sum(dmsic8Lu$GeneName %in% ref_val$V1)
dw7Ru<- sum(dmsic7Ru$GeneName %in% ref_val$V1)
dw7Lu<- sum(dmsic7Lu$GeneName %in% ref_val$V1)
dw6Ru<- sum(dmsic6Ru$GeneName %in% ref_val$V1)
dw6Lu<- sum(dmsic6Lu$GeneName %in% ref_val$V1)

aw10u<-sum(entbic10u$GeneName %in% ref_val$V1)
aw9Ru<-sum(entbic9Ru$GeneName %in% ref_val$V1)
aw9Lu<-sum(entbic9Lu$GeneName %in% ref_val$V1)
aw8Ru<-sum(entbic8Ru$GeneName %in% ref_val$V1)
aw8Lu<-sum(entbic8Lu$GeneName %in% ref_val$V1)
aw7Ru<-sum(entbic7Ru$GeneName %in% ref_val$V1)
aw7Lu<-sum(entbic7Lu$GeneName %in% ref_val$V1)
aw6Ru<-sum(entbic6Ru$GeneName %in% ref_val$V1)
aw6Lu<-sum(entbic6Lu$GeneName %in% ref_val$V1)


#unique hits by name and cleavage site
entbic10uh<-distinct(entbic10,GeneName,Hit,.keep_all = T)
entbic9Ruh<-distinct(entbic9R,GeneName,Hit,.keep_all = T)
entbic9Luh<-distinct(entbic9L,GeneName,Hit,.keep_all = T)
entbic8Ruh<-distinct(entbic8R,GeneName,Hit,.keep_all = T)
entbic8Luh<-distinct(entbic8L,GeneName,Hit,.keep_all = T)
entbic7Ruh<-distinct(entbic7R,GeneName,Hit,.keep_all = T)
entbic7Luh<-distinct(entbic7L,GeneName,Hit,.keep_all = T)
entbic6Ruh<-distinct(entbic6R,GeneName,Hit,.keep_all = T)
entbic6Luh<-distinct(entbic6L,GeneName,Hit,.keep_all = T)

dmsic10uh<-distinct(dmsic10,GeneName,Hit,.keep_all = T)
dmsic9Ruh<-distinct(dmsic9R,GeneName,Hit,.keep_all = T)
dmsic9Luh<-distinct(dmsic9L,GeneName,Hit,.keep_all = T)
dmsic8Ruh<-distinct(dmsic8R,GeneName,Hit,.keep_all = T)
dmsic8Luh<-distinct(dmsic8L,GeneName,Hit,.keep_all = T)
dmsic7Ruh<-distinct(dmsic7R,GeneName,Hit,.keep_all = T)
dmsic7Luh<-distinct(dmsic7L,GeneName,Hit,.keep_all = T)
dmsic6Ruh<-distinct(dmsic6R,GeneName,Hit,.keep_all = T)
dmsic6Luh<-distinct(dmsic6L,GeneName,Hit,.keep_all = T)



aw10uh<-sum(entbic10uh$GeneName %in% ref_val$V1)
aw9Ruh<-sum(entbic9Ruh$GeneName %in% ref_val$V1)
aw9Luh<-sum(entbic9Luh$GeneName %in% ref_val$V1)
aw8Ruh<-sum(entbic8Ruh$GeneName %in% ref_val$V1)
aw8Luh<-sum(entbic8Luh$GeneName %in% ref_val$V1)
aw7Ruh<-sum(entbic7Ruh$GeneName %in% ref_val$V1)
aw7Luh<-sum(entbic7Luh$GeneName %in% ref_val$V1)
aw6Ruh<-sum(entbic6Ruh$GeneName %in% ref_val$V1)
aw6Luh<-sum(entbic6Luh$GeneName %in% ref_val$V1)

dw10uh<-sum(dmsic10uh$GeneName %in% ref_val$V1)
dw9Ruh<-sum(dmsic9Ruh$GeneName %in% ref_val$V1)
dw9Luh<-sum(dmsic9Luh$GeneName %in% ref_val$V1)
dw8Ruh<-sum(dmsic8Ruh$GeneName %in% ref_val$V1)
dw8Luh<-sum(dmsic8Luh$GeneName %in% ref_val$V1)
dw7Ruh<-sum(dmsic7Ruh$GeneName %in% ref_val$V1)
dw7Luh<-sum(dmsic7Luh$GeneName %in% ref_val$V1)
dw6Ruh<-sum(dmsic6Ruh$GeneName %in% ref_val$V1)
dw6Luh<-sum(dmsic6Luh$GeneName %in% ref_val$V1)


df1 <- data.frame(Sequences=rep(c("alignments", "DMS"), each=9),window=c('10aa', '9R', '8R', '7R', '6R', '9L', '8L', '7L', '6L'),
                 cant=c(aw10,aw9R, aw8R, aw7R, aw6R, aw9L, aw8L, aw7L, aw6L,dw10, dw9R,dw8R,dw7R, dw6R, dw9L, dw8L, dw7L, dw6L))

ggplot(df1, aes(window, cant, label=cant)) + ggtitle("Raw hits")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('10aa', '9L', '8L', '7L', '6L','6R','7R', '8R', '9R')) +
 labs(y="Number of Reference Hits", x= "Peptide Window")+
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()

df2<- data.frame(Sequences=rep(c("alignments", "DMS"), each=9),window=c('10aa', '9R', '8R', '7R', '6R', '9L', '8L', '7L', '6L'),
                 cant=c(aw10u,aw9Ru, aw8Ru, aw7Ru, aw6Ru, aw9Lu, aw8Lu, aw7Lu, aw6Lu,dw10u, dw9Ru,dw8Ru,dw7Ru, dw6Ru, dw9Lu, dw8Lu, dw7Lu, dw6Lu))

ggplot(df2, aes(window, cant, label=cant)) + ggtitle("Unique gene names")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('10aa', '9L', '8L', '7L', '6L','6R','7R', '8R', '9R')) +
  labs(y="Number of Reference Hits", x= "Peptide Window")+
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()

df3 <- data.frame(Sequences=rep(c("alignments", "DMS"), each=9),window=c('10aa', '9R', '8R', '7R', '6R', '9L', '8L', '7L', '6L'),
                  cant=c(aw10uh,aw9Ruh, aw8Ruh, aw7Ruh, aw6Ruh, aw9Luh, aw8Luh, aw7Luh, aw6Luh,dw10uh, dw9Ruh,dw8Ruh,dw7Ruh, dw6Ruh, dw9Luh, dw8Luh, dw7Luh, dw6Luh))
                                                                            
ggplot(df3, aes(window, cant, label=cant)) + ggtitle("Unique cleavage sites")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('10aa', '9L', '8L', '7L', '6L','6R','7R', '8R', '9R')) +
  labs(y="Number of Reference Hits", x= "Peptide Window")+
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()                                                                          


ratio_a10  <-aw10/total10_a
ratio_a10u <-aw10u/total10_a
ratio_a10uh<-aw10uh/total10_a
ratio_a9R  <-aw9R/total9R_a
ratio_a9Ru <-aw9Ru/total9R_a
ratio_a9Ruh<-aw9Ruh/total9R_a
ratio_a9L  <-aw9L/total9L_a
ratio_a9Lu <-aw9Lu/total9L_a
ratio_a9Luh<-aw9Luh/total9L_a
ratio_a8R  <-aw8R/total8R_a
ratio_a8Ru <-aw8Ru/total8R_a
ratio_a8Ruh<-aw8Ruh/total8R_a
ratio_a8L  <-aw9L/total8L_a
ratio_a8Lu <-aw9Lu/total8L_a
ratio_a8Luh<-aw9Luh/total8L_a
ratio_a7R  <-aw7R/total7R_a
ratio_a7Ru <-aw7Ru/total7R_a
ratio_a7Ruh<-aw7Ruh/total7R_a
ratio_a7L  <-aw7L/total7L_a
ratio_a7Lu <-aw7Lu/total7L_a
ratio_a7Luh<-aw7Luh/total7L_a
ratio_a6R  <-aw6R/total6R_a
ratio_a6Ru <-aw6Ru/total6R_a
ratio_a6Ruh<-aw6Ruh/total6R_a
ratio_a6L  <-aw6L/total6L_a
ratio_a6Lu <-aw6Lu/total6L_a
ratio_a6Luh<-aw6Luh/total6L_a
  
ratio_d10  <-dw10/total10_d
ratio_d10u <-dw10u/total10_d
ratio_d10uh<-dw10uh/total10_d
ratio_d9R  <-dw9R/total9R_d
ratio_d9Ru <-dw9Ru/total9R_d
ratio_d9Ruh<-dw9Ruh/total9R_d
ratio_d9L  <-dw9L/total9L_d
ratio_d9Lu <-dw9Lu/total9L_d
ratio_d9Luh<-dw9Luh/total9L_d
ratio_d8R  <-dw8R/total8R_d
ratio_d8Ru <-dw8Ru/total8R_d
ratio_d8Ruh<-dw8Ruh/total8R_d
ratio_d8L  <-dw9L/total8L_d
ratio_d8Lu <-dw9Lu/total8L_d
ratio_d8Luh<-dw9Luh/total8L_d
ratio_d7R  <-dw7R/total7R_d
ratio_d7Ru <-dw7Ru/total7R_d
ratio_d7Ruh<-dw7Ruh/total7R_d
ratio_d7L  <-dw7L/total7L_d
ratio_d7Lu <-dw7Lu/total7L_d
ratio_d7Luh<-dw7Luh/total7L_d
ratio_d6R  <-dw6R/total6R_d
ratio_d6Ru <-dw6Ru/total6R_d
ratio_d6Ruh<-dw6Ruh/total6R_d
ratio_d6L  <-dw6L/total6L_d
ratio_d6Lu <-dw6Lu/total6L_d
ratio_d6Luh<-dw6Luh/total6L_d

ratplot<- data.frame(Sequences=rep(c("alignments", "DMS"), each=9),window=c('10aa', '9R', '8R', '7R', '6R', '9L', '8L', '7L', '6L'),
                  rat=c(ratio_a10,ratio_a9R, ratio_a8R, ratio_a7R, ratio_a6R, ratio_a9L, ratio_a8L, ratio_a7L, ratio_a6L,ratio_d10, ratio_d9R,ratio_d8R,ratio_d7R, ratio_d6R, ratio_d9L, ratio_d8L, ratio_d7L, ratio_d6L), cant=c(aw10,aw9R, aw8R, aw7R, aw6R, aw9L, aw8L, aw7L, aw6L,dw10, dw9R,dw8R,dw7R, dw6R, dw9L, dw8L, dw7L, dw6L))

ggplot(ratplot, aes(window, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background (raw)")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('10aa', '9L', '8L', '7L', '6L','6R','7R', '8R', '9R')) +
  labs(y="Ratio", x= "Peptide Window")+
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()

ratplot_u<- data.frame(Sequences=rep(c("alignments", "DMS"), each=9),window=c('10aa', '9R', '8R', '7R', '6R', '9L', '8L', '7L', '6L'),
                       rat=c(ratio_a10u,ratio_a9Ru, ratio_a8Ru, ratio_a7Ru, ratio_a6Ru, ratio_a9Lu, ratio_a8Lu, ratio_a7Lu, ratio_a6Lu,ratio_d10u, ratio_d9Ru,ratio_d8Ru,ratio_d7Ru, ratio_d6Ru, ratio_d9Lu, ratio_d8Lu, ratio_d7Lu, ratio_d6Lu), cant=c(aw10u,aw9Ru, aw8Ru, aw7Ru, aw6Ru, aw9Lu, aw8Lu, aw7Lu, aw6Lu,dw10u, dw9Ru,dw8Ru,dw7Ru, dw6Ru, dw9Lu, dw8Lu, dw7Lu, dw6Lu))

ggplot(ratplot_u, aes(window, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background (unique gene names)")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('10aa', '9L', '8L', '7L', '6L','6R','7R', '8R', '9R')) +
  labs(y="Ratio", x= "Peptide Window")+
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()

ratplot_uh<- data.frame(Sequences=rep(c("alignments", "DMS"), each=9),window=c('10aa', '9R', '8R', '7R', '6R', '9L', '8L', '7L', '6L'),
                        rat=c(ratio_a10uh,ratio_a9Ruh, ratio_a8Ruh, ratio_a7Ruh, ratio_a6Ruh, ratio_a9Luh, ratio_a8Luh, ratio_a7Luh, ratio_a6Luh,ratio_d10uh, ratio_d9Ruh,ratio_d8Ruh,ratio_d7Ruh, ratio_d6Ruh, ratio_d9Luh, ratio_d8Luh, ratio_d7Luh, ratio_d6Luh), cant=c(aw10uh,aw9Ruh, aw8Ruh, aw7Ruh, aw6Ruh, aw9Luh, aw8Luh, aw7Luh, aw6Luh,dw10uh, dw9Ruh,dw8Ruh,dw7Ruh, dw6Ruh, dw9Luh, dw8Luh, dw7Luh, dw6Luh))

ggplot(ratplot_uh, aes(window, rat, cant, label=Sequences)) + ggtitle("Ratio hits-background (unique cleavage sites)")+
  geom_bar(aes(fill = Sequences), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('10aa', '9L', '8L', '7L', '6L','6R','7R', '8R', '9R')) +
  labs(y="Ratio", x= "Peptide Window")+
  geom_text(aes(label= cant, group = Sequences),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()


#joining forces (Alignments and DMS)
both10<-rbind(entbic10,dmsic10,by="ProteinAcc")
both10u<-distinct(both10,GeneName,.keep_all=T)
both10uh<-distinct(both10,GeneName,Hit,.keep_all=T)
both9R<-rbind(entbic9R,dmsic9R,by="ProteinAcc")
both9Ru<-distinct(both9R,GeneName,.keep_all=T)
both9Ruh<-distinct(both9R,GeneName,Hit,.keep_all=T)
both9L<-rbind(entbic9L,dmsic9L,by="ProteinAcc")
both9Lu<-distinct(both9L,GeneName,.keep_all=T)
both9Luh<-distinct(both9L,GeneName,Hit,.keep_all=T)
both8R<-rbind(entbic8R,dmsic8R,by="ProteinAcc")
both8Ru<-distinct(both8R,GeneName,.keep_all=T)
both8Ruh<-distinct(both8R,GeneName,Hit,.keep_all=T)
both8L<-rbind(entbic8L,dmsic8L,by="ProteinAcc")
both8Lu<-distinct(both8L,GeneName,.keep_all=T)
both8Luh<-distinct(both8L,GeneName,Hit,.keep_all=T)
both7R<-rbind(entbic7R,dmsic7R,by="ProteinAcc")
both7Ru<-distinct(both7R,GeneName,.keep_all=T)
both7Ruh<-distinct(both7R,GeneName,Hit,.keep_all=T)
both7L<-rbind(entbic7L,dmsic7L,by="ProteinAcc")
both7Lu<-distinct(both7L,GeneName,.keep_all=T)
both7Luh<-distinct(both7L,GeneName,Hit,.keep_all=T)
both6R<-rbind(entbic6R,dmsic6R,by="ProteinAcc")
both6Ru<-distinct(both6R,GeneName,.keep_all=T)
both6Ruh<-distinct(both6R,GeneName,Hit,.keep_all=T)
both6L<-rbind(entbic6L,dmsic6L,by="ProteinAcc")
both6Lu<-distinct(both6L,GeneName,.keep_all=T)
both6Luh<-distinct(both6L,GeneName,Hit,.keep_all=T)

total10b<-nrow(both10)
total10bu<-nrow(both10u)
total10buh<-nrow(both10uh)
total9Rb<-nrow(both9R)
total9Rbu<-nrow(both9Ru)
total9Rbuh<-nrow(both9Ruh)
total9Lb<-nrow(both9L)
total9Lbu<-nrow(both9Lu)
total9Lbuh<-nrow(both9Luh)
total8Rb<-nrow(both8R)
total8Rbu<-nrow(both8Ru)
total8Rbuh<-nrow(both8Ruh)
total8Lb<-nrow(both8L)
total8Lbu<-nrow(both8Lu)
total8Lbuh<-nrow(both8Luh)
total7Rb<-nrow(both7R)
total7Rbu<-nrow(both7Ru)
total7Rbuh<-nrow(both7Ruh)
total7Lb<-nrow(both7L)
total7Lbu<-nrow(both7Lu)
total7Lbuh<-nrow(both7Luh)
total6Rb<-nrow(both6R)
total6Rbu<-nrow(both6Ru)
total6Rbuh<-nrow(both6Ruh)
total6Lb<-nrow(both6L)
total6Lbu<-nrow(both6Lu)
total6Lbuh<-nrow(both6Luh)


b10<-sum(both10$GeneName %in% ref_val$V1)
b10u<-sum(both10u$GeneName %in% ref_val$V1)
b10uh<-sum(both10uh$GeneName %in% ref_val$V1)
b9R<-sum(both9R$GeneName %in% ref_val$V1)
b9Ru<-sum(both9Ru$GeneName %in% ref_val$V1)
b9Ruh<-sum(both9Ruh$GeneName %in% ref_val$V1)
b9L<-sum(both9L$GeneName %in% ref_val$V1)
b9Lu<-sum(both9Lu$GeneName %in% ref_val$V1)
b9Luh<-sum(both9Luh$GeneName %in% ref_val$V1)
b8R<-sum(both8R$GeneName %in% ref_val$V1)
b8Ru<-sum(both8Ru$GeneName %in% ref_val$V1)
b8Ruh<-sum(both8Ruh$GeneName %in% ref_val$V1)
b8L<-sum(both8L$GeneName %in% ref_val$V1)
b8Lu<-sum(both8Lu$GeneName %in% ref_val$V1)
b8Luh<-sum(both8Luh$GeneName %in% ref_val$V1)
b7R<-sum(both7R$GeneName %in% ref_val$V1)
b7Ru<-sum(both7Ru$GeneName %in% ref_val$V1)
b7Ruh<-sum(both7Ruh$GeneName %in% ref_val$V1)
b7L<-sum(both7L$GeneName %in% ref_val$V1)
b7Lu<-sum(both7Lu$GeneName %in% ref_val$V1)
b7Luh<-sum(both7Luh$GeneName %in% ref_val$V1)
b6R<-sum(both6R$GeneName %in% ref_val$V1)
b6Ru<-sum(both6Ru$GeneName %in% ref_val$V1)
b6Ruh<-sum(both6Ruh$GeneName %in% ref_val$V1)
b6L<-sum(both6L$GeneName %in% ref_val$V1)
b6Lu<-sum(both6Lu$GeneName %in% ref_val$V1)
b6Luh<-sum(both6Luh$GeneName %in% ref_val$V1)


ratio_b10  <-b10/total10b
ratio_b10u <-b10u/total10b
ratio_b10uh<-b10uh/total10b
ratio_b9R  <-b9R/total9Rb
ratio_b9Ru <-b9Ru/total9Rb
ratio_b9Ruh<-b9Ruh/total9Rb
ratio_b9L  <-b9L/total9Lb
ratio_b9Lu <-b9Lu/total9Lb
ratio_b9Luh<-b9Luh/total9Lb
ratio_b8R  <-b8R/total8Rb
ratio_b8Ru <-b8Ru/total8Rb
ratio_b8Ruh<-b8Ruh/total8Rb
ratio_b8L  <-b9L/total8Lb
ratio_b8Lu <-b9Lu/total8Lb
ratio_b8Luh<-b9Luh/total8Lb
ratio_b7R  <-b7R/total7Rb
ratio_b7Ru <-b7Ru/total7Rb
ratio_b7Ruh<-b7Ruh/total7Rb
ratio_b7L  <-b7L/total7Lb
ratio_b7Lu <-b7Lu/total7Lb
ratio_b7Luh<-b7Luh/total7Lb
ratio_b6R  <-b6R/total6Rb
ratio_b6Ru <-b6Ru/total6Rb
ratio_b6Ruh<-b6Ruh/total6Rb
ratio_b6L  <-b6L/total6Lb
ratio_b6Lu <-b6Lu/total6Lb
ratio_b6Luh<-b6Luh/total6Lb


ratplotb<- data.frame(Source=rep(c("raw", "unique gene names", "unique cleavage sites"), each=9),window=c('10aa', '9R', '8R', '7R', '6R', '9L', '8L', '7L', '6L'),
                     rat=c(ratio_b10,ratio_b9R, ratio_b8R, ratio_b7R, ratio_b6R, ratio_b9L, ratio_b8L, ratio_b7L, ratio_b6L, 
                           ratio_b10u,ratio_b9Ru, ratio_b8Ru, ratio_b7Ru, ratio_b6Ru, ratio_b9Lu, ratio_b8Lu, ratio_b7Lu, ratio_b6Lu, 
                           ratio_b10uh,ratio_b9Ruh, ratio_b8Ruh, ratio_b7Ruh, ratio_b6Ruh, ratio_b9Luh, ratio_b8Luh, ratio_b7Luh, ratio_b6Luh), 
                            cant=c(b10,b9R, b8R, b7R, b6R, b9L, b8L, b7L, b6L,b10u,b9Ru, b8Ru, b7Ru, b6Ru, b9Lu, b8Lu, b7Lu, b6Lu,
                                   b10uh,b9Ruh, b8Ruh, b7Ruh, b6Ruh, b9Luh, b8Luh, b7Luh, b6Luh ))

ggplot(ratplotb, aes(window, rat, cant, label=Source)) + ggtitle("Ratio hits-background (Alignments-DMS)")+
  geom_bar(aes(fill = Source), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('10aa', '9L', '8L', '7L', '6L','6R','7R', '8R', '9R')) +
  labs(y="Ratio", x= "Peptide Window")+
  geom_text(aes(label= cant, group = Source),color = "white", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal(base_size = 19)



common8R <- both8Ruh[both8Ruh$GeneName %in% ref_val$V1,]


totals <- data.frame(Source=rep(c("alignments", "DMS", "both", "both unique genes", "both unique sites"), each=9),window=c('10aa','9L', '8L', '7L', '6L', '6R','7R', '8R', '9R' ),
                  cant=c(total10_a, total9L_a, total8L_a, total7L_a, total6L_a, total6R_a, total7R_a, total8R_a, total9R_a, total10_d, total9L_d, total8L_d, total7L_d, total6L_d, total6R_d, total7R_d, total8R_d, total9R_d, total10b, total9Lb, total8Lb, total7Lb, total6Lb, total6Rb, total7Rb, total8Rb, total9Rb, total10bu, total9Lbu, total8Lbu, total7Lbu, total6Lbu, total6Rbu, total7Rbu, total8Rbu, total9Rbu, total10buh,total9Lbuh, total8Lbuh, total7Lbuh, total6Lbuh, total6Rbuh, total7Rbuh, total8Rbuh, total9Rbuh )
                  )
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               

ggplot(totals, aes(window, cant, label=Source)) + ggtitle("Number of proteins PSSMSearch output")+
  geom_bar(aes(fill = Source), stat = "identity", position = "dodge")+
  scale_x_discrete(limits = c('10aa', '9L', '8L', '7L', '6L','6R','7R', '8R', '9R')) +
  labs(y="Number of hits", x= "Peptide Window")+
  geom_text(aes(label= cant, group = Source),color = "black", position = position_dodge2(width = .8), show.legend = FALSE, hjust = .6,vjust=4, size = 4)+ 
  theme_minimal()



write.table(both8R, file="3C/8Rcombination.txt", quote=F, sep= "\t", row.names = F)
