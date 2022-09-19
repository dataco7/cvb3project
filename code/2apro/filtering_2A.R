library(tidyverse)
library(readr)
library(stringi)

max.col= max(count.fields("access_windowed_modified.txt", sep = ','))

dasa = read.table("access_windowed_modified.txt", header = FALSE, sep = ",", 
                  col.names = paste0("V",seq_len(max.col)), fill = TRUE)

dss = read.table("secondary_structure_modified.txt", header = FALSE, sep = ",", 
                 col.names = paste0("V",seq_len(max.col)), fill = TRUE)

pssm <- read_delim("2A/9Rcombination.txt", delim= "\t")
len<- read_tsv("3C/proteomehuman_length.txt")
validated2a <- read.delim("2A/provalidated_2a.txt", header= F )
reactome <- read.delim("cap_trans_reactome.tsv") %>% select(c("Identifier")) %>% rename(c("ProteinID"="Identifier"))
intf <-read.delim("interferome_results.txt",skip = 19) %>% filter(intf$Fold.Change >=2 )
ifn <- read.delim("3C/filesbase/3c_interferome.txt", skip=19, header=T)

pssm <- pssm[!grepl("Secreted|Lumenal|Extracellular",pssm$Warnings),]

pssm$Length <- NA

data_p <- pssm %>% select(ProteinAcc, GeneName, Length, SeqStart, SeqStop,  Hit)
data_p$SeqStart <- as.numeric(data_p$SeqStart)
data_p$SeqStop <- as.numeric(data_p$SeqStop)

sum(pssm$GeneName %in% validated2a$V1) #71

for (i in 1:nrow(data_p)){
  # i=1
  if (data_p$ProteinAcc[i] %in% len$Entry){
    ic= data_p$ProteinAcc[i]
    v = grep(data_p$ProteinAcc[i], len$Entry)
    long = len$Length[v]
    data_p$Length[i] <- long
  }
}

sum(is.na(data_p$Length)) #cuantos NA tenemos en la col Length
which(is.na(pssm$Length))#vamos a buscar las long y asignarlas a mano




data_p$windowasa <- NA
data_p$windowss <- NA



data_p_clean <- subset(data_p,!is.na(data_p$Length))



for (i in 1:nrow(data_p_clean)){
  
  #i=1
  
  name = data_p_clean$ProteinAcc[i]
  st = data_p_clean$SeqStart[i]
  end = data_p_clean$SeqStop[i]
  long = data_p_clean$Length[i]
  corrend = long - 9
  
  
  
  if (st > 9 & end < corrend){
    
    temp_asa <-dasa[dasa$V1==name,c((st-9):(end+11))]
    temp_ss <- dss[dss$V1==name,c((st-9):(end+11))]
    
    
  } 
  if (st< 9){
    
    temp_asa <-dasa[dasa$V1==name,c((st+1):(end+11))]
    temp_ss <- dss[dss$V1==name,c((st+1):(end+11))]
    
  }
  
  
  
  if (st< 9 & end>corrend){
    
    temp_asa <-dasa[dasa$V1==name,c((st+1):(end))]
    temp_ss <- dss[dss$V1==name,c((st+1):(end))]
    
    
  }
  
  #if (end>corrend){
  
  #temp_asa <-dasa[dasa$V1==name,c((st-9):(end))]
  #temp_ss <- dss[dss$V1==name,c((st-9):(end))]
  
  
  #}
  
  data_p_clean$windowasa[i] <- paste(temp_asa[,],collapse = ",")
  data_p_clean$windowss[i] <-  paste(temp_ss[,],collapse = ",")
}

write.table(data_p_clean, "2A/filtering/preresults_unfiltered.csv", quote= F, row.names= F, sep= " ")


##FILTERING
data_raw <- data_p_clean
sum(data_raw$GeneName %in% validated2a$V1) #71

data_raw$exp <- NA
data_raw$SSe <- NA
data_raw$SSl <- NA
data_raw$SSh <- NA

data_p <- data_raw

for (i in 1:nrow(data_p)){
  
  #i=1
  
  tot = lengths(gregexpr(",", data_p$windowasa[i])) + 1
  asa = as.numeric(unlist(stri_split(data_p$windowasa[i],fixed=',')))
  
  data_p$exp[i] <- ((sum(asa>0.5))*100)/tot
  
  data_p$SSe[i] <- str_count(data_p$windowss[i],pattern="E")
  data_p$SSl[i] <- str_count(data_p$windowss[i],pattern="-")
  data_p$SSh[i] <- str_count(data_p$windowss[i],pattern="H|I|G")
  
}

#checkpoint prefiltering

sum(data_p$GeneName %in% ifn$Gene.Name)
sum(data_p$GeneName %in% validated3c$V1)
sum(data_p$ProteinAcc %in% reactome$ProteinID)

sel<- data.frame()
sel=data_p[data_p$exp>=100,]
sel<- sel[sel$SSe==0 | sel$SSl> 30 | sel$SSh> 20,]
sum(sel$GeneName %in% validated2a$V1)
sum(data_p$GeneName %in% ifn$Gene.Name)

cellexpr<-read_tsv('proteinatlas_757cd1cc.tsv')

virusprot <- cellexpr[grepl('Host-virus interaction',cellexpr$`Biological process`),]
sum(sel$GeneName %in% virusprot$`Gene`)
viral <- sel[sel$GeneName %in% virusprot$Gene,]

sel_int1 <- sel[sel$GeneName %in% ifn$Gene.Name,]
#or intf
sel_int2 <- sel[sel$GeneName %in% intf$Gene.Name,]
sel_react <- sel[sel$ProteinAcc %in% reactome$ProteinID,]

sel_join<- rbind(sel_int1,sel_react, viral)
sel_join <- distinct(sel_join, GeneName, Hit, .keep_all = T)

sum(sel_join$GeneName %in% validated3c$V1) #44

#sel_join <- pssm[pssm$GeneName %in% sel_join$GeneName,]
#protset=sel_join[!grepl("Secreted|Lumenal|Extracellular",sel_join$Warnings),

protset <- pssm[pssm$GeneName %in% sel_join$GeneName,]
sum(protset$GeneName %in% validated3c$V1)

unique_2a <- distinct(sel_join, GeneName,.keep_all=T)
sum(unique_2a$GeneName %in% validated2a$V1)
validated2a <- unique(validated2a)



selrank <- protset[protset$Rank<= 3000,]
sum(selrank$GeneName %in% validated3c$V1)
#selhitG <- selrank[grepl("QG", selrank$Hit),]
#selhitN <- selrank[grepl("QN", selrank$Hit),]

#hitsjoin <- rbind(selhitG, selhitN)

#sum(selhitG$GeneName %in% validated3c$V1)
#sum(selhitN$GeneName %in% validated3c$V1)

unique <- distinct(protset, GeneName)


write.table(sel_join, file="2A/filtering/filteredset_2Apro.txt", quote=F, row.names = F)




#Analysis GO

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# load the library
library("biomaRt")

# I prefer ensembl so that the one I will query, but you can
# query other bases, try out: listMarts() 
ensembl=useMart("ensembl")

# as it seems that you are looking for human genes:
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# if you want other model organisms have a look at:
#listDatasets(ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

goids = getBM(
  
  #you want entrezgene so you know which is what, the GO ID and
  # name_1006 is actually the identifier of 'Go term name'
  attributes=c('entrezgene','go_id', 'name_1006'), 
  
  filters='entrezgene', 
  values=hitsjoin$ProteinAcc, 
  mart=ensembl)
