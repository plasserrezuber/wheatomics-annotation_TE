library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

#localisation donnees entree
dirin<-"Y:/ANALYSES/results/clariTE/"
setwd(dirin)

classi=read_tsv("TE_classification.txt", col_names=T)

d=read_tsv("Length_TE_Renan_SUPERFamLevel_Renan.tsv", col_names=F)
names(d)<-c("chrom", "fam", "length")
ls.str(d)

d=d%>%mutate(genome=str_extract(chrom, "[:upper:]"), chrom=gsub("chr", "", chrom))

d=left_join(d, classi, by="fam")

total_dna=14195643615
d_genome=d%>%group_by(fam)%>%summarise(cumul_length=sum(length), percent_DNA=round(sum(length)/total_dna*100,2))

d_subgenome=d%>%group_by(genome,fullname, fam)%>%summarise(cumul_length=sum(length), percent_DNA=round(sum(length)/total_dna*100,2))
d_subgenome=d_subgenome%>%mutate(fullname=paste(fullname, " (", name, ")", sep=""))

### RENAN sub genome
g=ggplot(d_subgenome[!is.na(d_subgenome$fullname),], aes(x=genome, y=percent_DNA, fill=fullname))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette="Paired")
plot(g)

### RENAN sub genome
g2=ggplot(d_subgenome[!is.na(d_subgenome$fullname),], aes(x=genome, y=cumul_length, fill=fullname))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette="Paired")
plot(g2)


#### COMPARAISON RENAN REFSEQV2
D=read_tsv("TE_SUPERfam_length_REFSEQV2_RENAN.txt", col_names=T)
D=left_join(D, classi, by="fam")



############################################################################
d=read_tsv("Nb_TE_Renan_SUPERFamLevel_Renan.tsv", col_names=F)
names(d)<-c("chrom", "fam", "Nb")


