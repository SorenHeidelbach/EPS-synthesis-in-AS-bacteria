library(ape)
library(RFLPtools)
library(stringr)
library(tidyverse)
library(dplyr)
library(reshape2)
library(gplots)
library(seqinr)
library(plyr)
library(writexl)
library(varhandle)
require('Gviz')
library(data.table)
require(stats)
library(gridExtra)
library(seqinr)
library(ggplot2)
library(hrbrthemes)
library(extrafont)
library(forcats)
loadfonts(device = "win")

#Set working directory and path to psiblast input
setwd("C:/Users/heide/OneDrive - Aalborg Universitet/BIOT/8. semester/Project/Data")

#Define parameters
Prokkadistance = 8
Intergenetic = 2000
Perc_genes = 0.5
Perc_id = 25
upstream = 2
downstream = 2

fasta_export = FALSE
psiblastoutput = TRUE

magstats<-read_tsv("magstats.tsv")
gff<-read_tsv("gff.tsv")
gff_extend<-gff
gff_extend$Target_label<-paste(gff_extend$ID, gff_extend$ProkkaNO, sep="_")


#### Loading PSI-BLAST data ####
psiblast_path = "psiblast/curdlan.txt"
title = "curdlan"
psiblast = read.csv(file = psiblast_path, header = FALSE, 
                    sep = "\t", col.names = c("Query_label", "Target_label", "Percent_identity", "align_length", "mismatch", 
                                              "gap_opens", "start_query", "end_query", "start_target", 
                                              "end_target", "E_value", "Bit score"))

#Separating columns
psiblast <- psiblast %>%
  separate(Target_label, c("ID", "ProkkaNO"), sep = -5, remove = FALSE, extra = "merge") %>%
  separate(ID, c("ID", "rm"), sep =-1, remove = TRUE) %>%
  filter(Percent_identity>Perc_id)
psiblast = psiblast[ , !(names(psiblast) %in% "rm")]

#Removing redundant Target_labels (genes of MAG-IDs), keeping highest Bit.score
psiblast <- psiblast[order(psiblast$Target_label, -abs(psiblast$Bit.score) ), ] 
psiblast <- psiblast[ !duplicated(psiblast$Target_label), ] 

#### Joining psiblast dataframe with data from gff data and magstats data ####
psiblast <- psiblast %>% 
  left_join(gff, by = c("ID", "ProkkaNO"), keep = FALSE) %>%
  left_join(magstats, by = "ID", keep = FALSE)



#### Using column ID and Query_label to generate levels in the dataframe #### 
IDvec = unique(psiblast$ID)
Queryvec = unique(psiblast$Query_label)
Queryvec <- sort(Queryvec) 

psiblast <- psiblast %>%
  mutate(ID = factor(ID, levels = IDvec)) %>%
  mutate(Query_label = factor(Query_label, levels = Queryvec))

#Resets empty dataframes, run before looping 
psiblast_subset <- as.data.frame(matrix(0, ncol = ncol(psiblast), nrow = 0))

for (i in seq_along(levels(psiblast$ID))){  ##For each MAG ID
  #i
  df_i <- filter(psiblast, ID == levels(psiblast$ID)[i])
  contigsd = unique(df_i$seqname)
  
  df_i <- df_i %>% 
    mutate(seqname = factor(seqname, levels = contigsd))
  
  
  for (j in seq_along(levels(df_i$seqname))){   ##For each contig within the MAG ID
    df_j <- df_i 
    df_j <- df_j %>% filter(seqname == levels(df_j$seqname)[j]) ##Subsets to genes in one contig 
    df_j <- df_j %>% mutate(OperonNO = c(replicate(nrow(df_j), 0)))
    df_j <- df_j[order(df_j$ProkkaNO),]
    
    
    for (k in 1:nrow(df_j)){          ##For each gene within the contig
      if (k == 1){
        df_j$OperonNO[k] = 1          ##Count the putative operons in a contig, and number them
      } else {
        if ((as.numeric(df_j$ProkkaNO[[k]])-as.numeric(df_j$ProkkaNO[[k-1]])) < Prokkadistance 
            | (as.numeric(df_j$start[[k]]) - as.numeric(df_j$end[[k-1]])) < Intergenetic 
        ){
          df_j$OperonNO[k] = as.numeric(df_j$OperonNO[[k-1]])}
        else {df_j$OperonNO[k] = (df_j$OperonNO[[k-1]] + 1)}
      }
      
    } #k

    operons = unique(df_j$OperonNO)
    df_j <- df_j %>% mutate(OperonNO = factor(OperonNO, levels = operons))
    for (l in seq_along(levels(df_j$OperonNO))){    #For each putative operon within the contig
      df_l <- df_j
      df_l <- df_l %>% filter(OperonNO == levels(df_l$OperonNO)[l]) 
      heatmapdata <- matrix(0, nrow = 1, ncol = length(Queryvec))
      if (length(unique(df_l$Query_label)) >= Perc_genes*length(Queryvec)) {
        psiblast_subset <- rbind(psiblast_subset, df_l)
      }
      
    }#l
    
  }#j
  
}#i

psiblast_alt_subset<-psiblast_subset
if (psiblastoutput){
write_tsv(psiblast_alt_subset, paste("psiblast_subset/", title, ".tsv", sep=""))}

if (fasta_export == TRUE){
  ###Making fasta files for Interproscan and BLAST
  fasta_import<-unique(psiblast_alt_subset$ID)
  dir.create(paste("fasta_output/", title, sep=""), showWarnings = FALSE)
  for (f in 1:(length(fasta_import))){
    fastafile<-read.fasta(file=paste(paste("MGP1000_HQMAG1083_prot_db_split/", fasta_import[f], sep=""), ".faa", sep=""), seqtype="AA", as.string=TRUE, set.attributes=FALSE)
    psiblast_fasta_subset<-psiblast_alt_subset %>% filter(ID == fasta_import[f])
    Prokka_subset<-as.numeric(psiblast_fasta_subset$ProkkaNO)
    Prokka_subset<-(min(Prokka_subset)-2):(max(Prokka_subset+2))
    Prokka_subset<-str_pad(Prokka_subset, 5, pad = "0")
    subset_vector<-paste(fasta_import[f], Prokka_subset, sep="_")
    fasta_subset<-fastafile[names(fastafile) %in% subset_vector]
    write.fasta(fasta_subset, names=names(fasta_subset), file.out=paste(paste(paste(paste("fasta_output/", title, sep=""), title, sep="/"), fasta_import[f], sep = "_"), ".faa", sep=""))
  }
}





###Getting information on surrounding genes from gff file, identify contig borders
###and annotating multiple operons within a MAG###
psiblast_extended<-psiblast_alt_subset
IDs<-unique(psiblast_extended$ID)
psiblast_extended <- psiblast_extended %>% mutate(ID = factor(ID, levels = IDs))
psiblast_surr<-psiblast_extended[0,]
letters<-c("A", "B", "C", "D", "E", "F", "G", "H", "I")
gff_split <- split(gff, f = gff$ID)
for (o in 1:length(unique(psiblast_extended$ID))){
  gff_subset <- gff_split[[unique(psiblast_extended$ID)[o]]] %>% subset(ProkkaNO!="units")
  count<-0
  df_o <- filter(psiblast_extended, ID == levels(psiblast_extended$ID)[o])
  contigs2<-unique(df_o$seqname)
  df_o<-df_o %>% mutate(seqname = factor(seqname, levels = contigs2))
  for (t in 1:length(unique(df_o$seqname))){
  df_t<-filter(df_o, seqname == levels(df_o$seqname)[t])
  operon1<-unique(df_t$OperonNO)
  df_t<-df_t %>% mutate(OperonNO = factor(OperonNO, levels = operon1))
  for (u in 1:length(unique(df_t$OperonNO))){
    count<-count+1
    df_u<- filter(df_t, OperonNO == levels(df_t$OperonNO)[u])
  Prokka_extend<-(min(as.numeric(df_u$ProkkaNO))-upstream):(max(as.numeric(df_u$ProkkaNO))+downstream)
  Prokka_extend<-str_pad(Prokka_extend, 5, pad = "0")
  subset_vector<-paste(IDs[o], Prokka_extend, sep="_")
  gff_extend_subset<-gff_extend %>% filter(Target_label %in% subset_vector)
  
  if (length(unique(df_o$OperonNO))>1 | length(unique(df_o$seqname))>1){
    df_u$ID2<-paste(df_u$ID, letters[count], sep="_")
    gff_extend_subset$ID2<-paste(gff_extend_subset$ID, letters[count], sep="_")
  }else{df_u$ID2<-df_u$ID
        gff_extend_subset$ID2<-gff_extend_subset$ID}
  min_subset<-filter(gff_extend_subset, ProkkaNO==str_pad(as.character(min(as.numeric(df_u$ProkkaNO))), 5, pad="0"))
  if (min(as.numeric(df_u$ProkkaNO)) == 1) {
    min_compare <- min_subset
  } else {
    min_compare<-filter(gff_extend_subset, ProkkaNO==str_pad(as.character(min(as.numeric(df_u$ProkkaNO))-1), 5, pad="0"))
  }
  max_subset<-filter(gff_extend_subset, ProkkaNO==str_pad(as.character(max(as.numeric(df_u$ProkkaNO))), 5, pad="0"))
  if (max(as.numeric(df_u$ProkkaNO)) == max(as.numeric(gff_subset$ProkkaNO))) {
    max_compare <- max_subset
  } else {
    max_compare<-filter(gff_extend_subset, ProkkaNO==str_pad(as.character(max(as.numeric(df_u$ProkkaNO))+1), 5, pad="0"))
  }
  if(max_subset$seqname != max_compare$seqname | max_subset$Target_label == max_compare$Target_label){
    max_bind<-max_subset
    max_bind$Query_label<-"CB"
    max_bind$start<-max_subset$end + 100
    max_bind$end<-max_subset$end + 250
    max_bind$strand<-"*"
    gff_extend_subset<-bind_rows(gff_extend_subset, max_bind)
  }
  
  if(min_subset$seqname != min_compare$seqname | min_subset$Target_label == min_compare$Target_label){
    min_bind<-min_subset
    min_bind$Query_label<-"CB"
    min_bind$start<-min_subset$start - 250
    min_bind$end<-min_subset$start - 100
    min_bind$strand<-"*"
    gff_extend_subset<-bind_rows(gff_extend_subset, min_bind)
  }
  gff_extend_subset<-gff_extend_subset %>% filter(seqname == df_u$seqname[1]) 
  gff_extend_subset<-bind_rows(df_u, gff_extend_subset)
  gff_extend_subset<-gff_extend_subset[order(gff_extend_subset$ProkkaNO, gff_extend_subset$Query_label ), ] 
  gff_extend_subset<-gff_extend_subset[ !duplicated(gff_extend_subset$Target_label), ] 
  psiblast_surr<-bind_rows(psiblast_surr, gff_extend_subset)
  }
}}


#### Making heatmap dataframe####
psiblast_plot<-psiblast_surr %>% drop_na(Query_label) %>% 
  filter(Query_label != "CB") %>%
  mutate(ID2=factor(ID2,levels=rev(sort(unique(ID2)))))

psiblast_plot<-psiblast_plot[ !duplicated(psiblast_plot[c("Query_label", "ID2")]), ] 







#### Generate a list of annotationtracks for plotting ####
op_list <- list()
for (u in unique(psiblast_surr$ID2)){
  opvisdata <- psiblast_surr %>% filter(ID2 ==u)
  minstart <- min(opvisdata$start)
  opvisdata$start <- opvisdata$start - minstart
  opvisdata$end <- opvisdata$end - minstart
  opvisdata$width <- opvisdata$end - opvisdata$start
  op_list[u] <-AnnotationTrack(start = opvisdata$start, 
                               width = opvisdata$width, 
                               strand= opvisdata$strand, 
                               id    = opvisdata$Query_label,
                               group = opvisdata$Target_label,
                               mergeGroups = FALSE,
                               feature = opvisdata$Query_label,
                               name  = as.character(opvisdata$ID2)[1],
                               background.title="azure4",
                               stacking = "dense",
                               chromosome="chrX")
  
}


ref=GenomeAxisTrack(GRanges('TEST',IRanges(start =seq(0, 6000, by=200),
                                           end = c(seq(200, 6200, by=200)))),
                    lwd=20, fontsize=9)



#### Making a dataframe to the pdf ####
psiblast_alt_subset_bind <- psiblast_alt_subset[order(psiblast_alt_subset$Query_label, -abs(psiblast_alt_subset$start) ), ] 
psiblast_alt_subset_bind <- psiblast_alt_subset_bind[ !duplicated(psiblast_alt_subset_bind$ID), ]
df_pdf <- psiblast_alt_subset_bind %>% select(ID, GTDBTax, seqname, start, Query_label, ProkkaNO) %>% 
  arrange(ID)

n_pages_operon<-ceiling(length(op_list)/8)
n_pages_heatmap<-ceiling(length(unique(psiblast_surr$ID2))/15)

#### Make a pdf with the heatmap and operon structures ####
pdf(paste("r_output/", paste(title, "_1.pdf", sep=""), sep=""), width = 15, height = 9)
##Plot of heatmap data
for(p in 1:n_pages_heatmap){
  if (p == n_pages_heatmap){
    m<-15-length(unique(psiblast_surr$ID2))%%15
  } else {m<-0}
  subset_vector1<-unique(psiblast_plot$ID2)[(p*15-14):(p*15-m)]
  psiblast_plot1<-psiblast_plot %>% filter(ID2 %in% subset_vector1)
  plot<-ggplot(psiblast_plot1, aes(Query_label, ID2, fill = Percent_identity)) +
    geom_tile(colour="white",size=0.2)+
    geom_text(aes(label = round(Percent_identity, 0)), cex = 2)+
    scale_fill_gradient(low="lightblue", high="red", name="Percent Identity") +
    xlab("") +
    ylab("MAG")+
    ggtitle(paste(title, "Heatmap", sep=" ")) +
    theme(plot.title = element_text(size = 10, face = "bold"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    coord_fixed(ratio = 1L)+
    theme(legend.title = element_text(size=9, 
                                      face="bold"), legend.text = element_text(size=8))
  
  print(plot)
}


for (p in 1:n_pages_operon){
  if (p == n_pages_operon){
    m<-8-length(op_list)%%8
  } else {m<-0}
  plotTracks(c(op_list[(p*8-7):(p*8-m)], ref),
           showFeatureId=TRUE,
           #stacking = "dense",
           collapse=T, 
           fontcolor.item='black',
           fontcolor.group ="grey1",
           fontsize.group=25,
           title.width = 8.4,
           cex.main = 1,
           main = title,
           arrowHeadWidth= 10,
           arrowHeadMaxWidth=20,
           cex=0.3,
           cex.title=0.7,
           lty=1,
           lwd=1,
           GT="gray30",
           pgaA="#C8DBC8",
           pgaB="#82ACB9",
           pgaC="#BAB4D8",
           epsK="#B9DBF4",
           CB="red",
           rotation.title=1,
           trackLabelHAlign = "center")
}

grid.newpage()
theme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.6)),
  colhead = list(fg_params=list(cex = 0.6)),
  rowhead = list(fg_params=list(cex = 0.6)),
  base_size = 6)
a <- tableGrob(df_pdf, cols=colnames(df_pdf), theme = theme)
grid.draw(a)

dev.off(which = dev.cur())
#### ####
