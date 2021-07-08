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

#Set working directory and path to psiblast input
setwd("C:/Users/heide/OneDrive - Aalborg Universitet/BIOT/8. semester/Project/in_silico/Data")

#Define parameters
Prokkadistance = 3
Intergenetic = 2000
Perc_genes = 1
Perc_id = 20
upstream = 0
downstream = 0

fasta_export = TRUE

magstats <- read_tsv("magstats.tsv")
gff <- read_tsv("gff.tsv")
gff_extend <- gff
gff_extend$Target_label <- paste(gff_extend$ID, gff_extend$ProkkaNO, sep="_")


#### Loading PSI-BLAST data ####
title = "cellulose1"
psiblast_path = paste("psiblast/",title,".txt",sep="")
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
heatmapdf <- as.data.frame(matrix(0, ncol = length(Queryvec)+1, nrow = 0)) 
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
            & (as.numeric(df_j$start[[k]]) - as.numeric(df_j$end[[k-1]])) < Intergenetic 
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
        IDdf_l = unfactor(df_l$ID[1])
        heatmapdata[1] = IDdf_l
        psiblast_subset <- rbind(psiblast_subset, df_l)
        for (p in 1:length(Queryvec)) {             ##For each gene in the query, search the putative operon
          df_p <- df_l %>% filter(Query_label == Queryvec[p])
          if (nrow(df_p) > 0) {
            heatmapdata[p+1] = df_p$Percent_identity[1]    #If the query gene is present in an operon, print Percent_identity to a dataframe
          } else {
            heatmapdata[p+1] = NA
          }
        }
        heatmapdf = rbind(heatmapdf, heatmapdata)
        
      }
      
    }#l
    
  }#j
  
}#i

psiblast_alt_subset<-psiblast_subset
write.csv(psiblast_alt_subset, file = paste(title,".csv",sep=""))

if (fasta_export == TRUE){
  ###Making fasta files for Interproscan and BLAST
  fasta_import<-unique(psiblast_alt_subset$ID)
  dir.create(paste("fasta_output/", title, sep=""))
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



#### 

#Run till here, now run Interproscan on outputted FASTA files

####



###Getting information on surrounding genes from gff file
psiblast_extended<-psiblast_alt_subset
IDs<-unique(psiblast_extended$ID)
psiblast_extended <- psiblast_extended %>% mutate(ID = factor(ID, levels = IDs))
psiblast_surr<-psiblast_extended[0,]
for (o in 1:length(unique(psiblast_extended$ID))){
  df_o <- filter(psiblast_extended, ID == levels(psiblast_extended$ID)[o])
  Prokka_extend<-(min(as.numeric(df_o$ProkkaNO))-upstream):(max(as.numeric(df_o$ProkkaNO))+downstream)
  Prokka_extend<-str_pad(Prokka_extend, 5, pad = "0")
  subset_vector<-paste(IDs[o], Prokka_extend, sep="_")
  gff_extend_subset<-gff_extend %>% filter(Target_label %in% subset_vector) 
  gff_extend_subset<-bind_rows(df_o, gff_extend_subset)
  gff_extend_subset<-gff_extend_subset[order(gff_extend_subset$ProkkaNO, gff_extend_subset$Query_label ), ] 
  gff_extend_subset<-gff_extend_subset[ !duplicated(gff_extend_subset$Target_label), ] 
  psiblast_surr<-rbind(psiblast_surr, gff_extend_subset)
}

### Import tsv files for domain annotation



#### Making a dataframe to the pdf ####
psiblast_alt_subset_bind <- psiblast_alt_subset[order(psiblast_alt_subset$Query_label, -abs(psiblast_alt_subset$start) ), ] 
psiblast_alt_subset_bind <- psiblast_alt_subset_bind[ !duplicated(psiblast_alt_subset_bind$ID), ]
df_pdf <- psiblast_alt_subset_bind %>% select(ID, GTDBTax, seqname, start, Query_label, ProkkaNO) %>% 
  arrange(ID)

#### Generate a list of annotationtracks for plotting ####
op_list <- list()
for (u in unique(psiblast_alt_subset$ID)){
  opvisdata <- psiblast_surr %>% filter(ID ==u)
  minstart <- min(opvisdata$start)
  opvisdata$start <- opvisdata$start - minstart
  opvisdata$end <- opvisdata$end - minstart
  opvisdata$width <- opvisdata$end - opvisdata$start
  op_list[u] <-AnnotationTrack(start = opvisdata$start, 
                               width = opvisdata$width, 
                               strand= opvisdata$strand, 
                               id    = paste(opvisdata$Query_label, opvisdata$ProkkaNO),
                               group = opvisdata$Target_label,
                               mergeGroups = FALSE,
                               feature = opvisdata$Query_label,
                               name  = as.character(opvisdata$ID)[1],
                               background.title="azure4",
                               stacking = "dense",
                               chromosome="chrX")
  
}


ref=GenomeAxisTrack(GRanges('TEST',IRanges(start =seq(0, 6000, by=200),
                                           end = c(seq(200, 6200, by=200)))),
                    lwd=20, fontsize=9)

#### Making heatmap dataframe suitable for plot ####
colnames(heatmapdf) <- c("ID", Queryvec) 
heatmapdf[is.na(heatmapdf)] <- 0

heatmapdf <- as.matrix(heatmapdf)
heatmapdf_plot <- heatmapdf[, 2:(length(Queryvec)+1)]
mode(heatmapdf_plot) <- "numeric"
heatmapdf_plot <- round(heatmapdf_plot, digits = 0)
rownames(heatmapdf_plot) <- heatmapdf[, 1]
colnames(heatmapdf_plot) <- Queryvec

legend_values <- c()
for (r in 1:(length(Queryvec))){
  query_count <- psiblast %>% filter(Query_label == Queryvec[r]) 
  legend_values[r] <- paste(Queryvec[r], nrow(query_count))
}




#### Make a pdf with the heatmap ####
pdf(paste(title, ".pdf"))
##Plot of heatmap data
heatmap.2(heatmapdf_plot, 
          col = c("White", "Yellow", "darkorange","red", "darkred"),
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram = "none",
          density.info = "none",
          trace = "none",
          key = TRUE,
          cexCol = 1,
          cexRow = 0.5,
          margins = c(4, 12),
          lwid = c(2, 7),
          lhei = c(2, 7),
          cellnote = heatmapdf_plot,
          notecex = 0.4,
          notecol = "Darkblue",
          key.xlab = "Identity to Query (%)",
          scale = "none",
          breaks = seq(0, 100, length.out = 6),
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=1:ncol(heatmapdf_plot),
          rowsep=1:nrow(heatmapdf_plot)
)
title(title, cex.main = 1)
legend("topright", legend = legend_values, cex = 0.8)

plotTracks(c(op_list[1:length(op_list)], ref),
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
           col="gray30",
           pgaA="#C8DBC8",
           pgaB="#82ACB9",
           pgaC="#BAB4D8",
           epsK="#B9DBF4",
           rotation.title=1,
           trackLabelHAlign = "center")

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
