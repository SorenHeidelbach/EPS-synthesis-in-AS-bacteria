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
library(ggplot2)
library(hrbrthemes)
library(extrafont)
library(forcats)
library(readxl) 
#library(ShortRead)
loadfonts(device = "win")

#Set working directory and path to psiblast input
setwd("C:/Users/heide/OneDrive - Aalborg Universitet/BIOT/8. semester/Project/Data/")

outputMSA=TRUE

upstream = 5
downstream = 5

polysaccharide = "pnag_pga"
ips_path<-paste("ips/", polysaccharide, "_fasta_removed/", sep="") ##Path to interproscan results
ipsq_path<-paste("ipsq/PNAG_fasta_removed/pga_operon_Escherichia_coli_strain_K12.gff3", sep="")

psiblast_surr<-read_tsv(file=paste("psiblast_subset_tsv/", polysaccharide, ".tsv", sep="")) %>%
  separate(MiDAS3_6Tax, c("rm1", "MIDAS"), "Bacteria,") %>%
  separate(MIDAS, c("M_phylum", "M_class", "M_order", "M_family", "M_genus", "M_species"), ",") 
psiblast_surr$Function<-replace_na(psiblast_surr$Function, "N.A.")
psiblast_surr$plot_tax<-paste(psiblast_surr$M_phylum, psiblast_surr$M_order, psiblast_surr$M_genus, sep = ",")
psiblast_surr$Query_label<-replace_na(psiblast_surr$Query_label, " ")

psiblast_surr[replace_na(psiblast_surr$Query_label=="CB", FALSE), "Function"] <- "CB"
psiblast_surr$Query_label <- recode(psiblast_surr$Query_label, CB = ".") 
col_pasra<-c("#d4a0a0", "#d7ceb3", "#b6d4a4", "#94b1b7", "#899da9")

sheet_names = excel_sheets(path = "Query_figur.xlsx")
list_all <- lapply(sheet_names, function(x) read_excel(path = "Query_figur.xlsx", sheet = x, skip = 1))
genedata <- rbind.fill(list_all)
genedata<-genedata[c("Genename", "Function", "Psiblast", "Start", "End", "Strand")]
colnames(genedata)[1]<- "Query_label"
subset_MAGs<-c("Bjer_18-Q3-R1-45_MAXAC.014" , "Ega_18-Q3-R5-49_BATAC.239"  , "EsbE_18-Q3-R3-47_BAT3C.110" ,
               "EsbW_18-Q3-R4-48_MAXAC.044" , "Fred_18-Q3-R57-64_BATAC.644", "Fred_18-Q3-R57-64_MAXAC.027",
               "Fred_18-Q3-R57-64_MAXAC.292", "Hade_18-Q3-R52-61_BATAC.14",  "Hade_18-Q3-R52-61_BATAC.726",
                "Kalu_18-Q3-R12-55_BAT3C.208" ,"Lyne_18-Q3-R50-59_MAXAC.318", "OdNE_18-Q3-R46-58_BAT3C.415",
                "Vibo_18-Q3-R45-57_BATAC.199")

psiblast_surr<-psiblast_surr[psiblast_surr$ID2 %in% subset_MAGs, ]
psiblast_surr<-arrange(psiblast_surr, ID2, ProkkaNO)

psiblast_surr$remove<-c(0, 0, 1, 1, 1, 1, 1, 0,
                        0, 1, 1, 1, 0, 0,
                        0, 0, 1, 1, 1, 1, 0,
                        0, 1, 1, 1, 1, 0, 0,
                        0, 1, 1, 1, 1, 1, 0, 0,
                        0, 0, 1, 1, 1, 1, 0,
                        0, 1, 1, 1, 1, 0, 0,
                        0, 1, 1, 1, 1, 1, 1, 0, 0,
                        0, 1, 1, 1, 1, 0, 0,
                        0, 1, 1, 1, 1, 0, 0,
                        0, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                        0, 0, 1, 1, 1, 1, 0,
                        0, 1, 1, 1, 1, 0, 0)

psiblast_surr<-psiblast_surr %>% filter(psiblast_surr$remove == 1)

#Genedata
genedata1<-filter(genedata, Psiblast == polysaccharide)
Queryvec = genedata1$Query_label
Queryvec <- sort(Queryvec) 


#### Making heatmap dataframe####
psiblast_plot<-psiblast_surr %>% drop_na(Query_label) %>% 
  filter(Query_label != "CB") %>%
  mutate(ID2=factor(ID2,levels=rev(sort(unique(ID2)))))

psiblast_plot<-psiblast_plot[ !duplicated(psiblast_plot[c("Query_label", "ID2")]), ] 

#### Query AnnotationTrack ####


Query_track<-AnnotationTrack(start=genedata1$Start,
                             width=as.numeric(genedata1$End)-as.numeric(genedata1$Start),
                             strand=genedata1$Strand,
                             id=genedata1$Query_label,
                             #group=genedata1$Psiblast,
                             mergeGroups = FALSE,
                             feature = genedata1$Function,
                             name  = "Query",
                             background.title="gray90", 
                             stacking = "dense",
                             chromosome="chrX")
displayPars(Query_track)<-list(showFeatureId=TRUE, 
                               #background.panel="lightgrey",
                               collapse=T, 
                               fontcolor.item='black',
                               fontcolor.group ="grey1",
                               col = "gray75",
                               fontsize.group=45,
                               arrowHeadWidth= 5,
                               arrowHeadMaxWidth=10,
                               cex=2.5,
                               cex.title=3,
                               size = 2.2,
                               lty=1,
                               lwd=1,
                               rotation.title=1,
                               stackHeight=0.8,
                               trackLabelHAlign = "center")

ipsq<-read.csv2(file=ipsq_path, header = FALSE, sep = ";", comment.char = "#")
ipsq<-ipsq[grep("Pfam", ipsq[,1]), ] %>% 
  separate(V2, c("rm2", "rm1"), "Target=") %>%
  separate(rm1, c("Query_label", "start1", "end1"), " ") %>%
  separate(V4, c("rm3", "Domain"), "signature_desc=") %>%
  separate(V5, c("rm4", "PF"), "Name=")

drops = c("V1", "rm3", "rm4", "rm2", "V3", "V6", "V7", "rm5")
ipsq<-ipsq[ , !(names(ipsq) %in% drops)]
ipsq<-left_join(genedata1, ipsq, by = "Query_label")

ipsq2<-ipsq[0,]
for (t in 1:nrow(ipsq)){
  ipsq1<-ipsq[t, ]
  if (ipsq1$Strand == "+"){
    ipsq1$start_d<-as.numeric(ipsq1$Start)+3*as.numeric(ipsq1$start1)
    ipsq1$end_d<-as.numeric(ipsq1$Start)+3*as.numeric(ipsq1$end1)
    ipsq1$width_d<-ipsq1$end_d - ipsq1$start_d
  }else{
    ipsq1$end_d<-as.numeric(ipsq1$End)-3*as.numeric(ipsq1$start1)
    ipsq1$start_d<-as.numeric(ipsq1$End)-3*as.numeric(ipsq1$end1)
    ipsq1$width_d<-ipsq1$end_d - ipsq1$start_d
  }
  ipsq2<-bind_rows(ipsq2, ipsq1)
}
ipsq<-ipsq2
ipsq<-ipsq[!is.na(ipsq$Domain), ]
ipsq$Domain[grep("Glycosyl transferase family 2", ipsq$Domain)] <- "GT2"
ipsq$Domain[grep("Glycosyl transferase Family 4", ipsq$Domain)] <- "GT4"
ipsq$Domain[grep("Glycosyltransferase Family 4", ipsq$Domain)] <- "GT4"
ipsq$Domain[grep("Glycosyl transferase 4-like", ipsq$Domain)] <- "GT4_like"
ipsq$Domain[grep("Glycosyl transferase family 41", ipsq$Domain)] <- "GT41"
ipsq$Domain[grep("Glycosyltransferase like family 2", ipsq$Domain)] <- "GT2_like"
ipsq$Domain[grep("Glycosyl transferases group 1", ipsq$Domain)] <- "GT_g1"
ipsq$Domain[grep("Glycosyl transferase family 2", ipsq$Domain)] <- "GT2"
ipsq$Domain[grep("Glycosyl transferase family 2", ipsq$Domain)] <- "GT 2"
ipsq$Domain[grep("Glycosyl transferase family 2", ipsq$Domain)] <- "GT 2"
ipsq$Domain[grep("Glycosyl transferase family 2", ipsq$Domain)] <- "GT 2"
ipsq$Domain[grep("Glycosyl hydrolases family 8", ipsq$Domain)] <- "GH8"
ipsq$Domain[grep("Glycosyl transferase family 2", ipsq$Domain)] <- "GT 2"
ipsq$Domain[grep("Domain of unknown function", ipsq$Domain)] <- "NA"
ipsq$Domain[grep("deacetylase", ipsq$Domain)] <- "PDA"
ipsq$Domain[grep("Tetratrico", ipsq$Domain)] <- "T"
ipsq$Domain[grep("TPR repeat", ipsq$Domain)] <- "T"
ipsq$Domain[grep("glycosyl hydrolase family 13", ipsq$Domain)] <- "GH13"
ipsq$Domain[grep("glycosyl hydrolase family", ipsq$Domain)] <- "GH"
ipsq$Domain[grep("PilZ domain", ipsq$Domain)] <- "PilZ"
ipsq$Domain[grep("Cellulose synthase", ipsq$Domain)] <- "CS"
ipsq$Domain[grep("Cellulose synthase operon protein C C???terminus (BCSC_C)", ipsq$Domain)] <- "BCSC_C"
ipsq$Domain[grep("Bacterial cellulose synthase subunit", ipsq$Domain)] <- "BCS"
ipsq$Domain[grep("Cellulose biosynthesis protein BcsQ", ipsq$Domain)] <- "BcsQ"
ipsq$Domain[grep("Cellulose biosynthesis protein BcsG", ipsq$Domain)] <- "BcsG"
ipsq$Domain[grep("Cellulose biosynthesis", ipsq$Domain)] <- "CS"
ipsq$Domain[grep("TPR repeat", ipsq$Domain)] <- "T"
ipsq$Domain[grep("Acyltransferase", ipsq$Domain)] <- "ATF"
ipsq$Domain[grep("ABC", ipsq$Domain)] <- "ABC"
ipsq$Domain[grep("Polysaccharide biosynthesis protein", ipsq$Domain)] <- "Pbp"
ipsq$Domain[grep("Polysaccharide biosynthesis C-terminal domain", ipsq$Domain)] <- "Pb_C"
ipsq$Domain[grep("Polysaccharide biosynthesis/export protein", ipsq$Domain)] <- "Pbep"
ipsq$Domain[grep("PgaD", ipsq$Domain)] <- "pgaD"


Query_track_domain<-AnnotationTrack(start = ipsq$start_d, 
                                    width = ipsq$width_d, 
                                    #strand= ips$strand, 
                                    id    = ipsq$Domain,
                                    #group = ips$Target_label,
                                    mergeGroups = FALSE,
                                    feature = ipsq$Domain,
                                    name  = " ",
                                    background.title="gray90", 
                                    col.title="gray90",
                                    stacking = "dense",
                                    chromosome="chrX")


displayPars(Query_track_domain)<-list(
  showFeatureId=TRUE, 
  collapse=T, 
  fontcolor.item='black',
  fontcolor.group ="grey1",
  fontsize.group=25,
  cex.main = 1,
  main = "Operon Structure",
  cex=1.4,
  cex.title=0.55,
  lty=5,
  lwd=0.1,
  col="transparent",
  trackLabelHAlign = "center",
  size=0.9,
  rotation.item=0
)


#### Generate a list of annotationtracks for plotting ####
op_list <- list()
loop_genes<-seq(1, length(unique(psiblast_surr$ID2))*3, by = 3)
loop_domain<-seq(2, length(unique(psiblast_surr$ID2))*3+1, by = 3)
loop_empty<-seq(3, length(unique(psiblast_surr$ID2))*3+2, by = 3)
IDs2<-unique(psiblast_surr$ID2)
psiblast_surr <- psiblast_surr %>%
  mutate(ID = factor(ID, levels = IDs2))

for (u in 1:length(IDs2)){
  opvisdata <- psiblast_surr %>% filter(ID2 == IDs2[u])
  minstart <- min(opvisdata$start)
  opvisdata$start <- opvisdata$start - minstart
  opvisdata$end <- opvisdata$end - minstart
  opvisdata$width <- opvisdata$end - opvisdata$start
  op_list[loop_genes[u]]<-AnnotationTrack(start = opvisdata$start, 
                                          width = opvisdata$width, 
                                          strand= opvisdata$strand, 
                                          id    = opvisdata$Query_label,
                                          group = opvisdata$Target_label,
                                          mergeGroups = FALSE,
                                          feature = opvisdata$Function,
                                          name  = as.character(paste(opvisdata$ID2[1], rev(sort(opvisdata$plot_tax))[1], sep="\n")),
                                          background.title="azure4",
                                          stacking = "dense",
                                          chromosome="chrX")
  
  
  ips<-read.csv2(file=paste(ips_path, polysaccharide, "_", IDs2[u], ".gff3", sep=""), header = FALSE, sep = ";", comment.char = "#")
  ips<-ips[grep("Pfam", ips[,1]), ] %>% 
    separate(V2, c("rm2", "rm1"), "Target=") %>%
    separate(rm1, c("Target_label", "start1", "end1"), " ") %>%
    separate(V4, c("rm3", "Domain"), "signature_desc=") %>%
    separate(V5, c("rm4", "PF"), "Name=") %>%
    mutate(Target_label1=Target_label) %>%
    separate(Target_label1, c("ID", "ProkkaNO"), sep=-5) %>%
    separate(ID, c("ID", "rm5"), sep=-1)
  
  drops = c("V1", "rm3", "rm4", "rm2", "V3", "V6", "V7", "rm5")
  ips<-ips[ , !(names(ips) %in% drops)]
  ips<-ips[ips$ProkkaNO %in% opvisdata$ProkkaNO, ]
  ips<-left_join(opvisdata, ips, by = "Target_label")
  ips2<-ips[0,]
  for (t in 1:nrow(ips)){
  ips1<-ips[t, ]
    if (ips1$strand == "+"){
  ips1$start_d<-as.numeric(ips1$start)+3*as.numeric(ips1$start1)
  ips1$end_d<-as.numeric(ips1$start)+3*as.numeric(ips1$end1)
  ips1$width_d<-ips1$end_d - ips1$start_d
  }else{
  ips1$end_d<-as.numeric(ips1$end)-3*as.numeric(ips1$start1)
  ips1$start_d<-as.numeric(ips1$end)-3*as.numeric(ips1$end1)
  ips1$width_d<-ips1$end_d - ips1$start_d
  }
  ips2<-bind_rows(ips2, ips1)
  }
  ips<-ips2
  ips<-ips[!is.na(ips$Domain), ]
  ips$Domain[grep("Glycosyl transferase family 2", ips$Domain)] <- "GT2"
  ips$Domain[grep("Glycosyl transferase Family 4", ips$Domain)] <- "GT4"
  ips$Domain[grep("Glycosyl transferase 4-like", ips$Domain)] <- "GT4_like"
  ips$Domain[grep("Glycosyl transferase family 41", ips$Domain)] <- "GT41"
  ips$Domain[grep("Glycosyltransferase like family 2", ips$Domain)] <- "GT2_like"
  ips$Domain[grep("Glycosyltransferase Family 4", ips$Domain)] <- "GT4"
  ips$Domain[grep("Glycosyl transferase Family 4", ips$Domain)] <- "GT4"
  ips$Domain[grep("Glycosyl transferase family 2", ips$Domain)] <- "GT2"
  ips$Domain[grep("Glycosyl transferase family group 2", ips$Domain)] <- "GT2"
  ips$Domain[grep("Glycosyl transferase family 2", ips$Domain)] <- "GT 2"
  ips$Domain[grep("Glycosyl hydrolases family 8", ips$Domain)] <- "GH8"
  ips$Domain[grep("Glycosyl transferase family 2", ips$Domain)] <- "GT 2"
  ips$Domain[grep("Domain of unknown function", ips$Domain)] <- "NA"
  ips$Domain[grep("deacetylase", ips$Domain)] <- "PDA"
  ips$Domain[grep("Tetratrico", ips$Domain)] <- "T"
  ips$Domain[grep("glycosyl hydrolase family 13", ips$Domain)] <- "GH13"
  ips$Domain[grep("glycosyl hydrolase", ips$Domain)] <- "GH"
  ips$Domain[grep("PilZ domain", ips$Domain)] <- "PilZ"
  ips$Domain[grep("Cellulose synthase", ips$Domain)] <- "CS"
  ips$Domain[grep("Cellulose synthase operon protein C C???terminus (BCSC_C)", ips$Domain)] <- "BCSC_C"
  ips$Domain[grep("Bacterial cellulose synthase subunit", ips$Domain)] <- "BCS"
  ips$Domain[grep("Cellulose biosynthesis protein BcsQ", ips$Domain)] <- "BcsQ"
  ips$Domain[grep("Cellulose biosynthesis protein BcsG", ips$Domain)] <- "BcsG"
  ips$Domain[grep("Cellulose biosynthesis", ips$Domain)] <- "CS"
  ips$Domain[grep("Tetratrico", ips$Domain)] <- "T"
  ips$Domain[grep("TPR repeat", ips$Domain)] <- "T"
  ips$Domain[grep("Acyltransferase", ips$Domain)] <- "ATF"
  ips$Domain[grep("Glycosyl transferases group 1", ips$Domain)] <- "GTg1"
  ips$Domain[grep("Polysaccharide biosynthesis protein", ips$Domain)] <- "Pbp"
  ips$Domain[grep("Polysaccharide biosynthesis/export protein", ips$Domain)] <- "Pbep"
  ips$Domain[grep("Polysaccharide biosynthesis C-terminal domain", ips$Domain)] <- "Pb_C"
  ips$Domain[grep("PgaD", ips$Domain)] <- "pgaD"
  
  
  
  op_list[loop_domain[u]]<-AnnotationTrack(start = ips$start_d, 
                                           width = ips$width_d, 
                                           #strand= ips$strand, 
                                           id    = ips$Domain,
                                           #group = ips$Target_label,
                                           mergeGroups = FALSE,
                                           feature = ips$Domain,
                                           name  = " ",
                                           background.title="transparent", 
                                           col.title ="transparent",
                                           stacking = "dense",
                                           chromosome="chrX")
  
  op_list[loop_empty[u]] <- GenomeAxisTrack(
    range=IRanges(
      start = c(0),
      end = c(7000)
    ),
    scale=1, 
    distFromAxis=0,
    col="gray30"
  )
  
}

ref <- GenomeAxisTrack(
  range=IRanges(
    start = c(0),
    end = c(7000)
  ),
  distFromAxis=0,
  col="gray5",
  cex=1.8
)


legend_genes<-c("GT - Glycosyltransferase", "PDA - Polysaccharide deacetylase", "TPR - Tetratricopeptide repeat", "GH - Glycosyl hydrolase")

#### Defining display parameters for gene and domain annotationtracks####
for (i in 1:length(IDs2)){
  displayPars(op_list[[loop_genes[i]]])<-list(showFeatureId=TRUE, 
                                              #background.panel="lightgrey",
                                              collapse=T, 
                                              fontcolor.item='black',
                                              fontcolor.group ="grey1",
                                              col = "gray75",
                                              fontsize.group=45,
                                              arrowHeadWidth= 5,
                                              arrowHeadMaxWidth=10,
                                              cex=2.3,
                                              cex.title=1.6,
                                              size = 1.8,
                                              lty=1,
                                              lwd=1,
                                              rotation.title=1,
                                              stackHeight=0.8,
                                              trackLabelHAlign = "center")
  
  displayPars(op_list[[loop_domain[i]]])<-list(
    showFeatureId=TRUE, 
    collapse=T, 
    fontcolor.item='black',
    fontcolor.group ="grey1",
    fontsize.group=25,
    cex.main = 1,
    main = "Operon Structure",
    cex=1.4,
    cex.title=0.55,
    lty=5,
    lwd=0.1,
    col="transparent",
    trackLabelHAlign = "center",
    size=0.9,
    rotation.item=0
    
  )
  displayPars(op_list[[loop_empty[i]]])<-list(
    showFeatureId=TRUE, 
    collapse=T, 
    fontsize.group=25,
    size=1
  )
}

#### Making a dataframe to the pdf ####
psiblast_alt_subset_bind <- psiblast_surr[order(psiblast_surr$Query_label, -abs(psiblast_surr$start) ), ] 
psiblast_alt_subset_bind <- psiblast_alt_subset_bind[ !duplicated(psiblast_alt_subset_bind$ID2), ]
df_pdf <- psiblast_alt_subset_bind %>% select(ID2, plot_tax, seqname, start, Query_label, ProkkaNO) %>% 
  arrange(ID2)


n_pages_operon<-ceiling(length(op_list)/12)
n_pages_heatmap<-ceiling(length(unique(psiblast_surr$ID2))/15)
n_pages_table<-ceiling(nrow(df_pdf)/25)

#### Make a pdf with the heatmap and operon structures ####

pdf(paste("r_output/", polysaccharide, "_operon.pdf", sep=""), width = 27, height = 27)

  plotTracks(c(Query_track, Query_track_domain, op_list[3], op_list, ref),
             showFeatureId=TRUE,
             background.title="transparent",
             col.title="black",
             #from = -100, 
             #to=14400,
             #stacking = "dense",
             collapse=T, 
             fontcolor.item='black',
             fontcolor.group ="grey1",
             title.width = 14,
             cex.main = 5,
             #main = "Manually Verified Cellulose Structures",
             arrowHeadWidth= 3,
             arrowHeadMaxWidth=10,
             #cex=0.3,
             #cex.title=0.9,
             rotation.title=0,
             trackLabelHAlign = "left",
             extend.right = 0,
             extend.left = 0,
             innermargin = 0,
             #col="gray70",
             
             
             pgaD = "#B3CFD0",
             PDA= '#f5f5c6',
             TPR="burlywood",
             GH8='#fcdbbb',
             GH13='#fcdbbb',
             CS="#B3CFD0",
             PilZ='#f5f5c6',
             ATF=col_pasra[4],
             BCS='#B3CFD0',
             BcsQ='#deb887',
             BcsG='#f5f5c6',
             CBS=col_pasra[3],
             
             SY="#B3CFD0",
             MOD="#deb887",
             Export="#EDE9EB",
             pgaA="#C8DBC8",
             pgaB="#82ACB9",
             pgaC="BAB4D8",
             epsK="#B9DBF4",
             GT="#82ACB9",
             GT2="#82ACB9",
             GT4="#82ACB9",
             GT41="#82ACB9",
             GT2_like="#82ACB9",
             GH="brown",
             N.A.="white",
             CB="red",
             ccpA="#C8DBC8",
             # 
             ABC = '#dccede',
             PE = '#C8DBC8',
             Export = '#EDE9EB',
             Branch = '#ffcb85',
             MOD = '#deb887',
             Hydrolase = '#fcdbbb', 
             Lyase = '#FBDADA',  
             GT = '#82ACB9',
             cmc = '#FBDADA',
             Polymerization = '#a1c2a1',
             Prim_GT = '#ECDED5',
             SY = '#B3CFD0',
             PS = '#D68A8A',  #Purple
             REG ='#f5f5c6',
             `T` = "gray90",
             .="red4",
             
             
             
             rotation.title=0,
             trackLabelHAlign = "center")

dev.off(which = dev.cur())
#### ####