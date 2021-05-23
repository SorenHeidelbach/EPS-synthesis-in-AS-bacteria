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
setwd("C:/Users/morte/OneDrive - Aalborg Universitet/8. semester/P8/Data")

#### Defining the gffRead function (Maaike's) and loading gff data ####
gffRead <- function(gffFile){
  gff <- read.delim(gffFile, header=F, comment.char="#")
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  gff <- gff[complete.cases(gff), ]
  
  gff <- gff %>%
    separate(attributes, c("prokkaID", "attributes"), sep = ";", extra = "merge") %>%
    separate(prokkaID, c("ID1", "ProkkaNO"), sep = -5)
  
  
  return(gff)}

#Using the gffRead function to import gff files
gff <- tibble(
  file = list.files("gff_files_reduced_out"),
  location = list.files("gff_files_reduced_out", full.names = TRUE),
  ID = str_remove(file, ".gff"),
  df = map(location, gffRead)) %>% 
  unnest(df)

#Removing columns redundant to psiblast data in gff dataframe
drops = c("score", "feature", "file", "location", "frame", "attributes", "MAG", "source", "ID1")
gff <- gff[ , !(names(gff) %in% drops)]



#### Loading MAG statistics data ####
magstats = read.csv(file = "MAG_statistics_STABLEX_20200213.tsv", header = TRUE, sep = "\t")
magstats <- magstats %>%
  separate(MAG, c("ID", "rm"), sep =-3, remove = TRUE)
magstats = magstats[ , !(names(magstats) %in% "rm")]
drops = c("MaxContigBP", "AvContigBP", "HQMAG", "HQdRep", "HQdRep99ANI", "TotBP", "NumContigs", 
          "HQSpRep", "Comp", "Cont", "Circ", "StrHet", "FLSSU", "FLLSU", "ilmcov", "npcov",
          "mm27f", "mm534r", "ttl_polymorphic_rate", "polymut_rate")
magstats <- magstats[ , !(names(magstats) %in% drops)]

write_tsv(magstats, file="magstats.tsv")
write_tsv(gff, file="gff.tsv")