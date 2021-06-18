library(dplyr)
library(tidyverse)
library(readxl)
library(readr)
library(here)
library("data.table")
library(here)
d <- function(x, raw = TRUE){paste0(here(),"/data/",ifelse(raw, "raw/", "processed/"), x)}
o <- function(x){paste(here(), x, sep = "/output/")}
f <- function(x){paste(here(), x, sep = "/figures/")}
s <- function(x){paste(here(), x, sep = "/scripts/")}

# Before running:
## Change working dictionary path
## Check the psiblast result path
## Check the query_figure excel file location

# Set working directory and path to psiblast input
setwd(d(""))

# Read query information
query_metadata <- excel_sheets("Query_figur.xlsx") %>% 
  sapply(function(X) read_xlsx("Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(names(.)[!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))]) %>% 
  rbindlist(fill=T)

# Get all MAG ids
MAGs <- read_tsv("magstats.tsv") %>% 
  as.data.frame() %>% 
  `[`(,1)
iTol_all <- data.frame(MAG = MAGs)

for (i in list.files(d("psiblast_subset_tsv/", raw=FALSE))) {
  psiblast_file_name <- tools::file_path_sans_ext(i)
  psiblast_subset_path <- d(paste("psiblast_subset_tsv/", psiblast_file_name, ".tsv", sep=""), raw=FALSE)
  
  Query <- subset(query_metadata, Psiblast==psiblast_file_name)
  query_sheet <- unique(subset(query_metadata, Psiblast==psiblast_file_name)$Polysaccheride)
  n_query_gene <- length(unique(Query$Genename))
  
  # Loading PSI-BLAST data (Removing NA rows)
  psiblast <- read_tsv(file = psiblast_subset_path) 
  psiblast <- psiblast[!is.na(psiblast$Query_label),]
  psiblast <- subset(psiblast, Query_label != "CB")
  
  # Generating iTol formatted data frame
  iTol <- data.frame(MAG = MAGs)
  MAG_IDs <- unique(psiblast$ID)
  for (MAG in MAG_IDs) {
    df <- subset(psiblast, ID == MAG)
    iTol[iTol$MAG == MAG, query_sheet] <- (length(unique(df$Query_label))/n_query_gene)*100
    if (length(unique(df$Query_label)) > n_query_gene) {
      show(unique(df$Query_label))
    }
  }
  iTol[is.na(iTol)] <- 0
  
  # Merge with previous results
  iTol_all <- merge(iTol_all, iTol, by="MAG", all=TRUE)
}

write_tsv(iTol_all, paste0(here(), "/data/processed/iTol_percgenes_found_operons.tsv"))
show(unique(df$Query_label))
