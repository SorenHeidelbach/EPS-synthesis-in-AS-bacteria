library("tidyverse")
library("readxl")
library("here")
library("data.table")
# Set working directory and path to psiblast input
setwd(here())

# Read query information
query_metadata <- excel_sheets("./data/raw/Query_figur.xlsx") %>% 
  sapply(function(X) read_xlsx("./data/raw/Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(names(.)[!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))]) %>% 
  rbindlist(fill=T)

# Get all MAG ids
MAGs <- read_tsv("./data/raw/magstats.tsv") %>% 
  as.data.frame() %>% 
  `[`(,1)
iTol_all <- data.frame(MAG = MAGs)

for (i in list.files("./output/psi_operon_filt/")) {
  psiblast_file_name <- tools::file_path_sans_ext(i)
  psiblast_subset_path <- paste("./output/psi_operon_filt/", i, sep="")
  psiblast_file_name = "pnag_ica"
  query <- subset(query_metadata, str_detect(query_metadata$Psiblast, psiblast_file_name))
  n_query_gene <- length(unique(query$Genename))
  
  # Loading PSI-BLAST data (Removing NA rows)
  psiblast <- read_tsv(file = psiblast_subset_path) %>% 
    filter(!is.na(Query_label)) %>% 
    subset(Query_label != "CB")
  
  # Generating iTol formatted data frame
  iTol <- data.frame(MAG = MAGs)
  MAG_IDs <- unique(psiblast$ID)
  for (MAG in MAG_IDs) {
    df <- subset(psiblast, ID == MAG)
    iTol[iTol$MAG == MAG, psiblast_file_name] <- (length(unique(df$Query_label))/n_query_gene)*100
    if (length(unique(df$Query_label)) > n_query_gene) {
      show(unique(df$Query_label))
    }
  }
  iTol[is.na(iTol)] <- 0
  
  # Merge with previous results
  iTol_all <- merge(iTol_all, iTol, by="MAG", all=TRUE)
}

write_tsv(iTol_all, "./output/iTol_operon/iTol_percgenes_found_operons.tsv")

