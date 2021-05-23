library(dplyr)
library(stringr)
library(tidyverse)
library(readxl)
library(readr)
# Before running:
## Change working dictionary path
## Check the psiblast result path
## Check the query_figure excel file location

# Set working directory and path to psiblast input
str_split(getwd(), pattern = "/") %>% 
  unlist() %>%  
  `[`(1:(length(.)-1)) %>% 
  paste(collapse = "/") %>% 
  paste0(., "/Data") %>% 
  setwd()
psiblast_name <- "HA_Pasteurella"
psiblast_path <- paste("psiblast/", psiblast_name, ".txt", sep="")
query_sheet <- "HA_Pasteurella multocida"

# Read query information
query_metadata <- excel_sheets("Query_figur.xlsx") %>% 
  sapply(function(X) read_xlsx("Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(!names(.) %in% c("Abbreveations","HA_S_pyogenes")) %>% 
  rbindlist()
n_query_gene <- length(unique(query_metadata$Genename))

# Check for previous processing and if this specific synthase has been processed
if ("iTol_percgenes_operons.tsv" %in% list.files("iTol")) {
  iTol_all <- read_tsv("iTol/iTol_percgenes_operons.tsv")
  if (query_sheet %in% names(iTol_all)) {
    continue <- "n"
    while (continue!="y") {
      print(names(iTol_all))
      continue <-  readline(
        prompt = paste("The synthase operon has already been processed.", 
                       "\nContinuing lead to troublesome formating of column names (.x and .y names)",
                       "\ncontinue? y/n", sep=""))
      if (continue=="n") {
        invokeRestart("abort")
      }
    }
  } 
}

# Loading PSI-BLAST data
psiblast <- read.csv(
  file = psiblast_path, header = FALSE,
  sep = "\t", col.names = c(
    "Query_label", "Target_label", "Percent_identity", "align_length", "mismatch",
    "gap_opens", "start_query", "end_query", "start_target",
    "end_target", "E_value", "Bit score"
  )
)

# Sort PSI-BLAST data
psiblast <- psiblast %>%
  separate(Target_label, c("ID", "ProkkaNO"), sep = -5, remove = FALSE, extra = "merge") %>%
  separate(ID, c("ID", "rm"), sep = -1, remove = TRUE) %>%
  filter(Percent_identity > 20)
psiblast <- psiblast[, !(names(psiblast) %in% "rm")]

# Generating iTol formatted data frame
MAG_IDs <- unique(psiblast$ID)
psi_list <- split(psiblast, f = psiblast$ID)
iTol <- data.frame()
for (MAG in MAG_IDs) {
  df <- psi_list[[MAG]]
  iTol[MAG,query_sheet] <- (length(unique(df$Query_label))/n_query_gene)*100
}

# Merge with previous results
if (!exists("iTol_all")) {
  iTol_all <- data.frame()
  iTol_all <- merge(iTol_all, iTol, by=0, all=TRUE)
} else {
  iTol$Row.names <- row.names(iTol)
  iTol_all <- merge(iTol_all, iTol, by="Row.names", all=TRUE)
}
iTol_all[is.na(iTol_all)] <- 0 

write_tsv(iTol_all, "iTol/iTol_percgenes_operons.tsv")
