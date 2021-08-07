library("data.table")
library("readxl")
library("glue")
library("gggenes")
library("ggtext")
library("broom")
library("here")
library("seqinr")
library("tidyverse")
needs::prioritize("dplyr")
setwd(here::here())

##----------------------------------------------------------------
##  Loading MAG statistic, Prokka annotations and query metadata  
##----------------------------------------------------------------
# Statistics on the HQ-MAGs from singleton et.al 2021
magstats. <- fread("./data/raw/magstats.tsv", sep = "\t")

# File with metadata on the query genes, e.g. function
query_metadata. <- excel_sheets("./data/raw/Query_figur.xlsx") %>%
  sapply(function(X) read_xlsx("./data/raw/Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))) %>%
  rbindlist() 

# Additional information of the prokka annotations
gff. <- fread("./data/raw/gff.tsv") %>%
  filter(!(ProkkaNO == "units")) %>% 
  mutate(Target_label = paste(ID, ProkkaNO, sep = "_"),
         ProkkaNO = as.numeric(ProkkaNO))



proximity_filtration <- function(filename_psiblast,
                                 perc_id = 20,
                                 max_dist_prok = 8,
                                 max_dist_gene = 2000,
                                 min_genes = 2,
                                 flanking_genes = 2,
                                 magstats = magstats.,
                                 query_metadata = query_metadata.,
                                 gff = gff.,
                                 exclude_gene = "",
                                 essential_genes = NA
                                 ){
  filename_psiblast_col <- paste(filename_psiblast, collapse = "_")
  unlink(glue("./output/psi_operon_full/{filename_psiblast_col}.tsv"))
  unlink(glue("./output/psi_proxi_filt/{filename_psiblast_col}.tsv"))
  unlink(glue("./output/psi_percID_filt/{filename_psiblast_col}.tsv"))
  ##---------------------------------------------------------------
  ##                Processing of PSI-BLAST results                
  ##---------------------------------------------------------------
  # Loading raw psiblast results
  test <- filename_psiblast %>%
    paste0("./data/raw/psiblast/", ., ".txt") %>%
    lapply(
      fread,
      header = FALSE, sep = "\t", fill = TRUE, 
      col.names = c(
        "Query_label", "Target_label", "Percent_identity", "align_length",
        "mismatch", "gap_opens", "start_query", "end_query",
        "start_target", "end_target", "E_value", "Bit score")
    ) %>%
    bind_rows() %>%
    filter(!(Query_label %in% exclude_gene)) %>% 
    # ASSIGNING psi_raw
    assign(x = "psi_raw", value = ., pos = 1) %>% 
  ##---------------------------------------------------------------
    # Filtering
    subset(!is.na(E_value)) %>%
    filter(Percent_identity > perc_id) %>%
    group_by(Target_label) %>%
    filter(`Bit score` == max(`Bit score`)) %>%
    filter(!duplicated(Target_label)) %>%
    # Cleaning
    separate(Target_label, c("ID", "ProkkaNO"),
             sep = -5, remove = FALSE, extra = "merge") %>%
    mutate(ProkkaNO = as.numeric(ProkkaNO)) %>%
    mutate(ID = str_sub(ID, start = 1, end = -2)) %>% 
    # Merging
    inner_join(gff) %>%
    left_join(magstats, by = "ID", keep = FALSE) %>%
    left_join(query_metadata[query_metadata$Psiblast %in% filename_psiblast, c("Genename", "Function")], 
              by = c("Query_label" = "Genename"), keep=FALSE)  %>%
    arrange(Target_label, ProkkaNO) %>%
    # Operon Grouping
    group_by(seqname, ID) %>%
    # define gene and prokka distance to posterior and prior psiblast hit
    mutate(prio_prok = ProkkaNO - shift(ProkkaNO) < max_dist_prok,
           prio_gene = start - shift(end, 1) < max_dist_gene,
           prio_prok = replace_na(prio_prok, FALSE),
           prio_gene = replace_na(prio_gene, FALSE)
           ) %>%
    ungroup() %>% 
    # If a gene only satisfy distance to prior hit and not posterior = start
    mutate(operon_place = ifelse(!(prio_prok | prio_gene), "start", NA),
           operon = ifelse(operon_place == "start", row_number(), NA)) %>%
    fill(operon, .direction = "down") %>%
    # ASSIGNING psi_filtered
    assign(x = "psi_perc_ID_filt", value = ., pos = 1) %>% 
  ##---------------------------------------------------------------
    # Operon number of genes filtering
    group_by(operon) %>%
    filter(length(unique(Query_label)) >= min_genes & (all(essential_genes %in% Query_label) | is.na(essential_genes))) %>%
    {if(nrow(.) == 0) stop("No results, try less strict filtration. Filter by less unique genes, remove essential genes, increase percent identity...") else .} %>% 
    select(-prio_prok, -prio_gene) %>% 
    ungroup() %>% 
    nest(cols = !c(ID, operon)) %>% 
    group_by(ID) %>% 
    mutate(letter = toupper(letters)[seq_along(unique(operon))],
           ID2 = case_when(
             length(unique(operon)) > 1 ~ paste(ID, letter, sep = "_"),
             TRUE ~ ID)
           ) %>% 
    select(-letter) %>% unnest(cols = c(cols)) %>% 
    # ->ASSIGNING "psiblast"<- 
    assign(x = "psi_proxi_filt", value = ., pos = 1) %>% 
  ##---------------------------------------------------------------
    ## Expansion of filtered psiblast hits with surrounding genes
    # Identify surrounding genes
    group_by(operon) %>%
    fill(ID2, .direction = "updown") %>% 
    summarise(
      ID = unique(ID),
      tig = unique(seqname),
      max = max(ProkkaNO),
      min = min(ProkkaNO)
      ) %>%
    full_join(gff) %>%
    group_by(operon) %>%
    filter(ProkkaNO <= (max + flanking_genes) &
           ProkkaNO >= (min - flanking_genes) & 
           seqname  == tig) %>%
    # Merging surrounding and psiblast genes
    full_join(psi_proxi_filt) %>%
  ##---------------------------------------------------------------
    # Prepare for plotting in gggenes - see plot_operon.R
    group_by(operon) %>%
    mutate(
      # gggenes require 1 and -1 as strand direction
      strand = case_when(strand == "+" ~ 1, strand == "-" ~ -1, TRUE ~ 0),
      # Relative distances (start always minimum value of gene posistion)
      end = end - min(start),
      start = start - min(start),
      # converting AA numbers to nucleotide numbers
      end_target_plot =   ifelse(is.na(end_target),   end,   end_target   * 3 + start),
      start_target_plot = ifelse(is.na(start_target), start, start_target * 3 + start),
      # Assigning artificial percent identity to domains
      Percent_identity = ifelse(is.na(Percent_identity), 40, Percent_identity),
      # Order operons from highest bit score to lowest
      operon = factor(operon, levels = unique(arrange(., -`Bit score`)$operon))
      ) %>%
    fill(MiDAS3_6Tax, .direction = "updown") %>%
    fill(SILVA138Tax, .direction = "updown") %>%
    fill(GTDBTax, .direction = "updown") %>%
    fill(ID2, .direction = "updown") %>% 
    # ->ASSIGNING "genes"<- 
    assign(x = "psi_operon_full", value = ., pos = 1)
  
  ##---------------------------------------------------------------
  ##                   Writing data to tsv files                   
  ##---------------------------------------------------------------
  dir.create("./output/psi_percID_filt", showWarnings = FALSE)
  write.table(psi_perc_ID_filt, file = glue("./output/psi_percID_filt/{filename_psiblast_col}.tsv"),#_qualfilt_percID{perc_id}.tsv"), 
              quote = F, sep = "\t", row.names = F)
  
  dir.create("./output/psi_proxi_filt", showWarnings = FALSE)
  #write.table(psi_proxi_filt, file = glue("./output/psi_proxi_filt/{filename_psiblast_col}_genes{min_genes}_dist{max_dist_prok}&{max_dist_gene}_percID{perc_id}.tsv"), 
  #            quote = F, sep = "\t", row.names = F)
  write.table(psi_proxi_filt, file = glue("./output/psi_proxi_filt/{filename_psiblast_col}.tsv"), 
              quote = F, sep = "\t", row.names = F)
  
  dir.create("./output/psi_operon_full", showWarnings = FALSE)
  #write.table(psi_operon_full, file = glue("./output/psi_operon_full/{filename_psiblast_col}_genes{min_genes}_dist{max_dist_prok}&{max_dist_gene}_percID{perc_id}_flank{flanking_genes}.tsv"), 
  #            quote = F, sep = "\t", row.names = F)
  write.table(psi_operon_full, file = glue("./output/psi_operon_full/{filename_psiblast_col}.tsv"), 
              quote = F, sep = "\t", row.names = F)
  
  ##------------------------------------------------------------------
  ##  Export fasta files of identified genes, with surrounding genes  
  ##------------------------------------------------------------------
  for (f in unique(psi_operon_full$ID)){
    dir.create(glue("./data/processed/fasta_output/{filename_psiblast_col}"), showWarnings = FALSE, recursive = TRUE)
    read.fasta(file = glue("./data/raw/MGP1000_HQMAG1083_prot_db_split/{f}.faa"), 
                            seqtype="AA", 
                            as.string=TRUE, 
                            set.attributes=FALSE) %>% 
      `[`(names(.) %in% psi_operon_full$Target_label) %>% 
      write.fasta(names=names(.), file.out=glue("./data/processed/fasta_output/{filename_psiblast_col}/{f}.faa"))
  }
}









# Exploratory code to check the effect of different operon parameters
# PC analasys used to check

# library("ggbiplot")
# Manually verified operons ()
# verified <- glue("./data/processed/manually_verified/{filename_psiblast}.txt") %>% 
#   file.info() %>% 
#   filter(size > 0) %>% 
#   rownames() %>% 
#   lapply(fread, header = FALSE) %>% 
#   rbindlist() %>% 
#   unique() %>% 
#   `[[`(1)
# 
# library(ggbiplot)
# psi_compare %>% 
#   filter(!is.na(E_value)) %>% 
#   ungroup() %>% 
#   mutate(coverage = (end_target - start_target) / (end - start)/3,) %>% 
#   select(`Bit score`, Percent_identity, coverage) %>% 
#   prcomp(center = TRUE, scale. = TRUE) %>% 
#   ggbiplot(groups = psi_compare$in_operon, ellipse = TRUE)
# 
# psi_PCA_mean <-  psiblast %>% 
#   ungroup() %>% 
#   mutate(coverage = (end_target - start_target) / (end - start)/3,) %>%
#   group_by(operon) %>% 
#   dplyr::summarize(
#     mean_bit = mean(`Bit score`),
#     mean_cov = mean(coverage),
#     mean_perc_id = mean(Percent_identity),
#     size = length(coverage)
#     ) %>%
#   inner_join(psiblast, by = "operon") %>% 
#   ungroup() %>% 
#   full_join(psi_filtered) %>%
#   mutate(coverage = (end_target - start_target) / (end - start)/3,) %>%
#   mutate(
#     in_operon = case_when(
#       Target_label %in% subset(psiblast, ID %in% verified)$Target_label ~ "Manual",
#       Target_label %in% psiblast$Target_label ~ "Automatic",
#       TRUE ~ "Removed"
#     ),
#     mean_bit = ifelse(is.na(mean_bit), yes = `Bit score`, no = mean_bit), 
#     mean_cov = ifelse(is.na(mean_cov), yes = coverage, no = mean_cov),
#     mean_perc_id = ifelse(is.na(mean_perc_id), yes = Percent_identity, no = mean_perc_id),
#     size = ifelse(is.na(size), 1, size),
#     size2 = as.character(size)
#     ) %>%
#   ungroup()
# psi_PCA_mean2 <- psi_PCA_mean %>%  
#  select(mean_bit, mean_cov, mean_perc_id, size) 
# 
# 
# prcomp(psi_PCA_mean2, scale. = TRUE, center = TRUE) %>% 
#   ggbiplot(groups = psi_PCA_mean$size2)
# 
