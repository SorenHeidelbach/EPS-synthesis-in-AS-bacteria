library("data.table")
library("tidyverse")
library("tidyr")
library("treeio")
library("ggtree")
library("readxl")
setwd(here::here())


##----------------------------------------------------------------
##                    Import phylum of HQ-MAGs                    
##----------------------------------------------------------------
phylum <- read_tsv("./data/raw/magstats.tsv") %>% 
  mutate(
    phylum =   str_extract(MiDAS3_6Tax, "p:[^,]*"),
    phylum =   str_replace(phylum, "p:",replacement = "")
  ) %>% 
  select(ID, phylum) %>% 
  setNames(c("label", "phylum"))


##------------------------------------------------------------------------
##  Import percentage of identified genes in each HQ-MAG for all queries  
##------------------------------------------------------------------------
# File with metadata on the query genes, e.g. function
query_metadata. <- excel_sheets("./data/raw/Query_figur.xlsx") %>%
  sapply(function(X) read_xlsx("./data/raw/Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))) %>%
  rbindlist()

psi_proxi_filt <- list.files("./output/psi_proxi_filt/") %>% 
  map(.f = function(x){
    name_query <- tools::file_path_sans_ext(x)
    fread(paste0("./output/psi_proxi_filt/", x)) %>% 
      group_by(ID) %>% 
      summarise(genes_percent = length(unique(Query_label))) %>% 
      mutate(
        genes_percent = genes_percent/(query_metadata. %>% filter(Psiblast == "alginate") %>% nrow)
      ) %>% 
      setNames(c("ID", name_query))
    }
  ) %>% 
  reduce(full_join, by = "ID") %>% 
  tibble::column_to_rownames(var = "ID")

##---------------------------------------------------------------
##            Import meta information for each HQ-MAG            
##---------------------------------------------------------------
# Import phylum of HQ-MAGs                    
phylum <- read_tsv("./data/raw/magstats.tsv") %>% 
  mutate(
    phylum =   str_extract(MiDAS3_6Tax, "p:[^,]*"),
    phylum =   str_remove(phylum, "p:")
  ) %>% 
  select(ID, phylum) %>% 
  setNames(c("label", "phylum"))

# combine information to each node
tree_data <- read.tree("./data/processed/iTol/MGP1000_bac_1080_maaike_ITOL.tree") %>% 
  as_tibble() %>% 
  full_join(phylum, by = "label") %>% 
  select(node, phylum) %>% 
  filter(!is.na(phylum) & !is.na(node)) %>% 
  as_tibble() 


tree_tib <- read.tree("./data/processed/iTol/MGP1000_bac_1080_maaike_ITOL.tree") %>% 
  tidytree::treedata(phylo = ., data = tree_data) %>% 
  as_tibble() %>% 
  filter(!is.na(node)) %>% 
# Add phylum to intermediate nodes (only thoose with unqiue offspring phylum)
  mutate(
    offspring_phylum = map(node, function(x) offspring(., x) %>% pull(phylum) %>% unique %>% `[`(!is.na(.))),
    phylum2 = case_when(
      offspring_phylum %>% map(length) <= 1 ~ as.character(offspring_phylum),
      TRUE ~ "NA"
    ),
    phylum_nodes = case_when(
      phylum2 == "character(0)" ~ phylum,
      TRUE ~ phylum2
    ),
    phylum_nodes = ifelse(phylum_nodes == "NA", NA, phylum_nodes)
  ) %>% 
  select(-phylum2, -offspring_phylum) 

tree <- tree_tib %>% 
  as.treedata()

#ggsave("./figures/test.pdf", width = 10, height = 20,
  tree_plot <- ggtree(tree, aes(color = phylum_nodes), layout = "circular", lwd = 0.15)
#)
ggsave("./figures/test.pdf", width = 10, height = 20, 
  gheatmap(tree_plot, data = psi_proxi_filt, 
           width = 1,
           colnames_angle=90, font.size = 2,
           hjust = 1) +
    scale_color_continuous()
)




