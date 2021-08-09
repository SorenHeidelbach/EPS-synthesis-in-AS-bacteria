library("data.table")
library("tidyverse")
library("tidyr")
library("treeio")
library("ggtree")
library("readxl")
library("dplyr")
setwd(here::here())


##------------------------------------------------------------------------
##  Import percentage of identified genes in each HQ-MAG for all queries  
##------------------------------------------------------------------------
# File with metadata on the query genes, mostly for total genes in query
query_metadata. <- excel_sheets("./data/raw/Query_figur.xlsx") %>%
  sapply(function(X) read_xlsx("./data/raw/Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))) %>%
  rbindlist()

# Percentage of percent identity filtrated genes in each HQ-MAG
psi_perc_filt <- list.files("./output/psi_percID_filt/") %>% 
  map(function(query){
    name_query <- tools::file_path_sans_ext(query)
    fread(paste0("./output/psi_percID_filt/", query)) %>% 
      group_by(ID) %>% 
      summarise(genes_percent = length(unique(Query_label))) %>% 
      mutate(
        genes_percent = genes_percent/(query_metadata. %>% filter(Psiblast == name_query) %>% nrow)
      ) %>% 
      setNames(c("ID", name_query))
    }) %>% 
  reduce(full_join, by = "ID") %>% 
  tibble::column_to_rownames(var = "ID") %>% 
  replace(., is.na(.), 0)

# Percentage of proximity filtrated genes in each HQ-MAG
psi_proxi_filt <- list.files("./output/psi_proxi_filt/") %>% 
  map(.f = function(query){
    name_query <- tools::file_path_sans_ext(query)
    fread(paste0("./output/psi_proxi_filt/", query)) %>% 
      group_by(ID) %>% 
      summarise(genes_percent = length(unique(Query_label))) %>% 
      mutate(
        genes_percent = genes_percent/(query_metadata. %>% filter(Psiblast == name_query) %>% nrow)
      ) %>% 
      setNames(c("ID", name_query))
    }
  ) %>% 
  reduce(full_join, by = "ID") %>% 
  tibble::column_to_rownames(var = "ID") %>% 
  replace(., is.na(.), 0)

colnames(psi_proxi_filt) <- psi_proxi_filt %>% colnames() %>% 
  str_replace("alginate", "Alginate") %>% 
  str_replace("cellulose1", "Cellulose I") %>% 
  str_replace("cellulose2", "Cellulose II") %>% 
  str_replace("HA_Pasteurella", "HA (pmHAS)") %>% 
  str_replace("HA_streptococcus", "HA (has)") %>% 
  str_replace("NulO_merged", "NulO") %>% 
  str_replace("pel_merged", "Pel") %>% 
  str_replace("pnag_pga", "PNAG (pga)") %>% 
  str_replace("xanthan", "Xanthan")  %>% 
  str_replace("psl", "Psl") 
  
  
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

##---------------------------------------------------------------
##                         Creating tree                         
##---------------------------------------------------------------
# Creating tibble with tree information
tree_tib <- read.tree("./data/processed/iTol/MGP1000_bac_1080_maaike_ITOL.tree") %>% 
  as_tibble() %>% 
  full_join(phylum, by = "label") %>% 
  filter(!is.na(node)) %>% 
  as.treedata() %>% as_tibble() %>% 
  # Add phylum to intermediate nodes (only those with unique offspring phylum)
  mutate(
    offspring_phylum = map(node, function(x) offspring(., x) %>% pull(phylum) %>% unique %>% `[`(!is.na(.))),
    phylum = case_when(
      offspring_phylum %>% map(length) == 0 ~ phylum,
      offspring_phylum %>% map(length) == 1 ~ as.character(offspring_phylum),
      TRUE ~ "NA"),
    phylum = ifelse(phylum == "NA", NA, phylum)
  ) %>% 
  select(-offspring_phylum) 

# Creating treedata object
tree <- tree_tib %>% 
  as.treedata()

##---------------------------------------------------------------
##                         Plotting tree                         
##---------------------------------------------------------------
tree_plot <- ggtree(tree, aes(color = phylum), layout = "rectangular", 
                    lwd = 0.1, open.angle = 20)

# # Percent identity plot
# ggsave("./figures/percID_filt_HQ_MAG.pdf", width = 10, height = 10, 
#        gheatmap(tree_plot, data = psi_perc_filt, 
#                 width = 2,colnames_offset_y = 0.5, 
#                 colnames_angle = -90, font.size = 2,
#                 hjust = 1,  border_color = NA) +
#          scale_fill_gradient(low = "white", high = "red", na.value = "white") 
# )

# Proximity filtrated plot
ggsave("./figures/trees/HQ_MAG_tree_1.pdf", width = 10, height = 15, 
  gheatmap(tree_plot, data = psi_proxi_filt, 
           width = 2,colnames_offset_y = 0.5,
           colnames_angle = -90, font.size = 2,
           hjust = 1) +
    scale_fill_gradient(low = "gray90", high = "red", na.value = "white") +
    theme(legend.position="top")
)


ggsave("./figures/trees/HQ_MAG_tree_2.pdf", width = 10, height = 30, 
       gheatmap(tree_plot, data = psi_proxi_filt, 
                width = 1.5, colnames_offset_y = 0.5,
                colnames_angle = -90, font.size = 2,
                hjust = 1) +
         scale_fill_gradient(low = "gray90", high = "red", na.value = "white") +
         theme(legend.position="top",
               line = element_line(size = 0),
               panel.grid = element_line(colour = "red"),
               rect = element_rect(colour = "blue")) +
         scale_y_continuous(expand = c(0,0))
)




## load example tree from package
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

##  circular tree
circ <- ggtree(tree, layout = "circular") +
  geom_tiplab(size=3, color="black", align=TRUE, offset = 1)

##  external data which should be attached
df <- data.frame(habitat=c("a", "b", "a", "c", "d", "d", "a", "b", "e", "e", "f", "c", "f"))
row.names(df) <- pull(circ$data[circ$data$isTip=="TRUE", "label"])

##  circular tree with heatmap without colname but remaining gap
p1 <- gheatmap(circ, df[,"habitat", drop=FALSE], offset = 0.8, width=0.1,
               colnames = FALSE) 
p1

##  circular tree with heatmap without gap
p2 <- gheatmap(circ, df[,"habitat", drop=FALSE], offset = 0.8, width=0.1,
               colnames = FALSE) +
  ##  close gap 
  scale_y_continuous(expand = c(0,0))
p2
