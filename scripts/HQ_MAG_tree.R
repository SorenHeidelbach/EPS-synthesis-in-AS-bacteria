library("data.table")
library("tidyverse")
library("tidyr")
library("treeio")
library("ggtree")
library("ggtreeExtra")
library("ggnewscale")
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
    phylum =   str_remove(phylum, "p:"),
    genus =  str_extract(MiDAS3_6Tax, "g:[^,]*"),
    genus =  str_remove(genus, "g:"),
    genus =  str_remove(genus, ";")
  ) %>% 
  select(ID, phylum, genus) %>% 
  setNames(c("label", "phylum", "genus"))

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
    parent_offspring_phylum = map(parent, function(x) {
        offspring(., x) %>% 
        pull(phylum) %>% 
        unique() %>% 
        `[`(!is.na(.))
      }),
    phylum = case_when(
      offspring_phylum %>% map(length) == 0 ~ phylum,
      offspring_phylum %>% map(length) == 1 ~ as.character(offspring_phylum),
      TRUE ~ "NA"),
    phylum = ifelse(phylum == "NA", NA, phylum),
    phylum_ancestor = case_when(
      parent_offspring_phylum %>% map(length) > 1 & offspring_phylum %>% map(length) == 1 ~ as.character(offspring_phylum),
      TRUE ~ "NA"
    )
  ) %>% 
  select(-offspring_phylum) 

# Creating treedata object
tree <- tree_tib %>%
  as.treedata()

##---------------------------------------------------------------
##                         Plotting tree                         
##---------------------------------------------------------------

tree_plot <- ggtree(tree, layout = "fan", 
                    lwd = 0.2, open.angle = 20, ) +
  geom_cladelab(data = filter(tree_tib, phylum_ancestor != "NA") %>% mutate(node2 = node), 
                  mapping = aes(node = node2, label = phylum, color = phylum),
                  geom = "text", fontsize = 2.5, horizontal = TRUE, align = FALSE, fontface  = 2) +
  scale_color_discrete(guide = "none") +
  geom_hilight(data = filter(tree_tib, phylum_ancestor != "NA") , aes(node = node, fill = phylum_ancestor)) +
  scale_fill_discrete(guide = "none") +
  new_scale_fill() +
  geom_fruit(data = psi_proxi_filt %>% mutate(label = row.names(.)) %>%  pivot_longer(cols = -label), 
             geom = geom_tile,
             pwidth = 0.75,
             mapping = aes(y = label, x = name, fill = value),
             col = "grey94",
             alpha = 1,
             lwd = 0.0001,
             axis.params = list(axis = "x", text.size = 5, text.angle = 90, hjust = 1, vjust = 0.5),
             grid.params = list(vline = TRUE, color = "grey94", hline = FALSE),
             offset = 0.1) +
  scale_fill_gradient(low = "grey94", high = "red", na.value = "white",  guide = guide_legend(title = "Percent identified Genes", order = 1)) +
  theme(legend.position = "bottom")



ggsave("./figures/trees/HQ_MAG_tree_fruit_fan.pdf", width = 12, height = 16, limitsize = FALSE,
  plot = tree_plot
)

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


ggsave("./figures/trees/HQ_MAG_tree_2.pdf", width = 10, height = 50, limitsize = FALSE,
       gheatmap(tree_plot, data = psi_proxi_filt, 
                width = 1.5, colnames_offset_y = 0.5,
                colnames_angle = -90, font.size = 2,
                hjust = 1) +
         scale_fill_gradient(low = "gray90", high = "red", na.value = "white") +
         theme(legend.position="top")
)



