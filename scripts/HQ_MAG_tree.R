library("data.table")
library("ggtree")
library("ggtreeExtra")
library("ggnewscale")
library("readxl")
library("treeio")
library("tidyverse")
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
        genes_percent = 100 * genes_percent/(query_metadata. %>% filter(Psiblast == name_query) %>% nrow)
      ) %>% 
      setNames(c("ID", name_query))
    }) %>% 
  reduce(full_join, by = "ID") %>% 
  tibble::column_to_rownames(var = "ID")


colnames(psi_perc_filt) <- psi_perc_filt %>% colnames() %>% 
  str_replace("alginate", "Alginate") %>% 
  str_replace("cellulose1", "Cellulose I") %>% 
  str_replace("cellulose2", "Cellulose II") %>% 
  str_replace("HA_Pasteurella", "HA (pmHAS)") %>% 
  str_replace("HA_streptococcus", "HA (has)") %>% 
  str_replace("NulO_merged", "NulO") %>% 
  str_replace("pel_merged", "Pel") %>% 
  str_replace("pnag_pga", "PNAG (pga)") %>% 
  str_replace("pnag_eps", "PNAG (eps)") %>% 
  str_replace("xanthan", "Xanthan") %>% 
  str_replace("psl", "Psl") 

# Percentage of proximity filtrated genes in each HQ-MAG
psi_proxi_filt <- list.files("./output/psi_proxi_filt/") %>% 
  map(.f = function(query){
    name_query <- tools::file_path_sans_ext(query)
    fread(paste0("./output/psi_proxi_filt/", query)) %>% 
      group_by(ID) %>% 
      summarise(genes_percent = length(unique(Query_label))) %>% 
      mutate(
        genes_percent = 100 * genes_percent/(query_metadata. %>% filter(Psiblast == name_query) %>% nrow)
      ) %>% 
      setNames(c("ID", name_query))
    }
  ) %>% 
  reduce(full_join, by = "ID") %>% 
  tibble::column_to_rownames(var = "ID")

colnames(psi_proxi_filt) <- psi_proxi_filt %>% colnames() %>% 
  str_replace("alginate", "Alginate") %>% 
  str_replace("cellulose1", "Cellulose I") %>% 
  str_replace("cellulose2", "Cellulose II") %>% 
  str_replace("HA_Pasteurella", "HA (pmHAS)") %>% 
  str_replace("HA_streptococcus", "HA (has)") %>% 
  str_replace("NulO_merged", "NulO") %>% 
  str_replace("pel_merged", "Pel") %>% 
  str_replace("pnag_pga", "PNAG (pga)") %>% 
  str_replace("pnag_eps", "PNAG (eps)") %>% 
  str_replace("xanthan", "Xanthan") %>% 
  str_replace("psl", "Psl") 

psi_proxi_filt <- psi_proxi_filt %>% 
  mutate(ID = row.names(.)) %>% 
  pivot_longer(-ID)
psi_perc_filt <- psi_perc_filt %>% 
  mutate(ID = row.names(.)) %>% 
  pivot_longer(-ID)

psi_proxi_filt <- psi_proxi_filt %>% 
  rbind(
psi_perc_filt %>% group_by(name) %>% head(n=16) %>% mutate(value = NA) %>% filter(!(name %in% unique(psi_proxi_filt$name)))
)
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
      TRUE ~ "NA",
    )
  ) %>% 
  select(-offspring_phylum) 

# Creating treedata object
tree <- tree_tib %>%
  as.treedata()

##---------------------------------------------------------------
##                         Plotting tree                         
##---------------------------------------------------------------
# Tree with inspiration from https://yulab-smu.top/treedata-book/chapter10.html
fruit_list <- geom_fruit_list(
    geom_fruit(
      data = psi_perc_filt, 
      geom = geom_tile,
      pwidth = 0.75,
      mapping = aes(y = ID, x = name, fill = value),
      axis.params = list(axis = "x", text.size = 2, text.angle = -90, hjust = 0, vjust = 0.5),               # axis.params = list(axis = "x", text.size = 2.5, text.angle = -90, hjust = 0, vjust = 0.5),
      grid.params = list(vline = FALSE, color = "gray60", alpha = 0.3),
      offset = 0.1
    ),
    scale_fill_stepsn(
      colors = c("#B7C3FB", "#9AABF9", "#6E87F7", "#516FF5"), 
      na.value = "transparent",  
      limits = c(0,100),
      guide = guide_colorbar(
        title = "Percentage of Query Genes Present\nBefore Proxi Filtration",
        title.position = "top",
        title.hjust = 0.5,
        title.vjust = 0.5,
        ticks.colour = "black", 
        frame.colour = "black")
    ),
    new_scale_fill(),
    geom_fruit(
      data = psi_proxi_filt, 
      geom = geom_tile,
      pwidth = 0.75,
      mapping = aes(y = ID, x = name, fill = value),
      axis.params = list(axis = "x", text.size = 2, text.angle = -90, hjust = 0, vjust = 0.5),               # axis.params = list(axis = "x", text.size = 2.5, text.angle = -90, hjust = 0, vjust = 0.5),
      grid.params = list(vline = FALSE, color = "gray60", alpha = 0.3),
      offset = 0.1
    ),
    scale_fill_stepsn(
      colours = c("#FF9999", "#FF6666", "#FF3333", "#FF0000"),
      na.value = "transparent",  
      limits = c(0,100),
      guide = guide_colorbar(
        title = "Percentage of Query Genes Present\nAfter Proxi Filtration",
        title.position = "top",
        title.hjust = 0.5,
        title.vjust = 0.5, 
        ticks.colour = "black", 
        frame.colour = "black")
    )
  )
# circular
phylum_displayed <- phylum$phylum %>% table %>% `[`(order(.)) %>% tail(8) %>% names()
tree_plot <- 
  ggtree(
    tree, 
    layout = "fan", 
    lwd = 0.1, 
    open.angle = 20
  ) +
  geom_hilight(
    data = filter(tree_tib, phylum_ancestor != "NA") %>% 
      filter(phylum %in% phylum_displayed), 
    aes(node = node, fill = phylum_ancestor),
    extendto = 2.58, alpha = 0.4) +
  geom_cladelab(
    data = filter(tree_tib, phylum_ancestor != "NA") %>% 
      mutate(phylum = str_replace_all(phylum, "Ca_", "Candidatus ")) %>% 
      filter(phylum %in% phylum_displayed), 
    mapping = aes(node = node, label = phylum),
    angle = "auto", barsize = NA,
    geom = "text", fontsize = 1.8, horizontal = TRUE, 
    align = TRUE, fontface  = 0.8,  hjust = 1, offset.text = 0.28
    ) +
  scale_color_manual(guide = "none") +
  scale_fill_manual(
    values = c(
      "#FFC125","#87CEFA","#7B68EE","#808080","#800080",
      "#006400","#800000","#B0171F","#191970", "#006400",
      "#800000","#B0171F","#191970"), 
    guide = "none") +
  new_scale_fill() +
  fruit_list +
  theme(legend.position = "bottom")


ggsave("./figures/trees/HQ_MAG_tree_fan.pdf", width = 10, height = 10, limitsize = FALSE,
       plot = tree_plot
)


##---------------------------------------------------------------
##                         Trees of genes in each system                         
##---------------------------------------------------------------

plot_genes_HQ_mag <- function(eps){
  data <- fread(glue("./output/psi_percID_filt/{eps}.tsv")) %>% 
    select(Query_label, ID, Percent_identity)
  tree_plot <- ggtree(tree, layout = "fan", 
                      lwd = 0.2, open.angle = 20, ) +
    geom_cladelab(data = filter(tree_tib, phylum_ancestor != "NA") %>% mutate(phylum = str_replace_all(phylum, "Ca_", "Candidatus ")), 
                  mapping = aes(node = node, label = phylum),
                  angle = "auto", barsize = NA,
                  geom = "text", fontsize = 2.3, horizontal = TRUE, 
                  align = TRUE, fontface  = 0.8, offset = 0.08,
                  hjust = 1
    ) +
    scale_color_manual(guide = "none") +
    geom_hilight(data = filter(tree_tib, phylum_ancestor != "NA") , 
                 aes(node = node, fill = phylum_ancestor),
                 extendto = 2.37, alpha = 0.3) +
    scale_fill_manual(values = c("#FFC125","#87CEFA","#7B68EE","#808080","#800080",
                                 "#9ACD32","#D15FEE","#FFC0CB","#EE6A50","#8DEEEE",
                                 "#006400","#800000","#B0171F","#191970", "#006400","#800000","#B0171F","#191970"), guide = "none") +
    new_scale_fill() +
    geom_fruit(data = data, 
               geom = geom_tile,
               pwidth = 0.75,
               mapping = aes(y = ID, x = Query_label, fill = Percent_identity),
               lwd = 0.0001,
               axis.params = list(axis = "x", text.size = 2.5, text.angle = -90, hjust = 0, vjust = 0.5),
               grid.params = list(vline = FALSE, color = "gray60", alpha = 0.3),
               offset = 0.1) +
    scale_fill_gradientn(colors = c("lightblue2", "slateblue", "red3"), 
                         na.value = "transparent",  
                         limits = c(20,40),
                         guide = "colourbar") +
    scale_alpha_continuous(range = c(0,1)) +
    theme(legend.position = "top",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(5, 0, -80, 0)) +
    guides(fill = guide_colorbar(title = "Percentage of Query Genes Present", override.aes = list(col = "black"), ), #, fill = c("mediumpurple4", "orange2", "mediumpurple4", "orange2", "mediumpurple4"))),
           alpha = guide_legend(title = "Filtration Step")) 
  
  ggsave(glue("./figures/trees/HQ_MAG_tree_fan_percid_{eps}.pdf"), width = 12, height = 10, limitsize = FALSE,
         plot = tree_plot
  )
}
c("alginate", "cellulose1", "cellulose2", 
  #"HA_Pasteurella", 
  "HA_streptococcus", "NulO_merged", "pel_merged", "pnag_pga", "pnag_eps", "xanthan", "psl") %>% map(plot_genes_HQ_mag)
plot_genes_HQ_mag("cellulose1")
