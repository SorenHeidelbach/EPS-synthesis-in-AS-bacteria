library("data.table")
library("tidyverse")
library("tidyr")
library("treeio")
library("ggtree")
setwd(here::here())

phylum <- read_tsv("./data/raw/magstats.tsv") %>% 
  mutate(
    phylum =   str_extract(MiDAS3_6Tax, "p:[^,]*"),
    phylum =   str_replace(phylum, "p:",replacement = "")
  ) %>% 
  select(ID, phylum) %>% 
  setNames(c("label", "phylum"))



psi_proxi_filt <- list.files("./output/psi_proxi_filt/") %>% 
  map(.f = function(x){
    fread(paste0("./output/psi_proxi_filt/", x)) %>% 
      group_by(ID) %>% 
      summarise(genes_percent = length(unique(Query_label))) %>% 
      setNames(c("ID", tools::file_path_sans_ext(x)))
    }
  ) %>% 
  reduce(full_join, by = "ID") %>% 
  tibble::column_to_rownames(var = "ID")

data <- read.tree("./data/processed/iTol/MGP1000_bac_1080_maaike_ITOL.tree") %>% 
  as_tibble() %>% 
  full_join(phylum, by = "label") %>% 
  select(node, phylum) %>% 
  as_tibble()

tree <- read.tree("./data/processed/iTol/MGP1000_bac_1080_maaike_ITOL.tree") %>% 
  tidytree::treedata(phylo = ., data = data)

tree_plot <- ggtree(tree, aes(color = phylum))
tree_plot <- open_tree(tree_plot, 20, )
ggsave("./figures/test.pdf", width = 10, height = 20, 
gheatmap(tree_plot, data = psi_proxi_filt, 
         width = 1,
         colnames_angle=90, font.size = 2,
         hjust = 1) +
  ylim(0, 1100) 
)




