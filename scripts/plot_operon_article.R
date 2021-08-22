library("data.table")
library("tidyverse")
library("here")
library("gggenes")
library("ggtext")
library("glue")
library("readxl")
setwd(here())


# Create dir for figure
dir.create("./figures/operon_article", showWarnings = F)


query_title_test <- 
plot_operon_article <-  function(filename_psiblast,
                                 name_addon = "", 
                                 query_title = "Query",
                                 mags = "all"){
  ##---------------------------------------------------------------
  ##      Parameter deifinition (must be same as in filtration)
  ##---------------------------------------------------------------
  # The name of the .txt files
  filename_psiblast_col <- paste(filename_psiblast, 
                                 collapse = "_")
  
  # Loading results from proximity filtration
  genes <- fread(glue("./output/psi_operon_full/{filename_psiblast_col}.tsv"))
  if(any("all" != mags)) {
    genes <- filter(genes, ID %in%  mags)
  }
  ##----------------------------------------------------------------------
  ##  Adding midas taxonomy names and modify appearance to be more neat   
  ##----------------------------------------------------------------------
  add_midas_tax <- function(data) {
    genes %>% 
      select(ID2, MiDAS3_6Tax) %>% 
      separate(MiDAS3_6Tax, into = c("drop", "MiDAS3_6Tax"), sep = "=") %>% 
      select(-drop) %>% 
      mutate(
        MiDAS3_6Tax = str_remove_all(MiDAS3_6Tax, ".[\\:\\;]")
      ) %>% 
      separate(MiDAS3_6Tax, into = c("mi_domain","mi_phylum", "mi_class", "mi_order", "mi_family", "mi_genus", "mi_species"), sep = ",") %>% 
      distinct() %>% 
      right_join(data) %>% 
      mutate(
        mi_species = str_remove(mi_species, paste0(mi_genus, "_")),
        title = paste0("**", ID2, ", ",
                       mi_phylum, ", ",
                       mi_class, ", ",
                       mi_order, ", ",
                       mi_family, ", ",
                       "<em>", mi_genus, "</em>", ", ",
                       "<em>", mi_species, "</em>",
                       "**"),
        # Formating taxa names to be more inline with recommended guidelines
        title = str_replace_all(title, 
                                pattern = "<em>Ca_([^/]*)</em>", 
                                replacement = "*Candidatus* \\1"),
        title = str_replace_all(title, 
                                pattern = "\\Ca_(.*) ", 
                                replacement = "Candidatus \\1"),
        title = str_replace_all(title, 
                                pattern = "\\*(.*)\\_marine\\_group\\*", 
                                replacement = "\\1 marine group"),
        title = str_replace_all(title, 
                                pattern = "\\*(.*)\\_marine\\_group\\*", 
                                replacement = "\\1 marine group"),
        title = str_replace_all(title, 
                                pattern = "_Subgroup_(.)", 
                                replacement = " (Subgroup \\1)"),
      )
  }
  genes <- add_midas_tax(genes) %>% 
    mutate(Query_label = replace_na(Query_label, " "))
  
  
  ##----------------------------------------------------------------
  ##                    Reversing strand direction                    
  ##----------------------------------------------------------------
  operon_reversed <- genes %>%
    group_by(operon) %>% 
    summarize(direction = sum(strand)) %>% 
    filter(direction < 0) %>% 
    pull(operon)
  
  
  genes <- genes %>% group_by(operon) %>% 
    mutate(
      reverse = operon %in% operon_reversed,
      operon_middle = (min(c(start, end)) + max(c(start, end)))/2,
      start = ifelse(reverse, 2*operon_middle - start, start),
      end = ifelse(reverse, 2*operon_middle - end, end),
      min_operon = min(c(start, end)),
      start_target_plot = ifelse(reverse, 2*operon_middle - start_target_plot, start_target_plot) - min_operon,
      end_target_plot = ifelse(reverse, 2*operon_middle - end_target_plot, end_target_plot) - min_operon,
      end = end - min_operon,
      start = start - min_operon
      #strand = ifelse(reverse, strand, strand)
    ) %>% 
    select(-reverse, -operon_middle)

  ##----------------------------------------------------------------
  ##                    Adding query information                    
  ##----------------------------------------------------------------
  # File with metadata on the query genes, e.g. function
  query_metadata. <- excel_sheets("./data/raw/Query_figur.xlsx") %>%
    sapply(function(X) read_xlsx("./data/raw/Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
    lapply(as.data.frame) %>% 
    `[`(!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))) %>%
    rbindlist() 
  genes <- query_metadata. %>% 
    mutate(
      start = Start,
      start_target_plot = Start,
      end = End,
      end_target_plot = End,
      title = glue("**{query_title}**"),
      ID2 = "Query",
      mi_phylum = "a",
      operon = 0,
      strand = ifelse(Strand == "+", 1, -1),
      Query_label = Genename,
    ) %>% 
    filter(Psiblast %in% filename_psiblast) %>% 
    full_join(genes) %>% 
    mutate(
      Function = str_replace(Function, "MOD", "Modification"),
      Function = str_replace(Function, "PE", "Polymerization & Export"),
      Function = str_replace(Function, "GT", "Glycosyl Transferase"),
      Function = str_replace(Function, "ABC", "ABC transporter"),
      Function = str_replace(Function, "SY", "Synthase"),
      Function = str_replace(Function, "PS", "Precourser Synthesis"),
      Function = str_replace(Function, "REG", "Regulatory")
    )

  ##---------------------------------------------------------------
  ##            Plotting operons with domain annotation            
  ##---------------------------------------------------------------
  assign("genes", genes, pos = 1)
  
  gene_height <- 5
  n_operon <- length(unique(genes$ID2))
  ggsave(
    glue("./figures/operon_article/operon_", paste(filename_psiblast, collapse = "_"), "{name_addon}.pdf"), 
    width = unit(13, "mm"),
    height = unit(0.9 * n_operon + 1, "mm"),
    limitsize = FALSE,
    ggplot(
      genes, aes(xmin = start, xmax = end, y = title, forward = strand)
      ) +
      # Empty gene arrows
      geom_gene_arrow(
        arrowhead_height = unit(gene_height, "mm"),
        arrow_body_height = unit(gene_height, "mm"),
        arrowhead_width = unit(5, "mm")
      ) +
     # Colored gene arrows (match in psiblast)
      geom_subgene_arrow(
       data = genes,
       mapping = aes(xmin = start, xmax = end, y = title,
                     xsubmin = start_target_plot,
                     xsubmax = end_target_plot,
                     fill = Function,
                     forward = strand
       ),
       arrowhead_height = unit(gene_height, "mm"),
       arrow_body_height = unit(gene_height, "mm"),
       arrowhead_width = unit(5, "mm")
       #position = position_nudge(y = 0.3)
      ) +
      geom_text(
        data = genes %>% mutate(start = (start + end)/2),
        mapping = aes(x = start, label = Query_label),
        size = 4.5,
        nudge_y = -0.31, 
        angle = 0,
        hjust = 0.5
      ) +
      geom_richtext(
        data = genes %>% group_by(operon) %>% slice(1),
        mapping = aes(x = 0, label = title),
        size = 4,
        nudge_y = 0.33,
        hjust = 0,
        fill = NA, label.color = NA
      ) +
      facet_wrap(
        ~ mi_phylum + mi_class + mi_order + mi_family + mi_genus + mi_species + title, 
        scales = "free_y", 
        ncol = 1
      ) +
      guides(alpha = guide_legend(override.aes = list(fill = "black"))) +
      theme_genes() +
      theme(
       legend.position = "top",
       legend.spacing.x = unit(6, "mm"),
       legend.text = element_text(margin = margin(b  = -15)),
       axis.text.y = element_blank(),
       axis.title.x = element_blank(),
       axis.title.y = element_blank(), plot.margin = unit(c(0,0,0,0), "mm")
      ) + 
      scale_fill_brewer(palette = "Set3", na.value = "white") +
      guides(
        fill = guide_legend(
          title = "Function of Matched Query Gene", 
          nrow = 1, 
          title.position = "top", 
          title.hjust = 0.5, 
          label.vjust = 1,
          label.theme = element_text(size = 10),
          keyheight = 0.9,
          label.position = "bottom")
      ) +
      scale_x_continuous(expand = rep(max(genes$end, genes$start) * 0.0000003, 2))
  )
}


