library("tidyverse")
library("here")
library("gggenes")
library("ggtext")
library("glue")
setwd(here())


# Create dir for figure
dir.create("./figures/operon_w_domains", showWarnings = F)


plot_operon <-  function(filename_psiblast, name_addon = ""){
##---------------------------------------------------------------
##      Parameter deifinition (must be same as in filtration)
##---------------------------------------------------------------
# The name of the .txt files
filename_psiblast_col <- paste(filename_psiblast, collapse = "_")

# Loading results from proximity filtration
genes <- fread(glue("./output/psi_operon_full/{filename_psiblast_col}.tsv"))

##---------------------------------------------------------------
##  Loading interproscan data and merging with psiblast operons  
##---------------------------------------------------------------
## Loading
filename_psiblast %>% 
  lapply(function(query) {
    # Use this function on each polysaccharide name (e.g. "cellulose1)
    list.files(paste0("./data/raw/ips/", query, "_fasta_removed/")) %>%
      subset(grepl(., pattern = query)) %>% 
      lapply(function(mag) {
        # Use this function on each mag ID (e.g. "Ega_18-Q3-R5-49_MAXAC.199")
        mag_path = paste0("./data/raw/ips/", query, "_fasta_removed/", mag)
        if (file.info(mag_path)$size > 0) {
          read.table(
            mag_path,
            sep = "\t",
            fill = TRUE,
            comment.char = "#",
            quote = "\""
          )
        }
      }) %>%
      bind_rows()
  }) %>%
  bind_rows() %>%
  # Extraction of relevant annotation information
  filter(!str_detect("polypeptide", V3)) %>%
  mutate(
    domain = str_extract(V9, "signature_desc=[^;]*;"),
    domain = str_sub(domain, 1, -2),
    domain = gsub("signature_desc=", "", x = domain)) %>%
  subset(select = -c(V2, V3, V7, V8, V9)) %>%
  setNames(c("Target_label", "start1", "end1", "e_value", "Domain")) %>%
  ## Combining information from genes to domains
  full_join(genes) %>%
  mutate(start2 = start + as.numeric(start1) * 3,
         end2 = start + as.numeric(end1) * 3) %>%
  subset(select = c(
    "start", "end", "start2", "end2", "Domain", "operon", "ID",
    "ID2", "Function", "strand")) %>%
  mutate(Percent_identity = 50) %>%
  filter(!is.na(Domain)) %>%
  distinct() %>%
  # Formating of domain names (a bit of confusing regular expressions)
  mutate(
    Domain = str_replace(Domain, "[gG]lycosyl.*transferase.*[fF]amily ", "GT f"),
    Domain = str_replace(Domain, "[gG]lycosyl.*transferase.*[gG]roup ", "GT g"),
    Domain = str_replace(Domain, "[gG]lycosyl.*transferase.[lL]ike.[fF]amily", "GT like f"),
    Domain = str_replace(Domain, "[gG]lycosyl.*transferase", "GT"),
    Domain = str_replace(Domain, "[gG]lycosyl [hH]ydrolase.*[fF]amily " , "GH"),
    Domain = str_replace(Domain, "[gG]lycosyl [hH]ydrolase" , "GH"),
    Domain = str_replace(Domain, ".*[cC]ellulose.*synth.*protein[^' ']*" , "CS "),
    Domain = str_replace(Domain, ".*[cC]ellulose.*synth[^' ']*", "CS "),
    Domain = str_replace(Domain, " N.terminal domain" , "_N"),
    Domain = str_replace(Domain, " C.terminal domain" , "_C"),
    Domain = str_replace(Domain, ".*BCSC_C.*", "bcsC"),
    Domain = str_replace(Domain, "[iI]nitation [fF]actor" , "IF"),
    Domain = str_replace(Domain, "[eE]longation [fF]actor" , "EF"),
    Domain = str_replace(Domain, ".*Tetratrico.*|.*TPR.*", "T"),
    Domain = str_replace(Domain, ".*[dD]omain of unknown function.*", "NA"),
  ) %>% 
  # -> ASSIGNING "domains"
  assign(x = "domains", value = ., pos = 1)

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
      title = paste0(ID2, "<br>",
                     mi_phylum, "<br>",
                     mi_class, "<br>",
                     mi_order, "<br>",
                     mi_family, "<br>",
                     "*", mi_genus, "*", "<br>",
                     "*", mi_species, "*"),
      # Formating taxa names to be more inline with recommended guidelines
      title = str_replace_all(title, 
                              pattern = "\\*Ca_([^*]*)\\*", 
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
domains <- add_midas_tax(domains) %>% 
  filter(!is.na(ID2))



##---------------------------------------------------------------
##            Plotting operons with domain annotation            
##---------------------------------------------------------------
gene_height <- 4
ggsave(glue("./figures/operon_w_domains/gggenes_", paste(filename_psiblast, collapse = "_"), "{name_addon}.pdf"), 
       width = unit(13, "mm"),
       height = unit(1.5*length(unique(genes$ID2)), "mm"),
       limitsize = FALSE,
       ggplot(genes, aes(xmin = start, xmax = end, y = title, forward = strand)) +
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
                         forward = strand,
                         alpha = Percent_identity
           ),
           arrowhead_height = unit(gene_height, "mm"),
           arrow_body_height = unit(gene_height, "mm"),
           arrowhead_width = unit(5, "mm"),
           #position = position_nudge(y = 0.3)
         ) +
         geom_text(data = genes %>% mutate(start = (start_target_plot + end_target_plot)/2),
                   aes(x = start, label = Query_label)) +
         geom_text(data = genes %>% mutate(Percent_identity = ifelse(Percent_identity == 40,
                                                                     yes = " ",
                                                                     no = paste0(
                                                                       signif(Percent_identity, digits = 2),
                                                                       "%"))),
                   aes(x = start_target_plot, label = Percent_identity),
                    nudge_y = 0.16, size = 2) +
         # Domains boxes
         geom_gene_arrow(
           data = domains,
           mapping = aes(
             xmin = start2,
             xmax = end2,
             y = title,
             forward = strand,
             fill = Function,
             alpha = Percent_identity
           ),
           arrowhead_height = unit(gene_height - 1, "mm"),
           arrow_body_height = unit(gene_height - 1, "mm"),
           arrowhead_width = unit(0, "mm"),
           position = position_nudge(y = -0.35)
         ) +
         geom_text(
           data = domains %>% mutate(start = (start2 + end2) / 2),
           aes(x = start, label = Domain,  y = title),
           nudge_y = -0.35,
           size = 2,
           angle = 8
         ) +
         facet_wrap( ~ operon, scales = "free", ncol = 1) +
         scale_alpha_continuous(range = c(0.1, 1), limits = c(20, 40)) +
         guides(alpha = guide_legend(override.aes = list(fill = "black"))) +
         theme_genes() +
         theme(
           legend.position = "top",
           axis.text.y = element_markdown(),
           axis.title.y = element_blank()
         ) +
         scale_fill_brewer(palette = "Set3")
)
}


