library("data.table")
library("tidyverse")
library("here")
library("gggenes")
library("ggtext")
library("glue")
library("readxl")
setwd(here())


# Create dir for figure
dir.create("./figures/operon_w_domains", showWarnings = F)


plot_operon <-  function(filename_psiblast,
                         name_addon = "",
                         query_title = "Query",
                         article = FALSE,
                         article_plot_domain = FALSE,
                         mags = "all"){
##---------------------------------------------------------------
##      Parameter deifinition (must be same as in filtration)
##---------------------------------------------------------------
# The name of the .txt files
filename_psiblast_col <- paste(filename_psiblast, collapse = "_")

# Loading results from proximity filtration
genes <- fread(glue("./output/psi_operon_full/{filename_psiblast_col}.tsv"))
if(any("all" != mags)) {
  genes <- filter(genes, ID %in%  mags)
}
# File with metadata on the query genes, e.g. function
query_metadata. <- excel_sheets("./data/raw/Query_figur.xlsx") %>%
  sapply(function(X) read_xlsx("./data/raw/Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))) %>%
  rbindlist() 
##---------------------------------------------------------------
##  Loading interproscan data and merging with psiblast operons  
##---------------------------------------------------------------
clean_gff3 <- function(df){
  # Extraction of relevant annotation information
  df %>% filter(!str_detect("polypeptide", V3)) %>%
    mutate(
      domain = str_extract(V9, "signature_desc=[^;]*;"),
      domain = str_sub(domain, 1, -2),
      domain = gsub("signature_desc=", "", x = domain)) %>%
    subset(select = -c(V2, V3, V7, V8, V9)) %>%
    setNames(c("Target_label", "start1", "end1", "e_value", "Domain")) %>%
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
      Domain = str_replace(Domain, "Phosphoglucomutase/phosphomannomutase, alpha/beta/alpha domain II", "PGM_PMM_II"),
      Domain = str_replace(Domain, "Phosphoglucomutase/phosphomannomutase, alpha/beta/alpha domain I", "PGM_PMM_I"),
      Domain = str_replace(Domain, "Phosphoglucomutase/phosphomannomutase, alpha/beta/alpha domain III", "PGM_PMM_III"),
      Domain = str_replace(Domain, "Phosphoglucomutase/phosphomannomutase,_C", "PGM_PMM_IV"),
      Domain = str_replace(Domain, "RTX calcium-binding nonapeptide repeat.*", "Hemolysn_Ca-bd"),
      Domain = str_replace(Domain, "UDP-glucose/GDP-mannose dehydrogenase family, UDP binding domain", "UDPG_MGDP_dh_C"),
      Domain = str_replace(Domain, "UDP-glucose/GDP-mannose dehydrogenase family, central domain", "UDP-Glc/GDP-Man_DH_dimer"),
      Domain = str_replace(Domain, "UDP-glucose/GDP-mannose dehydrogenase family, NAD binding domain", "UDPG_MGDP_dh_N"),
      Domain = str_replace(Domain, "SGNH hydrolase-like domain, acetyltransferase AlgX", "ALGX/ALGJ_SGNH-like"),
      Domain = str_replace(Domain, "MBOAT, membrane-bound O-acyltransferase family", "MBOAT_fam"),
      Domain = str_replace(Domain, "Periplasmic copper-binding protein.*", "NosD_dom"),
      Domain = str_replace(Domain, "UDP-N-acetylglucosamine 2-epimerase", "UDP-GlcNAc_Epase")
      
    ) 
}
## Loading
test <- filename_psiblast %>% 
  lapply(function(query) {
    # Use this function on each polysaccharide name (e.g. "cellulose1)
    list.files(paste0("./data/raw/ips/", query, "_fasta_removed/")) %>%
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
  clean_gff3() %>% 
  ## Combining information from genes to domains
  full_join(genes) %>%
  mutate(start2 = start + as.numeric(start1) * 3,
         end2 = start + as.numeric(end1) * 3) %>%
  subset(select = c(
    "start", "end", "start2", "end2", "Domain", "operon", "ID",
    "ID2", "Function", "strand")) %>%
  mutate(Percent_identity = 50) %>%
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


##----------------------------------------------------------------
##                    Reversing strand direction                    
##----------------------------------------------------------------
operon_reversed <- genes %>%
  group_by(operon) %>% 
  summarize(direction = sum(strand)) %>% 
  filter(direction < 0) %>% 
  pull(operon)
  
  
genes <- genes %>% 
  group_by(operon) %>% 
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

domains <- domains %>% 
  group_by(operon) %>% 
  mutate(
    reverse = operon %in% operon_reversed,
    operon_middle = (min(c(start, end)) + max(c(start, end)))/2,
    start = ifelse(reverse, 2*operon_middle - start, start),
    end = ifelse(reverse, 2*operon_middle - end, end),
    min_operon = min(c(start, end)),
    gene_middle = (start + end)/2,
    start2 = ifelse(reverse, 2*operon_middle - start2, start2) -min_operon,
    start2 = ifelse(reverse, 2*gene_middle - start2, start2)-min_operon,
    end2 = ifelse(reverse, 2*operon_middle - end2, end2)-min_operon,
    end2 = ifelse(reverse, 2*gene_middle - end2, end2)-min_operon
    
    #strand = ifelse(reverse, strand, strand)
  ) %>% 
  select(-reverse, -operon_middle)

  

##----------------------------------------------------------------
##                    Adding query information                    
##----------------------------------------------------------------
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

# Extraction of relevant annotation information
domains <- read.table( 
    glue("./data/raw/ipsq/{filename_psiblast}.gff3"),
    sep = "\t", fill = TRUE, comment.char = "#", quote = "\""
    ) %>%
  clean_gff3() %>% 
  mutate(Genename = Target_label) %>% 
  # Adding gene information about query genes
  left_join(query_metadata.) %>% 
  mutate(
    start2 = Start + as.numeric(start1) * 3,
    end2 = Start + as.numeric(end1) * 3,
    #ID = "Query",
    ID2 = "Query",
    title = glue("**{query_title}**"),
    mi_phylum = "a",
    operon = 0,
    strand = ifelse(Strand == "+", 1, -1),
    start = Start,
    end = End
  ) %>% 
  select(c("start2", "end2", "Domain", "title",
    "ID2", "Function", "strand", "operon", "start", "end", "mi_phylum")) %>%
  mutate(Percent_identity = 50) %>%
  filter(!is.na(Domain)) %>% 
  full_join(domains) %>% 
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
assign("domains", domains, pos = 1)


if(article == FALSE){
gene_height <- 4
ggsave(
  glue("./figures/operon_w_domains/gggenes_", paste(filename_psiblast, collapse = "_"), "{name_addon}.pdf"), 
  width = unit(13, "mm"),
  height = unit(1.5 * length(unique(genes$ID2)), "mm"),
  limitsize = FALSE,
  ggplot(
    genes, aes(xmin = start, xmax = end, y = title, forward = strand)) +
    # Empty gene arrows
    geom_gene_arrow(
      arrowhead_height = unit(gene_height, "mm"),
      arrow_body_height = unit(gene_height, "mm"),
      arrowhead_width = unit(5, "mm")
    ) +
    # Colored gene arrows (match in psiblast)
    geom_subgene_arrow(
      data = genes,
      mapping = aes(
        xmin = start, xmax = end, y = title,
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
    geom_text(
      data = genes %>% mutate(start = (start_target_plot + end_target_plot)/2),
      aes(x = start, label = Query_label)
    ) +
    geom_text(
      data = genes %>% 
      mutate(
        Percent_identity = ifelse(Percent_identity == 40,
        yes = " ",
        no = paste0(signif(Percent_identity, digits = 2),"%"))),
      aes(x = (start_target_plot+end_target_plot)/2, label = Percent_identity),
      nudge_y = 0.16, size = 2
    ) +
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
      size = 1.6,
      angle = 20
    ) +
    facet_wrap( 
      ~ mi_phylum + mi_family + mi_genus + mi_species + title, 
      scales = "free", ncol = 1
    ) +
    guides(alpha = guide_legend(override.aes = list(fill = "black"))) +
    theme_genes() +
    theme(
      legend.position = "top",
      axis.text.y = element_markdown(),
      axis.title.y = element_blank()
    ) +
    scale_fill_brewer(palette = "Set3") +
    guides(
      fill = guide_legend(
        title = "Function of Matched Query Gene", 
        nrow = 1, 
        title.position = "top",
        title.hjust = 0.5)
    ) +
    scale_alpha_continuous(range = c(0.1, 1), limits = c(20, 40))
)
} else{
  ##---------------------------------------------------------------
  ##            Article plotting            
  ##---------------------------------------------------------------
  update_title = function(x){
    y <- x %>% 
      mutate(
        title = ifelse(
          ID2 == "Query", 
          title,
          glue("**{ID2}, {mi_phylum}, {mi_class}, {mi_order}, {mi_family},
                <em>{mi_genus}</em>, <em>{mi_species}</em>**")))
    return(y)
  }
  genes <- update_title(genes)
  domains <- update_title(domains)

  gene_height <- 5
  n_operon <- length(unique(genes$ID2))
  operon_plot <- ggplot(
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

      
      facet_wrap(
        ~ mi_phylum + mi_class + mi_order + mi_family + mi_genus + mi_species + title, 
        scales = "free_y", 
        ncol = 1
      ) +
      theme_genes() +
      theme(
        legend.position = "bottom",
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
  if (article_plot_domain) {
    ##---------------------------------------------------------------
    ##            Article plotting with domain            
    ##---------------------------------------------------------------
    operon_plot <- operon_plot + 
      # Domains boxes
      geom_text(
        data = genes %>% mutate(start = (start + end)/2),
        mapping = aes(x = start, label = Query_label),
        size = 4.5,
        nudge_y = -0.25, 
        angle = 0,
        hjust = 0.5
      ) +
      geom_richtext(
        data = genes %>% group_by(operon) %>% slice(1),
        mapping = aes(x = 0, label = title),
        size = 4,
        nudge_y = 0.45,
        hjust = 0,
        fill = NA, label.color = NA
      ) +
      geom_gene_arrow(
        data = domains,
        mapping = aes(
          xmin = start2,
          xmax = end2,
          y = title,
          forward = strand,
          fill = Function
        ),
        arrowhead_height = unit(gene_height - 2.5, "mm"),
        arrow_body_height = unit(gene_height - 2.5, "mm"),
        arrowhead_width = unit(0, "mm"),
        position = position_nudge(y = +0.22)
      ) +
      geom_text(
        data = domains %>% mutate(start = (start2 + end2) / 2),
        aes(x = start, label = Domain,  y = title),
        nudge_y = +0.22,
        size = 1.6,
        angle = 0
      )
    
    ggsave(
      plot = operon_plot,
      glue("./figures/operon_article/operon_", paste(filename_psiblast, collapse = "_"), "{name_addon}.pdf"), 
      width = unit(13, "mm"),
      height = unit(1 * n_operon + 1, "mm"),
      limitsize = FALSE)
    
    
  } else{
    operon_plot <- operon_plot + 
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
      )
    ggsave(
      plot = operon_plot,
      glue("./figures/operon_article/operon_", paste(filename_psiblast, collapse = "_"), "{name_addon}.pdf"), 
      width = unit(13, "mm"),
      height = unit(0.9 * n_operon + 1, "mm"),
      limitsize = FALSE)
  }
}
}


