setwd(here::here())

dir.create("./data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("./data/processed", recursive = TRUE, showWarnings = FALSE)


# Checking if the correct files are present, and how to place them if not
if(!"gff.tsv" %in% list.files("./data/raw/")) {
  if(!"gff_files_reduced_out" %in% list.files("./data/raw/")){
    print("Download/move 'gff_files_reduced_out' folder to 'data/raw/'")
  } else {
    print('Please wait while gff.tsv is generated')
    source("./scripts/generate_gff.R")
    rm(gff)
  }
}
if(!"magstats.tsv" %in% list.files("./data/raw/")) {
  if(!"MAG_statistics_STABLEX_20200213.tsv" %in% list.files("./data/raw/")){
    print("Download/move 'MAG_statistics_STABLEX_20200213.tsv' to 'data/raw/'")
  } else {
    source("./scripts/generate_magstats.R")
    rm(magstats)
  }
}
if(!"Query_figur.xlsx" %in% list.files("./data/raw/")) "Download/move 'Query_figur.xlsx' to data/raw/"
if(!"MGP1000_HQMAG1083_prot_db_split" %in% list.files("./data/raw/")) "Download/move 'MGP1000_HQMAG1083_prot_db_split' folder to data/raw/"
if(!"psiblast" %in% list.files("./data/raw/")) "Download 'psiblast' with psiblast results to data/raw/"








# The main proximity filtration
source("./scripts/proximity_filtration.R")
# Default: min_genes=2 and perc_id=20
proximity_filtration("alginate", 
                     min_genes = 6)
proximity_filtration("psl", 
                     min_genes = 7)
proximity_filtration("pel_merged", 
                     min_genes = 6,
                     essential_genes = "pelF")
proximity_filtration("cellulose1", 
                     min_genes = 2,
                     essential_genes = "bcsAI")
proximity_filtration("cellulose2",
                     min_genes = 2,
                     essential_genes = "bcsABII-A")
proximity_filtration("succinoglycan", 
                     min_genes = 9)
proximity_filtration("xanthan", 
                     min_genes = 6)
proximity_filtration("curdlan", 
                     min_genes = 2, 
                     essential_genes = "crdS")
proximity_filtration("pnag_pga",
                     min_genes = 3,
                     essential_genes = "pgaC")
proximity_filtration("pnag_ica",
                     min_genes = 3,
                     essential_genes = "icaA")
proximity_filtration("pnag_eps", 
                     min_genes = 3,
                     essential_genes = c("epsH", "epsJ"))
proximity_filtration("diutan", 
                     min_genes = 10, 
                     exclude_gene = c("rmlA", "rmlB", "rmlC", "rmlD"))
# proximity_filtration("gellan", 
#                      min_genes = 10,
#                      exclude_gene = c("rmlA", "rmlB", "rmlC", "rmlD",
#                                       "gelG", "gelS", "gelR", "gelA"))
proximity_filtration("S88", 
                     min_genes = 9, 
                     exclude_gene = c("rmlA", "rmlB", "rmlC", "rmlD",
                                      "rhsB", "rhsA", "rhsC", "rhsB"))
proximity_filtration("NulO_merged", 
                     min_genes = 3,
                     essential_genes = c("neuA", "neuB"))
proximity_filtration("HA_Pasteurella", 
                     min_genes = 1, 
                     perc_id = 33)
proximity_filtration("HA_streptococcus", 
                     min_genes = 3, 
                     perc_id = 20, 
                     essential_genes = "hasA",
                     exclude_gene = c("glmU", "pgi"))










# Plotting of operons from results, remember to run ips berfore this
if(!"ips" %in% list.files("./data/raw/")) "Download/move 'ips' folder with interproscan results to data/raw/"
source("./scripts/plot_operon.R")

plot_operon("alginate")
plot_operon("psl")
plot_operon("pel_merged")
plot_operon("cellulose1")
plot_operon("cellulose2")
plot_operon("succinoglycan")
plot_operon("xanthan")
plot_operon("curdlan")
plot_operon("pnag_pga")
plot_operon("pnag_ica")
plot_operon("pnag_eps")
plot_operon("diutan")
# plot_operon("gellan")
plot_operon("S88")
plot_operon("NulO_merged")
plot_operon("HA_Pasteurella")
plot_operon("HA_streptococcus")


source("./scripts/plot_operon_article.R")

plot_operon_article(
  "alginate",
  query_title = "Query, 
                Proteobacteria, Gammaproteobacteria, 
                Pseudomonadales, Pseudomonadaceae, 
                <em>Pseudomonas</em>, <em>aeruginosa</em>",
  mags = c("Aved_18−Q3−R54−62_MAXAC.392",
           "Rand_18−Q3−R56−63_BAT3C.326"))
#plot_operon_article("psl")
plot_operon_article(
  "pel_merged",
  query_title = "Query, 
                Proteobacteria, Gammaproteobacteria, 
                Pseudomonadales, Pseudomonadaceae, 
                <em>Pseudomonas</em>, <em>aeruginosa</em>",
  mags = c("Hade_18−Q3−R52−61_MAXAC.304",
           "Hade_18−Q3−R52−61_BATAC.311",
           "Kalu_18−Q3−R12−55_BATAC.288",
           "Aved_18−Q3−R54−62_BAT3C.540"))
plot_operon_article(
  "cellulose1", 
  query_title = "Query, 
                Proteobacteria, Alphaproteobacteria, 
                Rhodospirillales, Acetobacteraceae, 
                <em>Komagataeibacter</em>, <em>xylinus</em>",
  mags = c("Hjor_18−Q3−R7−51_BAT3C.155",
           "Kalu_18−Q3−R12−55_BAT3C.261",
           "Skiv_18−Q3−R9−52_BATAC.396",
           "Aved_18−Q3−R54−62_BAT3C.394",
           "Hade_18−Q3−R52−61_BATAC.364",
           "EsbW_18−Q3−R4−48_BATAC.453",
           "OdNE_18−Q3−R46−58_BATAC.187"))
#plot_operon_article("cellulose2")
#plot_operon_article("succinoglycan")
#plot_operon_article("xanthan")
#plot_operon_article("curdlan")
plot_operon_article(
  "pnag_pga",
  query_title = "Query,
                Proteobacteria,
                Gammaproteobacteria,
                Enterobacterales,
                Enterobacteriaceae,
                <em>Escherichia</em>,
                <em>coli</em>",
  mags = c("Kalu_18−Q3−R12−55_BAT3C.208",
           "OdNE_18−Q3−R46−58_BAT3C.415"))
#plot_operon_article("pnag_ica")
#plot_operon_article("pnag_eps")
#plot_operon_article("diutan")
## plot_oper_articleon("gellan")
#plot_operon_article("S88")
plot_operon_article(
  "NulO_merged",
  mags = c("EsbW_18−Q3−R4−48_MAXAC.279_cln",
           "Ega_18−Q3−R5−49_MAXAC.062",
           "EsbW_18−Q3−R4−48_MAXAC.012",
           "OdNE_18-Q3-R46-58_MAXAC.071"))
#plot_operon_article("HA_Pasteurella")
#plot_operon_article("HA_streptococcus")

