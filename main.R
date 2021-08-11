setwd(here::here())

dir.create("./data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("./data/processed", recursive = TRUE, showWarnings = FALSE)

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


# The main proximity filtration pipeline
source("./scripts/proximity_filtration.R")
# Default: min_genes=2 and perc_id=20
proximity_filtration("alginate", 
                     min_genes = 6)
proximity_filtration("psl", 
                     min_genes = 7)
proximity_filtration("pel_merged", 
                     min_genes = 6)
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



