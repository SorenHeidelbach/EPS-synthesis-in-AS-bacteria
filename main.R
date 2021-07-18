library("here")
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
proximity_filtration("alginate", min_genes = 6)
proximity_filtration("psl", min_genes = 5)
proximity_filtration("pel_merged", min_genes = 3)
proximity_filtration("cellulose1")
proximity_filtration("cellulose2")
proximity_filtration("succinoglycan", min_genes = 6)
proximity_filtration("xanthan", min_genes = 4, )
proximity_filtration("curdlan", min_genes = 1, perc_id = 30)
proximity_filtration("pnag_pga")
proximity_filtration("pnag_ica")
proximity_filtration("pnag_eps", min_genes = 3) # changed from 2
proximity_filtration("diutan", min_genes = 6, exclude_gene = c("rmlB", "rmlC", "rmlD"))
proximity_filtration("gellan", min_genes = 6,exclude_gene = c("rmlB", "rmlC", "rmlD"))
proximity_filtration("S88", min_genes = 6, exclude_gene = c("rmlB", "rmlC", "rmlD"))
proximity_filtration("NulO_merged", min_genes = 2)
proximity_filtration("HA_Pasteurella", min_genes = 1, perc_id = 33)
proximity_filtration("HA_streptococcus", min_genes = 3, exclude_gene = c("glmU", "pgi"))

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
plot_operon("gellan")
plot_operon("S88")
plot_operon("NulO_merged")
plot_operon("HA_Pasteurella")
plot_operon("HA_streptococcus")



