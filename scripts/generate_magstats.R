library("here")
library("dplyr")
setwd(here())


#----------------------------------------------------------------
#  Loading, cleaning, filtering and saving HQ-MAG statistics data   
#----------------------------------------------------------------
magstats = read.csv(file = "data/raw/MAG_statistics_STABLEX_20200213.tsv", header = TRUE, sep = "\t")
magstats <- magstats %>%
  separate(MAG, c("ID", "rm"), sep =-3, remove = TRUE) %>% 
  select(!rm) %>% 
  select(!c("MaxContigBP", "AvContigBP", "HQMAG", "HQdRep", "HQdRep99ANI", "TotBP", "NumContigs", 
          "HQSpRep", "Comp", "Cont", "Circ", "StrHet", "FLSSU", "FLLSU", "ilmcov", "npcov",
          "mm27f", "mm534r", "ttl_polymorphic_rate", "polymut_rate"))
write_tsv(magstats, file="data/raw/magstats.tsv")
