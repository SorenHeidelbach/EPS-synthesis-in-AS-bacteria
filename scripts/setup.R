library(here)
library(fs)

dir_create(here("data", "raw"))        # raw data
dir_create(here("data", "processed"))  # intermediate data
dir_create(here(c("scripts",           # storing all .R, .Rmd, .py files
                  "figures",           # output figures
                  "output",            # output files
                  "previous",          # previous versions keep for record
                  "paper",             # current versions of the document
                  "notes",             # random notes
                  "documents"          # supporting documents
)
)
)


# have this code chunk in the start of every file to refer to the folders
library(here)
d <- function(x, raw = TRUE){paste0(here(),"/data/",ifelse(raw, "raw/", "processed/"), x)}
o <- function(x){paste(here(), x, sep = "/output/")}
f <- function(x){paste(here(), x, sep = "/figures/")}
s <- function(x){paste(here(), x, sep = "/scripts/")}
