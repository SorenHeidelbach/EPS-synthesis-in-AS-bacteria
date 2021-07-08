knit_with_date <- function(input, ...) {
  dir.create(paste0(dirname(input), "/html_output"), showWarnings = FALSE)
  rmarkdown::render(
    input,
    output_file = paste0(
      dirname(input),
      "/html_output/",
      basename(xfun::sans_ext(input)), ' ', format(Sys.time(), "%d.%m.%y"), '.',
      "html"
    ),
    envir = globalenv()
  )
}
# Use this function by copying the below chunk:
#```{r , eval=FALSE, include=FALSE}
#source("utilities.R")
#library(rstudioapi)
#knit_with_date(rstudioapi::getActiveDocumentContext()$path)
#```