library(autonomics)
library(testthat)
library(magrittr)
library(autonomics.import)

# Define context
context("Test analyze_eset on Atkin hypoglycemia")

# Load data
result_dir <- tempdir() %>% stringi::stri_replace_all('/', fixed = '\\') %>%
              paste0('/test_analyze_eset')
pset <-  atkin.2014::soma
pset$subgroup
pset$block
contrast_defs <- atkin.2014::contrasts[1]
