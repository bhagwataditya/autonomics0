library(autonomics)
library(testthat)
library(magrittr)
library(autonomics.import)

# Define context
context("Test analyze_eset on Billing2016")

# Load data
result_dir <- tempdir() %>% stringi::stri_replace_all('/', fixed = '\\') %>%
              paste0('/test_analyze_contrasts')
pset <-  autonomics.data::billing2016

# Run tests
test_that("analyze_contrasts: default args", {
   result_dir %T>% message() %>% dir.create(showWarnings = FALSE)
   pset2 <- pset %>% analyze_eset(result_dir = result_dir, universe = NULL)
   assertive.types::assert_is_inherited_from(pset2, 'SummarizedExperiment')
   unlink(result_dir, recursive = TRUE)
})

