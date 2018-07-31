library(magrittr)
library(testthat)

context('run_ora_on_eset: test "universe" argument')

# Set parameters
(pset <- autonomics.data::billing2016 %>% 
         magrittr::extract(1:100, )   %>% 
         autonomics.find::add_limma_to_fdata())
(contrasts <- c(BM_E = 'BM - E'))
(result_dir <- paste0(tempdir(), '/test_run_ora_on_eset'))
(topdef <- 'p < 0.05')
(organism <- 'hsa')

# Test

test_that("run_ora_on_eset doesn't break for universe == NULL",
  {
    pset %>% run_ora_on_eset(
      contrasts      = contrasts, 
      result_dir     = result_dir,
      topdef         = topdef,
      universe       = NULL
    ) %>% expect_null()
  } 
)

test_that("run_ora_on_eset works for universe == 'detectome'",
  {
    pset %>% run_ora_on_eset(
                contrasts      = contrasts, 
                result_dir     = result_dir,
                topdef         = topdef,
                universe       = 'detectome', 
                gene_set_collections = 'gobp'
             )
    
    result_files <- list.files(paste0(result_dir, '/contrasts/BM_E/ora_in_detectome'))
    result_files %>% stringi::stri_detect_fixed('gobp') %>% any() %>% expect_true()
    unlink(result_dir, recursive = TRUE)
  } 
)

test_that("run_ora_on_eset works for universe == 'genome'",
  {
    pset %>% run_ora_on_eset(
                contrasts            = contrasts, 
                result_dir           = result_dir,
                topdef               = topdef,
                universe             = 'genome', 
                gene_set_collections = 'gobp'
             )
    
    result_files <- list.files(paste0(result_dir, '/contrasts/BM_E/ora_in_genome'))
    result_files %>% stringi::stri_detect_fixed('gobp') %>% any() %>% expect_true()
    unlink(result_dir, recursive = TRUE)
  }
)

