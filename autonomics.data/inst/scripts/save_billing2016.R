
#' Save eset with preprocessed Billing 2016 data
#' 
#' Billing 2016 is a stem cell comparison experiment, 
#' in which three different types of stem cells (E, EM, BM)
#' are compared to each other.
#' 
#' @noRd
#' @importFrom magrittr %<>%
save_billing2016 <- function(){
   autonomics.import::create_maxquant_design_file( 'inst/extdata/billing2016/proteinGroups.txt')
   billing2016 <- autonomics.import::load_proteingroups('inst/extdata/billing2016/proteinGroups.txt') %>% 
                  autonomics.preprocess::invert_ratios(
                     invert_subgroups = c('E_EM', 'E_BM', 'EM_BM'), 
                     channel_frac     = '/', 
                     subgroup_frac    = '_')
   billing2016 %<>% autonomics.import::annotate_proteingroups()
   autonomics.import::fdata(billing2016)$`Fasta headers` <- NULL
   billing2016 %>%  autonomics.import::fdata() %>% head()
   billing2016 %<>% autonomics.find::add_limma_to_fdata()
   save(billing2016, file = 'data/billing2016.RData', compress = 'xz')
}

require(magrittr)
save_billing2016()
