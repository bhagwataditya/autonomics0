require(magrittr)

# Stem cell comparison
#=====================

# PROTEINGROUPS
stemcomp.proteinratios <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>% 
                           system.file(package = 'autonomics.data') %>% 
                           autonomics.import::load_proteingroups(infer_design_from_sampleids = TRUE) %>% 
                           autonomics.preprocess::invert_ratios(
                              invert_subgroups = c('E_EM', 'E_BM', 'EM_BM'), 
                              subgroup_frac = '_')
stemcomp.proteinratios %>% autonomics.import::sdata() %>% str()
save(stemcomp.proteinratios, file = 'data/stemcomp.proteinratios.RData', compress = 'xz')

# SOMA
stemcomp.soma <- 'extdata/stemcomp/soma/stemcomp.adat'      %>% 
                  system.file(package = 'autonomics.data')  %>%
                  autonomics.import::load_soma()
stemcomp.soma %>% autonomics.import::sdata() %>% str()
save(stemcomp.soma, file = 'data/stemcomp.soma.RData', compress = 'xz')


# Stem cell differentiation
#==========================

# PROTEINGROUPS
stemdiff.proteinratios <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>% 
                           system.file(package = 'autonomics.data')     %>% 
                           autonomics.import::load_proteingroups(infer_design_from_sampleids = TRUE)           %>% 
                           autonomics.import::filter_samples(subgroup %>% stringi::stri_detect_fixed('_STD'))  %>% 
                              # ony interested in ./STD(L) ratios
                           autonomics.import::filter_samples(subgroup %>% stringi::stri_detect_fixed('BLANK_') %>% magrittr::not()) %>% 
                              # not interested in BLANK/STD ratios
                           autonomics.preprocess::invnorm()

autonomics.import::sdata(stemdiff.proteinratios) %>% str()
stemdiff.proteinratios$subgroup %<>% factor(c('EM00_STD', 'EM01_STD', 'EM02_STD', 'EM05_STD','EM15_STD', 'EM30_STD', 'BM00_STD'))
stemdiff.proteinratios %<>% autonomics.import::arrange_samples(subgroup)

save(stemdiff.proteinratios, file = 'data/stemdiff.proteinratios.RData', compress = 'xz')
stemdiff.proteinratios %>% autonomics.explore::plot_pca_samples()
stemdiff.proteinratios %>% autonomics.plot::plot_sample_distributions()
stemdiff.proteinratios %>% autonomics.plot::default_color_values(color_var = 'subgroup')

# Glutaminase (METABOLON)
#========================

glutaminase <- 'extdata/glutaminase/glutaminase.xlsx'     %>% 
                          system.file(package = 'autonomics.data')  %>% 
                          autonomics.import::load_metabolon()
save(glutaminase, file = 'data/glutaminase.RData', compress = 'xz')



#' Save eset with preprocessed Billing 2016 data
#' 
#' Billing 2016 is a stem cell comparison experiment, 
#' in which three different types of stem cells (E, EM, BM)
#' are compared to each other.
#' 
#' @noRd
#' @importFrom magrittr %<>%
save_billing2016 <- function(){
   autonomics.import::create_maxquant_sample_file( 'inst/extdata/billing2016/proteinGroups.txt')
   billing2016 <- autonomics.import::load_proteingroups('inst/extdata/billing2016/proteinGroups.txt') %>% 
      autonomics.preprocess::invert_ratios(
         invert_subgroups = c('E_EM', 'E_BM', 'EM_BM'), 
         subgroup_frac    = '_')
   billing2016 %<>% autonomics.import::annotate_proteingroups()
   autonomics.import::fdata(billing2016)$`Fasta headers` <- NULL
   billing2016 %>%  autonomics.import::fdata() %>% head()
   billing2016 %<>% autonomics.find::add_limma_to_fdata()
   save(billing2016, file = 'data/billing2016.RData', compress = 'xz')
}

#save_billing2016()
