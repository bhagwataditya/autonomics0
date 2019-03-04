require(magrittr)

# Stem cell comparison
#=====================

# PROTEINGROUPS
stemcomp.proteinratios <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>% 
                           system.file(package = 'autonomics.data')     %>% 
                           autonomics.import::read_proteingroups()      %>% 
                           autonomics::prepare_proteingroups(invert_subgroups = c('E_EM', 'E_BM', 'EM_BM'), 
                                                             deconvolution_fastafile = '../data/uniprot_hsa_20140515.fasta') %>% 
                           autonomics.import::set_contrastdefs(c( EM_E =  'EM_E', BM_E =  'BM_E', BM_EM = 'BM_EM')) %>% 
                           autonomics.find::add_limma()
save(stemcomp.proteinratios, file = 'data/stemcomp.proteinratios.RData', compress = 'xz')

# SOMA
stemcomp.soma <- 'extdata/stemcomp/soma/stemcomp.adat'      %>% 
                  system.file(package = 'autonomics.data')  %>%
                  autonomics.import::read_somascan()        %>%
                  autonomics::prepare_somascan()     %>% 
                  autonomics.find::add_limma(c(EM_E  = 'EM-E', BM_E  = 'BM-E', BM_EM = 'BM-EM'))
save(stemcomp.soma, file = 'data/stemcomp.soma.RData', compress = 'xz')


# Stem cell differentiation
#==========================

# PROTEINGROUPS
stemdiff.proteinratios <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>% 
                           system.file(package = 'autonomics.data')     %>% 
                           autonomics.import::read_proteingroups()      %>% 
                           autonomics::prepare_proteingroups(deconvolution_fastafile = '../data/uniprot_hsa_20140515.fasta') %>% 
                           autonomics.import::filter_samples(subgroup %>% stringi::stri_detect_fixed('_STD'))  %>% 
                              # ony interested in ./STD(L) ratios
                           autonomics.import::filter_samples(subgroup %>% stringi::stri_detect_fixed('BLANK_') %>% magrittr::not()) %>% 
                              # not interested in BLANK/STD ratios
                           autonomics.preprocess::invnorm() %>% 
                          (function(x){x$subgroup  %<>% stringi::stri_replace_first_fixed('_STD', '')
                                       x$sample_id %<>% stringi::stri_replace_first_fixed('_STD', '')
                                       autonomics.import::snames(x)    %<>% stringi::stri_replace_first_fixed('_STD', '')
                                       x$subgroup %<>% factor(c('EM00', 'EM01', 'EM02', 'EM05','EM15', 'EM30', 'BM00'))
                                       x %<>% autonomics.import::arrange_samples(subgroup)
                                       x
                           }) %>% 
                           autonomics.import::set_contrastdefs(autonomics.find::make_ref_contrasts(.)) %>% 
                           autonomics.find::add_limma()

save(stemdiff.proteinratios, file = 'data/stemdiff.proteinratios.RData', compress = 'xz')
stemdiff.proteinratios %>% autonomics.plot::plot_pca_samples()
stemdiff.proteinratios %>% autonomics.plot::plot_sample_distributions()
stemdiff.proteinratios %>% autonomics.plot::default_color_values(color_var = 'subgroup')

# Glutaminase (METABOLON)
#========================
require(magrittr)
glutaminase <- 'extdata/glutaminase/glutaminase.xlsx'     %>% 
                system.file(package = 'autonomics.data')  %>% 
                autonomics.import::read_metabolon()       %>% 
                autonomics::prepare_metabolon()
                #autonomics.import::set_contrastdefs(autonomics.find::make_ref_contrasts(.)) %>% 
                #autonomics.find::add_limma()
autonomics.import::sdata(glutaminase) %>% str()
glutaminase$subgroup %<>% factor(autonomics.support::vsprintf('%s_%s', c('UT', 'Veh', 'uM05', 'uM10'), 
                                                                       c('h10', 'h24', 'h48', 'h72')))
usethis::use_data(glutaminase, compress = 'xz', overwrite = TRUE)


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
