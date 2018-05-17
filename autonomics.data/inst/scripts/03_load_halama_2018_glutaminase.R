# Load functions
require(magrittr)
ggplot2::theme_set(ggplot2::theme_bw())

# Load data
glutaminase <- autonomics.import::load_metabolon('inst/extdata/glutaminase/glutaminase.xlsx')

# Rm blank samples
glutaminase$subgroup

# Explore
glutaminase %>% autonomics.explore::plot_projected_samples2(
   method    = c('pca', 'lda', 'sma', 'pls'),
   facet_var = c('pca', 'lda', 'sma', 'pls'),
   nrow = 2)


glutaminase %>% autonomics.explore::plot_pca_samples2()
glutaminase %>% autonomics.explore::plot_sma_samples2()
glutaminase %>% autonomics.explore::plot_pls_samples2()
glutaminase %>% autonomics.explore::plot_lda_samples2()
