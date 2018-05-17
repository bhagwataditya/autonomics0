# Load functions
require(magrittr)

# Load data
glutaminase <- 'extdata/glutaminase/glutaminase.xlsx'   %>% 
                system.file(package = 'autonomics.data') %>% 
                autonomics.import::load_metabolon()

# Rm blank samples
glutaminase$subgroup

# Explore
ggplot2::theme_set(ggplot2::theme_bw())
glutaminase %>% autonomics.explore::plot_projected_samples2(
   method    = c('pca', 'lda', 'sma', 'pls'),
   facet_var = c('pca', 'lda', 'sma', 'pls'),
   nrow = 2)