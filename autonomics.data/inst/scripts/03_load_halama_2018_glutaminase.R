# Load functions
require(magrittr)
ggplot2::theme_set(ggplot2::theme_bw())

# Load data
glutaminase <- 'extdata/glutaminase/glutaminase.xlsx'    %>% 
                system.file(package = 'autonomics.data') %>% 
                autonomics.import::load_metabolon()

# PCA: max var(samples)
glutaminase %>% autonomics::plot_pca_samples()
   # X1: time      (Right -> Left)
   # X2: treatment (Bottom -> Top)
   # X1 and X2 are correlated: effect of time is in fact effect of treatment becoming visible

# PLS: max cov(sample, subgroup)
glutaminase %>% autonomics::plot_pls_samples()
   # X1: time (Right -> Left)
   # X2: treatment (Bottom -> Top)

# LDA: max (within subgroups) / (between subgroups)
glutaminase %>% autonomics::plot_lda_samples()
   # X1: time (Left -> Right)
   # X2: treatment (Top -> Bottom)

# Compare: 
glutaminase %>% autonomics::plot_projected_samples(
   method    = c('pca', 'lda', 'sma', 'pls'),
   facet_var = c('pca', 'lda', 'sma', 'pls'),
   nrow = 2)
