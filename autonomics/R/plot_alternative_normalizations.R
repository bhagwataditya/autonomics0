#' Plot preprocessing alternatives
#' @param object           eset
#' @param result_dir         result dir
#' @param color_var          color var
#' @param invert_subgroups   vector of subgroup levels which require inversion
#' @importFrom magrittr  %>%  %<>%
#' @export
plot_preprocessing_alternatives <- function(
   object,
   result_dir,
   color_var = autonomics.plot::default_color_var(object),
   invert_subgroups = NULL
){

   # Perform alternative normalizations
   eset_original <- object
   autonomics.import::exprs(eset_original) %<>% (function(x) 2^x)

   # Original
   esets <- list(original = eset_original)

   # Log2
   esets %<>% c(list(object))
   names(esets)[length(esets)] <- 'log2'

   # Invert subgroups (labeled proteomics)
   if (!is.null(invert_subgroups)){
      esets %<>% c(list(autonomics.preprocess::invert_ratios(object, invert_subgroups)))
      names(esets)[length(esets)] <- sprintf('%s.flipratios', names(esets)[length(esets)-1])
   }
   cur_name <- names(esets)[length(esets)]

   # Center
   if (autonomics.import::prepro(object)$quantity == 'ratio'){
      esets %<>% c(list(autonomics.preprocess::mode_center(object)))
      names(esets)[length(esets)] <- sprintf('%s.center', cur_name)
   }

   # Normalize within samples
   esets %<>% c(list(object %>% autonomics.preprocess::invnorm()))
   names(esets)[length(esets)] <- sprintf('%s.inv', cur_name)

   # Normalize between samples
   esets %<>% c(list(object %>% autonomics.preprocess::quantnorm()))
   names(esets)[length(esets)] <- sprintf('%s.qnorm', cur_name)

   # Quantnorm within subgroups
   if (autonomics.import::prepro(object)$quantity == 'ratio'){
       esets %<>% c(list(object %>% autonomics.preprocess::quantnorm_within_subgroups()))
      names(esets)[length(esets)] <- sprintf('%s.qgroups', cur_name)
   }

   # Z score
   esets %<>% c(list(object %>% autonomics.preprocess::z_score_samples()))
   names(esets)[length(esets)] <- sprintf('%s.zscore', cur_name)

   # Invnorm and qnorm
   esets %<>% c(list(object %>% autonomics.preprocess::invnorm() %>% autonomics.preprocess::quantnorm()))
   names(esets)[length(esets)] <- sprintf('%s.inv.qnorm', cur_name)

   # Invnorm and qnorm_within_subgroups]
   esets %<>% c(list(object %>% autonomics.preprocess::invnorm() %>% autonomics.preprocess::quantnorm_within_subgroups()))
   names(esets)[length(esets)] <- sprintf('%s.inv.qgroups', cur_name)

   # Center, invnorm, and qnorm
   esets %<>% c(list(object %>% autonomics.preprocess::mode_center() %>% autonomics.preprocess::invnorm() %>% autonomics.preprocess::quantnorm()))
   names(esets)[length(esets)] <- sprintf('%s.center.inv.qnorm', cur_name)

   # Center, invnorm, and qnorm within subgroups
   esets %<>% c(list(object %>% autonomics.preprocess::mode_center() %>% autonomics.preprocess::invnorm() %>% autonomics.preprocess::quantnorm_within_subgroups()))
   names(esets)[length(esets)] <- sprintf('%s.center.inv.qgroups', cur_name)

   # Plot
   pcas                 <- autonomics.explore::plot_pca_samples %>%
                           mapply(object = esets, title = names(esets), MoreArgs = list(color_var = color_var), SIMPLIFY = FALSE)
   sample_distributions <- autonomics.plot::plot_sample_distributions %>%
                           mapply(object = esets, title = names(esets), MoreArgs = list(color_var = color_var), SIMPLIFY = FALSE)
   plotlist <- c(sample_distributions, pcas)
   nplots <- length(esets)  + (4 - length(esets) %% 4)

   grDevices::pdf(paste0(result_dir, '/preprocessing_distributions.pdf'), width = 21, height = 14)
   autonomics.plot::multiplot(plotlist = sample_distributions, layout = matrix(1:nplots, byrow = TRUE, ncol = 4))
   grDevices::dev.off()

   grDevices::pdf(paste0(result_dir, '/preprocessing_pca_plots.pdf'), width = 21, height = 14)
   autonomics.plot::multiplot(plotlist = pcas, layout = matrix(1:nplots, byrow = TRUE, ncol = 4))
   grDevices::dev.off()
}
