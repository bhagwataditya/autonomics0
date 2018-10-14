#' Cluster features
#' @param object SummarizedExperiment
#' @examples 
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::exiqon
#'    object %>% autonomics.plot::plot_overlayed_sample_distributions()
#' }
#' @importFrom magrittr %>% 
#' @export
cluster_features <- function(object){
   
   design <- object %>% autonomics.import::create_design_matrix(intercept = TRUE)
   contrastdefs <- colnames(design)[-1] %>% magrittr::set_names(., .)
   object %<>% autonomics.find::add_limma(contrastdefs = contrastdefs, design = design)
   object %<>% autonomics.import::filter_features_('F.p < 0.05')
   
   #object %<>% autonomics.import::filter_features(F.p < 0.05, verbose = TRUE)
   effectmat <- object %>% autonomics.import::limma()  %>% magrittr::extract(, , 'effect')
   tmat      <- object %>% autonomics.import::limma()  %>% magrittr::extract(, , 't')
   pmat      <- object %>% autonomics.import::limma()  %>% magrittr::extract(, , 'p')
      
   # How to deal with NAs - no NAS left in exiqon data
   #idx <- tmat %>% is.na() %>% matrixStats::rowAnys()
   #tmat[idx, ]
   #effectmat[idx, ]
   #pmat[idx, ]
   cormat <- tmat %>% t() %>% propagate::bigcor()
   cormat %<>% magrittr::extract(1:nrow(.), 1:ncol(.))
   
   # Cluster
   apres <- apcluster::apcluster(s = cormat, details = TRUE, q=0)
   
   # Plot exemplars
   object[apres@exemplars, ] %>% autonomics.plot::plot_features(x = 'time', color_var = 'subgroup', line = FALSE, group_var = 'condition')
   
   # Plot clusters
   object[apres@exemplars, ] %>% autonomics.plot::plot_features(x = 'time', color_var = 'condition', line = TRUE, group_var = 'condition')   
   i <- 1
   for (i in 1:length(apres)){
      object[apres@clusters[[i]], ] %>% autonomics.plot::plot_features(
                                           color_var = 'condition', 
                                           line      = FALSE, 
                                           group_var = 'condition', 
                                           file = sprintf('~/Projects/datasets/subramanian.dio/results/cluster_%d.pdf', i), 
                                           width = 13, 
                                           height = 9)
   }
}