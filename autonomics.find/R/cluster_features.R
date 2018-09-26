#' Cluster features
#' @param object SummarizedExperiment
#' @examples 
#' require(magrittr)
#' 
#' if (require(autonomics.data)){
#' 
#'    # Stem Cell Comparison (LCMSPROT)
#'    #--------------------------------
#'    object <- autonomics.data::stemcomp.proteinratios
#'    
#'    # Stem Cell differentiation (LCMSPROT)
#'    #-------------------------------------
#'    object <- autonomics.data::stemdiff.proteinratios
#' 
#'    # Glutaminase inhibitor (METABOLON)
#'    #----------------------------------
#'    object <- autonomics.data::glutaminase
#'    
#' }
#' @importFrom magrittr
#' @export
cluster_features <- function(object){
   
   design <- object %>% autonomics.import::create_design_matrix(intercept = TRUE)
   contrastdefs <- colnames(design)[-1] %>% magrittr::set_names(., .)
   object %<>% autonomics.find::add_limma(contrastdefs = contrastdefs, design = design)
   #object %<>% autonomics.import::filter_features(F.p < 0.05, verbose = TRUE)
   effectmat <- object %>% autonomics.import::limma()  %>% magrittr::extract(, , 'effect')
   tmat      <- object %>% autonomics.import::limma()  %>% magrittr::extract(, , 't')
   pmat      <- object %>% autonomics.import::limma()  %>% magrittr::extract(, , 'p')
      
   idx <- tmat %>% is.na() %>% matrixStats::rowAnys()
   tmat[idx, ]
   effectmat[idx, ]
   pmat[idx, ]
   
   cormat <- tmar %>% t() %>% propagate::bigcor()
   cormat
   
   # How to deal with NAs
   
}