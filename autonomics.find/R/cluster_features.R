#' @examples 
#' require(magrittr)
#' object <- autonomics.data::stemdiff.proteinratios
#' 
cluster_features <- function(object){
   
   design <- object %>% autonomics.import::create_design_matrix(intercept = TRUE)
   contrastdefs <- colnames(design)[-1] %>% magrittr::set_names(., .)
   object %<>% autonomics.find::add_limma(contrastdefs = contrastdefs, design = design)
   object %<>% autonomics.import::filter_features(F.p < 0.05)
   cormat <- object %>% autonomics.import::limma()  %>% 
                        magrittr::extract(, , 't')  %>% 
                        cor()
   
   # How to deal with NAs
   
}