#' Eliminate confounders
#' @param object      exprs object
#' @param confounders confounding svars (character)
#' @return object with confounder effects removed
#' @importFrom magrittr   %>%   %<>%
#' @export
eliminate_confounders <- function(
   object, 
   confounders = character()
){
   for (cur_confounder in confounders){
      values  <- autonomics.import::sdata(object) %>% magrittr::extract2(cur_confounder)
      design  <- object %>%  autonomics.find::create_design_matrix(confounders = setdiff(confounders, cur_confounder))
      if (is.numeric(values))   autonomics.import::exprs(object) %<>% limma::removeBatchEffect(covariates = values, design = design)
      else                      autonomics.import::exprs(object) %<>% limma::removeBatchEffect(batch      = values, design = design)
      object
   }
   object
}
