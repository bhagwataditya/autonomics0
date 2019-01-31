#' Create design matrix for statistical analysis
#' @param object      SummarizedExperiment
#' @param intercept   logical(1): whether to include an intercept in the design
#' @param confounders confounder svars (character)
#' @return design matrix
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    object %>% autonomics.import::create_design_matrix()
#'    object %>% autonomics.import::create_design_matrix(intercept = TRUE)
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>%
#'    autonomics.import::create_design_matrix()      %>%
#'    magrittr::extract(1:3,)
#' }
#' @importFrom magrittr  %<>%
#' @export
create_design_matrix <- function(
   object,
   intercept   = length(unique(autonomics.import::sdata(object)$subgroup)) == 1,
   confounders = character(0)
){

   # Assert
   if (assertive.types::is_inherited_from(object, 'eSet')){ # Can also be EList
      autonomics.import::assert_is_valid_eset(object)
   }
   assertive.sets::assert_is_subset('subgroup', autonomics.import::svars(object))

   # Ensure that subgroup vector is a factor to preserve order of levels
   sdata1 <- autonomics.import::sdata(object)
   if(is.character(sdata1$subgroup)){
      sdata1$subgroup %<>% autonomics.support::factorify()
   }

   # Create formula
   formula <- if (intercept) '~ 1' else '~ 0'
   if (length(unique(autonomics.import::sdata(object)$subgroup)) > 1) formula %<>% paste0(' + subgroup')
   # formula <- if (length(unique(sdata1$subgroup)) < 2){   '~ 1'
   #            } else if (intercept){                      '~ 1 + subgroup'
   #            } else {                                    '~ 0 + subgroup'
   #            }
   if (length(confounders)>0){
      formula %<>% sprintf('%s + %s', ., paste0(confounders, collapse = ' + '))
   }
   formula %<>% stats::as.formula()

   # Create design matrix
   myDesign <- stats::model.matrix(formula,  data = sdata1)

   # Rename coefficients
   # ~ 1 + subgroup
   if (intercept){
      subgroup1 <- unique(autonomics.import::sdata(object)$subgroup)[1] %>% as.character()
      colnames(myDesign) %<>% gsub('(Intercept)', subgroup1, ., fixed = TRUE)
      colnames(myDesign) %<>% stringi::stri_replace_first_regex('subgroup(.+)', paste0('$1_', subgroup1))
   # ~ 0 + subgroup
   } else {
      colnames(myDesign) %<>% gsub('subgroup',    '',        ., fixed = TRUE)
   }

   # Rename confounders
   if (length(confounders) > 0){
      numeric_confounders <- autonomics.import::sdata(object)[, confounders, drop = FALSE]   %>%
                             vapply(is.numeric, logical(1)) %>%  which() %>% names()
      factor_confounders  <- confounders %>% setdiff(numeric_confounders)
   }

   # Validify names
   colnames(myDesign) %<>% gsub(':', '..', ., fixed = TRUE)
   colnames(myDesign) %<>% make.names()

   # Return
   return(myDesign)
}
