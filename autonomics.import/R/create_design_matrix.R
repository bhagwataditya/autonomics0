#' Create design matrix for statistical analysis
#' @param object      SummarizedExperiment
#' @param intercept   TRUE or FALSE: include an intercept in the design?
#' @param confounders confounder svars (character)
#' @return design matrix
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'
#'    # STEM CELL COMPARISON
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% create_design_matrix()
#'
#'    # GLUTAMINASE
#'    object <- autonomics.data::glutaminase
#'    object %>% create_design_matrix() %>% extract(1:10, 1:10)
#' }
#' @importFrom magrittr  %<>%
#' @export
create_design_matrix <- function(
   object,
   intercept   = length(unique(sdata(object)$subgroup)) == 1,
   confounders = character(0)
){

   # Assert
   assert_is_valid_object(object)
   assertive.sets::assert_is_subset('subgroup', svars(object))

   # Ensure that subgroup vector is a factor to preserve order of levels
   sdata1 <- sdata(object)
   if(is.character(sdata1$subgroup)){
      sdata1$subgroup %<>% autonomics.support::factorify()
   }

   # Create formula
   formula <- if (intercept) '~ 1' else '~ 0'
   if (length(unique(sdata(object)$subgroup)) > 1) formula %<>% paste0(' + subgroup')
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
      subgroup1 <- unique(sdata(object)$subgroup)[1] %>% as.character()
      colnames(myDesign) %<>% gsub('(Intercept)', subgroup1, ., fixed = TRUE)
      colnames(myDesign) %<>% stringi::stri_replace_first_regex('subgroup(.+)', paste0('$1_', subgroup1))
   # ~ 0 + subgroup
   } else {
      colnames(myDesign) %<>% gsub('subgroup',    '',        ., fixed = TRUE)
   }

   # Rename confounders
   if (length(confounders) > 0){
      numeric_confounders <- sdata(object)[, confounders, drop = FALSE]   %>%
                             vapply(is.numeric, logical(1)) %>%  which() %>% names()
      factor_confounders  <- confounders %>% setdiff(numeric_confounders)
   }

   # Validify names
   colnames(myDesign) %<>% gsub(':', '..', ., fixed = TRUE)
   colnames(myDesign) %<>% make.names()

   # Return
   return(myDesign)
}
