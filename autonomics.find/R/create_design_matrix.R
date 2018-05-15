#' Create design matrix for statistical analysis
#' @param object      eSet or eList
#' @param confounders confounder svars (character)
#' @return design matrix
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::billing2016            %>%
#'    autonomics.find::create_design_matrix() %>% 
#'    magrittr::extract(1:3,)
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>% 
#'    autonomics.find::create_design_matrix()      %>% 
#'    magrittr::extract(1:3,)
#' }
#' @importFrom magrittr  %<>%
#' @export
create_design_matrix <- function(object, confounders = character(0)){

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
   formula <- if (length(unique(sdata1$subgroup)) < 2){   '~ 1'
              } else {                                    '~ 0 + subgroup'
              }
   if (length(confounders)>0){
      formula %<>% sprintf('%s + %s', ., paste0(confounders, collapse = ' + '))
   }
   formula %<>% stats::as.formula()

   # Create design matrix
   myDesign <- stats::model.matrix(formula,  data = sdata1)

   # Rename intercepts
   subgroup1 <- unique(autonomics.import::sdata(object)$subgroup)[1] %>% as.character()
   colnames(myDesign) %<>% gsub('(Intercept)', subgroup1, ., fixed = TRUE) # ~ 1
   colnames(myDesign) %<>% gsub('subgroup',    '',        ., fixed = TRUE) # ~ 0 + subgroup

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
