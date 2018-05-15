#' Create contrast matrix
#' @param contrast_defs vector of contrast defs
#' @param design        design matrix
#' @return contrast matrix
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    autonomics.data::billing2016 %>% 
#'       autonomics.find::create_design_matrix() %>% 
#'       autonomics.find::create_contrast_matrix()
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>% 
#'       autonomics.find::create_design_matrix() %>% 
#'       autonomics.find::create_contrast_matrix(
#'          c(EM0.8_0 = 'EM0.8 - EM0.0'))
#' }
#' @importFrom magrittr  %>%
#' @export
create_contrast_matrix <- function(design, contrast_defs = colnames(design)){
   
   # Create contrast matrix
   if (!assertive.properties::has_names(contrast_defs)){
      names(contrast_defs) <- make.names(contrast_defs)
   }
   if (length(contrast_defs) == 0)  return(NULL)
   assertive.properties::assert_has_no_duplicates(names(contrast_defs))
   limma::makeContrasts(contrasts = contrast_defs, levels = design) %>%
      magrittr::set_colnames(names(contrast_defs))
}