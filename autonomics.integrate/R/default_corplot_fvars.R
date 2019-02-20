#' Default correlation plot fvars
#' @param object exprs object
#' @importFrom magrittr  %>%
#' @export
default_corplot_fvars <- function(object){
   autonomics.plot::default_fvars(object) %>%
   setdiff('feature_id')
}
