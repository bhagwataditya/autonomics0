#' Plot contrast venns
#' @param object         SummarizedExperiment
#' @param contrast_names character vector
#' @param topdef         top definition
#' @param title          plot title
#' @param euler          logical(1): whether to plot euler rather than venn diagram
#' @param file           file to which print to 
#' @param ...            autonomics.plot::plot_venn(...)
#' @return the vennlist used for venn digram plotting
#' @export
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::exiqon
#'    object %>% autonomics.find::plot_contrast_venns()
#' }
plot_contrast_venns <- function(
   object,
   contrast_names = object %>% autonomics.find::infer_contrast_names(),
   topdef         = autonomics.find::default_topdef(object), 
   title          = topdef,
   euler          = FALSE,
   file           = NULL, 
   ...
){
   vennlist <- mapply(autonomics.find::get_contrast_features,
                      contrast_name = contrast_names,
                      MoreArgs = list(object = object, topdef = topdef),
                      SIMPLIFY = FALSE)
   vennlist %>% autonomics.plot::plot_venn(
                   filename = file, 
                   title    = title, 
                   euler    = euler,
                   ...)
   return(vennlist)
}
