#' Plot venns showing top features
#' @param object  exprs object
#' @param topdef top definition
#' @param file file to which print to 
#' @return the vennlist used for venn digram plotting
#' @export
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::exiqon
#'    object %>% autonomics.find::plot_top_venns()
#' }
plot_top_venns <- function(
   object,
   topdef = autonomics.find::default_topdef(object), 
   file = NULL
){
   contrast_names <- object %>% autonomics.find::infer_contrast_names()
   if (length(contrast_names) > 4){
      autonomics.support::cmessage('\t\tTaking max first 4 contrasts')
      contrast_names %<>% magrittr::extract(1:4)
   }
   vennlist <- mapply(autonomics.find::get_top_features,
                      contrast_name = contrast_names,
                      MoreArgs = list(object = object, topdef = topdef),
                      SIMPLIFY = FALSE)
   vennlist %>% autonomics.plot::plot_venn(
                   filename = file, 
                   title = topdef)
   return(vennlist)
}
