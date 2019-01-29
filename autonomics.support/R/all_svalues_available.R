#' Are all svalues available?
#' @param object SummarizedExperiment
#' @param svar   character(1)
#' @return logical(1)
#' @importFrom magrittr %>%
#' @export
all_svalues_available <- function(object, svar){
   svalues1 <- object %>% autonomics.import::svalues(svar)
   if (is.null(svalues1))                                                return(FALSE)
   if (autonomics.support::all_are_missing_or_empty_character(svalues1)) return(FALSE)
   TRUE
}
