#' Complete top definition
#' @param top_definition top definition
#' @param contrast_name  contrast name (string)
#' @param direction      direction
#' @return completed top definition
#' @examples
#' contrast <- c(P3S = 'P3.S - P3.A')
#' autonomics.find::complete_top_definition(
#'    top_definition = "abs(quantile) < 0.05", 
#'    contrast_name  = names(contrast), 
#'    direction = 'neg')
#' @importFrom magrittr %<>%
#' @export
complete_top_definition <- function(top_definition, contrast_name, direction = c('both', 'neg', 'both')){
   assertive.types::assert_is_a_string(top_definition)
   assertive.types::assert_is_a_string(contrast_name)
   top_definition %<>% gsub('(coef|rank|p|fdr|bonf|quantile)', '\\1.xxx', .)
   top_definition %<>% gsub('xxx', contrast_name, .)
   top_definition %<>% sprintf('(%s) & (coef.%s %s 0)', .,
                               contrast_name,
                               switch(direction, neg = '<', pos = '>', both = '!='))
   top_definition
}
