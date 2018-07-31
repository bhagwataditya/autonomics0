#' Complete top definition
#' @param topdef top definition
#' @param contrast_name  contrast name (string)
#' @param direction      direction
#' @return completed top definition
#' @examples
#' contrast <- c(P3S = 'P3.S - P3.A')
#' autonomics.find::complete_topdef(
#'    topdef = "abs(quantile) < 0.05", 
#'    contrast_name  = names(contrast), 
#'    direction = 'neg')
#' @importFrom magrittr %<>%
#' @export
complete_topdef <- function(topdef, contrast_name, direction = c('both', 'neg', 'both')){
   assertive.types::assert_is_a_string(topdef)
   assertive.types::assert_is_a_string(contrast_name)
   topdef %<>% gsub('(coef|rank|p|fdr|bonf|quantile)', '\\1.xxx', .)
   topdef %<>% gsub('xxx', contrast_name, .)
   topdef %<>% sprintf('(%s) & (coef.%s %s 0)', .,
                               contrast_name,
                               switch(direction, neg = '<', pos = '>', both = '!='))
   topdef
}
