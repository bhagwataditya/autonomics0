#' Create dataframe for feature plots
#' @param object   SummarizedExperiment
#' @param fvars    fvars used to annotate each plot
#' @param verbose  logical
#' @return dataframe
#' @author Aditya Bhagwat
#' @examples
#' if (require(autonomics.data)){
#'   require(magrittr)
#'
#'   # STEM CELL COMPARISON (Max Quant)
#'   autonomics.data::stemcomp.proteinratios %>%
#'   create_feature_plot_df(c('feature_id', 'Gene names')) %>% head()
#'
#'   # GLUTAMINASE
#'   autonomics.data::glutaminase %>%
#'   create_feature_plot_df(c('feature_id', 'BIOCHEMICAL')) %>% head()
#'
#'
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
create_feature_plot_df <- function(object, fvars, verbose = FALSE){

   # Make sure fvars is unique (otherwise fvar.2 -> strange behaviour later on)
   fvars %<>% unique()

   # Replace -Inf with NA
   idx <- is.infinite(autonomics.import::exprs(object))
   if (sum(idx) > 0){
      if (verbose) autonomics.support::cmessage('\t\t %d -Inf -> NA', sum(idx) )
      autonomics.import::exprs(object)[idx] <- NA
   }

   # Create 'wide' plotDF
   plotDF <- data.frame(autonomics.import::exprs(object), check.names = FALSE, stringsAsFactors = FALSE)

   # Collapse fvars and add to plotDF
   fdata1 <- autonomics.import::fdata(object)[,fvars, drop = FALSE]
   numeric_fvars <- fdata1 %>% vapply(is.numeric, logical(1))
   fdata1[,!numeric_fvars] <- fdata1[,!numeric_fvars, drop = FALSE] %>% lapply(as.character)
   fdata1[, numeric_fvars] <- fdata1[, numeric_fvars, drop = FALSE] %>%
                              lapply(function(x) sprintf('%.2e', signif(x, digits=3)))
   feature_facet <- fdata1 %>%
                    Reduce(function(...) paste(..., sep = " :: "), .) %>% #https://stackoverflow.com/a/33015233/2894278
                    autonomics.support::uniquify() %>%
                    autonomics.support::factorify()
   plotDF <- cbind(feature_facet, fdata1, plotDF)

   # Gather columns
   plotDF2 <- plotDF %>% tidyr::gather_(key = 'sample', value = 'value', gather_cols = setdiff(colnames(plotDF), c('feature_facet', fvars)))
   plotDF2 %<>% magrittr::extract(!is.na(plotDF2$value), ) # https://github.com/hadley/ggplot2/issues/791

   # Add pData
   idx <- match(plotDF2$sample, autonomics.import::snames(object))
   plotDF3 <- cbind(plotDF2, autonomics.import::sdata(object)[idx, , drop = FALSE], row.names = NULL)
   return(plotDF3)
}
