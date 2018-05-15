#' Mean center each sample
#' @param object eSet
#' @return mean centered eSet
#' @importFrom magrittr %<>%
#' @export
mean_center <- function(object){
   col_means <- colMeans(autonomics.import::exprs(object), na.rm = TRUE)
   autonomics.import::exprs(object) %<>% sweep(2, col_means)
   object
}

#' Median center each sample
#' @param object eSet
#' @return median centered eSet
#' @importFrom magrittr %<>%
#' @export
median_center <- function(object){
   col_medians <- matrixStats::colMedians(autonomics.import::exprs(object), na.rm = TRUE)
   autonomics.import::exprs(object) %<>% sweep(2, col_medians)
   object
}

# Estimate the mode of a vector
# http://stackoverflow.com/a/18763500/2894278
#' @importFrom stats density
estimate_mode <- function(x) {
   d <- density(x, from = min(x, na.rm = TRUE), to = max(x, na.rm = TRUE), na.rm = TRUE)
   d$x[which.max(d$y)]
}

#' Mode center eset
#' @param object eSet
#' @param verbose whether to report progress message (logical)
#' @return mode centered eSet
#' @importFrom magrittr  %>%  %<>%
#' @export
mode_center <- function(object, verbose = TRUE){
   if (verbose) autonomics.support::cmessage('\t\tMode center sample distributions')
   col_modes <- autonomics.import::exprs(object) %>% apply(2, estimate_mode)
   autonomics.import::exprs(object) %<>% sweep(2, col_modes)
   object
}

#' Feature center eset
#' @param object eSet
#' @param feature_center Vector of numerical indexes or \code{feature_ids} (row names) to be used for centering
#' @param verbose whether to report progress message (logical)
#' @return feature centered eSet
#' @importFrom magrittr  %>%  %<>%
#' @export
feature_center <- function(object, feature_center, verbose = TRUE){
   if (verbose) autonomics.support::cmessage('\t\tFeature center sample distributions')
   col_centers <- object %>%
                  autonomics.import::exprs() %>%
                  magrittr::extract(feature_center, , drop = FALSE) %>%
                  matrixStats::colMedians(na.rm = TRUE)
   autonomics.import::exprs(object) %<>% sweep(2, col_centers)
   return(object)
}
