#' mean_and_sd
#'
#' Updates data with feature/sample (row/column)-wise mean,
#' standard deviation, as well as a completion fraction (number of \code{\link{NA}}/measurements.
#' Under the hood \code{\link[matrixStats]{rowMeans2}} and \code{\link[matrixStats]{rowSds}} or
#' the \code{col} equivalent are used, all with \code{na.rm = TRUE}.
#' 
#' For subgroups without replicates, \code{NA} is returned.
#'
#' @param object SummarizedExperiment, eSet, or EList
#' @param MARGIN a vector giving the subscripts which the calculations will
#' be done for. 1 indicates rows/features, 2 indicates columns/samples,
#' c(1, 2) indicates rows/features and columns/samples.
#' @param by_subgroup Whether or not to calculate summary statistics per subgroup
#' (using \code{autonomics.import::sdata(object)$subgroup}) - used for
#' \code{MARGIN == 1} only;
#' @return updated object
#' @importFrom magrittr %<>% %>%
#' @export
mean_and_sd <- function(object, MARGIN = c(1, 2), by_subgroup = TRUE)
{
   object %>%
      autonomics.import::assert_is_valid_eset()
   MARGIN %<>%
      as.character() %>%
      match.arg(
         choices    = c(1,2),
         several.ok = TRUE)
   by_subgroup %>%
      assertive.types::assert_is_a_bool()
   if(MARGIN == 2 && by_subgroup)
   {
      warning('\'by_subgroup\' is disregarded when MARGIN == 2 (column-wise operation)')
   }

   if(1 %in% MARGIN)
   {
      if(by_subgroup)
      {
         subgroups <- object %>%
            autonomics.import::sdata() %>%
            magrittr::extract2('subgroup')

         autonomics.import::fdata(object) %<>% dplyr::bind_cols(
            subgroups %>%
               unique() %>%
               lapply(
                  function(x)
                  {
                     sub_object_exprs <- object %>%
                        autonomics.import::exprs() %>%
                        magrittr::extract(, subgroups == x)
                     if(length(dim(sub_object_exprs)) < 2) 
                     {
                       data.frame(
                         sub_mean = rep(x, times = ncol(sub_object_exprs)),
                         sub_sd   = rep(x, times = ncol(sub_object_exprs)),
                         sub_cmpl = rep(x, times = ncol(sub_object_exprs))) %>%
                         magrittr::set_colnames(
                           paste0(x, c('__mean', '__sd', '__cmpl')))
                     } else {
                       data.frame(
                         sub_mean = sub_object_exprs %>%
                           matrixStats::rowMeans2(na.rm = TRUE),
                         sub_sd   = sub_object_exprs %>%
                           matrixStats::rowSds(na.rm = TRUE),
                         sub_cmpl = sub_object_exprs %>%
                           is.na() %>%
                           magrittr::not() %>%
                           matrixStats::rowSums2() %>%
                           magrittr::divide_by(
                             sub_object_exprs %>%
                               ncol())) %>%
                         magrittr::set_colnames(
                           paste0(x, c('__mean', '__sd', '__cmpl')))
                     }
                  }
               ) %>%
               dplyr::bind_cols())
      } else {
         autonomics.import::fdata(object) %<>%
            dplyr::mutate(
               mean = object %>%
                  autonomics.import::exprs() %>%
                  matrixStats::rowMeans2(na.rm = TRUE),
               sd = object %>%
                  autonomics.import::exprs() %>%
                  matrixStats::rowSds(na.rm = TRUE),
               cmpl = object %>%
                  autonomics.import::exprs() %>%
                  is.na() %>%
                  magrittr::not() %>%
                  matrixStats::rowSums2() %>%
                  magrittr::divide_by(
                     object %>%
                        autonomics.import::exprs() %>%
                        ncol()))
      }
   }
   if(2 %in% MARGIN)
   {
      autonomics.import::sdata(object) %<>%
         dplyr::mutate(
            mean = object %>%
               autonomics.import::exprs() %>%
               matrixStats::colMeans2(na.rm = TRUE),
            sd = object %>%
               autonomics.import::exprs() %>%
               matrixStats::colSds(na.rm = TRUE),
            cmpl = object %>%
               autonomics.import::exprs() %>%
               is.na() %>%
               magrittr::not() %>%
               matrixStats::colSums2() %>%
               magrittr::divide_by(
                  object %>%
                     autonomics.import::exprs() %>%
                     nrow()))
   }

   object %>%
      invisible()
}
