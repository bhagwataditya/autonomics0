#' @rdname missify_inconsistent_zeroes
#' @export
zero_to_na_if_not_all_zero <- function(
   no_zero = c('fail', 'warn_passthrough', 'passthrough')[1],
   ...
){
   missify_inconsistent_zeroes(...)
}

#' Missify inconsistent zeroes
#'
#' Convert inconsistent zeroes into NAs. \cr
#' Inconsistent zeroes are those not replicated across samples of the same subgroup.
#' Such zeroes are likely caused by identification/quantification stochasticity.
#' Their proper representation is NA rather than zero.
#' Failure to make this distinction diminishes statistical power.
#'
#' Inconsistent zeroes are common in LCMS proteomics and metabolomics, due to the
#' stochastic nature of the identification process.
#'
#' @param object SummarizedExperiment
#' @param no_zero \code{\link{character}} defining how to react of there are no \code{0}s (zeros) present in the data set
#' @return updated object
#' @importFrom magrittr %<>% %>%
#' @export
missify_inconsistent_zeroes <- function(
  object,
  verbose = FALSE
){
# Check prerequisites -----------------------------------------------------
  autonomics.import::assert_is_valid_eset(object)
  if (!autonomics.import::has_complete_subgroup_values(object)){
     autonomics.support::cmessage('Return unmodified - object lacks complete subgroup values')
     return(object)
  }
  if(object %>% autonomics.import::exprs() %>% magrittr::equals(0) %>% any() %>% magrittr::not()){
    autonomics.support::cmessage("'object' contains no '0's (zeros).")
    return(object)
  }

# Processing --------------------------------------------------------------
  # Split by subgroup
  data_by_subgroup <- object %>%
    autonomics.import::subgroup_levels() %>%
    lapply(
      function(sg)
      {
        object %>%
          autonomics.import::filter_samples(subgroup == sg)
      }) %>%
    set_names(autonomics.import::subgroup_levels(object)) %>%
    ## Isolate expressions
    lapply(autonomics.import::exprs)

  # Count zeros per row
  zeros_per_row_by_subgroup <- data_by_subgroup %>%
    ## Check for zero-ness
    lapply(magrittr::equals, 0) %>%
    ## Check for zero-counts
    lapply(matrixStats::rowSums2)

  # Count replicates by subgroup
  replicate_counts <- object %>%
    autonomics.import::sdata() %>%
    magrittr::extract2('subgroup') %>%
    table()

  # Which rows contain 1 to replicates-1 zeros (which should be replaced by NA)?
  rows_with_missing_data_by_subgroup <- zeros_per_row_by_subgroup %>%
    names() %>%
    lapply(
      function(sg)
      {
        zeros_per_row_by_subgroup[[sg]] > 0 & zeros_per_row_by_subgroup[[sg]] < replicate_counts[sg]
      }) %>%
    magrittr::set_names(
      zeros_per_row_by_subgroup %>%
        names())

  # Subset zeros in identified rows
  subset_data_by_subroup <- data_by_subgroup %>%
    names() %>%
    lapply(
      function(sg)
      {
        ## Expand row indicators to matrix size
        rpl <- replicate_counts %>%
          magrittr::extract(sg)
        in_row_with_missing_data <- rows_with_missing_data_by_subgroup %>%
          magrittr::extract2(sg) %>%
          matrix(
            ncol = rpl,
            nrow = rows_with_missing_data_by_subgroup %>%
              magrittr::extract2(sg) %>%
              length())
        ## Subset data
        subsetter <- data_by_subgroup %>%
          magrittr::extract2(sg) %>%
          magrittr::equals(0) %>%
          magrittr::and(
            in_row_with_missing_data)
        data_by_subgroup[[sg]][subsetter] <- NA
        ## Return
        return(data_by_subgroup[[sg]])
      }) %>%
    magrittr::set_names(
      data_by_subgroup %>%
        names())

  # Reassemble
  autonomics.import::exprs(object) <- do.call(cbind, subset_data_by_subroup) %>%
    magrittr::extract(,colnames(autonomics.import::exprs(object)))
  autonomics.import::assert_is_valid_eset(object)

  # Return
  if (verbose)  autonomics.support::cmessage('\t\tConvert inconsistent zeroes into NAs')
  return(object)
}
