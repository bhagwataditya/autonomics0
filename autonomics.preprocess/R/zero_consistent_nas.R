#' NA inconsistent zeroes
#'
#' These functions allow to differentiate between two types of non-detects: consistent (which are 0 rather than NA) and
#' inconsistent (NA rather than 0). This distinction is required in LCMS platforms (proteomics and metabolomics),
#' where a high degree of inconsistent non-detects are generated due to platform stochasticity. It is not required in
#' transcriptomics platforms (RNAseq, qPCR, microarrays), where non-detecs arise due to low abundance rather than platform
#' stochasticity.
#'
#' @param object SummarizedExperiment
#' @param no_zero \code{\link{character}} defining how to react of there are no \code{0}s (zeros) present in the data set
#' @return updated object
#' @return updated object
#' @importFrom data.table data.table :=
#' @importFrom magrittr %<>% %>%
#' @export
zero_to_na_if_not_all_zero <- function(
   object,
   no_zero = c('fail', 'warn_passthrough', 'passthrough')[1]
){

# Check prerequisites -----------------------------------------------------
  subgroup <- NULL
  autonomics.import::assert_is_valid_object(object)
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
    magrittr::set_names(autonomics.import::subgroup_levels(object)) %>%
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
  autonomics.import::assert_is_valid_object(object)

  # Return
  autonomics.support::cmessage('\t\tConvert inconsistent zeroes into NAs')
  return(object)
}

