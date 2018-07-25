#' Zero to NA
#'
#' Prominently due to dynamic range limitations in mass spectrometry, data from proteomics
#' and metabolomics LCMS experiments is plagued by missing/dropout values, which occasionaly
#' (e.g. by MaxQuant LFQ) are reported as \code{0} (zero).
#' 
#' Based on the assumption that features truly absent will produce \code{0} (zero) measurements
#' across replicate samples, this function evaluates such data and replcaces \code{0} (zero)
#' with \code{NA} (not available) \bold{only} if any other replicate for the same feature is
#' \bold{not} \code{0} (zero) within the same subgroup.
#'
#' @param object SummarizedExperiment, eSet, or EList
#' @param no_zero \code{\link{character}} defining how to react of there are no \code{0}s (zeros) present in the data set
#' @return updated object
#' @importFrom magrittr %<>% %>%
#' @export
zero_to_na_if_not_all_zero <- function(
  object,
  no_zero = c('fail', 'warn_passthrough', 'passthrough')[1])
{
# Check prerequisites -----------------------------------------------------
  autonomics.import::assert_is_valid_eset(object)
  no_zero %<>%
    match.arg(
      choices = c('fail', 'warn_passthrough', 'passthrough'),
      several.ok = FALSE)
  
  if(
    object %>%
      autonomics.import::exprs() %>%
      magrittr::equals(0) %>%
      any() %>%
      magrittr::not())
  {
    msg <- "'object' contains no '0's (zeros)."
    if(no_zero == 'fail'){ stop(msg) }
    if(no_zero == 'warn_passthrough'){ warning(msg) }
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
  return(object)
}