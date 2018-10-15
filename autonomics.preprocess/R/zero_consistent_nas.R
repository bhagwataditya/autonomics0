#' @rdname zero_consistent_nas
#' @export
zero_to_na_if_not_all_zero <- function(
   no_zero = c('fail', 'warn_passthrough', 'passthrough')[1],
   ...
){
   na_inconsistent_zeroes(...)
}


#' @rdname zero_consistent_nas
#' @importFrom magrittr %<>% %>%
#' @export
na_inconsistent_zeroes <- function(
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


#' Zero consistent NAs - NA inconsistent zeroes
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
#' @examples
#' if (require(subramanian.2016)){
#'    object <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                           package = 'subramanian.2016') %>%
#'              autonomics.import::load_exiqon(infer_design_from_sampleids = TRUE)
#'    object %>% autonomics.preprocess::zero_consistent_nas(verbose = TRUE)
#'
#'    object <- system.file('extdata/metabolon/subramanian.2016.metabolon.xlsx',
#'                           package = 'subramanian.2016') %>%
#'              autonomics.import::load_metabolon(sheet = 5, infer_design_from_sampleids = TRUE)
#'    object %>% autonomics.preprocess::zero_consistent_nas(verbose = TRUE)
#' }
#' @return updated object
#' @importFrom data.table data.table :=
#' @importFrom magrittr %>%
#' @export
zero_consistent_nas <- function(object, verbose = FALSE){

   # Melt
   dt <- object %>% autonomics.import::sumexp_to_long_dt(svars = 'subgroup')

   # Zero consistent NAs
   dt %>% magrittr::extract(, consistent_na := all(is.na(value)), by = c('feature_id', 'subgroup')) %>%
          magrittr::extract(consistent_na==TRUE, value := 0)

   # Report
   n_consistentfeature  <- dt[consistent_na==TRUE, length(unique(feature_id))]
   n_consistentsubgroup <- dt[consistent_na==TRUE, length(unique(subgroup  ))]
   n_subgroup           <- dt[, length(unique(subgroup  ))]
   if (verbose) autonomics.support::cmessage("\t\tZero consistent NAS for %d/%d subgroups in %d/%d features",
                                             n_consistentsubgroup, n_subgroup, n_consistentfeature,  nrow(object))

   # Cast into exprs
   autonomics.import::exprs(object) <- dt %>% data.table::dcast.data.table(feature_id ~ sample_id, value.var = 'value') %>%
                                              autonomics.support::matrixify() %>%
                                              magrittr::extract(rownames(object), ) %>%
                                              magrittr::extract(, colnames(object))

   # Return
   object
}

# Thought of re-writing na_inconsistent_zeroes
# in the same style as zero_consistent_nas(), using data.table,
#
# Started working on this function while working on rnaseq data.
# But then noticed that rnaseq data does not need this.
# So then stopped working
# Will pick up once another need arises.
#
# na_inconsistent_zeroes <- function(object, verbose = FALSE){
#
#    # Melt
#    fid_var <- autonomics.import::fid_var(object)
#    dt <- object %>% autonomics.import::sumexp_to_long_dt(fid = fid_var, svars = 'subgroup')
#    dt %<>% magrittr::extract(order(gene_id))
#
#    # NA inconsistent zeroes
#    dt %>%  magrittr::extract(, inconsistent_0 := value==0 & any(value!=0), by = c(fid_var, 'subgroup'))
#    dt %>%  magrittr::extract(, any_inconsistent_0 := any(inconsistent_0), by='gene_id')
#    dt %<>% magrittr::extract(any_inconsistent_0 == TRUE)
#    mat <- dt %>%  data.table::dcast.data.table(gene_id ~ sample_id, value.var = 'value')
#    mat
# }

