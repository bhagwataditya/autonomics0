# Impute ProtSetSet expressions
#
# Imputes missing values in the "exprs" element of a ProtSet
# @param pset A ProtSet.
# @param scale A string denoting whether or not ratios and intensities in the
# ProtSet are on a natural or logarithmic scale.
# @return The input ProtSet, but with no missing values in the \code{exprs}
# element.
# @details This algorithm uses the \code{exprs}, \code{intensities_num}, and
# \code{intensities_den} elements of the ProtSet.
# 1. If the ratios and intensities are on a natural scale, they are converted
# to a log scale.
# 2. \code{intensities_num} that are -Inf (i.e., that were 0 on a natural
# scale) are replaced by half the smallest finite values of
# \code{intensities_num}.  (That is, an estimate of half the limit of
# detection.)
# 3. Likewise, \code{intensities_den}  that are -Inf are replaced by half the
# smallest finite values of \code{intensities_den}.
# 4. NaN ratios are replaced by \code{intensities_num / intensities_den}.
# 5. If the original scale was natural, this is restored.
# 6. The new imputed ratios are merged back into the \code{exprs} element of
# the ProtSet.
impute_exprs <- function(pset, scale = c("log2", "natural")) {
   
   # Satisfy CHECK
   # These ProtSet/PhosSet methods are now obsolete. They refered to the numerator
   # and denominator intensities of a labeled LCMS proteomics experiment
   # The impute_exprs function has been retained in this package
   # To make it work, however, new implementations of intensities_num etc must be developed.
   intensities_num    <-  intensities_den    <- NULL
  `intensities_num<-` <- `intensities_den<-` <- NULL
   
  scale <- match.arg(scale)
  exprs_pset <- Biobase::exprs(pset)
  # For convenience of reasoning about the algorithm, always use a log-scale.
  # This is computationally inefficient, so it might need changing at some point.
  if(scale == "natural"){
    exprs_pset <- log2(exprs_pset)
    intensities_num(pset) <- log2(intensities_num(pset))
    intensities_den(pset) <- log2(intensities_den(pset))
  }
  # Get the ratios that are NaN, along with the corresponding intensity
  # numerators and denominators, and their sample IDs
  isnan_ij <- which(is.nan(exprs_pset), arr.ind = TRUE)
  intensity_parts <- dplyr::data_frame(
    feature = rownames(exprs_pset)[isnan_ij[, "row"]],
    sample = colnames(exprs_pset)[isnan_ij[, "col"]],
    num = intensities_num(pset)[isnan_ij],
    den = intensities_den(pset)[isnan_ij]
  )
  
  # In the case of missing values, the numerator or denominator is sometimes
  # `-Inf` (on a log scale, `0` on a natural scale).  We substitute half the
  # smallest observed value for these cases.  The ratio is now calculated as the
  # numerator divided by the denominator.
  half_smallest_num <- 0.5 * intensity_parts %>%
    dplyr::filter_(~ is.finite(num)) %>%
    dplyr::summarize_(~ min(num)) %>%
    magrittr::extract2(1)
  
  half_smallest_den <- 0.5 * intensity_parts %>%
    dplyr::filter_(~ is.finite(den)) %>%
    dplyr::summarize_(~ min(den)) %>%
    magrittr::extract2(1)
  
  new_exprs <- intensity_parts %>%
    dplyr::mutate_(
      new_value = ~ ifelse(is.finite(num), num, half_smallest_num) /
        ifelse(is.finite(den), den, half_smallest_den)
    ) %>%
    # Convert back to matrix form
    reshape2::acast(feature ~ sample, value.var = "new_value")
  
  # Fold back into the ExpressionSet
  combined_exprs <- ifelse(is.nan(exprs_pset), new_exprs, exprs_pset)
  # Undo the log-scale, if necessary
  if(scale == "natural"){
    combined_exprs <- 2 ^ combined_exprs
  }
  Biobase::exprs(pset) <- combined_exprs
  pset
}
