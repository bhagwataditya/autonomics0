

#==================================================================================================================

#' @rdname impute_consistent_nas
#' @importFrom magrittr %<>%
#' @export
impute_common_nas <- function(
   object,
   imputefun = function(x) imputeLCMD::impute.QRILC(x)[[1]],
   verbose = FALSE
){

   # Assert no NaN or -Inf (nondetects should be represented as NA)
   assertive.numbers::assert_any_are_not_nan(autonomics.import::exprs(object))
   assertive.numbers::assert_all_are_finite(autonomics.import::exprs(object) %>% c() %>% stats::na.exclude())

   # Typify nondetects
   is_na <- is.na(autonomics.import::exprs(object))
   is_consistent_na <- matrixStats::rowAlls(is_na)
   if (verbose)  autonomics.support::cmessage("\t\tImpute %d/%d features with NA value in all %d samples",
                                               sum(is_consistent_na),
                                               length(is_consistent_na),
                                               ncol(object))
   is_inconsistent_na <- is_na & !matrix(is_consistent_na, nrow(object), ncol(object))

   # Replace inconsistent nondetects by median
   exprs1 <- autonomics.import::exprs(object)
   for (sample in 1:ncol(object)){
      exprs1[, sample] %<>% (function(x){x[is_inconsistent_na[, sample]] <- stats::median(x, na.rm=TRUE); x})
   }

   # Impute consistent nondetects
   autonomics.import::is_imputed(object) <- is.na(exprs1)
   exprs1 %<>% imputefun()

   # Re-NA inconsistent nondetects
   exprs1[is_inconsistent_na] <- NA_real_

   # Update object
   autonomics.import::exprs(object) <- exprs1

   # Return
   object
}


#' Impute consistent NA values
#'
#' Impute values missing in all (subgroup) samples.
#'
#' @param object     SummarizedExperiment
#' @param imputefun  imputation function
#' @param svar character(1)
#' @param verbose    logical(1)
#' @examples
#' require(magrittr)
#' if (require(graumann.lfq)){
#'
#'    # Read object
#'    object <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq') %>%
#'              autonomics::read_proteingroups(impute_consistent_nas = FALSE, plot = FALSE)
#'
#'    # Common NA values - missing in all samples
#'    object %>% autonomics.import::split_by_svar() %>%
#'               magrittr::extract2(2) %>%
#'               autonomics::impute_common_nas(verbose = TRUE)
#'
#'    # Consistent NA values - missing in all subgroup samples
#'    object %>% autonomics::impute_consistent_nas(verbose = TRUE, plot = TRUE)
#' }
#' @return SummarizedExperiment with updated exprs
#' @importFrom magrittr %>%
#' @export
impute_consistent_nas = function(
   object,
   imputefun    = function(x) imputeLCMD::impute.QRILC(x)[[1]],
   svar         = 'subgroup',
   verbose      = FALSE,
   plot         = FALSE
){

   # Check
   if (!svar %in% autonomics.import::svars(object)){
      autonomics.support::cmessage("\t\tNo svar '%s' => no imputation performed", svar)
      return(object)
   }
   if (any(autonomics.support::is_missing_or_empty_character(autonomics.import::svalues(object, svar)))){
      autonomics.support::cmessage("\t\tSome '%s' values are missing => no imputation performed", svar)
      return(object)
   }

   # Impute
   imputed_object <- object %>% autonomics.import::split_by_svar(svar) %>%
                                lapply(autonomics::impute_common_nas, imputefun=imputefun) %>% do.call(SummarizedExperiment::cbind, .)

   # Message
   if (verbose) autonomics.support::cmessage('\t\tImpute consistent NA values among %d/%d features', 
                                              autonomics.import::is_imputed(imputed_object) %>% matrixStats::rowAnys() %>% sum(), 
                                              nrow(object))
   # Plot
   if (plot) imputed_object %>% autonomics::plot_detects_per_subgroup(svar = svar) %>% print()

   # Return
   return(imputed_object)
}


