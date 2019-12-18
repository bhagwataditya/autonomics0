

#==================================================================================================================

#' @rdname impute_consistent_nas
#' @importFrom magrittr %<>%
#' @export
impute_common_nas <- function(
   object,
   imputefun = function(x) impute_around_zero(x),
   verbose = FALSE
){

   # Assert no NaN or -Inf (nondetects should be represented as NA)
   assertive.numbers::assert_any_are_not_nan(exprs(object))
   assertive.numbers::assert_all_are_finite(exprs(object) %>% c() %>% stats::na.exclude())

   # Typify nondetects
   is_na <- is.na(exprs(object))
   is_consistent_na <- matrixStats::rowAlls(is_na)
   if (verbose)  autonomics.support::cmessage("\t\tImpute %d/%d features with NA value in all %d samples",
                                               sum(is_consistent_na),
                                               length(is_consistent_na),
                                               ncol(object))
   is_inconsistent_na <- is_na & !matrix(is_consistent_na, nrow(object), ncol(object))

   # Replace inconsistent nondetects by median
   exprs1 <- exprs(object)
   for (sample in 1:ncol(object)){
      exprs1[, sample] %<>% (function(x){x[is_inconsistent_na[, sample]] <- stats::median(x, na.rm=TRUE); x})
   }

   # Impute consistent nondetects
   is_imputed(object) <- is.na(exprs1)
   exprs1 %<>% imputefun()

   # Re-NA inconsistent nondetects
   exprs1[is_inconsistent_na] <- NA_real_

   # Update object
   exprs(object) <- exprs1

   # Return
   object
}


#' @rdname filter_samples_available_for_some_feature
#' @importFrom magrittr %>%
#' @export
is_available_for_some_feature <- function(object){
   subsetter <- (!is.na(autonomics.import::exprs(object))) & (autonomics.import::exprs(object) != 0)
   subsetter %>% matrixStats::colAnys() %>% magrittr::set_names(autonomics.import::snames(object))
}

#' Filter samples available for some feature
#' @param object SummarizedExperiment
#' @param verbose TRUE or FALSE
#' @return SummarizedExperiment
#' @export
filter_samples_available_for_some_feature <- function(object, verbose = FALSE){
   subsetter <- object %>% is_available_for_some_feature()
   if (any(!subsetter)){
      if (verbose) autonomics.support::cmessage('\t\tRetain %d/%d samples with a value available for some feature',
                                                sum(subsetter), length(subsetter))
      object %<>% magrittr::extract(, subsetter)
   }
   object
}


#' Impute consistent NA values
#'
#' Impute values missing in all (subgroup) samples.
#'
#' @param object     SummarizedExperiment
#' @param imputefun  imputation function
#' @param svar       string
#' @param verbose    TRUE or FALSE
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'
#'    # Read object
#'    object <- 'extdata/glutaminase/glutaminase.xlsx'     %>%
#'               system.file(package = 'autonomics.data')  %>%
#'               read_metabolon()
#'
#'    # Common NA values - missing in all samples
#'    object %>% split_by_svar() %>%
#'               magrittr::extract2(2) %>%
#'               impute_common_nas(verbose = TRUE)
#'
#'    # Consistent NA values - missing in all subgroup samples
#'    object %>% impute_consistent_nas(verbose = TRUE)
#'    # Use a different imputation method
#'    object %>% impute_consistent_nas(
#'       verbose = TRUE,
#'       imputefun = function(x) imputeLCMD::impute.QRILC(x)[[1]])
#' }
#' @return SummarizedExperiment with updated exprs
#' @importFrom magrittr %>%
#' @export
impute_consistent_nas <- function(
   object,
   imputefun    = function(x) impute_around_zero(x),
   svar         = 'subgroup',
   verbose      = FALSE
){

   # Check
   if (!svar %in% svars(object)){
      autonomics.support::cmessage("\t\tNo svar '%s' => no imputation performed", svar)
      return(object)
   }
   if (any(autonomics.support::is_missing_or_empty_character(svalues(object, svar)))){
      autonomics.support::cmessage("\t\tSome '%s' values are missing => no imputation performed", svar)
      return(object)
   }
   selector <- is_available_for_some_feature(object)
   if (!all(selector)) stop('First rm ', names(selector)[!selector] %>% collapse_words(), ', which contain only NA/0 values')

   # Impute
   imputed_object <- object %>% split_by_svar(svar) %>%
                                lapply(impute_common_nas, imputefun=imputefun) %>% do.call(SummarizedExperiment::cbind, .)

   # Message
   if (verbose) autonomics.support::cmessage('\t\tImpute consistent NA values among %d/%d features',
                                              is_imputed(imputed_object) %>% matrixStats::rowAnys() %>% sum(),
                                              nrow(object))
   # Return
   return(imputed_object)
}


# Linguistic collapse
# x <- c('a', 'b', 'c')
# x %>% collapse_words()
collapse_words <- function(x, collapsor = 'and'){
   if (length(x) == 1) return(x)
   if (length(x) == 2) return(sprintf('%s %s %s', x[[1]], collapsor, x[[2]]))
   paste0(x[-length(x)], collapse = ', ') %>% paste0(' and ', x[[3]])
}

#' Impute around zero
#' @param x exprs matrix
#' @return exprs matrix
#' @examples
#' require(magrittr)
#' x <- matrix(c(20, 21, 22,
#'               30, 31, 32,
#'               NA, NA, NA,
#'               40, NA, 41,
#'               NA, NA, NA,
#'               NA, 50, NA), ncol=3, byrow=TRUE)
#' x %>% impute_around_zero()
#' object <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(exprs = x))
#' object %>% autonomics.import::impute_common_nas(
#'   imputefun = impute_around_zero) %>%
#'   autonomics.import::exprs()
#' object$subgroup <- rep('EV', 3)
#' object %>% autonomics.import::impute_consistent_nas(
#'   imputefun = impute_around_zero) %>%
#'   autonomics.import::exprs()
#' @importFrom magrittr %>% %<>%
#' @export
impute_around_zero <- function(x){
   meansd <- x %>%
      matrixStats::rowSds() %>%
      mean(na.rm = TRUE)
   is_common_na_feature <- x %>% is.na() %>% matrixStats::rowAlls()
   x[is_common_na_feature, ] %<>% apply(
      1,
      function(y) impute_row(y, meansd)) %>% t()
   x
}

impute_row <- function(x, sd){
   abs(rnorm(length(x), sd = sd * 2))
}
