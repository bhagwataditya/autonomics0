#======================================================
# MERGE SDATA
#======================================================

#' Merge sdata into SummarizedExperiment
#'
#' @details Common variables (other than by) are removed from newdata prior to merge.
#' @param object         SummarizedExperiment
#' @param newdata        dataframe
#' @param by             variable by which to merge
#' @param newdata_first  logical: whether to have newdata first in sdata
#' @param verbose        logical
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- system.file('extdata/stemcomp/soma/stemcomp.adat',
#'                           package = 'autonomics.data') %>%
#'              autonomics::read_somascan()
#'    newdata <- autonomics.import::sdata(object) %>%
#'               magrittr::extract(, 'sample_id', drop = FALSE) %>%
#'               cbind(Letter = letters[1:nrow(.)])
#'    object %>% autonomics.import::sdata()
#'    object %>% autonomics::merge_sdata(newdata) %>%
#'               autonomics.import::sdata()
#'    object %>% autonomics.import::merge_sdata(newdata, newdata_first = TRUE) %>% autonomics.import::sdata()
#' }
#' @importFrom magrittr %>%
#' @export
merge_sdata <- function(object, newdata, by = 'sample_id', newdata_first = FALSE, verbose = TRUE){

   # Assert
   assertive.sets::assert_is_subset(by, autonomics.import::svars(object))
   assertive.sets::assert_is_subset(by, names(newdata))
   assertive.strings::assert_all_are_non_missing_nor_empty_character(object    %>% autonomics.import::svalues(by) %>% as.character())
   assertive.strings::assert_all_are_non_missing_nor_empty_character(newdata %>% magrittr::extract2(by) %>% as.character())
   assertive.sets::assert_are_set_equal(object %>% autonomics.import::svalues(by), newdata %>% magrittr::extract2(by))

   # Common variables
   common_vars <- names(newdata) %>% setdiff(by) %>% intersect(autonomics.import::svars(object))
   if (length(common_vars)>0){
      if (verbose) autonomics.support::cmessage("Ignore from design file (already in sumexp): %s",
                                                 sprintf("'%s'", common_vars) %>% paste0(collapse = ', '))
      newdata %<>% magrittr::extract(, setdiff(names(.), common_vars), drop = FALSE)
   }

   # Merge into sdata
   if (ncol(newdata)>1){
      autonomics.import::sdata(object) %<>% autonomics.support::left_join_keeping_rownames(newdata, by = by)
      if (newdata_first) autonomics.import::sdata(object) %<>%
                         autonomics.support::pull_columns(names(newdata) %>% setdiff(by))
   }

   # Return
   object
}

#======================================================
# MERGE FDATA
#======================================================

#' Merge fdata into SummarizedExperiment
#'
#' @details Common variables (other than by) are removed from newdata prior to merge.
#' @param object         SummarizedExperiment
#' @param newdata        dataframe
#' @param by             variable by which to merge
#' @param newdata_first  logical: whether to have newdata first in sdata
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- system.file('extdata/stemcomp/soma/stemcomp.adat',
#'                           package = 'autonomics.data') %>%
#'              autonomics::read_somascan()
#'    newdata <- autonomics.import::fdata(object)           %>%
#'               magrittr::extract(, 'feature_id', drop = FALSE) %>%
#'               cbind(Letter = sample(letters, nrow(.), TRUE))
#'    object %>% autonomics.import::fdata() %>% head()
#'    object %>% autonomics::merge_fdata(newdata) %>%
#'               autonomics.import::fdata() %>% head()
#'    object %>% autonomics.import::merge_fdata(newdata, newdata_first = TRUE) %>%
#'               autonomics.import::fdata() %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
merge_fdata <- function(object, newdata, by = 'feature_id', newdata_first = FALSE){

   # Assert
   assertive.sets::assert_is_subset(by, autonomics.import::fvars(object))
   assertive.sets::assert_is_subset(by, names(newdata))
   assertive.strings::assert_all_are_non_missing_nor_empty_character(object  %>% autonomics.import::fvalues(by) %>% as.character())
   assertive.strings::assert_all_are_non_missing_nor_empty_character(newdata %>% magrittr::extract2(by) %>% as.character())
   assertive.sets::assert_are_set_equal(object %>% autonomics.import::fvalues(by), newdata %>% magrittr::extract2(by))

   # Common variables
   common_vars <- names(newdata) %>% setdiff(by) %>% intersect(autonomics.import::fvars(object))
   if (length(common_vars)>0){
      autonomics.support::cmessage("Ignore from design file (already in sumexp): %s",
                                   sprintf("'%s'", common_vars) %>% paste0(collapse = ', '))
      newdata %<>% magrittr::extract(, setdiff(names(.), common_vars), drop = FALSE)
   }

   # Merge into sdata
   if (ncol(newdata)>1){
      autonomics.import::fdata(object) %<>% autonomics.support::left_join_keeping_rownames(newdata, by = by)
      if (newdata_first) autonomics.import::fdata(object) %<>%
                         autonomics.support::pull_columns(names(newdata) %>% setdiff(by))
   }

   # Return
   object
}

