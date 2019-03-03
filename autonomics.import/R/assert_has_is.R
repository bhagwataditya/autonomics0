
#=====================================================================================
# is

is_elist <- function(x, .xname = assertive.base::get_name_in_parent(x)){
   if (!methods::is(x, 'EList')){
      return(assertive.base::false('%s is not an EList', .xname))
   }
   TRUE
}

is_eSet <- function(x, .xname = assertive.base::get_name_in_parent(x)){
   if (!methods::is(x, 'eSet')){
      return(assertive.base::false('%s is not an eSet', .xname))
   }
   TRUE
}

is_summarized_experiment <- function(x, .xname = assertive.base::get_name_in_parent(x)){
   if (!methods::is(x, 'SummarizedExperiment')){
      return(assertive.base::false('%s is not a SummarizedExperiment', .xname))
   }
   TRUE
}

#=====================================================================================
# has/contains

has_valid_featureNames <- function(x, .xname = assertive.base::get_name_in_parent(x)){
   # SummarizedExperiments do not allow row naming of fdata
   # if (!all(autonomics.import::fnames(x) == rownames(autonomics.import::fdata(x)))){
   #   return(assertive.base::false('fnames(%s) differ from rownames(fdata(%s))', .xname, .xname))
   # }
   if (!all(autonomics.import::fnames(x) == rownames(autonomics.import::exprs(x)))){
      return(assertive.base::false('fnames(%s) differ from rownames(exprs(%s))', .xname, .xname))
   }
   TRUE
}

has_valid_sampleNames <- function(x, .xname = assertive.base::get_name_in_parent(x)){
   if (!all(autonomics.import::snames(x) == rownames(autonomics.import::sdata(x)))){
      return(assertive.base::false('snames(%s) differ from rownames(sdata(%s))', .xname, .xname))
   }
   if (!all(autonomics.import::snames(x) == colnames(autonomics.import::exprs(x)))){
      return(assertive.base::false('snames(%s) differ from colnames(exprs(%s))', .xname, .xname))
   }
   TRUE
}

#' Does object have complete svalues
#' @param object SummarizedExperiment
#' @param svar   sample var
#' @return logical
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% autonomics.import::has_complete_subgroup_values()
#'    object %>% autonomics.import::has_complete_block_values()
#' }
#' @importFrom magrittr %>%
#' @export
has_complete_svalues <- function(object, svar){

   # svar missing
   var_present <- svar %in% autonomics.import::svars(object)
   if (!var_present) return(FALSE)

   # svalues missing
   values_present <- object %>% autonomics.import::svalues(svar) %>%
                                assertive.strings::is_empty_character() %>%
                                any()
   if (values_present) return(FALSE)

   # svar and svalues both present
   return(TRUE)
}

#' @rdname has_complete_svalues
#' @importFrom magrittr %>%
#' @export
has_complete_subgroup_values <- function(object){
   object %>% autonomics.import::has_complete_svalues('subgroup')
}

#' @rdname has_complete_svalues
#' @importFrom magrittr %>%
#' @export
has_complete_block_values <- function(object){
   object %>% autonomics.import::has_complete_svalues('block')
}

#' @title Does object contain prepro?
#' @description Does object contain preprocessing info?
#' @param object SummarizedExperiment
#' @return logical
#' @export
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>% contains_prepro()
#' }
contains_prepro <- function(object){
   autonomics.import::assert_is_valid_object(object)
   length(autonomics.import::prepro(object))!=0
}

#' Does object contain ratio values?
#' @param object SummarizedExperiment
#' @return logical
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>%
#'    autonomics.import::contains_ratios()
#'
#'    autonomics.data::stemcomp.soma %>%
#'    autonomics.import::contains_ratios()
#' }
#'
#' # STEM CELL DIFFERENTIATION
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::protein.ratios %>%
#'       autonomics.import::contains_ratios()
#'    billing.differentiation.data::rna.voomcounts %>%
#'       autonomics.import::contains_ratios()
#' }
#' @export
contains_ratios <- function(object){
   autonomics.import::assert_is_valid_object(object)
   if (autonomics.import::contains_prepro(object)){
      autonomics.import::prepro(object)$quantity %in% c('Ratio', 'Ratio normalized', 'occupancy')
   } else {
      FALSE
   }
}

#=========================================================================================
# Assert

#' Is valid eset?
#' @param x eset
#' @param .xname see assertive.base::get_name_in_parent
#' @return logical
#' @export
is_valid_object <- function(x, .xname = assertive.base::get_name_in_parent(x)){
   if (!(ok <- is_eSet(x, .xname = .xname)  |  is_elist(x, .xname = .xname)  |  is_summarized_experiment(x, .xname = .xname))){   return(ok)}
   if (!(ok <- has_valid_featureNames(x, .xname = .xname)))                  {   return(ok)}
   if (!(ok <- has_valid_sampleNames(x,  .xname = .xname)))                  {   return(ok)}
   TRUE
}

#' Assert that x is a valid eSet
#' @param x eset
#' @return error if not true
#' @export
assert_is_valid_object <- function(x){
  assertive.base::assert_engine(is_valid_object, x, .xname = assertive.base::get_name_in_parent(x))
}

#' Assert that features are valid
#' @param features \code{\link{numeric}} index or \code{\link{character}} names
#' of \code{\link[autonomics.import]{fdata}} entries
#' @param object eSet
#' @importFrom magrittr %>%
#' @export
assert_all_are_valid_features <- function(features, object){
   if (is.character(features)) {
      features %>%
      assertive.sets::assert_is_subset(
         object %>%
         autonomics.import::fdata() %>%
         magrittr::extract2("feature_id"))
   } else if (is.numeric(features)) {
      features %>%
      assertive.numbers::assert_all_are_whole_numbers() %>%
      assertive.numbers::assert_all_are_in_closed_range(1, nrow(object))
   } else {
      stop("Features must be numeric indexes or character objects.")
   }
   features %>%
   invisible()
}


