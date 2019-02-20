#===================================================================================================
#' Guess separator
#' @param x                   character vector or SummarizedExperiment
#' @param var                 svar or fvar
#' @param possible_separators character vector with possible separators to look for
#' @param verbose             logical
#' @param ...                 used for proper S3 method dispatch
#' @return separator (string) or NULL (if no separator could be identified)
#' @examples
#' require(magrittr)
#'
#' # charactervector
#'    x <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]', 'PERM_NON.R4[H/L]')
#'    x %>% autonomics.import::guess_sep()
#'
#'    x <- c('WT untreated 1', 'WT untreated 2', 'WT treated 1')
#'    x %>% autonomics.import::guess_sep()
#'
#'    x <- c('group1', 'group2', 'group3.R1')
#'    x %>% autonomics.import::guess_sep()
#'
#' # SummarizedExperiment
#'    if (require(autonomics.data))   autonomics.data::glutaminase %>%
#'                                    autonomics.import::guess_sep()
#'
#'    if (require(autonomics.data))   autonomics.data::stemcomp.proteinratios %>%
#'                                    autonomics.import::guess_sep()
#'
#'    if (require(graumann.lfq))      graumann.lfq::lfq.intensities %>%
#'                                    autonomics.import::guess_sep()
#' @export
guess_sep <- function (x, ...) {
   UseMethod("guess_sep", x)
}


#' @rdname guess_sep
#' @importFrom magrittr %>%
#' @export
guess_sep.character <- function(
   x,
   possible_separators = c('.', ' ', '_'),
   verbose = FALSE,
   ...
){
   . <- NULL
   sep_freqs <- Map(function(y) stringi::stri_split_fixed(x, y), possible_separators)        %>%
                lapply(function(y) y %>% vapply(length, integer(1)))                                  %>%
                magrittr::extract( vapply(., autonomics.support::has_identical_values, logical(1)))   %>%
                vapply(unique, integer(1))

   # No separator detected - return NULL
   if (all(sep_freqs==1)){
      if (verbose)  autonomics.support::cmessage('%s: no (consistent) separator. Returning NULL', x[1])
      return(NULL)   # no separator detected
   }

   # Find best separator
   best_sep <- sep_freqs %>%
               magrittr::extract(.!=1)  %>%
               magrittr::extract(autonomics.support::is_max(vapply(., magrittr::extract, integer(1), 1)))   %>%
               names()

   # Ambiguous separator - take first from tail
   if (length(best_sep)>1){
      pattern <- best_sep %>% paste0(collapse='') %>% paste0('[', ., ']')
      best_sep <- x[1] %>% stringi::stri_extract_last_regex(pattern)
   }

   # Separator identified - return
   if (verbose) autonomics.support::cmessage("\t\tGuess sep: '%s'", best_sep)
   return(best_sep)
}

#' @rdname guess_sep
#' @importFrom magrittr %>%
#' @export
guess_sep.factor <- function(x, ...) x %>% levels %>% guess_sep.character()


#' @rdname guess_sep
#' @importFrom magrittr %>%
#' @export
guess_sep.SummarizedExperiment <- function(
   x,
   var = 'sample_id',
   possible_separators = c('.', '_', ' '),# if (autonomics.import::contains_ratios(x)) c('.', ' ') else c('.', '_', ' '),
   verbose = FALSE,
   ...
){
   assertive.sets::assert_is_subset(var, c(autonomics.import::svars(x), autonomics.import::fvars(x)))
  (if (var %in% autonomics.import::svars(x)) autonomics.import::slevels(x, var) else autonomics.import::flevels(x, var)) %>%
   guess_sep(possible_separators = possible_separators,
             verbose             = verbose)
}

infer_design_sep <- function(...){
   .Deprecated('guess_sep')
   guess_sep(...)
}

#' @rdname guess_sep
#' @importFrom magrittr %>%
#' @export
ssep <- function(...){
   .Deprecated('guess_sep')
   guess_sep(...)
}

#' @rdname guess_sep
#' @importFrom  magrittr %>%
#' @export
subgroup_sep <- function(...){
   .Deprecated('guess_sep')
   guess_sep(...)
}



#=================================================================================
#' Guess subgroup values
#' @param x             charactervector, SummarizedExperiment
#' @param sep           character(1)
#' @param verbose       logical(1)
#' @param ...           used for proper S3 method dispatch
#' @return character(n)
#' @examples
#' require(magrittr)
#'
#' # charactervector
#'    # No sep: subgroup = x
#'       x <- c("EM00", "EM01", "EM02")
#'       x %>% guess_subgroup_values()
#'
#'    # Sep: subgroup = head components of x
#'       x <- c("UT_10h_R1", "UT_10h_R2", "UT_10h_R3")
#'       x %>% guess_subgroup_values()
#'
#'       x <- c("EM00_STD.R1", "EM01_STD.R1", "EM01_EM00.R1")
#'       x %>% guess_subgroup_values()
#'
#' @export
guess_subgroup_values <- function (x, ...) {
   UseMethod("guess_subgroup_values", x)
}

#' @rdname guess_subgroup_values
#' @importFrom magrittr %>%
#' @export
guess_subgroup_values.character <- function(
   x,
   sep     = x %>% autonomics.import::guess_sep(),
   verbose = FALSE,
   ...
){
   # Guess
   subgroup_values <- if (is.null(sep)){ x
                      } else {           x %>% stringi::stri_split_fixed(sep) %>%
                                               vapply(function(y) y %>%
                                                                  magrittr::extract(1:(length(y)-1)) %>%
                                                                  paste0(collapse = sep), character(1))
                      }
   # Inform
   if (verbose)   autonomics.support::cmessage('\t\tGuess subgroup values: %s => %s', x[1], subgroup_values[1])

   # Return
   return(subgroup_values)
}

#' @rdname guess_subgroup_values
#' @importFrom magrittr %>%
#' @export
guess_subgroup_values.SummarizedExperiment <- function(
   x,
   sep          = x %>% autonomics.import::guess_sep(),
   verbose      = FALSE,
   ...
){

   # already in x
   if ('subgroup' %in% autonomics.import::svars(x)){
      if (verbose) autonomics.support::cmessage("\t\tUse 'subgroup' values in x ")
      return(autonomics.import::sdata(x)$subgroup)
   }

   # guess from sampleid values
   x %>% autonomics.import::sampleid_values() %>%
         autonomics.import::guess_subgroup_values(verbose = verbose)
}


#==========================================================

# guess_subject_values <- function (x, ...) {
#    UseMethod("guess_subject_values", x)
# }
#
# guess_subject_values.character(
#    x,
#    sep     = autonomics.import::guess_sep(x),
#    verbose = FALSE
# ){
#    NULL
# }


