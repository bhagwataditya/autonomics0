#' All elements of vector are identical
#' @param x vector
#' @return logical(1)
#' @examples
#' require(magrittr)
#' x <- c(2,2,1,2)
#' x %>% has_identical_values()
#' @export
has_identical_values <- function(x){
  length(unique(x))==1
}

#' Convenient equals operator
#' 
#' Performs x == y, but returns FALSE rather than NA for NA elements of x.
#' @param x numeric vector or scalar
#' @param y numeric scalar
#' @examples
#' require(magrittr)
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' y <- 3
#' x %>% equals(y)
#' @export
na_aware_equals <- function(x,y){
  result <- rep(FALSE, length(x)) %>% 
    magrittr::set_names(names(x))
  if (is.na(y)){
    result[ is.na(x)] <- TRUE
    result[!is.na(x)] <- FALSE
  } else {
    result[ is.na(x)] <- FALSE
    result[!is.na(x)] <- x[!is.na(x)] == y
  }
  result
}

#' Is maximal
#' @param x numeric vector
#' @examples 
#' require(magrittr)
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' x %>% is_max()
#' @export
is_max <- function(x){
  na_aware_equals(x, max(x, na.rm = TRUE))
}

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
#'    x %>% guess_sep()
#'
#'    x <- c('WT untreated 1', 'WT untreated 2', 'WT treated 1')
#'    x %>% guess_sep()
#'
#'    x <- c('group1', 'group2', 'group3.R1')
#'    x %>% guess_sep()
#'
#' # SummarizedExperiment
#'    if (require(autonomics.data))   autonomics.data::glutaminase %>%
#'                                    guess_sep()
#'
#'    if (require(autonomics.data))   autonomics.data::stemcomp.proteinratios %>%
#'                                    guess_sep()
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
                magrittr::extract( vapply(., has_identical_values, logical(1)))   %>%
                vapply(unique, integer(1))

   # No separator detected - return NULL
   if (all(sep_freqs==1)){
      if (verbose)  autonomics.support::cmessage('%s: no (consistent) separator. Returning NULL', x[1])
      return(NULL)   # no separator detected
   }

   # Find best separator
   best_sep <- sep_freqs %>%
               magrittr::extract(.!=1)  %>%
               magrittr::extract(is_max(vapply(., magrittr::extract, integer(1), 1)))   %>%
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
   possible_separators = c('.', '_', ' '),# if (contains_ratios(x)) c('.', ' ') else c('.', '_', ' '),
   verbose = FALSE,
   ...
){
   assertive.sets::assert_is_subset(var, c(svars(x), fvars(x)))
  (if (var %in% svars(x)) slevels(x, var) else flevels(x, var)) %>%
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
   sep     = x %>% guess_sep(),
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
   sep          = x %>% guess_sep(),
   verbose      = FALSE,
   ...
){

   # already in x
   if ('subgroup' %in% svars(x)){
      if (verbose) autonomics.support::cmessage("\t\tUse 'subgroup' values in x ")
      return(sdata(x)$subgroup)
   }

   # guess from sampleid values
   x %>% sampleid_values() %>%
         guess_subgroup_values(verbose = verbose)
}


#==========================================================

# guess_subject_values <- function (x, ...) {
#    UseMethod("guess_subject_values", x)
# }
#
# guess_subject_values.character(
#    x,
#    sep     = guess_sep(x),
#    verbose = FALSE
# ){
#    NULL
# }


