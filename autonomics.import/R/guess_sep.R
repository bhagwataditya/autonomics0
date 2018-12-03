#===================================================================================================


#' Guess separator
#' @param x                   character vector or SummarizedExperiment
#' @param svar                string
#' @param possible_separators character vector with possible separators to look for
#' @param verbose             logical
#' @return separator (string) or NULL (if no separator could be identified)
#' @examples
#' require(magrittr)
#'
#' x <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]', 'PERM_NON.R4[H/L]')
#' x %>% guess_sep()
#'
#' x <- c('WT untreated 1', 'WT untreated 2', 'WT treated 1')
#' x %>% guess_sep()
#'
#' x <- c('group1', 'group2', 'group3.R1')
#' x %>% guess_sep()
#'
#' if (require(autonomics.data))   autonomics.data::glutaminase %>%
#'                                 autonomics.import::guess_sep()
#'
#' if (require(autonomics.data))   autonomics.data::stemcomp.proteinratios %>%
#'                                 autonomics.import::guess_sep()
#'
#' if (require(subramanian.2016))  subramanian.2016::metabolon  %>%
#'                                 autonomics.import::guess_sep()
#'
#' if (require(graumann.lfq))      graumann.lfq::lfq.intensities %>%
#'                                 autonomics.import::guess_sep()
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
   verbose = FALSE
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

   # Ambiguous separator - return NULL
   if (length(best_sep)>1){
      if (verbose)   autonomics.support::cmessage('%s: %s separator? Returning NULL.',
                                                  x[1],
                                                  paste0(sprintf("'%s'", best_sep), collapse = ' or '))
      return(NULL)  # ambiguous separator
   }

   # Separator identified - return
   if (verbose) autonomics.support::cmessage("\t\tGuess sep: '%s'", best_sep)
   return(best_sep)
}


#' @rdname guess_sep
#' @importFrom magrittr %>%
#' @export
guess_sep.SummarizedExperiment <- function(
   x,
   svar = 'subgroup',
   possible_separators = if (autonomics.import::contains_ratios(x)) c('.', ' ') else c('.', '_', ' '),
   verbose = FALSE
){
   x %>%
   autonomics.import::svalues(svar) %>%
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



