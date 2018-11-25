#==========================================================================================

#' Count no of components
#' @param x character(n) or SummarizedExperiment
#' @param svar character(1)
#' @param sep character(1)
#' @return integer(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::glutaminase            %>% count_components()
#'    autonomics.data::glutaminase$subgroup   %>% count_components()
#' }
#' @export
count_components <- function (x, ...) {
   UseMethod("count_components", x)
}


#' @rdname count_components
#' @importFrom magrittr %>%
#' @export
count_components.SummarizedExperiment <- function(
   x,
   svar = 'subgroup',
   sep  = x %>% autonomics.import::slevels(svar) %>% autonomics.import::guess_component_sep(.)
){
   n <- x %>%
      autonomics.import::slevels(svar) %>%
      (function(x) if (is.null(sep)) x else x %>% stringi::stri_split_fixed(sep)) %>%
      vapply(length,integer(1)) %>%
      unique()
   assertive.properties::assert_is_scalar(n)
   n
}

#' @rdname count_components
#' @importFrom magrittr %>%
#' @export
count_components.character <- function (x,
                                        sep = x %>% as.character() %>% unique() %>% autonomics.import::guess_component_sep(.)
){
   subgroup_levels <- x %>% as.character() %>% unique()
   n <- subgroup_levels %>% (function(x) if (is.null(sep)) x else x %>% stringi::stri_split_fixed(sep)) %>% vapply(length, integer(1)) %>% unique()
   assertive.properties::assert_is_scalar(n)
   n
}

n_components <- function(...){
   .Deprecated('count_components')
   count_components(...)
}

#===================================================================================================


#' Guess component separator
#' @param x character vector with sampleids
#' @param possible_separators character vector with possible separators to look for
#' @param verbose logical
#' @return separator (string) or NULL (if no separator could be identified)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::glutaminase %>% subgroup_values() %>% guess_component_sep()
#'    autonomics.data::glutaminase %>%                       guess_component_sep()
#' }
#'
#' x <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]', 'PERM_NON.R4[H/L]')
#' x %>% guess_component_sep()
#'
#' x <- c('WT untreated 1', 'WT untreated 2', 'WT treated 1')
#' x %>% guess_component_sep()
#'
#' x <- c('group1', 'group2', 'group3.R1')
#' x %>% guess_component_sep()
#' @export
guess_component_sep <- function (x, ...) {
   UseMethod("guess_component_sep", x)
}


#' @rdname guess_component_sep
#' @importFrom magrittr %>%
#' @export
guess_component_sep.character <- function(x, possible_separators = c('.', ' ', '_'), verbose = FALSE){
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
   if (verbose) autonomics.support::cmessage("\t\tInfer design sep '%s'", best_sep)
   return(best_sep)
}

#' @rdname guess_component_sep
#' @importFrom magrittr %>%
#' @export
guess_component_sep.SummarizedExperiment <- function(x, svar = 'subgroup'){
   x %>% autonomics.import::svalues(svar) %>% guess_component_sep()
}

infer_design_sep <- function(...){
   .Deprecated('guess_component_sep')
   guess_component_sep(...)
}


#===================================================================================================


#' Split components
#'
#' Get data.table in which subgroup components have been decomposed
#'
#' @param x    character(n) or SummarizedExperiment
#' @param sep  character(1): separator
#' @param svar character
#' @return data.table with decomposed subgroup components
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::protein.ratios %>% subgroup_values() %>% split_components()
#'    billing.differentiation.data::protein.ratios %>%                       split_components()
#' }
#' if (require(autonomics.data)){
#'    autonomics.data::glutaminase %>% subgroup_values()  %>% split_components()
#'    autonomics.data::glutaminase                            split_components()
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% subgroup_values()   %>% split_components()
#'    subramanian.2016::metabolon %>%                         split_components()
#' }
#' if (require(graumann.lfq)){
#'    graumann.lfq::lfq.intensities %>% subgroup_values() %>% split_components()
#'    graumann.lfq::lfq.intensities %>%                       split_components()
#' }
#' if (require(atkin.2014)){
#'    object <- 'extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat' %>%
#'               system.file(package='atkin.2014') %>%
#'               autonomics.import::load_soma()
#'    subgroup_values <- object %>% subgroup_values() %>% as.character() %>% unique() %>% split_components()
#'    subgroup_values <- object %>%                                                       split_components()
#' }
#' @export
split_components <- function (x, ...) {
   UseMethod("split_components", x)
}

#' @rdname split_components
#' @importFrom magrittr %>%
#' @export
split_components.character <- function(
   values,
   sep = autonomics.import::guess_component_sep(values, verbose = FALSE)
){

   n <- values %>% count_components()

   dt <- data.table::data.table(subgroup = values)

   # Single component
   if (is.null(sep)){
      dt %>% cbind(V1 = values)
      # Multiple components
   } else {
      dt  %>% cbind(dt[, data.table::tstrsplit(subgroup, sep, fixed = TRUE)])
   }

   # Old approach
   # y <- values %>% stringi::stri_split_fixed(sep)
   # n.component <- length(y[[1]])
   # 1:n.component %>% lapply(function(z) vapply(y, magrittr::extract, character(1), z)) %>%
   #                   data.table::as.data.table()
}


#' @importFrom magrittr %>%
#' @export
split_components.SummarizedExperiment <- function(object, svar = 'subgroup'){
   object %>% autonomics.import::svalues(svar) %>% autonomics.import::split_components()
}


#' @rdname split_components
#' @export
subgroup_components <- function(...){
   .Deprecated('split_components')
   split_components(...)
}

#' @rdname split_components
#' @export
scomponents <- function(...){
   .Deprecated('split_components')
   split_components(...)
}

