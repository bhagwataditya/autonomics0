
#' Invert
#'
#' For character input: invert collapsed strings.
#' For SummarizedExperiment input: invert expression values, invert subgroup definitions, update sample id values.
#'
#' @param  x          character(.) or SummarizedExperiment
#' @param  sep        character(1): collapsed string separator
#' @param  subgroups  character(.): subgroup levels for which inversion is required
#' @param  ... to enable S3 method dispatch
#' @return character(.) or SummarizedExperiment
#' @examples
#' require(magrittr)
#'
#' # character
#' #----------
#'    x <- c('Ctrl_A', 'Ctrl_B')
#'    x %>% invert()
#'
#' # SummarizedExperiment
#' #---------------------
#' if (require(autonomics.data)){
#'    x <- system.file('extdata/stemcomp/maxquant/proteinGroups.txt',
#'                      package = 'autonomics.data') %>%
#'         autonomics::read_proteingroups()
#'
#'    x %>% autonomics.import::sdata()
#'
#'    x %>% autonomics::invert(subgroups = c('E_EM', 'E_BM', 'EM_BM')) %>%
#'          autonomics.import::sdata()
#' }
#' @importFrom magrittr %>%
#' @export
invert <- function(x, ...){
   UseMethod('invert', x)
}


#' @rdname invert
#' @export
invert.character <- function(
   x,
   sep = autonomics.import::guess_sep(x),
   ...
){
   x                              %>%
   stringi::stri_split_fixed(sep) %>%
   lapply(rev)                    %>%
   vapply(paste, character(1), collapse = sep)
}


#' @rdname invert
#' @importFrom magrittr  %>%   %<>%
#' @export
invert.SummarizedExperiment <- function(
   x,
   subgroups = x %>% autonomics.import::slevels('subgroup'),
   sep       = x %>% autonomics.import::guess_sep('subgroup'),
   ...
){

   if (is.null(invert_subgroups)){
      return(x)
   }

  # Assert
  autonomics.import::assert_is_valid_eset(x)
  assertive.sets::assert_is_subset('subgroup', autonomics.import::svars(x))
  assertive.sets::assert_is_subset(subgroups, autonomics.import::subgroup_levels(x))

  # Invert (log) ratios
  idx <- which(x$subgroup %in% subgroups)
  if (all(autonomics.import::exprs(x) > 0, na.rm = TRUE)){ # Ratios
     autonomics.support::cmessage('\t\tInvert subgroups %s: exprs = 1/exprs', paste0(subgroups, collapse = ', '))
     autonomics.import::exprs(x)[, idx] %<>% (function(x){1/x})
  } else {                                                 # Log Ratios
     autonomics.support::cmessage('\t\tInvert subgroups %s: exprs = -exprs',  paste0(subgroups, collapse = ', '))
     autonomics.import::exprs(x)[, idx] %<>% (function(x){-x})
  }

  # Invert subgroup and sampleid values
  for (i in idx){
     oldsubgroup <- autonomics.import::sdata(x)$subgroup[i]
     newsubgroup <- autonomics.import::sdata(x)$subgroup[i] %>% invert.character(sep = sep)
     autonomics.import::sdata(x)$subgroup[i] <- newsubgroup
     autonomics.import::sdata(x)$sample_id[i] %<>% stringi::stri_replace_first_fixed(oldsubgroup, newsubgroup)
     autonomics.import::snames(x)[i]          %<>% stringi::stri_replace_first_fixed(oldsubgroup, newsubgroup)
  }

  # Order on subgroup
  x %<>% autonomics.import::arrange_samples_('subgroup')

  # Return
  return(x)
}
