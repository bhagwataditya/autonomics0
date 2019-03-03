
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
#'         read_proteingroups()
#'    x %>% invert(subgroups = c('E_EM', 'E_BM', 'EM_BM'))
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

   if (length(subgroups)==0) return(x)

  # Assert
  autonomics.import::assert_is_valid_object(x)
  assertive.sets::assert_is_subset('subgroup', autonomics.import::svars(x))
  assertive.sets::assert_is_subset(subgroups, autonomics.import::subgroup_levels(x))

  # Initialize message
  idx <- which(x$subgroup %in% subgroups)
  first <- autonomics.import::exprs(x)[, idx[1]] %>% (function(y) which(!is.na(y))[[1]])
  oldvalue <- autonomics.import::exprs(x)[first, idx[1]] %>% round(2) %>% as.character()
  autonomics.support::cmessage('\t\tInvert subgroups %s', paste0(subgroups, collapse = ', '))

  # Invert (log) ratios
  if (all(autonomics.import::exprs(x) > 0, na.rm = TRUE)){ autonomics.import::exprs(x)[, idx] %<>% (function(x){1/x})
  } else {                                                 autonomics.import::exprs(x)[, idx] %<>% (function(x){ -x})}
  newvalue <- autonomics.import::exprs(x)[first, idx[1]] %>% round(2) %>% as.character()
  autonomics.support::cmessage('\t\t\texprs    : %s -> %s', as.character(oldvalue), as.character(newvalue))

  # Invert subgroup and sampleid values
  oldsubgroups <- autonomics.import::sdata(x)$subgroup[idx]
  newsubgroups <- autonomics.import::sdata(x)$subgroup[idx] %>% vapply(invert.character, character(1), sep = sep)
  oldsampleids <- autonomics.import::sdata(x)$sample_id[idx]
  for (i in seq_along(idx)){
     autonomics.import::sdata(x)$subgroup[ idx[i]] <- newsubgroups[i]
     autonomics.import::sdata(x)$sample_id[idx[i]] %<>% stringi::stri_replace_first_fixed(oldsubgroups[i], newsubgroups[i])
     autonomics.import::snames(x)[         idx[i]] %<>% stringi::stri_replace_first_fixed(oldsubgroups[i], newsubgroups[i])
  }
  newsampleids <- autonomics.import::sdata(x)$sample_id[idx]
  autonomics.support::cmessage('\t\t\tsubgroups: %s -> %s', oldsubgroups[1], newsubgroups[1])
  autonomics.support::cmessage('\t\t\tsampleids: %s -> %s', oldsampleids[1], newsampleids[1])

  # Order on subgroup
  x %<>% autonomics.import::arrange_samples_('subgroup')
  autonomics.support::cmessage('\t\tOrder on subgroup')

  # Return
  return(x)
}
