utils::globalVariables('.')

#' Invert collapsed strings
#' @param x vector of collapsed strings
#' @param sep separator used to collapse strings
#' @examples
#'  invert_collapsed_strings(c('Ctrl_A', 'Ctrl_B'), '_')
#' @importFrom magrittr %>%
#' @export
invert_collapsed_strings <- function(x, sep){

   nsep <- x %>% stringi::stri_count_regex(sep) %>% magrittr::equals(1) %>%
           magrittr::set_names(x)
   if(!all(nsep)){
      stop()
   }

   assertive.base::assert_all_are_true(nsep)

  as.character(x) %>%
    vapply(
      function(x){
        unlist(strsplit(x, sep)) %>%
          magrittr::extract(length(.):1)     %>%
          paste(collapse = sep)
      },
      character(1)) %>%
    unname()
}


#' Invert ratios for some subgroup levels
#' @param object        eset
#' @param invert_subgroups vector of subgroup levels that require inversion
#' @param channel_frac  fractionator symbol for channel values
#' @param subgroup_frac fractionator symbol for subgroup values
#' @return eset with inverted ratios for specified subgroups
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    proteingroups_file <- system.file('extdata/billing2016/proteinGroups.txt',
#'                                       package = 'autonomics.data')
#'    object <- autonomics.import::load_proteingroups(proteingroups_file)
#'    invert_subgroups <- unique(object$subgroup)
#'    object$subgroup;  subgroup_frac <- '_'
#'    object$channel;   channel_frac <- '/'
#'    object %>% autonomics.preprocess::invert_ratios(invert_subgroups, channel_frac, subgroup_frac)
#' }
#' @importFrom dplyr     n
#' @importFrom magrittr  %>%   %<>%
#' @export
invert_ratios <- function(
   object,
   invert_subgroups,
   channel_frac,
   subgroup_frac
){

   if (is.null(invert_subgroups)){
      return(object)
   }

  # Assert
  autonomics.import::assert_is_valid_eset(object)
  assertive.base::assert_is_identical_to_true('subgroup' %in% autonomics.import::svars(object))
  assertive.base::assert_all_are_true(invert_subgroups %in% object$subgroup)

  # Select
  selector <- object$subgroup %in% invert_subgroups

  # Invert (log) ratios
  if (all(autonomics.import::exprs(object) > 0, na.rm = TRUE)){
     autonomics.support::cmessage('\t\tInvert subgroups %s (ratios)',     paste0(invert_subgroups, collapse = ', '))
    autonomics.import::exprs(object)[, selector] %<>% (function(x){1/x})  # ratios
  } else {
     autonomics.support::cmessage('\t\tInvert subgroups %s (log ratios)', paste0(invert_subgroups, collapse = ', '))
     autonomics.import::exprs(object)[, selector] %<>% (function(x){-x})   # log ratios
  }

  # Invert labels and subgroup
  sdata1 <- autonomics.import::sdata(object)
  sdata1$channel  %<>% as.character()
  sdata1$subgroup %<>% as.character()
  sdata1$channel[ selector] %<>% autonomics.preprocess::invert_collapsed_strings(channel_frac)
  sdata1$subgroup[selector] %<>% autonomics.preprocess::invert_collapsed_strings(subgroup_frac)

  # Redefine replicates and sample names
  sdata1 %<>% dplyr::group_by_('subgroup') %>%    # must be character to allow mapping to shape in ggplot
                                dplyr::mutate(replicate = as.character(seq_len(n()))) %>%
                                dplyr::ungroup() %>%
                                as.data.frame()
  sdata1$sample_id <- paste0(object$subgroup, '.R', object$replicate)
  autonomics.import::sdata(object)  <- sdata1
  autonomics.import::snames(object) <- sdata1$sample_id

  # Order on subgroup
  object %<>% autonomics.import::arrange_samples_(c('subgroup', 'replicate'))

  # Return
  return(object)
}
