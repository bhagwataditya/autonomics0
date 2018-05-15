
#' Infer design separator
#' @param sample_ids character vector with sample ids
#' @param possible_separators character vector with possible separators to look for
#' @examples
#' require(magrittr)
#' sample_ids <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]', 'PERM_NON.R4[H/L]')
#' sample_ids %>% infer_design_sep()
#' sample_ids <- c('WT untreated 1', 'WT untreated 1', 'WT untreated 2')
#' sample_ids %>% infer_design_sep()
#' @importFrom magrittr %>%
#' @export
infer_design_sep <- function(sample_ids, possible_separators = c('.', ' ', '_')){
   . <- NULL
   Map(function(x) stringi::stri_split_fixed(sample_ids, x), possible_separators)             %>%
   lapply(function(x) x %>% vapply(length, integer(1)))                                       %>%
   magrittr::extract( vapply(., autonomics.support::has_identical_values, logical(1)   ))     %>%
   magrittr::extract(!vapply(., autonomics.support::all_have_value,       logical(1), 1))     %>%
   magrittr::extract(autonomics.support::is_max(vapply(., magrittr::extract, integer(1), 1))) %>%
   names() %>%
   (function(x) if (length(x)>0) x[1])
}

#' Convert sample IDs into sample design
#' @param sample_ids character vector with sample ids
#' @examples
#' require(magrittr)
#' sample_ids <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]', 'PERM_NON.R4[H/L]')
#' sample_ids %>% infer_design_from_sampleids()
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon$sample_id %>%
#'    autonomics.import::infer_design_from_sampleids()
#' }
#' @importFrom magrittr %>%
#' @export
infer_design_from_sampleids <- function(sample_ids){

   # Rm label tags from Max Quant ratios
   maxquant_ratios <- sample_ids %>% stringi::stri_detect_fixed('[H/L]') %>% any()
   if (maxquant_ratios) sample_ids %<>% stringi::stri_replace_first_fixed('[H/L]', '') %>%
                                       stringi::stri_replace_first_fixed('[H/M]', '') %>%
                                       stringi::stri_replace_first_fixed('[M/L]', '')
   # Infer sep
   possible_separators <- if (maxquant_ratios)  c('.', ' ')  else  c('.', ' ', '_')
   sep <- sample_ids %>% infer_design_sep(possible_separators)

   # Return dataframe with only sample ids if no separator could be infered
   if (is.null(sep))   return(data.frame(sample_id = sample_ids,
                                         subgroup  = '',
                                         replicate = '',
                                         block     = '',
                                         row.names = sample_ids))
   # Separator infered
  return(data.frame(sample_id = sample_ids,
                    subgroup  = sample_ids %>% stringi::stri_split_fixed(sep) %>% vapply(function(y) y %>% magrittr::extract(1:(length(y)-1)) %>% paste0(collapse = sep), character(1)),
                    replicate = sample_ids %>% stringi::stri_split_fixed(sep) %>% vapply(function(y) y %>% magrittr::extract(length(y)),                                  character(1)),
                    block     = '',
                    row.names = sample_ids))
}

#' Infer and add design to sdata
#' @param sdata sample dataframe
#' @return augmented sample dataframe
#' @importFrom magrittr %>%
#' @export
infer_add_design <- function(sdata){
   . <- NULL
   sdata %>%
   autonomics.support::left_join_keeping_rownames(
      autonomics.import::infer_design_from_sampleids(
         .$sample_id),
         .,
         by = 'sample_id')
}

