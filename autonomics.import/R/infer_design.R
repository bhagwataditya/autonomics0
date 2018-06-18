
#' Infer design separator
#' @param sample_ids character vector with sample ids
#' @param possible_separators character vector with possible separators to look for
#' @param verbose logical
#' @return separator (string) or NULL (if no separator could be identified)
#' @examples
#' require(magrittr)
#' sample_ids <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]', 'PERM_NON.R4[H/L]')
#' sample_ids %>% infer_design_sep()
#'
#' sample_ids <- c('WT untreated 1', 'WT untreated 2', 'WT treated 1')
#' sample_ids %>% infer_design_sep()
#'
#' sample_ids <- c('group1', 'group2', 'group3.R1')
#' sample_ids %>% infer_design_sep()
#' @importFrom magrittr %>%
#' @export
infer_design_sep <- function(sample_ids, possible_separators = c('.', ' ', '_'), verbose = TRUE){
   . <- NULL
   sep_freqs <- Map(function(x) stringi::stri_split_fixed(sample_ids, x), possible_separators)        %>%
                lapply(function(x) x %>% vapply(length, integer(1)))                                  %>%
                magrittr::extract( vapply(., autonomics.support::has_identical_values, logical(1)))   %>%
                vapply(unique, integer(1))

   # No separator detected - return NULL
   if (all(sep_freqs==1)){
      if (verbose)  autonomics.support::cmessage('No (consistent) separator. Returning NULL')
      return(NULL)   # no separator detected
   }

   # Find best separator
   best_sep <- sep_freqs %>%
               magrittr::extract(.!=1)  %>%
               magrittr::extract(autonomics.support::is_max(vapply(., magrittr::extract, integer(1), 1)))   %>%
               names()

   # Ambiguous separator - return NULL
   if (length(best_sep)>1){
      if (verbose)   autonomics.support::cmessage('No unambiguous separator (%s). Returning NULL',
                                                  paste0(sprintf("'%s'", best_sep), collapse = ' or '))
      return(NULL)  # ambiguous separator
   }

   # Separator identified - return
   return(best_sep)
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
#' sample_ids <- c("UT_10h_R1", "UT_10h_R2", "UT_10h_R3", "UT_10h_R4")
#' sample_ids %>% autonomics.import::infer_design_from_sampleids()
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
   sep <- sample_ids %>% autonomics.import::infer_design_sep(possible_separators)

   # Return dataframe with only sample ids if no separator could be infered
   if (is.null(sep))   return(data.frame(sample_id = sample_ids,
                                         subgroup  = '',
                                         replicate = '',
                                         row.names = sample_ids))

   # Extract subgroup and replicate
   subgroup_values  <- sample_ids %>% stringi::stri_split_fixed(sep) %>%
                                      vapply(function(y) y %>% magrittr::extract(1:(length(y)-1)) %>%
                                      paste0(collapse = sep), character(1))
   replicate_values <- sample_ids %>% stringi::stri_split_fixed(sep) %>%
                                      vapply(function(y) y %>% magrittr::extract(length(y)), character(1))

   # Return df
   return(data.frame(sample_id = sample_ids,
                     subgroup  = subgroup_values,
                     replicate = replicate_values,
                     row.names = sample_ids,
                     stringsAsFactors = FALSE))
}


#' Standardize design in sdata
#'
#' @param object        Summarizedexperiment
#' @param sampleid_var  sampleid svar
#' @param subgroup_var  subgroup svar
#' @param infer_design  logical
#' @return SummarizedExperiment
#' @importFrom magrittr %>%
#' @export
prepare_design <- function(
   object,
   sampleid_var,
   subgroup_var = NULL,
   infer_design = if (is.null(subgroup_var)) TRUE else FALSE
){

   design <- NULL

   # Assert
   assertive.sets::assert_is_subset(sampleid_var, autonomics.import::svars(object))
   if (!is.null(subgroup_var)){
      assertive.sets::assert_is_subset(subgroup_var, autonomics.import::svars(object))
   }

   # Either infer design from sampleids (if possible)
   if (infer_design){
      autonomics.support::cmessage('Infer design from sampleids')
      sep <- autonomics.import::infer_design_sep(autonomics.import::sdata(object)[[sampleid_var]])
      if (is.null(sep)){
         autonomics.support::cmessage('No consistent unambiguous separator - design could not be prepared')
      } else
         design <- autonomics.import::sdata(object)[[sampleid_var]] %>%
                   autonomics.import::infer_design_from_sampleids()
   }

   # Or extract from sdata
   if (!is.null(subgroup_var)){
      design <- data.table::data.table(sample_id = autonomics.import::sdata(object)[[sampleid_var]],
                                       subgroup  = autonomics.import::sdata(object)[[subgroup_var]])  %>%
                magrittr::extract(, replicate := sprintf('R%d', 1:.N), by = 'subgroup')               %>%
                data.frame(row.names = .$sample_id, stringsAsFactors = FALSE)
   }

   # Add to sdata in standardized format
   if (!is.null(design)){
      autonomics.import::sdata(object) %<>% cbind(design, .)                       %>%
                                            autonomics.support::dedupe_varnames()  %>%
                                            magrittr::set_rownames(.$sample_id)
   }

   # Return
   object
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

