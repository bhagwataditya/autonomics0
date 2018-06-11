
#' Write sample df to file
#' @param sample_df sample dataframe
#' @param sample_file sample file
#' @return path to sample file
#' @importFrom magrittr %>%
#' @export
write_sample_file <- function(sample_df, sample_file){

   # Abort
   if (file.exists(sample_file)) {
      autonomics.support::cmessage('\tAbort - file already exists: %s', sample_file)
      return(invisible(sample_file))
   }

   # Write
   sample_df %>% autonomics.support::print2txt(sample_file)
   autonomics.support::cmessage('\tWriting to: %s', sample_file)
   autonomics.support::cmessage('\tOpen this file in a Excel or LibrOffice, complete it manually and save.')

   # Open
   if(interactive() && assertive.reflection::is_windows()){
      tryCatch(shell.exec(sample_file), return(invisible(sample_file)))
   }

   # Return
   return(invisible(sample_file))
}



is_missing_or_empty_character <- function(x){
   x == '' | is.na(x)
}

#' @importFrom magrittr %>%
is_neither_missing_nor_empty_character <- function(x){
   x %>% is_missing_or_empty_character() %>% magrittr::not()
}

#' @importFrom magrittr %>%
all_are_missing_or_empty_character <- function(x){
   x %>% is_missing_or_empty_character() %>% all()
}

#' Load sample design
#' @param sample_file sample design file
#' @return sample design datatable
#' @examples
#' if (require(autonomics.data)){
#'    sample_file <- system.file('extdata/billing2016/sample_design.txt',
#'                                       package = 'autonomics.data')
#'    autonomics.import::load_sample_file(sample_file)
#' }
#' if (require(billing.differentiation.data)){
#'    sample_file <- system.file('extdata/maxquant/sample_design.txt',
#'                                       package = 'billing.differentiation.data')
#'    load_sample_file(sample_file)
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
load_sample_file <- function(sample_file){
   # Prevent CHECK warnings
   n <- NULL

   # Load
   #assertive.files::assert_all_are_readable_files(sample_file, warn_about_windows = FALSE)
   sample_design <- autonomics.support::cfread(sample_file,
                                               colClasses = c(sample_id = 'character',
                                                              subgroup  = 'character',
                                                              replicate = 'character',
                                                              block     = 'character'),
                                               data.table = FALSE)
   # Check contents
   assertive.strings::assert_all_are_non_empty_character(sample_design$sample_id)
   assertive.strings::assert_any_are_non_missing_nor_empty_character(sample_design$subgroup)

   # Remove variable block if all values empty
   if (all_are_missing_or_empty_character(sample_design$block))   sample_design$block <- NULL

   # Create replicates if all absent.
   if (all_are_missing_or_empty_character(sample_design$replicate)){
      sample_design %<>% dplyr::group_by_('subgroup')                   %>%
         dplyr::mutate(replicate = as.character(1:n())) %<>%
         dplyr::ungroup()                               %>%
         as.data.frame()
   }

   # Return
   sample_design %<>% data.table::data.table()
   return(sample_design)
}

load_maxquant_design <- function(...){
   .Deprecated('load_sample_file')
   autonomics.import::load_sample_file(...)
}


#' Nameify strings
#' @param x character vector
#' @param verbose character
#' @return validified subgroups
#' @examples
#' require(magrittr)
#' @importFrom magrittr   %>%    %<>%
#' @export
nameify_strings <- function(x, verbose = TRUE){

   # Satisfy CHECK
   . <- NULL

   old_values <- x
   new_values <- old_values %>% make.names()

   old_values %<>% setdiff(new_values)
   new_values %<>% setdiff(old_values)

   if (length(old_values) > 0){
      msg <- ''
      for (i in seq_along(old_values)){
         x   %<>% stringi::stri_replace_first_fixed(old_values[i], new_values[i])
         msg %<>%    paste0(sprintf('\n\t\t%s -> %s', old_values[i], new_values[i]))
      }
      msg %<>% substr(3, nchar(.))
      if (verbose) message(msg)
   }
   return(x)

}

#' Add sdata
#' @param object SummarizedExperiment
#' @param sample_file sample file
#' @return sumexp with sdata added
#' @importFrom magrittr  %>%   %<>%
#' @export
add_sdata <- function(object, sample_file){

   # Prevent CHECK notes
   sample_id <- NULL

   # Load
   sample_design <- load_sample_file(sample_file)
   sample_design %<>% as.data.frame()
   rownames(sample_design) <- sample_design$sample_id

   # Match and merge
   assertive.sets::assert_are_set_equal(autonomics.import::snames(object), sample_design$sample_id)
   idx <- match(autonomics.import::snames(object), sample_design$sample_id)
   sample_design %<>% magrittr::extract(idx, )
   autonomics.import::sdata(object) <- sample_design

   # Restrict to samples with subgroup annotation
   idx <- !is.null(object$subgroup) & object$subgroup!=''
   assertive.base::assert_any_are_true(idx)
   if (any(!idx)){
      autonomics.support::cmessage('\trm samples with missing subgroup annotation: \n\t\t%s',
                                   autonomics.import::snames(object)[!idx] %>% paste0(collapse = '\n\t\t'))
   }
   object %<>% magrittr::extract(, idx)

   # Validify subgroups
   message('\tValidify subgroups')
   autonomics.import::sdata(object)$subgroup  %<>%  nameify_strings(verbose = TRUE)

   # Intuify sample ids
   subgroups  <- object$subgroup %>% as.character()
   replicates <- object$replicate %>% as.character()
   new_sample_ids <- paste0(subgroups, '.R', replicates)
   new_ids_suited <- all(assertive.strings::is_non_empty_character(subgroups))  &
      all(assertive.strings::is_non_empty_character(replicates)) &
      assertive.properties::has_no_duplicates(new_sample_ids)
   if (new_ids_suited){
      autonomics.import::snames(object) <- new_sample_ids
      autonomics.import::sdata(object)$sample_id <- new_sample_ids
   }

   # Return
   return(object)
}

add_maxquant_sdata <- function(...){
   .Deprecated('add_sdata')
   autonomics.import::add_sdata(...)
}
