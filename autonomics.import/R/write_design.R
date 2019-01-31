#=======================================
# sampleid_varname & subgroup_varname
#=======================================

#' Get sampleid svar name
#' @param platform  'exiqon', 'maxquant', 'metabolonlipids', 'metabolon', 'soma'
#' @return string
#' @examples
#' sampleid_varname('metabolonlipids')
#' sampleid_varname('metabolon')
#' sampleid_varname('soma')
#' @export
sampleid_varname <- function(platform){
   switch(platform,
          exiqon          = 'sample_id',
          maxquant        = 'sample_id',
          metabolonlipids = 'Client Identifier',
          metabolon       = 'sample_id',
          soma            = 'sample_id',
          rnaseq          = 'sample_id')
}

#' Get subgroup svar name
#' @param platform   'exiqon', 'maxquant', 'metabolonlipids', 'metabolon', 'soma'
#' @return string
#' @examples
#' subgroup_varname('metabolonlipids')
#' subgroup_varname('metabolon')
#' subgroup_varname('soma')
#' @export
subgroup_varname <- function(platform){
   switch(platform,
          exiqon          = NULL,
          maxquant        = NULL,
          metabolonlipids = 'Group',
          metabolon       = 'subgroup',
          soma            = 'subgroup',
          rnaseq          = NULL)
}

#========================================
# WRITE DESIGN FILE
#========================================

#' Write sample df to file
#' @param design_df sample dataframe
#' @param design_file sample file
#' @return path to sample file
#' @importFrom magrittr %>%
#' @export
write_design_file <- function(design_df, design_file){

   # Abort
   if (file.exists(design_file)) {
      autonomics.support::cmessage('\tAbort - file already exists: %s', design_file)
      return(invisible(design_file))
   }

   # Write
   design_df %>% autonomics.support::print2txt(design_file)
   autonomics.support::cmessage('\tWriting to: %s', design_file)
   autonomics.support::cmessage('\tOpen this file in a Excel or LibrOffice, complete it manually and save.')

   # Open
   if(interactive() && assertive.reflection::is_windows()){
      tryCatch(shell.exec(design_file), return(invisible(design_file)))
   }

   # Return
   return(invisible(design_file))
}

#=====================
# ADD REPLICATE VALUES
#=======================

#' Add replicate values
#'
#' Replicate values are set to the unique sampleid tails within each subgroup
#'
#' @param design_df dataframe with columns 'sample_id' and 'subgroup'
#' @return dataframe
#' @examples
#' require(magrittr)
#' design_df <- data.frame(sample_id = c('E_1', 'E_2', 'E_3', 'EM_1', 'EM_2', 'EM_3'),
#'                         subgroup  = c('E',   'E',   'E',   'EM',   'EM',   'EM')) %>%
#'              magrittr::set_rownames(.$sample_id)
#' design_df
#' design_df %>% add_replicate_values()
#' @importFrom magrittr %>%
#' @export
add_replicate_values <- function(design_df){
   sample_id <- NULL
   design_df %>% data.table::data.table() %>%
      magrittr::extract(, replicate := autonomics.support::get_unique_tails(sample_id), by = 'subgroup') %>%
      data.frame(row.names = rownames(design_df), check.names = FALSE)
}

#========================================
# WRITE DESIGN
#========================================

#' Write design
#' @param object        SummarizedExperiment
#' @param subgroup_var  character(1)
#' @param sep           design separator
#' @param file          character(1)
#' @return design dataframe
#' @examples
#' require(magrittr)
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'    object <- 'extdata/glutaminase/glutaminase.xlsx'     %>%
#'               system.file(package = 'autonomics.data')  %>%
#'               autonomics.import::read_metabolon()
#'    object %>% autonomics.import::write_design() %>% head()
#'    object %>% autonomics.import::write_design(subgroup_var = "Group   HMDB_ID") %>% head()
#' }
#'
#' # SOMASCAN
#' if (require(autonomics.data)){
#'    object <- 'extdata/stemcomp/soma/stemcomp.adat'     %>%
#'               system.file(package = 'autonomics.data') %>%
#'               autonomics.import::read_somascan()
#'    object %>% autonomics.import::write_design()
#'    object %>% autonomics.import::write_design(subgroup_var = 'SampleGroup')
#' }
#'
#' # EXIQON
#' if (require(subramanian.2016)){
#'    object <- 'extdata/exiqon/subramanian.2016.exiqon.xlsx'  %>%
#'               system.file(package = 'subramanian.2016')     %>%
#'               autonomics.import::read_exiqon()
#'    object %>% write_design() %>% head()
#' }
#'
#' # PROTEINGROUPS
#' if (require(autonomics.data)){
#'    object <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'               system.file(package = 'autonomics.data')     %>%
#'               autonomics.import::read_proteingroups(simplify_snames = TRUE)
#'    object %>% write_design()
#' }
#'
#' @importFrom magrittr %>%
#' @export
write_design <- function(
   object,
   subgroup_var = NULL,
   sep          = object %>% autonomics.import::guess_sep(verbose = TRUE),
   file         = NULL
){

   # Create design dataframe
   sampleid_values <- object %>% autonomics.import::sampleid_values()
   subgroup_values <- object %>% autonomics.import::guess_subgroup_values(subgroup_var, verbose = TRUE)
   design_df <- data.frame(sample_id        = sampleid_values,
                           subgroup         = subgroup_values,
                           row.names        = sampleid_values,
                           check.names      = FALSE,
                           stringsAsFactors = FALSE)

   # Write to file
   if (!is.null(file))  design_df %>% autonomics.import::write_design_file(file)

   # Return
   design_df

}

