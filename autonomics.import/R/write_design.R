#=======================================
# sampleid_varname & subgroup_varname
#=======================================

#' Get sampleid svar name
#' @param platform 'metabolonlipids', 'metabolon', 'soma'
#' @return string
#' @examples
#' sampleid_varname('metabolonlipids')
#' sampleid_varname('metabolon')
#' sampleid_varname('soma')
#' @export
sampleid_varname <- function(platform){
   switch(platform,
          maxquant        = 'sample_id',
          metabolonlipids = 'Client Identifier',
          metabolon       = 'CLIENT_IDENTIFIER',
          soma            = 'SampleId')
}

#' Get subgroup svar name
#' @param platform 'metabolonlipids', 'metabolon', 'soma'
#' @return string
#' @examples
#' subgroup_varname('metabolonlipids')
#' subgroup_varname('metabolon')
#' subgroup_varname('soma')
#' @export
subgroup_varname <- function(platform){
   switch(platform,
          maxquant        = NULL,
          metabolonlipids = 'Group',
          metabolon       = 'Group',
          soma            = 'SampleGroup')
}

#==================================================
# INFER DESIGN FROM SAMPLEIDS
#==================================================

#' Infer design separator
#' @param sample_ids character vector with sampleids
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
      if (verbose)  autonomics.support::cmessage('%s: no (consistent) separator. Returning NULL', sample_ids[1])
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
                                                  sample_ids[1],
                                                  paste0(sprintf("'%s'", best_sep), collapse = ' or '))
      return(NULL)  # ambiguous separator
   }

   # Separator identified - return
   return(best_sep)
}


#' Infer design from (non-maxquant) sample ids
#' @param sample_ids  character vector
#' @param sep         string: design separator
#' @examples
#' require(magrittr)
#' sample_ids <- c("UT_10h_R1", "UT_10h_R2", "UT_10h_R3", "UT_10h_R4")
#' sample_ids %>% autonomics.import::infer_design_from_default_sampleids()
#' @importFrom magrittr %>%
#' @export
infer_design_from_default_sampleids <- function(
   sample_ids,
   sep = sample_ids %>% autonomics.import::infer_design_sep(c('.', ' ', '_'))
){

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


#' Infer design from maxquant sample ids
#' @param sample_ids  character vector
#' @param sep         string: design separator
#' @examples
#' require(magrittr)
#' sample_ids <- c("E(L).EM(M).BM(H).R1[M/L]", "E(L).EM(M).BM(H).R1[H/L]")
#' sample_ids %>% infer_design_from_maxquant_sampleids()
#' @importFrom magrittr %>%
#' @export
infer_design_from_maxquant_sampleids <- function(
   sample_ids,
   sep = sample_ids %>% autonomics.import::infer_design_sep(c('.', ' ', '_'))
){
   sample_ids %>%
   autonomics.import::designify_maxquant_sampleids(sep = sep) %>%
   autonomics.import::infer_design_from_sampleids(sep = sep) %>%
   (function(x){x$sample_id <- rownames(x) <- sample_ids; x})
}


#' Infer design from sample ids
#' @param sample_ids  character vector
#' @param sep         string: design separator
#' @param maxquant    logical: maxquant sample ids?
#' @examples
#' require(magrittr)
#' sample_ids <- c("UT_10h_R1", "UT_10h_R2", "UT_10h_R3", "UT_10h_R4")
#' sample_ids %>% autonomics.import::infer_design_from_sampleids()
#'
#' sample_ids <- c("E(L).EM(M).BM(H).R1[M/L]", "E(L).EM(M).BM(H).R1[H/L]")
#' sample_ids %>% infer_design_from_sampleids(maxquant=TRUE)
#' @importFrom magrittr %>%
#' @export
infer_design_from_sampleids <- function(
   sample_ids,
   sep = sample_ids %>% autonomics.import::infer_design_sep(c('.', ' ', '_')),
   maxquant = FALSE
){
   if (maxquant){ sample_ids %>% autonomics.import::infer_design_from_maxquant_sampleids(sep = sep)
   } else {       sample_ids %>% autonomics.import::infer_design_from_default_sampleids( sep = sep)
   }
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

#' Write design for metabolon, metabolonlipids, or soma data
#' @param file                   string: path to metabolon file
#' @param platform               'metabolon', 'metabolonlipids', or 'soma'
#' @param infer_design_from_sampleids   logical: whether to infer design from CLIENT_IDENTIFIER
#' @param design_sep             design separator
#' @param design_file            string: path to design file
#' @param sheet                  character(1) or numeric(1): xls sheet
#' @param quantity               either NULL, or any of: 'Ratio', 'Ratio normalized', 'Intensities', 'LFQ intensities', 'Reporter intensities'
#' @return design dataframe
#' @examples
#' require(magrittr)
#'
#' # MAXQUANT
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% write_design('maxquant')
#'    file %>% write_design('maxquant', infer_design_from_sampleids = TRUE)
#' }
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                                   package = 'autonomics.data')
#'    file %>% write_design('metabolon') %>% head()
#'    file %>% write_design('metabolon') %>% head()
#'    file %>% write_design('metabolon', infer_design_from_sampleids = TRUE) %>% head()
#' }
#'
#' # METABOLONLIPIDS
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% write_design('metabolonlipids') %>% head(3)
#'    file %>% write_design('metabolonlipids', infer_design_from_sampleids = TRUE) %>% head(3)
#' }
#'
#' # SOMA
#' if (require(autonomics.data)){
#'    soma_file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                              package = 'autonomics.data')
#'    soma_file %>% autonomics.import::write_design('soma')
#'    soma_file %>% autonomics.import::write_design('soma', infer_design_from_sampleids = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
write_design <- function(
   file,
   platform,
   infer_design_from_sampleids = FALSE,
   design_sep = NULL,
   design_file  = NULL,
   sheet = 2,
   quantity = if (platform=='maxquant') autonomics.import::infer_maxquant_quantity(file) else NULL
){

   # sdata
   sdata1 <- autonomics.import::load_sdata(file                        = file,
                                           platform                    = platform,
                                           sheet                       = sheet,
                                           quantity                    = quantity)
   sampleid_var <- autonomics.import::sampleid_varname(platform)
   subgroup_var <- autonomics.import::subgroup_varname(platform)

   # Construct design
   design_df <- data.frame(x           = sdata1[[sampleid_var]],
                           row.names   = sdata1[[sampleid_var]],
                           check.names = FALSE) %>%
                magrittr::set_names(names(.) %>% stringi::stri_replace_first_fixed('x', sampleid_var))

   # Infer subgroup from sampleids and add to design
   if (infer_design_from_sampleids){
      if (is.null(design_sep)) design_sep <- design_df[[sampleid_var]] %>% autonomics.import::infer_design_sep()
      inferred_design <- design_df[[sampleid_var]] %>% autonomics.import::infer_design_from_sampleids(sep = design_sep, maxquant = platform=='maxquant')
      design_df$sample_id <- NULL
      design_df %<>% cbind(., inferred_design)

   # Merge in design file
   } else {
      design_df$sample_id <- sdata1[[sampleid_var]]
      if (!is.null(subgroup_var)){
         design_df$subgroup  <- sdata1[[subgroup_var]]
         missing_subgroups <- any(autonomics.support::is_missing_or_empty_character(design_df$subgroup))
         if (!missing_subgroups){
            design_df$subgroup %<>% autonomics.support::nameify_strings()
            design_df %<>% add_replicate_values()
         }
      }
   }

   # Write to file
   if (!is.null(design_file))  design_df %>% autonomics.import::write_design_file(design_file)

   # Return
   design_df

}


#' Write maxquant design
#'
#' For automated infer from sampleids, use the following sample naming
#' scheme (before running max quant): WT(L).KD(M).OE(H).R1
#'
#' @param file string: full path to protein groups file
#' @param design_file        string: full path to sample design file
#' @param value_type         string: any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @param infer_design_from_sampleids       logical: should design be infered from sample ids (see details)
#' @param ... backward compatibility to deprecated functions
#' @return sample design dataframe
#' @examples
#' require(magrittr)
#'
#' # LABELED RATIOS and INTENSITIES
#' if (require(autonomics.data)){
#'    file <- system.file(
#'       'extdata/stemcell.differentiation/maxquant/proteinGroups.txt',
#'        package = 'autonomics.data')
#'    file %>% write_maxquant_design()
#'    file %>% write_maxquant_design(infer_design_from_sampleids = TRUE)
#'    file %>% write_maxquant_design(infer_design_from_sampleids = TRUE, value_type = 'raw.intensity')
#'    file %>% write_maxquant_design(infer_design_from_sampleids = TRUE, value_type = 'raw.ratio')
#' }
#'
#' # LFQ INTENSITIES
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::write_maxquant_design()
#'    file %>% autonomics.import::write_maxquant_design(infer_design_from_sampleids = TRUE)
#' }
#'
#' # UNLABELED INTENSITIES
#' if (require(alnoubi.2017)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'alnoubi.2017')
#'    file %>% autonomics.import::write_maxquant_design()
#'    file %>% autonomics.import::write_maxquant_design(infer_design_from_sampleids = TRUE)
#' }
#'
#' # REPORTER INTENSITIES
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::write_maxquant_design()
#'    file %>% autonomics.import::write_maxquant_design(infer_design_from_sampleids = TRUE)
#' }
#'
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
write_maxquant_design <- function(
   file,
   value_type   = autonomics.import::infer_maxquant_value_type(file),
   infer_design_from_sampleids = FALSE,
   design_file = NULL
){

   # Assert
   assertive.files::assert_all_are_dirs(dirname(file))
   #assertive.files::assert_all_are_readable_files(file, warn_about_windows = FALSE)

   # Read
   DT <- autonomics.support::cfread(file)

   # Extract
   sampleids <- DT %>% autonomics.import::get_maxquant_value_columns(value_type) %>% names()

   # Infer design and return
   design <- if (infer_design_from_sampleids){
      sampleids %>%
         autonomics.import::designify_maxquant_sampleids()   %>%
         autonomics.import::infer_design_from_sampleids('.') %>%
         (function(x){x$sample_id <- sampleids; x})
   } else {
      data.frame(sample_id = sampleids,
                 subgroup  = '',
                 replicate = '',
                 block     = '')
   } %>% magrittr::set_rownames(sampleids)


   # Sample file
   if (!is.null(design_file))  design %>% autonomics.import::write_design_file(design_file)
   design

}


#===================
# read_design
#===================


#' Read design
#' @param design_file string: path to design file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                 package = 'autonomics.data')
#'    design_file <- tempfile()
#'    autonomics.import::write_design(
#'       file, 'metabolon', infer_design_from_sampleids = TRUE, design_file = design_file)
#'    design_file %>% read_design() %>% head()
#' }
#'
#' # METABOLONLIPIDS
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    design_file <- tempfile()
#'    autonomics.import::write_design(
#'       file, 'metabolonlipids', design_file = design_file) %>% head()
#'    read_design(design_file) %>% head()
#' }
#'
#' # SOMA
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                 package = 'autonomics.data')
#'    design_file <- tempfile()
#'    write_design(file, 'soma', infer_design_from_sampleids = TRUE, design_file = design_file)
#'    read_design(design_file) %>% head()
#' }
#' @importFrom magrittr %<>%
#' @export
read_design <- function(design_file){
   design_df <- autonomics.support::cfread(design_file, data.table = FALSE)
   design_df %<>% magrittr::set_rownames(.[[1]])
   design_df
}










#==========================
# MAXQUANT
#==========================






#'@rdname write_maxquant_design
#'@export
create_maxquant_design_df <- function(...){
   .Deprecated('write_maxquant_design')
   write_maxquant_design(...)
}

#'@rdname write_maxquant_design
#'@export
create_sample_design_df <- function(...){
   .Deprecated('write_maxquant_design')
   write_maxquant_design(...)
}

#' @rdname write_maxquant_design
#' @export
create_maxquant_design_file <- function(...){
   .Deprecated('write_maxquant_design')
   write_maxquant_design(...)
}

#' @rdname write_maxquant_design
#' @export
create_sample_design_file <- function(...){
   .Deprecated('write_maxquant_design')
   write_maxquant_design(...)
}


#' Read maxquant design file
#' @param design_file sample design file
#' @return sample design datatable
#' @examples
#' if (require(autonomics.data)){
#'    design_file <- system.file('extdata/billing2016/sample_design.txt',
#'                                       package = 'autonomics.data')
#'    autonomics.import::read_maxquant_design(design_file)
#' }
#' if (require(billing.differentiation.data)){
#'    design_file <- system.file('extdata/maxquant/sample_design.txt',
#'                                       package = 'billing.differentiation.data')
#'    read_maxquant_design(design_file)
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
read_maxquant_design <- function(design_file){
   # Prevent CHECK warnings
   n <- NULL

   # Load
   #assertive.files::assert_all_are_readable_files(design_file, warn_about_windows = FALSE)
   sample_design <- autonomics.support::cfread(design_file, data.table = FALSE)

   # Check contents
   assertive.strings::assert_all_are_non_empty_character(sample_design$sample_id)
   assertive.strings::assert_any_are_non_missing_nor_empty_character(sample_design$subgroup)

   # Remove variable block if all values empty
   if ('block' %in% names(sample_design)){
      if (autonomics.support::all_are_missing_or_empty_character(sample_design$block))   sample_design$block <- NULL
   }

   # Create replicates if all absent.
   if (autonomics.support::all_are_missing_or_empty_character(sample_design$replicate)){
      sample_design %<>% dplyr::group_by_('subgroup')                   %>%
         dplyr::mutate(replicate = as.character(1:n())) %<>%
         dplyr::ungroup()                               %>%
         as.data.frame()
   }

   # Return
   sample_design %<>% data.table::data.table()
   return(sample_design)
}


#========================
# EXIQON
#========================


#' Write exiqon design
#' @param file            string: path to exiqon file
#' @param infer_design_from_sampleids   logical: whether to infer design from sampleids
#' @param design_file            string: path to sample file
#' @return string: path to design file
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                                package = 'subramanian.2016')
#'    file %>% autonomics.import::write_exiqon_design()
#'    file %>% autonomics.import::write_exiqon_design(infer_design_from_sampleids = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
write_exiqon_design <- function(
   file,
   infer_design_from_sampleids = FALSE,
   design_file = NULL
){
   sampleids <- autonomics.import::load_exiqon_sdata(file)  %>%
                magrittr::extract2('sample_id')
   design_df <- if (infer_design_from_sampleids){
                   sampleids %>% autonomics.import::infer_design_from_sampleids()
                } else {
                   return(data.frame(sample_id = sampleids,
                                     subgroup  = '',
                                     replicate = '',
                                     block     = ''))
                }
   if (!is.null(design_file)) design_df %>% autonomics.import::write_design_file(design_file)
   return(design_df)
}




#====================
# OBSOLETE FUNCTIONS
#====================

# Standardize design in sdata
#
# @param object        Summarizedexperiment
# @param sampleid_var  sampleid svar
# @param subgroup_var  subgroup svar
# @param infer_design  logical
# @return SummarizedExperiment
# @importFrom magrittr %>%
# @export
# prepare_design <- function(
#    object,
#    sampleid_var,
#    subgroup_var = NULL,
#    infer_design = if (is.null(subgroup_var)) TRUE else FALSE
# ){
#
#    design <- NULL
#
#    # Assert
#    assertive.sets::assert_is_subset(sampleid_var, autonomics.import::svars(object))
#    if (!is.null(subgroup_var)){
#       assertive.sets::assert_is_subset(subgroup_var, autonomics.import::svars(object))
#    }
#
#    # Either infer design from sampleids (if possible)
#    if (infer_design){
#       autonomics.support::cmessage('Infer design from sampleids')
#       sep <- autonomics.import::infer_design_sep(autonomics.import::sdata(object)[[sampleid_var]])
#       if (is.null(sep)){
#          autonomics.support::cmessage('No consistent unambiguous separator - design could not be prepared')
#       } else
#          design <- autonomics.import::sdata(object)[[sampleid_var]] %>%
#          autonomics.import::infer_design_from_sampleids()
#    }
#
#    # Or extract from sdata
#    if (!is.null(subgroup_var)){
#       design <- data.table::data.table(sample_id = autonomics.import::sdata(object)[[sampleid_var]],
#                                        subgroup  = autonomics.import::sdata(object)[[subgroup_var]])  %>%
#          magrittr::extract(, replicate := sprintf('R%d', 1:.N), by = 'subgroup')               %>%
#          data.frame(row.names = .$sample_id, stringsAsFactors = FALSE)
#    }
#
#    # Add to sdata in standardized format
#    if (!is.null(design)){
#       autonomics.import::sdata(object) %<>% cbind(design, .)                       %>%
#          autonomics.support::dedupe_varnames()  %>%
#          magrittr::set_rownames(.$sample_id)
#    }
#
#    # Return
#    object
# }


