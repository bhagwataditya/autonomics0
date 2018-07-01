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
#' @param sample_ids         character vector with sampleids
#' @param sep                string separator
#' @param drop_common_parts  logical
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
infer_design_from_sampleids <- function(
   sample_ids,
   sep = sample_ids %>% autonomics.import::infer_design_sep(c('.', ' ', '_')),
   drop_common_parts = FALSE
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
   design_df %>% data.table::data.table() %>%
                 magrittr::extract(, replicate := autonomics.support::get_unique_tails(sample_id), by = 'subgroup') %>%
                 data.frame(row.names = rownames(design_df))
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


#==========================
# MAXQUANT
#==========================

#' Designify maxquant sampleids
#' @param sampleids character vector
#' @examples
#' sampleids <- c("STD(L).EM00(M).EM01(H).R1[H/L]",
#'                "STD(L).EM00(M).EM01(H).R1[M/L]")
#' sampleids %>% autonomics.import::designify_maxquant_sampleids()
#'
#' sampleids <- c("Gel 11 1", "Gel 11 2", "Gel Ctrl 1")
#' sampleids %>% autonomics.import::designify_maxquant_sampleids()
#'
#' sampleids <- c("ESC(0).NCM(1).CM(2).MV(3).EX(4).KIT(5).R1[0]",
#'                "ESC(0).NCM(1).CM(2).MV(3).EX(4).KIT(5).R1[1]")
#' sampleids %>% autonomics.import::designify_maxquant_sampleids()
#' @importFrom magrittr %>%
#' @export
designify_maxquant_sampleids <- function(sampleids){

   # Break into parts         #pattern <- '(.*)\\.R\\((.+)\\)\\[(.+)\\]'
   sep <- sampleids %>% autonomics.import::infer_design_sep()
   parts <- strsplit(sampleids, split = sep, fixed = TRUE)
   n <- length(parts[[1]])

   # replicate values
   replicate_values <- parts %>% vapply(extract, character(1), n)
   is_labeled <- parts %>% vapply(extract, character(1), n) %>% stringi::stri_detect_fixed('[') %>% any()
   if (is_labeled){
      label_values <- parts %>% vapply(extract, character(1), n) %>%
         stringi::stri_extract_first_regex('(?<=\\[)(.+)(?=\\])')  %>%
         strsplit(split = '/', fixed = TRUE)
      replicate_values %<>% stringi::stri_replace_first_regex('(.+)(\\[.+\\])', '$1')
      replicate_values %<>% paste0('_', label_values %>% vapply(paste0, character(1), collapse = ''))
   }
   parts %<>% lapply(function(x)x %>% magrittr::extract(1:(length(x)-1)))

   # subgroup values
   subgroup_values <- if (is_labeled){
      parts %>% lapply(function(x){
         cursample <- x %>% stringi::stri_replace_first_regex('(.+)\\(([HML0-9]+)\\)', '$1')
         curname   <- x %>% stringi::stri_replace_first_regex('(.+)\\(([HML0-9]+)\\)', '$2')
         cursample %>% magrittr::set_names(curname)
      }) %>%
         autonomics.support::vextract(label_values) %>%
         vapply(paste0, character(1), collapse = '_')
   } else {
      parts %>% vapply(paste0, character(1), collapse = '.')
   }

   # sampleid values
   long_sampleid_values <- sprintf('%s.%s', subgroup_values, replicate_values)
   sampleid_values <- long_sampleid_values %>% stringi::stri_replace_first_regex('_[HML0-9]+', '')
   idx <- sampleid_values %>% autonomics.support::cduplicated()
   sampleid_values[idx] <- long_sampleid_values[idx]
   sampleid_values
}


#' Write maxquant design
#'
#' For automated infer from sampleids, use the following sample naming
#' scheme (before running max quant): WT(L).KD(M).OE(H).R1
#'
#' @param file string: full path to protein groups file
#' @param design_file        string: full path to sample design file
#' @param value_type         string: any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @param infer_from_sampleids       logical: should design be infered from sample ids (see details)
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
#'    file %>% write_maxquant_design(infer_from_sampleids = TRUE)
#'    file %>% write_maxquant_design(infer_from_sampleids = TRUE, value_type = 'raw.intensity')
#'    file %>% write_maxquant_design(infer_from_sampleids = TRUE, value_type = 'raw.ratio')
#' }
#'
#' # LFQ INTENSITIES
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::write_maxquant_design()
#'    file %>% autonomics.import::write_maxquant_design(infer_from_sampleids = TRUE)
#' }
#'
#' # UNLABELED INTENSITIES
#' if (require(alnoubi.2017)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'alnoubi.2017')
#'    file %>% autonomics.import::write_maxquant_design()
#'    file %>% autonomics.import::write_maxquant_design(infer_from_sampleids = TRUE)
#' }
#'
#' # REPORTER INTENSITIES
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::write_maxquant_design()
#'    file %>% autonomics.import::write_maxquant_design(infer_from_sampleids = TRUE)
#' }
#'
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
write_maxquant_design <- function(
   file,
   value_type   = autonomics.import::infer_maxquant_value_type(file),
   infer_from_sampleids = FALSE,
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
   design <- if (infer_from_sampleids){
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
#' @param infer_from_sampleids   logical: whether to infer design from sampleids
#' @param design_file            string: path to sample file
#' @return string: path to design file
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                                package = 'subramanian.2016')
#'    file %>% autonomics.import::write_exiqon_design()
#'    file %>% autonomics.import::write_exiqon_design(infer_from_sampleids = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
write_exiqon_design <- function(
   file,
   infer_from_sampleids = FALSE,
   design_file = NULL
){
   sampleids <- autonomics.import::load_exiqon_sdata(file)  %>%
                magrittr::extract2('sample_id')
   design_df <- if (infer_from_sampleids){
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


#===========================
# METABOLON
#===========================


#' Write metabolon design
#' @param  file                   string: path to metabolon file
#' @param  sheet                  character(1) or numeric(1): xls sheet
#' @param  infer_from_sampleids   logical: whether to infer design from CLIENT_IDENTIFIER
#' @param  design_file            string: path to design file
#' @return design dataframe
#' @examples
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                                   package = 'autonomics.data')
#'    file %>% autonomics.import::write_metabolon_design(sheet=2)
#'    file %>% autonomics.import::write_metabolon_design(infer_from_sampleids = TRUE, sheet=2)
#' }
#' @importFrom magrittr %>%
#' @export
write_metabolon_design <- function(
   file,
   sheet,
   infer_from_sampleids = FALSE,
   design_file  = NULL
){
   sdata1 <- file %>% autonomics.import::load_metabolon_sdata(sheet = sheet)
   design_df <- data.frame(CLIENT_IDENTIFIER = sdata1$CLIENT_IDENTIFIER,
                           row.names         = sdata1$CLIENT_IDENTIFIER)
   if (infer_from_sampleids){
      design_df %<>% cbind(autonomics.import::infer_design_from_sampleids(.$CLIENT_IDENTIFIER))
   } else {
      design_df$sample_id <- sdata1$CLIENT_IDENTIFIER
      design_df$subgroup  <- sdata1$Group
      missing_subgroups <- any(autonomics.support::is_missing_or_empty_character(design_df$subgroup))
      if (!missing_subgroups){
         design_df$subgroup %<>% autonomics.support::nameify_strings()
         design_df %<>% add_replicate_values()
      }
   }
   if (!is.null(design_file))  design_df %>% autonomics.import::write_design_file(design_file)
   design_df

}

#' Read metabolon design
#' @param design_file string: path to design file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    design_file <- tempfile()
#'    system.file('extdata/glutaminase/glutaminase.xlsx',
#'                 package = 'autonomics.data') %>%
#'    autonomics.import::write_metabolon_design(design_file, infer_from_sampleids = TRUE)
#'    design_file %>% read_metabolon_design()
#' }
#' @importFrom magrittr %<>%
#' @export
read_metabolon_design <- function(design_file){
   design_df <- autonomics.support::cfread(design_file, data.table = FALSE)
   assertive::assert_is_subset('CLIENT_IDENTIFIER', names(design_df))
   design_df %<>% magrittr::set_rownames(.$CLIENT_IDENTIFIER)
   design_df
}


#=======================
# SOMA
#=======================

#' Write soma design
#' @param soma_file             string: path to soma file
#' @param infer_from_sampleids  logical: whether to infer design from sampleids
#' @param design_file           string: path of sample file to which design is written
#' @return sample design dataframe
#' @examples
#' if (require(autonomics.data)){
#'    soma_file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                              package = 'autonomics.data')
#'    soma_file %>% autonomics.import::write_soma_design()
#'    soma_file %>% autonomics.import::write_soma_design(infer_from_sampleids = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
write_soma_design <- function(
   soma_file,
   infer_from_sampleids = FALSE,
   design_file = NULL
){

   sdata1 <- autonomics.import::load_soma_sdata(soma_file)
   design_df <- data.frame(SampleId  = sdata1$SampleId,
                           row.names = sdata1$SampleId)
   if (infer_from_sampleids){
      design_df %<>% cbind(autonomics.import::infer_design_from_sampleids(.$SampleId))
   } else {
      design_df$sample_id <- sdata1$SampleId
      design_df$subgroup  <- sdata1$SampleGroup
      missing_subgroups <- any(autonomics.support::is_missing_or_empty_character(design_df$subgroup))
      if (!missing_subgroups){
         design_df$subgroup %<>% autonomics.support::nameify_strings()
         design_df %<>% add_replicate_values()
      }
   }
   if (!is.null(design_file)) design_df %>% autonomics.import::write_design_file(design_file)
   design_df
}


#' Read soma design
#' @param design_file string: path to design file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    design_file <- tempfile()
#'    system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                 package = 'autonomics.data') %>%
#'    autonomics.import::write_soma_design(design_file, infer_from_sampleids = TRUE)
#'    design_file %>% read_soma_design()
#' }
#' @importFrom magrittr %<>%
#' @export
read_soma_design <- function(design_file){
   design_df <- autonomics.support::cfread(design_file, data.table = FALSE)
   assertive::assert_is_subset('SampleId', names(design_df))
   design_df %<>% magrittr::set_rownames(.$SampleId)
   design_df
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


