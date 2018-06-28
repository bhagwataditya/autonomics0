#==================================================
# INFER DESIGN FROM SAMPLEIDS
#==================================================

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
#' @param sample_ids         character vector with sample ids
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
#' @param sample_df dataframe with columns 'sample_id' and 'subgroup'
#' @return
#' @examples
#' require(magrittr)
#' sample_df <- data.frame(sample_id = c('E_1', 'E_2', 'E_3', 'EM_1', 'EM_2', 'EM_3'),
#'                         subgroup  = c('E',   'E',   'E',   'EM',   'EM',   'EM')) %>%
#'              magrittr::set_rownames(.$sample_id)
#' sample_df
#' sample_df %>% add_replicate_values()
#' @importFrom magrittr %>%
#' @export
add_replicate_values <- function(sample_df){
   sample_df %>% data.table::data.table() %>%
                 magrittr::extract(, replicate := autonomics.support::get_unique_tails(sample_id), by = 'subgroup') %>%
                 data.frame(row.names = rownames(sample_df))
}

#========================================
# WRITE SAMPLE FILE
#========================================

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

#' Create maxquant sample df
#'
#' For automated infer of design, use the following sample naming
#' scheme (before running max quant): WT(L).KD(M).OE(H).R1
#'
#' @param proteingroups_file string: full path to protein groups file
#' @param sample_file        string: full path to sample design file
#' @param value_type         string: any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @param infer_design       logical: should design be infered from sample ids (see details)
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
#'    file %>% create_maxquant_sample_df()
#'    file %>% create_maxquant_sample_df(infer_design = TRUE)
#'    file %>% create_maxquant_sample_df(infer_design = TRUE, value_type = 'raw.intensity')
#'    file %>% create_maxquant_sample_df(infer_design = TRUE, value_type = 'raw.ratio')
#' }
#'
#' # LFQ INTENSITIES
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::create_maxquant_sample_df()
#'    file %>% autonomics.import::create_maxquant_sample_df(infer_design = TRUE)
#' }
#'
#' # UNLABELED INTENSITIES
#' if (require(alnoubi.2017)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'alnoubi.2017')
#'    file %>% autonomics.import::create_maxquant_sample_df()
#'    file %>% autonomics.import::create_maxquant_sample_df(infer_design = TRUE)
#' }
#'
#' # REPORTER INTENSITIES
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::create_maxquant_sample_df()
#'    file %>% autonomics.import::create_maxquant_sample_df(infer_design = TRUE)
#' }
#'
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
create_maxquant_sample_df <- function(
   proteingroups_file,
   value_type   = autonomics.import::infer_maxquant_value_type(proteingroups_file),
   infer_design = FALSE
){

   # Assert
   assertive.files::assert_all_are_dirs(dirname(proteingroups_file))
   #assertive.files::assert_all_are_readable_files(proteingroups_file, warn_about_windows = FALSE)

   # Read
   DT <- autonomics.support::cfread(proteingroups_file)

   # Extract
   sampleids <- DT %>% autonomics.import::get_maxquant_value_columns(value_type) %>% names()

   # Infer design and return
   if (infer_design){
      design <- sampleids %>%
         autonomics.import::designify_maxquant_sampleids()   %>%
         autonomics.import::infer_design_from_sampleids('.')
      design$sample_id <- sampleids
      rownames(design) <- sampleids
      return(design)
   }

   # Return design template if not able to infer
   return(data.frame(
      sample_id = sampleids,
      subgroup  = '',
      replicate = '',
      block     = ''))

}


#'@rdname create_maxquant_sample_df
#'@export
create_maxquant_design_df <- function(...){
   .Deprecated('create_maxquant_sample_df')
   create_maxquant_sample_df(...)
}

#'@rdname create_maxquant_sample_df
#'@export
create_sample_design_df <- function(...){
   .Deprecated('create_maxquant_design_df')
   create_maxquant_design_df(...)
}

#' Create maxquant sample file
#'
#' For automated infer of design, use the following sample naming
#' scheme (before running max quant): WT(L).KD(M).OE(H).R1
#'
#' @param proteingroups_file string: full path to protein groups file
#' @param sample_file        string: full path to sample design file
#' @param value_type         string: any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @param infer_design       logical: should design be infered from sample ids (see details)
#' @param ... backward compatibility to deprecated functions
#' @return sample design dataframe
#' @examples
#' require(magrittr)
#'
#' # LABELED RATIOS and INTENSITIES
#' if (require(autonomics.data)){
#'    file <- system.file(
#'              'extdata/stemcell.differentiation/maxquant/proteinGroups.txt',
#'               package = 'autonomics.data')
#'    file %>% create_maxquant_sample_df()
#'    file %>% create_maxquant_sample_file(sample_file = tempfile())
#' }
#'
#' @importFrom magrittr %>%
#' @export
create_maxquant_sample_file <- function(
   proteingroups_file,
   sample_file  = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   value_type   = autonomics.import::infer_maxquant_value_type(proteingroups_file),
   infer_design = FALSE
){

   # Create
   autonomics.import::create_maxquant_sample_df(proteingroups_file = proteingroups_file,
                                                value_type         = value_type,
                                                infer_design       = infer_design) %>%
      autonomics.import::write_sample_file(sample_file)

   # Return
   return(invisible(sample_file))
}


#' @rdname create_maxquant_sample_file
#' @export
create_maxquant_design_file <- function(...){
   create_maxquant_sample_file(...)
}

#' @rdname create_maxquant_sample_file
#' @export
create_sample_design_file <- function(...){
   .Deprecated('write_maxquant_design')
   write_maxquant_design(...)
}


#========================
# EXIQON
#========================


#' Create exiqon sample dataframe
#' @param exiqon_file  string
#' @param sample_file  string
#' @param value_type   string: any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @return string: path to design file
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    exiqon_file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                                package = 'subramanian.2016')
#'    exiqon_file %>% autonomics.import::create_exiqon_sample_df()
#'    exiqon_file %>% autonomics.import::create_exiqon_sample_df(infer_design = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
create_exiqon_sample_df <- function(
   exiqon_file,
   infer_design = FALSE
){
   sampleids <- autonomics.import::load_exiqon_sdata(exiqon_file)  %>%
                magrittr::extract2('sample_id')
   if (infer_design){
      sampleids %>% autonomics.import::infer_design_from_sampleids()
   } else {
      return(data.frame(sample_id = sampleids,
                        subgroup  = '',
                        replicate = '',
                        block     = ''))
   }
}

#' Create exiqon sample file
#' @param exiqon_file  string: path to exiqon file
#' @param sample_file  string: path to sample file
#' @param infer_design logical: whether to infer design from sampleids
#' @return string: path to design file
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016){
#'    exiqon_file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                                package = 'subramanian.2016')
#'    exiqon_file %>% autonomics.import::create_exiqon_sample_file(sample_file = tempfile())
#'
#' }
#' @importFrom magrittr %>%
#' @export
create_exiqon_sample_file <- function(
   exiqon_file,
   sample_file  = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   infer_design = FALSE
){
   autonomics.import::create_exiqon_sample_df(exiqon_file,
                                              infer_design = infer_design) %>%
      autonomics.import::write_sample_file(sample_file)
}


#===========================
# METABOLON
#===========================


#' Create metabolon sample dataframe
#' @param  metabolon_file  string: path to metabolon file
#' @param  sample_file     string: path to sample file
#' @param  infer_design    logical: whether to infer design from CLIENT_IDENTIFIER
#' @return design dataframe
#' @examples
#' if (require(autonomics.data)){
#'    metabolon_file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                                   package = 'autonomics.data')
#'    metabolon_file %>% autonomics.import::create_metabolon_sample_df()
#'    metabolon_file %>% autonomics.import::create_metabolon_sample_df(infer_design = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
create_metabolon_sample_df <- function(
   metabolon_file,
   infer_design = FALSE
){
   sdata1 <- metabolon_file %>% autonomics.import::load_metabolon_sdata()
   design_df <- data.frame(CLIENT_IDENTIFIER = sdata1$CLIENT_IDENTIFIER,
                           row.names         = sdata1$CLIENT_IDENTIFIER)
   if (infer_design){
      design_df %<>% cbind(autonomics.import::infer_design_from_sampleids(.$CLIENT_IDENTIFIER))
   } else {
      design_df$sample_id <- sdata1$CLIENT_IDENTIFIER
      design_df$subgroup  <- sdata1$Group
      design_df %<>% add_replicate_values()
   }
   design_df
}

#' Create metabolon sample file
#' @param metabolon_file  string: path to metabolon file
#' @param sample_file     string: path to sample file
#' @param infer_design    logical: whether to infer design from sampleids
#' @return path to metabolon file
#' @examples
#' if (require(autonomics.data)){
#'    metabolon_file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                                   package = 'autonomics.data')
#'    metabolon_file %>% autonomics.import::create_metabolon_sample_file(
#'                          sample_file = tempfile())
#' }
#' @importFrom magrittr %>%
#' @export
create_metabolon_sample_file <- function(
   metabolon_file,
   sample_file     = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   infer_design    = FALSE
){
   df <- metabolon_file %>%
         autonomics.import::create_metabolon_sample_df(infer_design = infer_design)
   df %>% autonomics.import::write_sample_file(sample_file)
   df
}


#=======================
# SOMA
#=======================


#' Create soma sample dataframe
#' @param soma_file    string: path to soma file
#' @param infer_design logical: whether to infer design from sampleids
#' @return sample design dataframe
#' @examples
#' if (require(autonomics.data)){
#'    soma_file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                              package = 'autonomics.data')
#'    soma_file %>% autonomics.import::create_soma_sample_df()
#'    soma_file %>% autonomics.import::create_soma_sample_df(infer_design = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
create_soma_sample_df <- function(
   soma_file,
   infer_design = FALSE
){

   sdata1 <- autonomics.import::load_soma_sdata(soma_file)
   design_df <- data.frame(SampleId  = sdata1$SampleId,
                           row.names = sdata1$SampleId)
   if (infer_design){
      design_df %<>% cbind(autonomics.import::infer_design_from_sampleids(.$SampleId))
   } else {
      design_df$sample_id <- sdata1$SampleId
      design_df$subgroup  <- sdata1$SampleGroup
      design_df %<>% add_replicate_values()
   }
   design_df
}



#' Create soma sample file
#' @param soma_file     string: path to soma file
#' @param sample_file   string: path to sample file
#' @param infer_design  logical: whether to infer design from sampleids
#' @return path to soma design file
#' @examples
#'    soma_file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                              package = 'autonomics.data')
#'    soma_file %>% autonomics.import::create_soma_sample_file(sample_file = tempfile())
#' @importFrom magrittr %>%
#' @export
create_soma_sample_file <- function(
   soma_file,
   sample_file,
   infer_design = FALSE
){
   df <- soma_file %>%
      autonomics.import::create_soma_sample_df(infer_design = infer_design)
   df %>% autonomics.import::write_sample_file(sample_file)
   df
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


