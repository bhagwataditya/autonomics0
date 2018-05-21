################################################################################
#                                                                              #
#                        proteingroups -> sample design                        #
#                                                                              #
################################################################################

#' Extract normalized ratio columns
#' @param DT maxquant data.table
#' @return normalized ratios data.table
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/proteinGroups.txt',
#'                 package = 'billing.differentiation.data') %>%
#'       data.table::fread() %>%
#'       get_maxquant_normalized_ratios()
#' }
#' @export
get_maxquant_normalized_ratios <- function(DT){

   # Satisfy CHECK
   . <- NULL

   # Identify
   pattern <- '^Ratio ([HM]/[ML]) normalized (.+)'
   colnames <- names(DT) %>%
      magrittr::extract(stringi::stri_detect_regex(., pattern)) %>%
      magrittr::extract(!stringi::stri_detect_fixed(., '___'))  # phosphosites %>%

   # Extract
   if (length(colnames)==0) return(NULL)
   values <- DT[, colnames, with = FALSE]
   names(values) %<>% stringi::stri_replace_all_regex(pattern, '$2[$1]')
   values
}

#' Extract raw ratio columns
#' @param DT maxquant data.table
#' @return raw ratio data.table
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/proteinGroups.txt',
#'                 package = 'billing.differentiation.data') %>%
#'       data.table::fread() %>%
#'       get_maxquant_raw_ratios()
#' }
#' @export
get_maxquant_raw_ratios <- function(DT){

   # Satisfy CHECK
   . <- NULL

   # Identify
   pattern <-'^Ratio ([HM]/[ML]) (.+)$'
   colnames <- names(DT) %>%
      magrittr::extract( stringi::stri_detect_regex(., pattern))                %>%
      magrittr::extract(!stringi::stri_detect_regex(., '(normalized|variability|type|iso|count)')) %>%
      magrittr::extract(!stringi::stri_detect_regex(., '(localized|unmod. pep.|nmods)'))           %>%  # phosphosites
      magrittr::extract(!stringi::stri_detect_fixed(., '___')) # phosphosites

   # Extract
   if (length(colnames)==0) return(NULL)
   values <- DT[, colnames, with = FALSE]
   names(values) %<>% stringi::stri_replace_all_regex(pattern, '$2[$1]')
   values
}

#' Extract lfq intensity columns
#' @param DT maxquant data.table
#' @return LFQ intensity data.table
#' @examples
#' require(magrittr)
#' if (require(alnoubi.2017)){
#'    system.file('extdata/proteinGroups.txt',
#'                 package = 'alnoubi.2017') %>%
#'    data.table::fread() %>%
#'    autonomics.import::get_maxquant_lfq_intensities()
#' }
#' if (require(graumann.lfq)){
#'    system.file('extdata/proteinGroups.txt',
#'                 package = 'graumann.lfq') %>%
#'    data.table::fread() %>%
#'    get_maxquant_lfq_intensities()
#' }
#' @export
get_maxquant_lfq_intensities <- function(DT){
   # Satisfy CHECK
   . <- experiment <- colname <- label <- newname <- NULL

   # PATTERN
   # Unlabeled: 'LFQ intensity WT.R1'
   # Labeled:   'LFQ intensity H WT.R1' (no 'LFQ intensity WT.R1',
   #                            even though 'Intensity WT.R1' (= Intensity H + L) exists
   pattern <- '^LFQ intensity( [HML])? (.+)$'
   design <- names(DT) %>% magrittr::extract( stringi::stri_detect_regex(., pattern)) %>% data.table::data.table(colname=.)
   design %>% magrittr::extract(, experiment := colname %>% stringi::stri_replace_first_regex(pattern, '$2'))
   design %>% magrittr::extract(, label      := colname %>% stringi::stri_replace_first_regex(pattern, '$1') %>% trimws())
   design[,          newname := experiment]
   design[label!='', newname := sprintf('%s[%s]', experiment, label)]

   if (length(design$colname)>0){
      values <- DT %>%
         magrittr::extract(, design$colname, with = FALSE)
      names(values) <- design$newname
      return(values)
   }

   # Return NULL elsewise
   return(NULL)
}

#' Extract reporter intensities
#' @param DT maxquant data.table
#' @return reporter intensity data.table
#' @examples
#' require(magrittr)
#' if (require(billing.vesicles)){
#'    system.file('extdata/proteinGroups.txt',
#'                 package = 'billing.vesicles') %>%
#'    data.table::fread() %>%
#'    get_maxquant_reporter_intensities()
#' }
#' @export
get_maxquant_reporter_intensities <- function(DT){
   # Satisfy CHECK
   . <- NULL

   # Identify
   pattern  <- '^Reporter intensity corrected ([0-9]+) (.+)$'
   colnames <- names(DT) %>%
      magrittr::extract( stringi::stri_detect_regex(., pattern))

   # Extract
   if (length(colnames)==0) return(NULL)
   values <- DT[, colnames, with = FALSE]
   names(values) %<>% stringi::stri_replace_all_regex(pattern, '$2[$1]')
   values
}

#' Extract raw intensities
#' @param DT maxquant data.table
#' @return raw intensity data.table
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/proteinGroups.txt',
#'                 package = 'billing.differentiation.data') %>%
#'    data.table::fread() %>%
#'    autonomics.import::get_maxquant_raw_intensities()
#' }
#' if (require(graumann.lfq)){
#'    system.file('extdata/proteinGroups.txt',
#'                 package = 'graumann.lfq') %>%
#'    data.table::fread() %>%
#'    autonomics.import::get_maxquant_raw_intensities()
#' }
#' if (require(alnoubi.2017)){
#'    system.file('extdata/proteinGroups.txt',
#'                 package = 'alnoubi.2017') %>%
#'    data.table::fread() %>%
#'    autonomics.import::get_maxquant_raw_intensities()
#' }
#' @importFrom data.table   data.table
#' @importFrom magrittr     %>%   %<>%
#' @export
get_maxquant_raw_intensities <- function(DT){

   # Satisfy CHECK
   . <- experiment <- n.labels <- label <- NULL

   # Define label-specific intensity columns
   pattern <- '^Intensity ([HML]) (.+$)'
   labeled_columns     <- names(DT) %>% magrittr::extract( stringi::stri_detect_regex(., pattern))
   labeled_experiments <- labeled_columns %>% stringi::stri_replace_first_regex(pattern, '$2')
   labeled_labels      <- labeled_columns %>% stringi::stri_replace_first_regex(pattern, '$1')
   labeled_newnames    <- sprintf('%s[%s]', labeled_experiments, labeled_labels)

   # Define label-agnostic intensity columns
   unlabeled_columns  <- names(DT) %>% magrittr::extract(   stringi::stri_detect_regex(., '^Intensity .+$')) %>%
      setdiff(magrittr::extract(., stringi::stri_detect_regex(., '^Intensity [HML]')))
   unlabeled_experiments <- unlabeled_columns %>% substr(11, nchar(.))
   unlabeled_labels   <- rep('', length(unlabeled_columns))
   unlabeled_newnames <- unlabeled_experiments

   # Restrict label agnostic intensity columns to labelfree experiments
   design <- data.table::data.table(column      = c(labeled_columns,     unlabeled_columns),
                                    experiment  = c(labeled_experiments, unlabeled_experiments),
                                    label       = c(labeled_labels,      unlabeled_labels),
                                    newname     = c(labeled_newnames,    unlabeled_newnames))
   design %>% magrittr::extract(, n.labels := sum(label != ''), by = experiment)
   design %<>% magrittr::extract(n.labels==0 | (n.labels>0 & label != ""))

   # Extract relevant intensity columns
   if (length(colnames)>0){
      values <- DT %>%
         magrittr::extract(, design$column, with = FALSE)
      names(values) <- design$newname
      return(values)
   }
}

#' Infer value type
#'
#' @param file maxquant file
#' @return value type
#' @examples
#' require(magrittr)
#'
#' # normalized.ratios
#' #------------------
#' if (require(billing.differentiation.data)){
#'    file <- system.file('extdata/maxquant/proteinGroups.txt',
#'                       package = 'billing.differentiation.data')
#'    file %>% infer_maxquant_value_type()
#' }
#'
#' if (require(alnoubi.2017)){
#'    file <- system.file('extdata/proteinGroups.txt',
#'                       package = 'alnoubi.2017')
#'    file %>% infer_maxquant_value_type()
#' }
#'
#' # reporter.intensities
#' #---------------------
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt',
#'                         package = 'billing.vesicles')
#'    file %>% infer_maxquant_value_type()
#' }
#'
#' # lfq intensities
#' #----------------
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt',
#'                         package = 'graumann.lfq')
#'    file %>% infer_maxquant_value_type()
#' }
#'
#' # raw intensities
#' #----------------
#' if (require(alnoubi.2017)){
#'    file <- system.file('extdata/proteinGroups.txt',
#'                         package = 'alnoubi.2017')
#'    file %>% infer_maxquant_value_type()
#' }
#' @export
infer_maxquant_value_type <- function(file){

   # Read
   DT <- data.table::fread(file)

   # Normalized ratios
   x <- DT %>% get_maxquant_normalized_ratios()
   if (!is.null(x)) return('normalized.ratio')

   # Reporter intensities
   x <- DT %>% get_maxquant_reporter_intensities()
   if (!is.null(x)) return('reporter.intensity')

   # LFQ intensities
   x <- DT %>% get_maxquant_lfq_intensities()
   if (!is.null(x)) return('lfq.intensity')

   # Raw ratios
   x <- DT %>% get_maxquant_raw_ratios()
   if (!is.null(x)) return('raw.ratio')

   # Raw labeled intensities
   x <- DT %>% get_maxquant_raw_intensities()
   if (!is.null(x)) return('raw.intensity')

   stop('infer_maxquant_value_type: unable to infer value columns')
}

#' Max Quant value types
#'@export
MAXQUANT_VALUE_TYPES <- c('normalized.ratio', 'reporter.intensity', 'lfq.intensity', 'raw.ratio', 'raw.intensity')


#' Extract value columns
#'
#' @param DT maxquant data.table
#' @param value_type any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @return value type (infer_maxquant_value_type), value column names (extract_value_colnames),
#'         sample names (update_value_colnames)
#' @examples
#' require(magrittr)
#'
#' # normalized.ratios
#' #------------------
#' if (require(billing.differentiation.data)){
#'    file <- system.file('extdata/maxquant/proteinGroups.txt',
#'                       package = 'billing.differentiation.data')
#'    value_type <- infer_maxquant_value_type(file)
#'    DT <- data.table::fread(file)
#'    get_maxquant_value_columns(DT, value_type)
#' }
#'
#' if (require(alnoubi.2017)){
#'    file <- system.file('extdata/proteinGroups.txt',
#'                       package = 'alnoubi.2017')
#'    value_type <- infer_maxquant_value_type(file)
#'    DT <- data.table::fread(file)
#'    get_maxquant_value_columns(DT, value_type)
#' }
#'
#' # reporter.intensities
#' #---------------------
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt',
#'                         package = 'billing.vesicles')
#'    value_type <- infer_maxquant_value_type(file)
#'    DT <- data.table::fread(file)
#'    get_maxquant_value_columns(DT, value_type)
#' }
#'
#' # lfq intensities
#' #----------------
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt',
#'                         package = 'graumann.lfq')
#'    value_type <- infer_maxquant_value_type(file)
#'    DT <- data.table::fread(file)
#'    get_maxquant_value_columns(DT, value_type)
#' }
#'
#' # raw intensities
#' #----------------
#' if (require(alnoubi.2017)){
#'    file <- system.file('extdata/proteinGroups.txt',
#'                         package = 'alnoubi.2017')
#'    value_type <- infer_maxquant_value_type(file)
#'    DT <- data.table::fread(file)
#'    get_maxquant_value_columns(DT, value_type)
#' }
#' @export
#' @importFrom magrittr %>%
#' @export
get_maxquant_value_columns <- function(DT, value_type){

   # Assert
   assertive.sets::assert_is_subset(value_type, MAXQUANT_VALUE_TYPES)

   # Infer value columns
   switch(value_type,
          normalized.ratio        = DT %>% get_maxquant_normalized_ratios(),
          reporter.intensity      = DT %>% get_maxquant_reporter_intensities(),
          lfq.intensity           = DT %>% get_maxquant_lfq_intensities(),
          raw.ratio               = DT %>% get_maxquant_raw_ratios(),
          raw.intensity           = DT %>% get_maxquant_raw_intensities()
   )
}

make_sample_names <- function(injection, channel = rep(1, length(injection)) ){
   if (length(unique(channel))==1){
      injection
   } else {
      paste0(injection, '[', channel, ']')
   }
}


#'@rdname create_maxquant_design_file
#'@export
create_sample_design_df <- function(...){
   .Deprecated('create_maxquant_design_df')
   create_maxquant_design_df(...)
}

#'@rdname create_maxquant_design_file
#'@export
create_sample_design_file <- function(...){
   .Deprecated('create_maxquant_design_file')
   create_maxquant_design_file(...)
}

#' @rdname create_maxquant_design_file
#' @importFrom magrittr %>%
#' @export
create_maxquant_design_df <- function(
   proteingroups_file,
   value_type = autonomics.import::infer_maxquant_value_type(proteingroups_file)
){

   # Assert
   #-------
   assertive.files::assert_all_are_dirs(dirname(proteingroups_file))
   #assertive.files::assert_all_are_readable_files(proteingroups_file, warn_about_windows = FALSE)

   # Read
   #-----
   DT <- data.table::fread(proteingroups_file)

   # Extract channel & injection
   #----------------------------
   value_columns <- DT %>% get_maxquant_value_columns(value_type)
   injection <- names(value_columns) %>% stringi::stri_replace_first_regex('\\[.+\\]', '')   # injection [channel]
   channel   <- names(value_columns) %>% stringi::stri_extract_first_regex('\\[.+\\]') %>%
      stringi::stri_replace_first_fixed('[', '') %>%
      stringi::stri_replace_first_fixed(']', '')
   # Return
   #-------
   sample_design <- data.frame(injection = injection, channel = channel, subgroup  = '', replicate = '')
   if (all(is.na(channel))){
      sample_design$channel <- NULL
   }
   sample_design

}


#' Create maxquant design dataframe/file
#'
#' @param proteingroups_file full path to protein groups file
#' @param sample_design_file full path to sample design file
#' @param value_type         any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @param ... backward compatibility to deprecated functions create_sample_design_df and create_sample_design_file
#' @return sample design dataframe (create_maxquant_design_df) or file (create_maxquant_design_file)
#' @examples
#' require(magrittr)
#'
#' # normalized ratios
#' #------------------
#' if (require(billing.differentiation.data)){
#'    proteingroups_file <- system.file('extdata/maxquant/proteinGroups.txt',
#'                                      package = 'billing.differentiation.data')
#'    autonomics.import::create_maxquant_design_df(proteingroups_file)
#'    sample_design_file <- paste0(tempdir(), '/billing.differentiation/sample_design.txt')
#'    dir.create(dirname(sample_design_file))
#'    create_maxquant_design_file(proteingroups_file, sample_design_file)
#' }
#'
#' # reporter intensities
#' #---------------------
#' if (require(billing.vesicles)){
#'    proteingroups_file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    autonomics.import::create_maxquant_design_df(proteingroups_file)
#'    sample_design_file <- paste0(tempdir(), '/billing.vesicles/sample_design.txt')
#'    dir.create(dirname(sample_design_file))
#'    create_maxquant_design_file(proteingroups_file, sample_design_file)
#' }
#'
#' # lfq intensities
#' #----------------
#' if (require(graumann.lfq)){
#'    proteingroups_file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    create_maxquant_design_df(proteingroups_file)
#'    sample_design_file <- paste0(tempdir(), '/graumann.lfq/sample_design.txt')
#'    dir.create(dirname(sample_design_file))
#'    create_maxquant_design_file(proteingroups_file, sample_design_file)
#' }
#'
#' # raw labeled intensities
#' #---------------------
#' if (require(billing.differentiation.data)){
#'    proteingroups_file <- system.file('extdata/maxquant/proteinGroups.txt',
#'                                      package = 'billing.differentiation.data')
#'    create_maxquant_design_df(proteingroups_file, value_type = 'raw.intensity')
#' }
#'
#'
#' # raw unlabeled intensities
#' #--------------------------
#' if (require(alnoubi.2017)){
#'    proteingroups_file <- system.file('extdata/proteinGroups.txt', package = 'alnoubi.2017')
#'    create_maxquant_design_df(proteingroups_file)
#' }
#' @importFrom magrittr %>%
#' @export
create_maxquant_design_file <- function(
   proteingroups_file,
   sample_design_file = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   value_type = infer_maxquant_value_type(proteingroups_file)
){

   # Abort
   if (file.exists(sample_design_file)) {
      autonomics.support::cmessage('\tAbort - file already exists: %s', sample_design_file)
      return(invisible(sample_design_file))
   }

   # Create
   sample_design_df <- create_maxquant_design_df(proteingroups_file = proteingroups_file, value_type = value_type)
   sample_design_df %>% autonomics.support::print2txt(sample_design_file)
   autonomics.support::cmessage('\tWriting to: %s', sample_design_file)
   autonomics.support::cmessage('\tOpen this file in a Excel or LibrOffice, complete it manually and save.')

   # Open
   if(interactive() && assertive.reflection::is_windows()){
      tryCatch(shell.exec(sample_design_file), return(invisible(sample_design_file)))
   }

   # Return
   return(invisible(sample_design_file))
}


################################################################################
#                                                                              #
#                proteingroups   ->   data.table                               #
#                phosphosites                                                  #
#                                                                              #
################################################################################

#' Get autonomics summary attributes
#' @param x object
#' @return autonomics summary attribute of object
#' @export
get_summary_attr <- function(x)
{
   tmp_summary_attr <- attr(x, 'autonomics_summary')
   if (is.null(tmp_summary_attr))
   {
      tmp_summary_attr <- list()
   }
   return(tmp_summary_attr)
}

#' Set autonomics summary attributes
#' @param x object
#' @param y autonomics summary attribute
#' @return updated object
#' @export
set_summary_attr <- function(x, y)
{
   attr(x, 'autonomics_summary') <- y
   return(x)
}

PROTEIN_VARS_MAXQUANT   <- c('id',         'Majority protein IDs', 'Protein names', 'Gene names', 'Fasta headers')
PHOSPHO_VARS_MAXQUANT   <- c('id',         'Protein group IDs', 'Proteins',           'Gene names', 'Positions within proteins', 'Phospho (STY) Probabilities', 'Sequence window', 'Amino acid', 'Fasta headers')
PROTEIN_VARS_AUTONOMICS <- c('feature_id', 'Uniprot accessions',   'Protein names', 'Gene names', 'Fasta headers')
PHOSPHO_VARS_AUTONOMICS <- c('feature_id', 'Protein group IDs', 'Uniprot accessions', 'Gene names', 'Positions within proteins', 'Phospho (STY) Probabilities', 'Sequence window', 'Amino acid', 'Fasta headers')

setnames_if_available <- function(DT, old, new){
   if(old %in% names(DT)){
      DT %>% data.table::setnames(old, new)
   }
   DT
}

#' Rename annotation columns
#' @param DT maxquant datatable
#' @importFrom magrittr %<>%
#' @export
rename_maxquant_annotation_columns <- function(DT){
   # Don't use stri_replace_fixed (which replaces id in 'amino acid' as well!)
   DT %<>% setnames_if_available('id',                   'feature_id')          # proteingroups
   DT %<>% setnames_if_available('Majority protein IDs', 'Uniprot accessions')  # proteingroups
   DT %<>% setnames_if_available('Proteins',             'Uniprot accessions')  # phosphosites
   DT
}

#' @importFrom data.table   data.table   :=
log2_transform_ratios <- function(DT){
   ratio <- NULL
   DT[, ratio:=log2(ratio)]
   data.table::setnames(DT, 'ratio', 'l2r')
   DT
}

#' Remove reverse hits
#' @param PG maxquant dataframe
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %<>%   %>%
#' @export
rm_reverse <- function(PG){
   # Satisfy CHECK
   Reverse <- NULL

   # Replace NA with empty string
   # (when all values are missing, fread reads this in as a numeric column with all NAs)
   PG[, Reverse := as.character(Reverse)]
   PG[is.na(Reverse), Reverse:='']

   # Remove reverse
   selector <- PG$Reverse != '+'
   PG %<>% magrittr::extract(selector, )
   sprintf('\t\tRetain %d/%d features after removing reverse proteins', sum(selector), length(selector)) %>% message()
   return(PG)
}

#' @importFrom stringr regex
CONTAMINANT_FIELD_RX <- '(?i)(?:potential )?contaminant(?-i)'

#' Remove contaminants
#' @param PG maxquant dataframe
#' @importFrom magrittr  %<>%   %>%
#' @export
rm_contaminants <- function(PG){

   # What is the name of the contaminant variable
   fields <- colnames(PG)
   contaminant_var <- fields[fields %>% stringr::str_detect(CONTAMINANT_FIELD_RX)]

   # This code didn't generalize well towards phosphosites files
   # And frankly, it's equally easy to manually change '+' into '' in the proteinGroups table
   #  # Allow specified contaminants
   #  permitted_contaminants = c(GFPL1_ZOASP = 'Q9U6Y5')
   #  selector <- PG$`Majority protein IDs` %>% stringr::str_detect(permitted_contaminants)
   #  PG[[contaminant_var]][selector] <- ''

   # Replace NA with empty string
   # (when all values are missing, fread reads this in as a numeric column with all NAs)
   PG[[contaminant_var]] %<>% as.character()
   PG[[contaminant_var]][is.na(PG[[contaminant_var]])] <- ''

   # Remove remaining contaminants
   selector <- PG[[contaminant_var]]!= '+'
   PG %<>% magrittr::extract(selector, )
   sprintf('\t\tRetain %d/%d features after removing contaminants', sum(selector), length(selector)) %>% message()

   # Return
   return(PG)
}

#' Load proteingroups file into wide data.table
#' @param proteingroups_file path to proteinGroups file
#' @param value_type any value in names(VALUE_COLNAME_PATTERNS)
#' @examples
#' require(magrittr)
#'
#' # NORMALIZED RATIOS
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/proteinGroups.txt', package = 'billing.differentiation.data') %>%
#'       load_proteingroups_to_wide_dt() %>%
#'       str()
#' }
#'
#' # REPORTER INTENSITIES
#' if (require(billing.vesicles)){
#'    system.file('extdata/proteinGroups.txt', package = 'billing.vesicles') %>%
#'       load_proteingroups_to_wide_dt() %>%
#'       str()
#' }
#'
#' # RAW RATIOS
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/proteinGroups.txt', package = 'billing.differentiation.data') %>%
#'       load_proteingroups_to_wide_dt(value_type = 'raw.ratio') %>%
#'       str()
#' }
#'
#' # LFQ INTENSITIES
#' if (require(graumann.lfq)){
#'    system.file('extdata/proteinGroups.txt', package = 'graumann.lfq') %>%
#'       load_proteingroups_to_wide_dt() %>%
#'       str()
#' }
#'
#' # RAW INTENSITIES
#' if (require(alnoubi.2017)){
#'    system.file('extdata/proteinGroups.txt', package = 'alnoubi.2017') %>%
#'       load_proteingroups_to_wide_dt(value_type = 'raw.intensity') %>%
#'       str()
#'    system.file('extdata/proteinGroups.txt', package = 'alnoubi.2017') %>%
#'       load_proteingroups_to_wide_dt() %>%
#'       str()
#' }
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%   %<>%
#' @export
load_proteingroups_to_wide_dt <- function(
   proteingroups_file,
   value_type = infer_maxquant_value_type(proteingroups_file)
){
   autonomics.support::cmessage('\tLoad proteingroups')
   DT <- suppressWarnings(data.table::fread(proteingroups_file, verbose = FALSE, integer64 = 'numeric', colClasses = c(id='character')))
   tmp_summary_attr <- list(proteingroups_n = nrow(DT))

   DT %<>% autonomics.import::rm_contaminants()
   tmp_summary_attr %<>% c(proteingroups_n_rmcontaminants = nrow(DT))

   DT %<>% autonomics.import::rm_reverse()
   tmp_summary_attr %<>% c(proteingroups_n_rmreverse = nrow(DT))

   autonomics.support::cmessage('\t\tExtract %s values', value_type)
   tmp_summary_attr %<>% c(proteingroups_valuetype = value_type)

   value_columns <- DT %>% autonomics.import::get_maxquant_value_columns(value_type = value_type)
   annotation_columns <- DT %>% magrittr::extract(, intersect(PROTEIN_VARS_MAXQUANT, names(DT)), with = FALSE)
   cbind(annotation_columns, value_columns) %>%
      set_summary_attr(tmp_summary_attr) %>%
      return()
}

#' Load phosphosites file into wide data.table
#' @param phosphosites_file path to phosphoSites file
#' @param value_type any value in names(VALUE_COLNAME_PATTERNS)
#' @param min_loc_prob minimum localization probability
#' @param remove_multiple_pg whether to phosphosites mapping to multiple protein groups (logical)
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/phospho (STY)Sites.txt',
#'                 package = 'billing.differentiation.data') %>%
#'       load_phosphosites_to_wide_dt() %>%
#'       str()
#' }
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%
#' @export
load_phosphosites_to_wide_dt <- function(
   phosphosites_file,
   value_type         = infer_maxquant_value_type(phosphosites_file),
   min_loc_prob       = 0.75,
   remove_multiple_pg = FALSE
){
   # Satisfy CHECK
   `Protein group IDs` <- `Localization prob` <- NULL

   autonomics.support::cmessage('\tLoad phosphosite ratios')
   DT <- suppressWarnings(data.table::fread(phosphosites_file, verbose = FALSE, integer64 = 'numeric', colClasses = c(id = 'character')))
   tmp_summary_attr <- list(phosphosites_n = nrow(DT))

   DT %<>% rm_contaminants()
   tmp_summary_attr %<>% c(phosphosites_n_rmcontaminants = nrow(DT))
   DT %<>% rm_reverse()
   tmp_summary_attr %<>% c(phosphosites_n_rmreverse      = nrow(DT))

   DT %<>% magrittr::extract(`Localization prob` >= min_loc_prob)
   tmp_summary_attr %<>% c(list(phosphosites_n_locprob = nrow(DT), phosphosites_min_locprob = min_loc_prob))
   if(remove_multiple_pg){
      DT %<>% magrittr::extract(!stringi::stri_detect_fixed(`Protein group IDs`, ';'))
      tmp_summary_attr %<>% c(phosphosites_n_singleprotid = nrow(DT))
   }
   autonomics.support::cmessage('\t\tExtract %s values', value_type)
   tmp_summary_attr %<>% c(phosphosites_value_type = value_type)
   annotation_columns <- DT %>% magrittr::extract(, PHOSPHO_VARS_MAXQUANT, with = FALSE)
   value_columns      <- DT %>% get_maxquant_value_columns(value_type = value_type)
   cbind(annotation_columns, value_columns) %>%
      set_summary_attr(tmp_summary_attr) %>%
      return()
}

#' Tidy max quant table
#' @param DT maxquant datatable
#' @param log2_transform logical
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/proteinGroups.txt', package = 'billing.differentiation.data') %>%
#'       load_proteingroups_to_wide_dt() %>% melt_wide_maxquant_dt() %>% print()
#' }
#' if (require(graumann.lfq)){
#'   system.file('extdata/proteinGroups.txt', package = 'graumann.lfq') %>%
#'      load_proteingroups_to_wide_dt() %>% melt_wide_maxquant_dt() %>% print()
#' }
#' if (require(billing.vesicles)){
#'   system.file('extdata/proteinGroups.txt', package = 'billing.vesicles') %>%
#'      load_proteingroups_to_wide_dt() %>% melt_wide_maxquant_dt() %>% print()
#' }
#' @importFrom magrittr     %>%   %<>%
#' @importFrom data.table   data.table   :=
#' @export
melt_wide_maxquant_dt <- function(DT, log2_transform = TRUE){

   # Satisfy CHECK
   .SD <- feature_id <- value <- NULL

   # Extract value vars
   annotation_vars <- unique(c(PROTEIN_VARS_MAXQUANT, PHOSPHO_VARS_MAXQUANT))
   value_vars <- setdiff(names(DT), annotation_vars)

   # Rename annotation vars
   DT %<>% autonomics.import::rename_maxquant_annotation_columns()

   # Convert integer columns to numeric to avoid warnings when melting
   # https://stackoverflow.com/questions/32940580/convert-some-column-classes-in-data-table
   integer_columns <- DT %>% vapply(is.integer, logical(1)) %>% which() %>% names()
   if (length(integer_columns) > 0){
      DT[, (integer_columns) := lapply(.SD, as.numeric), .SDcols = integer_columns]
   }

   # Melt
   n0 <- DT[, length(unique(feature_id))]
   tmp_summary_attr <- get_summary_attr(DT)
   DT %<>% data.table::melt(measure.vars = value_vars, variable.name = 'sample', value.name = 'value', na.rm = TRUE)
   DT %<>% set_summary_attr(tmp_summary_attr)

   # DT %<>% data.table::melt(id.vars = PROTEIN_VARS_MAXQUANT, variable.name = 'sample', value.name = 'value', na.rm = TRUE)
   n1 <- DT[, length(unique(feature_id))]
   autonomics.support::cmessage('\t\tRetain %d/%d features with a value for at least one sample', n1, n0)

   # Log2 transform
   if (log2_transform){
      autonomics.support::cmessage('\t\tLog2 transform.')
      # pseudoval <- if (quantity %in% c('lfq.intensity'))   0.005 else 0
      pseudoval <- 0
      DT[, value:=log2(value + pseudoval)]
   }
   DT
}

#' Load protein table
#' @param proteingroups_file     proteingroups file
#' @param value_type       any value in names(VALUE_COLNAME_PATTERNS)
#' @param log2_transform   whether to log2 transform
#' @return proteingroups data.table
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'   system.file('extdata/maxquant/proteinGroups.txt', package = 'billing.differentiation.data') %>%
#'      autonomics.import::load_proteingroups_to_long_dt() %>% str()
#' }
#' if (require(graumann.lfq)){
#'   system.file('extdata/proteinGroups.txt', package = 'graumann.lfq') %>%
#'      autonomics.import::load_proteingroups_to_long_dt() %>% str()
#' }
#' if (require(billing.vesicles)){
#'   system.file('extdata/proteinGroups.txt', package = 'billing.vesicles') %>%
#'      autonomics.import::load_proteingroups_to_long_dt() %>% str()
#' }
#' @importFrom magrittr %>%
#' @export
load_proteingroups_to_long_dt <- function(
   proteingroups_file,
   value_type     = infer_maxquant_value_type(proteingroups_file),
   log2_transform = TRUE
){
   output <- load_proteingroups_to_wide_dt(proteingroups_file = proteingroups_file, value_type = value_type)
   tmp_summary_attr <- output %>% get_summary_attr()
   long_dt <- output %>% autonomics.import::melt_wide_maxquant_dt(log2_transform = log2_transform)
   long_dt %>% set_summary_attr(tmp_summary_attr %>%
                                c(proteingroups_n_quantified = long_dt %>%
                                                               magrittr::extract2('feature_id') %>%
                                                               unique() %>%
                                                               length()))
}

#' Load datatable with phosphosite ratios
#' @param phosphosites_file  phosphosites file
#' @param value_type    any value in names(VALUE_COLNAME_PATTERNS)
#' @param min_loc_prob  mimum localization probability
#' @param remove_multiple_pg whether to remove phosphosites mapping to multiple protein groups (logical)
#' @return phosophosites data.table
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'   system.file('extdata/maxquant/phospho (STY)Sites.txt',
#'                package = 'billing.differentiation.data') %>%
#'   load_phosphosite_ratios_to_long_dt() %>%
#'   str()
#' }
#' @importFrom magrittr %>%
#' @export
load_phosphosite_ratios_to_long_dt <- function(
   phosphosites_file,
   value_type   = infer_maxquant_value_type(phosphosites_file),
   min_loc_prob = 0.75,
   remove_multiple_pg = FALSE
){
   output <- load_phosphosites_to_wide_dt(phosphosites_file = phosphosites_file,
                                          value_type   = value_type,
                                          min_loc_prob = min_loc_prob,
                                          remove_multiple_pg = remove_multiple_pg)
                                          tmp_summary_attr <- get_summary_attr(output)
   long_dt <- output %>% melt_wide_maxquant_dt()
   long_dt %>% set_summary_attr(tmp_summary_attr %>%
                                c(phosphosites_n_quantified = long_dt %>%
                                                              magrittr::extract2('feature_id') %>%
                                                              unique() %>%
                                                              length()))
}

# tidy_phospho_datatable <- function(DT){
#    DT %<>% melt(measure.vars = renamed_value_columns, variable.name = 'sample', value.name = 'value', na.rm = TRUE)
#    #DT %<>% melt(id.vars = PHOSPHO_VARS_MAXQUANT, variable.name = 'sample', value.name = 'ratio', na.rm = TRUE)
#
#    DT[, value:=log2(value)]
#    DT %<>% rename_maxquant_annotation_columns()
#    # DT %>% setnames('id', 'feature_id')
#    # DT %>% setnames('Proteins', 'Uniprot accessions')
#    DT[, feature_id := paste0('PS', formatC(feature_id, digits = max(floor(log10(feature_id))), flag = '0'))]
# }

#' Load datatable with phosphosite occupancies
#' @param phosphosites_file  phosphosites  file
#' @param proteingroups_file  proteingroups file
#' @param value_type    any value in names(VALUE_COLNAME_PATTERNS)
#' @param min_loc_prob  minimum localization probability (numeric)
#' @return datatable with phosphosite occupancies
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'   system.file('extdata/maxquant/phospho (STY)Sites.txt',
#'                package = 'billing.differentiation.data') %>%
#'   load_phosphosite_occupancies_to_long_dt() %>%
#'   str()
#' }
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%
#' @export
load_phosphosite_occupancies_to_long_dt <- function(
   phosphosites_file,
   proteingroups_file = paste0(dirname(phosphosites_file), '/proteinGroups.txt'),
   value_type   = infer_maxquant_value_type(proteingroups_file),
   min_loc_prob = 0.75
){
   # Satisfy CHECK
   `Protein group IDs` <- feature_id <- sample <- value <- protein.value <- phospho.value <- NULL

   # Load protein data
   proteins <- load_proteingroups_to_long_dt(proteingroups_file = proteingroups_file, value_type = value_type)
   tmp_summary_attr_prot <- get_summary_attr(proteins)
   proteins %<>% magrittr::extract(, list(feature_id, sample, value))
   data.table::setnames(proteins, 'feature_id', 'Protein group IDs')
   data.table::setnames(proteins, 'value', 'protein.value')
   proteins[, `Protein group IDs`:= as.character(`Protein group IDs`)]

   # Load phosphosite data
   phosphosites <- load_phosphosite_ratios_to_long_dt(
                      phosphosites_file = phosphosites_file,
                      value_type = value_type,
                      min_loc_prob = min_loc_prob,
                      remove_multiple_pg = TRUE)
   tmp_summary_attr_phos <- get_summary_attr(phosphosites)
   data.table::setnames(phosphosites, 'value', 'phospho.value')

   # Merge and compute occupancies
   autonomics.support::cmessage('\tCompute phosphosite occupancies')
   data.table::setkey(proteins,      `Protein group IDs`, sample)
   data.table::setkey(phosphosites,  `Protein group IDs`, sample)
   occupancies <- merge(phosphosites, proteins)
   occupancies[,  value := phospho.value - protein.value]
   occupancies %>% set_summary_attr(c(tmp_summary_attr_prot, tmp_summary_attr_phos))
}


# # Cast ratios into wide data table
# cast_ratios <- function(DT){
#    DT %>% dcast(feature_id ~ sample, value.var = 'l2r')
# }

get_autonomics_fvars <- function(entity){
   switch(
      entity,
      proteingroup = PROTEIN_VARS_AUTONOMICS,
      phosphosite  = PHOSPHO_VARS_AUTONOMICS
   )
}

backquote <- function(x)  paste0("`", x, "`")

#' @importFrom magrittr %>% %<>%
create_dcast_formula <- function(feature_vars){
   feature_vars %<>% backquote()
   feature_vars %>% paste0(collapse = ' + ') %>%
      paste0(' ~ sample')
}

################################################################################
#                                                                              #
#                   data.table -> SummarizedExperiment                         #
#                                                                              #
################################################################################


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
#'    load_maxquant_design(sample_file)
#' }
#' if (require(billing.differentiation.data)){
#'    sample_file <- system.file('extdata/maxquant/sample_design.txt',
#'                                       package = 'billing.differentiation.data')
#'    load_maxquant_design(sample_file)
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
load_maxquant_design <- function(sample_file){
   # Prevent CHECK warnings
   n <- NULL

   # Load
   #assertive.files::assert_all_are_readable_files(sample_file, warn_about_windows = FALSE)
   sample_design <- suppressWarnings(data.table::fread(sample_file,
                                                       verbose    = FALSE,
                                                       integer64  = 'numeric',
                                                       colClasses = c(injection = 'character',
                                                                      #channel   = 'character',
                                                                      subgroup  = 'character',
                                                                      replicate = 'character'),
                                                       data.table = FALSE))
   # Check contents
   assertive.strings::assert_all_are_non_empty_character(sample_design$injection)
   #assertive.strings::assert_all_are_non_empty_character(sample_design$channel)
   assertive.strings::assert_any_are_non_missing_nor_empty_character(sample_design$subgroup)

   # Create replicates if all absent.
   if (all_are_missing_or_empty_character(sample_design$replicate)){
      sample_design %<>% dplyr::group_by_('subgroup') %>%
         dplyr::mutate(replicate = as.character(1:n())) %<>%
         dplyr::ungroup() %>%
         as.data.frame()
   }

   # Return
   sample_design %<>% data.table::data.table()
   return(sample_design)
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
         x %<>% stringi::stri_replace_first_fixed(old_values[i], new_values[i])
         msg             %<>%    paste0(sprintf('\n\t\t%s -> %s', old_values[i], new_values[i]))
      }
      msg %<>% substr(3, nchar(.))
      if (verbose) message(msg)
   }
   return(x)

}


#' @importFrom magrittr  %>%   %<>%
add_maxquant_sdata <- function(object, sample_file){

   # Prevent CHECK notes
   sample_id <- injection <- channel <- NULL

   # Load
   sample_design <- load_maxquant_design(sample_file)
   sample_design %>% magrittr::extract(, sample_id := injection)
   if ('channel' %in% names(sample_design)){
      sample_design %>% magrittr::extract(channel!='', sample_id := sprintf('%s[%s]', injection, channel))
   }
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


#' @importFrom magrittr   %>%
load_maxquant_parameters <- function(parameter_file){
   data.table::fread(parameter_file) %>%
      (function(x) as.list(x$Value) %>% magrittr::set_names(x$Parameter))
}

add_maxquant_prepro <- function(object, entity, quantity, parameter_file){
   parameters <- load_maxquant_parameters(parameter_file)
   autonomics.import::prepro(object) <- autonomics.import::create_prepro_list(
      assay      = 'lcms',
      entity     = entity,
      quantity   = quantity,
      software   = 'maxquant',
      parameters = parameters)
   autonomics.import::annotation(object) <- basename(parameters$`Fasta file`)
   object
}

# Extract canonical accessions
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%  %<>%
extract_canonical_accessions <- function(object, verbose = FALSE){

   # Satisfy CHECK
   `Fasta headers` <- n.majority.accessions <- `Uniprot accessions` <- NULL
   db <- n.in.fasta.header <- pg.contains.swissprot <- pg.contains.canonical <- NULL
   swissprot.accessions <- canonical.accessions <- best.accessions <- NULL
   which.swissprot <- which.canonical <- which.best <- `Positions within proteins` <- NULL
   positions.within.best.proteins <- NULL

   # Extract fdata
   PG <- autonomics.import::fdata(object) %>% data.table::as.data.table()
   autonomics.support::cmessage('\t\tExtract canonical swissprot accessions')

   # Fasta headers column absent
   if (!'Fasta headers' %in% names(PG)){
      autonomics.support::cmessage('\t\t\tFasta headers column missing - abort')
      return(object)
   }

   # Fasta headers all NA
   if (PG[, all(is.na(`Fasta headers`))]){
      autonomics.support::cmessage('\t\t\tFasta header values missing - abort')
      return(object)
   }

   # Fasta headers corrupt
   if (all(! PG[, `Fasta headers`] %>% stringi::stri_detect_regex('^[>]?(sp|tr)') )){
      autonomics.support::cmessage('\t\t\tFasta header values corrupt - abort')
      return(object)
   }

   # Swissprot or trembl?
   # Limit to majority accessions.
   if (verbose) autonomics.support::cmessage('\t\t\tAre proteins in proteingroup from swissprot or tremble?')
   PG %>% magrittr::extract(,
                            db:= `Fasta headers` %>% data.table::tstrsplit(';')  %>%
                               stringi::stri_replace_all_fixed('>', '')      %>%
                               stringi::stri_extract_first_regex('^(tr|sp)') %>%
                               paste0(collapse = ';'),
                            by = 'feature_id')
   PG %>% magrittr::extract(,
                            n.majority.accessions := `Uniprot accessions` %>% data.table::tstrsplit(';') %>% length(),
                            by = 'feature_id')
   PG %>% magrittr::extract(,
                            n.in.fasta.header := db %>% data.table::tstrsplit(';') %>% length(),
                            by = 'feature_id')
   PG %>% magrittr::extract(,
                            db := db %>%
                               data.table::tstrsplit(';') %>%
                               magrittr::extract(seq(1, min(n.majority.accessions, n.in.fasta.header))) %>%
                               paste0(collapse = ';'),
                            by = 'feature_id')

   # Contains swissprot accessions?
   if (verbose) autonomics.support::cmessage('\t\t\tDoes protein group contain swissprot accessions?')
   PG %>% magrittr::extract(,
                            pg.contains.swissprot := db %>% stringi::stri_detect_fixed('sp'),
                            by = 'feature_id')

   # Which accessions are from swissprot?
   if (verbose) autonomics.support::cmessage('\t\t\tWhich uniprot accessions are from swissprot?')
   PG %>% magrittr::extract(pg.contains.swissprot==TRUE,
                            which.swissprot := data.table::tstrsplit(db, ';') %>%
                               magrittr::equals('sp') %>%
                               which() %>%
                               paste0(collapse = ';'),
                            by = 'feature_id')

   # Extract swissprot accessions
   if (verbose) autonomics.support::cmessage('\t\t\tExtract swissprot accessions')
   PG %>% magrittr::extract(pg.contains.swissprot==TRUE,
                            swissprot.accessions := data.table::tstrsplit(`Uniprot accessions`, ';') %>%
                               magrittr::extract(as.integer(unlist(data.table::tstrsplit(which.swissprot, ';')))) %>%
                               paste0(collapse = ';'),
                            by = 'feature_id')


   # Infer canonical accessions
   if (verbose) autonomics.support::cmessage('\t\t\tExtract canonical accessions')
   PG %>% magrittr::extract(pg.contains.swissprot==TRUE,
                            canonical.accessions := swissprot.accessions %>%
                               data.table::tstrsplit(';') %>%
                               unlist() %>%
                               stringi::stri_replace_all_regex('([-][0-9]+)', '') %>%
                               unique() %>%
                               paste0(collapse = ';'),
                            by = 'feature_id')
   PG[, pg.contains.canonical := FALSE]
   PG %>% magrittr::extract(pg.contains.swissprot==TRUE,
                            pg.contains.canonical := intersect(data.table::tstrsplit(`Uniprot accessions`, ';'),
                                                               data.table::tstrsplit(canonical.accessions, ';')) %>%
                               length() %>% magrittr::is_greater_than(0),
                            by = 'feature_id')
   PG %>% magrittr::extract(pg.contains.canonical==TRUE,
                            which.canonical    := which(data.table::tstrsplit(`Uniprot accessions`, ';') == canonical.accessions),
                            by = 'feature_id')

   # Infer best.accession
   # Note: use paste0(collapse=';') rather than magrittr::extract2(1)
   # The latter breaks when 'Uniprot accessions' is empty
   # This is the case for phosphosite '24023' in lendal.ovarian, which is a REV sequence without REV tag.
   PG %>% magrittr::extract(pg.contains.canonical==TRUE,
                            best.accessions := canonical.accessions %>% data.table::tstrsplit(';') %>% paste0(collapse = ';'),
                            by = 'feature_id')
   PG %>% magrittr::extract(pg.contains.canonical==FALSE & pg.contains.swissprot==TRUE,
                            best.accessions := swissprot.accessions %>% data.table::tstrsplit(';') %>% paste0(collapse = ';'),
                            by = 'feature_id')
   PG %>% magrittr::extract(pg.contains.canonical==FALSE & pg.contains.swissprot==FALSE,
                            best.accessions := `Uniprot accessions`  %>% data.table::tstrsplit(';') %>% paste0(collapse = ';'),
                            by = 'feature_id')
   PG %>% magrittr::extract(, which.best := which(data.table::tstrsplit(`Uniprot accessions`, ';') == best.accessions),  by = 'feature_id')

   # Positions within swissprot and canonical proteins
   if (verbose) autonomics.support::cmessage('\t\t\tExtract positions within canonical/swissprot accessions')
   if ('Positions within proteins' %in% names(PG)){ # for phospho esets
      PG %>% magrittr::extract(,
                               positions.within.best.proteins := data.table::tstrsplit(`Positions within proteins`, ';') %>%
                                  magrittr::extract2(which.best),
                               by = 'feature_id')
      #PG %>% magrittr::extract(pg.contains.swissprot==TRUE,
      #                         positions.within.swissprot.proteins := data.table::tstrsplit(`Positions within proteins`, ';') %>%
      #                                                                magrittr::extract(data.table::tstrsplit(`Uniprot accessions`, ';') %in%
      #                                                                                  data.table::tstrsplit(swissprot.accessions, ';')) %>%
      #                                                                paste0(collapse = ';'),
      #                         by = 'feature_id')
      #PG %>% magrittr::extract(pg.contains.canonical==TRUE,
      #                         positions.within.canonical.proteins := data.table::tstrsplit(`Positions within proteins`, ';') %>%
      #                                                                magrittr::extract(which.canonical) %>%
      #                                                                paste0(collapse = ';'),
      #                         by = 'feature_id')
   }


   # Cleanup
   PG[, n.majority.accessions := NULL]
   PG[, db                    := NULL]
   PG[, which.swissprot       := NULL]
   PG[, which.canonical       := NULL]
   PG[, which.best            := NULL]
   PG[, pg.contains.swissprot := NULL]
   PG[, pg.contains.canonical := NULL]

   # Update fdata and return
   autonomics.import::fdata(object) <- as.data.frame(PG) %>%
      magrittr::set_rownames(autonomics.import::fnames(object))
   return(object)
}

parse_gene_names_from_uniprot_fasta_hdrs <- function(x, sep = ';'){
   x %>% stringi::stri_extract_all_regex('(?<=GN=)[A-Za-z0-9_]+') %>%
      vapply(paste0, character(1), collapse = ';')
}

parse_protein_names_from_uniprot_fasta_hdrs <- function(x, sep = ';'){
   x %>% stringi::stri_extract_all_regex('[^|]+(?= OS=[A-Za-z ])') %>%
      vapply(paste0, character(1), collapse = ';') %>%
      stringi::stri_replace_all_regex('[A-Z0-9-]+_[A-Z]+ ', '')
}

#' Convert long data table into eset
#' @param DT                         long datatable
#' @param entity                     entity
#' @param quantity                   quantity
#' @param sample_file                sample design file
#' @param parameter_file             parameter file
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    maxquant_dir       <- system.file('extdata/maxquant', package = 'billing.differentiation.data')
#'    sample_file        <- paste0(maxquant_dir, '/sample_design.txt')
#'    proteingroups_file <- paste0(maxquant_dir, '/proteinGroups.txt')
#'    phosphosites_file  <- paste0(maxquant_dir, '/phospho (STY)Sites.txt')
#'    parameter_file     <- paste0(maxquant_dir, '/parameters.txt')
#'
#'    autonomics.import::load_proteingroups_to_long_dt(proteingroups_file) %>%
#'    autonomics.import::esetise_maxquant_dt(entity = 'proteingroup', quantity = 'normalized.ratio',
#'                           sample_file, parameter_file)
#'    autonomics.import::load_phosphosite_ratios_to_long_dt(phosphosites_file) %>%
#'    autonomics.import::esetise_maxquant_dt(entity = 'phosphosite', quantity = 'normalized.ratio',
#'                           sample_file, parameter_file)
#'    autonomics.import::load_phosphosite_occupancies_to_long_dt(phosphosites_file) %>%
#'    autonomics.import::esetise_maxquant_dt(entity = 'phosphosite', quantity = 'occupancy',
#'                           sample_file, parameter_file)
#' }
#' if (require(alnoubi.2017)){
#'    maxquant_dir <- system.file('extdata', package = 'alnoubi.2017')
#'    proteingroups_file <- paste0(maxquant_dir, '/proteinGroups.txt')
#'    sample_file    <- paste0(maxquant_dir, '/sample_design.txt')
#'    parameter_file <- paste0(maxquant_dir, '/parameters.txt')
#'    load_proteingroups_to_long_dt(proteingroups_file, value_type = 'lfq.intensity') %>%
#'       esetise_maxquant_dt(entity = 'proteingroup', quantity = 'lfq.intensity',
#'                           sample_file, parameter_file)
#' }
#' if (require(graumann.lfq)){
#'    maxquant_dir <- system.file('extdata', package = 'graumann.lfq')
#'    proteingroups_file <- paste0(maxquant_dir, '/proteinGroups.txt')
#'    sample_file <- paste0(maxquant_dir,    '/sample_design.txt')
#'    parameter_file <- paste0(maxquant_dir, '/parameters.txt')
#'    DT <- autonomics.import::load_proteingroups_to_long_dt(proteingroups_file)
#'    DT %>% autonomics.import::esetise_maxquant_dt(
#'              entity = 'proteingroup', quantity = 'lfq.intensity', sample_file, parameter_file)
#' }
#' @importFrom data.table  data.table
#' @importFrom magrittr    %<>%   %>%
#' @export
esetise_maxquant_dt <- function(
   DT,
   entity,
   quantity,
   sample_file,
   parameter_file
){
   # Prevent CHECK notes
   . <- NULL

   # Save custom attributes
   tmp_summary_attr <- autonomics.import::get_summary_attr(DT)

   # Dcast wide table
   annotation_cols <- get_autonomics_fvars(entity) %>% intersect(names(DT))
   dcast_formula   <- create_dcast_formula(annotation_cols)
   #value.var <- switch(quantity, ratio='value', occupancy='occupancy', intensity='intensity')
   DT %<>% data.table::dcast(dcast_formula, value.var = 'value')

   # Make eset
   feature_names <- DT$feature_id
   data_matrix   <- DT %>%
                    magrittr::extract(, !names(.) %in% annotation_cols, with = FALSE) %>%
                    data.matrix() %>%
                    magrittr::set_rownames(feature_names)
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = data_matrix))
   object %<>% add_maxquant_sdata(sample_file)
   autonomics.import::fdata(object) <- DT %>% magrittr::extract(, annotation_cols, with = FALSE) %>%
                                       as.data.frame(check.names = FALSE) %>%
                                       magrittr::set_rownames(feature_names)
   # Add prepro details
   object %<>% add_maxquant_prepro(entity, quantity, parameter_file)

   # Restore custom attributes
   autonomics.import::analysis(object) <- tmp_summary_attr

   # Assert validity
   autonomics.import::assert_is_valid_eset(object)

   # Return
   return(object)
}

################################################################################
#                                                                              #
#                  proteingroups  -> ->  SummarizedExperiment                  #
#                  phosphosites                                                #
#                                                                              #
################################################################################

#' Load proteingroup ratios/intensities
#'
#' @param proteingroups_file  full path to "proteinGroups.txt"
#' @param sample_file         full path to "sample_design.txt"
#' @param parameter_file      full path to "parameters.txt"
#' @param log2_transform      logical
#' @param value_type          any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @examples
#' library(magrittr)
#' if (require(autonomics.data)){
#'    system.file('extdata/billing2016/proteinGroups.txt', package = 'autonomics.data') %>%
#'    autonomics.import::load_proteingroups()
#' }
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/proteinGroups.txt',
#'                 package = 'billing.differentiation.data') %>%
#'    autonomics.import::load_proteingroups()
#' }
#' if (require(billing.vesicles)){
#'    system.file('extdata/proteinGroups.txt', package = 'billing.vesicles') %>%
#'    autonomics.import::load_proteingroups()
#' }
#' if (require(graumann.lfq)){
#'    system.file('extdata/proteinGroups.txt', package = 'graumann.lfq') %>%
#'    autonomics.import::load_proteingroups()
#' }
#' if (require(graumann.zebra)){
#'    system.file('extdata/proteinGroups.txt', package = 'graumann.zebra') %>%
#'    autonomics.import::load_proteingroups()
#' }
#' if (require(alnoubi.2017)){
#'    system.file('extdata/proteinGroups.txt', package = 'alnoubi.2017') %>%
#'    autonomics.import::load_proteingroups(value_type = 'raw.intensity')
#' }
#' @return data.table or ExpressionSet
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%
#' @export
load_proteingroups <- function(
   proteingroups_file,
   sample_file    = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   parameter_file = paste0(dirname(proteingroups_file), '/parameters.txt'),
   log2_transform = TRUE,
   value_type     = autonomics.import::infer_maxquant_value_type(proteingroups_file)
){
   # Satisfy CHECK
   feature_id <- NULL

   # Assert
   assertive.files::assert_all_are_existing_files(c(proteingroups_file, sample_file, parameter_file))

   # Load
   DT <- autonomics.import::load_proteingroups_to_long_dt(proteingroups_file, value_type = value_type, log2_transform = log2_transform)
   DT[, feature_id := paste0('PG', formatC(feature_id, digits = max(floor(log10(feature_id))), flag = '0'))]
   DT %>% autonomics.import::esetise_maxquant_dt(
             entity = 'proteingroup',
             quantity = value_type,
             sample_file = sample_file,
             parameter_file = parameter_file)
}

#' Load phosphosite ratios
#'
#' @param phosphosites_file   full path to "Phospho (STY)Sites.txt"
#' @param sample_file         full path to "sample_design.txt"
#' @param parameter_file      full path to "parameters.txt"
#' @param value_type          any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    'extdata/maxquant/Phospho (STY)Sites.txt' %>%
#'    system.file(package = 'billing.differentiation.data') %>%
#'    load_phosphosites()
#' }
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%
#' @export
load_phosphosites <- function(
   phosphosites_file,
   sample_file       = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   parameter_file    = paste0(dirname(proteingroups_file), '/parameters.txt'),
   value_type        = autonomics.import::infer_maxquant_value_type(phosphosites_file)
){
   # Satisfy CHECK
   feature_id <- NULL

   assertive.files::assert_all_are_existing_files(c(phosphosites_file, sample_file, parameter_file))
   DT <- autonomics.import::load_phosphosite_ratios_to_long_dt(phosphosites_file = phosphosites_file, value_type = value_type)
   DT[, feature_id := paste0('PS', formatC(feature_id, digits = max(floor(log10(feature_id))), flag = '0'))]
   DT %>% esetise_maxquant_dt(
             entity = 'phosphosite',
             quantity = value_type,
             sample_file,
             parameter_file)
}

#' Load phosphosite occupancies
#'
#' @param phosphosites_file   full path to "Phospho (STY)Sites.txt"
#' @param proteingroups_file  full path to "proteinGroups.txt"
#' @param sample_file         full path to "sample_design.txt"
#' @param parameter_file      full path to "parameters.txt"
#' @param value_type          any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    'extdata/maxquant/Phospho (STY)Sites.txt' %>%
#'     system.file(package = 'billing.differentiation.data') %>%
#'     load_phosphosite_occupanciess()
#' }
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%
#' @export
load_phosphosite_occupancies <- function(
   phosphosites_file,
   proteingroups_file = paste0(dirname(proteingroups_file), '/proteinGroups.txt'),
   sample_file        = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   parameter_file     = paste0(dirname(proteingroups_file), '/parameters.txt'),
   value_type         = infer_maxquant_value_type(proteingroups_file)
){
   assertive.files::assert_all_are_existing_files(c(proteingroups_file, phosphosites_file, sample_file, parameter_file))
   DT <- load_phosphosite_occupancies_to_long_dt(proteingroups_file = proteingroups_file, phosphosites_file  = phosphosites_file)
   DT %>% esetise_maxquant_dt(entity = 'phosphosite',
                              quantity = 'occupancy',
                              sample_file = sample_file,
                              parameter_file = parameter_file)
}

utils::globalVariables(c('.SD', 'goid', 'interpro', 'ngoid', 'ninterpro', 'score', '.', 'Protein names', 'Gene names'))

##############################################
# Annotate and deconvolute
##############################################
utils::globalVariables(c('keggid', 'reviewed', 'ngene', 'feature_id', 'existence', 'Canonical accessions'))

#' Annotate and deconvolute proteingroups SumExp
#'
#' Steps:
#' \enumerate{
#'    \item Separates proteingroups into uniprot accessions.
#'    \item Fetches up-to-date annotations for each from uniprot
#'    \item Keeps best annotated entries (per proteingroup)
#'    \item Collapses isoforms (per canonical accession)
#'    \item Keeps first of redundant uniprot entries
#'    \item Collapses paralogs (per protein group)
#'    \item Merges annotations into fdata(object) and returns
#' }
#'
#' @param object SummarizedExperiment
#' @examples
#' \dontrun{
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- system.file('extdata/billing2016/proteinGroups.txt',
#'                           package = 'autonomics.data') %>%
#'              autonomics.import::load_proteingroups() %>%
#'              magrittr::extract(1:10, )
#'    object %>% autonomics.import::fdata() %>% str()
#'    object %>% autonomics.import::annotate_proteingroups() %>%
#'               autonomics.import::fdata() %>% str()
#' }
#' }
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%
#' @export
annotate_proteingroups <- function(object){

   # rm existing annotation to avoid confusion
   autonomics.import::fdata(object)$`Gene names`    <- NULL
   autonomics.import::fdata(object)$`Protein names` <- NULL
   autonomics.import::fdata(object)$`Fasta headers` <- NULL

   # Extract
   uniprot_var <- autonomics.import::uniprot_var(object)
   fid_var     <- autonomics.import::fid_var(object)
   fdata1 <- autonomics.import::fdata(object)               %>%
             tidyr::separate_rows(uniprot_var, sep = ';') %>%
             data.table::data.table() %>%
             magrittr::extract(, ('Canonical accessions') := get(uniprot_var) %>% stringi::stri_replace_first_regex('[-][0-9]+$', ''))
   n0 <- nrow(object)
   n1 <- length(unique(fdata1$`Uniprot accessions`))
   n2 <- length(unique(fdata1$`Canonical accessions`))
   autonomics.support::cmessage('%d protein groups -> %d uniprot accessions -> %d canonical', n0, n1, n2)

   # Annotate
   autonomics.support::cmessage('Annotate with uniprot.ws')
   dt <- unique(fdata1$`Canonical accessions`) %>% autonomics.annotate::annotate_uniprot()
   fdata1 %<>% merge(dt, by.x = 'Canonical accessions', by.y = uniprot_var, all.x = TRUE)

   # Keep reviewed, drop unreviewed (per proteingroup)
   autonomics.support::cmessage('Simplify proteingroups')
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, .SD[reviewed == max(reviewed)], by = fid_var) %>% magrittr::extract(, reviewed:=NULL)
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: prefer reviewed (when available, per proteingroup)', n0, n1)

   # Keep best (annotation) score (per protein group)
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, .SD[score == max(score)],       by = fid_var) %>% magrittr::extract(, score:=NULL)
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: prefer best annotation score (per proteingroup)', n0, n1)

   # Keep best go annotation (per protein group)
   n0 <- nrow(fdata1)
   fdata1  %>% magrittr::extract( is.na(goid), ngoid := 0)
   fdata1  %>% magrittr::extract(!is.na(goid), ngoid := goid     %>% stringi::stri_count_fixed(';') %>% magrittr::add(1))
   fdata1 %<>% magrittr::extract(, .SD[ngoid  == max(ngoid)   ], by = fid_var)
   fdata1 %>%  magrittr::extract(, ngoid := NULL)
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: prefer best GO annotation (per protein group)', n0, n1)

   # Keep best interpro annotated (per protein group)
   n0 <- nrow(fdata1)
   fdata1  %>% magrittr::extract( is.na(interpro), ninterpro := 0)
   fdata1  %>% magrittr::extract(!is.na(interpro), ninterpro := interpro %>% stringi::stri_count_fixed(';') %>% magrittr::add(1))
   fdata1 %<>% magrittr::extract(, .SD[ninterpro == max(ninterpro)], by = fid_var)
   fdata1 %>%  magrittr::extract(, ninterpro := NULL)
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: prefer best interpro annotation (per protein group):', n0, n1)

   # Collapse isoforms (per canonical accession)
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, (uniprot_var):= paste0(get(uniprot_var), collapse = ';'), by = c(fid_var, 'Canonical accessions')) %>% unique()
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Collapse isoforms (per canonical accession)', n0, n1)

   # Keep those with "gene name" annotation (per protein group)
   n0 <- nrow(fdata1)
   fdata1[, ngene:=0]
   fdata1[`Gene names` != '', ngene:=1]
   fdata1 %<>% magrittr::extract(, .SD[ngene==max(ngene)], by = feature_id)
   fdata1[, ngene:=NULL]
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Prefer entries with gene names (per protein group)', n0, n1)

   # Keep best existence values (per protein group)
   n0 <- nrow(fdata1)
   fdata1[, existence:=as.numeric(existence)]
   fdata1 %<>% magrittr::extract(,.SD[existence==min(existence)], by = c(fid_var))
   fdata1[, existence := NULL]
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Keep entries with best uniprot existence score (per protein group)', n0, n1)

   # Keep only first of redundant entries (per gene name)
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, .SD[1], by = c(fid_var, 'Gene names'))
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Keep only first of redundant uniprot entries  (per protein group and gene name)', n0, n1)

   # Collapse paralogs (per protein group)
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, (uniprot_var)          :=  get(uniprot_var)      %>% paste0(collapse = '; '),                             by = fid_var)
   fdata1 %<>% magrittr::extract(, `Canonical accessions` := `Canonical accessions` %>% paste0(collapse = '; '),                             by = fid_var)
   fdata1 %<>% magrittr::extract(, `Gene names`           := `Gene names`           %>% paste0(collapse = '; '),                             by = fid_var)
   fdata1 %<>% magrittr::extract(,  keggid                :=  keggid                %>% paste0(collapse = '; '),                             by = fid_var)
   fdata1 %<>% magrittr::extract(,  goid                  :=  goid                  %>% autonomics.support::uniquify_collapsed_strings(';'), by = fid_var)
   fdata1 %<>% magrittr::extract(,  interpro              :=  interpro              %>% autonomics.support::uniquify_collapsed_strings(';'), by = fid_var)
   fdata1 %<>% magrittr::extract(, `Protein names`        := `Protein names`        %>% autonomics.support::commonify_strings(),             by = fid_var)
   fdata1 %<>% unique()
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Collapse paralogs (per protein group) and commonify protein names', n0, n1)

   # Merge back into object and return
   autonomics.support::cmessage('Merge %d annotations into SummarizedExperiment with %d protein groups', n1, nrow(object))
   fdata1 %<>% magrittr::extract(match(autonomics.import::fdata(object)$feature_id, .$feature_id), )
   fdata1 %<>% data.frame(row.names = .$feature_id, check.names = FALSE)
   autonomics.import::fdata(object) <- fdata1
   object
}
