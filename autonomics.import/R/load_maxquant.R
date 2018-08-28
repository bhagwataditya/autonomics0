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
      values <- DT %>% magrittr::extract(, design$colname, with = FALSE)
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
   colnames <- names(DT) %>% magrittr::extract( stringi::stri_detect_regex(., pattern))

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
   x <- DT %>% autonomics.import::get_maxquant_normalized_ratios()
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

#' Get npeptides matrix
#'
#' Get feature x sample matrix with number of peptides
#' @param DT data.table
#' @param value_type maxquant value type
#' @examples
#' require(magrittr)
#' file <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'          system.file(package = 'autonomics.data')
#' DT <- file %>% data.table::fread()
#' value_type <- autonomics.import::infer_maxquant_value_type(file)
#' autonomics.import::get_maxquant_npeptides(DT, value_type) %>% str()
#' @importFrom magrittr %>%
#' @export
get_maxquant_npeptides <- function(DT, value_type){

   # Assert
   assertive.types::assert_is_data.table(DT)
   assertive.sets::assert_is_subset(value_type, autonomics.import::MAXQUANT_VALUE_TYPES)
   assertive.base::assert_all_are_not_na(DT$id)

   # Extract
   sample_ids <- DT %>% autonomics.import::get_maxquant_value_columns(value_type) %>% names()
   peptide_columns <- sprintf('Razor + unique peptides %s', sample_ids %>% stringi::stri_replace_first_regex('\\[.+\\]', ''))
   DT %>% magrittr::extract(, peptide_columns, with = FALSE)  %>%
          magrittr::set_names(sample_ids)                     %>%
          data.matrix()                                       %>%
          magrittr::set_rownames(paste0('PG', formatC(DT$id, digits = max(floor(log10(DT$id))), flag = '0')))
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
set_summary_attr <- function(x, y){
   attr(x, 'autonomics_summary') <- y
   return(x)
}

#' PROTEIN_VARS_MAXQUANT
#' @export
PROTEIN_VARS_MAXQUANT   <- c('id',         'Majority protein IDs', 'Protein names', 'Gene names', 'Fasta headers')


#' PHOSPHO_VARS_MAXQUANT
#'@export
PHOSPHO_VARS_MAXQUANT   <- c('id',         'Protein group IDs', 'Proteins',           'Gene names', 'Positions within proteins', 'Phospho (STY) Probabilities', 'Sequence window', 'Amino acid', 'Fasta headers')


#' PROTEIN_VARS_AUTONOMICS
#' @export
PROTEIN_VARS_AUTONOMICS <- c('feature_id', 'Uniprot accessions',   'Protein names', 'Gene names', 'Fasta headers')


#' PHOSPHO_VARS_AUTONOMICS
#'@export
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

#' Remove reverse hits
#' @param PG maxquant dataframe
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %<>%   %>%
#' @export
rm_reverse_dt <- function(PG){
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
rm_contaminants_dt <- function(PG){

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

#' Rm missing id values
#' @param PG maxquant dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    PG <- system.file('extdata/stemdiff/maxquant/proteinGroups.txt',
#'                       package = 'autonomics.data') %>%
#'          data.table::fread()
#'    PG %>% dim()
#'    PG %>% autonomics.import::rm_missing_id_values() %>% dim()
#' }
#' @importFrom magrittr %<>% %>%
#' @export
rm_missing_id_values <- function(PG){
   # rm missing id values (occurs in stemdiff)
   # If not removed, dcast after melt is not possible.
   n0 <- nrow(PG)
   PG %<>% magrittr::extract(!is.na(.$id), )
   PG %<>% magrittr::extract(.$id != '', )
   n1 <- nrow(PG)
   if (n1<n0) autonomics.support::cmessage('\t\tRetain %d/%d features with available id values', n1, n0)
   return(PG)
}

#' Load proteingroups file into wide data.table
#' @param proteingroups_file path to proteinGroups file
#' @param value_type any value in names(VALUE_COLNAME_PATTERNS)
#' @examples
#' #require(magrittr)
#' #
#' # NORMALIZED RATIOS
#' #if (require(autonomics.data)){
#' #   proteingroups_file <- system.file(
#' #      'extdata/stemdiff/maxquant/proteinGroups.txt',
#' #       package = 'autonomics.data')
#' #   proteingroups_file %>% load_proteingroups_to_wide_dt() %>% dim()
#' #}
#'
#' # REPORTER INTENSITIES
#' #if (require(billing.vesicles)){
#' #   proteingroups_file <- system.file(
#' #      'extdata/proteinGroups.txt', package = 'billing.vesicles')
#' #   proteingroups_file %>% load_proteingroups_to_wide_dt() %>% dim()
#' #}
#'
#' # RAW RATIOS
#' #if (require(billing.differentiation.data)){
#' #   proteingroups_file <- system.file(
#' #      'extdata/maxquant/proteinGroups.txt', package = 'billing.differentiation.data')
#' #   proteingroups_file %>% load_proteingroups_to_wide_dt(value_type = 'raw.ratio') %>% dim()
#' #}
#'
#' # LFQ INTENSITIES
#' #if (require(graumann.lfq)){
#' #   proteingroups_file <- system.file(
#' #      'extdata/proteinGroups.txt', package = 'graumann.lfq')
#' #   proteingroups_file %>% load_proteingroups_to_wide_dt() %>% dim()
#' #}
#'
#' # RAW INTENSITIES
#' #if (require(alnoubi.2017)){
#' #   proteingroups_file <- system.file(
#' #      'extdata/proteinGroups.txt', package = 'alnoubi.2017')
#' #   proteingroups_file %>% load_proteingroups_to_wide_dt(value_type = 'raw.intensity') %>% dim()
#'
#' #    proteingroups_file <- system.file(
#' #      'extdata/proteinGroups.txt', package = 'alnoubi.2017')
#' #   proteingroups_file %>% load_proteingroups_to_wide_dt() %>% dim()
#' #}
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%   %<>%
#' @export
load_proteingroups_to_wide_dt <- function(
   proteingroups_file,
   value_type = infer_maxquant_value_type(proteingroups_file)
){
   autonomics.support::cmessage('\tLoad proteingroups')
   DT <- autonomics.support::cfread(proteingroups_file, colClasses = c(id = 'character'))
   tmp_summary_attr <- list(proteingroups_n = nrow(DT))

   # rm missing id values (occurs instemdiff)
   # If not removed, dcast after melt is not possible.
   DT %<>% autonomics.import::rm_missing_id_values()
   tmp_summary_attr %<>% c(proteingroups_n_validid = nrow(DT))

   # rm contaminants
   DT %<>% autonomics.import::rm_contaminants_dt()
   tmp_summary_attr %<>% c(proteingroups_n_rmcontaminants = nrow(DT))

   # rm reverse
   DT %<>% autonomics.import::rm_reverse_dt()
   tmp_summary_attr %<>% c(proteingroups_n_rmreverse = nrow(DT))

   # Extract relevant value_type
   autonomics.support::cmessage('\t\tExtract %s values', value_type)
   tmp_summary_attr %<>% c(proteingroups_valuetype = value_type)
   value_columns <- DT %>% autonomics.import::get_maxquant_value_columns(value_type = value_type)
   annotation_columns <- DT %>% magrittr::extract(, intersect(PROTEIN_VARS_MAXQUANT, names(DT)), with = FALSE)

   # Return
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
#' # require(magrittr)
#' # if (require(autonomics.data)){
#' #   phosphosites_file <- system.file(
#' #      'extdata/stemdiff/maxquant/phospho (STY)Sites.txt',
#' #       package = 'autonomics.data')
#' #   phosphosites_file %>% load_phosphosites_to_wide_dt() %>% dim()
#' #}
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
   DT <- autonomics.support::cfread(phosphosites_file, colClasses = c(id = 'character'))
   tmp_summary_attr <- list(phosphosites_n = nrow(DT))

   # rm missing id values (occurs instemdiff)
   # If not removed, dcast after melt is not possible.
   DT %<>% autonomics.import::rm_missing_id_values()
   tmp_summary_attr %<>% c(phosphosites_n_validid = nrow(DT))

   # rm contaminants
   DT %<>% rm_contaminants_dt()
   tmp_summary_attr %<>% c(phosphosites_n_rmcontaminants = nrow(DT))

   # rm reverse
   DT %<>% rm_reverse_dt()
   tmp_summary_attr %<>% c(phosphosites_n_rmreverse      = nrow(DT))

   # rm unreliable localization probabilities
   DT %<>% magrittr::extract(`Localization prob` >= min_loc_prob)
   tmp_summary_attr %<>% c(list(phosphosites_n_locprob = nrow(DT), phosphosites_min_locprob = min_loc_prob))

   # rm sites mapping to multiple proteingroups
   if(remove_multiple_pg){
      DT %<>% magrittr::extract(!stringi::stri_detect_fixed(`Protein group IDs`, ';'))
      tmp_summary_attr %<>% c(phosphosites_n_singleprotid = nrow(DT))
   }

   # Extract relevant valuetypes
   autonomics.support::cmessage('\t\tExtract %s values', value_type)
   tmp_summary_attr %<>% c(phosphosites_value_type = value_type)
   annotation_columns <- DT %>% magrittr::extract(, PHOSPHO_VARS_MAXQUANT, with = FALSE)
   value_columns      <- DT %>% get_maxquant_value_columns(value_type = value_type)

   # Return
   cbind(annotation_columns, value_columns) %>%
      set_summary_attr(tmp_summary_attr) %>%
      return()
}

#' Melt wide maxquant datatable
#' @param DT maxquant datatable
#' @param log2_transform logical
#' @param log2_offset offset added to value before log2 transformation ()
#' @examples
#' #require(magrittr)
#' #if (require(autonomics.data)){
#' #   DT <- system.file('extdata/stemdiff/maxquant/proteinGroups.txt',
#' #                      package = 'autonomics.data') %>%
#' #         autonomics.import::load_proteingroups_to_wide_dt()
#' #   DT %>% melt_wide_maxquant_dt() %>% print()
#' #}
#' #if (require(graumann.lfq)){
#' #  system.file('extdata/proteinGroups.txt', package = 'graumann.lfq') %>%
#' #     load_proteingroups_to_wide_dt() %>% melt_wide_maxquant_dt() %>% print()
#' #}
#' #if (require(billing.vesicles)){
#' #  system.file('extdata/proteinGroups.txt', package = 'billing.vesicles') %>%
#' #     load_proteingroups_to_wide_dt() %>% melt_wide_maxquant_dt() %>% print()
#' #}
#' @importFrom magrittr     %>%   %<>%
#' @importFrom data.table   data.table   :=
#' @export
melt_wide_maxquant_dt <- function(DT, log2_transform = TRUE, log2_offset = 0){

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
   tmp_summary_attr <- autonomics.import::get_summary_attr(DT)
   long_dt <- DT %>% data.table::melt(
                            measure.vars  = value_vars,
                            variable.name = 'sample',  # na.rm = FALSE is a must !!!
                            value.name    = 'value',   # Avoids an "only NaN values sample" from being dropped
                            na.rm         = FALSE)     # Which has happened in a labeled IP study (where H and L differed massively)
   long_dt %<>% autonomics.import::set_summary_attr(tmp_summary_attr)

   # Log2 transform
   if (log2_transform){
      autonomics.support::cmessage('\t\tLog2 transform (offset %d).', log2_offset)
      long_dt[, value:=log2(value + log2_offset)]
   }
   long_dt
}

#' Load protein table
#' @param proteingroups_file     proteingroups file
#' @param value_type       any value in names(VALUE_COLNAME_PATTERNS)
#' @param log2_transform   whether to log2 transform
#' @return proteingroups data.table
#' @examples
#' # require(magrittr)
#' # if (require(autonomics.data)){
#' #  proteingroups_file <- system.file(
#' #     'extdata/stemdiff/maxquant/proteinGroups.txt',
#' #     package = 'autonomics.data')
#' #  proteingroups_file %>% autonomics.import::load_proteingroups_to_long_dt() %>% str()
#' #}
#' #if (require(graumann.lfq)){
#' #  system.file('extdata/proteinGroups.txt', package = 'graumann.lfq') %>%
#' #     autonomics.import::load_proteingroups_to_long_dt() %>% str()
#' #}
#' #if (require(billing.vesicles)){
#' #  system.file('extdata/proteinGroups.txt', package = 'billing.vesicles') %>%
#' #     autonomics.import::load_proteingroups_to_long_dt() %>% str()
#' #}
#' @importFrom magrittr %>%
#' @export
load_proteingroups_to_long_dt <- function(
   proteingroups_file,
   value_type     = autonomics.import::infer_maxquant_value_type(proteingroups_file),
   log2_transform = TRUE
){
   output <- autonomics.import::load_proteingroups_to_wide_dt(proteingroups_file = proteingroups_file, value_type = value_type)
   tmp_summary_attr <- output %>% get_summary_attr()
   log2_offset <- if (value_type %in% c('raw.intensity', 'lfq.intensity')) 1 else 0
   long_dt <- output %>% autonomics.import::melt_wide_maxquant_dt(log2_transform = log2_transform, log2_offset = log2_offset)
   long_dt %>% set_summary_attr(tmp_summary_attr %>%
                                c(proteingroups_n_quantified = long_dt %>%
                                                               magrittr::extract2('feature_id') %>%
                                                               unique() %>%
                                                               length()))
}

#' Load datatable with phosphosite ratios
#' @param phosphosites_file  phosphosites file
#' @param value_type    any value in names(VALUE_COLNAME_PATTERNS)
#' @param log2_transform   whether to log2 transform
#' @param min_loc_prob  mimum localization probability
#' @param remove_multiple_pg whether to remove phosphosites mapping to multiple protein groups (logical)
#' @return phosophosites data.table
#' @examples
#' # require(magrittr)
#' # if (require(billing.differentiation.data)){
#' #  system.file('extdata/maxquant/phospho (STY)Sites.txt',
#' #               package = 'billing.differentiation.data') %>%
#' #  load_phosphosite_ratios_to_long_dt() %>%
#' #  str()
#' #}
#' @importFrom magrittr %>%
#' @export
load_phosphosite_ratios_to_long_dt <- function(
   phosphosites_file,
   value_type         = infer_maxquant_value_type(phosphosites_file),
   log2_transform     = TRUE,
   min_loc_prob       = 0.75,
   remove_multiple_pg = FALSE
){
   output <- load_phosphosites_to_wide_dt(phosphosites_file = phosphosites_file,
                                          value_type   = value_type,
                                          min_loc_prob = min_loc_prob,
                                          remove_multiple_pg = remove_multiple_pg)
                                          tmp_summary_attr <- get_summary_attr(output)
   long_dt <- output %>% melt_wide_maxquant_dt(log2_transform = log2_transform, log2_offset = 0)
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
#' @param log2_transform   whether to log2 transform
#' @param min_loc_prob  minimum localization probability (numeric)
#' @return datatable with phosphosite occupancies
#' @examples
#' #require(magrittr)
#' #if (require(billing.differentiation.data)){
#' #  system.file('extdata/maxquant/phospho (STY)Sites.txt',
#' #               package = 'billing.differentiation.data') %>%
#' #  load_phosphosite_occupancies_to_long_dt() %>%
#' #  str()
#' #}
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%
#' @export
load_phosphosite_occupancies_to_long_dt <- function(
   phosphosites_file,
   proteingroups_file = paste0(dirname(phosphosites_file), '/proteinGroups.txt'),
   value_type         = infer_maxquant_value_type(proteingroups_file),
   log2_transform     = TRUE,
   min_loc_prob       = 0.75
){
   # Satisfy CHECK
   `Protein group IDs` <- feature_id <- sample <- value <- protein.value <- phospho.value <- NULL

   # Load protein data
   proteins <- load_proteingroups_to_long_dt(
                  proteingroups_file = proteingroups_file,
                  value_type         = value_type,
                  log2_transform     = log2_transform)
   tmp_summary_attr_prot <- get_summary_attr(proteins)
   proteins %<>% magrittr::extract(, list(feature_id, sample, value))
   data.table::setnames(proteins, 'feature_id', 'Protein group IDs')
   data.table::setnames(proteins, 'value', 'protein.value')
   proteins[, `Protein group IDs`:= as.character(`Protein group IDs`)]

   # Load phosphosite data
   phosphosites <- load_phosphosite_ratios_to_long_dt(
                      phosphosites_file  = phosphosites_file,
                      value_type         = value_type,
                      log2_transform     = log2_transform,
                      min_loc_prob       = min_loc_prob,
                      remove_multiple_pg = TRUE)
   tmp_summary_attr_phos <- get_summary_attr(phosphosites)
   data.table::setnames(phosphosites, 'value', 'phospho.value')

   # Merge and compute occupancies
   autonomics.support::cmessage('\tCompute phosphosite occupancies')
   data.table::setkey(proteins,      `Protein group IDs`, sample)
   data.table::setkey(phosphosites,  `Protein group IDs`, sample)
   occupancies <- merge(phosphosites, proteins, all = TRUE) # all TRUE is required to avoid "all-NA samples" from being dropped
   occupancies %<>% magrittr::extract(!is.na(feature_id))   # but NA fetaures should be removed
   occupancies %>% magrittr::extract(,  value := if (log2_transform) phospho.value - protein.value else phospho.value / protein.value)
   occupancies %>% set_summary_attr(c(tmp_summary_attr_prot, tmp_summary_attr_phos))
}


# # Cast ratios into wide data table
# cast_ratios <- function(DT){
#    DT %>% dcast(feature_id ~ sample, value.var = 'l2r')
# }

# get_autonomics_fvars <- function(entity){
#    switch(
#       entity,
#       proteingroup = PROTEIN_VARS_AUTONOMICS,
#       phosphosite  = PHOSPHO_VARS_AUTONOMICS
#    )
# }


################################################################################
#                                                                              #
#                   data.table -> SummarizedExperiment                         #
#                                                                              #
################################################################################




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


#' Add maxquant sdata
#' @param object SummarizedExperiment
#' @param design_file sample file
#' @return sumexp with sdata added
#' @importFrom magrittr  %>%   %<>%
#' @export
add_maxquant_sdata <- function(object, design_file){

   # Prevent CHECK notes
   sample_id <- NULL

   # Load
   sample_design <- read_maxquant_design(design_file)
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

   # Nameify subgroups
   message('\tValidify subgroups')
   autonomics.import::sdata(object)$subgroup %<>% autonomics.support::nameify_strings(verbose = TRUE)

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

#' Convert long data table into eset
#' @param DT                         long datatable
#' @param entity                     entity
#' @param quantity                   quantity
#' @param design_file                sample design file
#' @param parameter_file             parameter file
#' @examples
#' #require(magrittr)
#' #if (require(autonomics.data)){
#' #   maxquant_dir       <- system.file('extdata/stemdiff/maxquant',
#' #                                      package = 'autonomics.data')
#' #   proteingroups_file <- paste0(maxquant_dir, '/proteinGroups.txt')
#' #   phosphosites_file  <- paste0(maxquant_dir, '/phospho (STY)Sites.txt')
#' #   parameter_file     <- paste0(maxquant_dir, '/parameters.txt')
#' #   design_file <- tempfile()
#' #   proteingroups_file %>% autonomics.import::write_maxquant_design(infer_from_sampleids = TRUE,
#' #                                                                   design_file = design_file)
#' #   DT <- proteingroups_file %>% autonomics.import::load_proteingroups_to_long_dt()
#' #   DT %>% autonomics.import::maxquant_dt_to_sumexp(
#' #             entity = 'proteingroup', quantity = 'normalized.ratio',
#' #             design_file, parameter_file)
#' #   DT <- phosphosites_file %>% autonomics.import::load_phosphosite_ratios_to_long_dt()
#' #   DT %>% autonomics.import::maxquant_dt_to_sumexp(
#' #             entity = 'phosphosite', quantity = 'normalized.ratio',
#' #             design_file, parameter_file)
#' #   DT <- autonomics.import::load_phosphosite_occupancies_to_long_dt(phosphosites_file)
#' #   DT %>% autonomics.import::maxquant_dt_to_sumexp(
#' #             entity = 'phosphosite', quantity = 'occupancy',
#' #             design_file, parameter_file)
#' #}
#' #if (require(graumann.lfq)){
#' #   maxquant_dir <- system.file('extdata', package = 'graumann.lfq')
#' #   proteingroups_file <- paste0(maxquant_dir, '/proteinGroups.txt')
#' #   design_file <- paste0(maxquant_dir,    '/sample_design.txt')
#' #   parameter_file <- paste0(maxquant_dir, '/parameters.txt')
#' #   DT <- autonomics.import::load_proteingroups_to_long_dt(proteingroups_file)
#' #   DT %>% autonomics.import::maxquant_dt_to_sumexp(
#' #             entity = 'proteingroup', quantity = 'lfq.intensity', design_file, parameter_file)
#' #}
#' @importFrom data.table  data.table
#' @importFrom magrittr    %<>%   %>%
#' @export
maxquant_dt_to_sumexp <- function(
   DT,
   entity,
   quantity,
   design_file,
   parameter_file
){
   # Prevent CHECK notes
   . <- NULL

   # Save custom attributes

   tmp_summary_attr <- autonomics.import::get_summary_attr(DT)

   # Dcast wide table
   annotation_cols <- switch(entity, proteingroup = PROTEIN_VARS_AUTONOMICS, phosphosite  = PHOSPHO_VARS_AUTONOMICS) %>%
                      intersect(names(DT))
   dcast_formula   <- annotation_cols %>% (function(x)  paste0("`", x, "`")) %>% paste0(collapse = ' + ') %>% paste0(' ~ sample')

   #value.var <- switch(quantity, ratio='value', occupancy='occupancy', intensity='intensity')
   wide_dt <- DT %>% data.table::dcast(dcast_formula, value.var = 'value')
   wide_dt %<>% magrittr::extract(!is.na(feature_id))
           # This is required for phospho occupancy datasets
           # merge(proteingroups, phosphosites, all = TRUE) is used to avoid an all-NA-sample from being dropped.
           # But now full NA fetaures need to be dropped

   # Make SummarizedExperiment
   feature_names <- wide_dt$feature_id
   data_matrix   <- wide_dt %>%
                    magrittr::extract(, !names(.) %in% annotation_cols, with = FALSE) %>%
                    data.matrix() %>%
                    magrittr::set_rownames(feature_names)
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = data_matrix))
   object %<>% add_maxquant_sdata(design_file)
   autonomics.import::fdata(object) <- wide_dt %>% magrittr::extract(, annotation_cols, with = FALSE) %>%
                                       as.data.frame(check.names = FALSE) %>%
                                       magrittr::set_rownames(feature_names)
   # Add prepro details
   object %<>% add_maxquant_prepro(entity, quantity, parameter_file)

   # Restore custom attributes
   autonomics.import::analysis(object) <- tmp_summary_attr

   # Assert validity
   autonomics.import::assert_is_valid_eset(object)

   # Filter rows with no available value for any sample
   # This filtering is placed here because it
   #    (1) is generic to proteingroups and phosphosites
   #    (2) does not belong in the wideDT -> longDT step (where a full NaN sample gets completely removed from the data - see comments there)
   #    (3) should come after including sample data (and hence removing unannotated samples)
   autonomics.support::cmessage('\tFilter features')
   object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample()

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
#' @param design_file         full path to "sample_design.txt"
#' @param parameter_file      full path to "parameters.txt"
#' @param log2_transform      logical
#' @param value_type          any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @examples
#' #library(magrittr)
#' #if (require(autonomics.data)){
#' #
#' #   # STEM CELL COMPARISON
#' #   proteingroups_file <- system.file('extdata/stemcomp/maxquant/proteinGroups.txt',
#' #                                      package = 'autonomics.data')
#' #   design_file <- tempfile()
#' #   autonomics.import::write_maxquant_design(proteingroups_file,
#' #      infer_from_sampleids = TRUE, design_file = design_file)
#' #   proteingroups_file %>% autonomics.import::load_proteingroups(design_file = design_file)
#' #
#' #   # STEM CELL DIFFERENTIATION
#' #   proteingroups_file <- system.file(
#' #      'extdata/stemdiff/maxquant/proteinGroups.txt',
#' #      package = 'autonomics.data')
#' #   design_file <- tempfile()
#' #   autonomics.import::write_maxquant_design(proteingroups_file,
#' #      design_file = design_file, infer_from_sampleids = TRUE)
#' #   proteingroups_file %>% autonomics.import::load_proteingroups(design_file = design_file)
#' #}
#' #if (require(billing.vesicles)){
#' #   proteingroups_file <- system.file('extdata/proteinGroups.txt',
#' #                                      package = 'billing.vesicles')
#' #   proteingroups_file %>% autonomics.import::load_proteingroups()
#' #}
#' #if (require(graumann.lfq)){
#' #   proteingroups_file <- system.file('extdata/proteinGroups.txt',
#' #                                      package = 'graumann.lfq') %>%
#' #   proteingroups_file %>% autonomics.import::load_proteingroups()
#' #}
#' #if (require(graumann.zebra)){
#' #   proteingroups_file <- system.file('extdata/proteinGroups.txt',
#' #                                      package = 'graumann.zebra')
#' #   proteingroups_file %>% autonomics.import::load_proteingroups()
#' #}
#' @return data.table or ExpressionSet
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%
#' @export
old_load_proteingroups <- function(
   proteingroups_file,
   design_file    = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   parameter_file = paste0(dirname(proteingroups_file), '/parameters.txt'),
   log2_transform = TRUE,
   value_type     = autonomics.import::infer_maxquant_value_type(proteingroups_file)
){
   # Satisfy CHECK
   feature_id <- NULL

   # Assert
   assertive.files::assert_all_are_existing_files(c(proteingroups_file, design_file, parameter_file))

   # Load
   DT <- autonomics.import::load_proteingroups_to_long_dt(proteingroups_file,
                                                          log2_transform = log2_transform,
                                                          value_type     = value_type)
   DT[, feature_id := paste0('PG', formatC(feature_id, digits = max(floor(log10(feature_id))), flag = '0'))]
   DT %>% autonomics.import::maxquant_dt_to_sumexp(
             entity = 'proteingroup',
             quantity = value_type,
             design_file = design_file,
             parameter_file = parameter_file)
}

#' Load phosphosite ratios
#'
#' @param phosphosites_file   full path to "Phospho (STY)Sites.txt"
#' @param design_file         full path to "sample_design.txt"
#' @param parameter_file      full path to "parameters.txt"
#' @param log2_transform      logical
#' @param value_type          any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    'extdata/maxquant/Phospho (STY)Sites.txt'             %>%
#'    system.file(package = 'billing.differentiation.data') %>%
#'    autonomics.import::load_phosphosites()
#' }
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%
#' @export
load_phosphosites <- function(
   phosphosites_file,
   design_file       = paste0(dirname(phosphosites_file), '/sample_design.txt'),
   parameter_file    = paste0(dirname(phosphosites_file), '/parameters.txt'),
   log2_transform    = TRUE,
   value_type        = autonomics.import::infer_maxquant_value_type(phosphosites_file)
){
   # Satisfy CHECK
   feature_id <- NULL

   assertive.files::assert_all_are_existing_files(c(phosphosites_file, design_file, parameter_file))
   DT <- autonomics.import::load_phosphosite_ratios_to_long_dt(phosphosites_file = phosphosites_file,
                                                               log2_transform    = log2_transform,
                                                               value_type        = value_type)
   DT[, feature_id := paste0('PS', formatC(feature_id, digits = max(floor(log10(feature_id))), flag = '0'))]
   DT %>% maxquant_dt_to_sumexp(
             entity = 'phosphosite',
             quantity = value_type,
             design_file,
             parameter_file)
}

#' Load phosphosite occupancies
#'
#' @param phosphosites_file   full path to "Phospho (STY)Sites.txt"
#' @param proteingroups_file  full path to "proteinGroups.txt"
#' @param design_file         full path to "sample_design.txt"
#' @param parameter_file      full path to "parameters.txt"
#' @param log2_transform      logical
#' @param value_type          any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    'extdata/maxquant/Phospho (STY)Sites.txt'              %>%
#'     system.file(package = 'billing.differentiation.data') %>%
#'     autonomics.import::load_phosphosite_occupancies()
#' }
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%
#' @export
load_phosphosite_occupancies <- function(
   phosphosites_file,
   proteingroups_file = paste0(dirname(phosphosites_file), '/proteinGroups.txt'),
   design_file        = paste0(dirname(phosphosites_file), '/sample_design.txt'),
   parameter_file     = paste0(dirname(phosphosites_file), '/parameters.txt'),
   log2_transform     = TRUE,
   value_type         = autonomics.import::infer_maxquant_value_type(phosphosites_file)
){
   assertive.files::assert_all_are_existing_files(c(proteingroups_file, phosphosites_file, design_file, parameter_file))
   DT <- autonomics.import::load_phosphosite_occupancies_to_long_dt(
            phosphosites_file  = phosphosites_file,
            proteingroups_file = proteingroups_file,
            log2_transform     = log2_transform)
   DT[, feature_id := paste0('PS', formatC(feature_id, digits = max(floor(log10(feature_id))), flag = '0'))]
   DT %>% autonomics.import::maxquant_dt_to_sumexp(
             entity         = 'phosphosite',
             quantity       = 'occupancy',
             design_file    = design_file,
             parameter_file = parameter_file)
}


