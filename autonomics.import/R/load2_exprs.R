#============================
# MAXQUANT
#============================

#' Extract maxquant intensity colnames
#' @param file full path to proteinGroups.txt
#' @param quantity 'Intensity', 'LFQ intensity', 'Reporter intensity'
#' @return character vector with intensity column names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames() %>% head(3)
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames() %>% head(3)
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames('LFQ intensity') %>% head(3)
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames('Reporter intensity') %>% head(3)
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_intensity_colnames <- function(file, quantity = 'Intensity'){

   # Asser
   assertive.files::assert_all_are_existing_files(file)
   assertive.sets::assert_is_subset(quantity, c('Intensity', 'LFQ intensity', 'Reporter intensity'))

   # Deduce injections and channels from unambiguous peptide columns
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()

   # Construct intensity colnames
   intensity_colnames <- if (length(channels)==0){                      sprintf('%s %s',    quantity,           injections)
   } else {                  autonomics.support::vsprintf('%s %s %s', quantity, channels, injections) }

   # Ensure identical order as in actual file
   names(autonomics.support::cfread(file)) %>% magrittr::extract(. %in% intensity_colnames)
}


#' Extract maxquant ratio colnames
#' @param file full path to proteinGroups.txt
#' @param quantity 'Ratio' or 'Ratio normalized'
#' @return character vector with ratio column names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% extract_maxquant_ratio_colnames('Ratio') %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_ratio_colnames <- function(file, quantity){

   assertive.sets::assert_is_subset(quantity, c('Ratio', 'Ratio normalized'))

   # Deduce injections and channels from unambiguous peptide columns
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()

   # Construct possible ratio colnames
   possible_ratio_columns <- autonomics.support::vsprintf('Ratio %s/%s%s %s',
                                                          channels,
                                                          channels,
                                                          if(quantity == 'Ratio normalized') ' normalized' else '',
                                                          injections)
   # Return actual colnames in correct order
   file %>%
      autonomics.support::cfread() %>%
      names() %>%
      magrittr::extract(. %in% possible_ratio_columns)
}

#' Infer maxquant quantity
#' @param file path to maxquant file
#' @return string
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::infer_maxquant_quantity()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::infer_maxquant_quantity()
#' }
#' @importFrom magrittr %>%
#' @export
infer_maxquant_quantity <- function(file){
   x <- autonomics.support::cfread(file) %>% names()
   if (any(stringi::stri_detect_fixed(x, 'Ratio') & stringi::stri_detect_fixed(x, 'normalized')))   return('Ratio normalized')
   if (any(stringi::stri_detect_fixed(x, 'Reporter intensity')))                                    return('Reporter intensity')
   if (any(stringi::stri_detect_fixed(x, 'LFQ intensity')))                                         return('LFQ intensity')
   return('Intensity')
}


#' Load maxquant exprs
#' @param file full path to proteinGroups.txt
#' @param quantity 'Ratio', 'Ratio normalized', Intensity', 'LFQ intensity', 'Reporter intensity'
#' @return character vector with sample names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% load_exprs_maxquant('Ratio')            %>% extract(1:3, 1:3)
#'    file %>% load_exprs_maxquant('Ratio normalized') %>% extract(1:3, 1:3)
#'    file %>% load_exprs_maxquant('Intensity')        %>% extract(1:3, 1:3)
#'             extract(1:3, 1:3)
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% load_exprs_maxquant('LFQ intensity') %>% extract(1:3, 1:3)
#'    file %>% load_exprs_maxquant('Intensity')     %>% extract(1:3, 1:3)
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% load_exprs_maxquant('Reporter intensity') %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_exprs_maxquant <- function(file, quantity){

   # Assert
   assertive.sets::assert_is_subset(quantity, c('Ratio', 'Ratio normalized', 'Intensity', 'LFQ intensity', 'Reporter intensity'))

   # Construct maxquant colnames
   col_names <- if (quantity %in% c('Ratio', 'Ratio normalized')){ file %>% autonomics.import::extract_maxquant_ratio_colnames(quantity)
   } else {                                                        file %>% autonomics.import::extract_maxquant_intensity_colnames(quantity)
   }

   # Extract exprs matrix
   exprs_mat <- file %>%
                autonomics.support::cfread(select = col_names) %>%
                data.matrix() %>%
                magrittr::set_rownames(autonomics.import::load_fnames_maxquant(file))

   # Rename samples
   colnames(exprs_mat) <- autonomics.import::load_snames_maxquant(file, quantity = quantity, infer_design_from_sampleids = FALSE)
   # if (quantity %in% c('Ratio', 'Ratio normalized')){
   #    colnames(exprs_mat) %<>% stringi::stri_replace_first_regex('Ratio (./.) (?:normalized )?(.+)',   '$2[$1]')
   #
   # } else {
   #    # It is better to do it this way (rather than use a stri_replace_regex as for the ratios)
   #    # because ' ' separators in sample names are difficult to differentiate from ' H' constructs
   #    injections <- file %>% autonomics.import::extract_maxquant_injections()
   #    channels   <- file %>% autonomics.import::extract_maxquant_channels()
   #    colnames(exprs_mat) <- if (length(channels)==0) injections else autonomics.support::vsprintf('%s[%s]', injections, channels)
   # }

   # Return
   exprs_mat
}


#======================================================================
# METABOLON
#======================================================================

#' Load metabolon exprs
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_exprs_metabolon(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_exprs_metabolon <- function(file, sheet){

   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   sdata1 <- file %>% autonomics.import::load_sdata_metabolon(sheet = sheet)
   fdata1 <- file %>% autonomics.import::load_fdata_metabolon(sheet = sheet)

   df %>% magrittr::extract((fstart+1):nrow(.), (sstart+1):ncol(.)) %>%
      data.matrix() %>%
      magrittr::set_colnames(sdata1$CLIENT_IDENTIFIER) %>%
      magrittr::set_rownames(fdata1$MCOMP_ID)
}


#=========================================================================
# METABOLONLIPIDS
#=========================================================================

#' Load metabolonlipids exprs
#'
#' Load exprs from metabolon clp (complex lipid panel) file
#' @param file path to metabolon clp (complex lipid panel) file
#' @param sheet name of excel sheet (any value in METABOLONLIPIDS_SHEETS)
#' @return exprs matrix
#' @examples
#' require(magrittr)
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_exprs_metabolonlipids('Lipid Class Concentrations') %>% extract(1:3,1:3)
#'    file %>% load_exprs_metabolonlipids(    'Species Concentrations') %>% extract(1:3,1:3)
#'    file %>% load_exprs_metabolonlipids( 'Fatty Acid Concentrations') %>% extract(1:3,1:3)
#'    file %>% load_exprs_metabolonlipids('Lipid Class Compositions')   %>% extract(1:3,1:3)
#'    file %>% load_exprs_metabolonlipids(    'Species Compositions')   %>% extract(1:3,1:3)
#'    file %>% load_exprs_metabolonlipids( 'Fatty Acid Compositions')   %>% extract(1:3,1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_exprs_metabolonlipids <- function(file, sheet){
   all_sheets <- readxl::excel_sheets(file) %>% (function(x) x %>% magrittr::set_names(x))

   x <- file %>% readxl::read_excel(sheet = sheet)
   row1 <- which(x[[2]]=='Client Identifier')
   coln <- which(x[row1, ] == 'Unit')
   fnamesrow <- if (stringi::stri_detect_fixed(all_sheets[[sheet]], 'Lipid Class')) row1 else row1-1
   snames1 <- x[[2]][(1+row1):nrow(x)]

   x %>% magrittr::extract((1+row1):nrow(x), (coln+1):ncol(x))  %>%
      magrittr::set_names(x[fnamesrow, (coln+1):ncol(x)])    %>%
      (function(y) {y[y=='.'] <- NA;  y})                     %>%
      data.matrix()                                          %>%
      magrittr::set_rownames(snames1)                        %>%
      t()
}

#=======================================================================================
# SOMA
#========================================================================================

#' Load soma exprs
#' @param file string: path to adat file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::load_exprs_soma() %>% extract(1:3, 1:3)
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::load_exprs_soma() %>% extract(1:3, 1:3)
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
load_exprs_soma <- function(file){
   x      <- file %>% autonomics.import::identify_soma_structure()
   fdata1 <- file %>% autonomics.import::load_fdata_soma()
   sdata1 <- file %>% autonomics.import::load_sdata_soma()

   file %>%
      data.table::fread(header = FALSE, sep = '\t', fill = TRUE) %>%
      magrittr::extract(x$row:nrow(.), (x$col):ncol(.))          %>%
      t()                                                        %>%
      data.matrix()                                              %>%
      magrittr::set_rownames(fdata1$SeqId)                       %>%
      magrittr::set_colnames(sdata1$SampleId)                    %>%
      (function(x){class(x) <- 'numeric'; x})
}

#=======================================================================
# GENERIC
#=======================================================================

#' Load exprs
#' @param file path to omics data file
#' @param platform 'metabolon', 'metabolonlipids'
#' @param sheet excel sheet number or name (if applicable)
#' @param quantity string: quantity to be loaded as exprs (if applicable)
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'     file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'     file %>% load_exprs('metabolon', sheet = 2) %>% extract(1:3, 1:3)
#'  }
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_exprs('metabolonlipids', sheet = 'Lipid Class Concentrations') %>%
#'             extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_exprs <- function(
   file,
   platform,
   sheet = NULL,
   quantity = NULL
){
   switch(platform,
          maxquant        = file %>% load_exprs_maxquant(quantity = quantity),
          metabolonlipids = file %>% load_exprs_metabolonlipids(sheet = sheet),
          metabolon       = file %>% load_exprs_metabolon(      sheet = sheet),
          soma            = file %>% load_exprs_soma())
}

