#============================================
# MAXQUANT
#============================================

#' Load maxquant fnames
#' @param file full path to proteinGroups.txt
#' @return character vector with sample names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::load_fnames_maxquant() %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_fnames_maxquant <- function(file){
   file %>% autonomics.support::cfread() %>%
            magrittr::extract2('Majority protein IDs') %>%
            stringi::stri_split_fixed(';') %>%
            vapply(extract, character(1), 1)
}


#' Load maxquant fdata
#' @param file path to proteinGroups.txt file
#' @return dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::load_fdata_maxquant() %>% head()
#' }
#' @importFrom data.table  data.table :=
#' @importFrom magrittr    %>%
#' @export
load_fdata_maxquant <- function(file){

   `Majority protein IDs` <- Reverse <- Contaminant <- NULL

   # Read
   dt <- file %>% autonomics.support::cfread()

   # Establish fvar names
   pattern <- '(?i)(?:potential )?contaminant(?-i)' # Contaminant -> Potential contaminant
   contaminant_var <- names(dt) %>% magrittr::extract(stringi::stri_detect_regex(., pattern))
   fdata_cols <- c('Majority protein IDs', 'Gene names', 'Protein names', 'Reverse', contaminant_var)

   # Process
   dt %>% magrittr::extract(, intersect(names(.), fdata_cols), with = FALSE) %>%

          # Standardize Contaminant varname
          magrittr::set_names(names(.) %>% stringi::stri_replace_first_regex(pattern, 'Contaminant')) %>%

          # Replace NA with empty string
          # when all values are missing, fread reads this in as a numeric column with all NAs
          magrittr::extract(, Reverse     := Reverse     %>% (function(x){ x[is.na(x)] <- ''; x})) %>%
          magrittr::extract(, Contaminant := Contaminant %>% (function(x){ x[is.na(x)] <- ''; x})) %>%

          # Define feature_id
          magrittr::extract(, feature_id := `Majority protein IDs` %>% stringi::stri_split_fixed(';') %>% vapply(extract, character(1), 1)) %>%

          # Convert into dataframe
          data.frame(stringsAsFactors = FALSE, check.names = FALSE, row.names = .$feature_id) %>%
          autonomics.support::pull_columns('feature_id')
}


#============================================
# METABOLON
#============================================

#' Load metabolon fdata
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_fdata_metabolon(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_fdata_metabolon <- function(file, sheet){
   . <- NULL
   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   df %>% magrittr::extract((fstart+1):nrow(.), 1:sstart)                     %>%
      as.data.frame(stringsAsFactors = FALSE)                             %>%
      magrittr::set_names(df[fstart, 1:sstart] %>% unlist() %>% unname()) %>%
      magrittr::set_rownames(paste0('M', .$COMP_ID))                      %>%
      cbind(MCOMP_ID = rownames(.), .)
}


#============================================
# METABOLONLIPIDS
#============================================

#' Load metabolonlipids fdata
#'
#' Load fdata from metabolon clp (complex lipid panel) file
#' @param file path to metabolon clp (complex lipid panel) file
#' @param sheet name of excel sheet (any value in METABOLONLIPIDS_SHEETS)
#' @return feature dataframe
#' @examples
#' require(magrittr)
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_fdata_metabolonlipids('Lipid Class Concentrations') %>% head()
#'    file %>% load_fdata_metabolonlipids(    'Species Concentrations') %>% head()
#'    file %>% load_fdata_metabolonlipids( 'Fatty Acid Concentrations') %>% head()
#'    file %>% load_fdata_metabolonlipids('Lipid Class Compositions')   %>% head()
#'    file %>% load_fdata_metabolonlipids(    'Species Compositions')   %>% head()
#'    file %>% load_fdata_metabolonlipids( 'Fatty Acid Compositions')   %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_fdata_metabolonlipids <- function(file, sheet){
   all_sheets <- readxl::excel_sheets(file) %>% (function(x) x %>% magrittr::set_names(x))

   x <- file %>% readxl::read_excel(sheet = sheet)
   row1 <- which(x[[2]]=='Client Identifier')
   coln <- which(x[row1, ] == 'Unit')
   fnamesrow <- if (stringi::stri_detect_fixed(all_sheets[[sheet]], 'Lipid Class')) row1 else row1-1

   data.frame(feature_id = x %>% magrittr::extract(fnamesrow, (coln+1):ncol(x)) %>% unname() %>% unlist(),
              lipidclass = x %>% magrittr::extract(row1,      (coln+1):ncol(x)) %>% unname() %>% unlist()) %>%
      magrittr::set_rownames(.$feature_id)
}

#============================================
# SOMA
#============================================

#' Load soma fdata
#' @param file string: path to adat file
#' @return feature dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcomp/soma/stemcomp.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::load_fdata_soma() %>% head()
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::load_fdata_soma() %>% head()
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
load_fdata_soma <- function(file){
   x <- file %>% autonomics.import::identify_soma_structure()

   file %>%
      data.table::fread(header = FALSE, sep = '\t', fill = TRUE)                  %>%
      magrittr::extract((x$row-length(x$fvars)-1):(x$row-2),  (x$col-1):ncol(.))  %>%
      t()                                                                         %>%
      data.table::data.table()                                                    %>%
      magrittr::set_names(unlist(unname(.[1,])))                                  %>%
      magrittr::extract(-1, )                                                     %>%
      data.frame(row.names = .$SeqId)
}

#============================================
# GENERIC
#============================================

#' Load fdata
#' @param file path to omics data file
#' @param platform 'maxquant', 'metabolon', 'metabolonlipids', 'soma'
#' @param sheet excel sheet number or name if applicable
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'     file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'     file %>% load_fdata(platform = 'metabolon', sheet = 2) %>% extract(1:3, 1:3)
#'  }
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_fdata(platform = 'metabolonlipids',
#'                        sheet = 'Lipid Class Concentrations') %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_fdata <- function(file, platform, sheet = NULL){
   switch(platform,
          maxquant        = file %>% load_fdata_maxquant(),
          metabolonlipids = file %>% load_fdata_metabolonlipids(sheet = sheet),
          metabolon       = file %>% load_fdata_metabolon(sheet = sheet),
          soma            = file %>% load_fdata_soma())
}

