#============================================
# RNASEQ
#============================================

# Load rnaseq fdata
# @param file path to rnaseq data file
# @return feature dataframe
# @examples
# require(magrittr)
# if (require(subramanian.2016)){
#    file <- system.file('extdata/rnaseq/gene_counts.txt',
#                         package = 'subramanian.2016')
#    file %>% autonomics.import::load_fdata_exiqon()
# }
# @importFrom magrittr %>%
# @export
#load_fdata_rnaseq <- function(
#   file,
#   fvars = c('gene_id', 'locus', 'gene_name', 'gene_type')
#){
#
#   file %>% autonomics.support::cfread()             %>%
#            magrittr::extract(, fvars, with = FALSE) %>%
#            data.frame(stringsAsFactors = FALSE, check.names = FALSE, row.names)
#
#}


#============================================
# EXIQON
#============================================

#' Load exiqon fdata
#' @param file path to exiqon data file
#' @return feature dataframe
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                         package = 'subramanian.2016')
#'    file %>% autonomics.import::load_fdata_exiqon()
#' }
#' @importFrom magrittr %>%
#' @export
load_fdata_exiqon <- function(file){

   Exiqon <- NULL

   file %>% readxl::read_excel()                                        %>%
            dplyr::select(tidyselect::matches('^[^#]'))                 %>%
            dplyr::filter(Exiqon %>% stringi::stri_detect_regex('^#'))  %>%
            t() %>%
            (function(x)
               x[-1, ] %>%
                magrittr::set_colnames(x[1,]) %>%
                data.frame(feature_id = rownames(.), ., stringsAsFactors = FALSE, check.names = FALSE) %>%
                magrittr::set_rownames(.$feature_id)
            )

}


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

   # Satisfy CHECK
   `Majority protein IDs` <- Reverse <- Contaminant <- feature_id <- NULL

   # Read
   dt <- file %>% autonomics.support::cfread()

   # Establish fvar names
   pattern <- '(?i)(?:potential )?contaminant(?-i)' # Contaminant -> Potential contaminant
   contaminant_var <- names(dt) %>% magrittr::extract(stringi::stri_detect_regex(., pattern))
   fdata_cols <- c('Majority protein IDs', 'id', 'Gene names', 'Protein names', 'Reverse', contaminant_var)

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
          data.table::setnames('Majority protein IDs', 'Uniprot accessions') %>%

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
          cbind(feature_id = rownames(.), .)
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
  (function(x){names(x)[1] <- 'feature_id'; x})                                %>%
   data.frame(row.names = .$SeqId)
}

#========
# RNASeq
#========

#' Load rnaseq fdata
#' @param file      path to exiqon xls file
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'  if (require(subramanian.2016)){
#'     file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'     file %>% autonomics.import::load_fdata_rnaseq() %>% head()
#'  }
#' @importFrom magrittr %>%
#' @export
load_fdata_rnaseq <- function(file){
   file %>%
   data.table::fread(header = TRUE) %>%
  (function(x) x %>% magrittr::extract(, vapply(x, is.numeric, logical(1)) %>% magrittr::not(), with = FALSE)) %>%
  (function(x) x %>% data.frame(row.names = x[[1]],
                                stringsAsFactors = FALSE,
                                check.names = FALSE)) %>%
  (function(x){names(x)[1] <- 'feature_id'; x})
}

#============================================
# GENERIC
#============================================

#' Load fdata
#' @param file path to omics data file
#' @param platform 'exiqon', 'maxquant', 'metabolon', 'metabolonlipids', 'soma'
#' @param sheet excel sheet number or name if applicable
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'
#' # EXIQON
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                         package = 'subramanian.2016')
#'    file %>% autonomics.import::load_fdata('exiqon')
#' }
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_fdata(platform = 'metabolon', sheet = 2) %>% extract(1:3, 1:3)
#' }
#'
#' # METABOLONLIPIDS
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_fdata(platform = 'metabolonlipids',
#'                        sheet = 'Lipid Class Concentrations') %>% head()
#' }
#'
#' # RNASEQ
#'  if (require(subramanian.2016)){
#'     file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'     file %>% autonomics.import::load_fdata('rnaseq')
#'  }
#' @importFrom magrittr %>%
#' @export
load_fdata <- function(file, platform, sheet = NULL){
   switch(platform,
          exiqon          = file %>% load_fdata_exiqon(),
          maxquant        = file %>% load_fdata_maxquant(),
          metabolonlipids = file %>% load_fdata_metabolonlipids(sheet = sheet),
          metabolon       = file %>% load_fdata_metabolon(sheet = sheet),
          soma            = file %>% load_fdata_soma(),
          rnaseq          = file %>% load_fdata_rnaseq())
}

