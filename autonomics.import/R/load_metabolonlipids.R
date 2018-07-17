#' Possible sheet values for metabolonlipids data
#' @export
METABOLONLIPIDS_SHEETS <- c('Lipid Class Concentrations',
                            'Lipid Class Compositions',
                            'Species Concentrations',
                            'Species Compositions',
                            'Fatty Acid Concentrations',
                            'Fatty Acid Compositions')


#' Load metabolonlipids sdata
#'
#' Load sdata from metabolon clp (complex lipid panel) file
#' @param file  path to metabolon clp (complex lipid panel) file
#' @param sheet name of excel sheet (any value in METABOLONLIPIDS_SHEETS)
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_metabolonlipids_sdata('Lipid Class Concentrations') %>% head()
#'    file %>% load_metabolonlipids_sdata(    'Species Concentrations') %>% head()
#'    file %>% load_metabolonlipids_sdata( 'Fatty Acid Concentrations') %>% head()
#'    file %>% load_metabolonlipids_sdata('Lipid Class Compositions')   %>% head()
#'    file %>% load_metabolonlipids_sdata(    'Species Compositions')   %>% head()
#'    file %>% load_metabolonlipids_sdata( 'Fatty Acid Compositions')   %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolonlipids_sdata <- function(file, sheet){
   x <- file %>% readxl::read_excel(sheet = sheet)
   row1 <- which(x[[2]]=='Client Identifier')
   coln <- which(x[row1, ] == 'Unit')

   x %>% magrittr::extract((1+row1):nrow(x), 1:coln) %>%
         data.frame() %>%
         magrittr::set_names(x[row1, 1:coln] %>% unname() %>% unlist()) %>%
         magrittr::set_rownames(.$`Client Identifier`)
}


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
#'    file %>% load_metabolonlipids_fdata('Lipid Class Concentrations') %>% head()
#'    file %>% load_metabolonlipids_fdata(    'Species Concentrations') %>% head()
#'    file %>% load_metabolonlipids_fdata( 'Fatty Acid Concentrations') %>% head()
#'    file %>% load_metabolonlipids_fdata('Lipid Class Compositions')   %>% head()
#'    file %>% load_metabolonlipids_fdata(    'Species Compositions')   %>% head()
#'    file %>% load_metabolonlipids_fdata( 'Fatty Acid Compositions')   %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolonlipids_fdata <- function(file, sheet){
   all_sheets <- readxl::excel_sheets(file) %>% (function(x) x %>% magrittr::set_names(x))

   x <- file %>% readxl::read_excel(sheet = sheet)
   row1 <- which(x[[2]]=='Client Identifier')
   coln <- which(x[row1, ] == 'Unit')
   fnamesrow <- if (stringi::stri_detect_fixed(all_sheets[[sheet]], 'Lipid Class')) row1 else row1-1

   data.frame(feature_id = x %>% magrittr::extract(fnamesrow, (coln+1):ncol(x)) %>% unname() %>% unlist(),
              lipidclass = x %>% magrittr::extract(row1,      (coln+1):ncol(x)) %>% unname() %>% unlist()) %>%
                                 magrittr::set_rownames(.$feature_id)
}


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
#'    file %>% load_metabolonlipids_exprs('Lipid Class Concentrations') %>% extract(1:3,1:3)
#'    file %>% load_metabolonlipids_exprs(    'Species Concentrations') %>% extract(1:3,1:3)
#'    file %>% load_metabolonlipids_exprs( 'Fatty Acid Concentrations') %>% extract(1:3,1:3)
#'    file %>% load_metabolonlipids_exprs('Lipid Class Compositions')   %>% extract(1:3,1:3)
#'    file %>% load_metabolonlipids_exprs(    'Species Compositions')   %>% extract(1:3,1:3)
#'    file %>% load_metabolonlipids_exprs( 'Fatty Acid Compositions')   %>% extract(1:3,1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolonlipids_exprs <- function(file, sheet){
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


#' Load metabolonlipids
#'
#' Load data from metabolon complex lipid panel (clp) file
#'
#' @param file         path to metabolon lipids file
#' @param sheet name of excel sheet (any value in METABOLONLIPIDS_SHEETS)
#' @param design_file  path to sample design file
#' @param infer_from_sampleids logical
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% autonomics.import::load_metabolonlipids('Lipid Class Concentrations')
#'    file %>% autonomics.import::load_metabolonlipids('Lipid Class Compositions')
#'    file %>% autonomics.import::load_metabolonlipids(    'Species Concentrations')
#'    file %>% autonomics.import::load_metabolonlipids(    'Species Compositions')
#'    file %>% autonomics.import::load_metabolonlipids( 'Fatty Acid Concentrations')
#'    file %>% autonomics.import::load_metabolonlipids( 'Fatty Acid Compositions')
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolonlipids <- function(
   file,
   sheet,
   log2_transform = TRUE,
   design_file = NULL,
   infer_design_from_sampleids = FALSE
){

   # Load and Create
   object <- autonomics.import::load_omics(file                        = file,
                                           sheet                       = sheet,
                                           platform                    = 'metabolonlipids',
                                           log2_transform              = log2_transform,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           sampleid_var                = 'Client Identifier')

   # Return
   object
}

