
#' Validify sample/feature names
#' @param x character vector with sample ids
#' @return character vector
#' @importFrom magrittr %<>%
#' @export
validify_sample_ids <- function(x){
   . <- NULL
   selector <- substr(x,1,1) %in% 0:9
   x[selector] %<>% paste0('s', .)
   x
}


#' Load sdata
#' @param file path to omics data file
#' @param sheet excel sheet number or name if applicable
#' @param platform 'metabolon', 'metabolonlipids'
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'     file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'     file %>% load_sdata(2, 'metabolon') %>% extract(1:3, 1:3)
#'  }
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_sdata('Lipid Class Concentrations', 'metabolonlipids') %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_sdata <- function(file, sheet = NULL, platform){
   switch(platform,
          metabolonlipids = load_metabolonlipids_sdata(file = file, sheet = sheet),
          metabolon       = load_metabolon_sdata(file = file, sheet = sheet))
}


#' Load fdata
#' @param file path to omics data file
#' @param sheet excel sheet number or name if applicable
#' @param platform 'metabolon', 'metabolonlipids'
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'     file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'     file %>% load_fdata(2, 'metabolon') %>% extract(1:3, 1:3)
#'  }
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_fdata('Lipid Class Concentrations', 'metabolonlipids') %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_fdata <- function(file, sheet = NULL, platform){
   switch(platform,
          metabolonlipids = load_metabolonlipids_fdata(file = file, sheet = sheet),
          metabolon       = load_metabolon_fdata(      file = file, sheet = sheet))
}


#' Load exprs
#' @param file path to omics data file
#' @param sheet excel sheet number or name if applicable
#' @param platform 'metabolon', 'metabolonlipids'
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'     file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'     file %>% load_exprs(2, 'metabolon') %>% extract(1:3, 1:3)
#'  }
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_exprs('Lipid Class Concentrations', 'metabolonlipids') %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_exprs <- function(file, sheet = NULL, platform){
   switch(platform,
          metabolonlipids = load_metabolonlipids_exprs(file = file, sheet = sheet),
          metabolon       = load_metabolon_exprs(      file = file, sheet = sheet))
}


#' Load omics
#' @param file path to omics data file
#' @param sheet excel sheet number or name if applicable
#' @param platform 'metabolon', 'metabolonlipids'
#' @param design_file path to design file
#' @param log2_transform logical
#' @param infer_design_from_sampleids logical
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'     file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'     file %>% load_omics(2, 'metabolon')
#'  }
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_omics('Lipid Class Concentrations', 'metabolonlipids')
#' }
#' @importFrom magrittr %>%
#' @export
load_omics <- function(
   file,
   sheet = 2,
   platform,
   log2_transform              = TRUE,
   design_file                 = NULL,
   infer_design_from_sampleids = FALSE,
   sampleid_var
){
   # Satisfy CHECK
   . <- NULL

   # Load components
   sdata1 <- file %>% autonomics.import::load_sdata(sheet = sheet, platform = platform)
   fdata1 <- file %>% autonomics.import::load_fdata(sheet = sheet, platform = platform)
   exprs1 <- file %>% autonomics.import::load_exprs(sheet = sheet, platform = platform)

   # Wrap into SummarizedExperiment
   object <- SummarizedExperiment::SummarizedExperiment(assays=list(exprs = exprs1))
   if (log2_transform)   autonomics.import::exprs(object) %<>% log2()
   autonomics.import::sdata(object)  <- sdata1
   autonomics.import::fdata(object)  <- fdata1

   # Merge in design
   design_df <- autonomics.import::write_metabolon_design(file, sheet = sheet, infer_from_sampleids = infer_design_from_sampleids)
   object %<>% autonomics.import::merge_sdata(design_df, by = sampleid_var)
   if (!is.null(design_file)){
      file_df <- autonomics.import::read_metabolon_design(design_file)
      object %<>% autonomics.import::merge_sdata(file_df, by = sampleid_var)
   }

   # Return
   object
}
