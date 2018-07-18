
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


#==========================================
# load_sdata
#==========================================

#' Load metabolon sdata
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_metabolon_sdata(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon_sdata <- function(file, sheet){

   # Load sample data
   . <- NULL
   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   df %<>% magrittr::extract(1:fstart, (sstart+1):ncol(.)) %>%
      t() %>%
      data.frame(stringsAsFactors = FALSE, check.names = FALSE)            %>%
      magrittr::set_names(df[[sstart]][1:fstart])     %>%
      magrittr::set_names(names(.) %>% stringi::stri_replace_first_fixed( 'Group   HMDB_ID', 'Group') %>% # recent metabolon files
                             stringi::stri_replace_first_fixed('Sample   HMDB_ID', 'Group'))    # older metabolon files
   df %<>% magrittr::set_rownames(.$CLIENT_IDENTIFIER)

   # Return
   df
}


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


#============================================
# load_fdata
#============================================

#' Load metabolon fdata
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_metabolon_fdata(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon_fdata <- function(file, sheet){
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


#=============================================
# load_exprs
#=============================================

#' Load metabolon exprs
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_metabolon_exprs(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon_exprs <- function(file, sheet){

   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   sdata1 <- file %>% autonomics.import::load_metabolon_sdata(sheet = sheet)
   fdata1 <- file %>% autonomics.import::load_metabolon_fdata(sheet = sheet)

   df %>% magrittr::extract((fstart+1):nrow(.), (sstart+1):ncol(.)) %>%
      data.matrix() %>%
      magrittr::set_colnames(sdata1$CLIENT_IDENTIFIER) %>%
      magrittr::set_rownames(fdata1$MCOMP_ID)
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

#===========================================
# load_omics
#===========================================

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


#' Load metabolon data
#' @param file                metabolon xlsx file
#' @param design_file         NULL or character (sample design file)
#' @param sheet               xls sheet name  or number
#' @param log2_transform      logical: whether to log2 transform
#' @param infer_design_from_sampleids        logical: whether to infer design from sample ids
#' @param add_kegg_pathways   logical: whether to add KEGG pathways to fdata
#' @param add_smiles          logical: whether to add SMILES to fdata
#' @param ... (backward compatibility)
#' @return SummarizedExperiment (load_metabolon) or dataframe (load_metabolon_sdata, load_metabolon_fdata)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'
#'    # Loading metabolon file is easy
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% autonomics.import::load_metabolon()
#'
#'    # Three ways to specify sample design
#'       # Use Group definition in metabolon file
#'       file %>% autonomics.import::load_metabolon() %>%
#'                          autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#'       # Infer from sample id values
#'       file %>% autonomics.import::load_metabolon(infer_design_from_sampleids = TRUE) %>%
#'                          autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#'       # Merge in from sample file
#'       design_file <- tempfile()
#'       file %>% autonomics.import::write_metabolon_design(design_file = design_file, infer_from_sampleids = TRUE)
#'       file %>% autonomics.import::load_metabolon(design_file = design_file) %>%
#'                          autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon <- function(
   file,
   sheet = 2,
   design_file                 = NULL,
   log2_transform              = TRUE,
   infer_design_from_sampleids = FALSE,
   add_kegg_pathways           = FALSE,
   add_smiles                  = FALSE
){
   # Satisfy CHECK
   . <- NULL

   # Sheet
   all_sheets <- readxl::excel_sheets(file)
   cur_sheet <- all_sheets %>% (function(x){ names(x) <- x; x}) %>% magrittr::extract2(sheet)
   autonomics.support::cmessage('Load  %s  %s', basename(file), cur_sheet)

   # Load sumexp
   object <- autonomics.import::load_omics(file                        = file,
                                           sheet                       = sheet,
                                           platform                    = 'metabolon',
                                           log2_transform              = log2_transform,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           sampleid_var                = 'CLIENT_IDENTIFIER')
   autonomics.import::prepro(object) <- list(assay='lcms', entity='metabolite', quantity='intensities', software='metabolon')
   autonomics.import::annotation(object) <- ''

   # Annotate
   if (add_kegg_pathways){
      autonomics.support::cmessage('\tAdd KEGG pathways to fdata')
      object %<>% autonomics.import::add_kegg_pathways_to_fdata()
   }
   if (add_smiles){
      autonomics.support::cmessage('\tAdd SMILES to fdata')
      object %<>% autonomics.import::add_smiles_to_fdata()
   }

   # Return
   object
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

