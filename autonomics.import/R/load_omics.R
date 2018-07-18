#=======================================
# sampleid_varname & subgroup_varname
#=======================================

#' Get sampleid svar name
#' @param platform 'metabolonlipids', 'metabolon', 'soma'
#' @return string
#' @examples
#' sampleid_varname('metabolonlipids')
#' sampleid_varname('metabolon')
#' sampleid_varname('soma')
#' @export
sampleid_varname <- function(platform){
   switch(platform,
          metabolonlipids = 'Client Identifier',
          metabolon       = 'CLIENT_IDENTIFIER',
          soma            = 'SampleId')
}

#' Get subgroup svar name
#' @param platform 'metabolonlipids', 'metabolon', 'soma'
#' @return string
#' @examples
#' subgroup_varname('metabolonlipids')
#' subgroup_varname('metabolon')
#' subgroup_varname('soma')
#' @export
subgroup_varname <- function(platform){
   switch(platform,
          metabolonlipids = 'Group',
          metabolon       = 'Group',
          soma            = 'SampleGroup')
}



#==================================================
# identify_structure
#==================================================

#' Identify soma structure
#'
#' Identify structure of soma file
#'
#' @param file adat file
#' @return list(row, col, fvars, svars)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::identify_soma_structure()
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::identify_soma_structure()
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
identify_soma_structure <- function(file){
   DT <- file %>% data.table::fread(header = FALSE, sep = '\t', fill = TRUE)
   fvars <- DT[1 + which(DT$V1 == '^COL_DATA')] %>% unlist() %>% unname() %>% magrittr::extract(.!='') %>% magrittr::extract(2:length(.))
   svars <- DT[1 + which(DT$V1 == '^ROW_DATA')] %>% unlist() %>% unname() %>% magrittr::extract(.!='') %>% magrittr::extract(2:length(.))
   row <- length(fvars) + 2 + which(DT$V1 == '^TABLE_BEGIN')
   col <- length(svars) + 2
   list(row = row, col = col, fvars = fvars, svars = svars)
}


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
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_sdata_metabolon(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_sdata_metabolon <- function(file, sheet){

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
#'    file %>% load_sdata_metabolonlipids('Lipid Class Concentrations') %>% head()
#'    file %>% load_sdata_metabolonlipids(    'Species Concentrations') %>% head()
#'    file %>% load_sdata_metabolonlipids( 'Fatty Acid Concentrations') %>% head()
#'    file %>% load_sdata_metabolonlipids('Lipid Class Compositions')   %>% head()
#'    file %>% load_sdata_metabolonlipids(    'Species Compositions')   %>% head()
#'    file %>% load_sdata_metabolonlipids( 'Fatty Acid Compositions')   %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_sdata_metabolonlipids <- function(file, sheet){
   x <- file %>% readxl::read_excel(sheet = sheet)
   row1 <- which(x[[2]]=='Client Identifier')
   coln <- which(x[row1, ] == 'Unit')

   x %>% magrittr::extract((1+row1):nrow(x), 1:coln) %>%
      data.frame() %>%
      magrittr::set_names(x[row1, 1:coln] %>% unname() %>% unlist()) %>%
      magrittr::set_rownames(.$`Client Identifier`)
}


#' Load soma sdata
#' @param file string: path to adat file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::load_sdata_soma() %>% head()
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::load_sdata_soma() %>% head()
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
load_sdata_soma <- function(file){
   x <- file %>% autonomics.import::identify_soma_structure()

   file %>%
      data.table::fread(header = FALSE, sep = '\t', fill = TRUE) %>%
      magrittr::extract((x$row-1):nrow(.),  1:(x$col-2), with = FALSE)    %>%
      magrittr::set_names(unlist(unname(.[1,])))                          %>%
      magrittr::extract(-1, )                                             %>%
      data.frame(row.names = .$SampleId)
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
          metabolonlipids = load_sdata_metabolonlipids(file = file, sheet = sheet),
          metabolon       = load_sdata_metabolon(file = file, sheet = sheet),
          soma            = load_sdata_soma(file))
}


#============================================
# load_fdata
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


#' Load soma fdata
#' @param file string: path to adat file
#' @return feature dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
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
          metabolonlipids = load_fdata_metabolonlipids(file = file, sheet = sheet),
          metabolon       = load_fdata_metabolon(      file = file, sheet = sheet),
          soma            = load_fdata_soma(file))
}


#=============================================
# load_exprs
#=============================================

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
          metabolonlipids = load_exprs_metabolonlipids(file = file, sheet = sheet),
          metabolon       = load_exprs_metabolon(      file = file, sheet = sheet),
          soma            = load_exprs_soma(file))
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
#'     file %>% load_omics(sheet=2, platform = 'metabolon')
#'  }
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_omics(sheet = 'Lipid Class Concentrations', platform = 'metabolonlipids')
#' }
#' @importFrom magrittr %>%
#' @export
load_omics <- function(
   file,
   sheet = 2,
   platform,
   log2_transform              = TRUE,
   design_file                 = NULL,
   infer_design_from_sampleids = FALSE
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
   design_df <- autonomics.import::write_design(file, platform = platform, infer_from_sampleids = infer_design_from_sampleids, sheet = sheet)
   object %<>% autonomics.import::merge_sdata(design_df, by = sampleid_varname(platform))
   if (!is.null(design_file)){
      file_df <- autonomics.import::read_design(design_file)
      object %<>% autonomics.import::merge_sdata(file_df, by = sampleid_varname(platform))
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
#' @return SummarizedExperiment (load_metabolon) or dataframe (load_sdata_metabolon, load_fdata_metabolon)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'
#'    # Loading metabolon file is easy
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_metabolon()
#'
#'    # Three ways to specify sample design
#'       # Use Group definition in metabolon file
#'       file %>% load_metabolon() %>%
#'                          autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#'       # Infer from sample id values
#'       file %>% load_metabolon(infer_design_from_sampleids = TRUE) %>%
#'                autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#'       # Merge in from sample file
#'       design_file <- tempfile()
#'       file %>% write_design('metabolon', infer_from_sampleids = TRUE, design_file = design_file)
#'       file %>% load_metabolon(design_file = design_file) %>%
#'                autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
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
                                           infer_design_from_sampleids = infer_design_from_sampleids)
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
#' @param file            path to metabolon lipids file
#' @param sheet           name of excel sheet (any value in METABOLONLIPIDS_SHEETS)
#' @param log2_transform  logical: whether to log2 transform
#' @param design_file     path to sample design file
#' @param infer_design_from_sampleids  logical: whether to infer design from sampleids
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
                                           infer_design_from_sampleids = infer_design_from_sampleids)

   # Return
   object
}



#' Load soma
#'
#' Loads data from soma file.
#' Extracts sample_id definitions from 'SampleId'    column
#' Extracts subgroup definitions  from 'SampleGroup' column when available.
#' When not available, attempts to extract subgroup definitions from 'sample_id'column.
#'
#' @param file                         string: path to adat file
#' @param design_file                  NULL or string:  path to sample file
#' @param infer_design_from_sampleids  logical: whether to infer design from SampleId values
#' @param log2_transform               logical: whether to log2 transform
#' @param rm_sample_type               character vector: sample  types to be removed. Probably a subset of c('Sample', 'QC', 'Buffer', 'Calibrator').
#' @param rm_feature_type              character vector: feature types to be removed. Probably a subset of c('Protein', 'Hybridization Control Elution', 'Rat Protein').
#' @param rm_sample_quality            character vector: sample  qualities to be removed. Probably a subset of c('PASS', 'FLAG', 'FAIL')
#' @param rm_feature_quality           character vector: feature qualities to be removed. Probably a subset of c('PASS', 'FLAG', 'FAIL')
#' @param rm_na_svars                  logical: whether to rm NA svars
#' @param rm_single_value_svars        logical: whether to rm single value svars
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'
#'    # Loading soma file is simple
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::load_soma()
#'
#'    # Three ways to specify sample design
#'       # Taken from 'SampleId' column in soma file
#'       file %>% autonomics.import::load_soma() %>% autonomics.import::sdata() %>% head()
#'
#'       # Inferred from sample ids
#'       file %>% autonomics.import::load_soma(infer_design = TRUE) %>%
#'                autonomics.import::sdata() %>% head()
#'
#'       # Specified through sample file
#'       design_file <- tempfile()
#'       file %>% write_design('soma', infer_from_sampleids = TRUE, design_file = design_file)
#'       file %>% load_soma(design_file = design_file) %>%
#'                autonomics.import::sdata() %>% head()
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::load_soma()
#'    file %>% autonomics.import::load_soma(rm_feature_quality = 'FAIL')
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
load_soma <- function(
   file,
   design_file                 = NULL,
   infer_design_from_sampleids = FALSE,
   log2_transform              = TRUE,
   rm_sample_type              = character(0),
   rm_feature_type             = character(0),
   rm_sample_quality           = character(0),
   rm_feature_quality          = character(0),
   rm_na_svars                 = TRUE,
   rm_single_value_svars       = TRUE
){
   SampleType <- RowCheck <- Type <- ColCheck <- NULL

   # Load sumexp
   object <- load_omics(file                        = file,
                        platform                    = 'soma',
                        log2_transform              = log2_transform,
                        design_file                 = design_file,
                        infer_design_from_sampleids = infer_design_from_sampleids)

   autonomics.import::prepro(object) <- list(assay      = "somascan",
                                             entity     = "epitope",
                                             quantity   = "abundance",
                                             software   = "somalogic",
                                             parameters = list())
   # Filter
   # On Sample Type
   if ('\nSampleType' %in% autonomics.import::svars(object)){
      message('')
      autonomics.support::cmessage_df('%s', table(`Sample types` = autonomics.import::sdata(object)$SampleType))
      if (length(rm_sample_type) > 0)      object %<>% autonomics.import::filter_samples(!SampleType %in% rm_sample_type, verbose = TRUE)
   }
   # On sample quality
   if ('RowCheck'   %in% autonomics.import::svars(object)){
      message('')
      autonomics.support::cmessage_df('%s', table(`Sample qualities ("RowCheck")` = autonomics.import::sdata(object)$RowCheck))
      if (length(rm_sample_quality)  > 0)  object %<>% autonomics.import::filter_samples(!RowCheck %in% rm_sample_quality, verbose = TRUE)
   }
   # On feature type
   if ('Type'       %in% autonomics.import::fvars(object)){
      message('')
      autonomics.support::cmessage_df('%s', table(`Type` = autonomics.import::fdata(object)$Type))
      if (length(rm_feature_type)    > 0)  object %<>% autonomics.import::filter_features(!Type %in% rm_feature_type, verbose = TRUE)
   }
   # On feature quality
   if ('ColCheck'   %in% autonomics.import::fvars(object)){
      message('')
      autonomics.support::cmessage_df('%s', table(`Feature qualities ("ColCheck")` = autonomics.import::fdata(object)$ColCheck))
      if (length(rm_feature_quality) > 0)  object %<>% autonomics.import::filter_features(!ColCheck %in% rm_feature_quality, verbose = TRUE)
   }

   # Select vars
   if (rm_na_svars)            autonomics.import::sdata(object) %<>% autonomics.support::rm_na_columns()
   if (rm_single_value_svars)  autonomics.import::sdata(object) %<>% autonomics.support::rm_single_value_columns()

   # Return
   object

}
