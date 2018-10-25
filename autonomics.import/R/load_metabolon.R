
#======================
# METABOLON
#======================

#' Load metabolon data
#' @param file                metabolon xlsx file
#' @param design_file         NULL or character (sample design file)
#' @param sheet               xls sheet name  or number
#' @param log2_transform      logical: whether to log2 transform
#' @param log2_offset         offset in mapping x -> log2(offset + x)
#' @param infer_design_from_sampleids        logical: whether to infer design from sample ids
#' @param design_sep          string: sample id separator
#' @param add_kegg_pathways   logical: whether to add KEGG pathways to fdata
#' @param add_smiles          logical: whether to add SMILES to fdata
#' @return SummarizedExperiment
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
#'       file %>% write_design('metabolon', infer_design_from_sampleids = TRUE,
#'                              design_file = design_file)
#'       file %>% load_metabolon(design_file = design_file) %>%
#'                autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#' }
#' if (require(subramanian.2016)){
#'    file <- 'extdata/metabolon/subramanian.2016.metabolon.xlsx' %>%
#'            system.file(package='subramanian.2016')
#'    object <- file %>% autonomics.import::load_metabolon(sheet = 5)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon <- function(
   file,
   sheet = 2,
   design_file                 = NULL,
   log2_transform              = TRUE,
   log2_offset                 = 0,
   infer_design_from_sampleids = FALSE,
   design_sep                  = NULL,
   add_kegg_pathways           = FALSE,
   add_smiles                  = FALSE
){
   # Satisfy CHECK
   . <- NULL

   # Sheet
   all_sheets <- readxl::excel_sheets(file)
   cur_sheet <- all_sheets %>% (function(x){ names(x) <- x; x}) %>% magrittr::extract2(sheet)
   autonomics.support::cmessage('\t\tLoad  %s  %s', basename(file), cur_sheet)

   # Load sumexp
   object <- autonomics.import::load_omics(file                        = file,
                                           sheet                       = sheet,
                                           platform                    = 'metabolon',
                                           log2_transform              = log2_transform,
                                           log2_offset                 = log2_offset,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep)
   autonomics.import::prepro(object) <- list(assay='lcms', entity='metabolite', quantity='intensities', software='metabolon')
   autonomics.import::annotation(object) <- ''

   # Zero consistent nas
   #object %<>% autonomics.preprocess::zero_consistent_nas(verbose = TRUE)

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


#============================
# METABOLONLIPIDS
#============================

#' Load metabolonlipids
#'
#' Load data from metabolon complex lipid panel (clp) file
#'
#' @param file            path to metabolon lipids file
#' @param sheet           name of excel sheet (any value in METABOLONLIPIDS_SHEETS)
#' @param log2_transform  logical: whether to log2 transform
#' @param log2_offset     numeric: offset in mapping x -> log2(offset + x)
#' @param design_file     path to sample design file
#' @param infer_design_from_sampleids  logical: whether to infer design from sampleids
#' @param design_sep      string: sample id separator
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
   log2_offset    = 0,
   design_file = NULL,
   infer_design_from_sampleids = FALSE,
   design_sep = NULL
){

   # Load and Create
   object <- autonomics.import::load_omics(file                        = file,
                                           sheet                       = sheet,
                                           platform                    = 'metabolonlipids',
                                           log2_transform              = log2_transform,
                                           log2_offset                 = log2_offset,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep)

   # Return
   object
}


