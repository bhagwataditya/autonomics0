#' Identify soma structure
#'
#' Identify structure of soma file
#'
#' @param file adat file
#' @return list(row, col, fvars, svars)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcomp/soma/stemcomp.adat',
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
#' @param design_sep                   string: sample id separator
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
#'    file <- system.file('extdata/stemcomp/soma/stemcomp.adat',
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
#'       file %>% write_design('soma', infer_design_from_sampleids = TRUE, design_file = design_file)
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
   design_sep                  = NULL,
   log2_transform              = TRUE,
   filter_sample_type          = 'Sample',
   filter_feature_type         = 'Protein',
   filter_sample_quality       = c('FLAG', 'PASS'),
   filter_feature_quality      = c('FLAG', 'PASS'),
   rm_na_svars                 = TRUE,
   rm_single_value_svars       = TRUE
){
   SampleType <- RowCheck <- Type <- ColCheck <- NULL

   # Load sumexp
   object <- autonomics.import::load_omics(file                        = file,
                                           platform                    = 'soma',
                                           log2_transform              = log2_transform,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep)

   autonomics.import::prepro(object) <- list(assay      = "somascan",
                                             entity     = "epitope",
                                             quantity   = "abundance",
                                             software   = "somalogic",
                                             parameters = list())

   # Filter on sample type
   if ('SampleType' %in% autonomics.import::svars(object)){ # older versions don't have it
      message('\t\t========================================================================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Sample types` = autonomics.import::sdata(object)$SampleType))
      object %<>% autonomics.import::filter_samples(SampleType %in% filter_sample_type, verbose = TRUE)
   }
   # Filter on sample quality
   if ('RowCheck'   %in% autonomics.import::svars(object)){
      message('\t\t========================================================================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Sample qualities ("RowCheck")` = autonomics.import::sdata(object)$RowCheck))
      object %<>% autonomics.import::filter_samples(RowCheck %in% filter_sample_quality, verbose = TRUE)
   }
   # Filter on feature type
   if ('Type'       %in% autonomics.import::fvars(object)){
      message('\t\t========================================================================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Type` = autonomics.import::fdata(object)$Type))
      object %<>% autonomics.import::filter_features(Type %in% filter_feature_type, verbose = TRUE)
   }
   # Filter on feature quality
   if ('ColCheck'   %in% autonomics.import::fvars(object)){
      message('\t\t========================================================================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Feature qualities ("ColCheck")` = autonomics.import::fdata(object)$ColCheck))
      object %<>% autonomics.import::filter_features(ColCheck %in% filter_feature_quality, verbose = TRUE)
      message('\t\t========================================================================================================================')
   }

   # Select vars
   if (rm_na_svars)            autonomics.import::sdata(object) %<>% autonomics.support::rm_na_columns()
   if (rm_single_value_svars)  autonomics.import::sdata(object) %<>% autonomics.support::rm_single_value_columns()

   # Return
   object

}

