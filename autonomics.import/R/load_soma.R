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


#' Load soma fdata
#' @param file string: path to adat file
#' @return feature dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::load_soma_fdata() %>% head()
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::load_soma_fdata() %>% head()
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
load_soma_fdata <- function(file){
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


#' Load soma sdata
#' @param file string: path to adat file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::load_soma_sdata() %>% head()
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::load_soma_sdata() %>% head()
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
load_soma_sdata <- function(file){
   x <- file %>% autonomics.import::identify_soma_structure()

   file %>%
   data.table::fread(header = FALSE, sep = '\t', fill = TRUE) %>%
   magrittr::extract((x$row-1):nrow(.),  1:(x$col-2), with = FALSE)    %>%
   magrittr::set_names(unlist(unname(.[1,])))                          %>%
   magrittr::extract(-1, )                                             %>%
   data.frame(row.names = .$SampleId)
}


#' Load soma exprs
#' @param file string: path to adat file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::load_soma_exprs()
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::load_soma_exprs() %>% head()
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
load_soma_exprs <- function(file){
   x      <- file %>% autonomics.import::identify_soma_structure()
   fdata1 <- file %>% autonomics.import::load_soma_fdata()
   sdata1 <- file %>% autonomics.import::load_soma_sdata()

   file %>%
   data.table::fread(header = FALSE, sep = '\t', fill = TRUE) %>%
   magrittr::extract(x$row:nrow(.), (x$col):ncol(.))          %>%
   t()                                                        %>%
   data.matrix()                                              %>%
   magrittr::set_rownames(fdata1$SeqId)                       %>%
   magrittr::set_colnames(sdata1$SampleId)                    %>%
  (function(x){class(x) <- 'numeric'; x})
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
#'       file %>% autonomics.import::write_soma_design(design_file = design_file, infer_from_sampleids = TRUE)
#'       file %>% autonomics.import::load_soma(design_file = design_file) %>%
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

   # Assemble components
   fdata1 <- autonomics.import::load_soma_fdata(file)
   sdata1 <- autonomics.import::load_soma_sdata(file)
   exprs1 <- autonomics.import::load_soma_exprs(file)

   # Forge sumexp
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs1))
   autonomics.import::sdata(object) <- sdata1
   autonomics.import::fdata(object) <- fdata1
   autonomics.import::prepro(object) <- list(assay     = "somascan",
                                             entity     = "epitope",
                                             quantity   = "abundance",
                                             software   = "somalogic",
                                             parameters = list())

   # Merge in design
   design_df <- autonomics.import::write_soma_design(file, infer_from_sampleids = infer_design_from_sampleids)
   object %<>% autonomics.import::merge_sdata(design_df, by = 'SampleId')
   if (!is.null(design_file)){
      file_df <- autonomics.import::read_soma_design(design_file)
      object %<>% autonomics.import::merge_sdata(file_df, by = 'SampleId')
   }

   # Preprocess
   if (log2_transform)         autonomics.import::exprs(object) %<>% log2()

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
      autonomics.support::cmessage_df('%s', table(`Sample qualities` = autonomics.import::sdata(object)$RowCheck))
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
      autonomics.support::cmessage_df('%s', table(`Feature qualities` = autonomics.import::fdata(object)$ColCheck))
      if (length(rm_feature_quality) > 0)  object %<>% autonomics.import::filter_features(!ColCheck %in% rm_feature_quality, verbose = TRUE)
   }

   # Select vars
   if (rm_na_svars)            autonomics.import::sdata(object) %<>% autonomics.support::rm_na_columns()
   if (rm_single_value_svars)  autonomics.import::sdata(object) %<>% autonomics.support::rm_single_value_columns()

   # Return
   object

}

#=============================================================================
# The following functions are deprecated, but still kept for now until
# the new function has been thoroughly validated.
#=============================================================================

#' Convert WideSomaLogic object to SummarizedExperiment
#' @param soma WideSomaLogic object
#' @param log2_transform logical
#' @export
soma_to_sumexp <- function (soma, log2_transform = TRUE){

   # Exprs
   myIntensities <- readat::getIntensities(soma, rowsContain = "sequences", reorder = TRUE)
   if (log2_transform) myIntensities %<>% log2()
   #myMeta <- readat::getMetadata(soma)
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = myIntensities))

   # Features
   featureDF <- readat::getSequenceData(soma)
   featureDF <- data.frame(featureDF, row.names = featureDF$SeqId)

   # Samples
   sampleDF <- readat::getSampleData(soma)
   sampleDF <- data.frame(sampleDF, row.names = soma$ExtIdentifier)
   autonomics.import::sdata(object) <- sampleDF
   autonomics.import::fdata(object) <- featureDF

   # Metadata
   parameters <- attributes(soma)$Metadata
   autonomics.import::prepro(object) <- list(assay    = "somascan",
                                           entity     = "epitope",
                                           quantity   = "abundance",
                                           software   = "somalogic",
                                           parameters = parameters)
   if ('StudyOrganism' %in% parameters){
      autonomics.import::annotation(object) <- parameters$StudyOrganism
   }

   # Return
   return(object)
}


#' Old function to load soma data
#'
#' Loads data from soma file.
#' Extracts sample_id definitions from 'SampleId'    column
#' Extracts subgroup definitions  from 'SampleGroup' column when available.
#' When not available, attempts to extract subgroup definitions from 'sample_id'column.
#'
#' @param soma_file              adat file
#' @param sample_file            NULL or character (file with sample design information)
#' @param log2_transform         logical: whether to log2 transform (logical)
#' @param keep_only_passes       logical: whether to keep only features and samples that pass QC
#' @param keep_only_samples      logical: whether to keep only biological samples (and rm QC, buffer, and calibrator samples).
#' @param infer_design           logical: whethe to infer design from sample_ids
#' @param rm_na_svars            logical: whether to rm na svars
#' @param rm_single_value_svars  logical: whether to rm single value svars
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'
#'    # Loading soma is easy
#'    soma_file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                              package = 'autonomics.data')
#'    soma_file %>% autonomics.import::old_load_soma()
#'
#'    # Three ways to specify sample design
#'       # Taken from from 'SampleId' column in soma file
#'       soma_file %>% autonomics.import::old_load_soma() %>% autonomics.import::sdata() %>% head()
#'
#'       # Inferred from sample ids
#'       soma_file %>% autonomics.import::old_load_soma(infer_design = TRUE) %>%
#'                     autonomics.import::sdata() %>% head()
#'
#'       # Specified through sample file
#'       sample_file <- tempfile()
#'       soma_file %>% autonomics.import::create_soma_sample_file(sample_file = sample_file, infer_design = TRUE)
#'       soma_file %>% autonomics.import::old_load_soma(sample_file = sample_file) %>%
#'                     autonomics.import::sdata() %>% head()
#' }
#'
#' @importFrom magrittr %>%
#' @export
old_load_soma <- function(
   soma_file,
   sample_file           = NULL,
   log2_transform        = TRUE,
   keep_only_passes      = TRUE,
   keep_only_samples     = TRUE,
   infer_design          = FALSE,
   rm_na_svars           = TRUE,
   rm_single_value_svars = TRUE
){
   # Read
   object <- suppressWarnings(readat::readAdat(soma_file,
                              keepOnlyPasses  = keep_only_passes,
                              keepOnlySamples = keep_only_samples,
                              verbose = FALSE)) %>%
             autonomics.import::soma_to_sumexp(log2_transform = log2_transform)

   # Add sdata
   sample_df <- if (is.null(sample_file)){
                   autonomics.import::create_soma_sample_df(soma_file, infer_design = infer_design)
                } else {
                   autonomics.support::cfread(sample_file, data.table = FALSE) %>%
                   magrittr::set_rownames(.$SampleId)
                }
   autonomics.import::sdata(object) %<>% (function(x) autonomics.support::left_join_keeping_rownames(sample_df, x, by = 'SampleId'))  %>%
                                          autonomics.support::dedupe_varnames()

   # Cleanup sdata
   autonomics.import::sdata(object) %<>% autonomics.support::rm_na_columns() %>%
                                         autonomics.support::rm_single_value_columns()

   # Return
   object
}

