



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

