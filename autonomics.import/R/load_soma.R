
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
   autonomics.import::prepro(object) <- list(assay      = "somascan",
                                           entity     = "epitope",
                                           quantity   = "abundance",
                                           software   = "somalogic",
                                           parameters = parameters)
   autonomics.import::annotation(object) <- parameters$StudyOrganism

   # Return
   return(object)
}

#' Load soma data
#' @param file adat file
#' @param keepOnlyPasses  A logical value indicating whether or not to keep only the
#'                        rows and columns where the data quality was considered to be passable.
#' @param keepOnlySamples A logical value indicating whether or not to keep only the rows
#'                        containing actual samples (as opposed to QC, buffer, and calibrator
#'                        samples).
#' @param verbose         Logical value indicating whether (lots of) diagnostic messages should be shown.
#' @examples
#' require(magrittr)
#' if (require(atkin.2014)){
#'    object <- system.file(
#'                 'extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                  package = 'atkin.2014') %>%
#'              load_soma()
#' }
#' @importFrom magrittr %>%
#' @export
load_soma <- function(file, keepOnlyPasses = TRUE, keepOnlySamples = TRUE, verbose = FALSE){
   object <- readat::readAdat(
                file,
                keepOnlyPasses  = keepOnlyPasses,
                keepOnlySamples = keepOnlySamples) %>%
             autonomics.import::soma_to_sumexp()
   #autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('SeqId',            'feature_id' )
   object
}
