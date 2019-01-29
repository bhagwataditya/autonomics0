#=============================
# MAXQUANT
#=============================

#' Load proteingroups
#' @param file                          path to proteinGroups.txt
#' @param quantity                     'Ratio normalized', 'Ratio', 'Intensity', 'LFQ intensity', 'Reporter intensity'
#' @param infer_design_from_sampleids   logical: whether to infer design from sampleids
#' @param design_sep                    string: design separator
#' @param design_file                   path to design file (created with write_maxquant_design)
#' @param fasta_file                    path to uniprot fasta database
#' @param log2_transform                logical: whether to log2 transform
#' @param rm_reverse_features           logical: whether to rm reverse features
#' @param rm_contaminant_features       logical: whether to rm contaminant features
#' @param rm_na_features                logical: whether to rm features which are NA, NaN or 0 in all samples
#' @param impute_consistent_nondetects: logical(1)
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcomp/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'
#'    # The sample design can be inferred from the sample ids
#'       object <- file %>% autonomics.import::load_proteingroups()
#'       object %>% autonomics.import::sdata()
#'
#'    # Or it can be loaded through a sample file
#'       design_file <- tempfile()
#'       write_design(file, platform = 'maxquant', infer_design_from_sampleids = TRUE,
#'                    design_file = design_file)
#'       object <- load_proteingroups(file, design_file = design_file)
#'       sdata(object)
#' }
#'
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    object <- file %>% autonomics.import::load_proteingroups()
#'    sdata(object) %>% str()
#' }
#'
#' # LFQ
#' if (require(graumann.lfq)){
#'    file <- 'extdata/proteinGroups.txt' %>% system.file(package = 'graumann.lfq')
#'    file %>% autonomics.import::load_proteingroups()
#' }
#'
#' @importFrom magrittr %>%
#' @export
load_proteingroups <- function(
   file,
   quantity                     = autonomics.import::infer_maxquant_quantity(file),
   infer_design_from_sampleids  = TRUE,
   design_sep                   = NULL,
   design_file                  = NULL,
   fasta_file                   = NULL,
   log2_transform               = TRUE,
   rm_reverse_features          = TRUE,
   rm_contaminant_features      = TRUE,
   rm_na_features               = TRUE,
   impute_consistent_nondetects = quantity %>% stringi::stri_detect_fixed('intensity')
){
   # Satisfy CHECK
   Reverse <- Contaminant <- feature_id <- NULL

   # Load
   object <- autonomics.import::load_omics(file                        = file,
                                           platform                    = 'maxquant',
                                           quantity                    = quantity,
                                           log2_transform              = FALSE,
                                           log2_offset                 = 0,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep,
                                           design_file                 = design_file)

   # NA zeroes (in intensity data)
   autonomics.support::cmessage('\t\tNA zeroes')
   autonomics.import::exprs(object) %<>% (function(x){x[x==0] <- NA; x})

   # Log2 transform
   if (log2_transform){
      autonomics.support::cmessage('\t\tLog2 transform')
      autonomics.import::exprs(object) %<>% log2()
   }

   # Filter features
   object %<>% autonomics.import::filter_features(!is.na(feature_id), verbose = TRUE) # Max Quant earlier version had bug that created corrupted lines without feature_id columns
   if (rm_reverse_features){
      object %<>% autonomics.import::filter_features(Reverse != '+',     verbose = TRUE)
      autonomics.import::fdata(object)$Reverse     <- NULL
   }
   if (rm_contaminant_features){
      object %<>% autonomics.import::filter_features(Contaminant != '+', verbose = TRUE)
      autonomics.import::fdata(object)$Contaminant <- NULL
   }
   if (rm_na_features){
      object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample(verbose = TRUE)
   }

   # Impute consistent nondetects (in intensity data)
   if (impute_consistent_nondetects) object %<>% autonomics.preprocess::qrilc_consistent_nondetects()

   # Intuify snames
   subgroup_values  <- object %>% autonomics.import::svalues('subgroup')
   replicate_values <- object %>% autonomics.import::svalues('replicate')
   # If not specified as argument: infer sep from sampleids
   design_sep %<>% (function(x) if (is.null(x)) object %>% autonomics.import::snames() %>% autonomics.import::guess_sep() else x) %>%
      # If not inferrable from sampleids: set sep to '.'
      (function(x) if (is.null(x)) '.' else x)
   ok <- !is.null( subgroup_values) &
      !is.null(replicate_values) &
      all(assertive.strings::is_non_missing_nor_empty_character(as.character(subgroup_values)))  &
      all(assertive.strings::is_non_missing_nor_empty_character(as.character(replicate_values))) &
      assertive.properties::has_no_duplicates(paste0(subgroup_values, '.', replicate_values))
   if (ok){
      new_sampleids <- paste0(subgroup_values, '.', replicate_values)
      autonomics.import::snames(object) <- new_sampleids
      autonomics.import::sdata(object)$sample_id <- new_sampleids
   }

   # Add prepro
   autonomics.import::prepro(object) <- list(assay    = 'lcms',
                                             entity   = 'proteingroup',
                                             quantity = quantity,
                                             software = 'maxquant')

   # Return
   return(object)

}

#' Open uniprot connection
#' @param object SummarizedExperiment
#' @return uniprot webservice connection
#' @examples
#' \dontrun{
#'    require(magrittr)
#'    if (require(autonomics.data)){
#'       object <- autonomics.data::stemcomp.proteinratios
#'       object %>% autonomics.import::open_uniprot_connection()
#'    }
#' }
#' @importFrom magrittr %>%
#' @export
open_uniprot_connection <- function(object){
   object                                                %>%
      autonomics.import::uniprot_values(first_only = TRUE)  %>%
      magrittr::extract(1:10)                               %>%
      autonomics.annotate::connect_to_uniprot()
}

#' Annotate proteingroups through uniprot.ws
#' @param object      SummarizedExperiment
#' @param connection  uniprot webservice connection
#' @param columns     uniprot webservice columns
#' @return Annotated SummarizedExperiment
#' @examples
#' require(magrittr)
#' \dontrun{
#'    if (require(autonomics.data)){
#'       object <- autonomics.data::stemcomp.proteinratios
#'       connection <- autonomics.import::open_uniprot_connection(object)
#'       object[1:10, ] %>% annotate_proteingroups(connection, c('SUBCELLULAR-LOCATIONS', 'GO-ID'))
#'    }
#' }
#' @importFrom magrittr %>%  %<>%
#' @export
annotate_proteingroups <- function(
   object,
   connection = object %>% autonomics.import::open_uniprot_connection(),
   columns    = c('SUBCELLULAR-LOCATIONS', 'INTERPRO', 'GO-ID')
){
   # Assert
   assertive.base::assert_is_identical_to_true(class(object) == 'SummarizedExperiment')

   # Restrict to first uniprot accession
   autonomics.import::fdata(object)$`Uniprot accessions` %<>% stringi::stri_split_fixed(';') %>% vapply(extract, character(1), 1)

   # Fetch annotations from uniprot
   annotations <- autonomics.import::fdata(object)$`Uniprot accessions` %>%
      autonomics.annotate::annotate_uniprot_with_webservice(connection = connection, columns = columns)

   # Merge in annotations
   autonomics.import::fdata(object) %<>% merge(annotations, by.x = 'Uniprot accessions', by.y = 'UNIPROTKB', sort = FALSE)

   # Return
   return(object)

}

