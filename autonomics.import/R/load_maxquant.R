#=============================
# MAXQUANT
#=============================

#' Load proteingroups
#' @param file            path to proteinGroups.txt
#' @param quantity       'Ratio normalized', 'Ratio', 'Intensity', 'LFQ intensity', 'Reporter intensity'
#' @param infer_design_from_sampleids  logical: whether to infer design from sampleids
#' @param design_sep      string: design separator
#' @param design_file     path to design file (created with write_maxquant_design)
#' @param fasta_file      path to uniprot fasta database
#' @param log2_transform  logical: whether to log2 transform
#' @param log2_offset     numeric: offset used in mapping: x -> log2(offset + x)
#' @param rm_reverse_features      logical: whether to rm reverse features
#' @param rm_contaminant_features logical: whether to rm contaminant features
#' @param rm_na_features  logical: whether to rm features which are NA, NaN or 0 in all samples
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcomp/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'
#'    # The sample design can be inferred from the sample ids
#'       object <- file %>% load_proteingroups(infer_design_from_sampleids = TRUE)
#'       object %>% sdata()
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
#'    object <- file %>% autonomics.import::load_proteingroups(
#'                          infer_design_from_sampleids = TRUE)
#'    sdata(object) %>% str()
#' }
#'
#' @importFrom magrittr %>%
#' @export
load_proteingroups <- function(
   file,
   quantity                    = autonomics.import::infer_maxquant_quantity(file),
   infer_design_from_sampleids = FALSE,
   design_sep                  = NULL,
   design_file                 = NULL,
   fasta_file                  = NULL,
   log2_transform              = TRUE,
   log2_offset                 = 0,
   rm_reverse_features         = TRUE,
   rm_contaminant_features     = TRUE,
   rm_na_features              = TRUE
){
   # Satisfy CHECK
   Reverse <- Contaminant <- feature_id <- NULL

   # Load
   object <- autonomics.import::load_omics(file                        = file,
                                           platform                    = 'maxquant',
                                           quantity                    = quantity,
                                           log2_transform              = log2_transform,
                                           log2_offset                 = log2_offset,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep,
                                           design_file                 = design_file)

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

   # Intuify snames
   subgroup_values  <- object %>% autonomics.import::svalues('subgroup')
   replicate_values <- object %>% autonomics.import::svalues('replicate')
   # If not specified as argument: infer sep from sampleids
   design_sep %<>% (function(x) if (is.null(x)) object %>% autonomics.import::snames() %>% autonomics.import::infer_design_sep() else x) %>%
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

#==============================================================

#' Deconvolute proteingroups
#' @param object             SummerizedExperiment with proteinGroups data
#' @param fastafile          path to fastafile
#' @param fastafields        character vector: fields to load from fastafile
#' @param drop_isoform_info  logical: whether to drop isoform info
#' @return deconvoluted and annotated SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    fastafile <- '../data/uniprot_hsa_20140515.fasta'
#'    if (file.exists(fastafile)){
#'       object <- 'extdata/stemcomp/maxquant/proteinGroups.txt'            %>%
#'                  system.file(package='autonomics.data')                  %>%
#'                  load_proteingroups(infer_design_from_sampleids = TRUE)
#'       autonomics.import::fdata(object) %>% head()
#'
#'       object %<>% magrittr::extract(1:100, )                             %>%
#'                   deconvolute_proteingroups(fastafile = fastafile)
#'       autonomics.import::fdata(object) %>% head()
#'    }
#' }
#' @importFrom magrittr %>% %<>%
#' @export
deconvolute_proteingroups <- function(
   object,
   fastafile,
   fastafields = c('GENES', 'PROTEIN-NAMES', 'EXISTENCE', 'REVIEWED'),
   drop_isoform_info = FALSE
){
   # Satisfy CHECK
   EXISTENCE <- GENES <- IS.FRAGMENT <- ISOFORM <- N <- NGENE <- NISOFORMS <- NPERACCESSION <- NULL
   `PROTEIN-NAMES` <- REVIEWED <- `Uniprot accessions` <- ngene <- nprotein <- nseq <- .SD <- NULL

   # Load fasta annotations
   autonomics.support::cmessage('Load fasta file')
   fasta_annotations <- fastafile %>% autonomics.annotate::load_uniprot_fasta_annotations(fastafields)

   # Uncollapse
   fdata1 <- object %>%
      autonomics.import::fdata() %>%
      magrittr::extract(, c('feature_id', 'Uniprot accessions'), drop = FALSE) %>%
      tidyr::separate_rows(`Uniprot accessions`, sep = ';') %>%
      data.table::data.table()

   # Split into CANONICAL and isoform
   fdata1 %<>% magrittr::extract(, ISOFORM              := `Uniprot accessions`) %>%
      magrittr::extract(, `Uniprot accessions` := `Uniprot accessions` %>% stringi::stri_replace_first_regex('[-][0-9]+',     '')) # %>%
   #magrittr::extract(, ISOFORM := ISOFORM %>% sort() %>% unique() %>% paste0(collapse=';'), by = c('feature_id', 'Uniprot accessions')) %>%
   #unique()

   # Merge in uniprot fasta annotations
   nunmapped <- fdata1 %>%
      magrittr::extract(ISOFORM %in% setdiff(ISOFORM, fasta_annotations$UNIPROTKB)) %>%
      magrittr::extract(, .SD[1], by = 'feature_id') %>%
      nrow()
   autonomics.support::cmessage('\nFastafile misses sequences of %d/%d proteingroups',
                                nunmapped, nrow(object))
   fdata1 %<>% merge(fasta_annotations, by.x = 'Uniprot accessions', by.y = 'UNIPROTKB', sort = FALSE)
   report_n <- function(dt, prefix='', suffix=''){
      n <- dt %>% magrittr::extract(, .SD[, list(nseq = .N,
                                                 ngene = length(unique(GENES)),
                                                 nprotein = length(unique(`PROTEIN-NAMES`)))],
                                    by = 'feature_id') %>%
         magrittr::extract(, list(ngroups = .N,
                                  nsinglegene    = sum(ngene==1),
                                  nsingleprotein = sum(nprotein==1),
                                  nsingleseq     = sum(nseq==1)))
      autonomics.support::cmessage('%s%d proteingroups -> %d singlegene -> %d singleprotein -> %d singleseq%s',
                                   stringi::stri_pad_right(prefix, width = 60), n$ngroups,          n$nsinglegene,   n$nsingleprotein,   n$nsingleseq, suffix)
   }
   message('')
   fdata1 %>% report_n(prefix = 'All proteingroups')

   # Prefer best existence
   fdata1 %<>% magrittr::extract(, .SD[EXISTENCE == min(EXISTENCE)], by = 'feature_id')
   fdata1 %>% report_n(prefix = 'Per group: drop inferior existences')

   # Drop trembl entries from swissprot groups
   fdata1 %<>% magrittr::extract(, .SD[REVIEWED == max(REVIEWED)], by = 'feature_id')
   swissprot <- fdata1[REVIEWED==1]
   trembl    <- fdata1[REVIEWED==0]
   fdata1 %>% report_n(prefix = 'Per group: drop trembl when swissprot available')

   # trembl groups
   #--------------
   if (nrow(trembl)>0){
      message('')
      trembl  %>% report_n(prefix = 'Trembl groups')

      # Drop fragments when full sequences available
      trembl  %>% magrittr::extract(, IS.FRAGMENT := 0)
      trembl  %>% magrittr::extract(, IS.FRAGMENT:= `PROTEIN-NAMES` %>% stringi::stri_detect_fixed('(Fragment)') %>% as.numeric())
      trembl %<>% magrittr::extract(, .SD[IS.FRAGMENT == min(IS.FRAGMENT)], by = 'feature_id')
      trembl[, IS.FRAGMENT:=NULL]
      trembl  %>% report_n(prefix = 'Per group: drop fragments when full available')

      # Use first sequence per gene
      trembl  %>% magrittr::extract(, N     := .N,                    by = 'feature_id')
      trembl  %>% magrittr::extract(, NGENE := length(unique(GENES)), by = 'feature_id')
      trembl %<>% magrittr::extract(, .SD[1],                         by = c('feature_id', 'GENES'))
      trembl  %>% report_n(prefix = 'Per group/gene: use first accession')
   }

   # swissprot groups
   #-----------------
   if (nrow(swissprot)>0){
      message('')
      swissprot  %>% report_n(prefix = 'Swissprot groups')

      # swissprot groups: collapse similar isoforms (shared accession)
      swissprot  %>% magrittr::extract(, ISOFORM := ISOFORM %>% paste0(collapse = ';'), by = c('feature_id', 'Uniprot accessions'))
      swissprot %<>% unique()
      swissprot  %>% report_n(prefix = 'Per group/accession: collapse spliceforms')

      # swissprot groups: collapse dissimilar isoforms: retain accession with maximum isoforms
      swissprot  %>% magrittr::extract(, NPERACCESSION := ISOFORM %>% stringi::stri_count_fixed(';'), by = c('feature_id', 'GENES'))
      swissprot  %>% magrittr::extract(, ISOFORM       := ISOFORM %>% paste0(collapse = ';'),         by = c('feature_id', 'GENES'))
      swissprot  %>% magrittr::extract(, `PROTEIN-NAMES` %>% autonomics.support::commonify_strings(), by = c('feature_id', 'GENES'))
      swissprot %<>% magrittr::extract(, .SD[NPERACCESSION==max(NPERACCESSION)],                      by = c('feature_id', 'GENES'))
      swissprot  %>% report_n(prefix = 'Per group/gene: use spliceform with most accessions')
      swissprot %<>% magrittr::extract(, .SD[1],                                                      by = c('feature_id', 'GENES'))
      swissprot  %>% report_n(prefix = 'Per group/gene: use first spliceform')
   }

   # paralogs
   #---------
   # Collapse paralogs: choose gene with most isoforms
   message('')
   swissprot[, NPERACCESSION:=NULL]
   trembl[, N:=NULL]
   trembl[, NGENE:=NULL]
   fdata1 <- rbind(swissprot, trembl)
   monologs <- fdata1[, .SD[length(unique(GENES))==1], by = c('feature_id')]
   paralogs <- fdata1[, .SD[length(unique(GENES))>1], by = c('feature_id')]
   if (nrow(paralogs)>0){
      paralogs  %>% report_n(prefix = 'Paralog groups')
      paralogs  %>% magrittr::extract(,  NISOFORMS := 1+ISOFORM %>% stringi::stri_count_fixed(';'))
      paralogs  %>% magrittr::extract(,  GENES          := GENES           %>% paste0(collapse = ';'),                  by = 'feature_id')
      paralogs  %>% magrittr::extract(,  ISOFORM        := ISOFORM         %>% paste0(collapse = ';'),                  by = 'feature_id')
      paralogs  %>% magrittr::extract(, `PROTEIN-NAMES` := `PROTEIN-NAMES` %>% autonomics.support::commonify_strings(), by = 'feature_id')
      paralogs %<>% magrittr::extract(, .SD[NISOFORMS == max(NISOFORMS)], by = 'feature_id')
      paralogs  %>% report_n(prefix = '   Per group: use paralog with most spliceforms')
      paralogs %<>% magrittr::extract(, .SD[1], by = 'feature_id')
      paralogs  %>% report_n(prefix = '   Per group: use first paralog')
      paralogs  %>% magrittr::extract(, NISOFORMS := NULL)
   }
   fdata1 <- rbind(monologs, paralogs)

   message('')
   fdata1 %>% report_n(prefix = 'All groups (deconvoluted)')

   # Merge back
   nullify_fvars <- function(object, fvars){
      for (curfvar in fvars)   autonomics.import::fdata(object)[[curfvar]] <- NULL
      return(object)
   }
   object %<>% nullify_fvars(fvars = c('Uniprot accessions', 'Protein names', 'Gene names'))
   autonomics.import::fdata(object) %<>% merge(fdata1, by = 'feature_id', sort = FALSE, all.x = TRUE)

   # Rename (MaxQuant style)
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('GENES',         'Gene names')
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('PROTEIN-NAMES', 'Protein names')
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('ISOFORM',       'Isoforms')

   # Remove unimportant fvars
   autonomics.import::fdata(object)$REVIEWED  <- NULL
   autonomics.import::fdata(object)$EXISTENCE <- NULL
   if (drop_isoform_info) autonomics.import::fdata(object)$Isoforms <- NULL

   # Return
   object

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

