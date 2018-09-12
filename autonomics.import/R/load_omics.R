


#===========================================
# GENERIC
#===========================================

#' Load omics
#' @param file        path to omics data file
#' @param platform   'exiqon', 'maxquant', 'metabolon', 'metabolonlipids', 'soma'
#' @param sheet       excel sheet number or name if applicable
#' @param quantity    string: which quantity to extract into exprs
#' @param design_file path to design file
#' @param log2_transform  logical
#' @param log2_offset     offset in mapping x -> log2(x+offset)
#' @param infer_design_from_sampleids logical
#' @param design_sep                  string: design separator
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'
#' # EXIQON
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                         package = 'subramanian.2016')
#'    file %>% autonomics.import::load_omics(
#'                platform = 'exiqon',
#'                infer_design_from_sampleids = TRUE)
#' }
#'
#' # PROTEINGROUPS
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemdiff/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    object <- file %>% autonomics.import::load_omics(
#'                          platform = 'maxquant',
#'                          quantity = 'Ratio normalized',
#'                          infer_design_from_sampleids = TRUE)
#' }
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_omics(sheet=2, platform = 'metabolon')
#' }
#'
#' # METABOLONLIPIDS
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_omics(sheet = 'Lipid Class Concentrations', platform = 'metabolonlipids')
#' }
#'
#' # RNASEQ
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    infer_design_from_sampleids <- TRUE
#' }
#' @importFrom magrittr %>%
#' @export
load_omics <- function(
   file,
   platform,
   sheet                       = 2,
   quantity                    = NULL,
   log2_transform              = TRUE,
   log2_offset                 = 0,
   design_file                 = NULL,
   infer_design_from_sampleids = FALSE,
   design_sep                  = NULL
){
   # Satisfy CHECK
   . <- NULL

   # Load components
   sdata1 <- file %>% autonomics.import::load_sdata(platform = platform, sheet = sheet, quantity = quantity)
   fdata1 <- file %>% autonomics.import::load_fdata(platform = platform, sheet = sheet)
   exprs1 <- file %>% autonomics.import::load_exprs(platform = platform, sheet = sheet, quantity = quantity)

   # Wrap into SummarizedExperiment
   object <- SummarizedExperiment::SummarizedExperiment(assays=list(exprs = exprs1))
   if (log2_transform)   autonomics.import::exprs(object) %<>% (function(x)log2(x+log2_offset))
   autonomics.import::sdata(object)  <- sdata1
   autonomics.import::fdata(object)  <- fdata1

   # Merge in design
   design_df <- autonomics.import::write_design(file, platform                    = platform,
                                                      infer_design_from_sampleids = infer_design_from_sampleids,
                                                      quantity                    = quantity,
                                                      design_sep                  = design_sep,
                                                      sheet                       = sheet) %>%
                autonomics.support::rm_empty_vars()

   object %<>% autonomics.import::merge_sdata(design_df, by = autonomics.import::sampleid_varname(platform))
   if (!is.null(design_file)){
      file_df <- autonomics.import::read_design(design_file)
      object %<>% autonomics.import::merge_sdata(file_df, by = autonomics.import::sampleid_varname(platform))
   }

   # Order on subgroup (and replicate)
   subgroup_values  <- object %>% autonomics.import::svalues('subgroup')
   replicate_values <- object %>% autonomics.import::svalues('replicate')
   if (!is.null(subgroup_values)){
      if (!is.null(replicate_values)){ object %<>% magrittr::extract(, order(subgroup_values, replicate_values))
      } else {                         object %<>% magrittr::extract(, order(subgroup_values))}
   }

   # Return
   object
}

#=============================
# EXIQON
#=============================

#' Load exiqon data
#' @param file                exiqon xlsx file
#' @param design_file         NULL or character (sample design file)
#' @param log2_transform      logical: whether to log2 transform
#' @param log2_offset         offset in mapping x -> log2(offset + x)
#' @param infer_design_from_sampleids  logical: whether to infer design from sample ids
#' @param design_sep          string: sample id separator
#' @param rm_ref_features     logical
#' @param rm_spike_features   logical
#' @param mean_center         logical: whether to mean_center exprs
#' @param flip_sign           logical: whether to flip sign
#' @param subtract_refgroup   logical
#' @param refgroup            character(1): reference group to subtract away
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'
#'    # Load exiqon
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                         package = 'subramanian.2016')
#'    file %>% load_exiqon(infer_design_from_sampleids = TRUE)
#'
#'    # Load exiqon with design file
#'    design_file <- tempfile()
#'    file %>% write_design('exiqon', infer_design_from_sampleids = TRUE,
#'                                    design_file = design_file)
#'    file %>% load_exiqon(design_file = design_file) %>%
#'             autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#' }
#' @importFrom magrittr %>%
#' @export
load_exiqon <- function(
   file,
   design_file                 = NULL,
   log2_transform              = NULL,
   log2_offset                 = NULL,
   infer_design_from_sampleids = FALSE,
   design_sep                  = NULL,
   rm_ref_features             = TRUE,
   rm_spike_features           = TRUE,
   mean_center                 = TRUE,
   flip_sign                   = TRUE,
   subtract_refgroup           = FALSE,
   refgroup                    = NULL
){
   # Satisfy CHECK
   . <- `#RefGenes` <- `#Spike` <- NULL

   # Load sumexp
   object <- autonomics.import::load_omics(file                        = file,
                                           platform                    = 'exiqon',
                                           log2_transform              = FALSE,
                                           log2_offset                 = log2_offset,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep)
   autonomics.import::prepro(object) <- list(assay = 'exiqon', entity = 'mirna', quantity = 'ct', software = 'genex')

   # Filter features
   if (rm_ref_features){
      object %<>% autonomics.import::filter_features(`#RefGenes`==0, verbose = TRUE)
      autonomics.import::fdata(object)$`#RefGenes` <- NULL
   }
   if (rm_spike_features){
      object %<>% autonomics.import::filter_features(`#Spike`   ==0, verbose = TRUE)
      autonomics.import::fdata(object)$`#Spike` <- NULL
   }

   # Mean center
   if (mean_center){
      autonomics.support::cmessage('\t\tMean center')
      autonomics.import::exprs(object) %<>% (function(x){
                                                sample_means <- x %>% (function(y){y[y>32] <- NA; y}) %>% colMeans(na.rm = TRUE)
                                                x %>% sweep(2, sample_means)
                                             })
   }

   # Flip sign
   if (flip_sign){
      autonomics.support::cmessage('\t\tFlip sign')
      object %<>% (function(x){autonomics.import::exprs(x) %<>% magrittr::multiply_by(-1); x})
   }

   # Subtract median reference exprs
   if (subtract_refgroup){
      if (is.null(refgroup))  refgroup <- autonomics.import::sdata(object)$subgroup[1]
      autonomics.support::cmessage("\t\tSubtract '%s' median", refgroup)
      object %<>% autonomics.preprocess::subtract_median_ref_exprs(ref = refgroup)
   }

   # Return
   object
}


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
#' @param object       SummerizedExperiment with proteinGroups data
#' @param fastafile    path to fastafile
#' @param fastafields  character vector: fields to load from fastafile
#' @return deconvoluted and annotated SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios[1:100, ]
#'    fastafile <- '../data/uniprot_hsa_20140515.fasta'
#'    if (file.exists(fastafile)){
#'       object %>% deconvolute_proteingroups(fastafile = fastafile)
#'    }
#' }
#' @importFrom magrittr %>% %<>%
#' @export
deconvolute_proteingroups <- function(
   object,
   fastafile,
   fastafields = c('GENES', 'PROTEIN-NAMES', 'EXISTENCE', 'REVIEWED')
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
   autonomics.support::cmessage('Load  %s  %s', basename(file), cur_sheet)

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



#==========================
# SOMA
#==========================

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
                        infer_design_from_sampleids = infer_design_from_sampleids,
                        design_sep                  = design_sep)

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

#=======================================================

#' Compute effective library sizes (using TMM)
#' @param counts counts matrix
#' @return vector (nsample)
#' @export
libsizes <- function(counts){
   colSums(counts) * edgeR::calcNormFactors(counts)
}

#' Convert counts into cpm (counts per million reads)
#' @param counts    count matrix
#' @param lib.size  scaled library sizes (vector)
#' @return cpm matrix
#' @export
counts_to_cpm <- function(counts, lib.size = libsizes(counts)){
   t(t(counts + 0.5)/(lib.size + 1) * 1e+06)
}

#' Invert cpm (counts per million reads) into counts
#'
#' @param cpm       cpm matrix
#' @param lib.size  scaled library sizes (vector)
#' @return count matrix
#' @export
cpm_to_counts <- function(cpm, lib.size){
   1e-06 * t(t(cpm) * (lib.size + 1))
}

#' @importFrom magrittr %>%
compute_precision_weights_once <- function(
   object,
   design = autonomics.import::create_design_matrix(object),
   plot = TRUE,
   ...
){

   # Extract
   log2cpm  <- autonomics.import::exprs(object)
   lib.size <- autonomics.import::sdata(object)$libsize

   # Assert
   n <- nrow(log2cpm)
   if (n < 2L) stop("Need at least two genes to fit a mean-variance trend")

   # Fit linear model
   fit <- limma::lmFit(log2cpm, design=design, ...)

   # Predict
   if (is.null(fit$Amean)) fit$Amean <- rowMeans(log2cpm, na.rm = TRUE)
   if (fit$rank < ncol(design)) {
      j <- fit$pivot[1:fit$rank]
      fitted.log2cpm <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,j, drop = FALSE])
   } else {
      fitted.log2cpm <- fit$coef %*% t(fit$design)
   }
   fitted.log2count <- (2^fitted.log2cpm) %>% cpm_to_counts(lib.size) %>% log2()

   # Fit mean-variance trend
   mean.log2count <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)  # mean log2 count
   sdrt.log2count <- sqrt(fit$sigma)                                     # sqrtsd(resid)
   all.identical <- matrixStats::rowVars(log2cpm)==0
   if (any(all.identical)) {
      mean.log2count <- mean.log2count[!all.identical]
      sdrt.log2count <- sdrt.log2count[!all.identical]
   }
   l <- lowess(mean.log2count, sdrt.log2count, f = 0.5)
   f <- approxfun(l, rule = 2)

   # Compute precision weights
   w <- 1/f(fitted.log2count)^4     # f(.) = sqrt(sd(.)) --> f(.)^4 = var(.)
   dim(w) <- dim(fitted.log2count)

   # Plot
   if (plot) {
      plot(mean.log2count, sdrt.log2count, xlab = "mean log2count", ylab = "sdrt log2count",
           pch = 16, cex = 0.25)
      title("voom: Mean-variance trend")
      lines(l, col = "red")
   }

   # Return
   return(w)
}


#' @rdname add_precision_weights
#' @importFrom magrittr %>%
#' @export
compute_precision_weights <- function(object, design = autonomics.import::create_design_matrix(object), plot = TRUE){

   # Estimate precision weights
   has_block <- autonomics.import::contains_block(object)
   weights <- object %>% compute_precision_weights_once(design = design, plot = !has_block)

   # Update precision weights using block correlation
   if (has_block){
      correlation <- log2cpm() %>%
         limma::duplicateCorrelation(design = design, block = object$block, weights = weights) %>%
         magrittr::extract2('consensus')
      weights <- object %>% compute_precision_weights_once(design = design,
                                                           block = object$block,
                                                           correlation = correlation,
                                                           plot = TRUE)
   }

   # Return
   dimnames(weights) <- dimnames(object)
   weights
}

#' Compute (and add) precision weights
#'
#' Compute (and add) precision weights for SummarizedExperiment with log2cpm values.
#'
#' Refactored version of limma::voom() which operates on a SummarizedExperiment with log2cpm values.
#' Created to allow a separation of the two steps in limma::voom():
#'    1. raw counts -> log2cpm (autonomics.import::counts_to_cpm)
#'    2. log2cpm    -> weights (autonomics.import::compute_precision_weights)
#' If a block factor is present, precision weights are updated using the block correlation information,
#' as suggested by Gordon Brown on BioC support (https://support.bioconductor.org/p/59700/)
#'
#' @param object    SummarizedExperiment
#' @param design    design matrix
#' @param plot      whether to plot mean-var trend (logical)
#' @param ...       passed to limma::lmFit
#' @return weights vector
#' @author Charity Law, Gordon Smyth, Aditya Bhagwat
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- 'extdata/rnaseq/gene_counts.txt'           %>%
#'               system.file(package = 'subramanian.2016') %>%
#'               autonomics.import::load_rnaseq(infer_design_from_sampleids = TRUE)
#'    object %>% autonomics.import::compute_precision_weights() %>% str()
#'    object %>% autonomics.import::add_precision_weights()
#'    object %>% autonomics.import::load_rnaseq()
#' }
#'
#' @importFrom magrittr  %>%
#' @export
add_precision_weights <- function(
   object,
   design = autonomics.import::create_design_matrix(object),
   plot = TRUE
){
   weights <- object %>% compute_precision_weights(design = design, plot = plot)
   SummarizedExperiment::assays(object)$weights <- weights
   object
}



#' Load RNAseq cpm
#' @param file                          rnaseq counts file
#' @param design_file                   sample design file
#' @param infer_design_from_sampleids   logical
#' @param design_sep                    character
#' @examples
#' require(magrittr)
#' if (require(autonomics.data){
#'    file <- system.file('extdata/stemdiff/rnaseq/gene_counts.txt', package = 'autonomics.data')
#' })
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    infer_design_from_sampleids <- TRUE
#'    object <- load_rnaseq(file, infer_design_from_sampleids = TRUE)
#'    object
#' }
#' @export
load_rnaseq <- function(
   file,
   design_file                 = NULL,
   infer_design_from_sampleids = FALSE,
   design_sep                  = NULL
){

   # Load sumexp
   autonomics.support::cmessage('\t\tload counts')
   object <- autonomics.import::load_omics(file                        = file,
                                           platform                    = 'rnaseq',
                                           log2_transform              = FALSE,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep)
   object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample()

   # Compute effective library size
   message('\t\tstore libsizes in sdata')
   autonomics.import::sdata(object)$libsize <- autonomics.import::exprs(object) %>% libsizes()

   # TMM-normalize counts
   message('\t\tcounts -> cpm')
   SummarizedExperiment::assays(object)$exprs %<>% counts_to_cpm()

   # Log2 transform
   message('\t\tcpm -> log2cpm')
   SummarizedExperiment::assays(object)$exprs %<>% log2()

   # Quantile normalize
   # Gordon:  prefer TMM over quantile normalization for most cases
   #          Use quantile normalization for extreme cases
   #          Don't use both
   # object %<>% limma::normalizeBetweenArrays(method = normalize.method)
   # Don't quantile normalize when using TMM
   # https://support.bioconductor.org/p/77664/

   # Add voom precision weights
   autonomics.support::cmessage('\t\tadd voom precision weights')
   object %<>% autonomics.import::add_precision_weights(plot = TRUE)

   # Add metadata
   autonomics.import::prepro(object) <- list(assay      = 'rnaseq',
                                             entity     = 'rna',
                                             quantity   = 'log2cpm',
                                             software   = 'limma',
                                             parameters = list())

   # Return
   object

}
