#==========================================================================================================
#' maxquant patterns
#' @export
maxquant_patterns <- c(`Ratio normalized`             =  '^Ratio ([HM]/[ML]) normalized (.+)$',
                       `Ratio`                        =  '^Ratio ([HM]/[ML]) (?!count|type|variability|iso-count|normalized)(.+)',
                       `LFQ intensity`                =  '^LFQ intensity ([HML])? ?(.+)$',
                       `Reporter intensity corrected` =  '^Reporter intensity corrected ([0-9]+) (.+)$',
                       `Reporter intensity`           =  '^Reporter intensity ([0-9]+) (.+)$',
                       `Intensity labeled`            =  '^Intensity ([HML]) (.+)$',
                       `Intensity`                    =  '^Intensity (.+)$')


#==========================================================================================================
#' Guess maxquant quantity from snames
#'
#' charactervector, dataframe, or SummarizedExperiment.
#'
#' @param x charactervector, dataframe, or SummarizedExperiment
#' @return  character(1): value from names(autonomics.import::maxquant_patterns)
#' @examples
#' require(magrittr)
#'
#' # charactervector
#'     autonomics.import::guess_maxquant_quantity("Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1")
#'     autonomics.import::guess_maxquant_quantity("Ratio M/L STD(L)_EM00(M)_EM01(H)_R1")
#'     autonomics.import::guess_maxquant_quantity("LFQ intensity EM00.R1")
#'     autonomics.import::guess_maxquant_quantity("Reporter intensity corrected 0 STD(0)EM00(1)EM01(2)_R1")
#'     autonomics.import::guess_maxquant_quantity("Reporter intensity 0 STD(0)EM00(1)EM01(2)_R1")
#'     autonomics.import::guess_maxquant_quantity("Intensity H STD(L)_EM00(M)_EM01(H)_R1")
#'
#' # dataframe
#'     if (require(autonomics.data)){
#'       x <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data') %>%
#'             data.table::fread()
#'       autonomics.import::guess_maxquant_quantity(x)
#'     }
#'
#' # SummarizedExperiment
#'     if (require(autonomics.data)){
#'       x <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data') %>%
#'             autonomics.import::read_proteingroups_asis()
#'       autonomics.import::guess_maxquant_quantity(x)
#'     }
#' @export
guess_maxquant_quantity <- function(x, ...){
   UseMethod("guess_maxquant_quantity", x)
}

#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.character <- function(x){      # x = character vector of maxquant colnames
   for (quantity in names(autonomics.import::maxquant_patterns)){
      pattern <- autonomics.import::maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.data.frame <- function(x){     # x = maxquant dataframe
   x <- names(x)
   for (quantity in names(autonomics.import::maxquant_patterns)){
      pattern <- autonomics.import::maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}

#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.SummarizedExperiment <- function(x){     # x = SummarizedExperiment
   x <- autonomics.import::snames(x)
   for (quantity in names(autonomics.import::maxquant_patterns)){
      pattern <- autonomics.import::maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}


#==========================================================================================================
#' Clean maxquant snames
#'
#' For charactervector or SummarizedExperiment
#'
#' Drop "Ratio normalized", "LFQ intensity" etc from maxquant snames & sample_id values
#'
#' @param x        character(.) or SummarizedExperiment
#' @param verbose  logical(1)
#' @examples
#' require(magrittr)
#'
#' # character vector
#'     autonomics.import::clean_maxquant_snames("Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1")
#'     autonomics.import::clean_maxquant_snames("Ratio M/L STD(L)_EM00(M)_EM01(H)_R1")
#'     autonomics.import::clean_maxquant_snames('LFQ intensity STD_R1')
#'     autonomics.import::clean_maxquant_snames('LFQ intensity L STD(L)_EM00(M)_EM01(H)_R1')
#'     autonomics.import::clean_maxquant_snames('Reporter intensity 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1')
#'     autonomics.import::clean_maxquant_snames('Reporter intensity corrected 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1')
#'
#' # SummarizedExperiment
#' if (require(autonomics.data)){
#'      x <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'            system.file(package = 'autonomics.data')     %>%
#'            autonomics.import::read_proteingroups_asis()
#'      x %>% autonomics.import::clean_maxquant_snames(verbose = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
clean_maxquant_snames <- function (x, ...) {
   UseMethod("clean_maxquant_snames", x)
}


#' @importFrom magrittr %>%
#' @export
#' @rdname clean_maxquant_snames
clean_maxquant_snames.character <- function(
   x,
   quantity = autonomics.import::guess_maxquant_quantity(x),
   verbose  = FALSE
){
   # x = mix + channel. Return mix if single channel.
   pattern <- autonomics.import::maxquant_patterns %>% magrittr::extract2(quantity)
   channel <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
   mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$2')
   if (all(channel=='')){ cleanx <- mix
   } else               { cleanx <- sprintf('%s{%s}', mix, channel)
   }
   autonomics.support::cmessage('\t\tClean snames: %s  ->  %s', x[1], cleanx[1])
   return(cleanx)
}

#' @importFrom magrittr %>%
#' @export
#' @rdname clean_maxquant_snames
clean_maxquant_snames.SummarizedExperiment <- function(
   x,
   quantity = autonomics.import::guess_maxquant_quantity(x),
   verbose  = FALSE
){
   newsnames <- autonomics.import::snames(x) %>% clean_maxquant_snames(quantity = quantity, verbose=verbose)
   autonomics.import::snames(x) <- autonomics.import::sdata(x)$sample_id <- newsnames
   x
}


#=========================================================================================================
#' Uniplex snames
#'
#' For charactervector or SummarizedExperiment
#'
#' @param x        character vector or SummarizedExperiment
#' @param verbose  logical(1)
#' @examples
#' require(magrittr)
#'
#' # character vector
#'     # Alternate multiplexing forms supported
#'      autonomics.import::uniplex_snames("STD(L)_EM00(M)_EM01(H)_R1{M/L}")     # Label Ratio
#'      autonomics.import::uniplex_snames('A(0)_B(1)_C(2)_D(3)_R1{0}'     )     # Reporter intensity
#'      autonomics.import::uniplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')       # Label Intensity
#'
#'    # Alternate separators supported
#'      autonomics.import::uniplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')       # underscore
#'      autonomics.import::uniplex_snames('STD(L).EM00(M).EM01(H).R1{L}')       # dot
#'      autonomics.import::uniplex_snames('STD(L)EM00(M)EM01(H).R1{L}')         # no separator
#'
#'    # Composite snames supported
#'      autonomics.import::uniplex_snames("WT.t0(L)_WT.t1(M)_WT.t2(H)_R1{H/L}") # composite snames
#'
#'    # Uniqueness ensured by appending labels when necessary
#'      autonomics.import::uniplex_snames(c("STD(L).BM00(M).BM00(H).R10{M/L}",  # implicit uniquification
#'                                          "STD(L).BM00(M).BM00(H).R10{H/L}"))
#'    # Uniplexed snames are returned unchanged
#'      autonomics.import::uniplex_snames(c('STD_R1', 'EM0_R1'))
#'
#' # SummarizedExperiment
#'   if (require(autonomics.data)){
#'      x <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'            system.file(package = 'autonomics.data')     %>%
#'            autonomics.import::read_proteingroups_asis() %>%
#'            autonomics.import::clean_maxquant_snames()
#'      x %>% autonomics.import::uniplex_snames(verbose = TRUE)
#'   }
#'
#' @export
uniplex_snames <- function (x, ...) {
   UseMethod("uniplex_snames", x)
}

#' @rdname uniplex_snames
#' @export
uniplex_snames.character <- function(x, verbose = FALSE){

   # Return unchanged if not multiplexed
   # KD(H)WT(L){H/L}
   pattern <- '(.+)\\{(.+)\\}'
   n_open   <- x %>% stringi::stri_count_fixed('(')
   n_closed <- x %>% stringi::stri_count_fixed(')')
   is_multiplexed <- all(stringi::stri_detect_regex(x, pattern) & (n_open==n_closed) & (n_open>0))
   if (!is_multiplexed) return(x)

   # Separate mix and channel
   mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
   channel <- x %>% stringi::stri_replace_first_regex(pattern, '$2')

   # Separate labels and samples
   pattern <- '\\(.+?\\)'
   labels  <- mix %>% stringi::stri_extract_all_regex(pattern) %>% lapply(stringi::stri_replace_first_fixed, '(', '') %>%
              lapply(stringi::stri_replace_first_fixed, ')', '')
   samples <- mix %>% stringi::stri_split_regex(pattern) %>%
                      # rm sep from samples (but not from replicate - needed to glue back later!)
                      lapply(function(y){y[1:length(labels[[1]])] %<>% stringi::stri_replace_first_regex('^[_. ]', ''); y})

   # Return unchanged if mixes differ in no of labels or samples
   are_all_identical <- function(y) if (length(y)==1) TRUE else all(y[-1] == y[1])
   n_samples <- vapply(samples,  length, integer(1))
   n_labels  <- vapply(labels, length, integer(1))
   if (!are_all_identical(n_samples) | !are_all_identical(n_labels)){
      autonomics.support::cmessage('\t\tCannot demultiplexing snames: mixes differ in number of samples or labels')
      return(x)
   }

   # Extract replicate
   n_samples %<>% unique()
   n_labels  %<>% unique()
   if (n_samples > n_labels){ replicate <- mix %>% stringi::stri_split_regex(pattern) %>% vapply((function(y) y %>% magrittr::extract(length(y))), character(1))
   samples %<>% lapply(extract, 1:(n_samples-1))
   } else {                   replicate <- rep('', length(samples))
   }

   # Extract channel samples from mix
   is_ratio <- channel %>% stringi::stri_detect_fixed('/') %>% all()
   samples %<>% mapply(set_names, ., labels, SIMPLIFY = FALSE)
   if (is_ratio){
      num_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(extract, character(1), 1)
      den_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(extract, character(1), 2)
      den_samples <- mapply(extract, samples, den_label)
      num_samples <- mapply(extract, samples, num_label)
      xdemultiplex <- sprintf('%s_%s%s', num_samples, den_samples, replicate)
   } else {
      samples %<>% mapply(extract, ., channel)
      xdemultiplex <- sprintf('%s%s', samples, replicate)
   }
   if (verbose) autonomics.support::cmessage('\t\tDemultiplex snames: %s  ->  %s', x[1], xdemultiplex[1])

   # Ensure uniqueness. Add labels if required.
   idx <- autonomics.support::cduplicated(xdemultiplex) %>% which()
   if (length(idx)>0){
      label_tags <- channel[idx] %>% stringi::stri_replace_first_fixed('/', '')
      if (verbose)   autonomics.support::cmessage('\t\tUniquify snames: %s -> %s%s (for %d/%d snames)',
                                                  xdemultiplex[idx][1], xdemultiplex[idx][1], label_tags[1],
                                                  length(idx), length(xdemultiplex))
      xdemultiplex[idx] %<>% paste0(label_tags)
   }

   # Return
   return(xdemultiplex)
}

#' @rdname uniplex_snames
#' @importFrom magrittr %>%
#' @export
uniplex_snames.SummarizedExperiment <- function(
   x,
   verbose  = FALSE
){
   newsnames <- autonomics.import::snames(x) %>% autonomics.import::uniplex_snames(verbose = verbose)
   autonomics.import::snames(x) <- autonomics.import::sdata(x)$sample_id <- newsnames
   x
}

#' proteingroups fvars
#' @export
proteingroups_fvars <- c(c('id', 'Majority protein IDs', 'Protein names', 'Gene names', 'Contaminant', 'Potential contaminant', 'Reverse', 'Phospho (STY) site IDs'))

#===========================================================================================================================================
#' @rdname read_proteingroups
#' @importFrom magrittr %>% %<>%
#' @export
read_proteingroups_asis <- function(
   file,
   quantity       = NULL,
   fvars          = autonomics.import::proteingroups_fvars,
   verbose        = TRUE
){

   # Initial Read
   assertive.files::assert_all_are_existing_files(file)
   dt <- data.table::fread(file, integer64 = 'numeric', header = TRUE)
   if (is.null(quantity)) quantity <- dt %>% autonomics.import::guess_maxquant_quantity()
   fvars %<>% intersect(names(dt))

   # Define components
   fid_rows   <- 2:nrow(dt)
   fid_cols   <- which(names(dt) == 'id')
   sid_rows   <- 1
   sid_cols   <- names(dt) %>% stringi::stri_detect_regex(autonomics.import::maxquant_patterns[[quantity]]) %>% which()
   expr_rows  <- 2:nrow(dt)
   expr_cols  <- sid_cols
   fvar_rows  <- 1
   fvar_cols  <- match(fvars, names(dt))
   fdata_rows <- 2:nrow(dt)
   fdata_cols <- fvar_cols

   # Read sumexp
   object <- file %>% autonomics.import::read_omics_asis(fid_rows   = fid_rows,     fid_cols   = fid_cols,
                                                         sid_rows   = sid_rows,     sid_cols   = sid_cols,
                                                         expr_rows  = expr_rows,    expr_cols  = expr_cols,
                                                         fvar_rows  = fvar_rows,    fvar_cols  = fvar_cols,
                                                         fdata_rows = fdata_rows,   fdata_cols = fdata_cols,
                                                         transpose  = FALSE,
                                                         verbose    = verbose)

   contaminant_var <- c('Contaminant', 'Potential contaminant') %>% intersect(autonomics.import::fvars(object))
   autonomics.import::fdata(object)[[contaminant_var]] %<>% (function(x){x[is.na(x)] <- ''; x})
   autonomics.import::fdata(object)[['Reverse'      ]] %<>% (function(x){x[is.na(x)] <- ''; x})

   # Return
   object

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
   autonomics.support::cmessage('\t\tLoad fasta file')
   fasta_annotations <- fastafile %>% autonomics.annotate::load_uniprot_fasta_annotations(fastafields)

   # Uncollapse
   fdata1 <- object %>%
             autonomics.import::fdata() %>%
             magrittr::extract(, c('feature_id', 'Uniprot accessions'), drop = FALSE) %>%
             tidyr::separate_rows(`Uniprot accessions`, sep = ';') %>%
             data.table::data.table()

   # Split into CANONICAL and isoform
   fdata1 %<>% magrittr::extract(, ISOFORM     := `Uniprot accessions`) %>%
      magrittr::extract(, `Uniprot accessions` := `Uniprot accessions` %>% stringi::stri_replace_first_regex('[-][0-9]+',     '')) # %>%
   #magrittr::extract(, ISOFORM := ISOFORM %>% sort() %>% unique() %>% paste0(collapse=';'), by = c('feature_id', 'Uniprot accessions')) %>%
   #unique()

   # Merge in uniprot fasta annotations
   nunmapped <- fdata1 %>%
                magrittr::extract(ISOFORM %in% setdiff(ISOFORM, fasta_annotations$UNIPROTKB)) %>%
                magrittr::extract(, .SD[1], by = 'feature_id') %>%
                nrow()
   autonomics.support::cmessage('\t\tDeconvolute %d/%d proteingroups with sequences from fastafile',
                                nrow(object) - nunmapped, nrow(object))
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
      autonomics.support::cmessage('\t\t\t%s%d proteingroups -> %d singlegene -> %d singleprotein -> %d singleseq%s',
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
      paralogs  %>% report_n(prefix = 'Per group: use paralog with most spliceforms')
      paralogs %<>% magrittr::extract(, .SD[1], by = 'feature_id')
      paralogs  %>% report_n(prefix = 'Per group: use first paralog')
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

#=========================================================================================================
#' Read proteingroups file
#' @param file                    character(1)
#' @param quantity                character(1): any value in names(maxquant_patterns)
#' @param fvars                   character(n): names of fvar columns
#' @param rm_reverse              logical(1)
#' @param rm_contaminants         logical(1)
#' @param rm_complete_nondetects  logical(1)
#' @param clean_snames            logical(1)
#' @param uniplex_snames          logical(1)
#' @param log2transform           logical(1)
#' @param fastafile               NULL or character(1): if provided, heuristic proteingroups deconvolution is attempted
#' @param verbose                 logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::read_proteingroups_asis()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% read_proteingroups()
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_proteingroups <- function(
   file,
   quantity               = NULL,
   fvars                  = autonomics.import::proteingroups_fvars,
   rm_reverse             = TRUE,
   rm_contaminants        = TRUE,
   rm_complete_nondetects = TRUE,
   clean_snames           = TRUE,
   uniplex_snames         = TRUE,
   log2transform          = TRUE,
   fastafile              = NULL,
   verbose                = TRUE
){

   # Read
   if (verbose) autonomics.support::cmessage('\tRead proteinGroups.txt')
   object <- file %>% autonomics.import::read_proteingroups_asis(quantity       = quantity,
                                                                 fvars          = fvars,
                                                                 verbose        = verbose)
   # Process exprs
   if (verbose) autonomics.support::cmessage('\texprs')
   if (rm_reverse)              object %<>% autonomics.import::filter_features_("Reverse != '+'", verbose = verbose)
   if (rm_contaminants){        contaminant_var <- c('Potential contaminant', 'Contaminant') %>% intersect(autonomics.import::fvars(object))
                                object %<>% autonomics.import::filter_features_(sprintf("`%s` != '+'", contaminant_var), verbose = verbose)
   }
   if (rm_complete_nondetects)  object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample(verbose = verbose)
   #object %>% autonomics.preprocess::qrilc_consistent_nondetects()
   if (log2transform){
      if (verbose) autonomics.support::cmessage('\t\tLog2 transform')
      autonomics.import::exprs(object) %<>% log2()
   }

   # Process sdata
   if (verbose) autonomics.support::cmessage('\tsdata')
   if (clean_snames)    object %<>% autonomics.import::clean_maxquant_snames(verbose = verbose)
   if (uniplex_snames)  object %<>% autonomics.import::uniplex_snames(verbose = verbose)
   object$subgroup <- object$sample_id %>% autonomics.import::guess_subgroup_values(verbose = verbose)
   #object$block    <- object$sample_id %>% autonomics.import::guess_subject_values( verbose = TRUE)

   # Process fdata
   if (verbose) autonomics.support::cmessage('\tfdata')
   if (verbose) autonomics.support::cmessage("\t\tRename: 'Majority protein IDs' -> 'Uniprot accessions'")
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('Majority protein IDs', 'Uniprot accessions')
   if (!is.null(fastafile)) object %>% autonomics.import::deconvolute_proteingroups(fastafile)

   # Return
   object
}


#======================================================================================================================================
#' Read phosphosites file
#' @param file                   character(1)
#' @param proteingroups_file     character(1)
#' @param quantity               NULL or value in names(autonomics.import::maxquant_patterns)
#' @param fvars                  character(n)
#' @param min_localization_prob  numeric(1)
#' @param add_occupancies        logical(1): whether to add occupancies(object) = phosphosite - proteingroups
#' @param clean_snames           logical(1)
#' @param uniplex_snames         logical(1)
#' @param verbose                logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemdiff/maxquant/phospho (STY)Sites.txt', package = 'autonomics.data')
#'    file %>% autonomics.import::read_phosphosites_asis()
#' }
#' @export
read_phosphosites_asis <- function(
   file,
   proteingroups_file     = file %>% stringi::stri_replace_first_fixed('phospho (STY)Sites.txt', 'proteinGroups.txt'),
   quantity               = NULL,
   fvars                  = c('id', 'Protein group IDs', 'Positions within proteins', 'Localization prob'),
   min_localization_prob  = 0.75,
   clean_snames           = TRUE,
   uniplex_snames     = TRUE
){
   # Initial Read
   assertive.files::assert_all_are_existing_files(c(file, proteingroups_file))
   dt <- data.table::fread(file, integer64 = 'numeric', header = TRUE)
   if (is.null(quantity)) quantity <- dt %>% autonomics.import::guess_maxquant_quantity()
   pattern <- autonomics.import::maxquant_patterns %>% magrittr::extract2(quantity)
   value_cols <- names(dt) %>% (function(x) stringi::stri_detect_regex(x,pattern) & !stringi::stri_detect_regex(x, '___[1-3]')) %>% which()
   fvar_cols  <- which(names(dt) %in% fvars)
   #dt %<>% magrittr::extract(, c(fvars, value_cols), with = FALSE)

   # Read phosphosites
   fid_rows   <- 2:nrow(dt)
   fid_cols   <- which(names(dt) == 'id')
   sid_rows   <- 1
   sid_cols   <- value_cols
   expr_rows  <- 2:nrow(dt)
   expr_cols  <- value_cols
   fvar_rows  <- 1
   fvar_cols  <- fvar_cols
   fdata_rows <- 2:nrow(dt)
   fdata_cols <- fvar_cols
   phosphosites  <- file  %>% autonomics.import::read_omics_asis(fid_rows   = fid_rows,     fid_cols   = fid_cols,
                                                                 sid_rows   = sid_rows,     sid_cols   = sid_cols,
                                                                 expr_rows  = expr_rows,    expr_cols  = expr_cols,
                                                                 fvar_rows  = fvar_rows,    fvar_cols  = fvar_cols,
                                                                 fdata_rows = fdata_rows,   fdata_cols = fdata_cols,
                                                                 transpose  = FALSE,
                                                                 verbose    = TRUE)
   # Filter phosphosites
   phosphosites %<>% autonomics.import::filter_features(!stringi::stri_detect_fixed(`Protein group IDs`, ';'), verbose = TRUE)
   phosphosites %<>% autonomics.import::filter_features(`Localization prob` >= min_localization_prob,          verbose = TRUE)

   # Calculate occupancies (allows to disentangle phosphorylation and protein expression)
   autonomics.support::cmessage('\t\tAdd occupancies(object) = phosphosites - proteingroups')
   proteingroups <- proteingroups_file %>%
      autonomics.import::read_proteingroups_asis(quantity = quantity, clean_snames = FALSE, uniplex_snames = FALSE, verbose = FALSE) %>%
      magrittr::extract(phosphosites %>% autonomics.import::fvalues("Protein group IDs"), )
   assertive.base::assert_all_are_true(autonomics.import::snames(proteingroups) == autonomics.import::snames(phosphosites))
   occupancies <- autonomics.import::exprs(phosphosites) %>% magrittr::subtract(autonomics.import::exprs(proteingroups))
   autonomics.import::occupancies(phosphosites) <- occupancies

   # Clean and demultiplex snames
   if (clean_snames)       autonomics.import::snames(phosphosites) %<>% autonomics.import::clean_maxquant_snames(verbose = TRUE)
   if (uniplex_snames) autonomics.import::snames(phosphosites) %<>% autonomics.import::uniplex_snames(verbose = TRUE)

   # Return
   phosphosites

}


