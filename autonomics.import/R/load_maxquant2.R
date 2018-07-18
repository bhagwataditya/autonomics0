#' Extract injection values from maxquant file
#' @param file 'proteinGroups.txt'
#' @return character vector: injection values
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_injections()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt',package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_injections()
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_injections()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_injections <- function(file){
   file %>%
   autonomics.support::cfread() %>%
   names() %>%
   magrittr::extract(stringi::stri_detect_fixed(., 'Razor + unique peptides ')) %>%
   stringi::stri_replace_first_fixed('Razor + unique peptides ', '')
}

#' Extract channel values from maxquant file
#' @param file string: path to proteinGroups.txt file (or other maxquant file)
#' @return character vector: maxquant channel values
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_channels()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt',package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_channels()
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_channels()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_channels <- function(file){
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   file %>% autonomics.support::cfread() %>%
            names()                      %>%
            # channel specific intensities available, except for reporter intensities
            magrittr::extract(autonomics.support::vstri_detect_fixed(., c('Intensity ', 'Reporter intensity corrected'))) %>%
            autonomics.support::vstri_replace_first_fixed(c('Intensity ', 'Reporter intensity corrected'), '')            %>%
            autonomics.support::vstri_replace_first_fixed(injections, '')                                                 %>%
            trimws()                     %>%
            unique()                     %>%
            setdiff("")
}


#' Extract maxquant fnames
#' @param file full path to proteinGroups.txt
#' @return character vector with sample names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_fnames() %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_fnames <- function(file){
   file %>%
   autonomics.support::cfread() %>%
   magrittr::extract2('Majority protein IDs') %>%
   stringi::stri_split_fixed(';') %>%
   vapply(extract, character(1), 1)
}


#' Extract maxquant intensity colnames
#' @param file full path to proteinGroups.txt
#' @param quantity 'Intensity', 'LFQ intensity', 'Reporter intensity'
#' @return character vector with intensity column names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames() %>% head(3)
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames() %>% head(3)
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames('LFQ intensity') %>% head(3)
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames('Reporter intensity') %>% head(3)
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_intensity_colnames <- function(file, quantity = 'Intensity'){

   # Asser
   assertive.files::assert_all_are_existing_files(file)
   assertive.sets::assert_is_subset(quantity, c('Intensity', 'LFQ intensity', 'Reporter intensity'))

   # Deduce injections and channels from unambiguous peptide columns
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()

   # Construct intensity colnames
   intensity_colnames <- if (length(channels)==0){                      sprintf('%s %s',    quantity,           injections)
                         } else {                  autonomics.support::vsprintf('%s %s %s', quantity, channels, injections) }

   # Ensure identical order as in actual file
   names(autonomics.support::cfread(file)) %>% magrittr::extract(. %in% intensity_colnames)
}

#' Extract maxquant ratio colnames
#' @param file full path to proteinGroups.txt
#' @param quantity 'Ratio' or 'Ratio normalized'
#' @return character vector with ratio column names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% extract_maxquant_ratio_colnames('Ratio') %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_ratio_colnames <- function(file, quantity){

   assertive.sets::assert_is_subset(quantity, c('Ratio', 'Ratio normalized'))

   # Deduce injections and channels from unambiguous peptide columns
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()

   # Construct possible ratio colnames
   possible_ratio_columns <- autonomics.support::vsprintf('Ratio %s/%s%s %s',
                                                          channels,
                                                          channels,
                                                          if(quantity == 'Ratio normalized') ' normalized' else '',
                                                          injections)
   # Return actual colnames in correct order
   file %>%
   autonomics.support::cfread() %>%
   names() %>%
   magrittr::extract(. %in% possible_ratio_columns)
}


#' Extract maxquant exprs
#' @param file full path to proteinGroups.txt
#' @param quantity 'Ratio', 'Ratio normalized', Intensity', 'LFQ intensity', 'Reporter intensity'
#' @return character vector with sample names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_exprs(quantity = 'Intensity') %>%
#'             magrittr::extract(1:3, 1:3)
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_exprs('Intensity')     %>%
#'             magrittr::extract(1:3, 1:3)
#'    file %>% autonomics.import::extract_maxquant_exprs('LFQ intensity') %>%
#'             magrittr::extract(1:3, 1:3)
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_exprs('Reporter intensity') %>%
#'             magrittr::extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_exprs <- function(file, quantity = 'Intensity'){

   # Assert
   assertive.sets::assert_is_subset(quantity, c('Ratio', 'Ratio normalized', 'Intensity', 'LFQ intensity', 'Reporter intensity'))

   # Construct maxquant colnames
   col_names <- if (quantity %in% c('Ratio', 'Ratio normalized')){ file %>% autonomics.import::extract_maxquant_ratio_colnames(quantity)
                } else {                                           file %>% autonomics.import::extract_maxquant_intensity_colnames(quantity)
                }

   # Extract exprs matrix
   exprs_mat <- file %>%
                autonomics.support::cfread(select = col_names) %>%
                data.matrix() %>%
                magrittr::set_rownames(autonomics.import::extract_maxquant_fnames(file))

   # Rename samples
   if (quantity %in% c('Ratio', 'Ratio normalized')){
     colnames(exprs_mat) %<>% stringi::stri_replace_first_regex('Ratio (./.) (?:normalized )?(.+)',   '$2[$1]')

   } else {
     # It is better to do it this way (rather than use a stri_replace_regex as for the ratios)
     # because ' ' separators in sample names are difficult to differentiate from ' H' constructs
     injections <- file %>% autonomics.import::extract_maxquant_injections()
     channels   <- file %>% autonomics.import::extract_maxquant_channels()
     colnames(exprs_mat) <- if (length(channels)==0) injections else autonomics.support::vsprintf('%s[%s]', injections, channels)
   }

   # Return
   exprs_mat
}

#' Extract maxquant fdata
#' @param file path to proteinGroups.txt file
#' @return dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_fdata() %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_fdata <- function(file){
   `Majority protein IDs` <- NULL

   file %>%
   autonomics.support::cfread(select = c('Majority protein IDs', 'Gene names', 'Protein names')) %>%
   magrittr::extract(, feature_id := `Majority protein IDs` %>% stringi::stri_split_fixed(';') %>% vapply(extract, character(1), 1)) %>%
   data.frame(stringsAsFactors = FALSE, check.names = FALSE, row.names = .$feature_id) %>%
   autonomics.support::pull_columns('feature_id')
}

#' Infer maxquant quantity
#' @param file path to maxquant file
#' @return string
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::infer_maxquant_quantity()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::infer_maxquant_quantity()
#' }
#' @importFrom magrittr %>%
#' @export
infer_maxquant_quantity <- function(file){
   x <- autonomics.support::cfread(file) %>% names()
   if (any(stringi::stri_detect_fixed(x, 'Ratio') & stringi::stri_detect_fixed(x, 'normalized')))   return('Ratio normalized')
   if (any(stringi::stri_detect_fixed(x, 'Reporter intensity')))                                    return('Reporter intensity')
   if (any(stringi::stri_detect_fixed(x, 'LFQ intensity')))                                         return('LFQ intensity')
                                                                                                    return('Intensity')
}


#' Load proteingroups
#' @param file            path to proteinGroups.txt
#' @param quantity       'Ratio normalized', 'Ratio', 'Intensity', 'LFQ intensity', 'Reporter intensity'
#' @param infer_design_from_sampleids  logical: whether to infer design from sampleids
#' @param design_file     path to design file (created with write_maxquant_design)
#' @param fasta_file      path to uniprot fasta database
#' @param log2_transform  logical: whether to log2 transform
#' @param log2_offset     numeric: offset used in mapping: x -> log2(offset + x)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    object <- file %>% autonomics.import::load_proteingroups2()
#' }
#' @importFrom magrittr %>%
#' @export
load_proteingroups2 <- function(
   file,
   quantity                    = autonomics.import::infer_maxquant_quantity(file),
   infer_design_from_sampleids = FALSE,
   design_file                 = NULL,
   fasta_file                  = NULL,
   log2_transform              = TRUE,
   log2_offset                 = 0
){

   # Load exprs
   exprs_mat <- file %>% autonomics.import::extract_maxquant_exprs(quantity)
   if (log2_transform) exprs_mat %<>% (function(x) log2(log2_offset + x))

   # Pack into Sumexp
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs_mat))
   autonomics.import::fdata(object) <- file %>% autonomics.import::extract_maxquant_fdata()
   autonomics.import::sdata(object) <- data.frame(sample_id = colnames(exprs_mat), row.names = colnames(exprs_mat))

   # Merge in sample design
   design_df <- autonomics.import::write_maxquant_design(file, infer_from_sampleids = infer_design_from_sampleids)
   object %<>% autonomics.import::merge_sdata(design_df, by = 'sample_id')
   if (!is.null(design_file)){
      file_df <- autonomics.import::read_maxquant_design(design_file)
      object %<>% autonomics.import::merge_sdata(file_df, by = 'sample_id')
   }

   # Return
   return(object)

}


#' Annotate and deconvolute proteingroups SumExp
#'
#' Steps:
#' \enumerate{
#'    \item Separates proteingroups into uniprot accessions.
#'    \item Fetches up-to-date annotations for each from uniprot
#'    \item Keeps best annotated entries (per proteingroup)
#'    \item Collapses isoforms (per canonical accession)
#'    \item Keeps first of redundant uniprot entries
#'    \item Collapses paralogs (per protein group)
#'    \item Merges annotations into fdata(object) and returns
#' }
#'
#' @param object SummarizedExperiment
#' @param fasta_file path to fasta file
#' @examples
#' \dontrun{
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- system.file('extdata/billing2016/proteinGroups.txt',
#'                           package = 'autonomics.data') %>%
#'              autonomics.import::load_proteingroups() %>%
#'              magrittr::extract(1:10, )
#'    object %>% autonomics.import::fdata() %>% str()
#'    object %>% autonomics.import::annotate_proteingroups() %>%
#'               autonomics.import::fdata() %>% str()
#' }
#' }
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%
#' @export
annotate_proteingroups2 <- function(object, fasta_file = NULL){
   Uniprot <- NULL

   # rm existing annotation to avoid confusion
   autonomics.import::fdata(object)$`Gene names`    <- NULL
   autonomics.import::fdata(object)$`Protein names` <- NULL
   autonomics.import::fdata(object)$`Fasta headers` <- NULL

   # Extract
   fdata1 <- autonomics.import::fdata(object)               %>%
             magrittr::set_colnames(colnames(.) %>% stringi::stri_replace_first_fixed('Majority protein IDs', 'Uniprot')) %>%
             tidyr::separate_rows('Uniprot', sep = ';') %>%
             data.table::data.table() %>%
             magrittr::extract(, ('Canonical') := Uniprot %>% stringi::stri_replace_first_regex('[-][0-9]+$', ''))
   n0 <- nrow(object)
   n1 <- length(unique(fdata1$Uniprot))
   n2 <- length(unique(fdata1$Canonical))
   autonomics.support::cmessage('%d protein groups -> %d uniprot accessions -> %d canonical', n0, n1, n2)

   # Annotate through uniprot.ws
   if (is.null(fasta_file)){
      autonomics.support::cmessage('Annotate with uniprot.ws')
      dt <- unique(fdata1$`Canonical accessions`) %>% autonomics.annotate::annotate_uniprot()
      fdata1 %<>% merge(dt, by.x = 'Canonical accessions', by.y = uniprot_var, all.x = TRUE)

   # Or annotate through fasta file
   } else {
      fasta_list <- seqinr::read.fasta(fasta_file)

      names <- fasta_list %>% (function(x)attr(x, 'name'))
      sp_or_tr           <- names %>% stringi::stri_split_fixed('|') %>% vapply(extract, character(1), 1)
      uniprot_accessions <- names %>% stringi::stri_split_fixed('|') %>% vapply(extract, character(1), 2)
      entry_names        <- names %>% stringi::stri_split_fixed('|') %>% vapply(extract, character(1), 3)

      annotations <- fasta_list %>% vapply(function(x) attr(x, 'Annot'), character(1))
      annotations[1:3] %>% stringi::stri_replace_first_fixed(paste0('>', names[1:3], ' '), '')

   }

   # Keep reviewed, drop unreviewed (per proteingroup)
   autonomics.support::cmessage('Simplify proteingroups')
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, .SD[reviewed == max(reviewed)], by = fid_var) %>% magrittr::extract(, reviewed:=NULL)
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: prefer reviewed (when available, per proteingroup)', n0, n1)

   # Keep best (annotation) score (per protein group)
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, .SD[score == max(score)],       by = fid_var) %>% magrittr::extract(, score:=NULL)
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: prefer best annotation score (per proteingroup)', n0, n1)

   # Keep best go annotation (per protein group)
   n0 <- nrow(fdata1)
   fdata1  %>% magrittr::extract( is.na(goid), ngoid := 0)
   fdata1  %>% magrittr::extract(!is.na(goid), ngoid := goid     %>% stringi::stri_count_fixed(';') %>% magrittr::add(1))
   fdata1 %<>% magrittr::extract(, .SD[ngoid  == max(ngoid)   ], by = fid_var)
   fdata1 %>%  magrittr::extract(, ngoid := NULL)
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: prefer best GO annotation (per protein group)', n0, n1)

   # Keep best interpro annotated (per protein group)
   n0 <- nrow(fdata1)
   fdata1  %>% magrittr::extract( is.na(interpro), ninterpro := 0)
   fdata1  %>% magrittr::extract(!is.na(interpro), ninterpro := interpro %>% stringi::stri_count_fixed(';') %>% magrittr::add(1))
   fdata1 %<>% magrittr::extract(, .SD[ninterpro == max(ninterpro)], by = fid_var)
   fdata1 %>%  magrittr::extract(, ninterpro := NULL)
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: prefer best interpro annotation (per protein group):', n0, n1)

   # Collapse isoforms (per canonical accession)
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, (uniprot_var):= paste0(get(uniprot_var), collapse = ';'), by = c(fid_var, 'Canonical accessions')) %>% unique()
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Collapse isoforms (per canonical accession)', n0, n1)

   # Keep those with "gene name" annotation (per protein group)
   n0 <- nrow(fdata1)
   fdata1[, ngene:=0]
   fdata1[`Gene names` != '', ngene:=1]
   fdata1 %<>% magrittr::extract(, .SD[ngene==max(ngene)], by = feature_id)
   fdata1[, ngene:=NULL]
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Prefer entries with gene names (per protein group)', n0, n1)

   # Keep best existence values (per protein group)
   n0 <- nrow(fdata1)
   fdata1[, existence:=as.numeric(existence)]
   fdata1 %<>% magrittr::extract(,.SD[existence==min(existence)], by = c(fid_var))
   fdata1[, existence := NULL]
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Keep entries with best uniprot existence score (per protein group)', n0, n1)

   # Keep only first of redundant entries (per gene name)
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, .SD[1], by = c(fid_var, 'Gene names'))
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Keep only first of redundant uniprot entries  (per protein group and gene name)', n0, n1)

   # Collapse paralogs (per protein group)
   n0 <- nrow(fdata1)
   fdata1 %<>% magrittr::extract(, (uniprot_var)          :=  get(uniprot_var)      %>% paste0(collapse = '; '),                             by = fid_var)
   fdata1 %<>% magrittr::extract(, `Canonical accessions` := `Canonical accessions` %>% paste0(collapse = '; '),                             by = fid_var)
   fdata1 %<>% magrittr::extract(, `Gene names`           := `Gene names`           %>% paste0(collapse = '; '),                             by = fid_var)
   fdata1 %<>% magrittr::extract(,  keggid                :=  keggid                %>% paste0(collapse = '; '),                             by = fid_var)
   fdata1 %<>% magrittr::extract(,  goid                  :=  goid                  %>% autonomics.support::uniquify_collapsed_strings(';'), by = fid_var)
   fdata1 %<>% magrittr::extract(,  interpro              :=  interpro              %>% autonomics.support::uniquify_collapsed_strings(';'), by = fid_var)
   fdata1 %<>% magrittr::extract(, `Protein names`        := `Protein names`        %>% autonomics.support::commonify_strings(),             by = fid_var)
   fdata1 %<>% unique()
   n1 <- nrow(fdata1)
   autonomics.support::cmessage('\t%d -> %d: Collapse paralogs (per protein group) and commonify protein names', n0, n1)

   # Merge back into object and return
   autonomics.support::cmessage('Merge %d annotations into SummarizedExperiment with %d protein groups', n1, nrow(object))
   fdata1 %<>% magrittr::extract(match(autonomics.import::fdata(object)$feature_id, .$feature_id), )
   fdata1 %<>% data.frame(row.names = .$feature_id, check.names = FALSE)
   autonomics.import::fdata(object) <- fdata1
   object
}
