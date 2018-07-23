


#===========================================
# GENERIC
#===========================================

#' Load omics
#' @param file path to omics data file
#' @param platform 'metabolon', 'metabolonlipids'
#' @param sheet excel sheet number or name if applicable
#' @param quantity string: which quantity to extract into exprs
#' @param design_file path to design file
#' @param log2_transform logical
#' @param log2_offset offset in mapping x -> log2(x+offset)
#' @param infer_design_from_sampleids logical
#' @param design_sep string: design separator
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'
#' # MAXQUANT: STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcomp/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    object <- file %>% load_omics(platform = 'maxquant',
#'                                  quantity = 'Ratio normalized',
#'                                  infer_design_from_sampleids = TRUE)
#' }
#'
#' # MAXQUANT: STEM CELL DIFFERENTIATION
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
                                                      design_sep                  = design_sep,
                                                      sheet                       = sheet)
   object %<>% autonomics.import::merge_sdata(design_df, by = sampleid_varname(platform))
   if (!is.null(design_file)){
      file_df <- autonomics.import::read_design(design_file)
      object %<>% autonomics.import::merge_sdata(file_df, by = sampleid_varname(platform))
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
#'    sdata(object) %>% head()
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
   log2_offset                 = 0
){

   # Load exprs
   object <- autonomics.import::load_omics(file                        = file,
                                           platform                    = 'maxquant',
                                           quantity                    = quantity,
                                           log2_transform              = log2_transform,
                                           log2_offset                 = log2_offset,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep,
                                           design_file                 = design_file)

   # Add prepro info
   autonomics.import::prepro(object) <- list(assay    = 'lcms',
                                             entity   = 'proteingroup',
                                             quantity = quantity,
                                             software = 'maxquant')

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
#'              autonomics.import::load_proteingroups()
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
#' @param ... (backward compatibility)
#' @return SummarizedExperiment (load_metabolon) or dataframe (load_sdata_metabolon, load_fdata_metabolon)
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
