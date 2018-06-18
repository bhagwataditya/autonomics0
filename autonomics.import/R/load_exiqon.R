#' @export
EXIQON_PATTERN <- '^(mmu|rno)[-](miR|let)[-]([0-9]+[a-z]?[-]?[0-9]?)([-][35]p)?$'

#' @export
EXIQON_FEATURE_ROWS <- c('#RefGenes', '#Spike',        '#MissDataFreq(%)')

#' @export
EXIQON_SAMPLE_COLS  <- c('Exiqon',    '#IPC-PlateID',  '#MissDataFreq(%)')

#' @rdname load_exiqon
#' @importFrom magrittr %>%
#' @export
load_exiqon_sdata <- function(file){

   # Satisfy CHECK
   . <- NULL

   # Load from file
   mir <- readxl::read_excel(file) %>% as.data.frame() %>% magrittr::set_rownames(.$Exiqon)

   # Define feature rows and sample cols
   feature_rows <- rownames(mir) %in% EXIQON_FEATURE_ROWS
   sample_cols  <- names(mir)    %in% EXIQON_SAMPLE_COLS

   # sdata
   sdata1 <- mir %>%
             magrittr::extract(!feature_rows, sample_cols) %>%
             magrittr::set_names(names(.) %>% stringi::stri_replace_first_fixed('Exiqon', 'sample_id'))
   sdata1
}


#' @rdname load_exiqon
#' @importFrom magrittr %>%
#' @export
load_exiqon_fdata <- function(file){

   # Satisfy CHECK
   . <- NULL

   # Load from file
   mir <- readxl::read_excel(file) %>% as.data.frame() %>% magrittr::set_rownames(.$Exiqon)

   # Define feature rows and sample cols
   feature_rows <- rownames(mir) %in% EXIQON_FEATURE_ROWS
   sample_cols  <- names(mir)    %in% EXIQON_SAMPLE_COLS

   # Extract fdata
   mir %>%
   magrittr::extract(feature_rows, !sample_cols) %>%
   t() %>%
   data.frame(check.names = FALSE, stringsAsFactors = FALSE) %>%
   cbind(miRNA = rownames(.), .)

   # %>%
   # dplyr::mutate(mir.org    = fid %>% stringi::stri_replace_all_regex(EXIQON_PATTERN, '$1'),
   #               mir.name   = fid %>% stringi::stri_replace_all_regex(EXIQON_PATTERN, '$2-$3'),
   #               mir.strand = fid %>% stringi::stri_replace_all_regex(EXIQON_PATTERN, '$4') %>% stringi::stri_replace_first_fixed('-', '')) %>%
   # magrittr::set_rownames(.$fid)
}

#' Filter mir features
#' @param object     eset
#' @param rm_ref_genes whether to remove ref genes (logical)
#' @param rm_spike     whether to remove spike-ins (logical)
#' @param verbose      whether to report progress
#' @return filtered eset
#' @importFrom magrittr   %<>%
#' @export
filter_exiqon_features <- function(object, rm_ref_genes = TRUE, rm_spike = TRUE, verbose = TRUE){

   if (rm_ref_genes){
      n0 <- nrow(object)
      object %<>% autonomics.import::filter_features_("`#RefGenes`==0")
      autonomics.import::fdata(object)$`#RefGenes` <- NULL
      n1 <- nrow(object)
      autonomics.support::cmessage('Retain %d/%d miRNAs after removing ref genes', n1, n0)
   }

   if (rm_spike){
      n0 <- nrow(object)
      object %<>% autonomics.import::filter_features_("`#Spike`==0")
      autonomics.import::fdata(object)$`#Spike` <- NULL
      n1 <- nrow(object)
      autonomics.support::cmessage('Retain %d/%d miRNAs after removing spike-ins', n1, n0)
   }
   object
}


#' Load exiqon file
#' @param file         path to exiqon txt file
#' @param infer_design logical: whether to infer design from samplesids
#' @param rm_ref_genes logical: whether to rm reference genes
#' @param rm_spike     logical: whether to rm spike-ins
#' @return eset
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                         package = 'subramanian.2016')
#'    file %>% autonomics.import::load_exiqon()
#'    file %>% autonomics.import::load_exiqon(infer_design = FALSE)
#'    file %>% autonomics.import::load_exiqon_sdata() %>% head(1)
#'    file %>% autonomics.import::load_exiqon_fdata() %>% head(1)
#' }
#' @importFrom magrittr %>%
#' @export
load_exiqon <- function(
   file,
   infer_design = TRUE,
   rm_ref_genes = TRUE,
   rm_spike = TRUE
){
   # Satisfy CHECK
   . <- NULL

  # Load from file
   mir <- readxl::read_excel(file) %>% as.data.frame() %>% magrittr::set_rownames(.$Exiqon)

  # Define feature rows and sample cols
  feature_rows <- rownames(mir) %in% autonomics.import::EXIQON_FEATURE_ROWS
  sample_cols  <- names(mir)    %in% autonomics.import::EXIQON_SAMPLE_COLS

  # Create components
  feature_df   <- autonomics.import::load_exiqon_fdata(file)
  sample_df    <- autonomics.import::load_exiqon_sdata(file)
  exprs_mat    <- mir %>% magrittr::extract(!feature_rows, !sample_cols) %>% data.matrix() %>% t()

  # Create eset
  my_eset <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs_mat))
  autonomics.import::sdata(my_eset) <- sample_df
  autonomics.import::fdata(my_eset) <- feature_df
  autonomics.import::assert_is_valid_eset(my_eset)
  autonomics.import::prepro(my_eset) <- autonomics.import::create_prepro_list(
     assay      = 'exiqon',
     entity     = 'mirna',
     quantity   = 'ct',
     software   = 'genex')

  # Standardize design
  my_eset %>% autonomics.import::prepare_design(sampleid_var = 'sample_id', infer_design = infer_design)

  # Filter features
  my_eset %<>% autonomics.import::filter_exiqon_features(rm_ref_genes = rm_ref_genes,
                                                         rm_spike     = rm_spike)
  # Return
  return(my_eset)
}


load_exiqon_file <- function(...){
   .Deprecated('load_exiqon')
   load_exiqon(...)
}
