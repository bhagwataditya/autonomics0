
#' Load metabolon sdata
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_metabolon_sdata(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon_sdata <- function(file, sheet){

   # Load sample data
   . <- NULL
   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   df %<>% magrittr::extract(1:fstart, (sstart+1):ncol(.)) %>%
           t() %>%
           data.frame(stringsAsFactors = FALSE, check.names = FALSE)            %>%
           magrittr::set_names(df[[sstart]][1:fstart])     %>%
           magrittr::set_names(names(.) %>% stringi::stri_replace_first_fixed( 'Group   HMDB_ID', 'Group') %>% # recent metabolon files
                                            stringi::stri_replace_first_fixed('Sample   HMDB_ID', 'Group'))    # older metabolon files
   df %<>% magrittr::set_rownames(.$CLIENT_IDENTIFIER)

   # Return
   df
}


#' Load metabolon fdata
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_metabolon_fdata(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon_fdata <- function(file, sheet){
   . <- NULL
   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   df %>% magrittr::extract((fstart+1):nrow(.), 1:sstart)                     %>%
          as.data.frame(stringsAsFactors = FALSE)                             %>%
          magrittr::set_names(df[fstart, 1:sstart] %>% unlist() %>% unname()) %>%
          magrittr::set_rownames(paste0('M', .$COMP_ID))                      %>%
          cbind(MCOMP_ID = rownames(.), .)
}


#' Load metabolon exprs
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_metabolon_exprs(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon_exprs <- function(file, sheet){

   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   sdata1 <- file %>% autonomics.import::load_metabolon_sdata(sheet = sheet)
   fdata1 <- file %>% autonomics.import::load_metabolon_fdata(sheet = sheet)

   df %>% magrittr::extract((fstart+1):nrow(.), (sstart+1):ncol(.)) %>%
          data.matrix() %>%
          magrittr::set_colnames(sdata1$CLIENT_IDENTIFIER) %>%
          magrittr::set_rownames(fdata1$MCOMP_ID)
}


#' Load metabolon data
#' @param file                metabolon xlsx file
#' @param design_file         NULL or character (sample design file)
#' @param sheet               xls sheet name  or number
#' @param log2_transform      logical: whether to log2 transform
#' @param infer_design_from_sampleids        logical: whether to infer design from sample ids
#' @param add_kegg_pathways   logical: whether to add KEGG pathways to fdata
#' @param add_smiles          logical: whether to add SMILES to fdata
#' @param ... (backward compatibility)
#' @return SummarizedExperiment (load_metabolon) or dataframe (load_metabolon_sdata, load_metabolon_fdata)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'
#'    # Loading metabolon file is easy
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% autonomics.import::load_metabolon()
#'
#'    # Three ways to specify sample design
#'       # Use Group definition in metabolon file
#'       file %>% autonomics.import::load_metabolon() %>%
#'                          autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#'       # Infer from sample id values
#'       file %>% autonomics.import::load_metabolon(infer_design_from_sampleids = TRUE) %>%
#'                          autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#'       # Merge in from sample file
#'       design_file <- tempfile()
#'       file %>% autonomics.import::write_metabolon_design(design_file = design_file, infer_from_sampleids = TRUE)
#'       file %>% autonomics.import::load_metabolon(design_file = design_file) %>%
#'                          autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon <- function(
   file,
   sheet = 2,
   design_file                 = NULL,
   log2_transform              = TRUE,
   infer_design_from_sampleids = FALSE,
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
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           sampleid_var                = 'CLIENT_IDENTIFIER')
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

#' @rdname load_metabolon
#' @export
load_metabolon_file <- function(...){
   .Deprecated('load_metabolon')
   load_metabolon(...)
}

