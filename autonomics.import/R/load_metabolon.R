
#' @rdname load_metabolon
#' @importFrom magrittr %>%
#' @export
load_metabolon_sdata <- function(file, sheet = 2){

   # Load sample data
   . <- NULL
   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   df %<>% magrittr::extract(1:fstart, (sstart+1):ncol(.)) %>%
           t() %>%
           data.frame(stringsAsFactors = FALSE)            %>%
           magrittr::set_names(df[[sstart]][1:fstart])     %>%
           magrittr::set_names(names(.) %>% stringi::stri_replace_first_fixed( 'Group   HMDB_ID', 'Group') %>% # recent metabolon files
                                            stringi::stri_replace_first_fixed('Sample   HMDB_ID', 'Group'))    # older metabolon files
   # Return
   df
}

#' @rdname load_metabolon
#' @importFrom magrittr %>%
#' @export
load_metabolon_fdata <- function(file, sheet=2){
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

#' Load metabolon data
#' @param file                         metabolon xlsx file
#' @param sheet                        xls sheet name  or number
#' @param log2_transform               logical: whether to log2 transform
#' @param infer_design  logical: whether to infer design from sample ids
#' @param add_kegg_pathways            logical: whether to add KEGG pathways to fdata
#' @param add_smiles                   logical: whether to add SMILES to fdata
#' @param ... (backward compatibility)
#' @return SummarizedExperiment (load_metabolon) or dataframe (load_metabolon_sdata, load_metabolon_fdata)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/metabolon/subramanian.2016.metabolon.xlsx',
#'                         package = 'subramanian.2016')
#'    file %>% autonomics.import::load_metabolon(sheet=5)
#'    file %>% autonomics.import::load_metabolon_sdata(sheet=5) %>% head(1)
#'    file %>% autonomics.import::load_metabolon_fdata(sheet=5) %>% head(1)
#' }
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% autonomics.import::load_metabolon()
#'    file %>% autonomics.import::load_metabolon(sheet=3)
#'    file %>% autonomics.import::load_metabolon() %>%
#'             autonomics.import::sdata() %>% magrittr::extract(1:2, 1:5)
#'    file %>% autonomics.import::load_metabolon(infer_design = TRUE) %>%
#'             autonomics.import::sdata() %>% magrittr::extract(1:2, 1:5)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon <- function(
   file,
   sheet=2,
   log2_transform     = TRUE,
   infer_design       = FALSE,
   add_kegg_pathways  = FALSE,
   add_smiles         = FALSE
){
   # Satisfy CHECK
   . <- NULL

   # Sheet
   cur_sheet <- readxl::excel_sheets(file) %>%
               (function(x){ names(x) <- x; x}) %>% magrittr::extract2(sheet)
   autonomics.support::cmessage('Load  %s  %s', basename(file), cur_sheet)

   # Get start points
   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]

   # Load components
   sdata1 <- load_metabolon_sdata(file, sheet=sheet)
   fdata1 <- load_metabolon_fdata(file, sheet=sheet)
   exprs1 <- df %>% magrittr::extract((fstart+1):nrow(.), (sstart+1):ncol(.)) %>%
                    data.matrix() %>%
                    magrittr::set_colnames(sdata1$sample_id) %>%
                    magrittr::set_rownames(fdata1$MCOMP_ID)

   # Wrap into SummarizedExperiment
   object <- SummarizedExperiment::SummarizedExperiment(assays=list(exprs = exprs1))
   if (log2_transform)   autonomics.import::exprs(object) %<>% log2()
   autonomics.import::sdata(object)  <- sdata1
   autonomics.import::fdata(object)  <- fdata1
   autonomics.import::prepro(object) <- list(assay='lcms', entity='metabolite', quantity='intensities', software='metabolon')
   autonomics.import::annotation(object) <- ''

   # Standardize design
   object %<>% autonomics.import::prepare_design(sampleid_var = 'CLIENT_IDENTIFIER',
                                                 subgroup_var = 'Group',
                                                 infer_design = infer_design)

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

