
#' @rdname load_metabolon
#' @importFrom magrittr %>%
#' @export
load_metabolon_sdata <- function(file, sheet = 2){
   . <- NULL
   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   df %<>% magrittr::extract(1:fstart, (sstart+1):ncol(.)) %>%
           t() %>%
           data.frame(stringsAsFactors = FALSE) %>%
           magrittr::set_names(df[[sstart]][1:fstart]) %>%
           magrittr::set_rownames(.$CLIENT_IDENTIFIER %>%
                                  autonomics.import::validify_sample_ids()) %>%
           magrittr::set_names(names(.) %>%
                               stringi::stri_replace_first_fixed('Group   HMDB_ID',  'subgroup') %>% # recent metabolon files
                               stringi::stri_replace_first_fixed('Sample   HMDB_ID', 'subgroup'))    # older metabolon files
   cbind(sample_id = df %>% rownames(),
         df %>% magrittr::extract(, 'subgroup', drop = FALSE),
         df %>% magrittr::extract(, setdiff(names(.), 'subgroup'), drop = FALSE))
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
#' @param file              metabolon xlsx file (character)
#' @param sheet             number or name
#' @param verbose           logical
#' @param add_kegg_pathways logical
#' @param add_smiles        logical
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
#' if (require(halama.2016)){
#'    file <- system.file('extdata/halama.2016.xlsx', package = 'halama.2016')
#'    file %>% autonomics.import::load_metabolon(sheet=4)
#'    file %>% autonomics.import::load_metabolon_sdata(sheet=4) %>% head(1)
#'    file %>% autonomics.import::load_metabolon_fdata(sheet=4) %>% head(1)
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon <- function(
   file,
   sheet=2,
   verbose            = TRUE,
   add_kegg_pathways  = FALSE,
   add_smiles         = FALSE
){
   # Satisfy CHECK
   . <- NULL

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
   sumexp1 <- SummarizedExperiment::SummarizedExperiment(assays=list(exprs = exprs1))
   autonomics.import::exprs(sumexp1) %<>% log2()
   autonomics.import::sdata(sumexp1)  <- sdata1
   autonomics.import::fdata(sumexp1)  <- fdata1
   autonomics.import::prepro(sumexp1) <- list(assay='lcms', entity='metabolite', quantity='intensities', software='metabolon')
   autonomics.import::annotation(sumexp1) <- ''

   # Annotate
   if (verbose)           autonomics.support::cmessage('\tAnnotate')
   if (add_kegg_pathways) sumexp1 %<>% autonomics.import::add_kegg_pathways_to_fdata()
   if (add_smiles)        sumexp1 %<>% autonomics.import::add_smiles_to_fdata()

   # Return
   sumexp1
}

#' @rdname load_metabolon
#' @export
load_metabolon_file <- function(...){
   .Deprecated('load_metabolon')
   load_metabolon(...)
}

