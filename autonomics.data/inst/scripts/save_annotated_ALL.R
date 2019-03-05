
#' Load an autonomics-ready version of the ALL data
#' 
#' Annotate the ALL eSet and prepare it for use by autonomics.
#' Save in /data/ALL.RData
#' 
#' @importFrom magrittr        %>%
#' @noRd
save_annotated_ALL <- function(){
   
   getFirstFeatureFromALL <- function(envir){
      AnnotationDbi::mget(autonomics.import::fnames(e$ALL), envir) %>% 
      sapply(getElement, 1) %>% 
      unname()
   }
   
   e <- new.env()
   utils::data(list = 'ALL', package = 'ALL', envir = e)
   e$ALL %<>% SummarizedExperiment::SummarizedExperiment()
   autonomics.import::fdata(e$ALL) <- data.frame(
      gene_symbols       = getFirstFeatureFromALL(hgu95av2.db::hgu95av2SYMBOL),
      gene_names         = getFirstFeatureFromALL(hgu95av2.db::hgu95av2GENENAME),
      ensembl_gene_ids   = getFirstFeatureFromALL(hgu95av2.db::hgu95av2ENSEMBL),
      entrez_gene_ids    = getFirstFeatureFromALL(hgu95av2.db::hgu95av2ENTREZID),
      uniprot_accessions = getFirstFeatureFromALL(hgu95av2.db::hgu95av2UNIPROT),
      row.names          = autonomics.import::fnames(e$ALL), 
      stringsAsFactors = FALSE
   )
   
   # Prepare pData
   e$ALL <- e$ALL[, !is.na(e$ALL$sex)]                                           # rm samples with NA sex
   e$ALL$subgroup <- paste0(gsub('([BT])[1234]?', '\\1', e$ALL$BT), e$ALL$sex)   # add var subgroup
   
   # Prepare fData
   autonomics.import::fdata(e$ALL)$feature_id <- rownames(autonomics.import::fdata(e$ALL))
   e$ALL <- e$ALL[!is.na(autonomics.import::fdata(e$ALL)$gene_symbols), ]                  # rm features with NA gene symbol
   e$ALL <- e$ALL[!duplicated(autonomics.import::fdata(e$ALL)$gene_symbols), ]             # rm features with duplicate gene symbols
   
   # Check and save
   ALL <- e$ALL
   autonomics.import::assert_is_valid_object(ALL)
   usethis::use_data(ALL, compress = 'xz', overwrite = TRUE)
 }
 
require(magrittr)
save_annotated_ALL()