
#' Load an autonomics-ready version of the ALL data
#' 
#' Annotate the ALL eSet and prepare it for use by autonomics.
#' Save in /data/ALL.RData
#' 
#' @importFrom magrittr        %>%
#' @noRd
save_annotated_ALL <- function(){
   
   getFirstFeatureFromALL <- function(envir){
      AnnotationDbi::mget(Biobase::featureNames(e$ALL), envir) %>% 
         sapply(getElement, 1) %>% 
         unname()
   }
   
   e <- new.env()
   utils::data(list = 'ALL', package = 'ALL', envir = e)
   Biobase::fData(e$ALL) <- data.frame(
      gene_symbols       = getFirstFeatureFromALL(hgu95av2.db::hgu95av2SYMBOL),
      gene_names         = getFirstFeatureFromALL(hgu95av2.db::hgu95av2GENENAME),
      ensembl_gene_ids   = getFirstFeatureFromALL(hgu95av2.db::hgu95av2ENSEMBL),
      entrez_gene_ids    = getFirstFeatureFromALL(hgu95av2.db::hgu95av2ENTREZID),
      uniprot_accessions = getFirstFeatureFromALL(hgu95av2.db::hgu95av2UNIPROT),
      row.names = Biobase::featureNames(e$ALL), 
      stringsAsFactors = FALSE
   )
   
   # Prepare pData
   e$ALL <- e$ALL[, !is.na(e$ALL$sex)]                                           # rm samples with NA sex
   e$ALL$subgroup <- paste0(gsub('([BT])[1234]?', '\\1', e$ALL$BT), e$ALL$sex)   # add var subgroup
   
   # Prepare fData
   Biobase::fData(e$ALL)$feature_id <- rownames(Biobase::fData(e$ALL))
   e$ALL <- e$ALL[!is.na(Biobase::fData(e$ALL)$gene_symbols), ]                  # rm features with NA gene symbol
   e$ALL <- e$ALL[!duplicated(Biobase::fData(e$ALL)$gene_symbols), ]             # rm features with duplicate gene symbols
   
   # Add preprocessing info
   autonomics.import::prepro(e$ALL) <- autonomics.import::create_prepro_list(
                                          assay     = 'microarray',
                                          entity    = 'rna',
                                          quantity  = 'abundance', 
                                          software  = 'affymetrix', 
                                          parameters = list())

   # Check and save
   ALL <- e$ALL
   autonomics.import::assert_is_valid_eset(ALL)
   save(ALL, file = 'data/ALL.RData', compress = TRUE)
 }
 
require(magrittr)
save_annotated_ALL()