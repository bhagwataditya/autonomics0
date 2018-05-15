
#' @importFrom magrittr %<>% extract
correct_for_universe <- function(pathways, universe, min_features_in_set = 5){
  pathways %<>% lapply(function(x){x[x %in% universe]})
  pathways %<>% magrittr::extract(lapply(., length) > min_features_in_set)
}

utils::globalVariables('TERM')
utils::globalVariables('GOID')

#' Load GO collection
#' @param collection_name   value in GO_COLLECTIONS
#' @param organism          value in ANNOTATED_ORGANISMS
#' @return Named list: names are GO IDs, each element contains a vector of entrez gene IDs
#' @examples 
#' gobp <- load_go_collection('GOBP', 'Homo sapiens')
#' @importFrom data.table      data.table   :=
#' @importFrom magrittr        %>%   %<>%
#' @export
load_go_collection <- function(collection_name, organism){
  
  # Check args
  assertive.sets::assert_is_subset(tolower(collection_name), GO_COLLECTIONS)
  assertive.sets::assert_is_subset(organism, names(ANNOTATED_ORGANISMS))
  
  # Select GOIDS: choose ontology
  GOIDS <- utils::getFromNamespace(sprintf('%sANCESTOR', toupper(collection_name)), 'GO.db') %>% AnnotationDbi::keys()
  
  # Select GOIDS: limit to terms in annotation map
  annotationPackage <- sprintf(
     'org.%s.%s.db',
     ANNOTATED_ORGANISMS[[organism]]['abrev'],
     ANNOTATED_ORGANISMS[[organism]]['source'])
  install_package_if_necessary(annotationPackage)
  GO2ALLIDS <- if(organism == 'yeast'){
     utils::getFromNamespace(
      sprintf(
         'org.%s.%sGO2ALLORFS',
         ANNOTATED_ORGANISMS[[organism]]['abrev'],
         ANNOTATED_ORGANISMS[[organism]]['source']),
     annotationPackage)
  } else {
     utils::getFromNamespace(
        sprintf(
           'org.%s.%sGO2ALLEGS',
           ANNOTATED_ORGANISMS[[organism]]['abrev'],
           ANNOTATED_ORGANISMS[[organism]]['source']),
        annotationPackage)
  }
  GOIDS %<>% magrittr::extract(. %in% AnnotationDbi::keys(GO2ALLIDS))
  
  # Select GOID: exclude top layer terms
  topLayerTerms <- c(BP = 'GO:0008150', CC = 'GO:0005575', MF = 'GO:0003674')
  GOIDS   %<>% magrittr::extract(!(. %in% topLayerTerms))
  suppressMessages(
    GODT <- AnnotationDbi::select(GO.db::GO.db, GOIDS, 'TERM') %>%
            data.table::data.table() %>% 
            magrittr::extract(, TERM := paste0(GOID, ' ', TERM))
  )
  
  # Load GO Collection
  GOCollection <- AnnotationDbi::mget(GODT$GOID, GO2ALLIDS) %>%
                  AnnotationDbi::as.list() %>%
                  magrittr::set_names(GODT$TERM)
  
  # Remove gene names (which are evidence codes) and duplicate entrez ids
  GOCollection <- lapply(GOCollection, function(geneVector){unique(unname(geneVector))})
  
  # Return
  GOCollection
}


#' @noRd
#' @examples
#' load_kegg_pathways_to_entrezg('hsa')
#' load_kegg_pathways_to_entrezg('mmu')
#' @importFrom magrittr %>%
load_kegg_pathways_to_entrezg <- function(organism){
  assertive.sets::assert_is_subset(organism, names(ANNOTATED_ORGANISMS))
  myFile <- sprintf('extdata/kegg_%s_pathways_to_entrezg.rds', ANNOTATED_ORGANISMS[[organism]]['abrev'])  %>%
            system.file(package = 'autonomics.ora')
  suppressWarnings(assertive.files::assert_all_are_readable_files(myFile))
  readRDS(myFile)
}

#' @noRd
#' @examples
#' load_kegg_pathways_to_entrezg('human')
#' load_kegg_pathways_to_entrezg('mouse')
#' @importFrom magrittr %>%
load_kegg_pathways_to_names <- function(organism){
  assertive.sets::assert_is_subset(organism, names(ANNOTATED_ORGANISMS))
  myFile <- sprintf('extdata/kegg_%s_pathways_to_names.rds', ANNOTATED_ORGANISMS[[organism]]['abrev'])  %>%
            system.file(package = 'autonomics.ora')
  suppressWarnings(assertive.files::assert_all_are_readable_files(myFile))
  readRDS(myFile)
}

load_kegg_collection <- function(organism){
  kegg_pathways <- load_kegg_pathways_to_entrezg(organism)
  kegg_names    <- load_kegg_pathways_to_names(organism)
  kegg_terms    <- unname(kegg_names[names(kegg_pathways)]) %>% 
                   stringi::stri_replace_first_fixed(ANNOTATED_ORGANISMS[[organism]]['kegg_pathway_suffix'], '')
  names(kegg_pathways) <- paste0(names(kegg_pathways), ' ', kegg_terms)
  kegg_pathways
}


#' @noRd
#' @importFrom magrittr %>%
load_uniprot_keywords <- function(organism){
  assertive.sets::assert_is_subset(organism, names(ANNOTATED_ORGANISMS))
  myFile <- sprintf('extdata/uniprot_%s_keywords_to_accessions.rds', ANNOTATED_ORGANISMS[[organism]]['kegg'])  %>%
            system.file(package = 'autonomics.ora')
  suppressWarnings(assertive.files::assert_all_are_readable_files(myFile))
  readRDS(myFile)
}

# get_pathways_old <- function(collection, organism, universe = NULL){
#   pathways <- switch(tolower(collection),
#                      gobp    = suppressMessages(EnrichmentBrowser::get.go.genesets(organism, 'BP')),
#                      gomf    = suppressMessages(EnrichmentBrowser::get.go.genesets(organism, 'MF')),
#                      gocc    = suppressMessages(EnrichmentBrowser::get.go.genesets(organism, 'CC')),
#                      kegg    = suppressMessages(EnrichmentBrowser::get.kegg.genesets(organism))
#   )
#   if (!is.null(universe)){
#     pathways %<>% correct_for_universe(universe)
#   }
#   return(pathways)
# 
# }

#' Get pathways
#' @param collection value from GENE_SET_COLLECTIONS
#' @param organism value from ANNOTATED_ORGANISMS
#' @param universe either NULL or a vector with feature IDs defining the universe
#' @param min_set_size minimum number of features in a set
#' @param max_set_size maximum number of features in a set
#' @importFrom magrittr   %<>%
#' @export
load_collection <- function(collection, organism, universe = NULL, min_set_size = 5, max_set_size = 2000){
  pathways <- switch(tolower(collection), 
                     gobp    = load_go_collection('GOBP', organism),
                     gomf    = load_go_collection('GOMF', organism),
                     gocc    = load_go_collection('GOCC', organism),
                     kegg    = load_kegg_collection(organism), 
                     uniprot = load_uniprot_keywords(organism))
  
  n_features <- vapply(pathways, length, integer(1))
  autonomics.support::cmessage('\t\tloading sets with %s-%s features', min_set_size, max_set_size)
  pathways %<>% magrittr::extract(n_features > min_set_size   &   n_features < max_set_size)
  
  if (!is.null(universe)){
    pathways %<>% correct_for_universe(universe)
  }
  return(pathways)
}

