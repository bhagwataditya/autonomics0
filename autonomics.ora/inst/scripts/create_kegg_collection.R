# Project: comics
# 
# Author: Aditya Bhagwat
###############################################################################

# Load packages
require(KEGGREST)
require(magrittr)

## Origin: https://stat.ethz.ch/pipermail/bioconductor/2013-September/054902.html
##         http://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r

#' Get KEGG pathway IDs
#' @examples
#' get_kegg_pathway_ids('Homo sapiens')
#' get_kegg_pathway_ids('Mus musculus')
get_kegg_pathway_ids <- function(organism){
  KEGGREST::keggList('pathway', ANNOTATED_ORGANISMS[[organism]]['kegg']) %>%
  (function(pathways.list){
    sub("path:", "", names(pathways.list))
  })
}


#'  Map single kegg pathway id to constituing entrez gene ids
#'  @examples
#'  single_kegg_pathway_to_entrezg('hsa00020'); # single_kegg_pathway_to_entrezg('hsaUnknown')
#'  single_kegg_pathway_to_entrezg('mmu05416'); # single_kegg_pathway_to_entrezg('hsaUnknown')
single_kegg_pathway_to_entrezg <- function(pwid, verbose = TRUE){
  pw <- tryCatch(
           { if(verbose) message(sprintf('%s: OK', pwid))
             pw <- KEGGREST::keggGet(pwid)
             return(pw[[1]]$GENE[c(TRUE, FALSE)]) # subsetting by c(TRUE, FALSE) -- which repeats as many times as needed, sorting
                                                  # through some unexpected packaging of geneIDs in the GENE element of each pw[[n]]
           }, 
           error = function(cond){
              if (verbose) message(sprintf('%s: %s', pwid, cond))
              return(NA)
           }, 
           warning = function(cond){
              if (verbose) message(sprintf('%s: %s', pwid, cond))
              return(NA)
           }
  )
  return(pw)
}


#' Map multiple kegg pathway ids to their entrez gene ids
#' @examples 
#' kegg_pathways_to_entrezg(c("hsa05414", "hsa05416"))
#' kegg_pathways_to_entrezg(c("mmu05414", "mmu05416"))
kegg_pathways_to_entrezg <- function(kegg_pathway_ids){
  message('\nMapping kegg pathway ids to entrez gene ids')
  dt <- system.time({
    map_of_kegg_pathways_to_entrezg <- sapply(kegg_pathway_ids, single_kegg_pathway_to_entrezg)
  })
  sprintf('Done in %s min', ceiling(dt[['elapsed']]/60))
  attr(map_of_kegg_pathways_to_entrezg, 'date')       <- lubridate::today()
  map_of_kegg_pathways_to_entrezg
}


#' Map single kegg pathway id to its name
#' @examples
#' single_kegg_pathway_to_name('hsa00020')
#' single_kegg_pathway_to_name('mmu05416')
single_kegg_pathway_to_name <- function(pwid, verbose = TRUE){
  pw <- tryCatch({ if(verbose) message(sprintf('%s: OK', pwid))
                   pw <- KEGGREST::keggGet(pwid)
                   return(pw[[1]]$NAME)
                 }, error = function(cond){
                               if (verbose) message(sprintf('%s: %s', pwid, cond))
                               return(NA)
                 }, warning = function(cond){
                                if (verbose) message(sprintf('%s: %s', pwid, cond))
                                return(NA)
                 }
  )
  return(pw)
}


#' Map multiple kegg pathway ids to their names
#' @examples
#' kegg_pathways_to_names(c("hsa05414", "hsa05416"))
#' kegg_pathways_to_names(c("mmu05414", "mmu05416"))
kegg_pathways_to_names <- function(kegg_pathway_ids){
  message('\nMapping kegg pathway ids to pathway names')
  dt <- system.time({
    map_of_kegg_pathways_to_names   <- sapply(kegg_pathway_ids, single_kegg_pathway_to_name)
  })
  sprintf('Done in %s min', ceiling(dt[['elapsed']]/60)) %>% message()
  attr(map_of_kegg_pathways_to_names, 'date') <- lubridate::today()
  map_of_kegg_pathways_to_names
}


#' Run KEGG mappings and save
run_kegg_mappings_and_save <- function(organism){
  kegg_pathway_ids <- get_kegg_pathway_ids(organism)
  kegg_pathways_to_entrezg(kegg_pathway_ids) %>%
     saveRDS(
        file = sprintf(
           'inst/extdata/kegg_%s_pathways_to_entrezg.rds',
           ANNOTATED_ORGANISMS[[organism]]['abrev']))
  kegg_pathways_to_names(kegg_pathway_ids) %>%
     saveRDS(
        file = sprintf(
           'inst/extdata/kegg_%s_pathways_to_names.rds',
           ANNOTATED_ORGANISMS[[organism]]['abrev']))
}

## Loop over relevant organisms
Map(run_kegg_mappings_and_save, names(ANNOTATED_ORGANISMS))
