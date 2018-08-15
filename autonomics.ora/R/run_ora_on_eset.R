
# # Version which just takes top 200
# #' @importFrom magrittr extract %<>% %>%
# which_query <- function(object, contrast_name, direction = 'both', top = 200){
#    # Order
#    pvar <- paste0('p.', contrast_name)
#    idx <- order(autonomics.import::fdata(object)[[pvar]])
#
#    # Restrict direction
#    diffvar <- paste0('diff.', contrast_name)
#    if (direction == 'neg'){
#       idx %<>% extract(autonomics.import::fdata(object)[[diffvar]] < 0
#    } else if (direction == 'pos'){
#      idx %<>% extract(autonomics.import::fdata(object)[[diffvar]] > 0)
#    }
#
#    # Take top
#    idx %<>% extract(1:top)
#    return(idx)
# }

#' Get entrez gene ids from uniprot
#' @param  object eset
#' @param  load_org.xx.xx.db   load_org.xx.xx.db
#' @return entrezg ids
#' @examples
#' require(magrittr)
#' if (require(autonomics.data) & require(org.Hs.eg.db)){
#'    object <- autonomics.data::billing2016
#'    object[1:10, ] %>% get_entrezg_from_uniprot(org.Hs.eg.db::org.Hs.eg.db)
#' }
#' if (require(atkin.2014) & require(org.Hs.eg.db)){
#'    object <- atkin.2014::soma
#'    object[1:10, ] %>% get_entrezg_from_uniprot(org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @importFrom magrittr   %>%
#' @export
get_entrezg_from_uniprot <- function(object, load_org.xx.xx.db){

  # Assert
  autonomics.import::assert_is_valid_eset(object)

  # Prepare uniprot accessions:
  #    - uncollapse
  #    - only uniprot accessions with a mapping in load_org.xx.xx.db
  #    - remove PG with no mappings
  uniprot_keys <- load_org.xx.xx.db %>% AnnotationDbi::keys(keytype = 'UNIPROT')
  my_uniprot <- autonomics.import::fdata(object) %>%
                magrittr::extract2(autonomics.import::uniprot_var(object))   %>%
                as.character()                                          %>%
                strsplit(autonomics.import::sep(object))           %>%
                lapply(function(x){x %>% magrittr::extract(. %in% uniprot_keys)}) %>%
                magrittr::extract(vapply(.,function(x) length(x) > 0, logical(1)))

  # Map uniprot accessions -> entrezg (for each PG)
  #   - isoforms map to a single gene.
  #   - paralogs map to multiple genes.
  #        simplest explanation is that each PG represents only one paralog
  #        which one is unknown
  #        we want to select only one to prevent false positive of ora results
  #           e.g. Anja's differentiation exp: histone modification (HIST1, HIST2, HIST3)
  #        simplest solution: take first paralog
  dbname <- load_org.xx.xx.db$packageName
  bimapname <- dbname %>% stringi::stri_replace_first_fixed('.db', 'UNIPROT')
  uniprot2eg <- utils::getFromNamespace(bimapname, dbname) %>% AnnotationDbi::revmap()
  my_entrezg <- my_uniprot %>%
                vapply(function(x){
                          AnnotationDbi::mget(x, uniprot2eg) %>%
                          unlist() %>%
                          unname() %>%
                          unique() %>%
                          magrittr::extract(1)},
                        character(1))
  unique(my_entrezg)

  # Old approach was to consider all possible uniprot -> entrezg mappings
  # Now discarded to prevent false positive ora results
  # autonomics.import::fdata(object)$uniprot_accessions %>%
  # uncollapse_and_unite(sep = sep) %>%
  # uniprot_set_to_entrezg_set(load_org.xx.xx.db)
}

#' Get entrez from gene symbols
#' @param  object       eset
#' @param  load_org.xx.xx.db   load_org.xx.xx.db
#' @return vector with entrez gene ids
#' @examples
#' require(magrittr)
#' if (require(autonomics.data) & require(org.Hs.eg.db)){
#'    object <- autonomics.data::billing2016
#'    load_org.xx.xx.db <- org.Hs.eg.db::org.Hs.eg.db
#'    object %>% get_entrezg_from_gsymbols(load_org.xx.xx.db)
#' }
#' if (require(atkin.2014) & require(org.Hs.eg.db)){
#'    atkin.2014::soma[1:10, ] %>%
#'    get_entrezg_from_gsymbols(org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @importFrom magrittr   %>%
#' @export
get_entrezg_from_gsymbols <- function(object, load_org.xx.xx.db){
  # Assert
  autonomics.import::assert_is_valid_eset(object)

  # Uncollapse gsymbols if required
  gsymbols <- autonomics.import::fdata(object) %>%
              magrittr::extract2(autonomics.import::fname_var(object)) %>%
              as.character()
  sep <- autonomics.import::sep(object)
  if (!is.null(sep)){
     gsymbols %<>% strsplit(sep) %>% vapply(extract, character(1), 1)
  }
  AnnotationDbi::mapIds(load_org.xx.xx.db, keys = gsymbols, column = 'ENTREZID', keytype = 'SYMBOL')

  # Map gsymbols to entrezg
  gsymbols %>% gsymbol_set_to_entrezg_set(load_org.xx.xx.db)
}

#' Get entrezg from eset
#' @param object eset
#' @param load_org.xx.xx.db   load_org.xx.xx.db
#' @return vector with entrez gene ids
#' @examples
#' require(magrittr)
#' if (require(autonomics.data) & require(org.Hs.eg.db)){
#'    object <- autonomics.data::billing2016[1:10, ]
#'    object %>% get_entrezg(org.Hs.eg.db::org.Hs.eg.db)
   #' }
#' if (require(atkin.2014) & require(org.Hs.eg.db)){
#'    atkin.2014::soma[1:10, ] %>%
#'    get_entrezg(org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @export
get_entrezg <- function(object, load_org.xx.xx.db){

  # Shotgun proteomics: map from uniprot accessions. (gene symbols not always available, cfr. Trypanosoma)
   if (autonomics.import::is_maxquant_eset(object)){
      object %>% get_entrezg_from_uniprot(load_org.xx.xx.db)
   }
   # somscan: gene symbols (available)
   else if (autonomics.import::is_soma_eset(object)){
      object %>% get_entrezg_from_gsymbols(load_org.xx.xx.db)
   }
  # RNA seq: map from gene symbols
  else if (Biobase::annotation(object) %in% 'cufflinks_gene_counts'){
     stop('implement autonomics.ora::get_entrezg for this type of eset')
    # object %>% get_entrezg_from_gsymbols(load_org.xx.xx.db)
  }
  # Microarrays: map from gene symbols
  else {
     stop('implement autonomics.ora::get_entrezg for this type of eset')
    # object %>% get_entrezg_from_gsymbols(load_org.xx.xx.db)
  }
}

# universe: 'detectome'  or 'genome'
#' @importFrom magrittr   %>%
write_ora_results <- function(ora_results, contrast, direction, collection, universe, result_dir){

  # Create result dir
  ora_dir <- sprintf('%s/ora_in_%s',
                     autonomics.find::get_contrast_subdir(result_dir, names(contrast)),
                     universe)
  dir.create(ora_dir, showWarnings = FALSE, recursive = TRUE)

  # Write ora results
  file <- sprintf('%s/%s_%s_%s_%s.txt', ora_dir, collection, names(contrast), direction, universe)
  ora_results %>% autonomics.support::print2txt(file)

  # Report
  autonomics.support::cmessage('\t\t%s %s 0   %s',
                               names(contrast),
                               c(neg='< ', pos=' >', both = '<>')[[direction]],
                               basename(file))
}

utils::globalVariables('annotate_ora')

#' Run ora on sumexp
#' @param object SummarizedExperiment
#' @param top_definition definition of 'top_features'
#' @param oraset_var ora set variable
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    autonomics.import::fvars(object)
#' }
run_ora_on_sumexp <- function(
  object, 
  top_definition = autonomics.find::default_top_definition(object), 
  oraset_var = 'goid'
){
  
  # Assert
  autonomics.import::assert_is_valid_eset(object)
  assertive.sets::assert_is_subset(oraset_var, autonomics.import::fvars(object))
  
  # Return if no limma info in sumexp
  if (!autonomics.find::contains_limma_in_fdata(object)){
    message('Abort - run object %<>% autonomics.find::add_limma_to_fdata() before run_ora_on_sumexp')
  }
  
  # Assert
  autonomics.import::assert_is_valid_eset(object)
  
  # Infer contrast names
  contrast_names <- object %>% autonomics.find::infer_contrast_names()
  
  # Run for each contrast
  universe <- object %>% autonomics.ora::extract_ora_universe()
  pathway_list <- object %>% autonomics.ora::extract_ora_sets(oraset_var)
  for (cur_contrast in contrast_names){
    for (direction in c('pos', 'neg', 'both')){
      query    <- object %>% autonomics.ora::extract_ora_query(cur_contrast, top_definition, direction)
      ora_res <- autonomics.ora::run_ora(
                    query                   = query, 
                    universe                = universe, 
                    pathway_list            = pathway_list, 
                    return_only_significant = TRUE)
      ora_res %>% annotate_ora
      utils::str(ora_res)
    }
  }
  
}


#' Run over representation analysis
#' @param object eset as returned by \code{\link[autonomics.find]{add_limma_to_fdata}}
#' @param contrasts                 contrasts vector
#' @param result_dir                directory to which to write results
#' @param gene_set_collections      gene set collections used for analysis (see \code{\link{GENE_SET_COLLECTIONS}} for implemented ones)
#' @param top_definition            definition of 'top features'.
#' @param universe                 'detectome' (all detected features), 'genome' (all known features),
#                                  or NULL (no ora)
#' @param min_set_size              minimum set size
#' @param max_set_size              maximum set size
#' @param ora_detectome_in_genome   whether to perform a detectome against genome ora (logical)
#' @importFrom magrittr  %>%
#' @export
run_ora_on_eset <- function(
  object,
  contrasts                    = autonomics.find::default_contrasts(object),
  result_dir,
  gene_set_collections         = GENE_SET_COLLECTIONS,
  top_definition               = autonomics.find::default_top_definition(object),
  universe                     = default_universe(object),
  min_set_size                 = 5,
  max_set_size                 = Inf,
  ora_detectome_in_genome = default_ora_detectome_in_genome()
){
   # Return NULL if no ora can be performed
   if (is.null(universe)){
     autonomics.support::cmessage('\t\tAbort ora: universe==NULL')
     return(NULL)
   }
   if (nrow(object)<50){
      autonomics.support::cmessage('\t\tAbort ora: nrow(object) < 50')
      return(NULL)
   }

   # Assert valid inputs
   autonomics.import::assert_is_valid_eset(object)
   #assertive.files::assert_all_are_writable_files(result_dir, warn_about_windows = FALSE)
   assertive.sets::assert_is_subset(gene_set_collections, GENE_SET_COLLECTIONS)
   assertive.sets::assert_is_subset(universe, c('detectome', 'genome'))
   assertive.types::assert_is_numeric(min_set_size)
   assertive.types::assert_is_numeric(max_set_size)

   # Constants
   organism <- object %>% infer_organism()
   autonomics.support::cmessage('\tUse %s gene set collections for ora', names(organism))
   load_org.xx.xx.db <- organism %>% load_org.xx.xx.db()
   entrezg_universe   <- get_entrezg(object, load_org.xx.xx.db = load_org.xx.xx.db)
   entrezg_multiverse <- AnnotationDbi::keys(load_org.xx.xx.db)

   # Loop across collection, contrast, direction
   return_list_a <- list()
   for (collection in gene_set_collections){
      autonomics.support::cmessage('\tora in %s %s', collection, universe)
      pathway_list <- load_collection(collection, organism = organism, min_set_size = min_set_size, max_set_size = max_set_size)

      # Is detectome universe enriched against genome multiverse?
      if (universe == 'detectome' & ora_detectome_in_genome){
        file <- sprintf('%s/ora_of_detectome_in_%s_genome.txt', result_dir, collection)
        autonomics.support::cmessage('\t\tdetectome   %s', basename(file))
        ora_res <- run_ora(entrezg_universe, entrezg_multiverse, pathway_list, min_set_size = min_set_size)
        # ora_res$gene_symbols <- ora_res$features %>% collapsed_entrezg_to_gsymbol(load_org.xx.xx.db)
        ora_res %>% autonomics.support::print2txt(file)
      }

      # Are contrasts enriched compared to universe?
      if (universe == 'genome')    entrezg_universe <- entrezg_multiverse
      return_list_b <- list()
      for (i in seq_along(contrasts)){
         return_list_c <- list()
         for (direction in c('neg', 'pos', 'both')){
           query_eset <- object %>% autonomics.find::filter_n_arrange_top_features(names(contrasts[i]), top_definition, direction)
           if (nrow(query_eset) < min_set_size)  break()
           entrezg_query <- query_eset %>% get_entrezg(load_org.xx.xx.db = load_org.xx.xx.db)
           if (length(entrezg_query) < min_set_size)  break()
           ora_results <- run_ora(entrezg_query, entrezg_universe, pathway_list, min_set_size = min_set_size)
           # write_ora_results(ora_results, query_eset, contrasts[i], direction, collection, result_dir)
           write_ora_results(ora_results, contrasts[i], direction, collection, universe, result_dir)
           return_list_c[[direction]] <- ora_results %>%
              dplyr::mutate(
                 direction =  direction,
                 contrast  =  contrasts[i],
                 collection = collection)
         }
         return_list_b[[contrasts[i]]] <- return_list_c
      }
      return_list_a[[collection]] <- return_list_b
   }
   
   # Gather into data.frame
   return_list_a %<>%
      unlist(recursive = FALSE) %>%
      unlist(recursive = FALSE) %>%
      dplyr::bind_rows()
   
   invisible(return_list_a)
}

