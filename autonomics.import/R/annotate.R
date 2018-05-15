#' Map a vector of KEGG IDs to KEGG Pathways
#' @param object SummarizedExperiment
#' @return vector with pathway IDs (collapsed)
#' @examples
#' require(magrittr)
#' \dontrun{
#'    if (require(subramanian.2016)){
#'       subramanian.2016::metabolon %>%
#'       autonomics.import::add_kegg_pathways_to_fdata()
#'    }
#'    if (require(halama.2016)){
#'       halama.2016::cell.metabolites %>%
#'       autonomics.import::add_kegg_pathways_to_fdata()
#'    }
#' }
#' @importFrom magrittr %>%
#' @export
add_kegg_pathways_to_fdata <- function(object){

   # Satisfy CHECK
   KEGG <- KEGGPATHWAY <- . <- NULL

   # KEGGIDS
   object %<>% autonomics.import::uncollapse_n_pick_first('KEGG', ',')
   keggids <-  autonomics.import::fdata(object)$KEGG %>%
      unique() %>%
      setdiff(NA_character_)

   # KEGGREST URL
   keggurl <-  keggids %>%
      paste0(collapse = '+') %>%
      sprintf('http://rest.kegg.jp/link/pathway/%s', .)

   # Download KEGGREST table
   cachefile <- tempfile()
   utils::download.file(keggurl, cachefile, quiet = TRUE)

   # Read
   keggdt <- cachefile %>%
      data.table::fread(header = FALSE, col.names = c("KEGG", "KEGGPATHWAY")) %>%
      magrittr::extract(, KEGG        := KEGG        %>% stringi::stri_replace_first_fixed('cpd:',  '')) %>%
      magrittr::extract(, KEGGPATHWAY := KEGGPATHWAY %>% stringi::stri_replace_first_fixed('path:', '')) %>%
      magrittr::extract(, list(KEGGPATHWAY = paste0(KEGGPATHWAY, collapse = ';')), by = 'KEGG')

   # Merge into fdata
   autonomics.import::fdata(object) %<>% autonomics.support::left_join_keeping_rownames(keggdt, by = 'KEGG')

   # Report
   autonomics.support::cmessage('\t\tAdd KEGG Pathways: %3d features -> %3d map to KEGG IDS -> %3d map to KEGG Pathways',
                                nrow(object),
                                sum(autonomics.import::fdata(object)$KEGG!=''),
                                sum(!is.na(autonomics.import::fdata(object)$KEGGPATHWAY)))

   # Return
   object
}


#' Map pubchem IDs to canonical smiles
#'
#' Uses \href{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial}{PUBCHEM's REST interface}
#'
#' @param object SummarizedExperiment
#' @return named character vector (elements = SMILES, names = PUBCHEMIDS)
#' @importFrom magrittr %>%
#' @export
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>%
#'       autonomics.import::add_smiles_to_fdata()
#' }
#'
add_smiles_to_fdata <- function(object){

   # Satify CHECK
   . <- NULL

   # Assert
   assertive.sets::assert_is_subset('PUBCHEM', autonomics.import::fvars(object))

   # Restrict to first mapped PUBCHEM ID
   object %<>% autonomics.import::uncollapse_n_pick_first('PUBCHEM', ';')

   # Map to smiles
   PUBCHEMIDS <- object %>% autonomics.import::flevels('PUBCHEM')   # NAS are never among levels
   SMILES <- PUBCHEMIDS %>%
            (function(x) split(x, ceiling(seq_along(x)/100))) %>%
             lapply(autonomics.annotate::pubchem_to_smiles) %>%
             data.table::rbindlist()

   # Merge into fdata
   autonomics.import::fdata(object) %<>% autonomics.support::left_join_keeping_rownames(SMILES, by = 'PUBCHEM')

   # Report
   autonomics.support::cmessage('\t\tAdd SMILES: %3d features -> %3d map to PUBCHEM -> %3d map to SMILES',
                                nrow(object),
                                sum(autonomics.import::fdata(object)$PUBCHEM!=''),
                                sum(!is.na(autonomics.import::fdata(object)$SMILES)))
   # Return
   object
}


