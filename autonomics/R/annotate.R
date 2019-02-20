#' Map miRNAs to target ENSGs
#' 
#' \code{\link[autonomics.annotate]{mir_to_ensg}}
#' @name mir_to_ensg
#' @importFrom autonomics.annotate   mir_to_ensg
#' @export  mir_to_ensg
NULL


#' Map miRNAs to gene symbols
#' 
#' \code{\link[autonomics.annotate]{mir_to_gsymbol}}
#' @name mir_to_gsymbol
#' @importFrom autonomics.annotate   mir_to_gsymbol
#' @export   mir_to_gsymbol
NULL


#===========================================================================
# KEGG
#===========================================================================

#' Map kegg entry values to kegg pathways
#' 
#' @param x           charactervector, factorvector, or SummarizedExperiment
#' @param entry_var   kegg entry fvar
#' @param pathway_var kegg pathway fvar
#' @param ...         provided to enable S3 dispatch
#' @return charactervector or SummarizedExperiment
#' @examples
#' # charactervector
#'    require(magrittr)
#'    x <- c("C07326", "C04742", "C18218", "C18218", NA_character_, NA_character_, "", "")
#'    x %>% autonomics::kegg_entry_to_pathways()
#'    
#' # SummarizedExperiment   
#' if (require(autonomics.data)){
#'    x <- 'extdata/glutaminase/glutaminase.xlsx'    %>%
#'          system.file(package = 'autonomics.data') %>%
#'          autonomics::read_metabolon()
#'    x %>% autonomics::kegg_entry_to_pathways()
#' }
#'    
#' @references Backend is \href{http://www.kegg.jp/kegg/rest/keggapi.html}{KEGG's REST API}
#' @export
kegg_entry_to_pathways <- function(x, ...){
   UseMethod('kegg_entry_to_pathways', x)
}


#' @rdname kegg_entry_to_pathways
#' @importFrom magrittr %>%
#' @export
kegg_entry_to_pathways.character <- function(x, ...){
   
   # Satisfy check
      Entry <- Pathway <- NULL

   # Map available values
      idx <- !is.na(x) & x!=''
      if (sum(idx)==0) return(character(length(x)))
      keggurl <- x[idx] %>% unique() %>% paste0(collapse = '+') %>% sprintf('http://rest.kegg.jp/link/pathway/%s', .)
      if (!RCurl::url.exists(keggurl)) return(character(length(x)))
      cachefile <- tempfile()
      utils::download.file(keggurl, cachefile, quiet = TRUE)
      if (readLines(cachefile, n=1)=='') return(character(length(x)))

   # Format and return
      data.table::fread(cachefile, header = FALSE, col.names = c("Entry", "Pathway"))            %>%
      magrittr::extract(, Entry   := Entry   %>% stringi::stri_replace_first_fixed('cpd:',  '')) %>%
      magrittr::extract(, Pathway := Pathway %>% stringi::stri_replace_first_fixed('path:', '')) %>%
      magrittr::extract(, list(Pathway = paste0(Pathway, collapse = ';')), by = 'Entry')         %>%
      merge(data.table::data.table(Entry = x), ., by = 'Entry', all.x = TRUE, sort = FALSE, )    %>%
      magrittr::extract2('Pathway')

}


#' @rdname kegg_entry_to_pathways
#' @importFrom magrittr %>%
#' @export
kegg_entry_to_pathways.factor <- function(x, ...){
   x %>% as.character() %>% kegg_entry_to_pathways.character()
}


#' @rdname kegg_entry_to_pathways
#' @importFrom magrittr %>%
#' @export
kegg_entry_to_pathways.SummarizedExperiment <- function(x, entry_var = 'KEGG', pathway_var = 'KEGGPATHWAYS', ...){

   # Add KEGG Pathways
   autonomics.import::fdata(x)[[pathway_var]] <- autonomics.import::fdata(x)[[entry_var]]     %>%
                                                 autonomics.import::extract_first_from_collapsed() %>%
                                                 autonomics::kegg_entry_to_pathways()
   # Report
   autonomics.support::cmessage('\t\tAdd KEGG Pathways: %3d features -> %3d map to KEGG IDS -> %3d map to KEGG Pathways',
                                nrow(x),
                                sum(!is.na(autonomics.import::fdata(x)[[entry_var  ]])),
                                sum(!is.na(autonomics.import::fdata(x)[[pathway_var]])))

   # Return
   x
}


#' Get KEGG compound data
#'
#' @param KEGGIDS KEGG compound ids
#' @param verbose logical
#' @return data.frame(compound.entry, compound.name, pathway.id, pathway.name)
#' @examples
#' \dontrun{
#'    require(magrittr)
#'    KEGGIDS <- c("C04100", "C14829", "C14775", "C00599", "C01035", "C14772")
#'    KEGGIDS %>% get_kegg_compound_data()
#' }
#' @importFrom magrittr  %>%
#' @noRd
get_kegg_compound_data <- function(KEGGIDS, verbose = FALSE){

   # Satisfy CHECK
   . <- NULL

   # Assert
   assertive.strings::assert_all_are_matching_regex(KEGGIDS, 'C[0-9]+')

   # GET KEGG DATA
   KEGG_GET <- KEGGIDS %>%
               unique() %>%
               paste0(collapse = '+') %>%
               paste0('http://rest.kegg.jp/get/', .) %>%
               xml2::read_html()

   # WRITE
   if (verbose) KEGG_GET %>% rvest::html_text() %>% writeLines()

   # PARSE
   KEGGRECORDS <- KEGG_GET %>% rvest::html_text() %>% strsplit('///\\\n') %>% unlist()
   ENTRIES <- KEGGRECORDS %>% stringi::stri_extract_first_regex('ENTRY[ ]+C[0-9]+[ ]+Compound') %>%
                              stringi::stri_replace_first_regex('ENTRY[ ]+', '') %>%
                              stringi::stri_replace_first_regex('[ ]+Compound', '')
   NAMES   <- KEGGRECORDS %>% stringi::stri_extract_first_regex('\\\nNAME.+;') %>%
                              stringi::stri_replace_first_regex('\\\nNAME[ ]+', '') %>%
                              stringi::stri_replace_first_regex(';', '')
   PATHWAYS <- KEGGRECORDS %>% stringi::stri_extract_first_regex('\\\nPATHWAY(.|\n)+ENZYME') %>%
                               stringi::stri_replace_first_regex('\\\nPATHWAY[ ]+', '') %>%
                               stringi::stri_replace_first_regex('\\\nENZYME', '') %>%
                               stringi::stri_replace_all_regex('\\\n[ ]+', ';')
   PATHWAY_IDS <- PATHWAYS %>% stringi::stri_extract_all_regex('map[0-9]+') %>%
                               stringi::stri_extract_all_regex('[0-9]+') %>%
                               vapply(paste0, character(1), collapse = ';') %>%
                               (function(x){x[x=="NA"] <- NA_character_; x})
   PATHWAY_NAMES <- PATHWAYS %>% stringi::stri_replace_all_regex('map[0-9]+', '') %>%
                                 trimws() %>%
                                 stringi::stri_replace_all_regex(';[ ]+', ';')

   # PACK INTO DF AND RETURN
   data.frame(compound.entry   = ENTRIES,
              compound.name    = NAMES,
              pathway.id       = PATHWAY_IDS,
              pathway.name     = PATHWAY_NAMES)
}


#===================================================================
# PUBCHEM
#===================================================================


#' Map a vector of pubchemids to (canonical) smiles
#' 
#' @param x            charactervector, factorvector, or SummarizedExperiment
#' @param pubchem_var  pubchem fvar
#' @param smiles_var   smiles fvar
#' @param ...          provided to enable S3 dispatch
#' @return charactervector or SummarizedExperiment
#' @examples
#' require(magrittr)
#' # charactervector
#'    x <- c(NA_character_, "11988421", "10236635", "5283147", "91477", NA_character_)
#'    x %>% autonomics::pubchem_to_smiles()
#'    
#' # SummarizedExperiment
#' if (require(autonomics.data)){
#'    x <- autonomics.data::glutaminase
#'    x %>% autonomics::pubchem_to_smiles()
#' }
#' @references Backend is \href{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial}{PUBCHEM's REST interface}
#' @export
pubchem_to_smiles <- function(x, ...){
   UseMethod('pubchem_to_smiles', x)
}


#' @rdname pubchem_to_smiles
#' @importFrom magrittr %>%
#' @export
pubchem_to_smiles.character <- function(x, ...){

   # Satisfy CHECK
   . <- NULL

   # Download pubchem smiles
   cachefile <- tempfile()
   resturl <- x %>% stats::na.exclude() %>% as.vector() %>% unique() %>% paste0(collapse = ',') %>%
              sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/CanonicalSMILES/CSV', .)
   utils::download.file(resturl, cachefile, quiet = TRUE)

   # Return
   data.table::data.table(CID = as.integer(x)) %>%
   merge(data.table::fread(cachefile), by = 'CID', all.x = TRUE, sort = FALSE) %>%
   magrittr::extract2('CanonicalSMILES') 
}

#' @rdname pubchem_to_smiles
#' @importFrom magrittr %>% 
#' @export
pubchem_to_smiles.factor <- function(x, ...){
   x %>% as.character() %>% pubchem_to_smiles.character()
}

#' @rdname pubchem_to_smiles
#' @importFrom magrittr %>% 
#' @export
pubchem_to_smiles.SummarizedExperiment <- function(x, pubchem_var = 'PUBCHEM', smiles_var = 'SMILES', ...){

   # Satify CHECK
   . <- NULL

   # Assert
   assertive.sets::assert_is_subset(pubchem_var, autonomics.import::fvars(x))

   # Map to smiles
   PUBCHEMIDS <- x %>% 
                 autonomics.import::fvalues(pubchem_var) %>% 
                 autonomics.import::extract_first_from_collapsed(';') 
   SMILES <- PUBCHEMIDS %>%
            (function(x) split(x, ceiling(seq_along(x)/100))) %>%
             lapply(pubchem_to_smiles) %>% 
             unlist() %>% 
             unname()
   autonomics.import::fdata(x)[[smiles_var]] <- SMILES

   # Report
   autonomics.support::cmessage('\t\tAdd SMILES: %3d features -> %3d map to PUBCHEM -> %3d map to SMILES',
                                nrow(x),
                                sum(!is.na(autonomics.import::fdata(x)[[pubchem_var]])),
                                sum(!is.na(autonomics.import::fdata(x)[[smiles_var ]])))
   # Return
   x
}
