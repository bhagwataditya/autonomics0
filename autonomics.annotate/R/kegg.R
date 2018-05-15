#' Get KEGG compound data
#'
#' Uses \href{http://www.kegg.jp/kegg/rest/keggapi.html}{KEGG's REST API}
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



