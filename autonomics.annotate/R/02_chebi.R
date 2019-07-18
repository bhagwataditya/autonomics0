#' Get canonical smiles from CHEBI
#'
#' Web scraping based, so non-vectorizable and slow.
#' No REST interface available (as per Dec 2017), only Perl and Java SOAP API.
#' Metabolon files do not give cheEBI annotations.
#'
#' @param CHEBIID chEBI ID (numeric, character, or factor)
#' @return canonical smiles string
#' @examples
#'    chebi_to_smiles('28733')
#' @export
chebi_to_smiles <- function(CHEBIID){

    # Satisfy CHECK
    . <- NULL

    # Assert
    assertive.properties::assert_is_scalar(CHEBIID)

    # Format
    CHEBIID %<>% as.character()

    # Map
    CHEBIID                                                               %>%
    sprintf('http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:%s', .) %>%
    xml2::read_html()                                                     %>%
    rvest::html_text()                                                    %>%
    stringi::stri_extract_first_regex('SMILES\\\n.+\\n')                  %>%
    stringi::stri_replace_first_regex('SMILES\\\n[ ]+', '')
}
