#' Map a vector of PUBCHEM IDS to canonical smiles
#' @param PUBCHEMIDS character vector of PUBCHEM IDS
#' @return character vector with canonical smiles
#' @examples
#' require(magrittr)
#' PUBCHEMIDS <- c("11988421", "10236635", "5283147",  "91477",    "500",      "5283142")
#' PUBCHEMIDS %>% autonomics.annotate::pubchem_to_smiles()
#' @importFrom magrittr %>%
#' @export
pubchem_to_smiles <- function(PUBCHEMIDS){

   # Satisfy CHECK
   . <- NULL

   # Download pubchem smiles
   cachefile <- tempfile()
   resturl <- PUBCHEMIDS %>%
      paste0(collapse = ',') %>%
      sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/CanonicalSMILES/CSV', .)
   utils::download.file(resturl, cachefile, quiet = TRUE)

   # Return
   data.table::fread(cachefile) %>% magrittr::set_names(c('PUBCHEM', 'SMILES'))
}

