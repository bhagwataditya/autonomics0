#' Accession-based mapping/interconversion for multiple ontologies
#'
#' @aliases pathway_accessions_to_names gene_accessions_to_pathway_accessions
#' @param accessions \code{\link{character}} object representing InterPro accessions
#' ('IPR*'), for which names are to be retrieved
#' @param ontology single \code{\link{character}} defining the ontology used or triggering
#' ontology-deduction (\code{auto})
#' @param mapping_object a mapping object as returned by e.g. \code{\link{fetch_interpro_maps}};
#' if \code{NULL} map fetching is handled internally
#' @param ... further arguments handed on to map fetching
#' @param verbose single \code{\link{logical}} indicating whether retreival progress etc.
#' are printed
#' @return Named \code{\link{character}} of InterPro names (using the corresponding
#' accessions as names)
#' @examples
#' library(magrittr)
#' # InterPro: Simple call - determine ontology and fetch data on the fly
#' c("IPR038348", "IPR011771", "IPR031559", "IPR007277") %>%
#'   pathway_accessions_to_names()
#'
#' # InterPro: Same thing - Format ERROR
#' \dontrun{
#' c("IPR038348", "IPR011771", "IPR031559", "IPR007277", "MyFavoriteAccession") %>%
#'   pathway_accessions_to_names()}
#'
#' # InterPro: Formatting trickery to demonstrate empty string returns
#' c("IPR038348", "IPR011771", "IPRMyFavoriteAccession", "IPR007277") %>%
#'   pathway_accessions_to_names()
#'
#' # InterPro: Pre-fetch map (reduced networking)
#' mo_ip <- fetch_interpro_maps()
#' c("IPR038348", "IPR011771", "IPR031559", "IPR007277") %>%
#'   pathway_accessions_to_names(mapping_object = mo_ip)
#'
#' # A KEGG workflow - fetching here is gene id dependent
#' ## A. Prefetch data for a bunch of gene/protein accessions
#' mo_kegg <- c("dre:793465", "dre:571876", "dre:798771", "dre:564511") %>%
#'   fetch_kegg_maps()
#' str(mo_kegg)
#' ## B. Retrieve pathways for a specific gene/protein accession
#' (pwa <- gene_accessions_to_pathway_accessions('dre:798771', mapping_object = mo_kegg))
#' ## C. Get names for retrieved pathways
#' pwa %>%
#'   strsplit(split = ';') %>%
#'   unlist(use.names = FALSE) %>%
#'   pathway_accessions_to_names(mapping_object = mo_kegg)
#' @author Johannes Graumann
#' @importFrom magrittr %>%
#' @rdname ontology_accession_mapping
#' @export
pathway_accessions_to_names <- function(
  accessions,
  ontology       = c('auto', 'interpro', 'kegg'),
  mapping_object = NULL,
  ...,
  verbose        = FALSE)
{
# Check prerequisites -----------------------------------------------------
  accessions %>%
    assertive.types::assert_is_character()

  ontology %<>%
    match.arg(
      choices = c('auto', 'interpro', 'kegg'),
      several.ok = FALSE)
  if(ontology == 'auto')
  {
    test_string <- accessions %>%
      magrittr::extract2(1)
    if(stringi::stri_startswith_fixed(test_string,'IPR'))
    {
      ontology <- 'interpro'
    } else if(stringi::stri_detect_regex(test_string,'^\\w{3,3}\\d+'))
    {
      ontology <- 'kegg'
    } else {
      stop('Ontology determination failed!')
    }
  }

  if(ontology == 'interpro')
  {
    accessions %>%
      assertive.strings::assert_all_are_matching_regex('^IPR')
  } else if(ontology == 'kegg')
  {
    accessions %>%
      assertive.strings::assert_all_are_matching_regex('^\\w{3,3}\\d+')
  }

  if(!is.null(mapping_object))
  {
    mapping_object %>%
      assertive.types::assert_is_all_of(c('ontology_map','list'))
    if(ontology == 'interpro')
    {
      mapping_object %>%
        assertive.types::assert_is_all_of('interpro_map')
    } else if(ontology == 'kegg')
    {
      mapping_object %>%
        assertive.types::assert_is_all_of('kegg_map')
    }
  }

  verbose %>%
    assertive.types::assert_is_a_bool()

# Process -----------------------------------------------------------------
  if(is.null(mapping_object))
  {
    if(ontology == 'interpro')
    {
      mapping_object <- fetch_interpro_maps(..., verbose = verbose)
      if (is.null(mapping_object)) return(NULL)
    } else if(ontology == 'kegg')
    {
      mapping_object <- fetch_kegg_maps(...)
    }
  }

  # Generate lookup table
  pm <- mapping_object %>%
    magrittr::extract2('pathway_accession_to_name')

  # QC
  if(accessions %>%
    assertive.sets::is_subset(
      pm %>%
        names()) %>%
    magrittr::not() &&
    verbose)
  {
    warning('Encountered InterPro accessions NOT present in the database.')
  }

  pm %>%
    magrittr::extract(accessions) %>%
    (function(x){ifelse(is.na(x),"", x)}) %>%
    magrittr::set_names(accessions) %>%
    return()
}

#' @rdname ontology_accession_mapping
#' @export
gene_accessions_to_pathway_accessions <- function(
  accessions,
  ontology       = c('auto', 'kegg'),
  mapping_object = NULL,
  ...,
  verbose        = FALSE)
{
# Check prerequisites -----------------------------------------------------
  accessions %>%
    assertive.types::assert_is_character()

  ontology %<>%
    match.arg(
      choices = c('auto', 'kegg'),
      several.ok = FALSE)
  if(ontology == 'auto')
  {
    test_string <- accessions %>%
      magrittr::extract2(1)
    if(stringi::stri_detect_regex(test_string,'^\\w{3,3}:'))
    {
      ontology <- 'kegg'
    } else {
      stop('Ontology determination failed or ontology not implemented!')
    }
  }

  if(ontology == 'kegg')
   {
     accessions %>%
       assertive.strings::assert_all_are_matching_regex('^\\w{3,3}:')
   }

   if(!is.null(mapping_object))
   {
     mapping_object %>%
       assertive.types::assert_is_all_of(c('ontology_map','list'))
     if(ontology == 'kegg')
     {
       mapping_object %>%
         assertive.types::assert_is_all_of('kegg_map')
     }
   }

# Process -----------------------------------------------------------------
   if(is.null(mapping_object))
   {
     if(ontology == 'kegg')
     {
       mapping_object <- fetch_kegg_maps(kegg_geneids = accessions, ...)
     }
   }

   # Generate lookup table
   pm <- mapping_object %>%
     magrittr::extract2('gene_accession_to_pathway_accession')
   # lut_mapping_object <-  pm %>%
     # magrittr::extract2('NAME') %>%
     # magrittr::set_names(
       # pm %>%
         # magrittr::extract2('ACCESSION'))

   # QC
   if(accessions %>%
      assertive.sets::is_subset(
        pm %>%
        names()) %>%
      magrittr::not() &&
      verbose)
   {
     warning('Encountered InterPro accessions NOT present in the database.')
   }

  pm %>%
     magrittr::extract(accessions) %>%
     (function(x){ifelse(is.na(x),"", x)}) %>%
     magrittr::set_names(accessions) %>%
     return()
}
utils::globalVariables(c('.', 'accessions', 'interpro_url', 'verbose'))

