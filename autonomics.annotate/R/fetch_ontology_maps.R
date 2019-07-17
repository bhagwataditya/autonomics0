#' Functions to generate maps for various ontology accession-based mapping operations
#'
#' @aliases fetch_kegg_maps fetch_interpro_maps
#' @param interpro_url single \code{\link{character}} representing the URL to the current
#' \code{names.dat} file from InterPro (default is \code{\link{paste0}}d to avoid line
#' length issues in the package documentation)
#' @param verbose      \code{\link{logical}} rendering the function silent (or not)
#' @param kegg_geneids \code{\link{character}} object representing KEGG gene/protein IDs,
#' for which information is retrieved via \code{\link[KEGGREST]{keggGet}} from package
#' \pkg{KEGGREST}
#' @param n_request_packet single \code{\link{numeric}} representing the number of accessions
#' send to KEGG in a single request (web service restricted to 10)
#' @param keep_raw         single \code{\link{logical}} indicating whether the complete
#' entire data returned by KEGG should be retained or solely preprocessed products thereoff
#' @examples
#' library(magrittr)
#' # Fetch InterPro map
#' fetch_interpro_maps() %>%
#'   str()
#'
#' # Fetch small KEGG map (gene/protein id specific, as whole data set not public)
#' c("dre:793465", "dre:571876", "dre:798771", "dre:564511") %>%
#'   fetch_kegg_maps() %>%
#'   str()
#' @author Johannes Graumann
#' @rdname fetch_ontology_maps
#' @export
fetch_interpro_maps <- function(
  interpro_url   = paste0(
     c('http:/',
     'ftp.ebi.ac.uk',
      'pub',
      'databases',
      'interpro',
      'current',
      'names.dat'), collapse = '/'),
  verbose        = FALSE)
{
  # Check prerequisites -----------------------------------------------------
  interpro_url %>%
    assertive.types::is_a_string()

  verbose %>%
    assertive.types::assert_is_a_bool()

  if (!RCurl::url.exists(interpro_url, timeout = 5)){
     autonomics.support::cmessage('Interpro FTP site not working - returning NULL\n%s', interpro_url)
     return(NULL)
  }

  # Processing --------------------------------------------------------------
  # Retrieve the matching table
  mapping_object <- interpro_url %>%
    data.table::fread(
      data.table = FALSE,
      header     = FALSE,
      verbose    = verbose,
      showProgress = verbose) %>%
    set_names(c('ACCESSION', 'NAME'))
  vectorized_mapping_object <- mapping_object %>%
    extract2('NAME') %>%
    set_names(
      mapping_object %>%
        extract2('ACCESSION'))

  output <- list(pathway_accession_to_name = vectorized_mapping_object)
  class(output) <- c('interpro_map', 'ontology_map', class(output))
  output %>%
    return()
}
#' @rdname fetch_ontology_maps
#' @export
fetch_kegg_maps <- function(
  kegg_geneids,
  n_request_packet = 10,
  keep_raw         = FALSE)
{
  # Check prerequisites -----------------------------------------------------
  kegg_geneids %<>%
    assertive.types::assert_is_character() %>%
    unique()

  n_request_packet %>%
    assertive.types::assert_is_a_number() %>%
    assertive.numbers::assert_all_are_whole_numbers() %>%
    assertive.numbers::assert_all_are_in_left_open_range(
      lower = 0,
      upper = 10)

  keep_raw %>%
    assertive.types::assert_is_a_bool()

  # Processing --------------------------------------------------------------
  # Split into packets (10 requests are KEGG REST limitation)
  list_kegg_geneids <- kegg_geneids %>%
    split(
      seq_along(.) %>%
        divide_by(n_request_packet) %>%
        ceiling)

  # Retreive and format
  retrieved_list <- list_kegg_geneids %>%
    # Retreive the data
    lapply(
      function(x)
      {
        x %>%
          KEGGREST::keggGet() %>%
          set_names(x) %>%
          return()
      }) %>%
    unlist(recursive = FALSE) %>%
    set_names(
      stringi::stri_replace_first_regex(names(.), '^\\d+\\.(.*)', '$1'))

  # Create patwayid/description map
  pathway_accession_to_name <- retrieved_list %>%
    lapply(extract2, 'PATHWAY') %>%
    lapply(as.list) %>%
    unname() %>%
    unlist(recursive = FALSE) %>%
    extract(!duplicated(names(.))) %>%
    unlist()

  # Create geneid/patwayid map
  gene_accession_to_pathway_accession <- retrieved_list %>%
    lapply(extract2, 'PATHWAY') %>%
    lapply(names) %>%
    lapply(paste0, collapse = ';') %>%
    unlist()

  # Assemble output & return
  output <- list(
    gene_accession_to_pathway_accession = gene_accession_to_pathway_accession,
    pathway_accession_to_name           = pathway_accession_to_name)
  if(keep_raw)
  {
    output[['raw_data']] <- retrieved_list
  }
  class(output) <- c('kegg_map', 'ontology_map', class(output))
  output %>%
    return()
}

