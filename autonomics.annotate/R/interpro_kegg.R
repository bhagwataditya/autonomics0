#' Functions to generate maps for various ontology accession-based mapping operations
#'
#' @aliases fetch_kegg_maps fetch_interpro_maps
#' @param url          single \code{\link{character}} representing the URL to the current
#'                     \code{names.dat} file from InterPro (default is \code{\link{paste0}}d to avoid line
#'                     length issues in the package documentation)
#' @param verbose      \code{\link{logical}} rendering the function silent (or not)
#' @param kegg_geneids \code{\link{character}} object representing KEGG gene/protein IDs,
#'                     for which information is retrieved via \code{\link[KEGGREST]{keggGet}} from package
#'                     \pkg{KEGGREST}
#' @param n_request_packet single \code{\link{numeric}} representing the number of accessions
#'                         sent to KEGG in a single request (web service restricted to 10)
#' @param keep_raw         single \code{\link{logical}} indicating whether the complete
#'                         entire data returned by KEGG should be retained or
#'                         solely preprocessed products thereoff
#' @examples
#' library(magrittr)
#'
#' # Fetch InterPro map
#'     fetch_interpro_maps() %>% str()
#'
#' # Fetch small KEGG map (gene/protein id specific, as whole data set not public)
#'     c("dre:793465", "dre:571876", "dre:798771", "dre:564511") %>%
#'     fetch_kegg_maps() %>%
#'     str()
#' @author Johannes Graumann
#' @rdname fetch_ontology_maps
#' @export
fetch_interpro_maps <- function(
    url      = 'http:/ftp.ebi.ac.uk/pub/databases/interpro/current/names.dat',
    verbose  = FALSE)
{
    # Check prerequisites -----------------------------------------------------
    assertive.types::is_a_string(url)
    assertive.types::assert_is_a_bool(verbose)
    if (!RCurl::url.exists(url, timeout = 5)){
        cmessage('Interpro FTP site not working - returning NULL\n', url)
        return(NULL)
    }

    # Processing --------------------------------------------------------------
    # Retrieve the matching table
    mapping_object <- data.table::fread(
                        url,
                        data.table = FALSE,
                        header     = FALSE,
                        verbose    = verbose,
                        showProgress = verbose) %>%
                        set_names(c('ACCESSION', 'NAME'))
    vectorized_mapping_object <-
        mapping_object %>%
        extract2('NAME') %>%
        set_names(mapping_object %>% extract2('ACCESSION'))

    output <- list(pathway_accession_to_name = vectorized_mapping_object)
    class(output) <- c('interpro_map', 'ontology_map', class(output))
    return(output)
}


#' @rdname fetch_ontology_maps
#' @export
fetch_kegg_maps <- function(
    kegg_geneids,
    n_request_packet = 10,
    keep_raw         = FALSE
){
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
    list_kegg_geneids <-    kegg_geneids %>%
                            split(  seq_along(.) %>%
                                    divide_by(n_request_packet) %>%
                                    ceiling)

    # Retreive and format
    retrieved_list <-
        list_kegg_geneids %>%
        # Retreive the data
        lapply(function(x){
                    x                    %>%
                    KEGGREST::keggGet()  %>%
                    set_names(x)         %>%
                    return()}
        ) %>%
        unlist(recursive = FALSE) %>%
        set_names(
            stringi::stri_replace_first_regex(names(.), '^\\d+\\.(.*)', '$1'))

    # Create patwayid/description map
    pathway_accession_to_name <-
        retrieved_list                 %>%
        lapply(extract2, 'PATHWAY')    %>%
        lapply(as.list)                %>%
        unname()                       %>%
        unlist(recursive = FALSE)      %>%
        extract(!duplicated(names(.))) %>%
        unlist()

    # Create geneid/patwayid map
    gene_accession_to_pathway_accession <-
        retrieved_list                 %>%
        lapply(extract2, 'PATHWAY')    %>%
        lapply(names)                  %>%
        lapply(paste0, collapse = ';') %>%
        unlist()

    # Assemble output & return
    output <- list(
        gene_accession_to_pathway_accession = gene_accession_to_pathway_accession,
        pathway_accession_to_name           = pathway_accession_to_name)
    if(keep_raw) output[['raw_data']] <- retrieved_list
    class(output) <- c('kegg_map', 'ontology_map', class(output))
    return(output)
}



#' Accession-based mapping/interconversion for multiple ontologies
#'
#' @aliases pathway_accessions_to_names gene_accessions_to_pathway_accessions
#' @param accessions \code{\link{character}} object representing InterPro accessions
#'                   ('IPR*'), for which names are to be retrieved
#' @param ontology single \code{\link{character}} defining the ontology used or triggering
#'                 ontology-deduction (\code{auto})
#' @param mapping_object a mapping object as returned by e.g. \code{\link{fetch_interpro_maps}};
#'                       if \code{NULL} map fetching is handled internally
#' @param ...            further arguments handed on to map fetching
#' @param verbose        single \code{\link{logical}} indicating whether retreival progress etc.
#'                       are printed
#' @return Named \code{\link{character}} of InterPro names
#'        (using corresponding accessions as names)
#' @examples
#' library(magrittr)
#'
#' # InterPro: Simple call - determine ontology and fetch data on the fly
#'     c("IPR038348", "IPR011771", "IPR031559") %>%
#'       pathway_accessions_to_names()
#'
#' # InterPro: Same thing - Format ERROR
#'     \dontrun{
#'         c("IPR038348", "IPR011771", "IPR031559", "MyFavoriteAccession") %>%
#'           pathway_accessions_to_names()
#'      }
#'
#' # InterPro: Formatting trickery to demonstrate empty string returns
#'     c("IPR038348", "IPR011771", "IPRMyFavoriteAccession", "IPR007277") %>%
#'     pathway_accessions_to_names()
#'
#' # InterPro: Pre-fetch map (reduced networking)
#'     mo_ip <- fetch_interpro_maps()
#'     c("IPR038348", "IPR011771", "IPR031559", "IPR007277") %>%
#'     pathway_accessions_to_names(mapping_object = mo_ip)
#'
#' # KEGG workflow - fetching here is gene id dependent
#' ## A. Prefetch data for a bunch of gene/protein accessions
#'     mo_kegg <-  c("dre:793465", "dre:571876", "dre:798771", "dre:564511") %>%
#'                 fetch_kegg_maps()
#'     str(mo_kegg)
#'
#' ## B. Retrieve pathways for a specific gene/protein accession
#'     (pwa <- gene_accessions_to_pathway_accessions(
#'                 'dre:798771', mapping_object = mo_kegg))
#' ## C. Get names for retrieved pathways
#'     pwa %>%
#'     strsplit(split = ';')     %>%
#'     unlist(use.names = FALSE) %>%
#'     pathway_accessions_to_names(mapping_object = mo_kegg)
#' @author Johannes Graumann
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
    assertive.types::assert_is_character(accessions)

    ontology %<>% match.arg(
                    choices = c('auto', 'interpro', 'kegg'),
                    several.ok = FALSE)

    if(ontology == 'auto'){
        test_string <- accessions %>% extract2(1)
        ontology <-
            if (stringi::stri_startswith_fixed(test_string,'IPR')){
                'interpro'
            } else if(stringi::stri_detect_regex(test_string,'^\\w{3,3}\\d+')){
                'kegg'
            } else {
                stop('Ontology determination failed!')
            }
    }

    if(ontology == 'interpro'){
        accessions %>%
        assertive.strings::assert_all_are_matching_regex('^IPR')
    } else if(ontology == 'kegg'){
        accessions %>%
        assertive.strings::assert_all_are_matching_regex('^\\w{3,3}\\d+')
    }

    if(!is.null(mapping_object)){
        mapping_object %>%
        assertive.types::assert_is_all_of(c('ontology_map','list'))
        map_type <- switch( ontology,
                            interpro = 'interpro_map',
                            kegg = 'kegg_map')
        mapping_object %>% assertive.types::assert_is_all_of(map_type)
    }

    assertive.types::assert_is_a_bool(verbose)

    # Process -----------------------------------------------------------------
    if(is.null(mapping_object)){
        mapping_object <-
        if(ontology == 'interpro'){
            mapping_object <- fetch_interpro_maps(..., verbose = verbose)
            if (is.null(mapping_object)) return(NULL)
        } else if(ontology == 'kegg'){
            mapping_object <- fetch_kegg_maps(...)
        }
    }

    # Generate lookup table
    pm <- mapping_object %>% extract2('pathway_accession_to_name')

    # QC
    if(!assertive.sets::is_subset(accessions, names(pm)) && verbose){
        warning('Encountered InterPro accessions NOT present in the database.')
    }

    pm %>%
    extract(accessions) %>%
    (function(x){ifelse(is.na(x),"", x)}) %>%
    set_names(accessions) %>%
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
    assertive.types::assert_is_character(accessions)

    ontology %<>% match.arg(choices = c('auto', 'kegg'), several.ok = FALSE)
    if(ontology == 'auto'){
        test_string <- accessions %>% extract2(1)
        ontology <-
            if(stringi::stri_detect_regex(test_string,'^\\w{3,3}:')){
                ontology <- 'kegg'
            } else {
                stop('Ontology determination failed or ontology not implemented!')
            }
    }

    if(ontology == 'kegg'){
        accessions %>%
        assertive.strings::assert_all_are_matching_regex('^\\w{3,3}:')
    }

    if(!is.null(mapping_object)){
        mapping_object %>%
            assertive.types::assert_is_all_of(c('ontology_map','list'))
        if(ontology == 'kegg'){
            mapping_object %>% assertive.types::assert_is_all_of('kegg_map')
        }
    }

    # Process -----------------------------------------------------------------
    if(is.null(mapping_object)){
        if(ontology == 'kegg'){
            mapping_object <- fetch_kegg_maps(kegg_geneids = accessions, ...)
        }
    }

    # Generate lookup table
    pm <- mapping_object %>% extract2('gene_accession_to_pathway_accession')
    # lut_mapping_object <-  pm %>%
    # extract2('NAME') %>%
    # set_names(
    # pm %>%
    # extract2('ACCESSION'))

    # QC
    if(!assertive.sets::is_subset(accessions, names(pm)) & verbose){
        warning('Encountered accessions NOT present in the database.')
    }

    pm %>%
    extract(accessions) %>%
    (function(x){ifelse(is.na(x),"", x)}) %>%
    set_names(accessions) %>%
    return()
}
utils::globalVariables(c('.', 'accessions', 'interpro_url', 'verbose'))

