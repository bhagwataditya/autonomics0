
#' Fetch interpro maps: interproid -> interproname
#'
#' @param verbose \code{\link{logical}} rendering the function silent (or not)
#' @return list (interproids_to_interpronames)
#' @examples
#' fetch_interpro_maps()
#' @author Johannes Graumann
#' @export
fetch_interpro_maps <- function(verbose  = FALSE){

    # Check prerequisites -----------------------------------------------------
    url <- 'http://ftp.ebi.ac.uk/pub/databases/interpro/current/names.dat'
    assertive.types::is_a_string(url)
    assertive.types::assert_is_a_bool(verbose)
    if (!RCurl::url.exists(url, timeout = 5)){
        message('Interpro FTP site not working - returning NULL\n', url)
        return(NULL)
    }

    # Processing --------------------------------------------------------------
    # Retrieve the matching table
    mapping_object <-
        data.table::fread(
            url,
            data.table   = FALSE,
            header       = FALSE,
            verbose      = verbose,
            showProgress = verbose) %>%
            set_names(c('ACCESSION', 'NAME'))
    vectorized_mapping_object <-    mapping_object$NAME %>%
                                    set_names(mapping_object$ACCESSION)

    output <- list(interproids_to_interpronames = vectorized_mapping_object)
    class(output) <- c('interpro_map', 'ontology_map', class(output))
    return(output)
}



#' Convert interproid to interproname
#'
#' @param interproids \code{\link{character}} with InterPro ids ('IPR*')
#' @param mapping_object mapping object as returned by \code{\link{fetch_interpro_maps}};
#' @param verbose        \code{\link{logical}(1)}
#' @return Named \code{\link{character}}: name = interproid, value = interproname
#' @examples
#' library(magrittr)
#'
#' # Fetch interpromap
#'     interpromap  <- fetch_interpro_maps()
#'
#' # Interproid -> Interproname
#'     interproids <-  c("IPR038348", "IPR011771", "IPRNotInDatabase")
#'     interproids %>% interproids_to_interpronames(interpromap)
#'     interproids %>% interproids_to_interpronames()
#'
#' @author Johannes Graumann
#' @export
interproids_to_interpronames <- function(
    interproids,
    mapping_object = fetch_interpro_maps(),
    verbose = FALSE
){
    # Check -------------------------------------------------------------------
    if (is.null(mapping_object)) return(NULL)
    interproids %>% assertive.types::assert_is_character()
    interproids %>% assertive.strings::assert_all_are_matching_regex('^IPR')
    mapping_object %>% assertive.types::assert_is_all_of('interpro_map')
    assertive.types::assert_is_a_bool(verbose)


    # Fetch maps and create lookup table---------------------------------------
    pm <- mapping_object$interproids_to_interpronames
    if(!assertive.sets::is_subset(interproids, names(pm)) && verbose){
        warning('InterPro ids NOT in database.')
    }

    # Return-------------------------------------------------------------------
    pm[interproids]                         %>%
    (function(y){ifelse(is.na(y),"", y)})   %>%
    set_names(interproids)
}

