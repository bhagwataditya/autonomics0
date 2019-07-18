
#' Fetch keggpathway mapps
#'
#' kegggeneid -> keggpathwayid -> keggpathwayname
#'
#' @param kegggeneids \code{\link{character}} object representing KEGG gene/protein IDs,
#'                     for which information is retrieved via \code{\link[KEGGREST]{keggGet}} from package
#'                     \pkg{KEGGREST}
#' @param n_request_packet single \code{\link{numeric}} representing the number of accessions
#'                         sent to KEGG in a single request (web service restricted to 10)
#' @param keep_raw         single \code{\link{logical}} indicating whether the complete
#'                         entire data returned by KEGG should be retained or
#'                         solely preprocessed products thereoff
#' @return list (keggpathwayid_to_keggpathwayname, keggeneid_to_keggpathwayid)
#' @examples
#' kegggeneids <- c("dre:793465", "dre:571876", "dre:798771", "dre:564511")
#' fetch_keggpathway_maps(kegggeneids)
#' @author Johannes Graumann
#' @export
fetch_keggpathway_maps <- function(
    kegggeneids,
    n_request_packet = 10,
    keep_raw         = FALSE
){
    # Check prerequisites -----------------------------------------------------
    kegggeneids  %<>%  assertive.types::assert_is_character() %>%
                        unique()
    n_request_packet %>% assertive.types::assert_is_a_number() %>%
                         assertive.numbers::assert_all_are_whole_numbers() %>%
                         assertive.numbers::assert_all_are_in_left_open_range(
                            lower = 0, upper = 10)
    assertive.types::assert_is_a_bool(keep_raw)
    . <- NULL

    # Processing --------------------------------------------------------------
    # Split into packets (10 requests are KEGG REST limitation)
    list_kegggeneids <- kegggeneids %>%
                        split( ceiling(seq_along(.) / n_request_packet))

    # Retreive and format
    retrieved_list <-   list_kegggeneids %>%
                        lapply( function(x){
                            x %>% KEGGREST::keggGet() %>% set_names(x)
                        }) %>%
                        unlist(recursive = FALSE) %>%
                        set_names(stringi::stri_replace_first_regex(
                                    names(.), '^\\d+\\.(.*)', '$1'))

    # Create patwayid/description map
    keggpathwayid_to_keggpathwayname <- retrieved_list   %>%
        lapply(extract2, 'PATHWAY')     %>%
        lapply(as.list)  %>%  unname()  %>%  unlist(recursive = FALSE)  %>%
        extract(!duplicated(names(.)))  %>%  unlist()

    # Create geneid/patwayid map
    keggeneid_to_keggpathwayid <- retrieved_list  %>%
        lapply(extract2, 'PATHWAY')     %>%  lapply(names)  %>%
        lapply(paste0, collapse = ';')  %>%  unlist()

    # Assemble output & return
    output <- list(
        keggeneid_to_keggpathwayid         = keggeneid_to_keggpathwayid,
        keggpathwayid_to_keggpathwayname   = keggpathwayid_to_keggpathwayname)
    if(keep_raw) output[['raw_data']] <- retrieved_list
    class(output) <- c('kegg_map', 'ontology_map', class(output))
    return(output)
}



#' @rdname keggpathwayids_to_keggpathwaynames
#' @export
kegggeneids_to_keggpathwayids <- function(
    kegggeneids,
    maps = fetch_keggpathway_maps(kegggeneids),
    verbose        = FALSE
){
    # Assert ------------------------------------------------------------------
    kegggeneids %>% assertive.strings::assert_all_are_matching_regex('^\\w{3,3}:')
    maps %>% assertive.types::assert_is_all_of('kegg_map')

    # Generate lookup table ---------------------------------------------------
    pm <- maps$keggeneid_to_keggpathwayid

    # QC ----------------------------------------------------------------------
    if(!assertive.sets::is_subset(kegggeneids, names(pm)) & verbose){
        warning('kegggeneid values NOT in database.')
    }

    pm[kegggeneids]                        %>%
    (function(y){ifelse(is.na(y),"", y)})  %>%
    set_names(kegggeneids)
}




#' Map between kegggeneid, keggpathwayid and keggpathwayname
#'
#' @param kegggeneids    \code{\link{character}} with kegggeneids
#' @param keggpathwayids \code{\link{character}} with keggpathwayids
#' @param maps object returned by \code{\link{fetch_keggpathway_maps}};
#' @param verbose        \code{\link{logical}(1)}
#' @return Named \code{\link{character}}: name = keggpathwayid, value = keggpathwayname
#' @examples
#' library(magrittr)
#'
#' # Fetch keggpathway maps
#'     kegggeneids <- c("dre:793465", "dre:571876", "dre:798771")
#'     kegg_map <- kegggeneids %>% fetch_keggpathway_maps()
#'     kegg_map %>% lapply(head, 1)
#'
#' # kegggeneid -> keggpathwayid
#'     keggpathwayids <- kegggeneids_to_keggpathwayids(kegggeneids, kegg_map)
#'     keggpathwayids
#'
#' ## keggpathwayid -> keggpathwayname
#'     keggpathwayids %>%
#'     strsplit(split = ';')     %>%
#'     unlist(use.names = FALSE) %>%
#'     keggpathwayids_to_keggpathwaynames(kegg_map)
#' @author Johannes Graumann
#' @export
keggpathwayids_to_keggpathwaynames <- function(
    keggpathwayids,
    maps = fetch_keggpathway_maps(keggpathwayids),
    verbose = FALSE
){

    # Assert-------------------------------------------------------------------
    keggpathwayids %>% assertive.types::assert_is_character()
    keggpathwayids %>% assertive.strings::assert_all_are_matching_regex('^\\w{3,3}\\d+')
    maps %>% assertive.types::assert_is_all_of('kegg_map')
    assertive.types::assert_is_a_bool(verbose)

    # Execute-------------------------s- --------------------------------------
    pm <- maps %>% extract2('keggpathwayid_to_keggpathwayname')
    if(!assertive.sets::is_subset(keggpathwayids, names(pm)) && verbose){
        warning('keggpathway ids NOT in database.')
    }

    pm[keggpathwayids]                                  %>%
    (function(y){ifelse(is.na(y),"", y)})  %>%
    set_names(keggpathwayids)

}
