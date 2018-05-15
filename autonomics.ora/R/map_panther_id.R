#' Map panther class_id to class_term
#' @param class_id panther class id
#' @return vector with uniprots
#' @examples 
#' map_panther_id_to_term('PC00193')
#' map_panther_id_to_term('PC00137')
#' @importFrom magrittr         %>%
#' @export
map_panther_id_to_term <- function(class_id){
   PANTHER.db::PANTHER.db %>% 
      PANTHER.db::select(keys     = class_id, 
                         columns  = "CLASS_TERM",
                         keytype = "CLASS_ID") %>% 
      magrittr::extract2('CLASS_TERM')
}

#' Map panther class_id to uniprots
#' @param class_id panther class id
#' @param organism organism name, formatted as 'Homo sapiens'. names(PANTHER_ORGANISMS) shows currently supported values.
#' @return vector with uniprots
#' @examples 
#' map_panther_id_to_uniprots(class_id = 'PC00193', organism = 'Homo sapiens')
#' map_panther_id_to_uniprots(class_id = 'PC00137', organism = 'Homo sapiens')
#' @importFrom magrittr     %>%
#' @export
map_panther_id_to_uniprots <- function(class_id, organism){
   assertive.sets::assert_is_subset(organism, names(PANTHER_ORGANISMS))
   
   DB <- PANTHER.db::PANTHER.db
   PANTHER.db::pthOrganisms(DB) <- PANTHER_ORGANISMS %>% magrittr::extract2(organism)
   
   DB %>% PANTHER.db::select(keys     = class_id, 
                             columns  = "UNIPROT",
                             keytype = "CLASS_ID") %>%
      magrittr::extract2('UNIPROT')
}