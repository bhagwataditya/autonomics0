
#' Load orgdb object
#' @param organism value in ANNOTATED_ORGANISMS
#' @examples
#' load_orgdb(organism = 'Homo sapiens')
#' load_orgdb(organism = 'Xenopus laevis')
#' @export
load_orgdb <- function(organism){
   orgdb_name <- get_orgdb_name(organism)
   install_package_if_necessary(orgdb_name)
   utils::getFromNamespace(orgdb_name, orgdb_name)
}
