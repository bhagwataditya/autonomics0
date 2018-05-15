#' Load annotation map
#' @param organism value in ANNOTATED_ORGANISMS
#' @examples 
#' load_org.xx.xx.db(organism = 'Homo sapiens')
#' load_org.xx.xx.db(organism = 'Xenopus laevis')
#' @importFrom magrittr %>% 
#' @export
load_org.xx.xx.db <- function(organism){
  organism <- ANNOTATED_ORGANISMS[[organism]]
  orgdb <- sprintf('org.%s.%s.db', organism['abrev'], organism['source'])
  install_package_if_necessary(orgdb)
  utils::getFromNamespace(orgdb, orgdb)
}