#' Infer organism from uniprot accessions
#' @param object eset
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    object %>% infer_organism_from_uniprot_accessions()
#'    autonomics.data::ALL %>% infer_organism_from_uniprot_accessions()
#' }
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::protein.ratios
#'    object %>% infer_organism_from_uniprot_accessions()
#' }
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma
#'    object %>% infer_organism_from_uniprot_accessions()
#' }
#' @importFrom magrittr   %>%
#' @export
infer_organism_from_uniprot_accessions <- function(object){

   sep <- autonomics.import::sep(object)
   eset_uniprots <- autonomics.import::fdata(object) %>%
                    magrittr::extract2(autonomics.import::uniprot_var(object)) %>%
                    as.character() %>%
                    (function(x) if (is.null(sep)) x else strsplit(x, sep)) %>%
                    unlist() %>%
                    unique()

    overlaps <- ANNOTATED_ORGANISMS %>%
                names() %>%
                lapply(function(x){
                         FEATURE_IDENTIFIERS %>%
                         magrittr::extract2("UNIPROT") %>%
                         magrittr::extract2(x) %>%
                         magrittr::extract2("keys") %>%
                         return()
                }) %>%
                vapply(function(x){
                          100*length(intersect(x, eset_uniprots))/length(eset_uniprots)
                       },
                       numeric(1))

    if (length(overlaps) == 0) return(NULL)

    organism_name <- names(ANNOTATED_ORGANISMS)[which.max(overlaps)]
    organism <- ANNOTATED_ORGANISMS[[organism_name]]
    organism_db <- sprintf('org.%s.%s.db', organism['abrev'], organism['source'])
    organism_db %>%
      install_package_if_necessary()
    if(FEATURE_IDENTIFIERS %>%
       magrittr::extract2("UNIPROT") %>%
       magrittr::extract2(organism_name) %>%
       magrittr::extract2("db_pkg_version") %>%
       magrittr::is_less_than(
          organism_db %>%
          utils::packageVersion()))
    {
      message("Dataset for inference of organism outdated.")
    }
    organism_name %>%
      return()
}

#' Infer organism from feature annotations
#' @param  object   eSet
#' @return A string corresponding to a name from \code{\link{ANNOTATED_ORGANISMS}}
#' @examples
#'    library(magrittr)
#'    if (require(autonomics.data)){
#'      autonomics.data::billing2016 %>% infer_organism()
#'    }
#'    if (require(atkin.2014)){
#'       atkin.2014::soma %>% infer_organism()
#'    }
#'    if (require(billing.differentiation.data)){
#'       require(magrittr)
#'       billing.differentiation.data::protein.ratios %>% infer_organism()
#'    }
#' @export
infer_organism <- function(object){
   if (autonomics.import::is_maxquant_eset(object) | autonomics.import::is_soma_eset(object)){
      infer_organism_from_uniprot_accessions(object)
   } else {
      stop('Implement autonomics.ora::infer_organism for other types of esets')
   }
}
# if (length(overlaps) == 0){
#   paste0(names(organisms), collapse = ', ') %>%
#     autonomics.support::cmessage('ora aborted - uniprot identifiers not from supported organisms (%s)', .)
#   return(NULL)
# }
#
# autonomics.support::cmessage('\t\tuniprot identifiers match %s genome', overlaps)

# Create a dataset for inferring organism
build_FEATURE_IDENTIFIERS <- function()
{
   key_types <- c('UNIPROT')
   FEATURE_IDENTIFIERS <- key_types %>%
      lapply(
         function(x){
            ANNOTATED_ORGANISMS                        %>%
              names()                                  %>%
              lapply(
                 function(y){
                    organism <- ANNOTATED_ORGANISMS[[y]]
                    org_db <- sprintf('org.%s.%s.db', organism['abrev'], organism['source'])
                    install_package_if_necessary(org_db)
                    list(
                       db_pkg = org_db,
                       db_pkg_version = org_db %>%
                          utils::packageVersion() %>%
                          as.character(),
                       keys = y %>%
                          load_org.xx.xx.db() %>%
                          AnnotationDbi::keys(keytype = x))
                 }) %>%
              magrittr::set_names(
                 ANNOTATED_ORGANISMS %>%
                   names())
            }) %>%
      magrittr::set_names(key_types)
   devtools::use_data(
      FEATURE_IDENTIFIERS,
      pkg = ".",
      internal = TRUE,
      compress = "xz",
      overwrite = TRUE)
   TRUE %>%
     invisible()
}
utils::globalVariables('FEATURE_IDENTIFIERS')
