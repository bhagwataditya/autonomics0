
#' Open uniprot connection
#' @param object SummarizedExperiment
#' @return uniprot webservice connection
#' @examples
#' \dontrun{
#'    require(magrittr)
#'    if (require(autonomics.data)){
#'       object <- autonomics.data::stemcomp.proteinratios
#'       object %>% autonomics.import::open_uniprot_connection()
#'    }
#' }
#' @importFrom magrittr %>%
#' @export
open_uniprot_connection <- function(object){
   object                                                %>%
      autonomics.import::uniprot_values(first_only = TRUE)  %>%
      magrittr::extract(1:10)                               %>%
      autonomics.annotate::connect_to_uniprot()
}

#' Annotate proteingroups through uniprot.ws
#' @param object      SummarizedExperiment
#' @param connection  uniprot webservice connection
#' @param columns     uniprot webservice columns
#' @return Annotated SummarizedExperiment
#' @examples
#' require(magrittr)
#' \dontrun{
#'    if (require(autonomics.data)){
#'       object <- autonomics.data::stemcomp.proteinratios
#'       connection <- autonomics.import::open_uniprot_connection(object)
#'       object[1:10, ] %>% annotate_proteingroups(connection, c('SUBCELLULAR-LOCATIONS', 'GO-ID'))
#'    }
#' }
#' @importFrom magrittr %>%  %<>%
#' @export
annotate_proteingroups <- function(
   object,
   connection = object %>% autonomics.import::open_uniprot_connection(),
   columns    = c('SUBCELLULAR-LOCATIONS', 'INTERPRO', 'GO-ID')
){
   # Assert
   assertive.base::assert_is_identical_to_true(class(object) == 'SummarizedExperiment')

   # Restrict to first uniprot accession
   autonomics.import::fdata(object)$`Uniprot accessions` %<>% stringi::stri_split_fixed(';') %>% vapply(extract, character(1), 1)

   # Fetch annotations from uniprot
   annotations <- autonomics.import::fdata(object)$`Uniprot accessions` %>%
      autonomics.annotate::annotate_uniprot_with_webservice(connection = connection, columns = columns)

   # Merge in annotations
   autonomics.import::fdata(object) %<>% merge(annotations, by.x = 'Uniprot accessions', by.y = 'UNIPROTKB', sort = FALSE)

   # Return
   return(object)

}

