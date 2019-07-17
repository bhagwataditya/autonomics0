#' Feature identifer types
#' @export
FEATURE_ID_TYPES <- c(uniprot = 'UNIPROT', ensg = 'ENSEMBL')


build_FEATURE_IDENTIFIERS <- function(){

   FEATURE_IDENTIFIERS <-
      Map(function(feature_id_type){
         Map(function(cur_organism){
            organism <- ANNOTATED_ORGANISMS[[cur_organism]]
            orgdb_name <- sprintf('org.%s.%s.db', organism['abrev'], organism['source'])
            install_package_if_necessary(orgdb_name)
            orgdb <- load_orgdb(cur_organism)
            keys <- if(feature_id_type %in% AnnotationDbi::keytypes(orgdb)){
               AnnotationDbi::keys(orgdb, keytype = feature_id_type)
            } else {
               character(0)
            }
            list(db_pkg         = orgdb_name,
                 db_pkg_version = orgdb_name %>% utils::packageVersion() %>% as.character(),
                 keys           = keys)
         }, names(ANNOTATED_ORGANISMS))
      }, FEATURE_ID_TYPES) %>%
         set_names(names(FEATURE_ID_TYPES))

   usethis::use_data(FEATURE_IDENTIFIERS, pkg = ".", internal = FALSE, compress = "xz", overwrite = TRUE)

   return(invisible(TRUE))
}


utils::globalVariables('FEATURE_IDENTIFIERS')
