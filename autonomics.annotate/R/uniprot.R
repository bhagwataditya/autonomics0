
#=========================

#' Annotate uniprot accessions
#' @param uniprot character vector with uniprot accessions
#' @param fastafile path to fasta file
#' @return data.table(uniprot, genename, proteinname, reviewed, existence)
#' @export
annotate_with_uniprot_fastafile <- function(uniprot, fastafile){

   # Assert
   assertive.types::assert_is_character(uniprot)
   assertive.files::assert_all_are_existing_files(fastafile)

   # Import
   `:=` <- data.table::`:=`

   # Load (relevant portion of) fasta
   fasta <- seqinr::read.fasta(fastafile)
   all_accessions <- fasta %>% names() %>% stringi::stri_split_fixed('|') %>% vapply(extract, character(1), 2)
   assertive.sets::assert_is_subset(uniprot, all_accessions)
   fasta %<>% magrittr::extract(all_accessions %in% uniprot)

   # Extract annotations
   dt <- data.table::data.table(
      uniprot    = fasta %>% names() %>% stringi::stri_split_fixed('|') %>% vapply(extract, character(1), 2),
      reviewed   = fasta %>% names() %>% stringi::stri_split_fixed('|') %>% vapply(extract, character(1), 1) %>% magrittr::equals('sp') %>% as.numeric(),
      #entryname  = fasta %>% names() %>% stringi::stri_split_fixed('|') %>% vapply(extract, character(1), 3),
      annotation = fasta %>% vapply(attr, character(1), 'Annot') %>% unname())

   # Sequence version - not of interest
   dt %>% magrittr::extract(, annotation := annotation %>% stringi::stri_replace_last_regex(' SV=[0-9]', ''))

   # Protein existence
   pattern <- ' PE=[0-9]'
   dt  %>% magrittr::extract(, existence  := annotation %>% stringi::stri_extract_last_regex(pattern) %>% substr(5,nchar(.)) %>% as.numeric())
   dt  %>% magrittr::extract(, annotation := annotation %>% stringi::stri_replace_last_regex(pattern, ''))

   # Gene names
   pattern <- ' GN=.+$'
   dt  %>% magrittr::extract(, genename   := annotation %>% stringi::stri_extract_last_regex(pattern) %>% substr(5,nchar(.)))
   dt  %>% magrittr::extract(, annotation := annotation %>% stringi::stri_replace_last_regex(pattern, ''))

   # Organism identifier
   pattern <- ' OX=[0-9]+'
   #dt %>% magrittr::extract(, orgid      := annotation %>% stringi::stri_extract_last_regex(pattern) %>% substr(5,nchar(.)))
   dt  %>% magrittr::extract(, annotation := annotation %>% stringi::stri_replace_last_regex(pattern, ''))

   # Organism name
   pattern <- ' OS=.+$'
   #dt %>% magrittr::extract(, orgname    := annotation %>% stringi::stri_extract_last_regex(pattern) %>% substr(5,nchar(.)))
   dt  %>% magrittr::extract(, annotation := annotation %>% stringi::stri_replace_last_regex(pattern, ''))

   # Protein names
   pattern <- ' .+$'
   dt  %>% magrittr::extract(, proteinname := annotation %>% stringi::stri_extract_last_regex(pattern) %>% substr(2,nchar(.)))
   dt  %>% magrittr::extract(, annotation := annotation %>% stringi::stri_replace_last_regex(pattern, ''))

   # Brushup
   dt  %>% magrittr::extract(, annotation := NULL)
   dt %<>% magrittr::extract(, c('uniprot', 'genename', 'proteinname', 'reviewed', 'existence'), with = FALSE)

   # Return
   return(dt)
}


#=======================================

#' Open a uniprot.ws connection
#'
#' Infers organism from uniprot ids and opens a
#' uniprot webservice connection.
#'
#' @param values uniprot ids
#' @return uniprot.ws connection
#' @examples
#' \dontrun{
#' require(magrittr)
#' values <- c("A0A024R4M0", "M0R210", "G3HSF3")
#' values %>% autonomics.annotate::connect_to_uniprot()
#' }
#' @importFrom magrittr %>%
#' @export
connect_to_uniprot <- function(values){
   . <- NULL
   values                                                          %>%
   autonomics.annotate::infer_organism('uniprot', verbose = FALSE) %>%
   magrittr::extract2(autonomics.annotate::ANNOTATED_ORGANISMS, .) %>%
   magrittr::extract2('taxonid') %>%
   as.numeric() %>%
   UniProt.ws::UniProt.ws()
}

#===============================================

#' Fetch uniprot annotations
#' @param values uniprot ids (character vector)
#' @param up uniprot.ws connection (output of autonomics.annotate::connect_to_uniprot)
#' @return dataframe with annotations
#' @examples
#' \dontrun{
#' require(magrittr)
#' values <- c("A0A024R4M0", "M0R210", "G3HSF3")
#' up <- values %>% autonomics.annotate::connect_to_uniprot()
#' values %>% autonomics.annotate::fetch_uniprot_annotations(up)
#' }
#' @importFrom magrittr %>%
#' @export
fetch_uniprot_annotations <- function(values, up){
   AnnotationDbi::select(
      up,
      keys = values,
      columns = c('SUBCELLULAR-LOCATIONS', 'KEGG', 'GO-ID', 'INTERPRO'))
}

#=====================================================

#' Clean uniprot "score" values
#' @param values character vector
#' @return character vectors
#' @examples
#' \dontrun{
#' require(magrittr)
#' values <- c("A0A024R4M0", "M0R210", "G3HSF3")
#' up <- values %>% autonomics.annotate::connect_to_uniprot()
#' annotations <- values %>% autonomics.annotate::fetch_uniprot_annotations(up)
#' annotations$SCORE %T>% print() %>% autonomics.annotate::clean_uniprot_score_values()
#' }
#' @importFrom magrittr %>%
#' @export
clean_uniprot_score_values <- function(values){
   values                              %>%
   substr(1, 1)                        %>%
  (function(x){x[is.na(x)] <- '0'; x}) %>%
   as.numeric()
}

#=============================================

#' Clean uniprot "protein existence" values
#' @param values uniprot protein existence values (character vector)
#' @return cleaned values (character vector)
#' @examples
#' \dontrun{
#' require(magrittr)
#' values <- c("A0A024R4M0", "M0R210", "G3HSF3")
#' up <- values %>% autonomics.annotate::connect_to_uniprot()
#' annotations <- values %>% autonomics.annotate::fetch_uniprot_annotations(up)
#' annotations$EXISTENCE %T>% print() %>% autonomics.annotate::clean_uniprot_existence_values()
#' }
#' @importFrom magrittr %>%
#' @export
clean_uniprot_existence_values <- function(values){
   values   %>%
   stringi::stri_replace_first_fixed("Evidence at protein level",    1)  %>%
   stringi::stri_replace_first_fixed("Evidence at transcript level", 2)  %>%
   stringi::stri_replace_first_fixed("Inferred from homology",       3)  %>%
   stringi::stri_replace_first_fixed("Predicted",                    4)  %>%
   stringi::stri_replace_first_fixed("Uncertain",                    5)  %>%
  (function(x){x[is.na(x)] <- '5'; x})                                   %>%
   as.numeric()
}

#===================================

#' Clean uniprot "reviewed" values
#' @param values uniprot "reviewed" values (character vector)
#' @return cleaned values (character vector)
#' @examples
#' \dontrun{
#' require(magrittr)
#' values <- c("A0A024R4M0", "M0R210", "G3HSF3")
#' up <- values %>% autonomics.annotate::connect_to_uniprot()
#' annotations <- values %>% autonomics.annotate::fetch_uniprot_annotations(up)
#' annotations$REVIEWED %T>% print() %>% autonomics.annotate::clean_uniprot_reviewed_values()
#' }
#' @importFrom magrittr %>%
#' @export
clean_uniprot_reviewed_values <- function(values){
   values %>%
   stringi::stri_replace_first_fixed('unreviewed', 0)  %>%
   stringi::stri_replace_first_fixed('reviewed',   1)  %>%
   (function(x){x[is.na(x)]<-'0'; x})                  %>%
   as.numeric()
}

#===================================
# Clean uniprot proteinname values
#====================================

PREFIX_RS_PATTERN <- '[(][0-9]+[RS][)]\\-'

#' @rdname rm_prefix_rs
#' @importFrom magrittr %>%
#' @export
has_prefix_rs <- function(values) values %>% stringi::stri_detect_regex(PREFIX_RS_PATTERN)

#' Detect or rm prefix (3R)-
#' @param values uniprot protein name values
#' @examples
#' require(magrittr)
#' values <- paste0("Very-long-chain (3R)-3-hydroxyacyl-CoA dehydratase 2 ",
#'                  "(EC 4.2.1.134) ",
#'                  "(3-hydroxyacyl-CoA dehydratase 2) ",
#'                  "(HACD2) ",
#'                  "(Protein-tyrosine phosphatase-like member B)")
#' values %>% has_prefix_rs()
#' values %>% rm_prefix_rs()
#' @importFrom magrittr %>%
#' @export
rm_prefix_rs <- function(values) values %>% gsub(PREFIX_RS_PATTERN, '', .)

#' Rm synonyms from uniprot proteinname values
#' @param values uniprot proteinname values
#' @return character vector with standard proteinname values
#' @examples
#' require(magrittr)
#' values <- paste0("Rho GTPase-activating protein 10 ",
#'                  "(GTPase regulator associated with focal adhesion kinase 2) ",
#'                  "(Graf-related protein 2) ",
#'                  "(Rho-type GTPase-activating protein 10)")
#' values %>% rm_protein_synonyms()
#' @importFrom magrittr %>%
#' @export
rm_protein_synonyms <- function(values){
   values                                              %>%
   gsub("\\((?>[^()]|(?R))*\\)", "", ., perl = TRUE)   %>%
   trimws()                                            %>%
   gsub(' +', ' ', .)                                  %>% # multiple white spaces -> single white space
   gsub(' ([,;])', '\\1', .)                           %>% # rm white spaces before "," or ";"
   gsub(' ]', ']', .)                                      # rm white spaces before "]". Occurs in [Includes: ] or [Cleaved into: ] constructs
}

INCLUDES_PATTERN <- ' \\[Includes:.+]'

#' @rdname rm_includes
#' @importFrom magrittr %>%
#' @export
has_includes <- function(values) values %>% stringi::stri_detect_regex(INCLUDES_PATTERN)

#' rm "includes" constructs from protein name values
#' @param values uniprot protein name values
#' @return updated values (include constructs removed)
#' @examples
#' require(magrittr)
#' values <- paste0("Bifunctional 3'-phosphoadenosine 5'-phosphosulfate synthase 1 ",
#'                  "(PAPS synthase 1) (PAPSS 1) (Sulfurylase kinase 1) (SK 1) (SK1) ",
#'                  "[Includes: Sulfate adenylyltransferase (EC 2.7.7.4) (ATP-sulfurylase) ",
#'                  "(Sulfate adenylate transferase) (SAT); Adenylyl-sulfate kinase (EC 2.7.1.25) ",
#'                  "(3'-phosphoadenosine-5'-phosphosulfate synthase) (APS kinase) ",
#'                  "(Adenosine-5'-phosphosulfate 3'-phosphotransferase) ",
#'                  "(Adenylylsulfate 3'-phosphotransferase)]")
#' values %>% has_includes()
#' values %>% rm_includes()
#' @importFrom magrittr %>%
#' @export
rm_includes <- function(values)   values %>% stringi::stri_replace_all_regex(INCLUDES_PATTERN, '')


CLEAVED_PATTERN <- ' \\[Cleaved into:.+]'

#' @rdname rm_cleaved
#' @importFrom magrittr %>%
#' @export
has_cleaved <- function(values) values %>% stringi::stri_detect_regex(CLEAVED_PATTERN)

#' rm "includes" constructs from protein name values
#' @param values uniprot protein name values
#' @return updated values (include constructs removed)
#' @examples
#' require(magrittr)
#' values <- paste0("Dystroglycan (Dystrophin-associated glycoprotein 1) ",
#'                  "[Cleaved into: Alpha-dystroglycan (Alpha-DG); Beta-dystroglycan (Beta-DG)]")
#' values %>% has_cleaved()
#' values %>% rm_cleaved()
#' @importFrom magrittr %>%
#' @export
rm_cleaved <- function(values)   values %>% stringi::stri_replace_all_regex(CLEAVED_PATTERN, '')


#' Clean uniprot "protein name" values
#' @param values uniprot protein name values (character vector)
#' @return cleaned uniprot protein name values (character vector)
#' @examples
#' protein1 <- paste0("Rho GTPase-activating protein 10 ",
#'                    "(GTPase regulator associated with focal adhesion kinase 2) ",
#'                    "(Graf-related protein 2) ",
#'                    "(Rho-type GTPase-activating protein 10)")
#' protein2 <- paste0("Dystroglycan (Dystrophin-associated glycoprotein 1) ",
#'                    "[Cleaved into: Alpha-dystroglycan (Alpha-DG); Beta-dystroglycan (Beta-DG)]")
#' protein3 <- paste0("Bifunctional 3'-phosphoadenosine 5'-phosphosulfate synthase 1 ",
#'                    "(PAPS synthase 1) (PAPSS 1) (Sulfurylase kinase 1) (SK 1) (SK1) ",
#'                    "[Includes: Sulfate adenylyltransferase (EC 2.7.7.4) (ATP-sulfurylase) ",
#'                    "(Sulfate adenylate transferase) (SAT); Adenylyl-sulfate kinase (EC 2.7.1.25) ",
#'                    "(3'-phosphoadenosine-5'-phosphosulfate synthase) (APS kinase) ",
#'                    "(Adenosine-5'-phosphosulfate 3'-phosphotransferase) ",
#'                    "(Adenylylsulfate 3'-phosphotransferase)]")
#' values <- c(protein1, protein2, protein3)
#' clean_uniprot_proteinname_values(values)
#' @importFrom magrittr %>%
#' @export
clean_uniprot_proteinname_values <- function(values){
   # Satisfy CHECK
   . <- NULL

   # Rm synonyms (...) in protein names
   # https://stackoverflow.com/questions/41749058/r-parse-nested-parentheses
   # 922
   # 1000
   # 415
   values                             %>%
   as.character()                     %>%
   rm_prefix_rs                       %>%
   rm_protein_synonyms()              %>%
   rm_includes()                      %>%
   rm_cleaved()                       %>%
   (function(x){x[is.na(x)]<-''; x})
}

#===============================

#' Clean uniprot "gene" values
#' @param values uniprot "gene" values (character vector)
#' @return cleaned uniprot gene values (character vector)
#' @examples
#' \dontrun{
#' require(magrittr)
#' values <- c("A0A024R4M0", "M0R210", "G3HSF3")
#' up <- values %>% autonomics.annotate::connect_to_uniprot()
#' annotations <- values %>% autonomics.annotate::fetch_uniprot_annotations(up)
#' annotations$GENE %>% print()
#' annotations$GENE %>% autonomics.annotate::clean_uniprot_gene_values() %>% print()
#' }
#' @importFrom magrittr %>%
#' @export
clean_uniprot_gene_values <- function(values){
   values                           %>%
   stringi::stri_split_fixed(' ')   %>%
   vapply(magrittr::extract, character(1), 1) %>%
  (function(x){x[is.na(x)]<-''; x})
}

#========================================================

#' Clean uniprot location values
#' @param values character vector
#' @return character vector
#' @examples
#' \dontrun{
#'    uniprot_vector <- c('P83940', 'E9Q4P1', 'E9Q3S4', 'Q9JKP5', 'Q8BVE3')
#'    up <- autonomics.annotate::connect_to_uniprot(uniprot_vector)
#'    dt <- AnnotationDbi::select(up, keys = uniprot_vector, columns = 'SUBCELLULAR-LOCATIONS')
#'    values <- dt$`SUBCELLULAR-LOCATIONS`
#' }
#' @importFrom magrittr %>%
#' @export
clean_uniprot_location_values <- function(values){
   values %>% stringi::stri_replace_last_regex(' Note=.+$', '') %>%
      stringi::stri_replace_all_fixed('SUBCELLULAR LOCATION: ', '') %>%
      stringi::stri_replace_all_regex(' \\{[^ ]+\\}', '')
}


#========================================================

#' Annotate uniprot ids
#' @param values uniprot accessions (character vector)
#' @param up uniprot.ws connection
#' @return data.table with uniprot annotations
#' @examples
#' \dontrun{
#' require(magrittr)
#' values <- c("A0A024R4M0", "M0R210", "G3HSF3", "P31941")
#' up <- values %>% autonomics.annotate::connect_to_uniprot()
#' autonomics.annotate::annotate_uniprot_through_ws(values, up) %>% print()
#' }
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%   %<>%
#' @export
annotate_with_uniprot_webservice <- function(
   values,
   up      = autonomics.annotate::connect_to_uniprot(values),
   columns = c('SUBCELLULAR-LOCATIONS', 'KEGG', 'GO-ID', 'INTERPRO')
){
   # Fetch uniprot annotations
   dt <- AnnotationDbi::select(up, keys = values, columns = columns)
   dt %<>% data.table::data.table()
   data.table::setnames('UNIPROTKB', 'Uniprot accessions')

   # Collapse redundant kegg ids
   if ('KEGG' %in% columns){
      n0 <- nrow(dt)
      dt %<>% magrittr::extract(, ('KEGG') := paste0(get('KEGG'), collapse = ';'), by = 'UNIPROTKB') %>% unique()
      n1 <- nrow(dt)
      autonomics.support::cmessage('\tCollapse %d KEGG ids mapping to same uniprot accession: %d -> %d features', n0-n1, n0, n1)
      data.table::setnames('KEGG',     'keggid')
   }

   # Clean values
   if ('SCORE'                 %in% columns)  dt %>% magrittr::extract(, SCORE                  := autonomics.annotate::clean_uniprot_score_values(       SCORE)                 )
   if ('EXISTENCE'             %in% columns)  dt %>% magrittr::extract(, EXISTENCE              := autonomics.annotate::clean_uniprot_existence_values(   EXISTENCE)             )
   if ('REVIEWED'              %in% columns)  dt %>% magrittr::extract(, REVIEWED               := autonomics.annotate::clean_uniprot_reviewed_values(    REVIEWED)              )
   if ('PROTEIN-NAMES'         %in% columns)  dt %>% magrittr::extract(,`PROTEIN-NAMES`         := autonomics.annotate::clean_uniprot_proteinname_values(`PROTEIN-NAMES`)        )
   if ('GENES'                 %in% columns)  dt %>% magrittr::extract(,`GENES`                 := autonomics.annotate::clean_uniprot_gene_values(        GENES)                 )
   if ('SUBCELLULAR-LOCATIONS' %in% columns)  dt %>% magrittr::extract(,`SUBCELLULAR-LOCATIONS` := autonomics.annotate::clean_uniprot_location_values(   `SUBCELLULAR-LOCATIONS`))

   # Return
   dt
}

#==============================================================

#' Annotate proteingroups
#' @param object SummerizedExperiment with proteinGroups data
#' @param fastafile path to fastafile
#' @param webservice NULL or character vector (annotations to fetch from uniprot webservice)
#' @return SummarizedExperiment with annotations
annotate_proteingroups <- function(
   object,
   fastafile,
   webservice_columns= NULL,
   webservice_connection  = NULL
){

   # Annotate by fasta file
   fdata1 <- object %>%
             autonomics.import::fdata() %>%
             magrittr::extract(, c('feature_id', 'Uniprot accessions'), drop = FALSE) %>%
             tidyr::separate_rows(`Uniprot accessions`, sep = ';') %>%
             data.table::data.table()
   annotations <- fdata1$`Uniprot accessions` %>% annotate_with_uniprot_fastafile(fastafile)
   fdata1 %<>% merge(annotations, by.x = 'Uniprot accessions', by.y = 'uniprot', sort = FALSE)

   # Prefer swissprot over trembl
   fdata1 %<>% magrittr::extract(, .SD[reviewed == max(reviewed)], by = 'feature_id')

   # Prefer canonical when present
   fdata1 %>%  magrittr::extract(, is.canonical := as.numeric(reviewed==TRUE & !stringi::stri_detect_fixed(`Uniprot accessions`, '-')))
   fdata1 %<>% magrittr::extract(, .SD[is.canonical == max(is.canonical)], by = 'feature_id')
   fdata1 %<>% magrittr::extract(, is.canonical := NULL)
   fdata1 %<>% magrittr::extract(, reviewed     := NULL)

   # Prefer best protein existence
   fdata1[is.na(existence), existence:=5]
   fdata1 %<>% magrittr::extract(, .SD[existence == min(existence)], by = 'feature_id')
   fdata1[, existence := NULL]

   # Annotate through webservice
   if (!is.null(annotate_through_webservice)){
      annotations2 <- fdata1$`Uniprot accessions` %>%
                      autonomics.annotate::annotate_with_uniprot_webservice(up = webservice_connection, columns = webservice_columns)
      # Merge in
   }

   # Collapse
   fdata1 %<>% magrittr::extract(, list(`Uniprot accessions` = `Uniprot accessions` %>% paste0(collapse = ';'),
                                        genename            =  unique(genename)     %>% paste0(collapse = ';'),
                                        proteinname         =  unique(proteinname)  %>% paste0(collapse = ';')), by = 'feature_id')

   # Merge back
   nullify_fvars <- function(object, fvars){
      for (curfvar in fvars)   autonomics.import::fdata(object)[[curfvar]] <- NULL
      return(object)
   }
   object %<>% nullify_fvars(fvars = c('Uniprot accessions', 'Protein names', 'Gene names'))
   autonomics.import::fdata(object) %<>% merge(fdata1, by = 'feature_id', sort = FALSE)

   # Return
   return(object)

}
