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
      columns = c('GENES', 'PROTEIN-NAMES',
                'SCORE', 'REVIEWED', 'EXISTENCE',
                'KEGG',
                'GO-ID', 'INTERPRO'))
}

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

#==========================
# Clean uniprot name values
#==========================

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

#' Annotate uniprot ids
#' @param values uniprot accessions (character vector)
#' @param up uniprot.ws connection
#' @return data.table with uniprot annotations
#' @examples
#' \dontrun{
#' require(magrittr)
#' values <- c("A0A024R4M0", "M0R210", "G3HSF3", "P31941")
#' up <- values %>% autonomics.annotate::connect_to_uniprot()
#' autonomics.annotate::annotate_uniprot(values, up) %>% print()
#' }
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%   %<>%
#' @export
annotate_uniprot <- function(
   values,
   up = autonomics.annotate::connect_to_uniprot(values)
){
   # Fetch uniprot annotations
   dt <- AnnotationDbi::select(up,
                               keys = values,
                               columns = c('GENES', 'PROTEIN-NAMES',
                                           'SCORE', 'REVIEWED', 'EXISTENCE',
                                           'KEGG',
                                           'GO-ID', 'INTERPRO'))
   # Collapse redundant kegg ids
   dt %<>% data.table::data.table()
   #n0 <- nrow(dt)
   dt %<>% magrittr::extract(, ('KEGG') := paste0(get('KEGG'), collapse = ';'), by = 'UNIPROTKB') %>% unique()
   #n1 <- nrow(dt)
   #autonomics.support::cmessage('\tCollapse %d KEGG ids mapping to same uniprot accession: %d -> %d features', n0-n1, n0, n1)

   # Clean values
   dt %<>% magrittr::extract(,  ('SCORE')          := autonomics.annotate::clean_uniprot_score_values(      get('SCORE')))         %>%
           magrittr::extract(,  ('EXISTENCE')      := autonomics.annotate::clean_uniprot_existence_values(  get('EXISTENCE')))     %>%
           magrittr::extract(,  ('REVIEWED')       := autonomics.annotate::clean_uniprot_reviewed_values(   get('REVIEWED')))      %>%
           magrittr::extract(,  ('PROTEIN-NAMES')  := autonomics.annotate::clean_uniprot_proteinname_values(get('PROTEIN-NAMES'))) %>%
           magrittr::extract(,  ('GENES')          := autonomics.annotate::clean_uniprot_gene_values(       get('GENES'))) %>%
           data.table::setnames(c('UNIPROTKB',            'GENES',      'PROTEIN-NAMES', 'SCORE', 'REVIEWED', 'EXISTENCE', 'KEGG',   'GO-ID',  'INTERPRO'),
                                c('Uniprot accessions',   'Gene names', 'Protein names', 'score', 'reviewed', 'existence', 'keggid', 'goid',   'interpro'))
   dt
}


