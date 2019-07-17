
#==============================================================================
# READ UNIPROT ANNOTATIONS FROM FASTAFILE
#==============================================================================

extract_accession <- function(fastahdrs){
    message('\t\t\tExtract ACCESSION')
    fastahdrs                         %>%
    stri_split_fixed('|')             %>%
    vapply(extract, character(1), 2)
}

extract_reviewed <- function(fastahdrs){
    message('\t\t\tExtract REVIEWED: 0=trembl, 1=swissprot')
    fastahdrs                         %>%
    stri_split_fixed('|')             %>%
    vapply(extract, character(1), 1)  %>%
    equals('sp')                      %>%
    as.numeric()
}

extract_entryname <- function(fastahdrs){
    message('\t\t\tExtract ENTRYNAME')
    fastahdrs                         %>%
    stri_split_fixed('|')             %>%
    vapply(extract, character(1), 3)
}


splitoff_VERSION <- function(dt){
    VERSION <- annotation <- NULL
    message('\t\t\tExtract VERSION')
    pattern <- ' SV=[0-9]'
    dt [ , VERSION :=   annotation                        %>%
                        stri_extract_last_regex(pattern)  %>%
                        substr(5,5)                       %>%
                        as.numeric() ]
    dt [ , annotation:= annotation %>%
                        stri_replace_last_regex(pattern, '') ]
}


splitoff_EXISTENCE <- function(dt){
    EXISTENCE <- annotation <- NULL
    pattern <- ' PE=[0-9]'
    message('\t\t\tExtract EXISTENCE: 1=protein, 2=transcript, ',
            '3=homology, 4=prediction, 5=uncertain, NA=isoform')
    dt  [ , EXISTENCE  :=   annotation %>%
                            stri_extract_last_regex(pattern) %>%
                            substr(5,5) %>%
                            as.numeric() ]
    dt  [ , annotation :=   annotation %>%
                            stri_replace_last_regex(pattern, '') ]
}

splitoff_GENES <- function(dt, gene_var){
    GENES <- annotation <- NULL
    pattern <- ' GN=.+$'
    autonomics.support::cmessage('\t\t\tExtract GENES')
    dt  [ , GENES :=    annotation %>%
                        stri_extract_last_regex(pattern) %>%
                        substr(5,nchar(.)) ]
    dt  [ , annotation :=   annotation %>%
                            stri_replace_last_regex(pattern, '') ]
}

splitoff_ORGID <- function(dt, orgid_var){
    ORGID <- annotation <- NULL
    pattern <- ' OX=[0-9]+'
    autonomics.support::cmessage('\t\t\tExtract ORGID')
    dt [ , ORGID :=     annotation %>%
                        stri_extract_last_regex(pattern) %>%
                        substr(5,nchar(.)) ]
    dt [ , annotation :=    annotation %>%
                            stri_replace_last_regex(pattern, '') ]
}

splitoff_ORGNAME <- function(dt, orgname_var){
    ORGNAME <- annotation <- NULL
    pattern <- ' OS=.+$'
    autonomics.support::cmessage('\t\t\tExtract ORGNAME')
    dt [ , ORGNAME :=   annotation %>%
                        stri_extract_last_regex(pattern) %>%
                        substr(5,nchar(.)) ]
    dt [ , annotation :=    annotation %>%
                            stri_replace_last_regex(pattern, '') ]
}

splitoff_PROTEINNAMES <- function(dt, proteinname_var){
    `PROTEIN-NAMES` <- annotation <- NULL
    pattern <- ' .+$'
    autonomics.support::cmessage('\t\t\tExtract PROTEINNAMES')
    dt [ , `PROTEIN-NAMES` :=   annotation %>%
                                stri_extract_last_regex(pattern) %>%
                                substr(2,nchar(.)) ]
    dt [ , annotation      :=   annotation %>%
                                stri_replace_last_regex(pattern, '') ]
}

#' Read uniprot annotations from fastafile
#'
#' @param fastafile    path to fasta file
#' @param fastafields  character vector
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'     fastafile <- system.file(
#'                      'extdata/fasta/uniprot_hsa_20140515.fasta',
#'                      package = 'autonomics.data')
#'     fastafile %>% read_uniprot_annotations()
#' }
#' @return data.table (uniprot, genename, proteinname, reviewed, existence)
#' @note EXISTENCE values are always those of canonical isoform
#' @importFrom data.table   data.table   :=
#' @import stringi
#' @export
read_uniprot_annotations <- function(
    fastafile,
    fastafields = c('GENES', 'PROTEIN-NAMES', 'REVIEWED', 'EXISTENCE',
                    'ENTRYNAME', 'ORGNAME', 'VERSION')
){

    # Assert & import
    assertive.files::assert_all_are_existing_files(fastafile)
    `:=` <- data.table::`:=`
    annotation <- NULL

    # Load (relevant portion of) fasta
    fasta <- seqinr::read.fasta(fastafile)
    all_accessions <-   names(fasta)          %>%
                        stri_split_fixed('|') %>%
                        vapply(extract, character(1), 2)

    # Extract annotations
     dt <- data.table::data.table(
            UNIPROTKB  = names(fasta) %>% extract_accession(),
            REVIEWED   = names(fasta) %>% extract_reviewed(),
            ENTRYNAME  = names(fasta) %>% extract_entryname(),
            annotation = seqinr::getAnnot(fasta) %>% unlist())

    dt %>% splitoff_VERSION()
    dt %>% splitoff_EXISTENCE()
    dt %>% splitoff_GENES()
    #dt %>% splitoff_ORGID() # annotations don't seem to contain this?!
    dt %>% splitoff_ORGNAME()
    dt %>% splitoff_PROTEINNAMES()
    dt  [ , annotation     := NULL ]
    dt %<>% extract(, c('UNIPROTKB', fastafields), with = FALSE)

    # Return
    return(dt)
}

#=============================================================================
# READ UNIPROT ANNOTATIONS FROM UNIPROT WEBSERVICE
# Functionality below works, but is slow.
# Recommended approach therefore is to use read_uniprot_annotations(fastafile)
#=============================================================================

#' Open a uniprot.ws connection
#'
#' Infers organism from uniprot ids and opens a
#' uniprot webservice connection.
#'
#' @param x uniprot ids
#' @return uniprot.ws connection
#' @examples
#' \dontrun{
#'    require(magrittr)
#'    x <- c("A0A024R4M0", "M0R210", "G3HSF3")
#'    up <- x %>% connect_to_uniprot()
#' }
#' @export
connect_to_uniprot <- function(x){
    . <- NULL
    x                                            %>%
    infer_organism('uniprot', verbose = FALSE)   %>%
    extract2(ANNOTATED_ORGANISMS, .)             %>%
    extract2('taxonid') %>%
    as.numeric() %>%
    UniProt.ws::UniProt.ws()
}

#----------------------------------------------

#' Fetch uniprot annotations
#' @param x uniprot ids (character vector)
#' @param up uniprot.ws connection (output of connect_to_uniprot)
#' @param columns character vector
#' @return dataframe with annotations
#' @examples
#' \dontrun{
#'     require(magrittr)
#'     x <- c("A0A024R4M0", "M0R210", "G3HSF3")
#'     up <- x %>% connect_to_uniprot()
#'     x %>% fetch_uniprot_annotations(up)
#' }
#' @export
fetch_uniprot_annotations <- function(
    x,
    up,
    columns = c('SUBCELLULAR-LOCATIONS', 'KEGG', 'GO-ID', 'INTERPRO')
){
    AnnotationDbi::select(up, keys = x, columns = columns)
}

#-------------------------------------------------

#' Clean uniprot "score" values
#' @param x character vector
#' @return character vectors
#' @examples
#' \dontrun{
#'     require(magrittr)
#'     x <- c("A0A024R4M0", "M0R210", "G3HSF3")
#'     up <- values %>% connect_to_uniprot()
#'     annotations <- x %>% fetch_uniprot_annotations(up)
#'     annotations$SCORE %T>% print() %>% clean_score()
#' }
#' @export
clean_score <- function(x){
    x                                     %>%
    substr(1, 1)                          %>%
    (function(x){x[is.na(x)] <- '0'; x})  %>%
    as.numeric()
}

#--------------------------------------------------

#' Clean uniprot "protein existence" values
#' @param x uniprot protein existence values (character vector)
#' @return cleaned values (character vector)
#' @examples
#' \dontrun{
#' require(magrittr)
#'     x <- c("A0A024R4M0", "M0R210", "G3HSF3")
#'     up <- x %>% connect_to_uniprot()
#'     annotations <- x %>% fetch_uniprot_annotations(up,'EXISTENCE')
#'     annotations$EXISTENCE %T>% print() %>% clean_existence()
#' }
#' @import stringi
#' @export
clean_existence <- function(x){
    x                                                             %>%
    stri_replace_first_fixed("Evidence at protein level",    1)   %>%
    stri_replace_first_fixed("Evidence at transcript level", 2)   %>%
    stri_replace_first_fixed("Inferred from homology",       3)   %>%
    stri_replace_first_fixed("Predicted",                    4)   %>%
    stri_replace_first_fixed("Uncertain",                    5)   %>%
    (function(y){y[is.na(y)] <- '5'; y})                          %>%
    as.numeric()
}

#------------------------------------------------------

#' Clean uniprot "reviewed" values
#' @param x uniprot "reviewed" values (character vector)
#' @return cleaned values (character vector)
#' @examples
#' \dontrun{
#' require(magrittr)
#' x <- c("A0A024R4M0", "M0R210", "G3HSF3")
#' up <- x %>% connect_to_uniprot()
#' annotations <- x %>% fetch_uniprot_annotations(up)
#' annotations$REVIEWED %T>% print() %>% clean_reviewed()
#' }
#' @import stringi
#' @export
clean_reviewed <- function(x){
    x                                          %>%
    stri_replace_first_fixed('unreviewed', 0)  %>%
    stri_replace_first_fixed('reviewed',   1)  %>%
    (function(x){x[is.na(x)]<-'0'; x})         %>%
    as.numeric()
}

#--------------------------------------------------------

# Clean uniprot proteinname values

PREFIX_RS_PATTERN <- '[(][0-9]+[RS][)]\\-'


#' @rdname rm_prefix_rs
#' @export
has_prefix_rs <- function(x) x %>% stri_detect_regex(PREFIX_RS_PATTERN)


#' Detect or rm prefix (3R)-
#' @param x uniprot protein name values
#' @examples
#' require(magrittr)
#' x <- paste0(
#'     "Very-long-chain (3R)-3-hydroxyacyl-CoA dehydratase 2 ",
#'     "(EC 4.2.1.134) ",
#'     "(3-hydroxyacyl-CoA dehydratase 2) ",
#'     "(HACD2) ",
#'     "(Protein-tyrosine phosphatase-like member B)")
#' x %>% has_prefix_rs()
#' x %>% rm_prefix_rs()
#' @export
rm_prefix_rs <- function(x) x %>% gsub(PREFIX_RS_PATTERN, '', .)


#' Rm synonyms from uniprot proteinname values
#' @param x uniprot proteinname values
#' @return character vector with standard proteinname values
#' @examples
#' require(magrittr)
#' x <- paste0(
#'     "Rho GTPase-activating protein 10 ",
#'     "(GTPase regulator associated with focal adhesion kinase 2) ",
#'     "(Graf-related protein 2) ",
#'     "(Rho-type GTPase-activating protein 10)")
#' x %>% rm_protein_synonyms()
#' @export
rm_protein_synonyms <- function(x){
    x                                                   %>%
    gsub("\\((?>[^()]|(?R))*\\)", "", ., perl = TRUE)   %>%
    trimws()                                            %>%
    # multiple white spaces -> single white space
    gsub(' +', ' ', .)                                  %>%
    # rm white spaces before "," or ";"
    gsub(' ([,;])', '\\1', .)                           %>%
    # rm white spaces before "]".
    # Occurs in [Includes: ] or [Cleaved into: ] constructs
    gsub(' ]', ']', .)
}


INCLUDES_PATTERN <- ' \\[Includes:.+]'

#' @rdname rm_includes
#' @export
has_includes <- function(x) x %>% stri_detect_regex(INCLUDES_PATTERN)


#' rm "includes" constructs from protein name values
#' @param x uniprot protein name values
#' @return updated values (include constructs removed)
#' @examples
#' require(magrittr)
#' values <- paste0(
#'     "Bifunctional 3'-phosphoadenosine 5'-phosphosulfate synthase 1 ",
#'     "(PAPS synthase 1) (PAPSS 1) (Sulfurylase kinase 1) (SK 1) (SK1) ",
#'     "[Includes: Sulfate adenase (EC 2.7.7.4) (ATP-sulfurylase) ",
#'     "(Sulfate adenase) (SAT); Adenylyl-sulfate kinase (EC 2.7.1.25) ",
#'     "(3'-phosphoadenosine-5'-phosphosulfate synthase) (APS kinase) ",
#'     "(Adenosine-5'-phosphosulfate 3'-phosphotransferase) ",
#'     "(Adenylylsulfate 3'-phosphotransferase)]")
#' values %>% has_includes()
#' values %>% rm_includes()
#' @export
rm_includes <- function(x)   x %>% stri_replace_all_regex(INCLUDES_PATTERN, '')


CLEAVED_PATTERN <- ' \\[Cleaved into:.+]'

#' @rdname rm_cleaved
#' @export
has_cleaved <- function(x) x %>% stri_detect_regex(CLEAVED_PATTERN)


#' rm "includes" constructs from protein name values
#' @param x uniprot protein name values
#' @return updated values (include constructs removed)
#' @examples
#' require(magrittr)
#' values <- paste0(
#'     "Dystroglycan (Dystrophin-associated glycoprotein 1) ",
#'     "[Cleaved into: Alpha-dglycan (Alpha-DG); Beta-dystroglycan (Beta-DG)]")
#' values %>% has_cleaved()
#' values %>% rm_cleaved()
#' @export
rm_cleaved <- function(x)  x %>% stri_replace_all_regex(CLEAVED_PATTERN, '')


#' Clean uniprot "protein name" values
#' @param x uniprot protein name values (character vector)
#' @return cleaned uniprot protein name values (character vector)
#' @examples
#' protein1 <- paste0(
#'     "Rho GTPase-activating protein 10 ",
#'     "(GTPase regulator associated with focal adhesion kinase 2) ",
#'     "(Graf-related protein 2) ",
#'     "(Rho-type GTPase-activating protein 10)")
#' protein2 <- paste0(
#'     "Dystroglycan (Dystrophin-associated glycoprotein 1) ",
#'     "[Cleaved into: Alpha-dglycan (Alpha-DG); Beta-dystroglycan (Beta-DG)]")
#' protein3 <- paste0(
#'     "Bifunctional 3'-phosphoadenosine 5'-phosphosulfate synthase 1 ",
#'     "(PAPS synthase 1) (PAPSS 1) (Sulfurylase kinase 1) (SK 1) (SK1) ",
#'     "[Includes: Sulfate adenylyltransferase (EC 2.7.7.4) (ATP-sulfurylase) ",
#'     "(Sulfate adenylase) (SAT); Adenylyl-sulfate kinase (EC 2.7.1.25) ",
#'     "(3'-phosphoadenosine-5'-phosphosulfate synthase) (APS kinase) ",
#'     "(Adenosine-5'-phosphosulfate 3'-phosphotransferase) ",
#'     "(Adenylylsulfate 3'-phosphotransferase)]")
#' x <- c(protein1, protein2, protein3)
#' clean_proteinnames(x)
#' @export
clean_proteinnames <- function(x){
    # Satisfy CHECK
    . <- NULL

    # Rm synonyms (...) in protein names
    # https://stackoverflow.com/questions/41749058/r-parse-nested-parentheses
    # 922
    # 1000
    # 415
    x                                  %>%
    as.character()                     %>%
    rm_prefix_rs                       %>%
    rm_protein_synonyms()              %>%
    rm_includes()                      %>%
    rm_cleaved()                       %>%
    (function(y){y[is.na(y)]<-''; y})
}

#------------------------------------------------------------

#' Clean uniprot "gene" values
#' @param x uniprot "gene" values (character vector)
#' @return cleaned uniprot gene values (character vector)
#' @examples
#' \dontrun{
#' require(magrittr)
#' x <- c("A0A024R4M0", "M0R210", "G3HSF3")
#' up <- x %>% connect_to_uniprot()
#' annotations <- x %>% fetch_uniprot_annotations(up)
#' annotations$GENE %>% print()
#' annotations$GENE %>% clean_genes() %>% print()
#' }
#' @import stringi
#' @export
clean_genes <- function(x){
    x                                 %>%
    stri_split_fixed(' ')             %>%
    vapply(extract, character(1), 1)  %>%
    (function(y){y[is.na(y)]<-''; y})
}

#---------------------------------------------------------------------

#' Clean uniprot location values
#' @param x character vector
#' @return character vector
#' @examples
#' require(magrittr)
#' \dontrun{
#'     uniprot_vector <- c('A0MZ66', 'A1KZ92', 'A1L020', 'A1XBS5',
#'                         'S4R350', 'A4D1P6')
#'     up <- connect_to_uniprot(uniprot_vector)
#'     dt <- AnnotationDbi::select(
#'             up,
#'             keys    = uniprot_vector,
#'             columns = 'SUBCELLULAR-LOCATIONS') %>%
#'           data.table::data.table()
#'    x <- dt$`SUBCELLULAR-LOCATIONS`
#'    x %>% clean_locations()
#' }
#' @import stringi
#' @export
clean_locations <- function(x){
    x                                                     %>%
    # rm Notes (not part of controlled vocabulary)
    stri_replace_last_regex(' Note=.+$', '')              %>%
    stri_replace_all_fixed('SUBCELLULAR LOCATION: ', '')  %>%
    # make pattern lazy by using [^.;] construct
    stri_replace_all_regex(' \\{[^;.]+\\}', '')
}


#---------------------------------------------------------------

#' Annotate uniprot ids
#' @param x      uniprot accessions (character vector)
#' @param connection  uniprot.ws connection
#' @param columns     character vector
#' @return data.table with uniprot annotations
#' @examples
#' \dontrun{
#'     require(magrittr)
#'     x <- c('A0MZ66', 'A1KZ92', 'A1L020', 'A1XBS5', 'S4R350', 'A4D1P6')
#'     up <- x %>% connect_to_uniprot()
#'     annotate_uniprot_with_webservice(
#'        x, up, columns =  'SUBCELLULAR-LOCATIONS') %>% print()
#' }
#' @importFrom data.table   data.table   :=
#' @export
annotate_uniprot_with_webservice <- function(
    x,
    connection = connect_to_uniprot(x),
    columns    = c('SUBCELLULAR-LOCATIONS', 'KEGG', 'GO-ID', 'INTERPRO')
){
    # Assert
    assertive.base::assert_is_identical_to_true(
        class(connection) == 'UniProt.ws')
    assertive.sets::assert_is_subset(columns, UniProt.ws::columns(connection))
    SCORE <- EXISTENCE <- REVIEWED <- `PROTEIN-NAMES` <- NULL
    GENES <- `SUBCELLULAR-LOCATIONS` <- NULL

    # Fetch uniprot annotations
    dt <- AnnotationDbi::select(connection, keys = x, columns = columns)
    dt %<>% data.table::data.table()

    # Collapse redundant kegg ids
    if ('KEGG' %in% columns){
        n0 <- nrow(dt)
        dt %<>% extract(, ('KEGG') := paste0(get('KEGG'), collapse = ';'),
                        by = 'UNIPROTKB') %>% unique()
        n1 <- nrow(dt)
        autonomics.support::cmessage(
            paste0('\tCollapse %d KEGG ids mapping to ',
                   'same uniprot accession: %d -> %d features'), n0-n1, n0, n1)
        #data.table::setnames('KEGG',     'keggid')
    }

    # Clean values
    if ('SCORE' %in% columns){
        dt [, SCORE           := clean_score(SCORE) ]
    }
    if ('EXISTENCE' %in% columns){
        dt [, EXISTENCE      := clean_existence(EXISTENCE) ]
    }
    if ('REVIEWED' %in% columns){
        dt [, REVIEWED       := clean_reviewed(REVIEWED) ]
    }
    if ('PROTEIN-NAMES' %in% columns){
        dt[, `PROTEIN-NAMES` := clean_proteinnames(`PROTEIN-NAMES`)]
    }
    if ('GENES' %in% columns){
        dt [,`GENES`         := clean_genes(GENES) ]
    }
    if ('SUBCELLULAR-LOCATIONS' %in% columns){
        dt [,`SUBCELLULAR-LOCATIONS` := clean_locations(`SUBCELLULAR-LOCATIONS`)]
    }
    # Return
    dt
}
