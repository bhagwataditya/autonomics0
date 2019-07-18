

#' Annotated organisms
#' @export
ANNOTATED_ORGANISMS <- list(
    `Homo sapiens`            = c(  abrev   = 'Hs',
                                    source  = 'eg',
                                    kegg    = 'hsa',
                                    panther = 'HUMAN',
                                    taxonid = '9606'),

    `Mus musculus`            = c(  abrev   = 'Mm',
                                    source  = 'eg',
                                    kegg    = 'mmu',
                                    panther = 'MOUSE',
                                    taxonid = '10090'),

    `Rattus norvegicus`       = c(  abrev   = 'Rn',
                                    source  = 'eg',
                                    kegg    = 'rno',
                                    panther = 'RAT',
                                    taxonid = '10116'),

    `Sacharomyces cerevisiae` = c(  abrev   = 'Sc',
                                    source  = 'sgd',
                                    kegg    = 'sce',
                                    panther = 'YEAST',
                                    taxonid = '559292'),

    `Danio rerio`             = c(  abrev   = 'Dr',
                                    source  = 'eg',
                                    kegg    = 'dre',
                                    panther = 'ZEBRAFISH',
                                    taxonid = '7955'),

    `Xenopus laevis`          = c(  abrev = 'Xl',
                                    source = 'eg',
                                    kegg = 'xla',
                                    panther = NA_character_,
                                    taxonid = '8355')
)


#' Get orgdb name
#' @param organism any value in \code{names(ANNOTATED_ORGANISMS)}
#' @return 'org.Hs.eg.db', org.Mm.eg.db, etc.
#' @examples
#' get_orgdb_name('Homo sapiens')
#' @export
get_orgdb_name <- function(organism){
    assertive.sets::assert_is_subset(organism, names(ANNOTATED_ORGANISMS))
    organism <- ANNOTATED_ORGANISMS %>% extract2(organism)
    orgdb_name <- sprintf('org.%s.%s.db', organism['abrev'], organism['source'])
    return(orgdb_name)
}


#' Load orgdb object
#' @param organism value in ANNOTATED_ORGANISMS
#' @return AnnotationDbi::OrgDb package
#' @examples
#' load_orgdb(organism = 'Homo sapiens')
#' load_orgdb(organism = 'Xenopus laevis')
#' @export
load_orgdb <- function(organism){
    orgdb_name <- get_orgdb_name(organism)
    install_package_if_necessary(orgdb_name)
    utils::getFromNamespace(orgdb_name, orgdb_name)
}

#' Map ensg to entrezg
#' @param x ensg vector
#' @param organism any value in ANNOTATED_ORGANISMS
#' @param verbose logical
#' @return unnamed character vector
#' @examples
#'    x <- c("ENSG00000066248", "ENSG00000066279", "ENSG99999", NA_character_)
#'    ensg_to_entrezg(x)
#' @export
ensg_to_entrezg <- function(
    x,
    organism = infer_organism(keys = x, keytype = 'ensg'),
    verbose = TRUE
){
    entrezg <- rep(NA_character_, length(x))
    idx <- !is.na(x)
    entrezg[idx] <- x[idx] %>%
                    (function(y){
                        idy <- y %in% AnnotationDbi::keys(
                                        load_orgdb(organism),
                                        keytype = 'ENSEMBL')
                        if (verbose){
                            autonomics.support::cmessage(
                                '%d ensg > %d available > %d to entrezg',
                                length(idx), length(idy), sum(idy))
                        }
                        AnnotationDbi::mapIds(
                            load_orgdb(organism),
                            keys    = y,
                            keytype = 'ENSEMBL',
                            column  = 'ENTREZID')
                    })
    entrezg
}

#' Map gsymbol to entrezg
#' @param x        gene symbol vector
#' @param organism any value in ANNOTATED_ORGANISMS
#' @param verbose  logical
#' @return unnamed character vector
#' @examples
#'    x <- c("NGEF", "ASPM", 'ABRACADABRA', NA_character_)
#'    organism <- 'Homo sapiens'
#'    gsymbol_to_entrezg(x, organism)
#' @export
gsymbol_to_entrezg <- function(
   x,
   organism,
   verbose = TRUE
){
    entrezg <- rep(NA_character_, length(x))
    idx <- !is.na(x)
    entrezg[idx] <-  x[idx] %>%
                    (function(y){
                        idy <- y %in% AnnotationDbi::keys(load_orgdb(organism),
                                                          keytype = 'SYMBOL')
                        if (verbose){
                            autonomics.support::cmessage(
                                '%d gsymbols > %d available > %d to entrezg',
                                length(idx), length(idy), sum(idy))
                        }
                        AnnotationDbi::mapIds(
                            load_orgdb(organism),
                            keys    = y,
                            keytype = 'SYMBOL',
                            column  = 'ENTREZID')
                    })
    entrezg
}


#' Load GO genes
#' @param goid       GO ID
#' @param organism  'Homo sapiens'
#' @param gene   'ENTREZID', 'ENSEMBL', 'UNIPROT', 'SYMBOL', 'GENENAME'
#' @return character vector
#' @examples
#' require(magrittr)
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'ENTREZID') %>% head()
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'ENSEMBL')  %>% head()
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'UNIPROT') %>% head()
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'SYMBOL') %>% head()
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'GENENAME') %>% head()
#' @export
load_go_genes <- function(goid, organism, gene = 'ENTREZID'){
   orgdb <- load_orgdb(organism)
   assertive.sets::assert_is_subset(gene, AnnotationDbi::columns(orgdb))
   orgdb %>% AnnotationDbi::select(keys = goid, keytype = 'GO', columns = gene) %>%
             extract2(gene)
}

#' Feature identifer types
#' @export
FEATURE_ID_TYPES <- c(uniprot = 'UNIPROT', ensg = 'ENSEMBL')

build_FEATURE_IDENTIFIERS <- function(){

    FEATURE_IDENTIFIERS <-
        Map(function(feature_id_type){
            Map(function(cur_organism){
                organism <- ANNOTATED_ORGANISMS[[cur_organism]]
                orgdb_name <- sprintf('org.%s.%s.db',
                                      organism['abrev'],
                                      organism['source'])
                install_package_if_necessary(orgdb_name)
                orgdb <- load_orgdb(cur_organism)
                keys <-
                    if(feature_id_type %in% AnnotationDbi::keytypes(orgdb)){
                        AnnotationDbi::keys(orgdb, keytype = feature_id_type)
                    } else {
                        character(0)
                    }
                list(db_pkg         =   orgdb_name,
                     db_pkg_version =   orgdb_name               %>%
                                        utils::packageVersion()  %>%
                                        as.character(),
                     keys           =   keys)
            }, names(ANNOTATED_ORGANISMS))
        }, FEATURE_ID_TYPES) %>%
        set_names(names(FEATURE_ID_TYPES))

    usethis::use_data(
        FEATURE_IDENTIFIERS,
        internal  = TRUE,
        compress  = "xz",
        overwrite = TRUE)

    return(invisible(TRUE))
}


#' Infer organism
#' @param keys    keys
#' @param keytype any value in names(FEATURE_IDENTIFIERS)
#' @param verbose whether to report message (logical)
#' @return organism name
#' @examples
#' require(magrittr)
#' keys <- c('ENSG00000000003',  'ENSG00000000005', 'ENSG00000000419')
#' infer_organism(keys, 'ensg')
#' @export
infer_organism <- function(keys, keytype, verbose = TRUE){

    overlaps <- FEATURE_IDENTIFIERS %>%
                extract2(keytype) %>%
                lapply(extract2, 'keys') %>%
                vapply( function(x){
                            selector <- keys %in% x
                            floor(100*sum(selector)/length(selector))
                        },
                        numeric(1))
    organism <- names(overlaps) %>% extract(which.max(overlaps))
    if (verbose){
        autonomics.support::cmessage('\t\t%s: %d %% %s match %s',
                                    organism,
                                    floor(overlaps[[organism]]),
                                    keytype,
                                    get_orgdb_name(organism))
    }
    return(organism)
}



