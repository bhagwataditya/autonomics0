
#' Constants in autonomics.ora
#' 
#' Which organisms have available annotations?
#' Which GO collections are supported?
#' Which gene set collections are supported?
#' 
#' @name Constants
#' @examples 
#' ANNOTATED_ORGANISMS
#' PANTHER_ORGANISMS
#' GO_COLLECTIONS
#' GENE_SET_COLLECTIONS
NULL

#' @rdname Constants
#' @export
ANNOTATED_ORGANISMS <- list(
   `Canis familiaris` = c(
      abrev               = 'Cf',
      source              = 'eg',
      kegg                = 'cfa',
      panther             = 'CANINE',
      taxonid             = '9615',
      kegg_pathway_suffix = ' - Canis familiaris (dog)'),
   `Danio rerio` = c(
      abrev               = 'Dr',
      source              = 'eg',
      kegg                = 'dre',
      panther             = 'ZEBRAFISH',
      taxonid             = '7955',
      kegg_pathway_suffix = ' - Danio rerio (zebra fish)'), 
   `Homo sapiens` = c(
      abrev               = 'Hs',
      source              = 'eg',
      kegg                = 'hsa',
      panther             = 'HUMAN',
      taxonid             = '9606',
      kegg_pathway_suffix = ' - Homo sapiens (human)'),
   `Mus musculus` = c(
      abrev               = 'Mm',
      source              = 'eg',
      kegg                = 'mmu',
      panther             = 'MOUSE',
      taxonid             = '10090',
      kegg_pathway_suffix = ' - Mus musculus (mouse)'),
   `Rattus norvegicus` = c(
      abrev               = 'Rn',
      source              = 'eg',
      kegg                = 'rno',
      panther             = 'RAT',
      taxonid             = '10116',
      kegg_pathway_suffix = ' - Rattus norvegicus (rat)'),
   `Sacharomyces cerevisiae` = c(
      abrev = 'Sc',
      source = 'sgd',
      kegg = 'sce',
      panther = 'YEAST',
      taxonid = '559292',
      kegg_pathway_suffix = ' - Saccharomyces cerevisiae (budding yeast)'),
   `Xenopus laevis` = c(
      abrev               = 'Xl',
      source              = 'eg',
      kegg                = 'xla',
      panther             = NA_character_,
      taxonid             = '8355',
      kegg_pathway_suffix = ' - Xenopus laevis (African clawed frog)')
)

#' @rdname Constants
#' @export
PANTHER_ORGANISMS <- ANNOTATED_ORGANISMS %>% 
                     vapply(magrittr::extract2, character(1), 'panther') %>% 
                     extract(!is.na(.))

# xenopus    = c(abrev = 'Xl', source = 'eg',  eg = 'Xl', kegg = 'xla', panther = 'XENOPUS',   taxonid = '8364'), 
   # AMB 22 Feb 2017: I removed  Xenopus of this list because this needs to be sorted out
   # Panther is based on Xenopus tropicalis (taxon id 8364)
   # Entrezg is based on Xenopus laevis     (taxon id 8355)
   # We often use concatenated FASTA dbs (X. tropicalis + X. laevis)


#' @rdname Constants
#' @export
GO_COLLECTIONS <- c('gobp', 'gocc', 'gomf')

#' @rdname Constants
#' @export
GENE_SET_COLLECTIONS <- c(GO_COLLECTIONS, 'kegg')
#GENE_SET_COLLECTIONS <- c(GO_COLLECTIONS, 'kegg', 'uniprot')
