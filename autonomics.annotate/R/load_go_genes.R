#' Load GO genes
#' @param goid       GO ID
#' @param organism  'Homo sapiens'
#' @param gene_id   'ENTREZID', 'GENENAME'
#' @return character vector
#' @examples
#' require(magrittr)
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'ENTREZID') %>% head()
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'ENSEMBL')  %>% head()
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'SYMBOL') %>% head()
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'GENENAME') %>% head()
#' load_go_genes(goid = 'GO:0006954', organism = 'Homo sapiens', gene = 'UNIPROT') %>% head()
#' @importFrom magrittr %>%
#' @export
load_go_genes <- function(goid, organism, gene = 'ENTREZID'){
   orgdb <- autonomics.annotate::load_orgdb(organism)
   assertive.sets::assert_is_subset(gene, AnnotationDbi::columns(orgdb))
   orgdb %>% AnnotationDbi::select(keys = goid, keytype = 'GO', columns = gene) %>%
             magrittr::extract2(gene)
}
