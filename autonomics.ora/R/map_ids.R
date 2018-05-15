
#' Map ids
#' 
#' Map a set of ids to another set of ids 
#' @examples
#' library(autonomics.ora)
#' \dontrun{
#'    annmap <- org.Hs.eg.db::org.Hs.eg.db
#'    map_ids(c('1', '2'),                     from = 'ENTREZID', to = 'SYMBOL',   annmap)
#'    map_ids(c('1', '2'),                     from = 'ENTREZID', to = 'UNIPROT',  annmap)
#'    map_ids(c("P04217", "V9HWD8", "P01023"), from = 'UNIPROT',  to = 'ENTREZID', annmap)
#'    map_ids(c("A1BG", "A2M"),                from = 'SYMBOL',   to = 'ENTREZID', annmap)
#' }
#' @param x vector with uniprot accessions / gene symbols
#' @param from 'UNIPROT', 'SYMBOL', 'ENTREZID'
#' @param to   'UNIPROT', 'SYMBOL', 'ENTREZID'
#' @param annotation_map org.xx.eg.db object
#' @importFrom magrittr   %>%    %<>%
#' @export
map_ids <- function(x, from, to, annotation_map){
   
   # Assert
   assertive.types::assert_is_character(x)
   assertive.sets::assert_is_subset(from, c('UNIPROT', 'SYMBOL', 'ENTREZID'))
   assertive.types::assert_is_inherited_from(annotation_map, 'OrgDb')
   
   # Retain those with existing mappings
   x %<>% magrittr::extract(. %in% AnnotationDbi::keys(annotation_map, keytype = from))
   
   # Convert
   if (length(x) == 0){
      return(character(0))
   } else {
      suppressMessages(AnnotationDbi::select(annotation_map, unique(x), to, from)) %>% 
      magrittr::extract(to) %>% 
      magrittr::extract(!is.na(.))  %>% 
      unique()
   }
}
#' Map uniprot accessions to entrezg ids
#' @param x vector with uniprot accessions
#' @param annotation_map org.xx.eg.db object
#' @return vector with entrezg ids
#' @examples
#' \dontrun{
#'    uniprot_set_to_entrezg_set(c("P04217", "V9HWD8", "P01023"), org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @export
uniprot_set_to_entrezg_set <- function(x, annotation_map){
   map_ids(x, from = 'UNIPROT', to = 'ENTREZID', annotation_map)
}

#' Map gene symbols to entrezg ids
#' @param x vector with gene symbols
#' @param annotation_map org.xx.eg.db object
#' @return vector with entrezg ids
#' @examples
#' \dontrun{
#'    gsymbol_set_to_entrezg_set(c("A1BG", "A2M"), org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @export
gsymbol_set_to_entrezg_set <- function(x, annotation_map){
   map_ids(x, from = 'SYMBOL', to = 'ENTREZID', annotation_map)
}


#' Map entrezg ids to uniprot accessions
#' @param x vector with entrezg ids
#' @param annotation_map org.xx.eg.db object
#' @return vector with uniprot accessions
#' @examples
#' \dontrun{
#'    entrezg_set_to_uniprot_set(c('1', '2'), org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @export
entrezg_set_to_uniprot_set <- function(x, annotation_map){
   map_ids(x, from = 'ENTREZID', to = 'UNIPROT', annotation_map = annotation_map)
#   entrezg_set %>%
#         AnnotationDbi::mget(org.Hs.eg.db::org.Hs.egUNIPROT) %>%
#         unlist()                                            %>%
#         unname()                                            %>%
#         na.exclude()
}

#' Map entrezg ids to gene symbols
#' @param x vector with entrezg ids
#' @param annotation_map org.xx.eg.db object
#' @return vector with gene symbols
#' @examples
#' \dontrun{
#'    entrezg_set_to_gsymbol_set(c('1', '2'), org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @export
entrezg_set_to_gsymbol_set <- function(x, annotation_map){
   map_ids(x, from = 'ENTREZID', to = 'SYMBOL', annotation_map = annotation_map)
#   entrezg_set %>%
#         AnnotationDbi::mget(org.Hs.eg.db::org.Hs.egSYMBOL)  %>%
#         unlist()                                            %>%
#         unname()                                            %>%
#         na.exclude()
}

#' Map collapsed entrezg ids to collapsed gene symbols
#' @param x vector with collapsed entrezg ids
#' @param org.xx.eg.db respective annotation map
#' @param split collapsing character to split on
#' @return vector with collapsed gene symbols
#' @examples 
#' \dontrun{
#'    collapsed_entrezg_to_gsymbol('1;2;3', org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @export
collapsed_entrezg_to_gsymbol <- function(x, org.xx.eg.db, split=';'){
  x %>% map_collapsed_ids(from = 'ENTREZID', to = 'SYMBOL', org.xx.eg.db = org.xx.eg.db, split = split)
}

#' Map collapsed ids from one identifer format to another
#' @param x vector with collapsed entrezg ids
#' @param from A string naming the type of ID.
#' @param to A string naming the type of ID.
#' @param org.xx.eg.db respective annotation map
#' @param split collapsing character to split on
#' @return vector with collapsed gene symbols
#' @examples
#' \dontrun{
#'    map_collapsed_ids('1;2;3', from = 'ENTREZID', to = 'SYMBOL', 
#'                      org.xx.eg.db = org.Hs.eg.db::org.Hs.eg.db)
#' }
#' @export
map_collapsed_ids <- function(
  x, 
  from         = c('SYMBOL', 'ENTREZID', 'UNIPROT'), 
  to           = c('SYMBOL', 'ENTREZID', 'UNIPROT'), 
  org.xx.eg.db, 
  split        = ';'
){
  from <- match.arg(from)
  to <- match.arg(to)
  x  %>%
  strsplit(split)  %>% 
  lapply(map_ids, from = from, to = to, annotation_map = org.xx.eg.db) %>% 
  vapply(paste0, character(1), collapse=split)
}
