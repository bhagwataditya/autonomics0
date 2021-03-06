TARGETSCAN_ORGANISMS <- c('H.sapiens', 'M.musculus', 'D.rerio')

SPECIES_TO_TAXONID <- c(
   D.rerio        = 7955,
   X.tropicalis   = 8364,
   G.gallus       = 9031,
   M.mulatta      = 9544,
   P.troglodytes  = 9598,
   H.sapiens      = 9606,
   C.familiaris   = 9615,
   B.taurus       = 9913,
   M.musculus     = 10090,
   R.norvegicus   = 10116,
   M.domestica    = 13616)

get_targetscan_cachefile <- function(organism){
   sprintf('~/.autonomics/targetscan/%s/target_predictions.txt', organism)
}

#' @importFrom   data.table   data.table
#' @importFrom   magrittr     %<>%   %>%
download_targetscan_fish <- function(cachefile){

   # Satisfy CHECK
   . <- NULL

   cachedir <- dirname(cachefile)
   remote <- 'http://www.targetscan.org//fish_62/fish_62_data_download/Summary_Counts.txt.zip'
   local <- paste0(cachedir, '/', basename(remote))

   utils::download.file(remote, local)
   utils::unzip(local, exdir = cachedir)
   local %<>% substr(1, nchar(.)-4)

   data.table::fread(local, select = c( 'Representative miRNA', 'Transcript ID', 'Gene Symbol', 'Species ID')) %>%
   data.table::setnames(c('Representative miRNA', 'Transcript ID', 'Gene Symbol', 'Species ID'),
                        c('mir',                   'ensg',         'gsymbol',     'taxonid'))
}

#' @importFrom  data.table   data.table   :=
#' @importFrom  magrittr     %>%          %<>%
download_targetscan_human_mouse <- function(organism = 'H.sapiens', cachefile){

   # Satisfy CHECK
   . <- `Transcript ID` <- NULL

   cachedir <- dirname(cachefile)

   # Target predictions
   subdir <- switch(organism, H.sapiens = 'vert_72/vert_72_data_download', M.musculus = 'mmu_72/mmu_72_data_download')
   remote <- sprintf('http://www.targetscan.org/%s/Summary_Counts.default_predictions.txt.zip', subdir)
   local  <-  paste0(cachedir, '/', basename(remote))
   utils::download.file(remote, local)
   utils::unzip(local, exdir = cachedir)
   local %<>% substr(1, nchar(.)-4)
   target_predictions <- data.table::fread(local, select = c('Representative miRNA', 'Transcript ID', 'Species ID' ))

   # Gene annotations
   remote <- sprintf('http://www.targetscan.org/%s/Gene_info.txt.zip', subdir)
   local  <- paste0(cachedir, '/', basename(remote))
   utils::download.file(remote, local)
   utils::unzip(local, exdir = cachedir)
   local %<>% substr(1, nchar(.)-4)
   gene_info <- data.table::fread(local, select = c('Transcript ID', 'Gene ID', 'Gene symbol'))

   # Merge
   target_predictions %<>% merge(gene_info, by = 'Transcript ID')
   target_predictions[, `Transcript ID` := NULL]
   target_predictions %>% data.table::setnames(c('Representative miRNA', 'Species ID', 'Gene ID', 'Gene symbol'),
                                               c('mir',                  'taxonid',    'ensg',    'gsymbol'))
}

#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%   %<>%
download_targetscan <- function(organism){

   # Abort if file already in cache
   cachefile <- get_targetscan_cachefile(organism)
   if(file.exists(cachefile)) return(cachefile)

   # Satisfy CHECK
   taxonid <- mir <- ensg <- . <- NULL

   # Assert validity
   assertive.sets::assert_is_subset(organism, c('H.sapiens', 'M.musculus', 'D.rerio'))

   # Create cache dir
   cachefile <- get_targetscan_cachefile(organism)
   cachedir <- dirname(cachefile)
   dir.create(cachedir, showWarnings = FALSE, recursive = TRUE)

   # Organism specific parts
   target_predictions <- switch(organism,
                                H.sapiens  = download_targetscan_human_mouse(organism, cachefile),
                                M.musculus = download_targetscan_human_mouse(organism, cachefile),
                                D.rerio    = download_targetscan_fish(cachefile))

   # Rename
   target_predictions[, mir := tolower(mir)]

   # Number
   target_predictions[, number := mir %>% substr(9, nchar(.))]

   # Filter
   mir_taxon_id = SPECIES_TO_TAXONID %>% magrittr::extract2(organism)
   target_predictions %<>% magrittr::extract(taxonid %in% mir_taxon_id)
   target_predictions %>%  magrittr::extract(, taxonid:=NULL)

   # Arrange
   target_predictions %<>% magrittr::extract(order(mir))

   # Trim ensgs (remove version suffix)
   target_predictions[, ensg := ensg %>% substr(1, stringi::stri_locate_first_fixed(., '.')[, 1] - 1)]

   # Save
   target_predictions %>% autonomics.support::print2txt(file = paste0(cachedir, '/target_predictions.txt'))

}




#' Load targetscan
#' @param organism Any value in TARGETSCAN_ORGANISMS
#' @examples
#' \dontrun{
#'    autonomics.annotate::load_targetscan('H.sapiens')
#'    autonomics.annotate::load_targetscan('M.musculus')
#'    autonomics.annotate::load_targetscan('D.rerio')
#' }
#' @export
load_targetscan <- function(organism){

   # Assert
   assertive.sets::assert_is_subset(organism, TARGETSCAN_ORGANISMS)

   # Download (only if required)
   download_targetscan(organism)

   # Read
   data.table::fread(get_targetscan_cachefile(organism))
}

#' Infer organism
#' @param mir vector of micrRNAs
#' @return name of infered organism (value in TARGETSCAN_ORGANISMS)
#' @examples
#' require(magrittr)
#' mir <- c("hsa-mir-7-5p", "hsa-mir-29a-3p", "hsa-mir-199a-5p ")
#' mir %>% infer_organism_from_mirs()
#' mir <- 'rna-blablabla'
#' infer_organism_from_mirs(mir)
#' @importFrom data.table   data.table
#' @importFrom magrittr     %>%
#' @export
infer_organism_from_mirs <- function(mir){
   # Satisfy CHECK
   J <- NULL
   mir %<>% tolower()

   TARGETSCAN_ORGANISMS %>%
   vapply(function(x){
           targetscan <- load_targetscan(x)
           targetscan %>% data.table::setkey('mir')
           # https://stackoverflow.com/questions/17331684/fast-exists-in-data-table
           query_mir <- mir
           .subset2(targetscan[J(query_mir), mult = 'first', nomatch=0], "mir") %>%
              length() %>%
              magrittr::is_greater_than(0)
        }, logical(1)) %>%
   (function(x){
     if (sum(x)==0) autonomics.support::cmessage('microRNAs do not match targetscan data of any organism - abort')
     if (sum(x) >1) autonomics.support::cmessage('microRNAs match targetscan data of multiple organisms - abort')
     names(x)[x]
   })
}

#' Map mir to genes
#'
#' Map vector of microRNAs to target genes
#' @param mir      vector of microRNAs
#' @param to       'ensg' or 'gsymbol'
#' @return named list that maps each mir to its target genes.
#' @examples
#' require(magrittr)
#' mir <- c('hsa-mir-7-5p', 'hsa-mir-29a-3p', 'hsa-mir-199a-5p ')
#' mir %>% mir_to('gsymbol')
#' @importFrom magrittr %>%
#' @export
mir_to <- function(mir, to){
   organism <- infer_organism_from_mirs(mir)
   targetscan <- autonomics.annotate::load_targetscan(organism)
   targetscan %>% data.table::setkey(mir)
   targetscan %<>% magrittr::extract(mir)
   targetscan %<>% split(by = 'mir', keep.by = FALSE) %>% lapply(magrittr::extract2, to)
   if (length(targetscan)==1)   targetscan %<>% magrittr::extract2(1)
   targetscan
}

#' Map vector of microRNAs to target gsymbols
#' @param mir vector of microRNAs
#' @return named list which maps each mir to its target gene symbols
#' @examples
#' require(magrittr)
#' mir <- c('hsa-mir-7-5p', 'hsa-mir-29a-3p', 'hsa-mir-199a-5p ')
#' mir %>% mir_to_gsymbol()
#' @export
mir_to_gsymbol <- function(mir){
   mir_to(mir, 'gsymbol')
}

#' Map vector of microRNAs to target ensgs
#' @param mir vector of microRNAs
#' @return named list which maps each mir to its target ensgs
#' @examples
#' require(magrittr)
#' mir <- c('hsa-mir-7-5p', 'hsa-mir-29a-3p', 'hsa-mir-199a-5p ')
#' mir %>% mir_to_ensg()
#' @importFrom magrittr %>%
#' @export
mir_to_ensg <- function(mir){
   mir_to(mir, 'ensg')
}

