#=====================================
# FID
#=====================================

#' Get feature id variable
#' @param object SummarizedExperiment
#' @importFrom magrittr   %<>%
#' @export
fid_var <- function(object) 'feature_id'


#' @rdname fid_var
#' @importFrom magrittr %>%
#' @export
fid_values <- function(object) object %>% fvalues('feature_id')


#================================================================
# FNAME
#================================================================

#' Get fname var/values
#' @param object   SummarizedExperiment
#' @return fvar with feature_name values
#' @examples
#' library(magrittr)
#' if (require(autonomics.data)){
#'    object <- 'extdata/stemcomp/soma/stemcomp.adat' %>%
#'               system.file(package = 'autonomics.data') %>%
#'               read_somascan()
#'    object %>% fname_var()
#'    object %>% fname_values()
#' }
#' @importFrom magrittr            %<>%
#' @export
fname_var <- function(object){
   if ('feature_name' %in% fvars(object)) 'feature_name' else 'feature_id'
}


#' @rdname fname_var
#' @importFrom magrittr %>%
#' @export
fname_values <- function(object){
   . <- NULL
   object %>% fvalues(fname_var(.))
}

#==================================================
# ENSG
#==================================================

#' Get ensg var
#' @param object SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- 'extdata/stemdiff/rnaseq/gene_counts.txt' %>%
#'               system.file(package = 'autonomics.data') %>%
#'               read_counts(fid_var = 'gene_id')
#'    fvars(object)
#' }
ensg_var <- function(object){
   fdata(object)[1:10, ] %>%
   lapply(stringi::stri_detect_fixed, 'ENSG') %>%
   vapply(any, logical(1)) %>%
   magrittr::extract(names(.), .) %>%
   magrittr::extract(1)
}


#===================================================
# UNIPROT
#===================================================

#' Get uniprot var
#' @param object  eset
#' @return fvar with uniprot accessions
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'               system.file(package = 'autonomics.data')     %>%
#'               read_proteingroups()
#'    fvars(object)
#'
#'
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% uniprot_var()
#' }
#' @importFrom magrittr            %<>%
#' @export
uniprot_var <- function(object) fvars(object) %>% magrittr::extract(stringi::stri_detect_regex(., '[uU]ni[pP]rot'))


#' Get uniprot values
#' @param object  eset
#' @param first_only logical(1)
#' @return character vector
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% uniprot_values() %>% head()
#'    object %>% uniprot_values(first_only = TRUE) %>% head()
#' }
#' @importFrom magrittr            %<>%
#' @export
uniprot_values <- function(object, first_only = FALSE){
   uniprot_var1 <- object %>% uniprot_var() %>% magrittr::extract(max(1, length(.)))
   object %>%
   fvalues(uniprot_var1) %>%
  (function (x) if (first_only) x %>% stringi::stri_split_fixed(';') %>% vapply(extract, character(1), 1) else x)
}


#=================================================================================
# ORAID
#=================================================================================

#' Get ora id var/values
#' @param object  eset
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- 'extdata/stemdiff/rnaseq/gene_counts.txt' %>%
#'               system.file(package='autonomics.data')   %>%
#'               read_counts(fid_var = 'gene_id')
#'    object %>% oraid_var()
#'    object %>% oraid_values() %>% head(3)
#' }
#' @export
oraid_var <- function(object) 'feature_id'

#' @rdname oraid_var
#' @importFrom magrittr %>%
#' @export
oraid_values <- function(object){
   object %>% fvalues(oraid_var(.))
}



#===============================
# SEP
#===============================

#' Get separator of collapsed fields
#' @param object eset
#' @return string with separator
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% sep()
#' }
#' @export
sep <- function(object){
   if      (is_maxquant_eset(object))    ';'
   else if (is_soma_eset(object))        ' '
   else if (is_rnaseq_eset(object))      NULL
   else                                  NULL
}


#======================================================
# SAMPLEID
#======================================================

#' Get sample id var
#' @param object   SummarizedExperiment
#' @return sample id var (character)
#' @importFrom magrittr           %<>%
#' @export
get_sample_id_var <- function(object) 'sample_id'

#' @rdname get_sample_id_var
#' @importFrom magrittr %>%
#' @export
get_sample_id_values <- function(object) object %>% svalues('sample_id')


