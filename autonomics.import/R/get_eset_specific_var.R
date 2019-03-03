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
fid_values <- function(object) object %>% autonomics.import::fvalues('feature_id')


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
#'               autonomics.import::read_somascan()
#'    object %>% autonomics.import::fname_var()
#'    object %>% autonomics.import::fname_values()
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
   object %>% autonomics.import::fvalues(autonomics.import::fname_var(.))
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
#'    autonomics.import::fvars(object)
#'
#'
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.import::uniprot_var()
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
#'    object %>% autonomics.import::uniprot_values() %>% head()
#'    object %>% autonomics.import::uniprot_values(first_only = TRUE) %>% head()
#' }
#' @importFrom magrittr            %<>%
#' @export
uniprot_values <- function(object, first_only = FALSE){
   object %>%
   autonomics.import::fvalues(autonomics.import::uniprot_var(.)) %>%
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
#'    object <- autonomics.data::billing2016
#'    object %>% autonomics.import::oraid_var()
#'    object %>% autonomics.import::oraid_values() %>% head(3)
#' }
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% autonomics.import::oraid_var()
#'    object %>% autonomics.import::oraid_values() %>% head(3)
#'
#'    object <- billing.differentiation.data::protein.ratios
#'    object %>% autonomics.import::oraid_var()
#'    object %>% autonomics.import::oraid_values()
#' }
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma
#'    object %>% autonomics.import::oraid_var()
#'    object %>% autonomics.import::oraid_values() %>% head(3)
#' }
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.import::oraid_var()
#'    object %>% autonomics.import::oraid_values() %>% head(3)
#'
#'    object <- subramanian.2016::exiqon
#'    object %>% autonomics.import::oraid_var()
#'    object %>% autonomics.import::oraid_values() %>% head(3)
#' }
#' @export
oraid_var <- function(object){

   # Assert
   autonomics.import::assert_is_valid_object(object)

   if (is_rnaseq_eset(object))           return(ensg_var(object))
   if (is_soma_eset(object))             return(uniprot_var(object))
   if (is_maxquant_eset(object))         return(uniprot_var(object))
   if (is_metabolon_eset(object))        return('BIOCHEMICAL')
   if (is_exiqon_eset(object))           return('feature_id')
   autonomics.support::cmessage('Abort - object must be a SummarizedExperiment with a format as created by one of the autonomics importers')
   return(NULL)
}

#' @rdname oraid_var
#' @importFrom magrittr %>%
#' @export
oraid_values <- function(object){
   object %>% autonomics.import::fvalues(autonomics.import::oraid_var(.))
}



#===============================
# SEP
#===============================

#' Get separator of collapsed fields
#' @param object eset
#' @return string with separator
#' @examples
#' if (require(billing.differentiation.data)){
#'    require(magrittr)
#'    object <- rna.voomcounts
#'    object %>% autonomics.import::sep()
#' }
#' @export
sep <- function(object){
   if      (is_maxquant_eset(object))    ';'
   else if (is_soma_eset(object))        ' '
   else if (is_rnaseq_eset(object))      NULL
   else                                  NULL
}

#' @rdname sep
#' @export
get_sep <- function(object){
   .Deprecated('sep')
   sep(object)
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
get_sample_id_values <- function(object) object %>% autonomics.import::svalues('sample_id')


