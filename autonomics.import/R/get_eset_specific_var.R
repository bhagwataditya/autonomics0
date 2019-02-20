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
#' @param object   eset
#' @return fvar with gene symbols
#' @examples
#' library(magrittr)
#' if (require(atkin.2014)){
#'    atkin.2014::soma %>% autonomics.import::fname_var()
#'    atkin.2014::soma %>% autonomics.import::fname_values()
#' }
#' @importFrom magrittr            %<>%
#' @export
fname_var <- function(object){

   # Keep default for esets without prepro
   fvar_name <- character(0)

   # Dispatch on eset type
   if (is_rnaseq_eset(object))    fvar_name <- 'gene_name'           # rnaseq
   if (is_exiqon_eset(object))    fvar_name <- 'feature_id'          # exiqon
   if (is_maxquant_eset(object))  fvar_name <- 'Gene names'          # max quant
   if (is_soma_eset(object))      fvar_name <- 'EntrezGeneSymbol'    # somascan
   if (is_metabolon_eset(object)) fvar_name <- 'BIOCHEMICAL'         # metabolon

   # Return
   return(fvar_name)
}

#' @rdname fname_var
#' @importFrom magrittr %>%
#' @export
fname_values <- function(object){
   . <- NULL
   object %>% autonomics.import::fvalues(autonomics.import::fname_var(.))
}

#' @rdname fname_var
#' @export
get_gene_symbol_var <- function(object){
   .Deprecated('fname_var')
   fname_var(object)
}


#' @rdname fname_var
#' @export
get_gene_symbol_values <- function(object){
   .Deprecated('fname_values')
   fname_values(object)
}

#================================================
# ENSG
#================================================

#' Get ensg var/values
#' @param object  eset
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% autonomics.import::ensg_var()
#'    object %>% autonomics.import::ensg_values()
#' }
#' @export
ensg_var <- function(object){
   if      (is_rnaseq_eset(object))      'gene_id'
   else if (is_soma_eset(object))        stop('no ensg in soma esets')
   else if (is_maxquant_eset(object))    stop('no ensg in maxquant esets')
   else                                  stop('eset of unknown type')
}


#' @rdname ensg_var
#' @importFrom magrittr %>%
#' @export
ensg_values <- function(object){
   . <- NULL
   object %>% autonomics.import::fvalues(autonomics.import::ensg_var(.))
}

#' @rdname ensg_var
#' @export
get_ensg_var <- function(object){
   .Deprecated('ensg_var')
   ensg_var(object)
}

#' @rdname ensg_var
#' @export
get_ensg_values <- function(object){
   .Deprecated('ensg_values')
   ensg_values(object)
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
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.import::uniprot_var()
#' }
#' @importFrom magrittr            %<>%
#' @export
uniprot_var <- function(object){

   # Keep default for esets without prepro
   fvar_name <- 'uniprot_accessions'

   if (is_maxquant_eset(object))    fvar_name <- 'Uniprot accessions'   # max quant
   if (is_soma_eset(object))        fvar_name <- 'UniProt'              # somascan

   assertive.sets::is_subset(fvar_name, autonomics.import::fvars(object))
   fvalues <- autonomics.import::fdata(object)[[fvar_name]]
   if (is.factor(fvalues))   fvalues %<>% as.character()
   assertive.strings::assert_any_are_non_empty_character(fvalues)
   return(fvar_name)
}

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
   autonomics.import::assert_is_valid_eset(object)

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


