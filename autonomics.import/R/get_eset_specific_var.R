#=====================================
# FID
#=====================================

#' Get fid var/values
#' @param object eset
#' @return feature id variable / values
#' @examples
#' library(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::ALL %>% autonomics.import::fid_var()
#'    autonomics.data::ALL %>% autonomics.import::fid_values() %>% head(1)
#' }
#' if (require(atkin.2014)){
#'    atkin.2014::soma %>% autonomics.import::fid_var()
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>% autonomics.import::fid_var()
#'    billing.differentiation.data::protein.ratios %>% autonomics.import::fid_var()
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% autonomics.import::fid_var()
#'    subramanian.2016::exiqon    %>% autonomics.import::fid_var()
#'    subramanian.2016::rnaseq    %>% autonomics.import::fid_var()
#'
#'    subramanian.2016::metabolon %>% autonomics.import::fid_values() %>% head(3)
#'    subramanian.2016::exiqon    %>% autonomics.import::fid_values() %>% head(3)
#'    subramanian.2016::rnaseq    %>% autonomics.import::fid_values() %>% head(3)
#' }
#' @importFrom magrittr   %<>%
#' @export
fid_var <- function(object){

   # Keep default for esets without prepro
   fvar_name <- 'feature_id'
   if (autonomics.import::is_rnaseq_eset(object))    fvar_name <- 'gene_id'     # rnaseq
   if (autonomics.import::is_exiqon_eset(object))    fvar_name <- 'feature_id'  # exiqon
   if (autonomics.import::is_maxquant_eset(object))  fvar_name <- 'feature_id'  # max quant
   if (autonomics.import::is_soma_eset(object))      fvar_name <- 'SeqId'       # somascan
   if (autonomics.import::is_metabolon_eset(object)) fvar_name <- 'MCOMP_ID'    # metabolon:  M + COMPOUNDID

   assertive.sets::is_subset(fvar_name, autonomics.import::fvars(object))
   fvalues <- autonomics.import::fdata(object)[[fvar_name]]
   if (is.factor(fvalues))   fvalues %<>% as.character()
   assertive.strings::assert_any_are_non_empty_character(fvalues)    # feature ids MUST be present
   return(fvar_name)
}

#' @rdname fid_var
#' @importFrom magrittr %>%
#' @export
fid_values <- function(object){
   . <- NULL
   object %>% autonomics.import::fvalues(autonomics.import::fid_var(.))
}

#' @rdname fid_var
#' @export
get_feature_id_var <- function(object){
   .Deprecated('fid_var')
   fid_var(object)
}

#' @rdname fid_var
#' @export
get_feature_id_values <- function(object){
   .Deprecated('fid_values')
   fid_values(object)
}


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

#' @rdname uniprot_var
#' @export
get_uniprot_var <- function(object){
   .Deprecated('uniprot_var')
   uniprot_var(object)
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
#' @param object   eset
#' @return sample id var (character)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::billing2016 %>% get_sample_id_var()
#'    autonomics.data::billing2016 %>% get_sample_id_values()
#' }
#' if (require(atkin.2014)){
#'    atkin.2014::soma %>% autonomics.import::get_sample_id_var()
#'    atkin.2014::soma %>% autonomics.import::get_sample_id_values()
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::rnaseq %>% autonomics.import::get_sample_id_var()
#'    subramanian.2016::rnaseq %>% autonomics.import::get_sample_id_values()
#' }
#' @importFrom magrittr           %<>%
#' @export
get_sample_id_var <- function(object){
   sample_id_var <- 'sample_id'           # metabolon
   if      (is_maxquant_eset(object)){  sample_id_var <- 'sample_id'
   } else if (is_soma_eset(object)){    sample_id_var <- 'SampleId'
   } else if (is_rnaseq_eset(object)){  sample_id_var <- 'sample_id'
   }
   assertive.sets::assert_is_subset(sample_id_var, autonomics.import::svars(object))
   sample_id_values <- autonomics.import::sdata(object)[[sample_id_var]]
   if (is.factor(sample_id_values))   sample_id_values %<>% as.character()
   assertive.strings::assert_all_are_non_empty_character(sample_id_values)
   return(sample_id_var)
}

#' @rdname get_sample_id_var
#' @importFrom magrittr %>%
#' @export
get_sample_id_values <- function(object){
   . <- NULL
   object %>% autonomics.import::svalues(autonomics.import::get_sample_id_var(.))
}

