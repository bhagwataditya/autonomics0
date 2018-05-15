

#' Supported assays
#' @export
"SUPPORTED_ASSAYS"

SUPPORTED_ASSAYS     <- c('microarray', 'rnaseq',       'lcms',        'somascan', 'exiqon')

#' Supported entities
#' @export
"SUPPORTED_ENTITIES"

SUPPORTED_ENTITIES   <- c('rna',        'proteingroup', 'phosphosite', 'epitope',  'mirna')


#' Supported quantities
#' @export
"SUPPORTED_QUANTITIES"

SUPPORTED_QUANTITIES <- c('abundance', 'voomcount', 'raw.ratio', 'normalized.ratio', 'occupancy',
                          'reporter.intensity', 'lfq.intensity', 'raw.intensity', 'ct')

#' Supported softwares
#' @export
"SUPPORTED_SOFTWARES"

SUPPORTED_SOFTWARES  <- c('maxquant', 'limma', 'somalogic', 'affymetrix', 'genex')


#' Create list with preprocessing details
#' @param assay    value in SUPPORTED_ASSAYS
#' @param entity   value in SUPPORTED_ENTITIES
#' @param quantity value in SUPPORTED_QUANTITIES
#' @param software value in SUPPORTED_SOFTWARES
#' @param parameters preprocessing parameters
#' @return list(assay, entity, quantity, software, parameters)
#' @export
create_prepro_list <- function(assay, entity, quantity, software, parameters = list()){
   
   # Assert
   assertive.sets::assert_is_subset(assay,      SUPPORTED_ASSAYS)
   assertive.sets::assert_is_subset(entity,     SUPPORTED_ENTITIES)
   assertive.sets::assert_is_subset(quantity,   SUPPORTED_QUANTITIES)
   assertive.sets::assert_is_subset(software,   SUPPORTED_SOFTWARES)
   assertive.types::assert_is_list(parameters)
   
   list(
      assay       = assay,
      entity      = entity,
      quantity    = quantity,
      software    = software,
      parameters  = parameters
   )
}



