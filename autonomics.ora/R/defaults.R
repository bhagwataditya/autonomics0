#' Default value for argument detectome_against_genome
#' 
#' Should ora of detectome versus genome be performed?
#' @return default value of detectome_against_genome
#' @export
default_ora_detectome_in_genome <- function(){
   FALSE
}

#' default universe
#' @param object eset
#' @return default value of universe
#' @export
default_universe <- function(object){
   if (autonomics.import::prepro(object)$entity == 'metabolite') NULL
   else                                                         'detectome'
}