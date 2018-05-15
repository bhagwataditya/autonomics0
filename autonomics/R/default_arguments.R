#' Default cluster features
#' @return default value for cluster_features
#' @export
default_cluster_features <- function(){
   FALSE
}

#' Default result dir
#' @param object eset
#' @noRd
#' @importFrom autonomics.import   prepro
default_result_dir <- function(object){
   tmp <- autonomics.import::prepro(object)
   sprintf('results/%s_%ss', tmp$entity, tmp$quantity)
}
