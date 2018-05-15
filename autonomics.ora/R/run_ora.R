utils::globalVariables('.')

#' Uncollapse and unite
#' @param x collapsed string
#' @param sep separator to strsplit on
#' @importFrom magrittr %>%
#' @export
uncollapse_and_unite <- function(x, sep = ';'){
  if (length(x) == 0){
    return(character(0))
  } else {
    x %>% strsplit(sep) %>% unlist() %>% unique()
  }
}

#' Run over representation analysis
#' 
#' Each pathway is analyzed for an enrichment in query genes, 
#' using a right-sided (over representation) Fisher Exact test.
#' 
#' First define:
#'    q = # query features in pathway
#'    m = # pathway genes
#'    n = # non pathway genes
#'    k = # query genes
#' 
#' The right-sided Fisher exact test investigates whether:
#'    H0: q = random draw from universe
#'    H1: q > random draw from universe
#' 
#' Under H0:
#'    q ~ phyper(m, n, k)
#'    
#' So the Fisher exact test examines whether 
#'    phyper(q, m, n, k, lower.tail = FALSE) <= alpha
#'    
#' @note 1. My original implementation used 'fisher.test', but that was not vectorizable.
#   2. I vectorized this using 'phyper', as used in the R package EnrichmentBrowser.
#      However, that approach gave different values than fisher.test.
#   3. Meng's blog gave me more insight into matter, explaining:
#         - how to properly use phyper
#         - why to decrement q (with one) for over representation, in order to get P(X>=x) rather than P(X>x).
#         - how to use phyper(..., lower.tail = FALSE)
#      \url{http://mengnote.blogspot.qa/2012/12/calculate-correct-hypergeometric-p.html}
#   4. My simulations showed that phyper(lower.tail = FALSE) and 1 - phyper(lower.tail = TRUE) 
#      give different results when the enrichment is strong but not absolute.
#         - phyper(lower.tail=TRUE): small, non-zero values, identical to those with fisher.test
#         - 1 - phyper: returns 0 for near-0 p values 
#   5. I switched to phyper(lower.tail = TRUE) for two reasons:
#         - I wish to differentiate among the strong enrichments
#         - I wish for the p values to be identical as when obtained by fisher.test.
#' @param query Character vector of feature IDs in the selection.
#' @param universe Character vector of feature IDs in the whole dataset (all 
#' detected features).
#' @param pathway_list A named list of character vectors. Names correspond to 
#' pathways. Each component of the list is a character vector of feature IDs.
#' @param min_set_size An integer. The minimum number of universe genes 
#' required for a pathway to be included in analysis.
#' @param return_only_significant A logical value. Should returned results be 
#' limited to those that are significant?
#' @param collapse_features_column Logical.  Should the \code{features} column 
#' be collapsed into a character vector?  See return value section.
#' @return A data frame with the following columns:
#' \describe{
#'   \item{rank}{Integer, 1 to the number of rows.  The rank of the significance 
#'   of the pathway.}
#'   \item{p}{Numeric, between 0 and 1. The p-value scoring the significance of
#'   the pathway.}
#'   \item{pathway}{Character. The ID of the pathway, taken from 
#'   \code{names(pathway_list)}.}
#'   \item{n_total}{Integer, positive. The number of features found in the 
#'   \code{pathway_list} with this pathway.}
#'   \item{n_detected}{Integer, positive. The number of features found in the 
#'   \code{universe} with this pathway.}
#'   \item{n_selected}{Integer, positive. The number of features found in the 
#'   \code{query} with this pathway.}
#'   \item{features}{Character or list of character.  If 
#'   \code{collapse_features_column} is \code{TRUE}, a character vector of semi-
#'   colon separated feature names.  Otherwise, a list of character vectors of 
#'   feature names.}
#' }
#' @importFrom magrittr   %<>%   %>%
#' @importFrom tibble     data_frame
#' @export
run_ora <- function(query, universe, pathway_list, min_set_size = 5, return_only_significant = FALSE, collapse_features_column = TRUE){
  
  # Remove duplicates
  query <- dedupe(query)
  universe <- dedupe(universe)
  pathway_list <- pathway_list %>% lapply(dedupe)
  
  n_total_features_per_pathway    <- pathway_list %>% lengths
  
  # Limit to universe of detected features
  features_in_universe_per_pathway <- pathway_list %>% lapply(intersect, universe)
  n_detected_features_per_pathway <- features_in_universe_per_pathway %>% lengths
  pathway_has_enough_detected_features <- n_detected_features_per_pathway >= min_set_size
  n_total_features_per_pathway     %<>% magrittr::extract(pathway_has_enough_detected_features)
  n_detected_features_per_pathway  %<>% magrittr::extract(pathway_has_enough_detected_features)
  features_in_universe_per_pathway %<>% magrittr::extract(pathway_has_enough_detected_features)
  
  # Collect q,m,n,k for hypergeometric test
  features_in_query_per_pathway <- features_in_universe_per_pathway %>% lapply(intersect, query)
  n_selected_features_per_pathway <- features_in_query_per_pathway %>% lengths
  n_selected_features                <- length(query)
  n_non_pathway_features_per_pathway <- length(universe) - n_detected_features_per_pathway
  
  # Perform hypergeometric test
  # See https://mengnote.blogspot.qa/2012/12/calculate-correct-hypergeometric-p.html
  # Why is there a -1 in the q argument?
  # > Because if lower.tail is TRUE (default), probabilities are P[X ≤ x], 
  # > otherwise, P[X > x]. We subtract x by 1, when P[X ≥ x] is needed.
  pvalues    <- stats::phyper(q = n_selected_features_per_pathway - 1,
                       m = n_detected_features_per_pathway, 
                       n = n_non_pathway_features_per_pathway, 
                       k = n_selected_features, 
                       lower.tail = FALSE)
  
  if(collapse_features_column)
  {
    features_in_query_per_pathway <- features_in_query_per_pathway %>% 
      vapply(paste, character(1), collapse = ';')
  }
  
  # Record results
  res_df <- tibble::data_frame(rank       = NA_integer_, # will be rank(pvalues), but saves computation to add it later
                       p          = pvalues,
                       pathway    = names(pvalues),
                       n_total    = n_total_features_per_pathway,
                       n_detected = n_detected_features_per_pathway,
                       n_selected = n_selected_features_per_pathway,
                       features   = features_in_query_per_pathway
  )
  res_df %<>% dplyr::arrange_(~ p)
  res_df$rank <- seq_len(nrow(res_df))
  
  # Limit to significant
  if (return_only_significant){
    res_df %<>% dplyr::filter_(~ p < 0.05)
  }
  
  # Return
  data.table::data.table(res_df)
}

