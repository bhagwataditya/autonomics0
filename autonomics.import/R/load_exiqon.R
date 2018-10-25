#' Load exiqon data
#' @param file                         exiqon xlsx file
#' @param design_file                  NULL or character (sample design file)
#' @param infer_design_from_sampleids  logical: whether to infer design from sample ids
#' @param design_sep                   string: sample id separator
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'
#'    # Load exiqon and infer design from sampleids
#'    #--------------------------------------------
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                         package = 'subramanian.2016')
#'    file %>% autonomics.import::load_exiqon(infer_design_from_sampleids = TRUE)
#'
#'    # Load exiqon with design file
#'    #-----------------------------
#'    design_file <- tempfile()
#'    file %>% write_design('exiqon', infer_design_from_sampleids = TRUE,
#'                                    design_file = design_file)
#'    file %>% load_exiqon(design_file = design_file) %>%
#'             autonomics.import::sdata() %>% magrittr::extract(1:3, 1:5)
#' }
#' @importFrom magrittr %>%
#' @export
load_exiqon <- function(
   file,
   design_file                 = NULL,
   infer_design_from_sampleids = FALSE,
   design_sep                  = NULL
){
   . <- `#RefGenes` <- `#Spike` <- NULL
   object <- autonomics.import::load_omics(file                        = file,
                                           platform                    = 'exiqon',
                                           log2_transform              = FALSE,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep)
   autonomics.import::prepro(object) <- list(assay = 'exiqon', entity = 'mirna', quantity = 'ct', software = 'genex')
   object
}

#==========================================================================

#' Preprocess Exiqon Ct
#' @param object                    SummarizedExperiment
#' @param filter_features           character(1) with feature filter condition
#' @param align_sample_means        logical: whether to align the sample maens
#' @param lod                       numeric(1): lod Ct value
#' @param zero_consistent_nas       logical(1)
#' @param filter_conserved_in_human logical: whether to filter for features consvered in human (useful for non-human microRNA profiles)
#' @param plot                      logical: whether to plot sample distributions
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- 'extdata/exiqon/subramanian.2016.exiqon.xlsx' %>%
#'               system.file(package = 'subramanian.2016') %>%
#'               autonomics.import::load_exiqon(infer_design_from_sampleids = TRUE)
#'    object %>% autonomics.import::prepro_exiqon(lod=36)
#' }
#' @importFrom magrittr %>%
#' @export
prepro_exiqon <- function(
   object,
   filter_features             = '`#RefGenes`==0 & `#Spike`   ==0',
   align_sample_means          = TRUE,
   lod                         = 37,
   zero_consistent_nas         = TRUE,
   filter_conserved_in_human   = FALSE,
   plot                        = TRUE
){
   # Plot function
   #--------------
   plotfun <- function(object, xlab, title){
      object %>%
         autonomics.plot::plot_overlayed_sample_distributions() +
         ggplot2::xlab(xlab) +
         ggplot2::ggtitle(title)
   }

   # Original
   #---------
   if (plot) object %>% plotfun('Ct', 'Load Ct values') %>% print()

   # Filter features
   #----------------
   object %<>% autonomics.import::filter_features_(filter_features, verbose = TRUE)

   # Align sample means
   #-------------------
   if (align_sample_means){
      autonomics.support::cmessage('\t\tAlign sample means')
      sample_means <- object %>% autonomics.import::exprs() %>% (function(y){y[y>32] <- NA; y}) %>% colMeans(na.rm = TRUE)
      sample_diffs <- sample_means - median(sample_means)
      autonomics.import::exprs(object) %<>% sweep(2, sample_diffs)
      if (plot) object %>% plotfun('Ct', 'Align sample means') %>% print()
   }

   # NA beyond lod
   #-----------------------------
   object %<>% autonomics.preprocess::na_exprs_weakly_gt(lod)
   if (plot) object %>% plotfun('Ct', 'NA Ct >= lod') %>% print()

   # Invert scale
   #-------------
   autonomics.support::cmessage('\t\tInvert exprs = %d - Ct', lod)
   autonomics.import::exprs(object) %<>% magrittr::subtract(lod, .)
   newmetric <- sprintf('%d - Ct', lod)
   if (plot) object %>% plotfun(newmetric, newmetric) %>% print()

   # Zero consistent NAs
   #--------------------
   if (zero_consistent_nas){
      object %<>% autonomics.preprocess::zero_consistent_nas(verbose = TRUE)
      if (plot) object %>% plotfun(newmetric, 'Zero consistent NA values') %>% print()
   }

   # Filter for mirs conserved in human
   #-----------------------------------
   if (filter_conserved_in_human){
      autonomics.import::fdata(object) %<>% (function(x){x$number <- x$feature_id %>% substr(9, nchar(.)); x})
      conserved_human_mirs <- autonomics.annotate::load_targetscan('H.sapiens') %>% magrittr::extract2('number')
      object %<>% autonomics.import::filter_features(number %in% conserved_human_mirs, verbose = TRUE)
   }

   # Return
   object
}
