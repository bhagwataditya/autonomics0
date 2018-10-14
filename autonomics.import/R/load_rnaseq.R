#' Load RNAseq cpm
#' @param file                          rnaseq counts file
#' @param design_file                   sample design file
#' @param infer_design_from_sampleids   logical
#' @param design_sep                    character
#' @examples
#' require(magrittr)
#' if (require(autonomics.data){
#'    file <- system.file('extdata/stemdiff/rnaseq/gene_counts.txt', package = 'autonomics.data')
#' })
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    load_rnaseq(file, infer_design_from_sampleids = TRUE)
#' }
#' @export
load_rnaseq <- function(
   file,
   design_file                 = NULL,
   infer_design_from_sampleids = FALSE,
   design_sep                  = NULL
){

   # Load sumexp
   autonomics.support::cmessage('\t\tload counts')
   object <- autonomics.import::load_omics(file                        = file,
                                           platform                    = 'rnaseq',
                                           log2_transform              = FALSE,
                                           design_file                 = design_file,
                                           infer_design_from_sampleids = infer_design_from_sampleids,
                                           design_sep                  = design_sep)
   object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample()

   # Add metadata
   autonomics.import::prepro(object) <- list(assay      = 'rnaseq',
                                             entity     = 'rna',
                                             quantity   = 'log2cpm',
                                             software   = 'limma',
                                             parameters = list())
   object
}


#' @rdname counts_to_cpm
#' @export
libsizes <- function(counts){
   colSums(counts) * edgeR::calcNormFactors(counts)
}

#' Convert between counts and cpm (counts per million scaled reads)
#' @param counts    count matrix
#' @param lib.size  scaled library sizes (vector)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    counts <- 'extdata/rnaseq/gene_counts.txt' %>%
#'               system.file(package = 'subramanian.2016') %>%
#'               autonomics.import::load_rnaseq(infer_design_from_sampleids = TRUE) %>%
#'               autonomics.import::exprs()
#'    lib.size <- counts %>% autonomics.import::libsizes()
#'    cpm      <- counts %>% autonomics.import::counts_to_cpm(lib.size)
#'    counts2  <- cpm %>% autonomics.import::cpm_to_counts(lib.size)
#'    sum(counts-counts2)
#' }
#' @return cpm matrix
#' @export
counts_to_cpm <- function(counts, lib.size = libsizes(counts)){
   t(t(counts + 0.5)/(lib.size + 1) * 1e+06)
}

#' @rdname counts_to_cpm
#' @export
cpm_to_counts <- function(cpm, lib.size){
   1e-06 * t(t(cpm) * (lib.size + 1)) - 0.5
}

#' @importFrom    magrittr %>%
compute_precision_weights_once <- function(
   object,
   design = object %>% autonomics.import::create_design_matrix(),
   plot   = TRUE,
   ...
){
   # Compute
   counts   <- object %>% autonomics.import::counts()
   lib.size <- counts %>% autonomics.import::libsizes()
   counts %>% limma::voom(design = design, lib.size = lib.size, plot = plot, ...) %>%
              magrittr::extract2('weights') %>%
              magrittr::set_rownames(rownames(object)) %>%
              magrittr::set_colnames(colnames(object))
}


#' Compute voom precision weights
#' @param object  SummarizedExperiment: exprs(.) with log2cpm, counts(.) with raw counts.
#' @param design  design matrix
#' @param plot    logical
#' @param ...     passed to limma::voom() -> limma::lmFit()
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- 'extdata/rnaseq/gene_counts.txt' %>%
#'               system.file(package = 'subramanian.2016') %>%
#'               autonomics.import::load_rnaseq(infer_design_from_sampleids = TRUE)
#'    object %>% compute_precision_weights()
#' }
#' @importFrom magrittr %>%
#' @export
compute_precision_weights <- function(
   object,
   design = autonomics.import::create_design_matrix(object),
   plot   = TRUE
){

   # Assert
   assertive.properties::assert_is_not_null(autonomics.import::counts(object))

   # Estimate precision weights
   has_block <- autonomics.import::has_complete_block_values(object)
   weights <- object %>% compute_precision_weights_once(design = design, plot = !has_block)

   # Update precision weights using block correlation
   if (has_block){
      correlation <- autonomics.import::counts(object) %>% counts_to_cpm() %>% log2() %>%
                     limma::duplicateCorrelation(design = design, block = object$block, weights = weights) %>%
                     magrittr::extract2('consensus')
      weights <- object %>% compute_precision_weights_once(design      = design,
                                                           block       = object$block,
                                                           correlation = correlation,
                                                           plot        = TRUE)
   }

   # Return
   dimnames(weights) <- dimnames(object)
   weights
}

#' Preprocess RNAseq counts
#' @param object SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- 'extdata/rnaseq/gene_counts.txt' %>%
#'               system.file(package = 'subramanian.2016') %>%
#'               autonomics.import::load_rnaseq(infer_design_from_sampleids = TRUE)
#'    object %<>% prepro_rnaseq()
#'    object %>% autonomics.plot::plot_overlayed_sample_distributions()
#' }
#' @importFrom magrittr %>%
#' @export
prepro_rnaseq <- function(object, plot = TRUE){

   # Store counts
   autonomics.support::cmessage('\t\tStore counts in autonomics.import::counts(object)')
   autonomics.import::counts(object) <- autonomics.import::exprs(object)

   # TMM-normalize
   message('\t\tTMM normalize exprs: counts -> cpm')
   autonomics.import::exprs(object) %<>% autonomics.import::counts_to_cpm()

   # Filter
   object %<>% autonomics.import::filter_exprs_replicated_in_some_subgroup('>', 1)

   # Log2 transform
   message('\t\tLog2 transform exprs: cpm -> log2cpm')
   autonomics.import::exprs(object) %<>% log2()

   # Add precision weights
   autonomics.support::cmessage('\t\tCompute and add precision weights')
   SummarizedExperiment::assays(object)$weights <- object %>% autonomics.import::compute_precision_weights(plot = plot)

   # Return
   object

   # Quantile normalize
   # Gordon:  prefer TMM over quantile normalization for most cases
   #          Use quantile normalization for extreme cases
   #          Don't use both
   # object %<>% limma::normalizeBetweenArrays(method = normalize.method)
   # Don't quantile normalize when using TMM
   # https://support.bioconductor.org/p/77664/

}



#' @importFrom magrittr %>%
explicitly_compute_precision_weights_once <- function(
   object,
   design = autonomics.import::create_design_matrix(object),
   plot = TRUE,
   ...
){

   # Extract
   log2cpm  <- autonomics.import::exprs(object)
   lib.size <- autonomics.import::sdata(object)$libsize

   # Assert
   n <- nrow(log2cpm)
   if (n < 2L) stop("Need at least two genes to fit a mean-variance trend")

   # Fit linear model
   fit <- limma::lmFit(log2cpm, design=design, ...)

   # Predict
   if (is.null(fit$Amean)) fit$Amean <- rowMeans(log2cpm, na.rm = TRUE)
   if (fit$rank < ncol(design)) {
      j <- fit$pivot[1:fit$rank]
      fitted.log2cpm <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,j, drop = FALSE])
   } else {
      fitted.log2cpm <- fit$coef %*% t(fit$design)
   }
   fitted.log2count <- (2^fitted.log2cpm) %>% cpm_to_counts(lib.size) %>% log2()

   # Fit mean-variance trend
   mean.log2count <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)  # mean log2 count
   sdrt.log2count <- sqrt(fit$sigma)                                     # sqrtsd(resid)
   all.identical <- matrixStats::rowVars(log2cpm)==0
   if (any(all.identical)) {
      mean.log2count <- mean.log2count[!all.identical]
      sdrt.log2count <- sdrt.log2count[!all.identical]
   }
   l <- lowess(mean.log2count, sdrt.log2count, f = 0.5)
   f <- approxfun(l, rule = 2)

   # Compute precision weights
   w <- 1/f(fitted.log2count)^4     # f(.) = sqrt(sd(.)) --> f(.)^4 = var(.)
   dim(w) <- dim(fitted.log2count)

   # Plot
   if (plot) {
      plot(mean.log2count, sdrt.log2count, xlab = "mean log2count", ylab = "sdrt log2count",
           pch = 16, cex = 0.25)
      title("voom: Mean-variance trend")
      lines(l, col = "red")
   }

   # Return
   return(w)
}

