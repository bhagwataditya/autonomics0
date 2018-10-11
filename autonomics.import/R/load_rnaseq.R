#' Compute effective library sizes (using TMM)
#' @param counts counts matrix
#' @return vector (nsample)
#' @export
libsizes <- function(counts){
   colSums(counts) * edgeR::calcNormFactors(counts)
}

#' Convert counts into cpm (counts per million reads)
#' @param counts    count matrix
#' @param lib.size  scaled library sizes (vector)
#' @return cpm matrix
#' @export
counts_to_cpm <- function(counts, lib.size = libsizes(counts)){
   t(t(counts + 0.5)/(lib.size + 1) * 1e+06)
}

#' Invert cpm (counts per million reads) into counts
#'
#' @param cpm       cpm matrix
#' @param lib.size  scaled library sizes (vector)
#' @return count matrix
#' @export
cpm_to_counts <- function(cpm, lib.size){
   1e-06 * t(t(cpm) * (lib.size + 1))
}

#' @importFrom magrittr %>%
compute_precision_weights_once <- function(
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


#' @rdname add_precision_weights
#' @importFrom magrittr %>%
#' @export
compute_precision_weights <- function(object, design = autonomics.import::create_design_matrix(object), plot = TRUE){

   # Estimate precision weights
   has_block <- autonomics.import::has_complete_block_values(object)
   weights <- object %>% compute_precision_weights_once(design = design, plot = !has_block)

   # Update precision weights using block correlation
   if (has_block){
      correlation <- log2cpm() %>%
         limma::duplicateCorrelation(design = design, block = object$block, weights = weights) %>%
         magrittr::extract2('consensus')
      weights <- object %>% compute_precision_weights_once(design = design,
                                                           block = object$block,
                                                           correlation = correlation,
                                                           plot = TRUE)
   }

   # Return
   dimnames(weights) <- dimnames(object)
   weights
}

#' Compute (and add) precision weights
#'
#' Compute (and add) precision weights for SummarizedExperiment with log2cpm values.
#'
#' Refactored version of limma::voom() which operates on a SummarizedExperiment with log2cpm values.
#' Created to allow a separation of the two steps in limma::voom():
#'    1. raw counts -> log2cpm (autonomics.import::counts_to_cpm)
#'    2. log2cpm    -> weights (autonomics.import::compute_precision_weights)
#' If a block factor is present, precision weights are updated using the block correlation information,
#' as suggested by Gordon Brown on BioC support (https://support.bioconductor.org/p/59700/)
#'
#' @param object    SummarizedExperiment
#' @param design    design matrix
#' @param plot      whether to plot mean-var trend (logical)
#' @param ...       passed to limma::lmFit
#' @return weights vector
#' @author Charity Law, Gordon Smyth, Aditya Bhagwat
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- 'extdata/rnaseq/gene_counts.txt'           %>%
#'               system.file(package = 'subramanian.2016') %>%
#'               autonomics.import::load_rnaseq(infer_design_from_sampleids = TRUE)
#'    object %>% autonomics.import::compute_precision_weights() %>% str()
#'    object %>% autonomics.import::add_precision_weights()
#'    object %>% autonomics.import::load_rnaseq()
#' }
#'
#' @importFrom magrittr  %>%
#' @export
add_precision_weights <- function(
   object,
   design = autonomics.import::create_design_matrix(object),
   plot = TRUE
){
   weights <- object %>% compute_precision_weights(design = design, plot = plot)
   SummarizedExperiment::assays(object)$weights <- weights
   object
}



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
#'    infer_design_from_sampleids <- TRUE
#'    object <- load_rnaseq(file, infer_design_from_sampleids = TRUE)
#'    object
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

   # Compute effective library size
   message('\t\tstore libsizes in sdata')
   autonomics.import::sdata(object)$libsize <- autonomics.import::exprs(object) %>% libsizes()

   # TMM-normalize counts
   message('\t\tcounts -> cpm')
   SummarizedExperiment::assays(object)$exprs %<>% counts_to_cpm()

   # Log2 transform
   message('\t\tcpm -> log2cpm')
   SummarizedExperiment::assays(object)$exprs %<>% log2()

   # Quantile normalize
   # Gordon:  prefer TMM over quantile normalization for most cases
   #          Use quantile normalization for extreme cases
   #          Don't use both
   # object %<>% limma::normalizeBetweenArrays(method = normalize.method)
   # Don't quantile normalize when using TMM
   # https://support.bioconductor.org/p/77664/

   # Add voom precision weights
   autonomics.support::cmessage('\t\tadd voom precision weights')
   object %<>% autonomics.import::add_precision_weights(plot = TRUE)

   # Add metadata
   autonomics.import::prepro(object) <- list(assay      = 'rnaseq',
                                             entity     = 'rna',
                                             quantity   = 'log2cpm',
                                             software   = 'limma',
                                             parameters = list())

   # Return
   object

}
