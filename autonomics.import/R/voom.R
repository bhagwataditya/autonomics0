#' Convert into cpm and compute precision weights for linear modeling
#'
#' @param object SummarizedExperiment
#' @param design
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- 'extdata/rnaseq/gene_counts.txt'            %>%
#'               system.file(package = 'subramanian.2016')  %>%
#'               autonomics.import::load_rnaseq(infer_design_from_sampleids = TRUE)
#' }
voom <- function (
   object,
   design,
   lib.size = NULL,
   normalize.method = "none",
   plot = FALSE,
   save.plot = FALSE, ...){

   # Fit model
   log2cpm <- counts %>% counts_to_cpm() %>% log2() %>% limma::normalizeBetweenArrays(method = normalize.method)

}

compute_precision_weights_core <- function(
   object,
   design = autonomics.find::create_design_matrix(object),
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

compute_precision_weights <- function(object, design = autonomics.find::create_design_matrix(object), plot = TRUE){

   # Estimate precision weights
   has_block <- autonomics.import::contains_block(object)
   weights <- object %>% compute_precision_weights(lib.size, design = design, plot = !has_block)

   # Update precision weights using block correlation
   if (has_block){
      correlation <- log2cpm() %>%
                     limma::duplicateCorrelation(design = design, block = object$block, weights = weights) %>%
                     magrittr::extract2('consensus')
      weights <- object %>% compute_precision_weights(design = design, block = object$block, correlation = correlation, plot = TRUE)
   }

   # Return
   weights

}

#' Compute precision weights (and add to sumexp)
#'
#' Compute precision weights for SummarizedExperiment with log2cpm exprs.
#'
#' Refactored and modularized version of limma::voom().
#' Created to separate limma::voom() into two separate steps:
#'    1. raw counts -> log2cpm (autonomics.import::counts_to_cpm)
#'    2. log2cpm    -> weights (autonomics.import::compute_precision_weights)
#'
#' @param object    SummarizedExperiment
#' @param design    design matrix
#' @param plot      whether to plot mean-var trend (logical)
#' @param ...       passed to limma::lmFit
#' @return weights vector
#' @author Charity Law, Gordon Smyth, Aditya Bhagwat
#' @importFrom magrittr  %>%
#' @export
add_precision_weights <- function(
   object,
   design = autonomics.find::create_design_matrix(object),
   plot = TRUE
){
   weights <- object %>% compute_precision_weights(design = design, plot = plot)
   SummarizedExperiment::assays(object)$weights <- weights
   object
}


