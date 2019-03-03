
#' Impute
#'
#' Data from proteomics and metabolomics LCMS experiments are missing not at random (MNAR).
#' They are well imputed with the QRILC (quantile regression imputation of left censored data).
#'
#' See \href{https://cran.r-project.org/package=imputeLCMD}{imputeLCMD package on CRAN}
#'
#' @param object SummarizedExperiment, eSet, or EList
#' @param method Imputation method employed. See Section 'Details'.
#' @param random_seed Integer used with \code{\link[R.utils]{withSeed}} for reproducible
#' imputation (the used \code{\link[imputeLCMD]{impute.QRILC}} implys random draws from a distribution)
#' @param ... Parameters handed on to the imputation methods.
#' @details The function under the performs no imutation if \code{method = 'none'}, uses
#'  \code{\link[missForest]{missForest}} if \code{method = 'missForest'} and the parameter-named
#'  methods from \pkg{imputeLCMD} in all other cases.
#' @return updated object
#' @importFrom imputeLCMD impute.MinDet impute.MinProb impute.QRILC
#' @importFrom magrittr %<>% %>%
#' @importFrom missForest missForest
#' @export
impute <- function(
  object,
  method = c('none', 'missForest', 'impute.minDet', 'impute.minProb', 'impute.QRILC'),
  random_seed = NULL,
  ...)
{
# Check prerequisites -----------------------------------------------------
  autonomics.import::assert_is_valid_object(object)
  method %<>%
    match.arg(
      choices = c('none', 'missForest', 'impute.MinDet', 'impute.MinProb', 'impute.QRILC'),
      several.ok = FALSE)
  if(!is.null(random_seed))
  {
    assertive.types::assert_is_a_number(random_seed)
  }

# Processing --------------------------------------------------------------
  # Eject if nothing is to be done
  if(method == 'none')
  {
    return(object)
  }

  # Retrieve the function to be used
  # https://stackoverflow.com/a/10022480/2103880
  pkg <- method %>%
    switch(
      missForest = 'missForest',
      'imputeLCMD')
  fn <- get(method, asNamespace(pkg))

  largs <- list(autonomics.import::exprs(object)) %>%
    magrittr::set_names(
      method %>%
        switch(
          missForest = 'xmis',
          'dataSet.mvs'))

  # Execute dependent on 'random_seed'
  ## Capture output is currently needed, as imputeLCMD::impute.MinProb contains a spurious print statement (v.2.0)
  if(is.null(random_seed))
  {
    utils::capture.output(
      impute_result <- do.call(fn, largs, ...))
  } else {
    utils::capture.output(
      impute_result <- R.utils::withSeed(
        { do.call(fn, autonomics.import::exprs(object), ...) },
        seed = random_seed))
  }

  # Extract dependent on method
  # missForest --> list[['ximp']]
  # MinDet --> matrix
  # MinProb --> matrix
  # QRILC --> list[[1]]
  index <- method %>%
    switch(
      missForest   = 'ximp',
      impute.QRILC = 1,
      NULL)
  if(!is.null(index))
  {
    impute_result %<>%
      magrittr::extract2(index)
  }

  # Reassemble & Return
  autonomics.import::exprs(object) <- impute_result
  return(object)
}


