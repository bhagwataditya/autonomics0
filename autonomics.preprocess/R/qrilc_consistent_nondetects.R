
#' Guess representation of nondetects
#' @param x matrix or SummarizedExperiment
#' @return NA or 0
#' @examples
#' require(magrittr)
#' if (require(graumann.lfq)){
#'    x <- graumann.lfq::lfq.intensities
#'    x %>% guess_nondetect_representation()
#'    x %>% autonomics.import::exprs() %>% guess_nondetect_representation()
#' }
#' @export
guess_nondetect_representation <- function(x, ...){
   UseMethod("guess_nondetect_representation", x)
}

#' @rdname guess_nondetect_representation
#' @export
guess_nondetect_representation.default <- function(x, verbose = FALSE){
   representation <- if        (any(is.na(x))){   NA
                     } else if (any(x==0)    ){    0
                     } else if (any(x==-Inf) ){ -Inf
                     } else                   { NULL
                     }
   if (!is.null(representation) & verbose){
      autonomics.support::cmessage('\t\tGuesss nondetect representation: %s', as.character(representation))
   }
   representation
}


#' @rdname guess_nondetect_representation
#' @export
guess_nondetect_representation.SummarizedExperiment <- function(x, verbose = FALSE){
   x %>%
   autonomics.import::exprs() %>%
   guess_nondetect_representation(verbose = verbose)
}


#=========================================================================================================
#' Is nondetect
#' @param x matrix or SummarizedExperiment
#' @param representation nondetect representation
#' @return logical matrix with same dims as x
#' @examples
#' require(magrittr)
#' if (require(graumann.lfq)){
#'    x <- graumann.lfq::lfq.intensities
#'    x %>% is_nondetect()
#'    x %>% autonomics.import::exprs() %>% is_nondetect()
#'    x %>% autonomics.import::exprs() %>% magrittr::extract(, 1) %>% is_nondetect()
#' }
#' @export
is_nondetect <- function(x, ...){
   UseMethod("is_nondetect", x)
}


#' @rdname is_nondetect
#' @export
is_nondetect.default <- function(
   x,
   representation = guess_nondetect_representation(x)
){
   if        (is.null(representation)){ x %>% (function(y){ y[] <- 0; mode(y) <- 'logical'})
   } else if (is.na(representation)  ){ is.na(x)
   } else if (representation == 0    ){ x == 0
   } else if (representation == -Inf ){ x == -Inf
   }
}


#' @rdname is_nondetect
#' @export
is_nondetect.SummarizedExperiment <- function(
   x,
   representation = guess_nondetect_representation(x)
){
   x %>% autonomics.import::exprs() %>% is_nondetect()
}


#==================================================================================================================
#' Is complete nondetect
#'
#' Is feature nondetect in all samples
#'
#' @param x matrix or SummarizedExperiment
#' @param representation nondetect representation
#' @return logical matrix with same dims as x
#' @examples
#' require(magrittr)
#' if (require(graumann.lfq)){
#'    x <- graumann.lfq::lfq.intensities
#'    x %>% is_common_nondetect()
#'    x %>% autonomics.import::exprs() %>% is_common_nondetect()
#' }
#' @export
is_common_nondetect <- function(
   x,
   representation = guess_nondetect_representation(x)
){
   x %>% autonomics.preprocess::is_nondetect() %>% matrixStats::rowAlls()
}


#==================================================================================================================


#' QRILC common nondetects
#'
#' QRILC-impute values missing in all samples.
#'
#' @param object          SummarizedExperiment
#' @param representation  representation of missing values. Currently supported are NA, 0, and -Inf
#' @param verbose         logical(1)
#' @examples
#' if require(graumann.lfq){
#'    object <- graumann.lfq::lfq.intensities %>%
#'              autonomics.import::split_by_svar() %>%
#'              magrittr::extract2(2)
#'    object %>% qrilc_common_nondetects()
#' }
#' @return SummarizedExperiment with updated exprs
#' @importFrom magrittr %<>%
#' @export
qrilc_common_nondetects <- function(
   object,
   representation = guess_nondetect_representation(object),
   verbose = FALSE
){

   # Typify exprs
   consistent_nas <- object %>% is_common_nondetect(representation = representation)
   if (verbose)  autonomics.support::cmessage("\t\tImpute (QRILC) %4d features with value '%s' in all %d samples",
                                               sum(consistent_nas),
                                               as.character(representation),
                                               ncol(object))
   inconsistent_nas <- is_nondetect(object, representation) & !matrix(consistent_nas, nrow(object), ncol(object))

   # Temporarily replace inconsistent misses by median
   exprs1 <- autonomics.import::exprs(object)
   for (sample in 1:ncol(object)){
      exprs1[, sample] %<>% (function(x){x[is_nondetect(x, representation) & !consistent_nas] <- median(x, na.rm=TRUE); x})
   }

   # QRILC remaining consistent misses
   exprs1 %<>% imputeLCMD::impute.QRILC() %>% magrittr::extract2(1)

   # Re-NA inconsistent misses
   exprs1[inconsistent_nas] <- representation

   # Update object
   autonomics.import::exprs(object) <- exprs1

   # Return
   object
}


#' QRILC consistent nondetects
#'
#' QRILC nondetects shared across all samples of the same subgroup.
#'
#' @param object SummarizedExperiment
#' @param svar character(1)
#' @examples
#' require(magrittr)
#' if (require(graumann.lfq)){
#'    object <- graumann.lfq::lfq.intensities
#'    object %>% autonomics.preprocess::qrilc_consistent_nondetects()
#' }
#' @importFrom magrittr %>%
#' @export
qrilc_consistent_nondetects = function(
   object,
   svar = 'subgroup'
){

   if (!svar %in% autonomics.import::svars(object)){
      autonomics.support::cmessage("\t\tNo svar '%s' => no QRILC imputation performed", svar)
      return(object)
   }
   if (any(autonomics.support::is_missing_or_empty_character(autonomics.import::svalues(object, svar)))){
      autonomics.support::cmessage("\t\tSome '%s' values are missing => no QRILC imputation performed", svar)
      return(object)
   }

   object_list <- object %>% autonomics.import::split_by_svar(svar)
   n_consistent_nondetects <- object_list %>% lapply(autonomics.preprocess::is_common_nondetect) %>% vapply(sum, integer(1))
   autonomics.support::cmessage('\t\tQRILC-impute consistent nondetects among %d features', nrow(object))
   autonomics.support::cmessage_df('\t\t%s', n_consistent_nondetects)

   Map(qrilc_common_nondetects, object_list) %>%
   do.call(SummarizedExperiment::cbind, .)
}
