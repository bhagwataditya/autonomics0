#=================
# FILTER EXPRS
#=================

utils::globalVariables('.')


#' Filter features with replicated expression in some subgroup
#' @param object      SummarizedExperiment
#' @param comparator  '>' or '!='
#' @param lod         numeric(1)
#' @return Filtered SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% filter_exprs_replicated_in_some_subgroup()
#'    object <- autonomics.data::glutaminase
#'    object %>% filter_exprs_replicated_in_some_subgroup()
#' }
#' @importFrom magrittr %>%
#' @export
filter_exprs_replicated_in_some_subgroup <- function(
    object,
    comparator = if (contains_ratios(object)) '!=' else '>',
    lod = 0
){

    # Return if no subgroup (filtering not possible) or no replicates (filtering leads to empty SumExp )
    if (!'subgroup' %in% svars(object))             return(object)
    if (all(!duplicated(sdata(object)$subgroup)))   return(object)

    # Datatablify
    replicated_in_its_subgroup <- replicated_in_any_subgroup <- value <- NULL
    dt <- object %>% sumexp_to_long_dt(svars = 'subgroup')

    # Is expr replicated in its subgroup?
    if (comparator ==  '>') dt %>% magrittr::extract(, replicated_in_its_subgroup := sum(value  > lod, na.rm=TRUE) > 1,   by = c('feature_id', 'subgroup'))
    if (comparator == '!=') dt %>% magrittr::extract(, replicated_in_its_subgroup := sum(value != lod, na.rm=TRUE) > 1,   by = c('feature_id', 'subgroup'))

    # Is it replicated in any subgroup
    dt %>% magrittr::extract(, replicated_in_any_subgroup := any(replicated_in_its_subgroup),    by =  'feature_id')

    # Keep only replicated features
    replicated_features <- dt %>% magrittr::extract(replicated_in_any_subgroup == TRUE, 'feature_id')
    idx <- fid_values(object) %in% replicated_features
    autonomics.support::cmessage('\t\tFilter %d/%d features: expr %s %s, for at least two samples in some subgroup', sum(idx), length(idx), comparator, as.character(lod))
    object %<>% extract_features(idx)
    # use this rather than magrittr::extract() directly
    # to ensure that the limma(.) object is properly taken care of
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>%
            c(structure(
                sum(idx),
                names = sprintf(
                    "expr %s %s, for at least two samples in some subgroup",
                    comparator, as.character(lod))))
    }
    object
}




#=================
# FILTER FEATURES
#=================

#' @rdname filter_features
#' @export
filter_features_ <- function(object, condition, verbose = FALSE){
    if (is.null(condition)) return(object)
    idx <- lazyeval::lazy_eval(condition, fdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx), ' features: ', if (class(condition)=='lazy') deparse(condition$expr) else condition)
    object %<>% extract_features(idx)
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>%
            c(structure(sum(idx), names = condition))
    }
    object
}

#' Filter features on condition
#' @param object SummarizedExperiment
#' @param condition filter condition
#' @param verbose logical
#' @return filtered eSet
#' @examples
#' require(magrittr)
#' \dontrun{
#' if (require(autonomics.data)){
#'    autonomics.data::ALL %>%
#'    filter_features(gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1') )
#'
#'    autonomics.data::ALL %>%
#'    filter_features(gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1'), verbose = TRUE)
#'
#'    autonomics.data::ALL %>%
#'    filter_features_("gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1')")
#'
#'    autonomics.data::ALL %>%
#'    filter_features_("gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1')", verbose = TRUE)
#'
#'    file <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'    system.file(package = 'autonomics.data')
#'    sumexp <- file %>% read_omics(fid_rows   = 2:9783,  fid_cols   = 383,
#'                        sid_rows   = 1,       sid_cols   = seq(124, 316, by = 6),
#'                        expr_rows  = 2:9783,  expr_cols  = seq(124, 316, by = 6),
#'                        fvar_rows  = 1,       fvar_cols  = c(2, 6, 7, 383),
#'                        fdata_rows = 2:9783,  fdata_cols = c(2, 6, 7, 383),
#'                        transpose  = FALSE)
#'    analysis(sumexp)
#'    sumexp %<>% filter_features_("c('UBA6', 'KRT8') %in% `Gene names`")
#'    analysis(sumexp)
#' }
#' }
#' @export
filter_features <- function(object, condition, verbose = FALSE){
    # earlier version based on lazyeval (still functional, but soft deprecated by Hadley and co)
    # filter_features_(object, lazyeval::lazy(condition), verbose = verbose)
    condition <- rlang::enquo(condition)
    idx <- rlang::eval_tidy(condition, fdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx), ' features: ', rlang::expr_text(condition))
    object %<>% magrittr::extract(idx,)
    fdata(object) %<>% droplevels()
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>%
            c(structure(sum(idx), names = rlang::as_name(condition)))
    }
    object
}

#' Identify reference features
#' @param object   eset
#' @param fvar       fvar
#' @param split      split
#' @param fvalues    fvalues
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    fvalues <- c("A0MZ66", "A0JLT2", "Q8NG66", "Q9NPI5")
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% is_fvalue_feature(fvar = 'Uniprot accessions',
#'                                 split = ';', fvalues = fvalues) %>%
#'               sum()
#' }
#' @importFrom magrittr  %>%
#' @export
is_fvalue_feature <- function(object, fvar, split, fvalues){
    fdata(object) %>%
        magrittr::extract2(fvar) %>%
        as.character() %>%
        strsplit(split) %>%
        vapply(function(x) any(x %in% fvalues), logical(1))
}


#' Filter subet for reference features
#' @param  object            eSet
#' @param  fvar                fvar
#' @param  split               split
#' @param  fvalues             fvalues
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    fvalues <- c("A0MZ66", "A0JLT2", "Q8NG66", "Q9NPI5")
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% filter_features_on_fvalues(fvar    = 'Uniprot accessions',
#'                                          split   = ';',
#'                                          fvalues = fvalues)
#'    if(is.null(analysis(object))) {
#'        analysis(object) <- list(nsamples = nrow(sdata(object)))
#'    }
#'    analysis(object)$nsamples
#'    object %<>% filter_features_on_fvalues(fvar    = 'Uniprot accessions',
#'                                          split   = ';',
#'                                          fvalues = fvalues)
#'    analysis(object)$nsamples
#' }
#' @return subset of object with reference features
#' @importFrom magrittr       %>%     %<>%
#' @export
filter_features_on_fvalues <- function(object, fvar, split, fvalues){
    fvalues %<>% unique()
    idx <- object %>% is_fvalue_feature(fvar, split, fvalues)
    object %<>% magrittr::extract(idx, )
    if (!is.null(analysis(object))) {
        analysis(object)$nsamples %<>%
            c(structure(
                sum(idx),
                names = sprintf("c('%s') %%in%% `%s`",
                                paste(fvalues, collapse = "', '"), fvar)))
    }
    object
}




#=========================
# FILTER SAMPLES
#=========================

#' @rdname filter_samples
#' @importFrom magrittr %<>%
#' @export
filter_samples_ <- function(object, condition, verbose = FALSE, record = TRUE){
    if (is.null(condition)) return(object)
    idx <- lazyeval::lazy_eval(condition, sdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx), ' samples: ', if (class(condition)=='lazy') deparse(condition$expr) else condition)
    object %<>% magrittr::extract(, idx)
    sdata(object) %<>% droplevels()
    if (recrod && !is.null(analysis(object))) {
        analysis(object)$nsamples %<>%
            c(structure(
                sum(idx),
                names = if (class(condition)=='lazy') deparse(condition$expr) else condition))
    }
    object
}

#' Filter samples on condition
#' @param object SummarizedExperiment
#' @param condition filter condition
#' @param verbose logical
#' @param record TRUE (default) or FALSE
#' @return filtered SummarizedExperiment
#' @examples
#' if (require(autonomics.data)){
#'    filter_samples( autonomics.data::ALL,  sex == 'M')
#'    filter_samples_(autonomics.data::ALL, "sex == 'M'")
#'
#'    sumexp <- autonomics.data::ALL
#'    if (is.null(analysis(sumexp))) {
#'        analysis(sumexp) <- list(nsamples = nrow(sdata(sumexp)))
#'    }
#'    analysis(sumexp)$nsamples
#'    sumexp <- filter_samples_(sumexp, "sex == 'M'")
#'    analysis(sumexp)$nsamples
#' }
#' @export
filter_samples <- function(object, condition, verbose = FALSE, record = TRUE){
    # earlier version based on lazyeval (still functional, but soft deprecated by Hadley and co)
    # filter_samples_(object, lazyeval::lazy(condition), verbose = verbose)
    condition <- rlang::enquo(condition)
    idx <- rlang::eval_tidy(condition, sdata(object))
    idx <- idx & !is.na(idx)
    if (verbose) if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx), ' samples: ', rlang::expr_text(condition))
    object %<>% magrittr::extract(, idx)
    sdata(object) %<>% droplevels()
    if (record && !is.null(analysis(object))) {
        analysis(object)$nsamples %<>%
            c(structure(
                sum(idx),
                names = rlang::expr_text(condition)))
    }
    object
}

