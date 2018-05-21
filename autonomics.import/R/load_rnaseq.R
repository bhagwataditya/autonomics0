REQUIRED_FVARS_RNASEQ <- c('gene_id', 'gene_name')
#' Required feature variables
#'
#' Required variables in annotated count file
"REQUIRED_FVARS_RNASEQ"

REQUIRED_SVARS_RNASEQ <- c('sample_id', 'subgroup', 'replicate')
#' Required sample variables
#'
#' Required variables in sample design file
"REQUIRED_SVARS_RNASEQ"

#' Load sample design file
#' @param sample_design_file sample design file
#' @return sample design dataframe
#' @examples
#' if (require(billing.differentiation.data)){
#'    sample_design_file <- system.file('extdata/rnaseq/sample_design.txt',
#'                             package = 'billing.differentiation.data')
#'    autonomics.import::load_sdata_rnaseq(sample_design_file)
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
load_sdata_rnaseq <- function(sample_design_file){

   # Check
   . <- NULL

   # Read
   sample_df <- sample_design_file %>% data.table::fread(data.table = FALSE)

   # sample_id has to be present
   sample_df$sample_id %<>% as.character()
   assertive.base::assert_all_are_not_na(sample_df$sample_id)                 # no NAs
   assertive.base::assert_all_are_true(stats::complete.cases(sample_df$sample_id)) # sample ids unique

   # subgroup and replicate are optional
   if ('subgroup'  %in% names(sample_df)){
      sample_df$subgroup  %<>% as.character()
      assertive.base::assert_all_are_not_na(sample_df$subgroup)
   }
   if ('replicate' %in% names(sample_df)){
      sample_df$replicate %<>% as.character() # should not be integer to allow mapping to shape
      assertive.base::assert_all_are_not_na(sample_df$replicate)
   }

   # Add rownames
   sample_df %<>% magrittr::set_rownames(.$sample_id)

   # Return
   sample_df

}



#' Load fdata
#' @param count_file  count file
#' @param map_entrezg whether to map features to their entrezg (logical)
#' @param verbose     whether to report messages (logical)
#' @examples
#' if (require(billing.differentiation.data)){
#'    count_file <- system.file('extdata/rnaseq/gene_counts.txt',
#'                               package = 'billing.differentiation.data')
#'    load_fdata_rnaseq(count_file)
#' }
#' if (require(subramanian.2016)){
#'    count_file <- system.file('extdata/rnaseq/gene_counts.txt',
#'                               package = 'subramanian.2016')
#' }
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%  %<>%
#' @export
load_fdata_rnaseq <- function(count_file, map_entrezg = FALSE, verbose = TRUE){

   # Satisfy CHECK
   . <- entrezg <- gene_id <- gene_name <- NULL

   fdata1 <- data.table::fread(count_file)  %>% magrittr::extract(, !sapply(., is.integer), with = FALSE)
   assertive.sets::assert_is_subset(REQUIRED_FVARS_RNASEQ, names(fdata1))
   fdata1 %<>% magrittr::extract(, gene_id := gene_id %>%
                                              stringi::stri_split_fixed('.') %>%
                                              vapply(magrittr::extract, character(1), 1)) #%>%
               #data.table::setnames('gene_id', 'feature_id') %>%
               #magrittr::extract(, ensg := feature_id)
   organism <- fdata1[, autonomics.annotate::infer_organism(gene_id, keytype = 'ensg', verbose = verbose)]

   # Add entrezg mappings
   if (map_entrezg){
      suppressMessages(fdata1[, entrezg := autonomics.annotate::ensg_to_entrezg(gene_id,      organism = organism, verbose = FALSE)])
      idx1 <- !is.na(fdata1[, entrezg])

      suppressMessages(fdata1[is.na(entrezg), entrezg := autonomics.annotate::gsymbol_to_entrezg(gene_name, organism = organism)])
      idx2 <- !is.na(fdata1[, entrezg])
      if (verbose){
      autonomics.support::cmessage('\t\t%d/%d features mapped to entrezg: %d through ensg, %d through gsymbol',
                                   sum(idx2), length(idx2), sum(idx1), sum(idx2)-sum(idx1))
      }
   }

   # Return
   fdata1 %<>% as.data.frame() %>% magrittr::set_rownames(.$gene_id)
   fdata1
}

#' Load exprs
#' @param count_file count file
#' @return count matrix
#' @importFrom magrittr %>%
#' @export
load_exprs_rnaseq <-    function(count_file){
   . <- NULL
   feature_ids <- load_fdata_rnaseq(count_file, map_entrezg = FALSE, verbose = FALSE) %>% magrittr::extract2('gene_id')
   exprs1 <- data.table::fread(count_file) %>%
             magrittr::extract(, sapply(., is.integer), with = FALSE) %>%
             data.matrix() %>%
             magrittr::set_rownames(feature_ids)
}

#' Load RNA seq counts
#' @param dir                directory with count file and sample design file
#' @param count_file         gene count file
#' @param sample_file        sample design file
#' @param map_entrezg        logical
#' @return eset
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    dir <- system.file('extdata/rnaseq', package = 'billing.differentiation.data')
#'    autonomics.import::load_rnaseq_counts(dir)
#' }
#' if (require(subramanian.2016)){
#'    dir <- system.file('extdata/rnaseq', package = 'subramanian.2016')
#'    autonomics.import::load_rnaseq_counts(dir, sample_file = NULL)
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
load_rnaseq_counts <- function(
   dir         = '.',
   count_file  = paste0(dir, '/gene_counts.txt'),
   sample_file = if (file.exists(paste0(dir, '/sample_design.txt'))) paste0(dir, '/sample_design.txt') else NULL,
   map_entrezg = FALSE
){

   # Assert
   assertive.files::assert_all_are_non_empty_files(count_file)
   if (!is.null(sample_file)) assertive.files::assert_all_are_non_empty_files(sample_file)

   # Load components
   fdata1 <- autonomics.import::load_fdata_rnaseq(count_file, map_entrezg = map_entrezg, verbose = TRUE)
   exprs1 <- autonomics.import::load_exprs_rnaseq(count_file)
   sdata1 <- if (is.null(sample_file)){
                data.frame(sample_id = colnames(exprs1), row.names = colnames(exprs1)) %>%
                autonomics.import::infer_add_design()
             } else {
                autonomics.import::load_sdata_rnaseq(sample_file)
             }
   prepro1 <- autonomics.import::create_prepro_list(
                 assay      = 'rnaseq',
                 entity     = 'rna',
                 quantity   = 'voomcount',
                 software   = 'limma',
                 parameters = list())

   # Forge eset
   eset1 <- SummarizedExperiment::SummarizedExperiment(list(exprs = exprs1))
   autonomics.import::fdata(eset1) <- fdata1
   autonomics.import::sdata(eset1) <- sdata1
   autonomics.import::prepro(eset1) <- prepro1
   autonomics.import::assert_is_valid_eset(eset1)

   # Filter
   eset1 %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample()

   # Return
   eset1
}


#' Convert raw counts into log cpm
#'
#' Convert raw counts into log cpm in exactly the same way as limma::voom does
#'
#' @param object eset
#' @return log cpm
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    eset1 <- billing.differentiation.data::rna.counts
#'    eset1 %>% autonomics.import::exprs()    %>% magrittr::extract(1:3, 1:3)
#'    eset1 %>% autonomics.import::logcpm() %>%
#'              autonomics.import::exprs() %>% magrittr::extract(1:3, 1:3)
#' }
#' @export
logcpm <- function(object){
   counts <- autonomics.import::exprs(object)
   lib.size <- colSums(counts)
   autonomics.import::exprs(object) <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
   object
}


#' Voom transform eset
#' @param object eset
#' @param normalize.method normalization method
#' @param plot whether to plot mean variance trend
#' @param verbose whether to report progress
#' @export
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'
#'    # eset with counts
#'       object <- billing.differentiation.data::rna.counts
#'       object %>% autonomics.import::exprs() %>% extract(1:3, 1:3)
#'
#'    # voom transformation
#'          object %>% autonomics.import::voom_transform() %>% autonomics.import::exprs() %>%
#'                     magrittr::extract(1:3, 1:3)
#'       # containing logcpm exprs
#'          object %>% autonomics.import::logcpm() %>% autonomics.import::exprs() %>%
#'                     extract(1:3, 1:3)
#'          object %>% autonomics.import::exprs() %>% magrittr::extract(1:3, 1:3)
#'       # and precision weights
#'          object %>% voom_transform() %>% autonomics.import::weights() %>%
#'                     extract(1:3, 1:3)
#' }
#' if (require(subramanian.2016)){
#'    dir <- system.file('extdata/rnaseq', package = 'subramanian.2016')
#'    object <- autonomics.import::load_rnaseq_counts(dir, sample_file = NULL)
#'    object %>% autonomics.import::voom_transform()
#' }
#' @export
voom_transform <- function(object, normalize.method = 'none', plot = TRUE, verbose = TRUE){

   autonomics.import::assert_is_valid_eset(object)
   assertive.sets::assert_is_subset('subgroup', autonomics.import::svars(object))

   # Normalize library size
   if (verbose)   autonomics.support::cmessage('\t\tNormalize library size')
   dge <- edgeR::DGEList(
      counts  = autonomics.import::exprs(object),
      samples = autonomics.import::sdata(object),
      group   = autonomics.import::sdata(object)$subgroup,
      genes   = autonomics.import::fdata(object)
   )
   dge <- edgeR::calcNormFactors(dge)

   # Voom transform
   if (verbose)   autonomics.support::cmessage('\t\tVoom transform')
   my_design <- stats::model.matrix(~ 0 + subgroup, data = autonomics.import::sdata(object))
   v <- limma::voom(dge, my_design, plot = plot, normalize.method = normalize.method, save.plot = TRUE)

   # Correct voom transformation for block effect
   # https://support.bioconductor.org/p/59700/
   if (autonomics.import::contains_block(object)){
      if (verbose)   autonomics.support::cmessage('\t\tAccount for block effect and rerun voom transformation')
      corfit <- limma::duplicateCorrelation(v, design = my_design, block = object$block)
      v <- limma::voom(dge, my_design, block = object$block, correlation = corfit$consensus,
                       plot = plot, normalize.method = normalize.method, save.plot = TRUE)
   }

   # Add voom transfo to SummarizedExperiment
   SummarizedExperiment::assays(object)$counts  <- SummarizedExperiment::assays(object)$exprs
   SummarizedExperiment::assays(object)$exprs   <- v$E
   SummarizedExperiment::assays(object)$weights <- v$weights
   S4Vectors::metadata(object)$voom.xy   <- v$voom.xy
   S4Vectors::metadata(object)$voom.line <- v$voom.line

   # Add prepro info
   autonomics.import::prepro(object) <- autonomics.import::prepro(object)

   # Return
   object
}


#' Plot voom mean variance trend
#' @param object SummarizedExperiment
#' @param normalize.method passed to limma::voom
#' @param title plot title
#' @return plot object
#' @importFrom magrittr %>%
#' @export
voom_plot <- function(object, normalize.method, title){
   x <- S4Vectors::metadata(object)$voom.xy$x
   y <- S4Vectors::metadata(object)$voom.xy$y
   xlab   <- S4Vectors::metadata(object)$voom.xy$xlab
   ylab   <- S4Vectors::metadata(object)$voom.xy$ylab
   graphics::plot(x, y, pch = 16, cex = 0.25, xlab = xlab, ylab = ylab)
   l <- stats::lowess(x, y, f=0.5)
   graphics::lines(l, col = 'red')
   title(title)
}

invlog2 <- function(x, pseudocount=0.5){
   function(x) 2^x - 0.5
}

#' Is mean var trend monotonic
#' @param object SummarizedExperiment
#' @param normalize.method normalization method
#' @return logical
#' @importFrom magrittr %>%
#' @export
mean_var_monotonic <- function(object, normalize.method){
   object %>%
   autonomics.import::voom_transform(normalize.method = normalize.method, verbose = FALSE) %>%
   (function(object){
      x <- S4Vectors::metadata(object)$voom.line$x
      y <- S4Vectors::metadata(object)$voom.line$y
      x[which.max(y)] <= x[1]
   })
}

#' Voom filter features
#'
#' Filter low count features (for all samples) to achieve monotonic mean variance trend,
#' as is recommended by Gordon Smyth & co. (see references)
#'
#' @param object exprs object
#' @param normalize.method passed to limma::voom
#' @param file file to print mean var curve to
#' @param width  fig width in inches
#' @param height fig height in inches
#' @return filtered eset
#' @references
#' https://support.bioconductor.org/p/64484/
#' https://support.bioconductor.org/p/80331/
#' @importFrom magrittr %>% %<>%
#' @export
voom_filter_n_transform <- function(object, normalize.method = 'none', file = NULL, width = 10, height = 5){

   # Filter until inflection point vanishes
   object0 <- object
   filter_threshold <- 0
   while(!mean_var_monotonic(object, normalize.method)){
      filter_threshold %<>% magrittr::add(1)
      object %<>% autonomics.import::filter_features_min_expr(filter_threshold + 1)
   }

   # Transform
   object %<>% autonomics.import::voom_transform(normalize.method = normalize.method)

   # Log
   summary_attr <- autonomics.import::get_summary_attr(object) %>%
                   c(list(voom_filter_n         = unname(nrow(object)),
                          voom_filter_threshold = filter_threshold))
   attr(object, 'autonomics_summary') <- summary_attr

   # Plot
   if (!is.null(file))   grDevices::pdf(file, width = width, height = height)
   graphics::layout(matrix(1:2, ncol = 2))
   object0 %>% voom_transform(plot = FALSE) %>%
               autonomics.import::voom_plot(normalize.method = normalize.method,
                                            title = sprintf('all (%d)', nrow(object0   )))
   object %>% autonomics.import::voom_plot(normalize.method = normalize.method,
                                           title = sprintf('expr < %0.0f (%d)', filter_threshold, nrow(object)))
   if (!is.null(file))   grDevices::dev.off()


   # Return
   return(object)
}

