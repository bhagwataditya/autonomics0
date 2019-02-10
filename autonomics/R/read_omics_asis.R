#=================================================
# GENERIC
#=================================================

#' Read omics data as is from file
#' @param file character(1)
#' @param sheet numeric(1) or character(1)
#' @return data.table
#' @importFrom magrittr %>%
#' @export
read_file <- function(file, sheet){
   assertive.files::assert_all_are_existing_files(file)
   if (tools::file_ext(file) %in% c('xls', 'xlsx')){
      file %>% readxl::read_excel(sheet = sheet, col_names = FALSE) %>% data.table::data.table()
   } else {
      file %>% data.table::fread(na.strings = "", header = FALSE, integer64 = 'numeric')
   }
}

# Extract row
extract_dt_row <- function(dt, i) dt %>% magrittr::extract(i,) %>% as.matrix() %>% magrittr::extract(1, ) %>% unname()
extract_dt_col <- function(dt, i) dt[[i]]


#' Extract rectangle from omics datatable
#' @param dt         datatable with raw omics data
#' @param rows       row selector: logical(n), numeric(n), or character(n)
#' @param cols       col selector: logical(n), numeric(n), or character(n)
#' @param transpose  logical(1)
#' @param drop       logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    dt <- file %>% read_file(sheet = 2)
#'    dt %>% extract_rectangle(rows = 11:401, cols = 15:86) %>% extract(1:3, 1:3)           # exprs
#'    dt %>% extract_rectangle(rows = 11:401, cols = 5    ) %>% extract(1:3, , drop=FALSE)  # fids
#'    dt %>% extract_rectangle(rows = 2,      cols = 15:86) %>% extract(, 1:3, drop=FALSE)  # sids
#'    dt %>% extract_rectangle(rows = 10:401, cols = 1:14 ) %>% extract(1:3, 1:3)           # fdata
#'    dt %>% extract_rectangle(rows = 1:10,   cols = 14:86, transpose = TRUE)  %>% extract(1:3, 1:3)
#'                                                                                          # sdata
#' }
#' @importFrom magrittr %>% %<>%
#' @export
extract_rectangle <- function(dt, rows, cols, transpose = FALSE, drop = FALSE){
   dt %<>% magrittr::extract(rows, cols, with = FALSE)
   rectangle <- dt %>% as.matrix()
   if (transpose) rectangle %<>% t()
   if (drop) if (nrow(rectangle)==1 | ncol(rectangle)==1) rectangle %<>% as.vector('character')
   rectangle
}

#' Read omics data
#' @param file       character(1): name of text (txt, csv, tsv) or excel (xls, xlsx) file
#' @param sheet      integer(1) or character(1)
#' @param fid_rows   featureid rows: logical(n) or numeric(n)
#' @param fid_cols   featureid cols: logical(n) or numeric(n)
#' @param sid_rows   sampleid  rows: logical(n) or numeric(n)
#' @param sid_cols   sampleid  cols: logical(n) or numeric(n)
#' @param expr_rows  expr      rows: logical(n) or numeric(n)
#' @param expr_cols  expr      rows: logical(n) or numeric(n)
#' @param fvar_rows  fvar      rows: logical(n) or numeric(n)
#' @param fvar_cols  fvar      cols: logical(n) or numeric(n)
#' @param svar_rows  svar      rows: logical(n) or numeric(n)
#' @param svar_cols  svar      cols: logical(n) or numeric(n)
#' @param fdata_rows fdata     rows: logical(n) or numeric(n)
#' @param fdata_cols fdata     cols: logical(n) or numeric(n)
#' @param sdata_rows sdata     rows: logical(n) or numeric(n)
#' @param sdata_cols sdata     cols: logical(n) or numeric(n)
#' @param transpose  logical(1)
#' @param verbose    logical(1)
#' @examples
#' require(magrittr)
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'   file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'   file %>% read_omics_asis(sheet      = 2,
#'                            fid_rows   = 11:401,    fid_cols   = 5,
#'                            sid_rows   = 3,         sid_cols   = 15:86,
#'                            expr_rows  = 11:401,    expr_cols  = 15:86,
#'                            fvar_rows  = 10,        fvar_cols  = 1:14,
#'                            svar_rows  = 1:10,      svar_cols  = 14,
#'                            fdata_rows = 11:401,    fdata_cols = 1:14,
#'                            sdata_rows = 1:10,      sdata_cols = 15:86,
#'                            transpose  = FALSE)
#' }
#'
#' # SOMASCAN
#' if (require(atkin.2014)){
#'    file <- 'extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat' %>%
#'             system.file(package = 'atkin.2014')
#'    file %>% read_omics_asis(fid_rows   = 73,       fid_cols   = 27:1343,
#'                             sid_rows   = 93:520,   sid_cols   = 6:6,
#'                             expr_rows  = 93:520,   expr_cols  = 27:1343,
#'                             fvar_rows  = 73:91,    fvar_cols  = 26,
#'                             svar_rows  = 92,       svar_cols  = 1:25,
#'                             fdata_rows = 73:91,    fdata_cols = 27:1343,
#'                             sdata_rows = 93:520,   sdata_cols = 1:25,
#'                             transpose  = TRUE)
#' }
#'
#' # EXIQON
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx', package = 'subramanian.2016')
#'    file %>% read_omics_asis(sheet      = 1,
#'                             fid_rows   = 1,        fid_cols   = 2:330,
#'                             sid_rows   = 2:73,     sid_cols   = 1,
#'                             expr_rows  = 2:73,     expr_cols  = 2:330,
#'                             fvar_rows  = 74:76,    fvar_cols  = 1,
#'                             svar_rows  = 1,        svar_cols  = 331:332,
#'                             fdata_rows = 74:76,    fdata_cols = 2:330,
#'                             sdata_rows = 2:73,     sdata_cols = 331:332,
#'                             transpose  = TRUE)
#' }
#'
#' # PROTEINGROUPS LABELFREE
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% read_omics_asis(fid_rows   = 2:5893,   fid_cols   = 275,
#'                             sid_rows   = 1,        sid_cols   = 227:248,
#'                             expr_rows  = 2:5893,   expr_cols  = 227:248,
#'                             fvar_rows  = 1,        fvar_cols  = 1:8,
#'                             fdata_rows = 2:5893,   fdata_cols = 1:8,
#'                             transpose  = FALSE)
#' }
#'
#' # PROTEINGROUPS LABELED
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemdiff/maxquant/proteinGroups.txt', package = 'autonomics.data')
#'    file %>% read_omics_asis(fid_rows   = 2:9783,  fid_cols   = 383,
#'                             sid_rows   = 1,       sid_cols   = seq(124, 316, by = 6),
#'                             expr_rows  = 2:9783,  expr_cols  = seq(124, 316, by = 6),
#'                             fvar_rows  = 1,       fvar_cols  = c(2, 6, 7, 383),
#'                             fdata_rows = 2:9783,  fdata_cols = c(2, 6, 7, 383),
#'                             transpose  = FALSE)
#' }
#'
#' # RNASEQ
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    file %>% read_omics_asis(fid_rows   = 2:48821,   fid_cols   = 1,
#'                             sid_rows   = 1,         sid_cols   = 5:86,
#'                             expr_rows  = 2:48821,   expr_cols  = 5:86,
#'                             fvar_rows  = 1,         fvar_cols  = 1:4,
#'                             fdata_rows = 2:48821,   fdata_cols = 1:4,
#'                             transpose  = FALSE)
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_omics_asis <- function(
   file,
   sheet = 1,
   fid_rows,           fid_cols,
   sid_rows,           sid_cols,
   expr_rows,          expr_cols,
   fvar_rows  = NULL,  fvar_cols  = NULL,
   svar_rows  = NULL,  svar_cols  = NULL,
   fdata_rows = NULL,  fdata_cols = NULL,
   sdata_rows = NULL,  sdata_cols = NULL,
   transpose  = FALSE,
   verbose    = FALSE
){
   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.types::assert_is_a_bool(transpose)
   
   # read
   dt <- read_file(file, sheet)
   
   # Extract exprs
   fids1  <- dt %>% extract_rectangle(fid_rows,  fid_cols,  transpose = transpose, drop = TRUE)
   sids1  <- dt %>% extract_rectangle(sid_rows,  sid_cols,  transpose = transpose, drop = TRUE)
   exprs1 <- dt %>% extract_rectangle(expr_rows, expr_cols, transpose = transpose) %>% (function(y){class(y) <- 'numeric'; y})
   
   # Extract feature annotations
   #    Leave rownames(fdata1) empty: fids1 may contain non-valid values
   #    This happens in MaxQuant files, which sometimes contain missing rows (I think after opening in excell)
   fdata1 <- data.frame(feature_id = fids1, stringsAsFactors = FALSE)
   fdata_available <- !is.null(fvar_rows) & !is.null(fvar_cols)
   if (fdata_available){
      fvars1 <- dt %>% extract_rectangle(fvar_rows,  fvar_cols,  transpose = transpose, drop = TRUE)
      fdata1 %<>% cbind(extract_rectangle(dt, fdata_rows, fdata_cols, transpose = transpose) %>%
                           magrittr::set_colnames(fvars1) %>%
                           data.frame(stringsAsFactors = FALSE, check.names = FALSE))
   }
   
   # Extract sample annotations
   #    Leave rownames(sdata1) empty: sids may contain non-valid values
   #    This happens in SOMA files, where CLIENT_IDENTIFIER is not unique for calibrator and buffer samples
   sdata1 <- data.frame(sample_id = sids1, stringsAsFactors = FALSE)
   sdata_available <- !is.null(svar_rows) & !is.null(svar_cols)
   if (sdata_available){
      svars1 <- dt %>% extract_rectangle(svar_rows,  svar_cols,  transpose =  transpose, drop = TRUE)
      sdata1 %<>% cbind(extract_rectangle(dt, sdata_rows, sdata_cols, transpose = !transpose) %>%
                           magrittr::set_colnames(svars1) %>%
                           data.frame(stringsAsFactors = FALSE, check.names = FALSE))
   }
   
   # Rm features with missing ids (happens in MaxQuant files due to empty interspersed lines)
   idx <- !is.na(fids1)
   if (any(!idx)){
      if (verbose) autonomics.support::cmessage("\t\tRm %d features with missing 'feature_id' values", sum(!idx))
      fids1  %<>% magrittr::extract(idx)
      fdata1 %<>% magrittr::extract(idx,)
      exprs1 %<>% magrittr::extract(idx, )
   }
   
   # Rm samples with missing ids
   idx <- !is.na(sids1)
   if (any(!idx)){
      if (verbose) autonomics.support::cmessage("\t\tRm %d samples with missing 'sample_id' values", sum(!idx))
      sids1  %<>% magrittr::extract(idx)
      sdata1 %<>% magrittr::extract(idx,)
      exprs1 %<>% magrittr::extract(, idx)
   }
   
   # Name features and samples
   fids1 %<>% autonomics.support::uniquify('make.unique')
   sids1 %<>% autonomics.support::uniquify('make.unique')
   rownames(exprs1) <- rownames(fdata1) <- fids1
   colnames(exprs1) <- rownames(sdata1) <- sids1
   
   # Wrap into Sumexp
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs1))
   autonomics.import::fdata(object) <- fdata1
   autonomics.import::sdata(object) <- sdata1
   
   # Return
   object
}



#=======================================================
# SOMASCAN
#=======================================================

#' @rdname read_somascan
#' @importFrom magrittr %>%
#' @export
read_somascan_asis <- function(
   file,
   fid_var = 'SeqId',
   sid_var = 'SampleId'
){
   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.types::assert_is_a_string(fid_var)
   assertive.types::assert_is_a_string(sid_var)
   
   # Peak
   dt <- data.table::fread(file, header = FALSE)
   fdata_row1 <- 1 + which(dt[[1]] == '^TABLE_BEGIN')[1]
   fvar_col <- which(dt %>% extract_dt_row(fdata_row1) != '')[1]
   svar_row <- which(dt[[1]] != '' &  (1:nrow(dt) > fdata_row1))[1]
   fvars <- dt %>% extract_dt_col(fvar_col) %>% magrittr::extract(fdata_row1:svar_row)
   svars <- dt %>% extract_dt_row(svar_row) %>% magrittr::extract(1:fvar_col-1)
   fid_row <- fdata_row1 + -1 + which(fvars == fid_var)
   sid_col <- which(svars == sid_var)
   
   # Read
   file %>% read_omics_asis(sheet      =  NULL,
                            fid_rows   =  fid_row,                   fid_cols   = (fvar_col+1):ncol(dt),
                            sid_rows   = (svar_row+1):nrow(dt),      sid_cols   =  sid_col,
                            expr_rows  = (svar_row+1):nrow(dt),      expr_cols  = (fvar_col+1):ncol(dt),
                            fvar_rows  =  fdata_row1:(svar_row-1),   fvar_cols  =  fvar_col,
                            fdata_rows =  fdata_row1:(svar_row-1),   fdata_cols = (fvar_col+1):ncol(dt),
                            svar_rows  =  svar_row,                  svar_cols  = 1:(fvar_col-1),
                            sdata_rows = (svar_row+1):nrow(dt),      sdata_cols = 1:(fvar_col-1),
                            transpose  = TRUE,
                            verbose    = TRUE)
}





#======================================================
# METABOLON
#======================================================

#' @rdname read_metabolon
#' @importFrom magrittr %>%
#' @export
read_metabolon_asis <- function(file, sheet = 2, fid_var = 'COMP_ID', sid_var = 'CLIENT_IDENTIFIER'){
   
   assertive.files::assert_all_are_existing_files(file)
   
   d_f <- readxl::read_excel(file, sheet, col_names = FALSE)
   
   fvar_rows <- which(!is.na(d_f %>% extract_dt_col(1))) %>% magrittr::extract(1)
   svar_cols <- which(!is.na(d_f %>% extract_dt_row(1))) %>% magrittr::extract(1)
   fvar_cols <- fdata_cols <- 1:svar_cols
   svar_rows <- sdata_rows <- 1:fvar_rows
   
   fid_rows  <- fdata_rows <- expr_rows <- (fvar_rows+1):nrow(d_f)
   sid_cols  <- sdata_cols <- expr_cols <- (svar_cols+1):ncol(d_f)
   fid_cols  <-  d_f %>% extract_dt_row(fvar_rows) %>% magrittr::extract(1:svar_cols) %>% magrittr::equals(fid_var) %>% which()
   sid_rows  <-  d_f %>% extract_dt_col(svar_cols) %>% magrittr::extract(1:fvar_rows) %>% magrittr::equals(sid_var) %>% which()
   
   file %>% read_omics_asis(sheet      = sheet,
                            fid_rows   = fid_rows,      fid_cols   = fid_cols,
                            sid_rows   = sid_rows,      sid_cols   = sid_cols,
                            expr_rows  = expr_rows,     expr_cols  = expr_cols,
                            fvar_rows  = fvar_rows,     fvar_cols  = fvar_cols,
                            svar_rows  = svar_rows,     svar_cols  = svar_cols,
                            fdata_rows = fdata_rows,    fdata_cols = fdata_cols,
                            sdata_rows = svar_rows,     sdata_cols = sdata_cols,
                            transpose  = FALSE,
                            verbose    = TRUE)
}


#=========================================================
# RNASEQ
#=========================================================

#' @rdname read_rnaseq
#' @importFrom magrittr %>%
#' @export
read_rnaseq_asis <- function(file, fid_var = 'gene_id'){
   
   assertive.files::assert_all_are_existing_files(file)
   
   dt <- data.table::fread(file, integer64='numeric')
   expr_cols   <- dt %>% vapply(is.integer, logical(1)) %>% unname() %>% which()
   fdata_cols  <- dt %>% vapply(is.integer, logical(1)) %>% magrittr::not() %>% unname() %>% which()
   fid_col     <- which(names(dt)==fid_var)
   
   file %>% read_omics_asis(sheet      = NULL,
                            fid_rows   = 2:nrow(dt),   fid_cols   = fid_col,
                            sid_rows   = 1,            sid_cols   = expr_cols,
                            expr_rows  = 2:nrow(dt),   expr_cols  = expr_cols,
                            fvar_rows  = 1,            fvar_cols  = fdata_cols,
                            fdata_rows = 2:nrow(dt),   fdata_cols = fdata_cols,
                            transpose  = FALSE,
                            verbose    = TRUE)
}

#' @rdname counts_to_cpm
#' @export
libsizes <- function(counts){
   colSums(counts) * edgeR::calcNormFactors(counts)
}

#' Convert between counts and cpm (counts per million scaled reads)
#' @param counts    count matrix
#' @param cpm       cpm matrix (counts per million scaled reads)
#' @param lib.size  scaled library sizes (vector)
#' @param verbose   logical(1)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    counts <- 'extdata/rnaseq/gene_counts.txt' %>%
#'               system.file(package = 'subramanian.2016') %>%
#'               read_rnaseq_asis()     %>%
#'               autonomics.import::exprs()
#'    lib.size <- counts %>% libsizes()
#'    cpm      <- counts %>% counts_to_cpm(lib.size)
#'    counts2  <- cpm    %>% cpm_to_counts(lib.size)
#'    sum(counts-counts2)
#' }
#' @return cpm matrix
#' @export
counts_to_cpm <- function(counts, lib.size = libsizes(counts), verbose = FALSE){
   if (verbose)  message('\t\tTMM normalize exprs: counts -> cpm')
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
   lib.size <- counts %>% libsizes()
   counts %>% limma::voom(design = design, lib.size = lib.size, plot = plot, ...) %>%
      magrittr::extract2('weights') %>%
      magrittr::set_rownames(rownames(object)) %>%
      magrittr::set_colnames(colnames(object))
}

#' Does object contain subgroup replicates?
#' @param object SummarizedExperiment
#' @return logical(1)
#' @importFrom magrittr %>%
#' @export
contains_replicates <- function(object){
   if (!'subgroup' %in% autonomics.import::svars(object)) return(FALSE)
   autonomics.import::sdata(object)$subgroup %>% duplicated() %>% any()
}

#' Create design for voom transformation
#' @param object SummarizedExperiment
#' @param verbose logical(1)
#' @return NULL (if no replicates) or design matrix
#' @export
create_voom_design <- function(object, verbose = TRUE){
   
   # Replicates
   if (contains_replicates(object)) return(autonomics.import::create_design_matrix(object))
   
   # No replicates
   autonomics.support::cmessage('\t\t\tsubgroup values not replicated: voom(design=NULL)')
   return(NULL)
   
}

#' Log2 transform
#' @param object SummarizedExperiment
#' @param verbose logical(1)
#' @return Updated SummarizedExperiment
#' @importFrom magrittr %<>%
#' @export
log2transform <- function(object, verbose = FALSE){
   if (verbose) message('\t\tLog2 transform')
   autonomics.import::exprs(object) %<>% log2()
   return(object)
}

#' Compute voom precision weights
#' @param object  SummarizedExperiment: exprs(.) with log2cpm, counts(.) with raw counts.
#' @param design  design matrix
#' @param plot    logical
#' @param verbose logical(1)
#' @param ...     passed to limma::voom() -> limma::lmFit()
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- 'extdata/rnaseq/gene_counts.txt'            %>%
#'               system.file(package = 'subramanian.2016')  %>%
#'               read_rnaseq_asis()
#'    autonomics.import::counts(object) <- autonomics.import::exprs(object)
#'    object$subgroup <- object %>% autonomics.import::guess_subgroup_values()
#'    object %>% compute_precision_weights()
#' }
#' @importFrom magrittr %>%
#' @export
compute_precision_weights <- function(
   object,
   design  = object %>% create_voom_design(),
   plot    = TRUE,
   verbose = FALSE
){
   
   # Assert & message
   assertive.properties::assert_is_not_null(autonomics.import::counts(object))
   if (verbose) autonomics.support::cmessage('\t\tCompute and add precision weights')
   
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

#' @importFrom magrittr %>%
explicitly_compute_precision_weights_once <- function(
   object,
   design = object %>% create_voom_design(),
   plot   = TRUE,
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
   l <- stats::lowess(mean.log2count, sdrt.log2count, f = 0.5)
   f <- stats::approxfun(l, rule = 2)
   
   # Compute precision weights
   w <- 1/f(fitted.log2count)^4     # f(.) = sqrt(sd(.)) --> f(.)^4 = var(.)
   dim(w) <- dim(fitted.log2count)
   
   # Plot
   if (plot) {
      plot(mean.log2count, sdrt.log2count, xlab = "mean log2count", ylab = "sdrt log2count",
           pch = 16, cex = 0.25)
      graphics::title("voom: Mean-variance trend")
      graphics::lines(l, col = "red")
   }
   
   # Return
   return(w)
}



#==========================================================
# EXIQON
#==========================================================

#' @rdname read_exiqon
#' @importFrom magrittr %>%
#' @export
read_exiqon_asis <- function(file){
   assertive.files::assert_all_are_existing_files(file)
   dt <- read_file(file, sheet=1)
   file %>% read_omics_asis(sheet = 1,
                            fid_rows   = 1,                       fid_cols   = 2:(ncol(dt)-2),
                            sid_rows   = 2:(nrow(dt)-3),          sid_cols   = 1,
                            expr_rows  = 2:(nrow(dt)-3),          expr_cols  = 2:(ncol(dt)-2),
                            fvar_rows  = (nrow(dt)-2):nrow(dt),   fvar_cols  = 1,
                            svar_rows  = 1,                       svar_cols  = (ncol(dt)-1):ncol(dt),
                            fdata_rows = (nrow(dt)-2):nrow(dt),   fdata_cols = 2:(ncol(dt)-2),
                            sdata_rows = 2:(nrow(dt)-3),          sdata_cols = (ncol(dt)-1):ncol(dt),
                            transpose  = TRUE,
                            verbose    = TRUE)
}

#====================================
# LCMS PROTEINGROUPS & PHOSPHOSITES
#====================================

#' maxquant patterns
#' @export
maxquant_patterns <- c(`Ratio normalized`             =  '^Ratio ([HM]/[ML]) normalized (.+)$',
                       `Ratio`                        =  '^Ratio ([HM]/[ML]) (?!count|type|variability|iso-count|normalized)(.+)',
                       `LFQ intensity`                =  '^LFQ intensity ([HML])? ?(.+)$',
                       `Reporter intensity corrected` =  '^Reporter intensity corrected ([0-9]+) (.+)$',
                       `Reporter intensity`           =  '^Reporter intensity ([0-9]+) (.+)$',
                       `Intensity labeled`            =  '^Intensity ([HML]) (.+)$',
                       `Intensity`                    =  '^Intensity (.+)$')


#' Guess maxquant quantity from snames
#'
#' charactervector, dataframe, or SummarizedExperiment.
#'
#' @param x charactervector, dataframe, or SummarizedExperiment
#' @param ... used for proper S3 method dispatch
#' @return  character(1): value from names(maxquant_patterns)
#' @examples
#' require(magrittr)
#'
#' # file
#'    x <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'          system.file(package = 'autonomics.data')
#'    guess_maxquant_quantity(x)
#'
#' # charactervector
#'     guess_maxquant_quantity("Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1")
#'     guess_maxquant_quantity("Ratio M/L STD(L)_EM00(M)_EM01(H)_R1")
#'     guess_maxquant_quantity("LFQ intensity EM00.R1")
#'     guess_maxquant_quantity("Reporter intensity corrected 0 STD(0)EM00(1)EM01(2)_R1")
#'     guess_maxquant_quantity("Reporter intensity 0 STD(0)EM00(1)EM01(2)_R1")
#'     guess_maxquant_quantity("Intensity H STD(L)_EM00(M)_EM01(H)_R1")
#'
#' # dataframe
#'     if (require(autonomics.data)){
#'       x <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data') %>%
#'             data.table::fread()
#'       guess_maxquant_quantity(x)
#'     }
#'
#' # SummarizedExperiment
#'     if (require(autonomics.data)){
#'       x <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data') %>%
#'             read_proteingroups()
#'       guess_maxquant_quantity(x)
#'     }
#' @export
guess_maxquant_quantity <- function(x, ...){
   UseMethod("guess_maxquant_quantity", x)
}

#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.character <- function(x, ...){      # x = character vector of maxquant colnames
   
   # read if x is filename
   if (assertive.files::is_existing_file(x)){
      x %<>% data.table::fread(header = TRUE, nrows = 1) %>% names()
   }
   
   # guess from character vector
   for (quantity in names(maxquant_patterns)){
      pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.data.frame <- function(x, ...){     # x = maxquant dataframe
   x <- names(x)
   for (quantity in names(maxquant_patterns)){
      pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}

#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.SummarizedExperiment <- function(x, ...){     # x = SummarizedExperiment
   x <- autonomics.import::snames(x)
   for (quantity in names(maxquant_patterns)){
      pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}




#' proteingroups fvars
#' @export
proteingroups_fvars <- c(c('id', 'Majority protein IDs', 'Protein names', 'Gene names', 'Contaminant', 'Potential contaminant', 'Reverse', 'Phospho (STY) site IDs'))

#' Clean maxquant snames
#'
#' For charactervector or SummarizedExperiment
#'
#' Drop "Ratio normalized", "LFQ intensity" etc from maxquant snames & sample_id values
#'
#' @param x        character(.) or SummarizedExperiment
#' @param quantity maxquant quantity
#' @param verbose  logical(1)
#' @param ...      allow for proper S3 method dispatch
#' @examples
#' require(magrittr)
#'
#' # character vector
#'     standardize_maxquant_snames("Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1")
#'     standardize_maxquant_snames("Ratio M/L STD(L)_EM00(M)_EM01(H)_R1")
#'     standardize_maxquant_snames('LFQ intensity STD_R1')
#'     standardize_maxquant_snames('LFQ intensity L STD(L)_EM00(M)_EM01(H)_R1')
#'     standardize_maxquant_snames('Reporter intensity 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1')
#'     standardize_maxquant_snames('Reporter intensity corrected 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1')
#'
#' # SummarizedExperiment
#' if (require(autonomics.data)){
#'      x <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'            system.file(package = 'autonomics.data')     %>%
#'            read_proteingroups(standardize_snames = FALSE, 
#'                               demultiplex_snames = FALSE)
#'      x %>% standardize_maxquant_snames(verbose = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
standardize_maxquant_snames <- function (x, ...) {
   UseMethod("standardize_maxquant_snames", x)
}


#' @importFrom magrittr %>%
#' @export
#' @rdname standardize_maxquant_snames
standardize_maxquant_snames.character <- function(
   x,
   quantity = guess_maxquant_quantity(x),
   verbose  = FALSE, 
   ...
){
   # x = mix + channel. Return mix if single channel.
   pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
   channel <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
   mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$2')
   if (all(channel=='')){ cleanx <- mix
   } else               { cleanx <- sprintf('%s{%s}', mix, channel)
   }
   autonomics.support::cmessage('\t\tClean snames: %s  ->  %s', x[1], cleanx[1])
   return(cleanx)
}

#' @importFrom magrittr %>%
#' @export
#' @rdname standardize_maxquant_snames
standardize_maxquant_snames.SummarizedExperiment <- function(
   x,
   quantity = guess_maxquant_quantity(x),
   verbose  = FALSE, 
   ...
){
   newsnames <- autonomics.import::snames(x) %>% standardize_maxquant_snames(quantity = quantity, verbose=verbose)
   autonomics.import::snames(x) <- autonomics.import::sdata(x)$sample_id <- newsnames
   x
}


#' Demultiplex snames
#'
#' For charactervector or SummarizedExperiment
#'
#' @param x        character vector or SummarizedExperiment
#' @param verbose  logical(1)
#' @examples
#' require(magrittr)
#'
#' # character vector
#'     # Alternate multiplexing forms supported
#'      demultiplex_snames("STD(L)_EM00(M)_EM01(H)_R1{M/L}")     # Label Ratio
#'      demultiplex_snames('A(0)_B(1)_C(2)_D(3)_R1{0}'     )     # Reporter intensity
#'      demultiplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')       # Label Intensity
#'
#'    # Alternate separators supported
#'      demultiplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')       # underscore
#'      demultiplex_snames('STD(L).EM00(M).EM01(H).R1{L}')       # dot
#'      demultiplex_snames('STD(L)EM00(M)EM01(H).R1{L}')         # no separator
#'
#'    # Composite snames supported
#'      demultiplex_snames("WT.t0(L)_WT.t1(M)_WT.t2(H)_R1{H/L}") # composite snames
#'
#'    # Uniqueness ensured by appending labels when necessary
#'      demultiplex_snames(c("STD(L).BM00(M).BM00(H).R10{M/L}",  # implicit uniquification
#'                                          "STD(L).BM00(M).BM00(H).R10{H/L}"))
#'    # Uniplexed snames are returned unchanged
#'      demultiplex_snames(c('STD_R1', 'EM0_R1'))
#'
#' # SummarizedExperiment
#'   if (require(autonomics.data)){
#'      x <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'            system.file(package = 'autonomics.data')     %>%
#'            read_proteingroups(standardize_snames = FALSE, demultiplex_snames = FALSE) %>%
#'            standardize_maxquant_snames()
#'      x %>% demultiplex_snames(verbose = TRUE)
#'   }
#'
#' @export
demultiplex_snames <- function (x, ...) {
   UseMethod("demultiplex_snames", x)
}

#' @rdname demultiplex_snames
#' @export
demultiplex_snames.character <- function(x, verbose = FALSE, ...){
   
   # Return unchanged if not multiplexed
   # KD(H)WT(L){H/L}
   pattern <- '(.+)\\{(.+)\\}'
   n_open   <- x %>% stringi::stri_count_fixed('(')
   n_closed <- x %>% stringi::stri_count_fixed(')')
   is_multiplexed <- all(stringi::stri_detect_regex(x, pattern) & (n_open==n_closed) & (n_open>0))
   if (!is_multiplexed) return(x)
   
   # Separate mix and channel
   mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
   channel <- x %>% stringi::stri_replace_first_regex(pattern, '$2')
   
   # Separate labels and samples
   pattern <- '\\(.+?\\)'
   labels  <- mix %>% stringi::stri_extract_all_regex(pattern) %>% lapply(stringi::stri_replace_first_fixed, '(', '') %>%
      lapply(stringi::stri_replace_first_fixed, ')', '')
   samples <- mix %>% stringi::stri_split_regex(pattern) %>%
      # rm sep from samples (but not from replicate - needed to glue back later!)
      lapply(function(y){y[1:length(labels[[1]])] %<>% stringi::stri_replace_first_regex('^[_. ]', ''); y})
   
   # Return unchanged if mixes differ in no of labels or samples
   are_all_identical <- function(y) if (length(y)==1) TRUE else all(y[-1] == y[1])
   n_samples <- vapply(samples,  length, integer(1))
   n_labels  <- vapply(labels, length, integer(1))
   if (!are_all_identical(n_samples) | !are_all_identical(n_labels)){
      autonomics.support::cmessage('\t\tCannot demultiplexing snames: mixes differ in number of samples or labels')
      return(x)
   }
   
   # Extract replicate
   n_samples %<>% unique()
   n_labels  %<>% unique()
   if (n_samples > n_labels){ replicate <- mix %>% stringi::stri_split_regex(pattern) %>% vapply((function(y) y %>% magrittr::extract(length(y))), character(1))
   samples %<>% lapply(extract, 1:(n_samples-1))
   } else {                   replicate <- rep('', length(samples))
   }
   
   # Extract channel samples from mix
   is_ratio <- channel %>% stringi::stri_detect_fixed('/') %>% all()
   samples %<>% mapply(set_names, ., labels, SIMPLIFY = FALSE)
   if (is_ratio){
      num_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(extract, character(1), 1)
      den_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(extract, character(1), 2)
      den_samples <- mapply(extract, samples, den_label)
      num_samples <- mapply(extract, samples, num_label)
      xdemultiplex <- sprintf('%s_%s%s', num_samples, den_samples, replicate)
   } else {
      samples %<>% mapply(extract, ., channel)
      xdemultiplex <- sprintf('%s%s', samples, replicate)
   }
   if (verbose) autonomics.support::cmessage('\t\tDemultiplex snames: %s  ->  %s', x[1], xdemultiplex[1])
   
   # Ensure uniqueness. Add labels if required.
   idx <- autonomics.support::cduplicated(xdemultiplex) %>% which()
   if (length(idx)>0){
      label_tags <- channel[idx] %>% stringi::stri_replace_first_fixed('/', '')
      if (verbose)   autonomics.support::cmessage('\t\tUniquify snames: %s -> %s%s (for %d/%d snames)',
                                                  xdemultiplex[idx][1], xdemultiplex[idx][1], label_tags[1],
                                                  length(idx), length(xdemultiplex))
      xdemultiplex[idx] %<>% paste0(label_tags)
   }
   
   # Return
   return(xdemultiplex)
}

#' @rdname demultiplex_snames
#' @importFrom magrittr %>%
#' @export
demultiplex_snames.SummarizedExperiment <- function(x,verbose  = FALSE, ...){
   newsnames <- autonomics.import::snames(x) %>% demultiplex_snames(verbose = verbose)
   autonomics.import::snames(x) <- autonomics.import::sdata(x)$sample_id <- newsnames
   x
}

#' Read proteingroups
#' @param file                character(1). Path to 'proteinGroups.txt'
#' @param quantity            character(1). Expression columns to extract from proteinGroups file: 'Ratio normalized', 'Ratio', 'LFQ intensity',
#'                                          'LFQ intensity', 'Reporter intensity corrected', 'Reporter intensity', 'Intensity labeled', or 'Intensity'.
#' @param fvars               character(n). Annotation columns to extract from proteinGroups file.
#' @param standardize_snames  logical(1).   Standardize maxquant snames: Ratio normalized H/L WT(L).KD(H).R1 -> WT(L).KD(H).R1{H/L} ?
#' @param demultiplex_snames  logical(1).   Demultiplex (standardized) maxquant snames: WT(L).KD(H).R1{H/L} -> KD(H)_WT(L).R1 ?
#' @param verbose             logical(1).   Message progress?
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'   file %>% read_proteingroups()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% read_proteingroups()
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_proteingroups <- function(
   file,
   quantity           = guess_maxquant_quantity(file),
   fvars              = proteingroups_fvars,
   standardize_snames = TRUE,
   demultiplex_snames = TRUE,
   verbose            = TRUE
){

   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.sets::assert_is_subset(quantity, names(maxquant_patterns))
   assertive.types::assert_is_a_bool(verbose)

   # Initial Read
   assertive.files::assert_all_are_existing_files(file)
   dt <- data.table::fread(file, integer64 = 'numeric', header = TRUE)
   fvars %<>% intersect(names(dt))

   # Define components
   fid_rows   <- 2:nrow(dt)
   fid_cols   <- which(names(dt) == 'id')
   sid_rows   <- 1
   sid_cols   <- names(dt) %>% stringi::stri_detect_regex(maxquant_patterns[[quantity]]) %>% which()
   expr_rows  <- 2:nrow(dt)
   expr_cols  <- sid_cols
   fvar_rows  <- 1
   fvar_cols  <- match(fvars, names(dt))
   fdata_rows <- 2:nrow(dt)
   fdata_cols <- fvar_cols

   # Read sumexp
   object <- file %>% read_omics_asis(fid_rows   = fid_rows,     fid_cols   = fid_cols,
                                      sid_rows   = sid_rows,     sid_cols   = sid_cols,
                                      expr_rows  = expr_rows,    expr_cols  = expr_cols,
                                      fvar_rows  = fvar_rows,    fvar_cols  = fvar_cols,
                                      fdata_rows = fdata_rows,   fdata_cols = fdata_cols,
                                      transpose  = FALSE,
                                      verbose    = verbose)

   contaminant_var <- c('Contaminant', 'Potential contaminant') %>% intersect(autonomics.import::fvars(object))
   autonomics.import::fdata(object)[[contaminant_var]] %<>% (function(x){x[is.na(x)] <- ''; x})
   autonomics.import::fdata(object)[['Reverse'      ]] %<>% (function(x){x[is.na(x)] <- ''; x})

   # Simplify snames
   if (standardize_snames) object %<>% standardize_maxquant_snames(verbose = verbose)
   if (demultiplex_snames) object %<>% demultiplex_snames(verbose = verbose)

   # Guess subgroup (to allow tweaking before prepro)
   object$subgroup <- object$sample_id %>% autonomics.import::guess_subgroup_values(verbose = verbose)
   #object$block    <- object$sample_id %>% autonomics.import::guess_subject_values( verbose = TRUE)
   
   # Return
   object

}

#' phosphosites fvars
#' @export
phosphosite_fvars <- c('id', 'Protein group IDs', 'Positions within proteins', 'Localization prob')


#' Read phosphosites
#' @param file                   character(1)
#' @param quantity               NULL or value in names(maxquant_patterns)
#' @param fvars                  character(n)
#' @examples
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemdiff/maxquant/phospho (STY)Sites.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% read_phosphosites()
#' }
#' @importFrom magrittr %>%
#' @export
read_phosphosites <- function(
   file,
   quantity = guess_maxquant_quantity(file),
   fvars    = phosphosite_fvars,
   verbose  = FALSE
){

   # Check
   `Protein group IDs` <- `Localization prob` <- NULL

   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.sets::assert_is_subset(quantity, names(maxquant_patterns))
   assertive.types::assert_is_character(fvars)

   # Initial Read
   dt <- data.table::fread(file, integer64 = 'numeric', header = TRUE)
   pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
   value_cols <- names(dt) %>% (function(x) stringi::stri_detect_regex(x,pattern) & !stringi::stri_detect_regex(x, '___[1-3]')) %>% which()
   fvar_cols  <- which(names(dt) %in% fvars)
   #dt %<>% magrittr::extract(, c(fvars, value_cols), with = FALSE)

   # Read phosphosites
   fid_rows   <- 2:nrow(dt)
   fid_cols   <- which(names(dt) == 'id')
   sid_rows   <- 1
   sid_cols   <- value_cols
   expr_rows  <- 2:nrow(dt)
   expr_cols  <- value_cols
   fvar_rows  <- 1
   fvar_cols  <- fvar_cols
   fdata_rows <- 2:nrow(dt)
   fdata_cols <- fvar_cols
   phosphosites  <- file  %>% read_omics_asis(fid_rows   = fid_rows,     fid_cols   = fid_cols,
                                                                 sid_rows   = sid_rows,     sid_cols   = sid_cols,
                                                                 expr_rows  = expr_rows,    expr_cols  = expr_cols,
                                                                 fvar_rows  = fvar_rows,    fvar_cols  = fvar_cols,
                                                                 fdata_rows = fdata_rows,   fdata_cols = fdata_cols,
                                                                 transpose  = FALSE,
                                                                 verbose    = verbose)

   # Calculate occupancies (allows to disentangle phosphorylation and protein expression)
   autonomics.support::cmessage('\t\toccupancies(phosphosites) = exprs(phosphosites) - exprs(proteingroups)')
   proteingroups <- read_proteingroups(proteingroups_file, quantity = quantity, verbose = FALSE) %>%
                    magrittr::extract(phosphosites %>% autonomics.import::fvalues("Protein group IDs"), ) %>%
                    magrittr::extract(, phosphosites %>% autonomics.import::snames())
   autonomics.import::occupancies(phosphosites) <- autonomics.import::exprs(phosphosites) - autonomics.import::exprs(proteingroups)


   # Simplify snames
   if (standardize_snames) phosphosites %<>% standardize_maxquant_snames(verbose = verbose)
   if (demultiplex_snames) phosphosites %<>% demultiplex_snames(verbose = verbose)

   # Return
   phosphosites

}

