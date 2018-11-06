# require(magrittr)
#
# # Metabolon
# file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
# sheet <- 2
# fdata_range <- c(10, 401,  1, 14)
# sdata_range <- c(1,   10, 14, 86)
# exprs_range <- c(11, 401, 15, 86)
# sid_var  <- 'CLIENT_IDENTIFIER'
# fid_var <- 'COMP_ID'
# features_in_rows <- TRUE
#
# # Somascan
# file <- system.file('extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat', package = 'atkin.2014')
# fdata_range <- c(73,  91,  26, 1343)
# sdata_range <- c(92, 520,  1,    25)
# exprs_range <- c(93, 520, 1343,  27)
# transpose <- TRUE
# sid_var <- 'SampleId'
# fid_var <- 'SeqId'

# Extract row
extract_dt_row <- function(dt, i) dt %>% magrittr::extract(i,) %>% as.matrix() %>% magrittr::extract(1, ) %>% unname()
extract_dt_col <- function(dt, i) dt[[i]]

#' Read rectangle from excel or txt file
#' @param file       character(1)
#' @param sheet      numeric(1), character(1), or NULL
#' @param range      integer(4): c(row1, rown, col1, coln)
#' @param transpose  logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% read_matrix(2, c(11, 401, 15, 86))                   %>% extract(1:3, 1:3)           # exprs
#'    file %>% read_matrix(2, c(11, 401,  5,  5))                   %>% extract(1:3, , drop=FALSE)  # fids
#'    file %>% read_matrix(2, c( 2,   2, 15, 86))                   %>% extract(, 1:3, drop=FALSE)  # sids
#'    file %>% read_matrix(2, c(10, 401,  1, 14))                   %>% extract(1:3, 1:3)           # fdata
#'    file %>% read_matrix(2, c(1,   10, 14, 86), transpose = TRUE) %>% extract(1:3, 1:3)           # sdata
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_matrix <- function(file, sheet, range, transpose = FALSE){
   if (tools::file_ext(file) %in% c('xls', 'xlsx')){
      xl_range <- sprintf('R%dC%d:R%dC%d', range[[1]], range[[3]], range[[2]], range[[4]])
      dt <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE, range = xl_range) %>%
               data.table::data.table()
   } else {
      skip <- range[[1]] - 1
      nrows <- range[[2]] - range[[1]] + 1
      select <- range[[3]]:range[[4]]
      dt <- file %>% data.table::fread(na.strings = "", header = FALSE, skip = skip, nrows = nrows, select = select, integer64 = 'numeric')
   }
   dt %<>% as.matrix()
   if (transpose) dt %<>% t()
   dt
}


#' Read dataframe block from excel or txt file
#' @param file       character(1)
#' @param sheet      numeric(1), character(1), or NULL
#' @param range      integer(4): c(row1, rown, col1, coln)
#' @param transpose  logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% read_dataframe(2, c(10, 401,  1, 14))                   %>% extract(1:3, 1:3)           # fdata
#'    file %>% read_dataframe(2, c(1,   10, 14, 86), transpose = TRUE) %>% extract(1:3, 1:3)           # sdata
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_dataframe <- function(file, sheet, range, transpose = FALSE){
   matrix1 <- file %>% read_matrix(sheet, range, transpose)
   if (nrow(matrix1)>1){
      vars <- matrix1[1, , drop = FALSE]
      matrix1 %<>% magrittr::extract(-1, , drop = FALSE) %>%
                   magrittr::set_colnames(vars)
   }
   matrix1 %>% data.frame(stringsAsFactors = FALSE, check.names = FALSE)
}



#' Read omics data from txt or xls file
#' @param file         character(1)
#' @param sheet        integer(1) or character(1)
#' @param exprs_range  numeric(4): c(row1, rown, col1, coln)
#' @param fid_range    numeric(4): c(row1, rown, col1, coln)
#' @param sid_range    numeric(4): c(row1, rown, col1, coln)
#' @param fdata_range  numeric(4): c(row1, rown, col1, coln)
#' @param sdata_range  numeric(4): c(row1, rown, col1, coln)
#' @param transpose    logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'   file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'   file %>% autonomics.import::read_omics(sheet       = 2,
#'                                          exprs_range = c(11, 401, 15, 86),
#'                                          fid_range   = c(11, 401,  5,  5),
#'                                          sid_range   = c( 3,   3, 15, 86),
#'                                          fdata_range = c(10, 401,  1, 14),
#'                                          sdata_range = c(1,   10, 14, 86),
#'                                          transpose   = FALSE)
#' }
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    fdata_range <- c(1, 48821, 1, 4)
#'    sdata_range <- c(1, 1, 5, 86)
#'    exprs_range <- c(2, 48821, 5, 86)
#'    transpose <- FALSE
#'    sid_var <- NULL
#'    fid_var <- 'gene_id'
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_omics <- function(
   file,
   sheet,
   exprs_range,
   fid_range   = NULL,
   sid_range   = NULL,
   fdata_range = NULL,
   sdata_range = NULL,
   transpose   = FALSE
){
   # Assert
   assertive.types::assert_is_numeric(exprs_range);  assertive.properties::assert_is_of_length(exprs_range, 4)
   assertive.types::assert_is_numeric(fid_range);    assertive.properties::assert_is_of_length(fid_range,   4)
   assertive.types::assert_is_numeric(sid_range);    assertive.properties::assert_is_of_length(sid_range,   4)
   assertive.types::assert_is_numeric(fdata_range);  assertive.properties::assert_is_of_length(fdata_range, 4)
   assertive.types::assert_is_numeric(sdata_range);  assertive.properties::assert_is_of_length(sdata_range, 4)
   assertive.types::assert_is_a_bool(transpose)

   # Read
   exprs1 <- file %>% autonomics.import::read_matrix(sheet = sheet, range = exprs_range, transpose =  transpose)
   fids   <- file %>% autonomics.import::read_matrix(sheet = sheet, range = fid_range,   transpose =  transpose) %>% as.vector('character')
   sids   <- file %>% autonomics.import::read_matrix(sheet = sheet, range = sid_range,   transpose =  transpose) %>% as.vector('character')
   fdata1 <- file %>% read_dataframe(sheet = sheet, range = fdata_range, transpose =  transpose)
   sdata1 <- file %>% read_dataframe(sheet = sheet, range = sdata_range, transpose = !transpose)
   rownames(exprs1) <- rownames(fdata1) <- fdata1$feature_id <- fids %>% autonomics.support::uniquify(method = 'make.unique', verbose = TRUE)
   colnames(exprs1) <- rownames(sdata1) <- sdata1$sample_id  <- sids %>% autonomics.support::uniquify(method = 'make.unique', verbose = TRUE)
   fdata1 %<>% autonomics.support::pull_columns('feature_id')
   sdata1 %<>% autonomics.support::pull_columns('sample_id' )

   # Pack into SumExp
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs1))
   autonomics.import::sdata(object) <- sdata1
   autonomics.import::fdata(object) <- fdata1
   autonomics.import::exprs(object) <- exprs1
   object
}

#' Read SOMAScan data
#' @param file character(1)
#' @param fid_var  character(1): feature id variable
#' @param sid_var  character(1): sample id variable
#' @return Summarizedexperiment
#' @examples
#' require(magrittr)
#' file <- system.file('extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat',
#'                      package = 'atkin.2014')
#' if (file.exists(file)){
#'    file %>% autonomics.import::read_somascan()
#' }
#' @importFrom magrittr %>%
#' @export
read_somascan <- function(file, fid_var = 'SeqId', sid_var = 'SampleId'){

   dt <- data.table::fread(file, header = FALSE)
   fdata_row <- 1 + which(dt[[1]] == '^TABLE_BEGIN')[1]
   fdata_col <- which(dt %>% extract_dt_row(fdata_row) != '')[1]
   sdata_row <- which(dt[[1]] != '' &  (1:nrow(dt) > fdata_row))[1]
   fvars <- dt %>% extract_dt_col(fdata_col) %>% magrittr::extract(fdata_row:sdata_row)
   svars <- dt %>% extract_dt_row(sdata_row) %>% magrittr::extract(1:fdata_col-1)
   fid_row <- fdata_row + -1 + which(fvars == fid_var)
   sid_col <- which(svars == sid_var)

   file %>% read_omics(sheet       = NULL,
                       exprs_range = c(sdata_row+1,  nrow(dt),       fdata_col+1,   ncol(dt)),
                       fid_range   = c(fid_row,      fid_row,        fdata_col+1,   ncol(dt)),
                       sid_range   = c(sdata_row+1,  nrow(dt),       sid_col,        sid_col),
                       fdata_range = c(fdata_row,    sdata_row-1,   fdata_col,     ncol(dt)),
                       sdata_range = c(sdata_row,    nrow(dt),       1,              fdata_col-1),
                       transpose   = TRUE)

}

#' Read data from metabolon file
#' @param file     character(1)
#' @param sheet    numeric(1) or character(1)
#' @param fid_var  character(1): feature id variable
#' @param sid_var  character(1): sample id variable
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% autonomics.import::read_metabolon()
#' }
#' @importFrom magrittr %>%
#' @export
read_metabolon <- function(file, sheet = 2, fid_var = 'COMP_ID', sid_var = 'CLIENT_IDENTIFIER'){

   d_f <- readxl::read_excel(file, sheet, col_names = FALSE)
   fdata_row <- which(!is.na(d_f %>% extract_dt_col(1))) %>% magrittr::extract(1)
   sdata_col <- which(!is.na(d_f %>% extract_dt_row(1)))        %>% magrittr::extract(1)
   fid_col   <- d_f %>% extract_dt_row(fdata_row) %>% magrittr::extract(1:sdata_col) %>% magrittr::equals(fid_var) %>% which()
   sid_row   <- d_f %>% extract_dt_col(sdata_col) %>% magrittr::extract(1:fdata_row) %>% magrittr::equals(sid_var) %>% which()

   file %>% read_omics(sheet       = sheet,
                       exprs_range = c(fdata_row+1,  nrow(d_f),  sdata_col+1,  ncol(d_f)),
                       fid_range   = c(fdata_row+1,  nrow(d_f),  fid_col,      fid_col),
                       sid_range   = c(sid_row,      sid_row,    sdata_col+1,  ncol(d_f)),
                       fdata_range = c(fdata_row,    nrow(d_f),  1,            sdata_col),
                       sdata_range = c(1,            fdata_row,  sdata_col,    ncol(d_f)),
                       transpose   = FALSE)
}

#' Read data from rnaseq file
#' @param file     character(1)
#' @param fid_var  character(1): feature id variable
#' @examples
#' file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
read_rnaseq <- function(file, fid_var = 'gene_id'){

   dt <- data.table::fread(file, integer64='numeric')
   expr_cols   <- dt %>% vapply(is.integer, logical(1)) %>% unname() %>% which()
   fid_col     <- which(names(dt)==fid_var)
   exprs_range <- c(2, nrow(dt)+1, min(expr_cols), max(expr_cols)) # Note: 2 and nrow(dt)+1 because dt is loaded WITH header
   fid_range   <- c(2, nrow(dt), fid_col, fid_col)
   sid_range   <- c(1, 1, min(expr_cols), max(expr_cols))

   file %>% read_omics(sheet       = NULL,
                       exprs_range = exprs_range,
                       fid_range   = fid_range,
                       sid_range   = sid_range)

}

