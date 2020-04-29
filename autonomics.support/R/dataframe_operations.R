#' Characterify a dataframe
#' 
#' Transform numeric variables into character variables in scientific notation
#' Limit length of character variables
#' 
#' @param x             dataframe
#' @param char_length   allowed length of strings
#' @return compactified dataframe
#' @examples characterify(data.frame(a = 1.09736, b = 4.902903))
#' @importFrom magrittr          %>%
#' @export
characterify <- function(x, char_length = 40){
  
    assertive.types::assert_is_data.frame(x)
    assertive.types::assert_is_numeric(char_length)
  
    # limit length of character variables
    selector <- x %>% vapply(is.character, logical(1))
    x[selector] %<>% apply(2, substring, 1, char_length)
  
    # numeric -> characters in scientific notation
    selector <- x %>% vapply(is.numeric, logical(1))
    x[selector] %<>% signif(2) %>% apply(2, function(x){sprintf('%1.2e', x)})
    x
}


#' Left join while preserving rownames
#' 
#' Convenient left_join function that preserves row names
#' and avoids warning messages when by values are of different type.
#' 
#' @param x first object
#' @param y second object
#' @param by by
#' @param ... passed to dplyr::left_join
#' @return merged dataframe
#' @seealso dplyr::left_join
#' @examples 
#' x <- data.frame(organism = c('Homo sapiens', 'Mus musculus'), 
#'         commonname   = c('human', 'mouse'), row.names = c('Hs', 'Mm'))
#' y <- data.frame(organism = c('Homo sapiens', 'Mus musculus'),
#'         chromosomes = c(23, 20), row.names = c('Hs', 'Mm'))
#' left_join_keeping_rownames(x, y, by = 'organism')
#' @importFrom magrittr %>% 
#' @export
left_join_keeping_rownames <- function(x, y, by, ...){
  
    # Avoid that by values are of different type
    # This throws a warning in dplyr::left_join
    x[[by]] %<>% as.character()
    y[[by]] %<>% as.character()
  
    # Left join  
    dplyr::left_join(x,y,by,...) %>% 
    magrittr::set_rownames(rownames(x))
}


#' Mutate while preserving rownames
#' @param .data dataframe
#' @param ... passed to dplyr::mutate
#' @return data.frame
#' @examples 
#' .data <- data.frame(organism = c('Homo sapiens', 'Mus musculus'), 
#'         commonname   = c('human', 'mouse'), row.names = c('Hs', 'Mm'))
#' .data %>% mutate_keeping_rownames(genus = rownames(.))
#' @importFrom magrittr  %>% 
#' @export
mutate_keeping_rownames <- function(.data, ...){
    dplyr::mutate(.data, ...) %>% 
    magrittr::set_rownames(rownames(.data))
}

#' Dedupe varnames
#' @param x list
#' @examples
#' require(magrittr)
#' x <- cbind(data.frame(a=1:3, b=1:3, c=1:3), data.frame(a=4:6))
#' x
#' x %>% dedupe_varnames()
#' @importFrom magrittr %>% 
#' @export
dedupe_varnames <- function(x){
    selector <- duplicated(names(x))
    if (any(selector)){
        cmessage('Rm repeated column "%s"', names(x)[selector] %>% 
            paste0(collapse = ', '))
        x %<>% magrittr::extract(!selector)
    }
    x
}

#' Pull columns in a dataframe to the front
#' @param df         data.frame
#' @param first_cols character vector: columns to be pulled to the front
#' @param verbose    TRUE (default) or FALSE
#' @return dataframe with re-ordered columns
#' @examples
#' require(magrittr)
#' df <- data.frame(
#'    symbol = c('A1BG', 'A2M'), 
#'    id     = c('1',    '2'),
#'    name   = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'), 
#'    type   = c('proteincoding', 'proteincoding'))
#' first_cols <- c('id', 'symbol', 'location', 'uniprot') 
#' df %>% autonomics.support::pull_columns(first_cols)
#' @importFrom magrittr %>% 
#' @export
pull_columns <- function(df, first_cols, verbose = TRUE){
  
    assertive.types::assert_is_data.frame(df)
    assertive.types::assert_is_character(first_cols)
    extract <- magrittr::extract
  
    idx <- first_cols %in% names(df)
    if (any(!idx)){
        if (verbose) cmessage(
            'pull_columns: ignore absent columns %s', 
            paste0(sprintf("'%s'", first_cols[!idx]), collapse = ', '))
        first_cols %<>% magrittr::extract(idx)
    }
  
  df %>% extract(, c(first_cols, setdiff(names(df), first_cols)), drop = FALSE)
}


#' Rm columns with only nas
#' @param df dataframe
#' @return dataframe with re-ordered columns
#' @examples 
#' require(magrittr)
#' df <- data.frame(
#'    symbol    = c('A1BG', 'A2M'), 
#'    id        = c('1',    '2'),
#'    name      = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'), 
#'    relevance = c(NA_character_, NA_character_),
#'    version   = c('NA', 'NA'), 
#'    type      = c('proteincoding', 'proteincoding'))
#' df %>% rm_na_columns()
#' @importFrom magrittr %>% 
#' @export
rm_na_columns <- function(df){
  Filter(function(x) !all(is.na(x)|x=='NA'), df) # 
}


#'Rm single value columns
#'@param df dataframe 
#'@return dataframe with informative columns
#'@examples
#' require(magrittr)
#' df <- data.frame(
#'    symbol    = c('A1BG', 'A2M'), 
#'    id        = c('1',    '2'),
#'    name      = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'), 
#'    relevance = c(NA_character_, NA_character_),
#'    type      = c('proteincoding', 'proteincoding'))
#' df %>% rm_single_value_columns()
#' @importFrom magrittr %>% 
#' @export
rm_single_value_columns <- function(df){
  Filter(function(x) length(unique(x))>1, df)
}


#' Matrixify datatable or dataframe
#' @param df dataframe or datatable (first column contains rownames)
#' @return matrix
#' @examples 
#' df <- data.frame(organims = c('Homo sapiens', 'Mus musculus'), 
#'                 chromosomes = c(23, 20))
#' matrixify(df)
#' @importFrom magrittr %>% 
#' @export
matrixify <- function(df){
  rownames1 <- df[[1]]
  df[[1]] <- NULL
  mat <- df %>% as.matrix() #data.matrix()
  rownames(mat) <- rownames1
  mat
}

