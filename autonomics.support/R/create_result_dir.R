utils::globalVariables('.')

#' Create result dir
#' 
#' Given projectdir 'root/analysis/collaborator/project', 
#' create resultdir 'root/results/collaborator/project/subdir', 
#' and return that path.
#' 
#' @param projectdir   project directory
#' @param subdir        results subdir
#' @param showWarnings  passed to dir.create (logical)
#' @importFrom magrittr  %<>%   %>%
#' @export
#' @return path to result dir
create_result_dir <- function(projectdir = getwd(), subdir = '', showWarnings = TRUE){
  project      <- projectdir %>% strsplit('/') %>% unlist() %>%  magrittr::extract(length(.))
  collaborator <- projectdir %>% strsplit('/') %>% unlist() %>%  magrittr::extract(length(.) - 1)
  platform     <- projectdir %>% strsplit('/') %>% unlist() %>%  magrittr::extract(length(.) - 3)
  root         <- projectdir %>% strsplit('/') %>% unlist() %>%  magrittr::extract(1:(length(.)-4)) %>% paste0(collapse='/')
  resultdir   <- file.path(root, platform, 'results', collaborator, project)
  if (subdir!='') resultdir %<>% paste0('/', subdir)
  dir.create(resultdir, showWarnings = showWarnings, recursive = TRUE)
  resultdir
}
