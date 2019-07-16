#' Get time stamp
#' @export
get_time_stamp <- function(){
   sprintf('%04d.%02d.%02d.%02dh.%02dm.%02ds', 
           lubridate::now() %>% lubridate::year(),
           lubridate::now() %>% lubridate::month(),
           lubridate::now() %>% lubridate::day(),
           lubridate::now() %>% lubridate::hour(),
           lubridate::now() %>% lubridate::minute(),
           lubridate::now() %>% lubridate::second() %>% round()
   )
}

#' Zip dir
#' @param dir             directory which needs to be zipped
#' @param add_time_stamp  whether to add a time stamp (logical)
#' @examples 
#' require(magrittr)
#' my_dir <- tempdir() %>% paste0('/test_zip')
#' dir.create(my_dir)
#' my_dir %>% paste0('/dirA') %>%  dir.create()
#' my_dir %>% paste0('/dirB') %>%  dir.create()
#' data.frame(a=1,b=2) %>% print2txt(my_dir %>% paste0('/dirA/fileA.txt'))
#' data.frame(a=3,b=4) %>% print2txt(my_dir %>% paste0('/dirB/fileB.txt'))
#' zip_dir(my_dir)
#' @export
zip_dir <- function(dir, add_time_stamp = TRUE){
  olddir <- getwd()
  
  rootdir <- dirname(dir)
  subdir <- basename(dir)
  setwd(dir)
  
  file_name <- paste0('../', subdir)
  if (add_time_stamp)   file_name %<>% paste0('.', get_time_stamp())
  file_name %<>% paste0('.zip')
  
  utils::zip(zipfile = file_name, files = list.files())
  setwd(olddir)
}