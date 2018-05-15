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