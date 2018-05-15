
#' Save session info
#' @param dir directory 
#' @importFrom utils capture.output   sessionInfo
#' @export
save_session_info <- function(dir){
  file_name <- sprintf("%s/sessionInfo.txt", dir)
  writeLines(capture.output(sessionInfo()), file_name)
}