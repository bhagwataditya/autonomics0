#' Deduplicate a vector
#' 
#' Deduplicates a vector, with a warning.
#' @param x A vector.
#' @return \code{x}, with duplicates removed.
#' @seealso \code{\link[base]{unique}}, \code{\link[base]{duplicated}}, 
#' \code{\link[assertive.properties]{has_no_duplicates}}
#' @examples
#' xx <- c(1, 2, 3, 4, 2, 3)
#' dedupe(xx)
#' @noRd
dedupe <- function(x)
{
  dupes <- duplicated(x)
  if(any(dupes))
  {
    dupe_index <- which(dupes)
    first_dupe_index <- utils::head(dupe_index, 10)
    dupe_data <- data.frame(
      Position = first_dupe_index,
      Value = x[first_dupe_index]
    )
    xname <- assertive.base::get_name_in_parent(x)
    warning(
      sprintf(
        "Removing %d duplicate values from %s:\n%s", 
        length(dupe_index), 
        xname, 
        print_and_capture(dupe_data)
      )
    )
    x[!dupes]
  } else 
  {
    x
  }
}

# Borrowed from assertive.base
#' @importFrom utils capture.output
print_and_capture <- function(x, ...)
{
  # call to enc2utf8 is a workaround for
  # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=16539
  enc2utf8(paste(capture.output(print(x, ...)), collapse = "\n"))
}