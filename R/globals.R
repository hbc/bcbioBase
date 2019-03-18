globalVariables(".")

packageVersion <- packageVersion("bcbioBase")

#' Cache URL
#' @export
#' @examples
#' bcbioBaseCacheURL
bcbioBaseCacheURL <- paste0(
    "http://bcbiobase.seq.cloud/",
    "v", packageVersion$major, ".", packageVersion$minor  # nolint
)

#' Project directory grep pattern
#' @export
#' @examples
#' projectDirPattern
projectDirPattern <- "^([[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2})_([^/]+)$"
