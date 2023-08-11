.version <- packageVersion(packageName())

#' bcbioBase test data URL
#' @export
#' @examples
#' bcbioBaseTestsURL
bcbioBaseTestsURL <- paste0(
    "https://r.acidgenomics.com/testdata/bcbiobase/",
    "v", .version$major, ".", .version$minor # nolint
)

#' Project directory grep pattern
#' @export
#' @examples
#' projectDirPattern
projectDirPattern <- "^([[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2})_([^/]+)$"
