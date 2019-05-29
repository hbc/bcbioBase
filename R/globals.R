globalVariables(".")

packageVersion <- packageVersion("bcbioBase")

#' bcbioBase test data URL
#' @export
#' @examples
#' bcbioBaseTestsURL
bcbioBaseTestsURL <- paste0(
    "http://tests.acidgenomics.com/bcbioBase/",
    "v", packageVersion$major, ".", packageVersion$minor  # nolint
)

#' Project directory grep pattern
#' @export
#' @examples
#' projectDirPattern
projectDirPattern <- "^([[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2})_([^/]+)$"
