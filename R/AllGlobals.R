globalVariables(".")

.version <- packageVersion("bcbioBase")

#' bcbioBase test data URL
#' @export
#' @examples
#' bcbioBaseTestsURL
bcbioBaseTestsURL <- paste0(
    "https://tests.acidgenomics.com/bcbioBase/",
    "v", .version$major, ".", .version$minor  # nolint
)

#' Project directory grep pattern
#' @export
#' @examples
#' projectDirPattern
projectDirPattern <- "^([[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2})_([^/]+)$"
