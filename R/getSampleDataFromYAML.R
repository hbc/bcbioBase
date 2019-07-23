#' Get sample data from YAML
#'
#' @author Michael Steinbaugh
#' @export
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "summary.yaml")
#' yaml <- basejump::import(file)
#' x <- getSampleDataFromYAML(yaml)
#' summary(x)
#' colnames(x)

## Updated 2019-07-23.
getSampleDataFromYAML <- function(yaml) {
    assert(is.list(yaml))
    message("Getting sample metadata from YAML.")
    data <- .sampleYAML(yaml, keys = "metadata")
    assert(.isSampleData(data))
    .makeSampleData(data)
}
