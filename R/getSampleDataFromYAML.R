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
#' file <- file.path(bcbioBaseCacheURL, "summary.yaml")
#' yaml <- basejump::import(file)
#' x <- getSampleDataFromYAML(yaml)
#' summary(x)
#' colnames(x)
getSampleDataFromYAML <- function(yaml) {
    message("Getting sample metadata from YAML.")
    data <- .sampleYAML(yaml, keys = "metadata")
    assert(.isSampleData(data))
    .makeSampleData(data)
}
