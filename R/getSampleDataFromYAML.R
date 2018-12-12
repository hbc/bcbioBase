#' Get Sample Data from YAML
#'
#' @author Michael Steinbaugh
#' @inheritParams basejump::params
#' @inheritParams params
#' @export
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "summary.yaml")
#' yaml <- import(file)
#' x <- getSampleDataFromYAML(yaml)
#' summary(x)
#' colnames(x)
getSampleDataFromYAML <- function(yaml) {
    message("Getting sample metadata from YAML.")
    data <- .sampleYAML(yaml, keys = "metadata")
    assert(.isSampleData(data))
    .makeSampleData(data)
}
