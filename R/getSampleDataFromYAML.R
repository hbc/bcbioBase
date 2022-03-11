#' Get sample data from YAML
#'
#' @author Michael Steinbaugh
#' @note Updated 2022-03-07.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "summary.yaml")
#' yaml <- import(file)
#' x <- getSampleDataFromYAML(yaml)
#' summary(x)
#' colnames(x)
getSampleDataFromYAML <- function(yaml) {
    assert(is.list(yaml))
    alert("Getting sample metadata from YAML.")
    data <- .sampleYAML(yaml, keys = "metadata")
    makeSampleData(data)
}
