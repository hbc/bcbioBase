#' Get quality control metrics from YAML
#'
#' Parse the summary YAML and return quality control metrics.
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#'
#' @author Michael Steinbaugh
#' @note Updated 2021-02-26.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `DataFrame`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "summary.yaml")
#' yaml <- basejump::import(file)
#' x <- getMetricsFromYAML(yaml)
#' summary(x)
#' colnames(x)
getMetricsFromYAML <- function(yaml) {
    assert(is.list(yaml))
    alert("Getting sample quality control metrics from YAML.")
    data <- .sampleYAML(yaml, keys = c("summary", "metrics"))
    ## Early return on empty metrics (e.g. fast mode).
    if (!hasLength(data)) {
        ## nocov start
        alertWarning("No quality control metrics were calculated.")
        return(NULL)
        ## nocov end
    }
    ## Drop any metadata columns. Note we're also dropping the duplicate `name`
    ## column present in the metrics YAML.
    yamlFlatCols <- c("description", "genome_build", "sam_ref")
    denylist <- camelCase(c(yamlFlatCols, "name"), strict = TRUE)
    ## Drop denylisted columns from the return.
    data <- data[, sort(setdiff(colnames(data), denylist)), drop = FALSE]
    data
}
