#' Get quality control metrics from YAML
#'
#' Parse the summary YAML and return quality control metrics.
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#'
#' @author Michael Steinbaugh
#' @note Updated 2020-01-17.
#' @export
#'
#' @inheritParams acidroxygen::params
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
    cli_alert("Getting sample quality control metrics from YAML.")
    data <- .sampleYAML(yaml, keys = c("summary", "metrics"))
    ## Early return on empty metrics (e.g. fast mode).
    if (!hasLength(data)) {
        ## nocov start
        cli_alert_warning("No quality control metrics were calculated.")
        return(NULL)
        ## nocov end
    }
    ## Drop any metadata columns. Note we're also dropping the duplicate `name`
    ## column present in the metrics YAML.
    yamlFlatCols <- c("description", "genome_build", "sam_ref")
    blacklist <- c(camelCase(yamlFlatCols), "name")
    ## Drop blacklisted columns from the return.
    data <- data[, sort(setdiff(colnames(data), blacklist)), drop = FALSE]
    data
}
