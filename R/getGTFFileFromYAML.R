#' Get GTF file path from YAML
#'
#' @author Michael Steinbaugh
#' @note Updated 2019-08-20.
#' @export
#'
#' @inheritParams acidroxygen::params
#'
#' @return `character(1)` or `NULL`.
#' File path if the file exists. `NULL` if the file does not exist.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "summary.yaml")
#' yaml <- basejump::import(file)
#' x <- getGTFFileFromYAML(yaml)
#' print(x)
getGTFFileFromYAML <- function(yaml) {
    assert(
        is.list(yaml),
        .isSummaryYAML(yaml)
    )
    ## Assume all samples are using the same GTF file.
    x <- yaml
    x <- x[["samples"]]
    x <- x[[1L]]
    x <- x[["genome_resources"]]
    x <- x[["rnaseq"]]
    x <- x[["transcripts"]]
    if (!isString(x)) {
        warning("bcbio GTF file is not defined in YAML.")
        return(NULL)
    }
    message(sprintf("bcbio GTF file: %s.", x))
    if (!file.exists(x)) {
        message("bcbio GTF file is not accessible.")
        return(NULL)
    }
    x
}
