#' Get GTF file path from YAML
#'
#' @author Michael Steinbaugh
#' @note Updated 2022-03-07.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return `character(1)` or `NULL`.
#' File path if the file exists. `NULL` if the file does not exist.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsUrl, "summary.yaml")
#' yaml <- import(file)
#' x <- getGtfFileFromYaml(yaml)
#' print(x)
getGtfFileFromYaml <- function(yaml) {
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
        alertWarning("bcbio GTF file is not defined in YAML.")
        return(NULL)
    }
    dl(c("bcbio GTF file" = x))
    if (!file.exists(x)) {
        alertWarning("bcbio GTF file is not accessible.")
        return(NULL)
    }
    x
}
