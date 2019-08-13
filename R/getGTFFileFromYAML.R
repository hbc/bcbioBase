#' Get GTF file path from YAML
#'
#' @author Michael Steinbaugh
#' @note Updated 2019-08-05.
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
    file <- yaml %>%
        .[["samples"]] %>%
        .[[1L]] %>%
        .[["genome_resources"]] %>%
        .[["rnaseq"]] %>%
        .[["transcripts"]]
    if (!isString(file)) {
        warning("bcbio GTF file is not defined in YAML.")
        return(NULL)
    }
    message(sprintf("bcbio GTF file: %s.", file))
    if (!file.exists(file)) {
        message("bcbio GTF file is not accessible.")
        NULL
    } else {
        file
    }
}
