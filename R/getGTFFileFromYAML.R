#' Get GTF File from YAML
#'
#' @author Michael Steinbaugh
#' @inheritParams basejump::params
#' @inheritParams params
#' @export
#'
#' @return `string`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "summary.yaml")
#' yaml <- import(file)
#' x <- getGTFFileFromYAML(yaml)
#' print(x)
getGTFFileFromYAML <- function(yaml) {
    .assertIsSummaryYAML(yaml)
    # Assume all samples are using the same GTF file.
    file <- yaml %>%
        .[["samples"]] %>%
        .[[1L]] %>%
        .[["genome_resources"]] %>%
        .[["rnaseq"]] %>%
        .[["transcripts"]]
    assert_is_a_string(file)
    assert_are_identical(basename(file), "ref-transcripts.gtf")
    file
}
