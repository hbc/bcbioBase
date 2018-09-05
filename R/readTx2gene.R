#' Transcript to Gene Annotations
#'
#' @note Doesn't attempt to strip transcript versions.
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `DataFrame`.
#' @export
#'
#' @examples
#' x <- readTx2gene("http://bcbiobase.seq.cloud/tx2gene.csv")
#' print(x)
readTx2gene <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file)
    data <- read_csv(file, col_names = c("transcriptID", "geneID"))
    # Arrange by transcript and check for duplicates.
    data <- arrange(data, !!sym("transcriptID"))
    assert_has_no_duplicates(data[["transcriptID"]])
    # Coerce to DataFrame.
    data <- as(data, "DataFrame")
    rownames(data) <- data[["transcriptID"]]
    assertIsTx2gene(data)
    data
}
