#' Transcript to Gene Annotations
#'
#' @note Doesn't attempt to strip transcript versions.
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `data.frame`.
#' @export
#'
#' @examples
#' readTx2gene("http://bcbiobase.seq.cloud/tx2gene.csv")
readTx2gene <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file)
    data <- read_csv(file, col_names = c("txID", "geneID"))
    data <- as.data.frame(data)
    assert_has_no_duplicates(data[["txID"]])
    rownames(data) <- data[["txID"]]
    assertIsTx2gene(data)
    data
}
