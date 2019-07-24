#' Read transcript-to-gene annotations
#'
#' Generates a `Tx2Gene` object containing `transcriptID` and `geneID` columns.
#'
#' @note Doesn't attempt to strip transcript versions.
#'
#' @author Michael Steinbaugh
#' @export
#' @inheritParams basejump::params
#'
#' @param organism `character(1)` or `NULL`.
#'   Full Latin organism name (e.g. `"Homo sapiens"`).
#' @param genomeBuild `character(1)` or `NULL`.
#'   Genome build assembly name (e.g. `"GRCh38"`).
#' @param ensemblRelease `integer(1)` or `NULL`.
#'   Ensembl release version (e.g. `90`).
#'
#' @return `Tx2Gene`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "tx2gene.csv")
#' x <- readTx2Gene(
#'     file = file,
#'     organism = "Mus musculus",
#'     genomeBuild = "GRCm38",
#'     ensemblRelease = 90L
#' )
#' print(x)

## Updated 2019-07-23.
readTx2Gene <- function(
    file,
    organism = NULL,
    genomeBuild = NULL,
    ensemblRelease = NULL
) {
    data <- read_csv(
        file = file,
        col_names = c("transcriptID", "geneID"),
        col_types = "cc"  # character
    )
    data <- as(data, "DataFrame")
    metadata(data) <- list(
        organism = as.character(organism),
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease)
    )
    Tx2Gene(data)
}
