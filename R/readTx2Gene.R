#' Read Transcript-to-Gene Annotations
#'
#' Generates a `Tx2Gene` object containing `transcriptID` and `geneID` columns.
#'
#' @note Doesn't attempt to strip transcript versions.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams basejump::params
#' @param organism `string` or `NULL`. Full Latin organism name
#'   (e.g. `"Homo sapiens"`).
#' @param genomeBuild `string` or `NULL`. Genome build assembly name
#'   (e.g. `"GRCh38"`).
#' @param ensemblRelease `scalar integer` or `NULL`. Ensembl release version
#'   (e.g. `90`).
#'
#' @return `Tx2Gene`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "tx2gene.csv")
#' x <- readTx2Gene(
#'     file = file,
#'     organism = "Mus musculus",
#'     genomeBuild = "GRCm38",
#'     ensemblRelease = 90L
#' )
#' print(x)
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
