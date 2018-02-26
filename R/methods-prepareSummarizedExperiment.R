#' Prepare Summarized Experiment
#'
#' This is a utility wrapper for `SummarizedExperiment()` that provides
#' automatic subsetting for `colData` and `rowData`.
#'
#' This function also provides automatic metadata slotting of multiple useful
#' environment parameters:
#'
#' - `date`: Today's date.
#' - `wd`: Working directory.
#' - `utilsSessionInfo`: [utils::sessionInfo()] return.
#' - `devtoolsSessionInfo`: [devtools::session_info()] return.
#'
#' @rdname prepareSummarizedExperiment
#' @name prepareSummarizedExperiment
#' @family bcbio Utilities
#'
#' @inheritParams general
#'
#' @param assays List containing RNA-seq count matrices with matching
#'   dimensions. Counts can be passed in either dense (`matrix`) or sparse
#'   (`dgCMatrix`, `dgTMatrix`) format.
#' @param rowData Object describing assay matrix rows. Must support
#'   [base::dim()].
#' @param colData Object describing assay matrix columns. Must support
#'   [base::dim()].
#' @param metadata *Optional*. Metadata list.
#'
#' @seealso
#' - [SummarizedExperiment::SummarizedExperiment()].
#' - [base::Sys.Date()].
#' - [base::getwd()].
#' - [utils::sessionInfo()].
#'
#' @return [SummarizedExperiment].
#' @export
#'
#' @examples
#' mat <- matrix(
#'     seq(1L:16L),
#'     nrow = 4L,
#'     ncol = 4L,
#'     dimnames = list(
#'         c(
#'             "ENSMUSG00000000001",
#'             "ENSMUSG00000000003",
#'             "ENSMUSG00000000028",
#'             "ENSMUSG00000000031"
#'         ),
#'         c(
#'             "sample_1",
#'             "sample_2",
#'             "sample_3",
#'             "sample_4"
#'         )
#'     )
#' )
#' rowData <- data.frame(
#'     ensgene = c(
#'         "ENSMUSG00000000001",
#'         "ENSMUSG00000000003",
#'         "ENSMUSG00000000028",
#'         "ENSMUSG00000000031"
#'     ),
#'     biotype = c(
#'         "coding",
#'         "coding",
#'         "coding",
#'         "coding"
#'     ),
#'     row.names = rownames(mat)
#' )
#' colData <- data.frame(
#'     genotype = c(
#'         "wildtype",
#'         "wildtype",
#'         "knockout",
#'         "knockout"
#'     ),
#'     age = c(3, 6, 3, 6),
#'     row.names = colnames(mat)
#' )
#' prepareSummarizedExperiment(
#'     assays = list(assay = mat),
#'     rowData = rowData,
#'     colData = colData
#' )
NULL



# Constructors =================================================================
#' @importFrom fs path_real
#' @importFrom scales percent
#' @importFrom sessioninfo session_info
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom tibble has_rownames
#' @importFrom utils head
.prepareSummarizedExperiment <- function(
    assays,
    rowData = NULL,
    colData = NULL,
    metadata = NULL) {
    assert_is_list(assays)
    validData <- c("data.frame", "DataFrame", "matrix", "NULL")
    assert_is_any_of(rowData, validData)
    assert_is_any_of(colData, validData)
    assert_is_any_of(metadata, c("list", "SimpleList", "NULL"))

    # Assays ===================================================================
    # Drop any `NULL` items from list.
    assays <- Filter(Negate(is.null), assays)

    assay <- assays[[1L]]
    assert_has_dimnames(assay)
    assert_has_rownames(assay)
    assert_has_colnames(assay)
    assert_has_no_duplicates(rownames(assay))
    assert_has_no_duplicates(colnames(assay))
    assert_are_identical(
        make.names(rownames(assay), unique = TRUE, allow_ = TRUE),
        rownames(assay)
    )
    assert_are_identical(
        make.names(colnames(assay), unique = TRUE, allow_ = TRUE),
        colnames(assay)
    )

    # Ensure that all slotted items have the same dimensions and names
    lapply(
        X = assays,
        FUN = function(x) {
            assert_are_identical(dim(x), dim(assay))
            assert_are_identical(dimnames(x), dimnames(assay))
    })

    # Row data =================================================================
    if (!is.null(rowData)) {
        rowData <- as.data.frame(rowData)
        rowData <- as(rowData, "DataFrame")
        assert_are_intersecting_sets(rownames(assay), rownames(rowData))

        # Check for unannotated genes not found in annotable. This typically
        # includes gene identifiers that are now deprecated on Ensembl and/or
        # FASTA spike-in identifiers. Warn on detection rather than stopping.
        if (!all(rownames(assay) %in% rownames(rowData))) {
            unannotatedGenes <- setdiff(rownames(assay), rownames(rowData)) %>%
                sort()
            warn(paste(
                "Unannotated genes detected in assay",
                paste0(
                    "(", percent(length(unannotatedGenes) / nrow(assay)), ")"
                )
            ))
        } else {
            unannotatedGenes <- NULL
        }

        # Allow for dynamic resizing of rows to match gene annotable input
        rowData <- rowData[rownames(assay), , drop = FALSE]
        rownames(rowData) <- rownames(assay)
    } else {
        warn(paste(
            "Summarizing experiment without row data",
            "(e.g. gene annotations)"
        ))
        rowData <- DataFrame(row.names = rownames(assay))
        unannotatedGenes <- NULL
    }

    # Column data ==============================================================
    if (!is.null(colData)) {
        colData <- as.data.frame(colData)
        colData <- as(colData, "DataFrame")
        assert_are_identical(colnames(assay), rownames(colData))
    } else {
        warn(paste(
            "Summarizing experiment without column data",
            "(e.g. sample metadata)"
        ))
        colData <- DataFrame(row.names = colnames(assay))
    }

    # Metadata =================================================================
    if (is.null(metadata)) {
        metadata <- list()
    }
    metadata[["date"]] <- Sys.Date()
    metadata[["wd"]] <- path_real(".")
    metadata[["sessionInfo"]] <- session_info()
    metadata[["unannotatedGenes"]] <- unannotatedGenes

    # Return ===================================================================
    SummarizedExperiment(
        assays = assays,
        rowData = rowData,
        colData = colData,
        metadata = metadata)
}



# Methods ======================================================================
#' @rdname prepareSummarizedExperiment
#' @export
setMethod(
    "prepareSummarizedExperiment",
    signature("list"),
    .prepareSummarizedExperiment)
