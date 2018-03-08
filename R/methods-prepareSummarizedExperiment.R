# TODO Add `isSpike` param
# By default drop genes that don't have an annotation. Maybe add this as
# `drop` argument, allowing user to keep unannotated genes that aren't spike-ins
# When slotting spike-ins, use `plasmid` in the GRanges seqnames



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
#' @name prepareSummarizedExperiment
#'
#' @inheritParams general
#'
#' @param assays List containing RNA-seq count matrices with matching
#'   dimensions. Counts can be passed in either dense (`matrix`) or sparse
#'   (`dgCMatrix`, `dgTMatrix`) format.
#' @param rowData Object describing assay rows. Must support [base::dim()].
#' @param colData Object describing assay columns. Must support [base::dim()].
#' @param metadata *Optional*. Metadata `list`.
#'
#' @seealso
#' - [SummarizedExperiment::SummarizedExperiment()].
#' - [base::Sys.Date()].
#' - [base::getwd()].
#' - [utils::sessionInfo()].
#'
#' @return `SummarizedExperiment`.
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
#' @importFrom basejump sanitizeColData
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom scales percent
#' @importFrom sessioninfo session_info
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom tibble has_rownames
#' @importFrom utils head sessionInfo
.prepareSummarizedExperiment <- function(
    assays,
    rowData = NULL,
    colData = NULL,
    metadata = NULL
) {
    assert_is_list(assays)
    assert_is_any_of(
        rowData,
        classes = c(
            "GRanges",
            "data.frame",
            "DataFrame",
            "matrix",
            "NULL"
        )
    )
    assert_is_any_of(
        x = colData,
        classes = c(
            "data.frame",
            "DataFrame",
            "matrix",
            "NULL"
        )
    )
    assert_is_any_of(
        x = metadata,
        classes = c("list", "SimpleList", "NULL")
    )

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
    invisible(lapply(
        X = assays,
        FUN = function(x) {
            assert_are_identical(dim(x), dim(assay))
            assert_are_identical(dimnames(x), dimnames(assay))
        }
    ))

    # Row data =================================================================
    # Check for unannotated genes not found in annotable. This typically
    # includes gene identifiers that are now deprecated on Ensembl and/or
    # FASTA spike-in identifiers. Warn on detection rather than stopping.
    unannotatedRows <- NULL
    if (is(rowData, "GRanges")) {
        assert_are_intersecting_sets(rownames(assay), names(rowData))
        unannotatedRows <- sort(setdiff(rownames(assay), names(rowData)))
        if (length(unannotatedRows)) {
            # TODO Improved method for stashing empty ranges?
            # Use `plasmid` here for spike-ins
            # Stash the missing rows at the first seqname (e.g. chromosome 1)
            # with the ranges 1-2, and no strand.
            vec <- paste(
                levels(seqnames(rowData))[[1L]],
                "1-2",
                sep = ":"
            )
            empty <- GRanges(
                replicate(
                    n = length(unannotatedRows),
                    expr = vec)
            )
            names(empty) <- unannotatedRows
            # Create the required empty metadata columns
            mcols(empty) <- matrix(
                nrow = length(unannotatedRows),
                ncol = ncol(mcols(rowData)),
                dimnames = list(
                    unannotatedRows,
                    colnames(mcols(rowData))
                )
            ) %>%
                as("DataFrame")
            rowData <- c(rowData, empty)
        }
        rowData <- rowData[rownames(assay)]
    } else if (!is.null(rowData)) {
        # Coerce to DataFrame, if necessary
        if (!is(rowData, "DataFrame")) {
            rowData <- rowData %>%
                as.data.frame() %>%
                as("DataFrame")
        }
        assert_are_intersecting_sets(rownames(assay), rownames(rowData))
        unannotatedRows <- sort(setdiff(rownames(assay), rownames(rowData)))
        # Allow for dynamic resizing of rows to match gene annotable input
        rowData <- rowData[rownames(assay), , drop = FALSE]
        rownames(rowData) <- rownames(assay)
    } else {
        warn("Summarizing experiment without row data")
        rowData <- DataFrame(row.names = rownames(assay))
    }

    # Warn the user on detection of unannotated rows
    if (length(unannotatedRows)) {
        warn(paste(
            length(unannotatedRows),
            "unannotated rows detected",
            paste0("(", percent(length(unannotatedRows) / nrow(assay)), "):"),
            toString(unannotatedRows)
        ))
    } else {
        unannotatedRows <- NULL
    }

    # Column data ==============================================================
    if (!is.null(colData)) {
        # Coerce to DataFrame, if necessary
        if (!is(colData, "DataFrame")) {
            colData <- colData %>%
                as.data.frame() %>%
                as("DataFrame")
        }
        assert_are_identical(colnames(assay), rownames(colData))
        colData <- sanitizeColData(colData)
    } else {
        warn("Summarizing experiment without column data")
        colData <- DataFrame(row.names = colnames(assay))
    }

    # Metadata =================================================================
    if (is.null(metadata)) {
        metadata <- list()
    }
    metadata[["date"]] <- Sys.Date()
    metadata[["wd"]] <- normalizePath(".")
    metadata[["utilsSessionInfo"]] <- sessionInfo()
    metadata[["devtoolsSessionInfo"]] <- session_info(include_base = TRUE)
    metadata[["unannotatedRows"]] <- unannotatedRows

    # Return ===================================================================
    SummarizedExperiment(
        assays = assays,
        rowData = rowData,
        colData = colData,
        metadata = metadata
    )
}



# Methods ======================================================================
#' @rdname prepareSummarizedExperiment
#' @export
setMethod(
    "prepareSummarizedExperiment",
    signature("list"),
    .prepareSummarizedExperiment
)
