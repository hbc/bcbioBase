#' Prepare Summarized Experiment
#'
#' This is a utility wrapper for [SummarizedExperiment::SummarizedExperiment()]
#' that provides automatic subsetting for row and column data.
#'
#' This function also provides automatic metadata slotting of multiple useful
#' environment parameters:
#'
#' - `date`: Today's date, returned from [Sys.Date()].
#' - `wd`: Working directory, returned from [getwd()].
#' - `utilsSessionInfo`: [utils::sessionInfo().
#' - `devtoolsSessionInfo`: [sessioninfo::session_info()].
#'
#' @importFrom basejump sanitizeColData
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom scales percent
#' @importFrom sessioninfo session_info
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom tibble has_rownames
#' @importFrom utils sessionInfo
#'
#' @inheritParams general
#' @param assays List containing RNA-seq count matrices with matching
#'   dimensions. Counts can be passed in either dense (`matrix`) or sparse
#'   (`dgCMatrix`, `dgTMatrix`) format.
#' @param rowRanges `GRanges` describing assay rows. Must contain genomic
#'   ranges. Can be left `NULL` if the genome is poorly annotated and/or ranges
#'   aren't available from AnnotationHub.
#' @param colData `DataFrame` `data.frame`, or `matrix` describing assay
#'   columns.
#' @param metadata *Optional*. Metadata `list`.
#' @param isSpike Character vector of spike-in sequence rownames.
#' @param dropRows Drop unannotated rows (e.g. deprecated genes/transcripts). Only
#'   applies when `rowRanges` is defined. Note that spike-ins declared in the
#'   `isSpike` argument will not be dropped.
#'
#' @return `RangedSummarizedExperiment`.
#' @export
#'
#' @examples
#' genes <- c(
#'     "EGFP",  # spike
#'     "gene_1",
#'     "gene_2",
#'     "gene_3",
#'     "dead_gene"
#' )
#' samples <- c(
#'     "sample_1",
#'     "sample_2",
#'     "sample_3",
#'     "sample_4"
#' )
#' mat <- matrix(
#'     seq(1L:20L),
#'     nrow = 5L,
#'     ncol = 4L,
#'     dimnames = list(genes, samples)
#' )
#' # Leave out the unannotated EGFP spike-in
#' rowRanges <- GRanges(
#'     seqnames = c("1", "1", "1"),
#'     ranges = IRanges(
#'         start = c(1L, 101L, 201L),
#'         end = c(100L, 200L, 300L)
#'     )
#' )
#' names(rowRanges) <- c("gene_1", "gene_2", "gene_3")
#' colData <- data.frame(
#'     "genotype" = c(
#'         "wildtype",
#'         "wildtype",
#'         "knockout",
#'         "knockout"
#'     ),
#'     "age" = c(3L, 6L, 3L, 6L),
#'     row.names = samples
#' )
#' prepareSummarizedExperiment(
#'     assays = list(assay = mat),
#'     rowRanges = rowRanges,
#'     colData = colData,
#'     isSpike = "EGFP"
#' )
prepareSummarizedExperiment <- function(
    assays,
    rowRanges = NULL,
    colData = NULL,
    metadata = NULL,
    isSpike = NULL
) {

    # Assays ===================================================================
    assert_is_list(assays)
    # Drop any `NULL` items from list
    assays <- Filter(Negate(is.null), assays)
    assay <- assays[[1L]]
    assert_has_dimnames(assay)
    assert_has_rownames(assay)
    assert_has_colnames(assay)
    assert_has_no_duplicates(rownames(assay))
    assert_has_no_duplicates(colnames(assay))
    assert_are_identical(
        x = make.names(rownames(assay), unique = TRUE, allow_ = TRUE),
        y = rownames(assay)
    )
    assert_are_identical(
        x = make.names(colnames(assay), unique = TRUE, allow_ = TRUE),
        y = colnames(assay)
    )
    # Ensure that all slotted items have the same dimensions and names
    invisible(lapply(
        X = assays,
        FUN = function(x) {
            assert_are_identical(dim(x), dim(assay))
            assert_are_identical(dimnames(x), dimnames(assay))
        }
    ))

    # Row ranges ===============================================================
    assert_is_any_of(rowRanges, c("GRanges", "NULL"))
    unannotatedRows <- character()
    if (is(rowRanges, "GRanges")) {
        assert_are_intersecting_sets(rownames(assay), names(rowRanges))
        unannotatedRows <- setdiff(rownames(assay), names(rowRanges))
        # Create placeholder ranges for spike-ins
        if (length(unannotatedRows) && is.character(isSpike)) {
            assert_is_subset(isSpike, unannotatedRows)
            unannotatedRows <- unannotatedRows %>%
                .[-match(isSpike, .)]
            vec <- paste("spike", "1-100", sep = ":")
            vec <- replicate(n = length(isSpike), expr = vec)
            spikes <- GRanges(vec)
            names(spikes) <- isSpike
            # Create the required empty metadata columns
            mcols(spikes) <- matrix(
                nrow = length(spikes),
                ncol = ncol(mcols(rowRanges)),
                dimnames = list(
                    isSpike,
                    colnames(mcols(rowRanges))
                )
            ) %>%
                as("DataFrame")
            # Warning about no sequence levels in common is expected here
            rowRanges <- suppressWarnings(c(spikes, rowRanges))
        }

        # Warn the user about dropping unannotated rows
        if (length(unannotatedRows)) {
            warn(paste(
                "Dropping", length(unannotatedRows), "unannotated rows",
                paste0(
                    "(",
                    percent(length(unannotatedRows) / nrow(assay)),
                    "):"
                ),
                toString(unannotatedRows)
            ))
            intersect <- intersect(rownames(assay), names(rowRanges))
            assays <- mapply(
                assay = assays,
                MoreArgs = list(intersect = intersect),
                FUN = function(assay, intersect) {
                    assay[intersect, , drop = FALSE]
                },
                USE.NAMES = TRUE,
                SIMPLIFY = FALSE
            )
            rowRanges <- rowRanges[intersect]
        }
    }

    # Column data ==============================================================
    assert_is_any_of(
        x = colData,
        classes = c("DataFrame", "data.frame", "matrix", "NULL")
    )
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
    assert_is_any_of(metadata, c("list", "NULL"))
    metadata <- as.list(metadata)
    metadata[["date"]] <- Sys.Date()
    metadata[["wd"]] <- normalizePath(".")
    metadata[["utilsSessionInfo"]] <- sessionInfo()
    metadata[["devtoolsSessionInfo"]] <- session_info(include_base = TRUE)
    metadata[["unannotatedRows"]] <- unannotatedRows

    # Return ===================================================================
    SummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata
    )
}
