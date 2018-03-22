context("prepareSummarizedExperiment")

genes <- c(
    "gene_1",
    "gene_2",
    "gene_3",
    "gene_4"
)
samples <- c(
    "sample_1",
    "sample_2",
    "sample_3",
    "sample_4"
)
mat <- matrix(
    seq(1L:16L),
    nrow = 4L,
    ncol = 4L,
    dimnames = list(genes, samples)
)
rowRanges <- GRanges(
    seqnames = replicate(n = 4L, expr = "1"),
    ranges = IRanges(
        start = c(1L, 101L, 201L, 301L),
        end = c(100L, 200L, 300L, 400L)
    )
)
names(rowRanges) <- genes
colData <- data.frame(
    "genotype" = c(
        "wildtype",
        "wildtype",
        "knockout",
        "knockout"
    ),
    "age" = c(3L, 6L, 3L, 6L),
    row.names = samples
)

test_that("RangedSummarizedExperiment", {
    rse <- prepareSummarizedExperiment(
        assays = list(assay = mat),
        rowRanges = rowRanges,
        colData = colData
    )
    expect_s4_class(rse, "RangedSummarizedExperiment")
    expect_identical(dim(rse), c(4L, 4L))
    expect_identical(names(rse), genes)
    expect_identical(
        lapply(metadata(rse), class),
        list(
            "date" = "Date",
            "wd" = "character",
            "utilsSessionInfo" = "sessionInfo",
            "devtoolsSessionInfo" = "session_info",
            "isSpike" = "character",
            "unannotatedRows" = "character"
        )
    )
})

test_that("Spike-in support", {
    rownames(mat)[[1L]] <- "EGFP"
    rowRanges <- rowRanges[2L:4L]
    rse <- prepareSummarizedExperiment(
        assays = list(mat),
        rowRanges = rowRanges,
        colData = colData,
        isSpike = "EGFP"
    )
    expect_identical(
        rownames(rse),
        c("EGFP", "gene_2", "gene_3", "gene_4")
    )
    expect_identical(
        metadata(rse)[["isSpike"]],
        "EGFP"
    )
})

test_that("Unannotated rows", {
    rowRanges <- rowRanges[seq_len(3L)]
    rse <- suppressWarnings(
        prepareSummarizedExperiment(
            assays = list(mat),
            rowRanges = rowRanges,
            colData = colData
        )
    )
    expect_identical(
        metadata(rse)[["unannotatedRows"]],
        "gene_4"
    )
})

test_that("Strict names", {
    # Don't allow any dashes and other illegal characters in names
    matBadRows <- mat
    rownames(matBadRows) <- paste0(rownames(matBadRows), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matBadRows),
            rowRanges = rowRanges,
            colData = colData
        ),
        "are_identical : makeNames\\(rownames\\(assay\\)"
    )
    matBadCols <- mat
    colnames(matBadCols) <- paste0(colnames(matBadCols), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matBadCols),
            rowRanges = rowRanges,
            colData = colData
        ),
        "are_identical : makeNames\\(colnames\\(assay\\)"
    )
})

test_that("Duplicate names", {
    matDupeRows <- mat
    rownames(matDupeRows) <- c(
        "gene_1",
        "gene_1",
        "gene_2",
        "gene_2"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matDupeRows),
            rowRanges = rowRanges,
            colData = colData
        ),
        paste(
            "has_no_duplicates :",
            "rownames\\(assay\\) has duplicates at positions 2, 4."
        )
    )
    matDupeCols <- mat
    colnames(matDupeCols) <- c(
        "sample_1",
        "sample_1",
        "sample_2",
        "sample_2"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matDupeCols),
            rowRanges = rowRanges,
            colData = colData
        ),
        paste(
            "has_no_duplicates :",
            "colnames\\(assay\\) has duplicates at positions 2, 4."
        )
    )
})

test_that("Column data pass-in failure in assays", {
    # Bad pass-in of objects not supporting `dimnames()`
    expect_error(
        prepareSummarizedExperiment(
            assays = list(c(xxx = "yyy")),
            rowRanges = rowRanges,
            colData = colData
        ),
        "has_dimnames : The dimension names of assay are NULL."
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowRanges = rowRanges,
            colData = c(xxx = "yyy")
        ),
        "is2 : colData"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowRanges = c(xxx = "yyy"),
            colData = colData
        ),
        "is2 : rowRanges"
    )
})

test_that("Invalid metadata", {
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowRanges = rowRanges,
            colData = colData,
            metadata = Sys.Date()
        ),
        "is2 : metadata"
    )
})
