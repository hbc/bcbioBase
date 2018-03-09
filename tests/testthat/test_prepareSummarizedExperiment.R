context("prepareSummarizedExperiment")

genes <- c(
    "EGFP",  # spike
    "gene_1",
    "gene_2",
    "gene_3",
    "dead_gene"
)
samples <- c(
    "sample_1",
    "sample_2",
    "sample_3",
    "sample_4"
)
mat <- matrix(
    seq(1L:20L),
    nrow = 5L,
    ncol = 4L,
    dimnames = list(genes, samples)
)
# Leave out the unannotated EGFP spike-in
rowRanges <- GRanges(
    seqnames = c("1", "1", "1"),
    ranges = IRanges(
        start = c(1L, 101L, 201L),
        end = c(100L, 200L, 300L)
    )
)
names(rowRanges) <- c("gene_1", "gene_2", "gene_3")
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
    rse <- suppressWarnings(prepareSummarizedExperiment(
        assays = list(assay = mat),
        rowRanges = rowRanges,
        colData = colData,
        isSpike = "EGFP"
    ))
    expect_s4_class(rse, "RangedSummarizedExperiment")
    expect_identical(dim(rse), c(4L, 4L))
    expect_identical(
        names(rse),
        c("EGFP", "gene_1", "gene_2", "gene_3")
    )
    expect_identical(
        metadata(rse)[["isSpike"]],
        "EGFP"
    )
    expect_identical(
        metadata(rse)[["unannotatedRows"]],
        "dead_gene"
    )
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

# This checks to see if there are any dashes in the names
test_that("Enforce strict names", {
    matBadRows <- mat
    rownames(matBadRows) <- paste0(rownames(matBadRows), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matBadRows),
            rowData = rowData,
            colData = colData
        ),
        paste(
            "are_identical :",
            "make.names\\(rownames\\(assay\\), unique = TRUE, allow_ = TRUE\\)",
            "and rownames\\(assay\\) are not identical."
        )
    )
    matBadCols <- mat
    colnames(matBadCols) <- paste0(colnames(matBadCols), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matBadCols),
            rowData = rowData,
            colData = colData
        ),
        paste(
            "are_identical :",
            "make.names\\(colnames\\(assay\\), unique = TRUE, allow_ = TRUE\\)",
            "and colnames\\(assay\\) are not identical."
        )
    )
})

test_that("Duplicate names", {
    matDupeRows <- mat
    rownames(matDupeRows) <- c(
        "ENSMUSG00000000001",
        "ENSMUSG00000000001",
        "ENSMUSG00000000003",
        "ENSMUSG00000000003"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matDupeRows),
            rowData = rowData,
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
            rowData = rowData,
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
            rowData = rowData,
            colData = colData
        ),
        "has_dimnames : The dimension names of assay are NULL."
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowData,
            colData = c(xxx = "yyy")
        ),
        "is2 : colData"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = c(xxx = "yyy"),
            colData = colData
        ),
        "is2 : rowData"
    )
})

test_that("Invalid metadata", {
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowData,
            colData = colData,
            metadata = Sys.Date()
        ),
        "is2 : metadata")
})

test_that("Dimension mismatch", {
    matExtraCol <- cbind(mat, "sample_5" = seq(17L, 20L))
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matExtraCol),
            rowData = rowData,
            colData = colData
        ),
        paste(
            "are_identical :",
            "colnames\\(assay\\) and rownames\\(colData\\) are not identical."
        )
    )
    matExtraRow <- rbind(mat, "ENSMUSG00000000037" = seq(17L, 20L))
    expect_warning(
        prepareSummarizedExperiment(
            assays = list(assay = matExtraRow),
            rowData = rowData,
            colData = colData
        ),
        "1 unannotated rows detected")
    seExtraRow <- suppressWarnings(
        prepareSummarizedExperiment(
            assays = list(assay = matExtraRow),
            rowData = rowData,
            colData = colData
        )
    )
    expect_identical(
        metadata(seExtraRow)[["unannotatedRows"]],
        "ENSMUSG00000000037"
    )
})

test_that("Missing rownames", {
    # Missing rownames
    matnorownames <- mat
    rownames(matnorownames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matnorownames),
            rowData = rowData,
            colData = colData
        ),
        "has_rownames : The row names of assay are NULL."
    )
    matnocolnames <- mat
    colnames(matnocolnames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = matnocolnames),
            rowData = rowData,
            colData = colData
        ),
        "has_colnames : The column names of assay are NULL."
    )
    rowdatanorownames <- rowData
    rownames(rowdatanorownames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowdatanorownames,
            colData = colData
        ),
        paste(
            "are_intersecting_sets :",
            "rownames\\(assay\\) and rownames\\(rowData\\)",
            "have no common elements."
        )
    )
    coldatanorownames <- colData
    rownames(coldatanorownames) <- NULL
    expect_error(
        prepareSummarizedExperiment(
            assays = list(assay = mat),
            rowData = rowData,
            colData = coldatanorownames
        ),
        paste(
            "are_identical :",
            "colnames\\(assay\\) and rownames\\(colData\\) are not identical."
        )
    )
})
