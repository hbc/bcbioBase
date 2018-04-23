context("Prepare Functions")



# prepareSampleData ============================================================
# FIXME This needs more coverage

test_that("prepareSampleData : Missing description column", {
    expect_error(
        prepareSampleData(mtcars),
        paste(
            "is_subset :",
            "The element 'description' in \"description\" is not in",
            "colnames\\(object\\)."
        )
    )
})



# prepareSummarizedExperiment ==================================================
test_that("prepareSummarizedExperiment : RangedSummarizedExperiment", {
    rse <- prepareSummarizedExperiment(
        assays = list("counts" = mat),
        rowRanges = rr,
        colData = cd
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
            "devtoolsSessionInfo" = "session_info"
        )
    )
})

test_that("prepareSummarizedExperiment : Super minimal", {
    rse <- suppressWarnings(prepareSummarizedExperiment(
        assays = list("counts" = mat),
        rowRanges = NULL,
        colData = NULL
    ))
    expect_s4_class(rse, "RangedSummarizedExperiment")
    expect_identical(levels(seqnames(rse)), "unknown")
})

test_that("prepareSummarizedExperiment : Spike-in support", {
    rownames(mat)[1L:2L] <- c("EGFP", "ERCC")
    rse <- prepareSummarizedExperiment(
        assays = list("counts" = mat),
        rowRanges = rr[3L:4L],
        colData = cd,
        transgeneNames = "EGFP",
        spikeNames = "ERCC"
    )
    expect_identical(
        rownames(rse),
        c("EGFP", "ERCC", genes[3L:4L])
    )
    expect_identical(
        levels(seqnames(rse)),
        c("spike", "transgene", "1")
    )
})

test_that("prepareSummarizedExperiment : Strict names", {
    # Don't allow any dashes and other illegal characters in names
    matBadRows <- mat
    rownames(matBadRows) <- paste0(rownames(matBadRows), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list("counts" = matBadRows),
            rowRanges = rr,
            colData = cd
        ),
        "are_identical : makeNames\\(rownames\\(assay\\)"
    )
    matBadCols <- mat
    colnames(matBadCols) <- paste0(colnames(matBadCols), "-XXX")
    expect_error(
        prepareSummarizedExperiment(
            assays = list("counts" = matBadCols),
            rowRanges = rr,
            colData = cd
        ),
        "are_identical : makeNames\\(colnames\\(assay\\)"
    )
})

test_that("prepareSummarizedExperiment : Duplicate names", {
    matDupeRows <- mat
    rownames(matDupeRows) <- c(
        "gene_1",
        "gene_1",
        "gene_2",
        "gene_2"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list("counts" = matDupeRows),
            rowRanges = rr,
            colData = cd
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
            assays = list("counts" = matDupeCols),
            rowRanges = rr,
            colData = cd
        ),
        paste(
            "has_no_duplicates :",
            "colnames\\(assay\\) has duplicates at positions 2, 4."
        )
    )
})

test_that("prepareSummarizedExperiment : Column data failure", {
    # Bad pass-in of objects not supporting `dimnames()`
    expect_error(
        prepareSummarizedExperiment(
            assays = list("counts" = "yyy"),
            rowRanges = rr,
            colData = cd
        ),
        "has_dimnames : The dimension names of assay are NULL."
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list("counts" = mat),
            rowRanges = rr,
            colData = c(xxx = "yyy")
        ),
        "is2 : colData"
    )
    expect_error(
        prepareSummarizedExperiment(
            assays = list("counts" = mat),
            rowRanges = c(xxx = "yyy"),
            colData = cd
        ),
        "is2 : rowRanges"
    )
})

test_that("prepareSummarizedExperiment : Invalid metadata", {
    expect_error(
        prepareSummarizedExperiment(
            assays = list("counts" = mat),
            rowRanges = rr,
            colData = cd,
            metadata = Sys.Date()
        ),
        "is2 : metadata"
    )
})



# prepareTemplate ==============================================================
test_that("prepareTemplate : All default shared files", {
    files <- c(
        "_footer.Rmd",
        "_header.Rmd",
        "_output.yaml",
        "_setup.R",
        "bibliography.bib"
    )
    expect_silent(prepareTemplate())
    expect_true(all(file.exists(files)))
    unlink(files)
})

test_that("prepareTemplate : Single file", {
    prepareTemplate("bibliography.bib")
    expect_true(file.exists("bibliography.bib"))
    unlink("bibliography.bib")
})

test_that("prepareTemplate : Missing source file", {
    expect_error(
        prepareTemplate("XXX.R"),
        "is_existing_file :"
    )
})
