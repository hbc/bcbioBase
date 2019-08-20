dir.create("cache", showWarnings = FALSE)
invisible(lapply(
    X = c(
        "bcbio-nextgen-commands.log",
        "bcbio-nextgen.log",
        "data-versions.csv",
        "programs.txt",
        "summary-invalid-metrics-mismatch.yaml",
        "summary-nested-metadata.yaml",
        "summary.yaml",
        "surecell-commands.log",
        "tx2gene.csv"
    ),
    FUN = function(file, url) {
        destfile <- file.path("cache", file)
        if (!file.exists(destfile)) {
            utils::download.file(
                url = paste(url, file, sep = "/"),
                destfile = destfile
            )
        }
    },
    url = bcbioBaseTestsURL
))
