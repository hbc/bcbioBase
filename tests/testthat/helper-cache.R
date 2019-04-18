dir.create("cache", showWarnings = FALSE)
invisible(lapply(
    X = c(
        "bcbio-nextgen.log",
        "data-versions.csv",
        "demultiplexed-invalid-duplicated.csv",
        "demultiplexed-invalid-legacy-samplename.csv",
        "demultiplexed-invalid-missing.csv",
        "demultiplexed-invalid-sample-id.csv",
        "demultiplexed.csv",
        "multiplexed-cellranger.csv",
        "multiplexed-indrops.csv",
        "multiplexed-invalid-duplicated.csv",
        "multiplexed-invalid-missing.csv",
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
