dir.create("cache", showWarnings = FALSE)
invisible(lapply(
    X = c(
        "bcbio-nextgen-commands.log",
        "bcbio-nextgen.log",
        "data-versions.csv",
        "metadata-demultiplexed-invalid-duplicated.csv",
        "metadata-demultiplexed-invalid-legacy-samplename.csv",
        "metadata-demultiplexed-invalid-missing-columns.csv",
        "metadata-demultiplexed-invalid-sample-id.csv",
        "metadata-demultiplexed.csv",
        "metadata-invalid-column-name.csv",
        "metadata-invalid-description.csv",
        "metadata-multiplexed-cellranger.csv",
        "metadata-multiplexed-indrops.csv",
        "metadata-multiplexed-invalid-duplicated.csv",
        "metadata-multiplexed-invalid-missing-columns.csv",
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
