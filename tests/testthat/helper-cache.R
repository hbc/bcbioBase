dir.create("cache", showWarnings = FALSE)
invisible(lapply(
    X = c(
        "bcbio-legacy-samplename.csv",
        "bcbio-nextgen.log",
        "cellranger-metadata.csv",
        "data-versions.csv",
        "demultiplexed.csv",
        "demultiplexed-duplicated-description.csv",
        "demultiplexed-missing-cols.csv",
        "multiplexed.csv",
        "multiplexed-duplicated-sampleName.csv",
        "multiplexed-missing-cols.csv",
        "project-summary.yaml",
        "project-summary-metrics-mismatch.yaml",
        "project-summary-nested-metadata.yaml",
        "programs.txt",
        "sampleID-column-defined.csv",
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
    url = "http://tests.acidgenomics.com/bcbioBase/v0.4"
))
