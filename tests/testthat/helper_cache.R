cacheURL <- "http://bcbiobase.seq.cloud"
files <- c(
    "bcbio_legacy_samplename.csv",
    "bcbio-nextgen.log",
    "cellranger_metadata.csv",
    "data_versions.csv",
    "demultiplexed.csv",
    "demultiplexed_duplicated_description.csv",
    "demultiplexed_missing_cols.csv",
    "multiplexed.csv",
    "multiplexed_duplicated_sampleName.csv",
    "multiplexed_missing_cols.csv",
    "project-summary.yaml",
    "project-summary-metrics-mismatch.yaml",
    "project-summary-nested-metadata.yaml",
    "programs.txt",
    "sampleID_column_defined.csv",
    "tx2gene.csv"
)
mapply(
    FUN = function(cacheURL, file, envir) {
        if (!file.exists(file)) {
            utils::download.file(
                url = paste(cacheURL, file, sep = "/"),
                destfile = file)
        }
    },
    file = files,
    MoreArgs = list(cacheURL = cacheURL, envir = environment())
)
