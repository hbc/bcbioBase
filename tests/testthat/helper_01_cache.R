cacheURL <- "http://bcbiobase.seq.cloud"
files <- c(
    "bcbio_legacy_samplename.csv",
    "bcbio-nextgen.log",
    "data_versions.csv",
    "demultiplexed.xlsx",
    "demultiplexed_duplicated_description.csv",
    "demultiplexed_missing_cols.csv",
    "demultiplexed_with_sampleName.csv",
    "multiplexed.xlsx",
    "multiplexed_duplicated_sampleName.csv",
    "multiplexed_missing_cols.csv",
    "project-summary.yaml",
    "project-summary-metrics-mismatch.yaml",
    "programs.txt",
    "sampleID_column_defined.xlsx"
)
mapply(
    FUN = function(cacheURL, file, envir) {
        if (!file.exists(file)) {
            utils::download.file(
                url = paste(cacheURL, file, sep = "/"),
                destfile = file)
        }
        # Load R Data file
        if (grepl("\\.rda$", file)) {
            message(paste("Loading", file))
            load(file, envir = envir)
        }
    },
    file = files,
    MoreArgs = list(cacheURL = cacheURL, envir = environment())
)
