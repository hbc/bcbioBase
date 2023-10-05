lst <- AcidDevTools::cacheTestFiles(
    pkg = .pkgName,
    files = c(
        "bcbio-nextgen-commands.log",
        "bcbio-nextgen.log",
        "data-versions.csv",
        "programs.txt",
        "summary-invalid-metrics-mismatch.yaml",
        "summary-nested-metadata.yaml",
        "summary.yaml",
        "surecell-commands.log",
        "tx2gene.csv"
    )
)
cacheDir <- lst[["cacheDir"]]
rm(lst)
