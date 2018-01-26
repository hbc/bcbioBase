.onAttach <- function(libname, pkgname) {
    # Attach imported packages
    imports <- c(
        "basejump",
        "SummarizedExperiment"
    )
    invisible(lapply(
        X = imports,
        FUN = function(package) {
            if (!package %in% (.packages())) {
                attachNamespace(package)
            }
        }
    ))

    # Attach suggested packages
    suggests <- c(
        "knitr",
        "magrittr",
        "rmarkdown",
        "tidyverse"
    )
    notInstalled <- setdiff(suggests, rownames(installed.packages()))
    if (length(notInstalled)) {
        source("https://bioconductor.org/biocLite.R")
        BiocInstaller::biocLite(pkgs = notInstalled)
    }
    invisible(lapply(
        X = suggests,
        FUN = function(package) {
            if (!package %in% (.packages())) {
                attachNamespace(package)
            }
        }
    ))
}
