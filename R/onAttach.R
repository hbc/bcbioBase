.onAttach <- function(libname, pkgname) {
    imports <- c(
        "basejump",
        "SummarizedExperiment"
    )

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
        X = c(imports, suggests),
        FUN = require,
        character.only = TRUE
    ))
}
