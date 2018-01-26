.onAttach <- function(libname, pkgname) {
    packages <- c(
        "basejump",
        "SummarizedExperiment"
    )
    lapply(packages, function(package) {
        if (!package %in% (.packages())) {
            attachNamespace(package)
        }
    })
    invisible()
}
