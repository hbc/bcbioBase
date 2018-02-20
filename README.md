# bcbioBase

[![Build Status](https://travis-ci.org/hbc/bcbioBase.svg?branch=master)](https://travis-ci.org/hbc/bcbioBase)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/hbc/bcbioBase/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioBase)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/r-bcbiobase/badges/version.svg)](https://anaconda.org/bioconda/r-bcbiobase)

Base functions and generics for [bcbio][] [R][] packages.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```r
source("https://bioconductor.org/biocLite.R")
biocLite(
    "hbc/bcbioBase",
    dependencies = c("Depends", "Imports", "Suggests")
)
```

### [conda][] method

```bash
conda install -c bioconda r-bcbiobase
```


[bcbio]: https://bcbio-nextgen.readthedocs.io
[Bioconductor]: https://bioconductor.org
[conda]: https://conda.io
[R]: https://www.r-project.org
