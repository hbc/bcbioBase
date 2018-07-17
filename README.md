# bcbioBase

[![Travis CI](https://travis-ci.org/hbc/bcbioBase.svg?branch=master)](https://travis-ci.org/hbc/bcbioBase)
[![AppVeyor CI](https://ci.appveyor.com/api/projects/status/j2o9aspoj8x4l9x7/branch/master?svg=true)](https://ci.appveyor.com/project/mjsteinbaugh/bcbiobase/branch/master)
[![Codecov](https://codecov.io/gh/hbc/bcbioBase/branch/master/graph/badge.svg)](https://codecov.io/gh/hbc/bcbioBase)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/r-bcbiobase/badges/version.svg)](https://anaconda.org/bioconda/r-bcbiobase)

Base functions and generics for [bcbio][] [R][] packages.


## Installation

This is an [R][] package.

### [Bioconductor][] method

We recommend using [R][] 3.5 / [Bioconductor][] 3.7.

#### R >= 3.5

```r
install.packages("BiocManager")
library(BiocManager)
install("devtools")
install("GenomeInfoDbData")
install("hbc/bcbioBase")
```

#### R < 3.5

Legacy support for [R][] 3.4 / [Bioconductor][] 3.6 is provided.

```r
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite("GenomeInfoDbData")
biocLite("hbc/bcbioBase")
```

### [conda][] method

```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -c bioconda r-bcbiobase
```


[bcbio]: https://bcbio-nextgen.readthedocs.io
[Bioconductor]: https://bioconductor.org
[conda]: https://conda.io
[R]: https://www.r-project.org
