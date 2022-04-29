# bcbioBase 0.7.1 (2022-04-29)

## Minor changes

- Bumped the minimum R dependency version to 4.2.

# bcbioBase 0.7.0 (2022-03-11)

## Major changes

- Removed `copyToDropbox` function, which is not commonly used and is
  difficult to properly unit test.

## Minor changes

- Updated `importDataVersions` and `importProgramVersions` to now always
  use base R engine, as this avoids some parsing issues that can pop up when
  using readr engine instead. Primarily this applies to inconsistent date
  formatting in data versions return, which can cause a POSIX date error to
  return when using readr engine.
- `sampleDirs` function now excludes nested pipeline directories such as
  `bcbioRNASeq`, which was added in 2021.
- Now reexporting `import` (from pipette), which is used in working examples.

# bcbioBase 0.6.22 (2021-09-08)

## Minor changes

- Updated NAMESPACE to import some functions from methods rather than using
  basejump reexports: `as`, `is`, `new`.

# bcbioBase 0.6.21 (2021-03-12)

## Minor changes

- Now importing `str_match` from basejump rather than stringr, reducing the
  number of required packages in imports.

# bcbioBase 0.6.20 (2021-02-26)

## Minor changes

- Renamed "blacklist" to "denylist" internally.
- Removed reexported functions that are no longer used in other packages.

# bcbioBase 0.6.19 (2021-02-21)

## Minor changes

- Updated internal YAML parser to use improved `rbindToDataFrame` function
  (see AcidPlyr and/or basejump for detail) instead of the now deprecated
  `unlistToDataFrame` approach. The new `rbindToDataFrame` function always
  returns 1:1 from nested list elements to rows, which is ideal for sample
  metadata and quality control metrics. Alternatively, `rbindlist` from
  data.table (with `fill = TRUE`) is also work a look, but it can introduce
  unwanted row expansion and there is no way currently to enforce 1:1 mapping.
  So we wrote our own function inside of basejump.

# bcbioBase 0.6.18 (2021-02-19)

## Minor changes

- `getSampleDataFromYAML`: Fix for call to `AcidPlyr::unlistToDataFrame`,
  which now returns "name" column for nested YAML processing. See internal
  `.sampleYAML` generator for details.

# bcbioBase 0.6.17 (2021-02-10)

## Minor changes

- NAMESPACE updates, following basejump v0.14 release.
- Hardened the YAML parser against nested metadata, and updated to use
  `unlistToDataFrame` function defined in AcidPlyr, rather than deprecated
  plyr `ldply` approach.

# bcbioBase 0.6.16 (2020-12-03)

## Minor changes

- `projectDir`: Improved error message on match failure.

# bcbioBase 0.6.15 (2020-10-08)

## Minor changes

- Updated package dependency version requirements.

# bcbioBase 0.6.14 (2020-07-24)

## Minor changes

- Maintenance release, increasing minimum R dependency to 4.0.

# bcbioBase 0.6.13 (2020-02-19)

## Minor changes

- Changed license from MIT to GPL-3.

# bcbioBase 0.6.12 (2020-01-20)

## Minor changes

- Now using cli package to improve console messages.

# bcbioBase 0.6.11 (2019-10-30)

## Minor changes

- Updated basejump dependencies.
- Rebuilt documentation and checked for Bioconductor 3.10 support.

# bcbioBase 0.6.10 (2019-08-27)

## Minor changes

- Removed deprecated and defunct functions declared prior to v0.6.
- Deprecated `read*` functions in favor of `import*`. This naming format is
  more consistent with other functions used in basejump.

# bcbioBase 0.6.9 (2019-08-20)

## Major changes

- Moved `readSampleData` and `readTxGene` functions to basejump. These are
  shared methods applicable outside of bcbio. They are still reexported here
  inside the package, maintaining full backward compatibility.

## Minor changes

- Reduced the number of imports and dependencies required.
- Improved handling of characters and numerics in YAML parsing.
- Lightened the package by removing dplyr and magrittr dependencies.
- `copyToDropbox` now looks for the optional rdrop2 package, rather than
  importing as a dependency. This keeps the package lighter.

# bcbioBase 0.6.8 (2019-08-05)

## Minor changes

- Miscellaneous documentation improvements. Now using AcidRoxygen as the shared
  parameter source, following the new conventions used in basejump.
- Updated dependency versions.

# bcbioBase 0.6.7 (2019-07-23)

## Minor changes

- Updated basejump dependency versions.

# bcbioBase 0.6.6 (2019-07-17)

## Minor changes

- Documentation improvements to pass BiocCheck. Need to have working examples in
  80% of documentation files.

# bcbioBase 0.6.5 (2019-07-17)

## Minor changes

- Updated basejump dependency.
- Improved Travis CI configuration.

# bcbioBase 0.6.4 (2019-05-29)

## Major changes

- Merging this release from Acid Genomics fork back to canonical HBC repo. The
  v0.6 release series will be used as the base for new updates to bcbioRNASeq
  and bcbioSingleCell packages.

## Minor changes

- Hardened user metadata input checks, particularly in `readSampleData`. Updated
  internal `.isSampleData` assert check function to return more informative
  error messages to the user on failure.

# bcbioBase 0.6.3 (2019-05-05)

## Major changes

- Now pinned to R >= 3.5.

## Minor changes

- Improved Travis and AppVeyor CI configuration.

# bcbioBase 0.6.2 (2019-04-22)

## Minor changes

- Improved testing to ensure compatibility with R 3.4.

# bcbioBase 0.6.1 (2019-04-18)

## Minor changes

- Switched Travis CI configuration to use `rnaseq` Docker image.
- Updated [basejump][] dependencies.
- Updated testthat cache to use `tests.acidgenomics.com` URL.

# bcbioBase 0.6.0 (2019-03-28)

Updated to work with basejump v0.10 release series.

## Deprecations

- Made the deprecated bcbio ggplot2 geom functions defunct. Applies to:
  `bcbio_geom_abline`, `bcbio_geom_label`, `bcbio_geom_label_average`,
  `bcbio_geom_label_repel`.
- Additional deprecations made defunct: `readLog`, `readTx2gene`,
  `readYAMLSampleData`, `readYAMLSampleMetrics`.

## Minor changes

- Consolidated unit test data into Acid Genomics S3 bucket
  (tests.acidgenomics.com).
- Example data used for unit testing, which gets cached into `tests/testthat/`
  are now consistently formatted using kebab case instead of snake case.

# bcbioBase 0.5.14 (2019-03-22)

## Minor changes

- Updated GitHub remotes to use acidgenomics instead of steinbaugh.

# bcbioBase 0.5.13 (2019-03-18)

## Minor changes

- Migrated `metadataBlacklist` to basejump package, for use in
  `makeSummarizedExperiment` and `makeSingleCellExperiment` functions. Global
  variable is still reexported here.
- Reduced the number of reexported functions, including pipe (`%>%`) and
  `import`.

# bcbioBase 0.5.12 (2019-02-11)

## Minor changes

- `runDate()`: Updated assert check to be compatible with R 3.4 / BioC 3.6.

# bcbioBase 0.5.11 (2019-01-22)

## Minor changes

- Documentation improvements.

# bcbioBase 0.5.10 (2019-01-17)

## Minor changes

- `sampleDirs`: Added an informative message for user regarding sample name
  sanitization of cellular barcodes.

# bcbioBase 0.5.9 (2019-01-13)

## Minor changes

- Added return value for `runDate` to documentation.
- Updated Travis CI and AppVeyor CI configuration.
- Added comment regarding top-level sample metadata parsing in internal YAML
  code.

# bcbioBase 0.5.8 (2019-01-07)

## Minor changes

- Miscellaneous improvements to documentation and goalie assert checks.

# bcbioBase 0.5.7 (2018-12-21)

## Minor changes

- Updated documentation to reflect new conventions used in [basejump][]. Switch
  to using `logical(1)` instead of `string`, for example. This approach was
  inspired by the conventions used in the checkmate package.
- Switched documentation titles to use sentence case instead of title case.

# bcbioBase 0.5.6 (2018-12-12)

## Major changes

- First release that has switched to using [goalie][] internally for assert
  checks, in place of assertive package.
- Reworked internal YAML parsing functions and improved the assert checks.

## Minor changes

- Reorganized and removed some deprecated functions.
- Updated unit tests to reflect switch to goalie.

# bcbioBase 0.5.5 (2018-11-29)

## Minor changes

- Miscellaneous documentation improvements.

# bcbioBase 0.5.4 (2018-11-26)

## Minor changes

- Improve GTF file assert checks and unit testing.

# bcbioBase 0.5.3 (2018-11-25)

## Minor changes

- `getBarcodeCutoffFromCommands`: Simplified flow in `cutoff` internal variable
  assignment.
- Reorganized deprecated functions. See `deprecated.R` file.
- Split out `projectDir` and `runDate` functions into separate R files, from
  `detect.R`.
- Split out YAML parsing functions into separate files from `yaml.R`.
- Tweaked `metadataBlacklist` global and added improved comments.
- Reworked `.assertIsSampleData` and `.makeSampleData` internal code.

# bcbioBase 0.5.2 (2018-11-19)

## Minor changes

- Updated [basejump][] and [goalie][] dependencies.
- Miscellaneous documentation improvements.

# bcbioBase 0.5.1 (2018-11-15)

## Minor changes

- Switched imports to simply [basejump][], instead of attempting to reference
  any basejump subpackages.
- Migrating to [goalie][] package for internal assert checks.

# bcbioBase 0.5.0 (2018-09-14)

Working towards a release candidate for [Bioconductor][] submission.

## New functions

- `projectDir`: Returns the latest dated directory path inside a [bcbio][]
  run. If [bcbio][] has been run multiple times to the same upload directory,
  the function will return the latest project directory and warn the user.

# Major changes

- `readTx2gene` now returns a `tx2gene` class object.
- Pinning to [basejump][] [Bioconductor][] release candidate.
- Migrated [ggplot2][] geoms to [basejump][] package.
- YAML functions now return `DataFrame` instead of `data.frame`.
- `sampleDirs` is now stricter about loading sample directories with names
  that are not valid in [R][]. This now checks for non-alphanumeric characters,
  including spaces, **dashes**, and samples that begin with a number. The
  [bcbio][] pipeline will be updated to enforce these rules, to avoid unexpected
  downstream behavior in [R][] due to invalid names.
- Improved camel case sanitization of metrics columns in the YAML parser. Now
  `readYAMLSampleMetrics` should return `percentGC` instead of `xGC`. 5'->3'
  bias should return as `x5x3Bias`. The `plyr::ldply` call used to coerce
  from a `list` to `data.frame` will sanitize names, so we need to apply our
  rules before this step.

# Minor changes

- Simplified NAMESPACE imports.
- Updated deprecations.
- Linking out to recommended guidelines for development.
- Moved global `lanePattern` variable to [basejump][] package.

# bcbioBase 0.4.1 (2018-08-19)

## Minor changes

- `prepareTemplate`: Migrated function to [basejump][] package, and simplifed
  to copy all files from `rmarkdown/shared` directory inside a package.
  Currently in use by [bcbioRNASeq][] and [bcbioSingleCell][].

# bcbioBase 0.4.0 (2018-08-08)

This is a maintenance release designed to simplify the package for long-term
stability. Here we are moving all of the current S4 method support to
[basejump][] for consistency, since the methods apply to `SummarizedExperiment`
and are not [bcbio][]-workflow specific. We are keeping the package simple by
exporting the `read*` family of functions here, which are designed to integrate
with the [bcbio][] output directories.

## Major changes

- Moved S4 method support to [basejump][]: `flatFiles`, `metrics`,
  `plotCorrelationHeatmap`, `plotDot`, `plotGene`, `plotHeatmap`,
  `plotQuantileHeatmap`, `plotViolin`.
- `prepareSummarizedExperiment` has been moved to [basejump][] and renamed
  to `makeSummarizedExperiment`, for consistency. The underlying code remains
  the same.
- Moved the `SummarizedExperiment` to `list` coercion method to [basejump][].
- Reduced the number of re-exported [ggplot2][] and [viridis][] functions, for
  simplicity.

# bcbioBase 0.3.2 (2018-07-31)

## Minor changes

- Offloaded `separatorBar` and `updateMsg` to [basejump][].
- Now using [roxygen2][] v6.1 for documentation.
- Miscellaneous documentation improvements.
- Simplified installation instructions, referring the user to [bcbioRNASeq][]
  and [bcbioSingleCell][] packages instead.
- [lintr][] check fixes.

# bcbioBase 0.3.1 (2018-07-24)

## Minor changes

- Added hexadecimal color support to heatmap functions.
- Improved documentation, specifying the supported types for each argument more
  clearly.
- Improved messages in `readYAML` family of functions.

# bcbioBase 0.3.0 (2018-07-17)

## Major changes

- Moved some generics that provide `SummarizedExperiment` method suport to
  [basejump][] package: `gene2symbol`, `interestingGroups`, `sampleData`,
  `sampleNames`, `sanitizeSampleData`, `selectSamples`,
  `uniteInterestingGroups`.
- Moved method support for `SummarizedExperiment` to [basejump][]:
  `convertGenesToSymbols`, `counts`, `gene2symbol`, `interestingGroups`,
  `sampleData`, etc.
- Heatmap functions now only provide `SummarizedExperiment` method support.
  The documentation for these heatmap functions has been improved (and
  simplified), showing the supported arguments more clearly.

## Minor changes

- Now using `transcriptID` instead of `txID` in tx2gene and GRanges metadata
  defined in `rowRanges` slot of object.
- Now using `aes` instead of `aes_string` for all internal [ggplot2][] code,
  which uses tidyeval (v3.0 update).
- Moved `rse_bcb` and `rse_bcb` example data to [basejump][], for unit testing.
  These datasets are still available in [bcbioBase][].
- Made `assertFormalAnnotationCol` function defunct. No longer using to check
  `annotationCol` argument integration in heatmap functions, since we're now
  providing only `SummarizedExperiment` method support, instead of matrix
  method support (which is too complicated and error prone for users).
- Moved `assertFormalInterestingGroups` function to [basejump][], where the
  other assert check functions are defined. This helps improve package
  consistency.
- Updated unit tests to reflect changes in the number of exported generics.

# bcbioBase 0.2.16 (2018-06-28)

## Minor changes

- Improved aggregate column handling in `bcbio_geom_abline`.
- Improved code coverage back to 100%.
- Reorganized R methods files to use `-methods.R` as a suffix, as recommended
  by [Bioconductor][], in preparation for package submission.

# bcbioBase 0.2.15 (2018-06-05)

## New functions

- `minimalSampleData` enables easy creation of a sample metadata `data.frame`
  by simplify specifying the sample names (e.g. "description" in [bcbio][]
  YAML). This function was added for easy metadata handling for 10X Cell Ranger
  in the [bcbioSingleCell][] package, but is generally applicable to other
  [bcbio][] datasets and may be incorporated into [bcbioRNASeq][] in a future
  update.

## Minor changes

- Improved internal `interestingGroups` handling inside `plotHeatmap` family
  of functions.
- Removed message about transgeneNames and spikeNames in
  `prepareSummarizedExperiment` if missing transcripts are present. Now the
  function simply lists the genes that don't have metadata in `rowRanges`.
- Reorganized assert checks imports in `bcbioBase-package.R` file.

# bcbioBase 0.2.14 (2018-05-18)

## Minor changes

- Improved formatting of package NEWS.
- Improved code coverage by adding ggplot2 unit testing.
- Added `sampleData<-` method support for standard `data.frame`.

# bcbioBase 0.2.13 (2018-05-07)

## Minor changes

- Added `overwrite` parameter to `prepareTemplate`, disabled by default.
- Removed validity check for `sampleData` accessor.
- Improved vector return in `sampleNames`.

# bcbioBase 0.2.12 (2018-05-08)

## Minor changes

- Migrated `sanitizeSampleData` from basejump to bcbioBase.
- Don't return `interestingGroups` column in `sampleData` return when
  `clean = TRUE`.
- Removed internal R Markdown shared files. These are already saved per bcbio R
  package.

# bcbioBase 0.2.11 (2018-05-07)

## New methods

- Added `SummarizedExperiment` method support for `sampleNames` generic.

# bcbioBase 0.2.10 (2018-05-03)

## New functions

- Exporting [ggplot2][] convenience functions for easier plotting in
  [bcbioRNASeq][] and [bcbioSingleCell][] packages: `bcbio_geom_abline`,
  `bcbio_geom_label`, `bcbio_geom_label_average`, and `bcbio_geom_label_repel`.

# bcbioBase 0.2.9 (2018-04-30)

## Minor changes

- Improved assert check messages for interesting groups.
- Switched back to internally using `make.unique` instead of `make.names`
  for `convertGenesToSymbols` coercion on SummarizedExperiment objects. We're
  using `make.names` only when coercing seurat objects
  (see [bcbioSingleCell][] code).
- Updated `RangedSummarizedExperiment` working example, based on `bcbioRNASeq`
  example dataset.

# bcbioBase 0.2.8 (2018-04-25)

## Minor changes

- Split out assertive imports so we can pin release to [bioconda][].

# bcbioBase 0.2.7 (2018-04-25)

## Minor changes

- Updated recommended Bioconductor installation method for 3.7 release.
- Attempt to add `description` to `metadataBlacklist`.
- `readSampleData`: Use `description` metadata column internally instead of
  `sampleID`.

# bcbioBase 0.2.6 (2018-04-25)

## Minor changes

- Using `make.names` instead of `make.unique` in `convertGenesToSymbols` method
  for `SummarizedExperiment` class.

# bcbioBase 0.2.5 (2018-04-24)

## Major changes

- `sampleData` now supports `clean = TRUE` argument, which will hide columns
  that contain metrics or other calculations from the user. This is useful when
  preparing an R Markdown report, where we only want to show relevant metadata
  (e.g. factor columns). This is enabled by default.
- Now using `metadataBlacklist` to hide specific columns in sample metadata.
- `sampleMetadata` is now deprecated in favor of `sampleData`. This function
  still works but will now warn the user, and should be removed from R Markdown
  template code.
- Added `SummarizedExperiment` to `list` coercion support, which uses the same
  code as `flatFiles`.
- `readSampleData` no longer requires or recommended `fileName` column. Only
  the `description` column is required for demultiplexed data. In the case of
  multiplexed samples (e.g. inDrops single-cell RNA-seq), then `sampleName`,
  `index` and `sequence` names are also required. Multiplexed cell ranger data
  only requires `sampleName` and `index`, since the index barcode isn't present
  in the counts matrix.

## Minor changes

- `sampleID` is no longer included in metadata priority columns. This is only
  to be used internally and should be hidden from the user where possible.
- `prepareSampleData` is no longer needed and is now defunct.

# bcbioBase 0.2.4 (2018-04-22)

## Minor changes

- Reworked internal code for `readYAML` family of functions.
- Reorganized and improved default arguments for heatmap functions.
- `readSampleData`: switched to using `merge` instead of `left_join` internally.
- Reexporting viridis family of color functions, including `inferno`, and
  the British variant `scale_colour_viridis`.

# bcbioBase 0.2.3 (2018-04-19)

- `plotHeatmap`, `plotQuantileHeatmap`: Always attempt to convert genes to
  symbols for heatmaps.

# bcbioBase 0.2.2 (2018-04-16)

## Minor changes

- Improved unit tests for `convertGenesToSymbols` and `gene2symbol`.
- Miscellaneous documentation fixes.

# bcbioBase 0.2.1 (2018-04-13)

## Major changes

- `convertGenesToSymbols` now has method support for `SummarizedExperiment`.
- Updated example datasets to use `rse_bcb` and `rse_dds`.
- Simplified internal S4 method code for heatmap functions.
- `prepareSummarizedExperiment` now supports transgenes with the
  `transgeneNames` argument (e.g. "EGFP"), and FASTA spike-ins with the
  `spikeNames` argument (e.g. "ERCCs"). Additionally, `rowRanges` and `colData`
  are no longer required and can be left `NULL`, although this isn't generally
  recommended.

## Minor changes

- Consolidated `sampleYAML` function code.
- Consolidated assert check imports into `bcbioBase-package.R` file.
- Consolidated globals into `globals.R` file.
- Improved internal [pheatmap][] color support.

# bcbioBase 0.2.0 (2018-03-22)

## Major changes

- Added `sampleData` generic. Providing method support for
  `SummarizedExperiment` class. This generic provides consistent sample
  metadata support in both the [bcbioRNASeq][] and [bcbioSingleCell][]
  packages. For [bcbioSingleCell][], `colData` returns information about the
  cells, not the samples. For [bcbioRNASeq][], `colData` and `sampleData`
  provide similar information on the samples.
- Added `uniteInterestingGroups` as a generic. Provides support for adding an
  `interestingGroups` column to a data frame.
- Added `SummarizedExperiment` method support for `gene2symbol`,
  `interestingGroups` generics.
- Converted `prepareSummarizedExperiment` from a generic to a standard function.
- `prepareSummarizedExperiment` now requires rowRanges and supports FASTA
  spike-ins with `isSpike` argument. This helps set up a `SingleCellExperiment`
  object in the [bcbioSingleCell][] package.
- `prepareTemplate` has been converted from a generic to a standard function.
- `readDataVersions` and `readProgramVersions` now return an empty tibble
  if the file is missing. bcbio doesn't always output these files, so we have
  changed the behavior from stopping on a missing file to simply warning and
  returning empty.
- Added `sampleDirs` function, that informs the user about the names of the
  sample directories in the bcbio run.

## Minor changes

- Improved internal factor sanitization of sample metadata for
  `readSampleMetadataFile` and `sampleYAMLMetadata`.
- Exporting `projectDirPattern` to match bcbio project directories
  (e.g. `2018-01-01_rnaseq`) across packages.
- Changed default of `rdsToken` in `copyToDropbox` to be `NULL` instead of `NA`.
- Reduced the number of imports from basejump package.
- Now using `sessioninfo::session_info` instead of
  `devtools::session_info`.
- Added AppVeyor CI support for code testing on Windows.
- Made Travis CI checking stricter and added `BiocCheck`.
- Simplified package dependencies in DESCRIPTION.

## Deprecations

- Deprecated `sampleMetadata` in favor of `sampleData`.
- Deprecated `flatFiles` in favor of using `as(object, "list")` coercion method
  instead.
- Deprecated `prepareSampleMetadata` in favor of `prepareSampleData`.
- Made `annotable` method on `SummarizedExperiment` objects defunct.

# bcbioBase 0.1.4 (2018-02-20)

- Don't include sample metadata in `summaryYAMLMetrics` return.

# bcbioBase 0.1.3 (2018-02-19)

- Now exporting `assertFormalInterestingGroups` in camel case.

# bcbioBase 0.1.2 (2018-02-18)

- Added back `checkInterestingGroups` since this code is still present on
  [bcbioRNASeq][] master branch.

# bcbioBase 0.1.1 (2018-02-17)

- Added internal assert checks.
- Now exporting `sampleMetadata<-` assignment generic.
- Updated encrypted token for rdrop2 working example.
- Added `assert_formal_interesting_groups` assert checks. Deprecated
  `checkInterestingGroups`.

# bcbioBase 0.1.0 (2018-02-13)

- Preparing version pinning for [bioconda][]. Relaxed [rlang][] dependency from
  v0.1.6 to v0.1.2 and [tidyr][] dependency from v0.7.2 to v0.7.1.
- Added `copyToDropbox` function, which enables input of a list of local files
  and returns [Dropbox][] paths using [rdrop2][].
- Updated [basejump][] dependency to v0.2.1.
- Added assertive checks for all functions.

# bcbioBase 0.0.3 (2018-01-27)

- NAMESPACE export fixes.
- Minor code cleanup to pass [lintr][] checks.
- Switched to [rlang][] methods for errors, messages, and warnings: `abort`,
  `inform`, and `warn`.
- Messages now consistently use backticks instead of apostrophes, as recommended
  by the [tidyverse style guide][].
- Removed [tibble][] rownames support from `prepareSampleMetadata`.
- Improved code coverage to 100%.

# bcbioBase 0.0.2 (2018-01-19)

- Re-export `assignAndSaveData`, `loadData`, `loadRemoteData`, `saveData` from
  [basejump][] package.
- Update `prepareTemplate` function to use internally stored data.

# bcbioBase 0.0.1 (2018-01-11)

- Initial release.

[basejump]: http://r.acidgenomics.com/packages/basejump/
[bcbio]: https://bcbio-nextgen.readthedocs.io/
[bcbiobase]: https://r.acidgenomics.com/packages/bcbiobase/
[bcbiornaseq]: https://r.acidgenomics.com/packages/bcbiornaseq/
[bcbiosinglecell]: http://r.acidgenomics.com/packages/bcbiosinglecell/
[bioconda]: https://bioconda.github.io/
[bioconductor]: http://bioconductor.org/
[dropbox]: https://www.dropbox.com/
[ggplot2]: https://ggplot2.tidyverse.org/
[goalie]: https://goalie.acidgenomics.com/
[lintr]: https://github.com/jimhester/lintr/
[pheatmap]: https://github.com/raivokolde/pheatmap/
[rdrop2]: https://github.com/karthik/rdrop2/
[rlang]: http://rlang.tidyverse.org/
[roxygen2]: https://cran.r-project.org/package=roxygen2
[tibble]: http://tibble.tidyverse.org/
[tidyr]: http://tidyr.tidyverse.org/
[tidyverse style guide]: http://style.tidyverse.org/
[viridis]: https://cran.r-project.org/package=viridis
