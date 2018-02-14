## bcbioBase 0.1.0 (2018-02-13)

- Preparing version pinning for [bioconda][]. Relaxed [rlang][] dependency from v0.1.6 to v0.1.2 and [tidyr][] dependency from v0.7.2 to v0.7.1.
- Added `copyToDropbox()` function, which enables input of a list of local files and returns [Dropbox][] paths using [rdrop2][].
- Updated [basejump][] dependency to v0.2.1.
- Added assertive checks for all functions.


## bcbioBase 0.0.3 (2018-01-27)

- NAMESPACE export fixes.
- Minor code cleanup to pass [lintr][] checks.
- Switched to [rlang][] methods for errors, messages, and warnings: `abort()`, `inform()`, and `warn()`.
- Messages now consistently use backticks instead of apostrophes, as recommended by the [tidyverse style guide][].
- Removed [tibble][] rownames support from `prepareSampleMetadata()`.
- Improved code coverage to 100%.


## bcbioBase 0.0.2 (2018-01-19)

- Re-export `assignAndSaveData()`, `loadData()`, `loadRemoteData()`, `saveData()` from [basejump][] package.
- Update `prepareTemplate()` function to use internally stored data.


## bcbioBase 0.0.1 (2018-01-11)

- Initial release.


[basejump]: http://steinbaugh.com/basejump
[bioconda]: https://bioconda.github.io
[Dropbox]: https://www.dropbox.com
[lintr]: https://github.com/jimhester/lintr
[rdrop2]: https://github.com/karthik/rdrop2
[rlang]: http://rlang.tidyverse.org
[tibble]: http://tibble.tidyverse.org
[tidyr]: http://tidyr.tidyverse.org
[tidyverse style guide]: http://style.tidyverse.org
