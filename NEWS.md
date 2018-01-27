## bcbioBase 0.0.3 (2018-01-27)

- NAMESPACE export fixes.
- Minor code cleanup to pass lintr checks.
- Switched to rlang methods for errors, messages, and warnings: `abort()`, `inform()`, and `warn()`.
- Messages now consistently use backticks instead of apostrophes, as recommended by the tidyverse style guide.
- Removed tibble rownames support from `prepareSampleMetadata()`.
- Improved code coverage to 100%.


## bcbioBase 0.0.2 (2018-01-19)

- Re-export `assignAndSaveData()`, `loadData()`, `loadRemoteData()`, `saveData()` from basejump package.
- Update `prepareTemplate()` function to use internally stored data.


## bcbioBase 0.0.1 (2018-01-11)

- Initial release.
