# Changelog

## eatModel 0.10.19

Changes based on available commit messages. The local Git history can be
traced back to the initial commit on 2018-10-30.

### 2025-08-25 to 2026-06-10

#### New features

- Added partial credit support for TAM, mirt, and ConQuest.
- Added DIF handling for partial credit models in TAM and ConQuest.
- Added partial credit support to
  [`plotICC()`](https://weirichs.github.io/eatModel/reference/plotICC.md).
- Added
  [`getTrafo()`](https://weirichs.github.io/eatModel/reference/getTrafo.md)
  to extract IQB transformation rules from a database.
- Added an equivalence table for partial credit models.

#### Improvements

- Improved WLE handling in multidimensional models by accounting for
  missing values on subdimensions.
- Updated `adaptSkeletonForAnchor()` so covariances are not constrained.
- Omit empty categories.
- Allow ConQuest variable names longer than 10 characters.
- Restrict `group.var` to `software = "conquest"`.

#### Bug fixes

- Fixed `getDevianceConquest()`.
- Fixed `getConquestDesc()`.
- Fixed
  [`runModel()`](https://weirichs.github.io/eatModel/reference/runModel.md)
  for partial credit models in mirt.
- Fixed `checkGenerateRefPop()`.
- Fixed issues in Thurstone threshold tests.

#### Documentation and tests

- Added documentation for
  [`getTrafo()`](https://weirichs.github.io/eatModel/reference/getTrafo.md).
- Updated the `defineModel` help page.
- Added and updated pkgdown documentation.
- Added `test_defineModel.R` with TAM and ConQuest examples from the
  `defineModel` documentation.
- Added tests for Thurstone thresholds and similarity tests for `it1`.
- Reworked [`stopifnot()`](https://rdrr.io/r/base/stopifnot.html) checks
  into `testthat` tests.
- Updated tests to ignore rows with missing values in comparisons.
- Added Edna Grewers as a package contributor.

### 2025-01-13 to 2025-05-13

- Added a message for
  [`equat1pl()`](https://weirichs.github.io/eatModel/reference/equat1pl.md).
- Increased Q3 estimation speed by incorporating the Rfast package.
- Adapted
  [`transformToBista()`](https://weirichs.github.io/eatModel/reference/transformToBista.md)
  for non-nested person groups.
- Rewrote and clarified
  [`multiEquatError()`](https://weirichs.github.io/eatModel/reference/multiEquatError.md).
- Added assertions.
- Removed selected assertions for `fixSlopeMat` and directories.
- Removed recursive calls.
- Fixed splitting by non-nested person groups.
- Fixed partial credit definitions.
- Fixed `scipen` handling for R \>= 4.5.0.
- Fixed and clarified help pages for
  [`defineModel()`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  and
  [`multiEquatError()`](https://weirichs.github.io/eatModel/reference/multiEquatError.md).

### 2024

#### Package structure and maintenance

- Restructured R files into focused modules for
  [`defineModel()`](https://weirichs.github.io/eatModel/reference/defineModel.md),
  [`runModel()`](https://weirichs.github.io/eatModel/reference/runModel.md),
  [`getResults()`](https://weirichs.github.io/eatModel/reference/getResults.md),
  [`equat1pl()`](https://weirichs.github.io/eatModel/reference/equat1pl.md),
  [`transformToBista()`](https://weirichs.github.io/eatModel/reference/transformToBista.md),
  [`plotICC()`](https://weirichs.github.io/eatModel/reference/plotICC.md),
  ConQuest helpers, TAM helpers, result extractors, Q3 helpers, and
  model checks.
- Added `cli` support.
- Replaced deprecated
  [`memory.limit()`](https://rdrr.io/r/utils/memory.size.html) usage
  with
  [`ps::ps_system_memory()`](https://ps.r-lib.org/reference/ps_system_memory.html).
- Updated the `eatTools` dependency.
- Removed outdated `getConquestVersion()`, `prepRep()`, the `date`
  dependency, and obsolete exports.
- Changed encoding in `getResults.Rd` from ANSI to UTF-8.
- Updated README content, logo, sticker, and R-CMD badge.

#### Assertions and checks

- Added assertions for `verbose`, `digits`, `isLetter()`,
  [`simEquiTable()`](https://weirichs.github.io/eatModel/reference/simEquiTable.md),
  `getConquestVersion()`,
  [`checkContextVars()`](https://weirichs.github.io/eatModel/reference/checkContextVars.md),
  `checkQmatrixConsistency()`, `checkPersonGroupsConsistency()`,
  [`splitModels()`](https://weirichs.github.io/eatModel/reference/splitModels.md),
  [`checkLinking()`](https://weirichs.github.io/eatModel/reference/checkLinking.md),
  [`defineModel()`](https://weirichs.github.io/eatModel/reference/defineModel.md),
  [`runModel()`](https://weirichs.github.io/eatModel/reference/runModel.md),
  `anker()`, result extractor functions, low-level ConQuest readers,
  [`plotDevianceConquest()`](https://weirichs.github.io/eatModel/reference/plotDevianceConquest.md),
  [`multiEquatError()`](https://weirichs.github.io/eatModel/reference/multiEquatError.md),
  [`transformToBista()`](https://weirichs.github.io/eatModel/reference/transformToBista.md),
  and logical arguments.
- Moved Q-matrix and person group consistency checks into shared check
  functions.
- Reworked
  [`checkLinking()`](https://weirichs.github.io/eatModel/reference/checkLinking.md)
  diagnostics and removed redundant checks.
- Updated warnings and messages in item consistency checks.
- Added and revised checks for `cutScores`, model name elements,
  `qMatrix`, `idVarName`, `logFile`, `weights`, `analysis.name`, and
  `irtmodel`.

#### Improvements and bug fixes

- Added customizable model statements and fixed related bugs.
- Fixed DIF identification and ETS DIF classification.
- Fixed `getTamEAPs()`,
  [`get.shw()`](https://weirichs.github.io/eatModel/reference/get.shw.md),
  [`equat1pl()`](https://weirichs.github.io/eatModel/reference/equat1pl.md),
  [`splitModels()`](https://weirichs.github.io/eatModel/reference/splitModels.md),
  [`simEquiTable()`](https://weirichs.github.io/eatModel/reference/simEquiTable.md),
  `getTamWLE()`, and item consistency checks.
- Enhanced
  [`get.dsc()`](https://weirichs.github.io/eatModel/reference/get.dsc.md).
- Implemented `stringr` in
  [`get.shw()`](https://weirichs.github.io/eatModel/reference/get.shw.md).
- Replaced `set.col.type` for better performance.
- Updated `getResults` documentation and the
  [`defineModel()`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  help file.
- Tried to resolve R CMD check failures.

### 2023

- Added
  [`compareDefineModelObjects()`](https://weirichs.github.io/eatModel/reference/compareDefineModelObjects.md).
- Added a seed argument in
  [`getResults()`](https://weirichs.github.io/eatModel/reference/getResults.md)
  for infit computation.
- Added warnings in
  [`equat1pl()`](https://weirichs.github.io/eatModel/reference/equat1pl.md).
- Applied rounding in
  [`simEquiTable()`](https://weirichs.github.io/eatModel/reference/simEquiTable.md).
- Fixed
  [`runModel()`](https://weirichs.github.io/eatModel/reference/runModel.md)
  with `show.output.on.console`.
- Fixed
  [`checkLinking()`](https://weirichs.github.io/eatModel/reference/checkLinking.md),
  [`equat1pl()`](https://weirichs.github.io/eatModel/reference/equat1pl.md),
  [`multiEquatError()`](https://weirichs.github.io/eatModel/reference/multiEquatError.md),
  `tripleEquatError()`, and DIF estimation.

### 2022

#### New features and improvements

- Added
  [`regcoefFromRes()`](https://weirichs.github.io/eatModel/reference/regcoefFromRes.md).
- Added
  [`checkContextVars()`](https://weirichs.github.io/eatModel/reference/checkContextVars.md).
- Added optional long-format output in
  [`q3FromRes()`](https://weirichs.github.io/eatModel/reference/q3FromRes.md).
- Added Q3 pairs to
  [`transformToBista()`](https://weirichs.github.io/eatModel/reference/transformToBista.md)
  output.
- Added diagnostics and alternatives for
  [`checkLinking()`](https://weirichs.github.io/eatModel/reference/checkLinking.md).
- Added error messages to
  [`simEquiTable()`](https://weirichs.github.io/eatModel/reference/simEquiTable.md)
  and `checkItemParLists()`.
- Updated exemplary data and GitHub Actions.
- Replaced `cat("Warning...")` with `warning("...")`.
- Tidied syntax and removed unnecessary lines and `assign(...)` calls.

#### Bug fixes

- Fixed
  [`itemFromRes()`](https://weirichs.github.io/eatModel/reference/itemFromRes.md)
  and TAM DIF analysis handling.
- Fixed
  [`transformToBista()`](https://weirichs.github.io/eatModel/reference/transformToBista.md).
- Fixed linking error calculations.
- Fixed `doAufb()`.
- Fixed `tripleEquatError()`, `equat1p()`, and
  [`equat1pl()`](https://weirichs.github.io/eatModel/reference/equat1pl.md).
- Fixed
  [`q3FromRes()`](https://weirichs.github.io/eatModel/reference/q3FromRes.md).
- Repaired `blockPositions` and updated position names in
  [`checkLinking()`](https://weirichs.github.io/eatModel/reference/checkLinking.md).
- Updated
  [`checkLinking()`](https://weirichs.github.io/eatModel/reference/checkLinking.md)
  documentation.

### 2021

#### New features and improvements

- Set up GitHub Actions for continuous integration.
- Added WLE reliability and EAP reliability to the results object.
- Added
  [`multiEquatError()`](https://weirichs.github.io/eatModel/reference/multiEquatError.md)
  for linking three measurements and updated it.
- Added empty `testletStr` handling and checks for testlet structures.
- Added examples for trend analyses with three measurement points.
- Added output information in
  [`defineModel()`](https://weirichs.github.io/eatModel/reference/defineModel.md).
- Asked for the ConQuest executable location.
- Prepared for more than two measurement points.
- Modified dependencies and removed `readxl` from help files.
- Removed console files and remotes from the repository metadata.
- Added and corrected help file notes, examples, and typos.

#### Bug fixes

- Fixed `getTamResults()`, `getTam2plDiscrim()`,
  `checkQmatrixConsistency()`, `doAufb()`,
  [`get.shw()`](https://weirichs.github.io/eatModel/reference/get.shw.md),
  `prepareAndCheckEatModelObject()`,
  [`itemFromRes()`](https://weirichs.github.io/eatModel/reference/itemFromRes.md),
  and
  [`transformToBista()`](https://weirichs.github.io/eatModel/reference/transformToBista.md).
- Eliminated notes and warnings.

### 2020

- Released versions 0.7.9 through 0.7.30.
- Updated `eatTools` and `eatRep` dependencies.
- Removed the `eatTools` remotes section from `DESCRIPTION`.
- Modified examples.
- Fixed bugs, including a bugfix after version 0.7.30.

### 2019

#### Releases and new functionality

- Released versions 0.6.27 through 0.7.9.
- Added `checkDesign()`.
- Added `quasiMontecarlo` for TAM methods.
- Added 2PL parameter handling from TAM when `omitRegr = TRUE`.
- Added remotes metadata.
- Updated
  [`plotICC()`](https://weirichs.github.io/eatModel/reference/plotICC.md)
  and defaults for
  [`plotICC()`](https://weirichs.github.io/eatModel/reference/plotICC.md).
- Added hints for TAM compatibility.

#### Bug fixes and maintenance

- Fixed
  [`runModel()`](https://weirichs.github.io/eatModel/reference/runModel.md)
  for 2PL models.
- Updated messages in `anker()` and
  [`transformToBista()`](https://weirichs.github.io/eatModel/reference/transformToBista.md).
- Fixed `checkDesign()`, `checkModel()`, and the set method for
  [`defineModel()`](https://weirichs.github.io/eatModel/reference/defineModel.md).
- Fixed model splitting.
- Updated `DESCRIPTION` and uploaded package files.

### 2018

- Initial commit on 2018-10-30.
- Released versions 0.6.23 through 0.6.26.
- Updated `eatModel.r`.
