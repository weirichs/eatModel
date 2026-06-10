# eatModel 0.10.19

Changes based on available commit messages. The local Git history can be traced
back to the initial commit on 2018-10-30.

## 2025-08-25 to 2026-06-10

### New features

* Added partial credit support for TAM, mirt, and ConQuest.
* Added DIF handling for partial credit models in TAM and ConQuest.
* Added partial credit support to `plotICC()`.
* Added `getTrafo()` to extract IQB transformation rules from a database.
* Added an equivalence table for partial credit models.

### Improvements

* Improved WLE handling in multidimensional models by accounting for missing
  values on subdimensions.
* Updated `adaptSkeletonForAnchor()` so covariances are not constrained.
* Omit empty categories.
* Allow ConQuest variable names longer than 10 characters.
* Restrict `group.var` to `software = "conquest"`.

### Bug fixes

* Fixed `getDevianceConquest()`.
* Fixed `getConquestDesc()`.
* Fixed `runModel()` for partial credit models in mirt.
* Fixed `checkGenerateRefPop()`.
* Fixed issues in Thurstone threshold tests.

### Documentation and tests

* Added documentation for `getTrafo()`.
* Updated the `defineModel` help page.
* Added and updated pkgdown documentation.
* Added `test_defineModel.R` with TAM and ConQuest examples from the
  `defineModel` documentation.
* Added tests for Thurstone thresholds and similarity tests for `it1`.
* Reworked `stopifnot()` checks into `testthat` tests.
* Updated tests to ignore rows with missing values in comparisons.
* Added Edna Grewers as a package contributor.

## 2025-01-13 to 2025-05-13

* Added a message for `equat1pl()`.
* Increased Q3 estimation speed by incorporating the Rfast package.
* Adapted `transformToBista()` for non-nested person groups.
* Rewrote and clarified `multiEquatError()`.
* Added assertions.
* Removed selected assertions for `fixSlopeMat` and directories.
* Removed recursive calls.
* Fixed splitting by non-nested person groups.
* Fixed partial credit definitions.
* Fixed `scipen` handling for R >= 4.5.0.
* Fixed and clarified help pages for `defineModel()` and `multiEquatError()`.

## 2024

### Package structure and maintenance

* Restructured R files into focused modules for `defineModel()`, `runModel()`,
  `getResults()`, `equat1pl()`, `transformToBista()`, `plotICC()`, ConQuest
  helpers, TAM helpers, result extractors, Q3 helpers, and model checks.
* Added `cli` support.
* Replaced deprecated `memory.limit()` usage with `ps::ps_system_memory()`.
* Updated the `eatTools` dependency.
* Removed outdated `getConquestVersion()`, `prepRep()`, the `date` dependency,
  and obsolete exports.
* Changed encoding in `getResults.Rd` from ANSI to UTF-8.
* Updated README content, logo, sticker, and R-CMD badge.

### Assertions and checks

* Added assertions for `verbose`, `digits`, `isLetter()`, `simEquiTable()`,
  `getConquestVersion()`, `checkContextVars()`, `checkQmatrixConsistency()`,
  `checkPersonGroupsConsistency()`, `splitModels()`, `checkLinking()`,
  `defineModel()`, `runModel()`, `anker()`, result extractor functions,
  low-level ConQuest readers, `plotDevianceConquest()`, `multiEquatError()`,
  `transformToBista()`, and logical arguments.
* Moved Q-matrix and person group consistency checks into shared check
  functions.
* Reworked `checkLinking()` diagnostics and removed redundant checks.
* Updated warnings and messages in item consistency checks.
* Added and revised checks for `cutScores`, model name elements, `qMatrix`,
  `idVarName`, `logFile`, `weights`, `analysis.name`, and `irtmodel`.

### Improvements and bug fixes

* Added customizable model statements and fixed related bugs.
* Fixed DIF identification and ETS DIF classification.
* Fixed `getTamEAPs()`, `get.shw()`, `equat1pl()`, `splitModels()`,
  `simEquiTable()`, `getTamWLE()`, and item consistency checks.
* Enhanced `get.dsc()`.
* Implemented `stringr` in `get.shw()`.
* Replaced `set.col.type` for better performance.
* Updated `getResults` documentation and the `defineModel()` help file.
* Tried to resolve R CMD check failures.

## 2023

* Added `compareDefineModelObjects()`.
* Added a seed argument in `getResults()` for infit computation.
* Added warnings in `equat1pl()`.
* Applied rounding in `simEquiTable()`.
* Fixed `runModel()` with `show.output.on.console`.
* Fixed `checkLinking()`, `equat1pl()`, `multiEquatError()`,
  `tripleEquatError()`, and DIF estimation.

## 2022

### New features and improvements

* Added `regcoefFromRes()`.
* Added `checkContextVars()`.
* Added optional long-format output in `q3FromRes()`.
* Added Q3 pairs to `transformToBista()` output.
* Added diagnostics and alternatives for `checkLinking()`.
* Added error messages to `simEquiTable()` and `checkItemParLists()`.
* Updated exemplary data and GitHub Actions.
* Replaced `cat("Warning...")` with `warning("...")`.
* Tidied syntax and removed unnecessary lines and `assign(...)` calls.

### Bug fixes

* Fixed `itemFromRes()` and TAM DIF analysis handling.
* Fixed `transformToBista()`.
* Fixed linking error calculations.
* Fixed `doAufb()`.
* Fixed `tripleEquatError()`, `equat1p()`, and `equat1pl()`.
* Fixed `q3FromRes()`.
* Repaired `blockPositions` and updated position names in `checkLinking()`.
* Updated `checkLinking()` documentation.

## 2021

### New features and improvements

* Set up GitHub Actions for continuous integration.
* Added WLE reliability and EAP reliability to the results object.
* Added `multiEquatError()` for linking three measurements and updated it.
* Added empty `testletStr` handling and checks for testlet structures.
* Added examples for trend analyses with three measurement points.
* Added output information in `defineModel()`.
* Asked for the ConQuest executable location.
* Prepared for more than two measurement points.
* Modified dependencies and removed `readxl` from help files.
* Removed console files and remotes from the repository metadata.
* Added and corrected help file notes, examples, and typos.

### Bug fixes

* Fixed `getTamResults()`, `getTam2plDiscrim()`,
  `checkQmatrixConsistency()`, `doAufb()`, `get.shw()`,
  `prepareAndCheckEatModelObject()`, `itemFromRes()`, and
  `transformToBista()`.
* Eliminated notes and warnings.

## 2020

* Released versions 0.7.9 through 0.7.30.
* Updated `eatTools` and `eatRep` dependencies.
* Removed the `eatTools` remotes section from `DESCRIPTION`.
* Modified examples.
* Fixed bugs, including a bugfix after version 0.7.30.

## 2019

### Releases and new functionality

* Released versions 0.6.27 through 0.7.9.
* Added `checkDesign()`.
* Added `quasiMontecarlo` for TAM methods.
* Added 2PL parameter handling from TAM when `omitRegr = TRUE`.
* Added remotes metadata.
* Updated `plotICC()` and defaults for `plotICC()`.
* Added hints for TAM compatibility.

### Bug fixes and maintenance

* Fixed `runModel()` for 2PL models.
* Updated messages in `anker()` and `transformToBista()`.
* Fixed `checkDesign()`, `checkModel()`, and the set method for
  `defineModel()`.
* Fixed model splitting.
* Updated `DESCRIPTION` and uploaded package files.

## 2018

* Initial commit on 2018-10-30.
* Released versions 0.6.23 through 0.6.26.
* Updated `eatModel.r`.
