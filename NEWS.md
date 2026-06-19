# eatModel 0.10.19 [Jun 2026]

* describe renam argument in `checkContextVars()`
* minor bug fix in `runModel()` for TAM: tam function selection between 1pl and 2pl is correct now
* Use `mirt::traditional2mirt()` for transformation of item parameters.
* implement 2pl transformation for item parameters according to a probability of 62.5% 

# eatModel 0.10.18 [Jun 2026]

* add help file for `getTrafo()`
* add partial credit to `plotICC()`

# eatModel 0.10.17 [May 2026]

* dif for partial credit in tam
* add function `getTrafo()` for extraction of IQB transformation rules from data base

# eatModel 0.10.10 [May 2026]

* added pkgdown page

# eatModel 0.10.9 [Apr 2026]

* updated tests to ignore NA rows in the comparison. 
* bugfix in `runModel()` for partial credit in mirt: column identification of 
  irtmodel via number instead colname
* defineModel.rd help page: more streamlined examples that run faster

# eatModel 0.10.8 [Mar 2026]

* add tests for thurstone thresholds
* added new test file test_defineModel.R, copied TAM and Conquest examples from 
  .Rd file (examples 8, 8.1), renamed objects for comparison
* changed the stopifnot lines to test_that tests

# eatModel 0.10.5 [Feb 2026]

* allow more than 10 characters in conquest variable names

# eatModel 0.10.4 [Dec 2025]

* add partial credit for conquest

# eatModel 0.10.3 [Nov 2025]

* do not constrain covariances in `adaptSkeletonForAnchor()`
* WLEs in multidimensional models for mirt: account for missings on sub dimensions

# eatModel 0.10.1 [Nov 2025]

* partial credit with TAM and mirt

# eatModel 0.9.5 [Aug 2025]

* use the benchmarkme package for `getDevianceConquest()`

# eatModel 0.9.4 [May 2025]

* bugfix in scipen (adapted for R version >= 4.5.0)

# eatModel 0.9.0 [Mar 2025]

* `multiEquatError()`: rewritten: tests added, help files adapted
* remove assertions in `fixSlopeMat(): In the assertions, the columns must not 
  be selected according to order, as this can vary.
* remove recursive calls in `defineModel()`: new internally needed function 
  `defineModelSingle()`

# eatModel 0.8.52 [Mar 2025]

* Fix in `multiEquatError()`: The old function lacked a check if the testlets are 
  constant. In addition, the comparison of the original and chained linking errors 
  was misleading. Both should now be corrected.

# eatModel 0.8.49 [Jan 2025]

* adapt `transformToBista()` for non-nested person groups

# eatModel 0.8.43 [Jan 2025]

* add message for `equat1pl()` that testlet column name and dimension names must be different
* increase Q3 estimation speed by incorporating Rfast package

# eatModel 0.8.37 [Nov 2024]

* `transformToBista()`: Bug fix in the weights argument

# eatModel 0.8.36 [Oct 2024]

* deleted `getConquestVersion()` (outdated) and `prepRep()` (now via eatGADS).
* deleted date dependency in DESCRIPTION
* `simEquiTable()`: added new checks for cutScores and applied feedback
* use `makeDateFrame()` for `checkPersonGroupsConsistency()`
* `anker()`: added column types for prm
* assertions for all `...FromRes()` functions
* assertions for all `get...()` functions
* assertions for `plotDevianceConquest()`
* `multiEquatError()`: bugfix testletStr

# eatModel 0.8.35 [Sep 2024]

* use cli for error/warning message presentation
* The functions are no longer contained in a single file (eatModel.r), but are 
  grouped according to their intended use. 
* `getConquestVersion()`, `simEquiTable()`: new assertions (input validation) 
* getResults.rd: change file encoding from ISO-8859-1 to UTF-8
* `checkQmatrixConsistency()`: moved checks from `splitModels()` in here (because 
  the function is also called by defineModel
* added checkmate checks for `defineModel()`
* remove problematic `eatTools::gsubAll()` from `checkLinking()` 

# eatModel 0.8.34 [Aug 2024]

* changing encoding from ANSI to UTF-8 in getResults.rd
* replace deprecated `memory.limit()` with `ps::ps_system_memory()`
* `get.shw()`: Conquest only allows variable names with a maximum of 11 characters. 
  If the variable name is exactly 11 characters long, the space between the variable 
  name and the item parameter value is missing in the ASCII file. In this case, the 
  columns could not be read correctly because this was done by identifying spaces 
  in the ASCII file. This bug has been fixed here. 

# eatModel 0.8.21 [Feb 2023]

* Bug fix for rounding in `simEquiTable()`: Depending on the number of decimal 
  places to be displayed, a less ambiguous rounding method is now selected, 
  which makes it possible to determine exactly which proficiency level an item 
  falls into when the transformed item parameter lies on the boundary between 
  proficiency levels. 
* add warnings in `equat1pl()` if dimension names do not correspond to refPop's 
  first column names
* add seed argument in `getResults()` (for infit computation)
* added assertions in `replaceLinkingError()`: verbose, digits
* add package sticker
* customize model statement in `defineModel()` for Conquest: The function checks 
  whether the variable names in the `model` statement are also included in the 
  dataset. However, this does not apply to reserved names such as "item,",  "step," etc. 

# eatModel 0.8.14 [Jan 2023]

* add `compareDefineModelObjects()` which facilitates quality checks by 
  comparing two objects returned by `defineModel()` whether there are 
  functionally equivalent.

# eatModel 0.8.13 [Oct 2022]

* add q3 pairs in `transformtobista()` output
* update exemplary data (trends.rda)

# eatModel 0.8.12 [Sep 2022]

* add exemplary conquest results files for demonstration purposes
* add appropriate error messages in `checkItemParLists()` when parameter 
  value column is not numeric
* provide optional long format return in `q3FromRes()`

# eatModel 0.8.10 [May 2022]

* add error messages to `simEquiTable()`  

# eatModel 0.8.1 [May 2022]

* The function `mergeAttr()` is being moved from the eatModel package 
  to the eatTools package. 

# eatModel 0.8.0 [Feb 2022]

* add function `checkContextVars()` which checks whether some types of 
  context vars, i.e. group, DIF and weighting variables, are consistent 
  with item variables. The function is mainly used for internal consistency checks.

# eatModel 0.7.57 [Feb 2022]

* add function `regcoefFromRes()` which lists the regression parameters 
  of the conditioning model from the results object. 

# eatModel 0.7.51 [Feb 2022]

* `transformToBista()` now gives linking errors in a data.frame which 
  is appropriate for eatRep

# eatModel 0.7.48 [Jan 2022]

* Bug fix in `transformToBista()`: drawing PVs failed when items are 
  only partially linked

# eatModel 0.7.47 [Dec 2021]

* add examples for trend analyses with 3 measurement occasions
* `defineModel()`: add output information about students with all items 
  solved or all items failed

# eatModel 0.7.45 [Nov 2021]

* add checks for testlet structure
* additional examples in `multiEquatError()` help file

# eatModel 0.7.39 [Oct 2021]

* add functionality for >2 measurement occasions
* Extracting the recursive function part from `defineModel()` which is
  now provided in a new function `doAufb()`
* add `multiEquatError()` for linking of three measurement  

# eatModel 0.7.35 [Jul 2021]

* add EAP reliability in the results object
* add checks in `checkQmatrixConsistency()`: Indicator columns must not 
  be consistently 0 (a consistent 1 would be allowed; this would then 
  constitute within-item multidimensionality). 

# eatModel 0.7.34 [Mar 2021]

* Bug fix in `getTam2plDiscrim()`: This internally used function reloads 
  the 2pl discrimination parameters from TAM. 
* `wleRelFromRes()` adds wle reliability to results object, also for 
  multidimensionality models
* help files for test design data provided as rda instead of xlsx  

# eatModel 0.7.33 [Mar 2021]

* If it cannot be found on the network drives, the user will be prompted 
  to specify the location of the Conquest console the first time it is needed. 

# eatModel 0.7.27 [Nov 2020]

* Add validity checks to verify that background variables are consistent 
  with themselves and with the item data. Also, determine the length 
  (number of characters required) for each variable. 

# eatModel 0.7.25 [Nov 2020]

* set compatibility to eatRep 0.13.0: replace `jk2.mean()` by `repMean()` 

# eatModel 0.7.24 [Aug 2020]

* Bug fix: Importing WLEs from Conquest now also works for multidimensional cases. 

# eatModel 0.7.21 [Jul 2020]

* set compatibility to eatTools 0.3.1

# eatModel 0.7.14 [Mar 2020]

* add ordinary least squares (OLS) to robust linking methods in 
  `equat1pl()`

# eatModel 0.7.13 [Feb 2020]

* `transformToBista()` calculates the mean and standard deviation of the 
  reference population using eatRep methods when the given sample is drawn 
  from the norm population and is not to be equated. 

# eatModel 0.7.11 [Feb 2020]

* Add Haberman linking to `equat1pl()` 

# eatModel 0.7.10 [Jan 2020]

* Add robust methods to `equat1pl()` linking procedures 

# eatModel 0.7.8 [Nov 2019]

* Based on an idea by Johannes Schult, the `simEquiTable()` function now 
  uses methods from Ivailo Partchev's irtoys package to generate the item 
  parameter equivalence table more quickly and efficiently. 

# eatModel 0.7.7 [Nov 2019]

* Temporary storage of console messages that cannot be displayed immediately 
  in the case of multicore. They are stored as strings and displayed on the 
  console once the multicore analysis is complete. 

# eatModel 0.7.6 [Nov 2019]

* Bug fix for recursive calls to `defineModel()` in cases involving multiple models 

# eatModel 0.7.4 [May 2019]

* compatibility issues: Add functionality that allow the ID variable to 
  be specified manually when output generated by earlier package versions 
  needs to be re-read by newer package versions. 

# eatModel 0.7.2 [Apr 2019]

* Extractor functions for Conquest. 

# eatModel 0.7.0 [Mar 2019]

* Added `quasiMontecarlo` for TAM methods.
* Fixed `runModel()` for 2PL models.

# eatModel 0.6.33 [Feb 2019]

* add `checkDesign()` for verifying the underlying test design

# eatModel 0.6.32 [Feb 2019]

* 2pl parameter extraction 

# eatModel 0.6.30 [Jan 2019]

* plotICC() for dichotomous 1pl item characteristic curves 

# eatModel 0.6.26 [Oct 2018]

* Initial release on 2018-10-30.
* Released versions 0.6.23 through 0.6.26.
* Updated `eatModel.r`.
