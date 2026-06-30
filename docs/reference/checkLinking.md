# check for linked design

Function checks whether all blocks in a complete or incomplete block
design are linked to each other.

## Usage

``` r
checkLinking(design, blocks=NULL, bookletColumn=NULL, verbose=FALSE)
```

## Arguments

- design:

  A data frame with the test design. All columns (except for the
  `booklet` identifier column) are expected to contain blocks.

- blocks:

  Optional. To check whether a subdomain is completely linked add here
  all the blocks that belong to this subdomain in a character vector.

- bookletColumn:

  Optional. Number or name of the booklet identifier column in the
  design data frame.

- verbose:

  Optional. If `TRUE` the function gives more messages.

## Details

## Value

A list. Containing information about linking diagnostics:

- completelyLinked:

  A logical value (TRUE/FALSE). Whether design is completely linked.

- occuringBlockCombinations:

  A data.frame containing the occureng frequencies of block pairs.

- blockPositions:

  A data.frame containing information about how often each block occurs
  at each position.

## Author

Sebastian Weirich and Karoline Sachse

## Examples

``` r
# 1. first examples
# a) design linked
des1a <- data.frame(booklet = paste0("B", 1:4),
       Pos1 = c("blockA", "blockB", "blockC", "blockD"),
       Pos2 = c("blockB", "blockC", "blockD", "blockE"),
       Pos3 = c("blockC", "blockD", "blockE", "blockF"))
test1 <- checkLinking(design = des1a, bookletColumn = "booklet", verbose=TRUE)
#> Dataset is completely linked.

# b) design not linked:
des1b <- data.frame(Pos1 = c("blockA", "blockH", "blockC", "blockF"),
       Pos2 = c("blockB", "blockC", "blockD", "blockA"),
       Pos3 = c("blockF", "blockG", "blockH", "blockB"))
test2 <- checkLinking(design = des1b, verbose=TRUE)
#> WARNING! Dataset is not completely linked.
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   24 cases unconnected. Following items are unconnected: 
#>                   blockH_1, blockH_2, blockH_3, blockC_1, blockC_2, blockC_3, blockG_1, blockG_2, blockG_3, blockD_1, blockD_2, blockD_3

# 2. second examples: use a 'IQB Bildungstrend 2016'-like design
data(des2)

# design contains three dimensions -- reading, listening and orthography
# linking check must be conducted for each dimension separately.

# a) domain reading (blocks contain "-L")
readblocks <- grep("-L", unique(unlist(des2[,-1])), value=TRUE)
test3 <- checkLinking(design = des2, blocks =readblocks, bookletColumn = "TH", verbose=TRUE)
#> Special characters and spaces will be removed from names in 'blocks'
#> WARNING! Dataset is not completely linked.
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   48 cases unconnected. Following items are unconnected: 
#>                   DL32_1, DL32_2, DL32_3

# b) domain listening (blocks contain "-H")
listenblocks <- grep("-H", unique(unlist(des2[,-1])), value=TRUE)
test4 <- checkLinking(design = des2[,-1], blocks =listenblocks, verbose=TRUE)
#> Special characters and spaces will be removed from names in 'blocks'
#> Dataset is completely linked.

# c) domain orthography (blocks contain "-R")
orthoblocks <- grep("-R", unique(unlist(des2[,-1])), value=TRUE)
test5 <- checkLinking(design = des2[,-1], blocks =orthoblocks, verbose=TRUE)
#> Special characters and spaces will be removed from names in 'blocks'
#> WARNING! Dataset is not completely linked.
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   156 cases unconnected. Following items are unconnected: 
#>                   DR13_1, DR13_2, DR13_3, DR01_1, DR01_2, DR01_3, DR23_1, DR23_2, DR23_3, DR0516S_1, DR0516S_2, DR0516S_3, DR0416S_1, DR0416S_2, DR0416S_3, DR02_1, DR02_2, DR02_3, DR12_1, DR12_2, DR12_3, DR33_1, DR33_2, DR33_3, DR21_1, DR21_2, DR21_3

# reconstruct test design from exemplary data, separately for each year and each domain
data(trends)
design<- by(data = trends, INDICES = trends[,c("year", "domain")], FUN = function (d) {
         rownames(d) <- NULL
         dw <- reshape2::dcast(unique(d[,c("booklet", "block", "pos")]), booklet~pos, value.var="block")
         message("\nCondition: \n", eatTools::print_and_capture(d[1,c("year", "domain")]))
         cl <- checkLinking(design=dw, bookletColumn ="booklet") } )
#> 
#> Condition: 
#>   year    domain
#> 1 2010 listening
#> Dataset is completely linked.
#> 
#> Condition: 
#>   year    domain
#> 1 2015 listening
#> WARNING! Dataset is not completely linked.
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   96 cases unconnected. Following items are unconnected: 
#>                   Bl12_1, Bl12_2, Bl12_3, Bl09_1, Bl09_2, Bl09_3
#> 
#> Condition: 
#>   year    domain
#> 1 2020 listening
#> Dataset is completely linked.
#> 
#> Condition: 
#>   year  domain
#> 1 2010 reading
#> Dataset is completely linked.
#> 
#> Condition: 
#>   year  domain
#> 1 2015 reading
#> WARNING! Dataset is not completely linked.
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   48 cases unconnected. Following items are unconnected: 
#>                   Bl25_1, Bl25_2, Bl25_3
#> 
#> Condition: 
#>   year  domain
#> 1 2020 reading
#> Dataset is completely linked.
```
