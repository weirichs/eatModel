# Dichotomous and partial credit responses for reading comprehension data obtained from large-scale assessment

This data set contains fictional achievement scores of 2595 students in
the long format along with item and person covariates. The data set is
based on an actual pilot survey from 2012, but has been anonymized for
use in the package.

## Usage

``` r
data(reading)
```

## Format

'data.frame': 64830 obs. of 22 variables:

- taskID:

  task ID

- bookletID:

  booklet ID

- idstud:

  unique student identifier

- sex:

  student's sex

- language:

  student's language spoken at home

- age:

  student's age in years

- country:

  the country the student stems from (anonymized)

- format:

  response format of the item

- type:

  item type: The items that were part of the pilot study are labeled
  'pilot'. The link items that provide the connection to the educational
  standards metric are called 'link'. Items from the IGLU study are
  called 'iglu'.

- taskLabel:

  short description of the task

- textType:

  brief description of the stimulus text on which the test tasks were
  based

- position:

  item position (i.e., block position)

- item:

  item identifier

- valueSum:

  The sum score for each student of all variables belonging to an item.
  Use this variable for partial credit model, i.e., polytomous
  responses.

- valueAgg:

  Aggregated item score. Use this variable for dichotomous responses.

- valueMax:

  maximum score which was possible for this item

## Source

pilot study of a reading comprehension test of the Institute of
Educational Progress
