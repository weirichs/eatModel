# Collect all results from Conquest, TAM, or mirt analyses into a common data frame

First the IRT model should be defined using
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md).
Afterwards, call
[`runModel`](https://weirichs.github.io/eatModel/reference/runModel.md)
with the argument returned by
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
to start the estimation. The last step then is to create a results frame
using `getResults`.

## Usage

``` r
getResults( runModelObj, overwrite = FALSE, Q3 = TRUE, q3theta = c("pv", "wle", "eap"), 
            q3MinObs = 0, q3MinType = c("singleObs", "marginalSum"), omitFit = FALSE, 
            omitRegr = FALSE, omitWle = FALSE, omitPV = FALSE, abs.dif.bound = 0.6, 
            sig.dif.bound = 0.3, p.value = 0.9, nplausible = NULL, ntheta = 2000,
            normal.approx = FALSE, samp.regr = FALSE, theta.model=FALSE, np.adj=8,
            group = NULL, beta_groups = TRUE, level = .95, n.iter = 1000,
            n.burnin = 500, adj_MH = .5, adj_change_MH = .05, refresh_MH = 50, 
            accrate_bound_MH = c(.45, .55),  sample_integers=FALSE, theta_init=NULL,
            print_iter = 20, verbose = TRUE, calc_ic=TRUE, omitUntil = 1, seed=NA)
```

## Arguments

- runModelObj:

  The object returned by
  [`runModel`](https://weirichs.github.io/eatModel/reference/runModel.md).

- overwrite:

  Logical. Should result files be overwritten if exist?

- Q3:

  Logical. Estimate the Q3 statistic according to Yen (1984)? Note: this
  is only possible for uni-dimensional models. If `software == "tam"`,
  Q3 statistic is estimated using the
  [`tam.modelfit`](https://rdrr.io/pkg/TAM/man/tam.modelfit.html)
  function. If `software == "Conquest"`, Q3 statistic is estimated using
  the [`Q3`](https://rdrr.io/pkg/sirt/man/Q3.html) function from the
  `sirt` package. Note that Q3 estimation does not work yet if
  `software == "mirt"`.

- q3theta:

  Specify whether the Q3 statistic should be estimated using PVs, WLEs
  or EAPs as the theta variable.

- q3MinObs:

  Q3 statistic might be untrustworthy if item covariance estimation is
  based on very few observations. Define the minimum number of
  observation which should be fulfilled for Q3 estimation.

- q3MinType:

  If `"singleObs"`, `q3MinObs` argument is based on the least number of
  observations in the \\2\times 2\\ 0/1 frequency table of item pairs.
  If `"marginalSum"`, `q3MinObs` argument is based on the sum of
  marginals in the \\2\times 2\\ 0/1 frequency table of item pairs.

- omitFit:

  Logical. Should item fit values be included into the results?

- omitRegr:

  Logical. Should regression parameters and their standard errors be
  included into the results?

- omitWle:

  Logical. Should WLE estimates be included into the results?

- omitPV:

  Logical. Should plausible values be included into the results?

- abs.dif.bound:

  Applies only if DIF analyses are performed before. When DIF-Parameter
  are evaluated, this specifies the critical value for absolute DIF. See
  the details section for further information.

- sig.dif.bound:

  Applies only if DIF analyses are performed before. When DIF-Parameter
  are evaluated, this specifies the critical value for confidence
  interval DIF. See the details section for further information.

- p.value:

  Applies only if DIF analyses are performed before. When DIF-Parameter
  are evaluated, this specifies the critical p-value for confidence
  interval DIF. See the details section for further information.

- nplausible:

  Applies only if `software = "tam"`: Number of plausible values to be
  drawn. Note: number of plausible values were already defined in
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md),
  because Conquest needs to know the number of PVs prior to estimation.
  In `TAM`, it is possible to redefine the number of plausible values
  and overwrite the definition that was given in
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md).

- ntheta:

  Applies only if `software = "tam"`. Following description is borrowed
  from the help file of
  [`tam.pv`](https://rdrr.io/pkg/TAM/man/tam.pv.html) from the `TAM`
  package: Number of ability nodes for plausible value imputation. Note
  that in this function ability nodes are simulated for the whole
  sample, not for every person (contrary to the software Conquest).

- normal.approx:

  Applies only if `software = "tam"`. Following description is borrowed
  from the help file of
  [`tam.pv`](https://rdrr.io/pkg/TAM/man/tam.pv.html) from the `TAM`
  package: An optional logical indicating whether the individual
  posterior distributions should be approximated by a normal
  distribution? The default is `FALSE`. In the case `normal.approx=TRUE`
  (normal distribution approximation), the number of ability nodes
  `ntheta` can be substantially smaller than 2000, say 200 or 500. The
  normal approximation is implemented for unidimensional and
  multidimensional models.

- samp.regr:

  Applies only if `software = "tam"`. Following description is borrowed
  from the help file of
  [`tam.pv`](https://rdrr.io/pkg/TAM/man/tam.pv.html) from the `TAM`
  package: An optional logical indicating whether regression
  coefficients should be fixed in the plausible value imputation or also
  sampled from their posterior distribution? The default is `FALSE`.
  Sampled regression coefficients are obtained by nonparametric
  bootstrap.

- theta.model:

  Applies only if `software = "tam"`. Following description is borrowed
  from the help file of
  [`tam.pv`](https://rdrr.io/pkg/TAM/man/tam.pv.html) from the `TAM`
  package: Logical indicating whether the theta grid from the `tamobj`
  object should be used for plausible value imputation. In case of
  `normal.approx=TRUE`, this should be sufficient in many applications.

- np.adj:

  Applies only if `software = "tam"`. Following description is borrowed
  from the help file of
  [`tam.pv`](https://rdrr.io/pkg/TAM/man/tam.pv.html) from the `TAM`
  package: This parameter defines the “spread” of the random theta
  values for drawing plausible values when `normal.approx=FALSE`. If
  \\s\_{EAP}\\ denotes the standard deviation of the posterior
  distribution of theta (in the one-dimensional case), then theta is
  simulated from a normal distribution with standard deviation `np.adj`
  times \\s\_{EAP}\\.

- group:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`.
  Optional vector of group identifiers. See the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- beta_groups:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. See
  the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- level:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`.
  Confidence level in bayesian approach. See the help page
  of[`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- n.iter:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. Number
  of iterations in the bayesian approach. See the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- n.burnin:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. Number
  of burn-in iterations in the bayesian approach. See the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- adj_MH:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. See
  the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- adj_change_MH:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. See
  the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- refresh_MH:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. See
  the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- accrate_bound_MH:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. See
  the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- sample_integers:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`.
  Logical indicating whether weights for complete cases should be
  sampled in bootstrap. See the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- theta_init:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`.
  Optional matrix with initial theta values. See the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- print_iter:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. See
  the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- verbose:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`. See
  the help page of
  [`tam.pv.mcmc`](https://rdrr.io/pkg/TAM/man/tam.pv.html) for further
  details.

- calc_ic:

  Applies only if `software = "tam"` and `pvMethod = "bayesian"`.
  Logical indicating whether information criteria should be computed.
  See the help page of `tam.pv.mcmc` for further details.

- omitUntil:

  Argument is passed to
  [`plotDevianceConquest`](https://weirichs.github.io/eatModel/reference/plotDevianceConquest.md):
  An optional value indicating number of iterations to be omitted for
  plotting.

- seed:

  Fixed simulation seed. This value is directly passed on to the
  [`tam.fit`](https://rdrr.io/pkg/TAM/man/tam.fit.html) function.

## Details

If
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
was run with software Conquest, a path argument (`'dir'`) is necessary.
The path argument is optional for software TAM. If `'dir'` was
specified, `getResults` additionally writes its output into the
specified folder, using the `analysis.name` argument for file naming.
Otherwise, `getResults` only returnes the result data frame.

If DIF analyses were performed before, the user can specify the criteria
according to which DIF should be interpreted or evaluated. By default,
the ETS criteria (Zieky, 1993) are used which classify DIF into three
distinct categories, "A", "B", or "C". Small DIF ("A") corresponds to
absolute DIF values below .43 (no significance test is performed here);
medium DIF ("B") corresponds to absolute DIF values between .43 and .64
which are significantly different from zero. High DIF ("C") corresponds
to absolute DIF values above .64 which are significantly different from
.43 (DeMars, 2011; Monahan et al. 2007). Alternatively, the three
arguments `abs.dif.bound`, `sig.dif.bound`, and `p.value` allow to
specify user-defined dichotomous criteria. If items should be flagged as
DIF, if the absolute value increases 0.5 and significantly exceeds 0.1
at a alpha level of 0.05, use `abs.dif.bound = 0.5` and
`sig.dif.bound = 0.1` and `p.value = 0.95`.

## Value

A data frame in the long format with ten columns.

- model:

  The name of the model (as specified by the user in `analysis.name`.

- source:

  The estimation software (i.e, conquest or TAM)

- var1:

  The variable name for which the corresponding value is given, i.e. its
  indicator.

- var2:

  Additional variable information if necessary.

- type:

  Type of coefficient (for example, random or fixed).

- indicator.group:

  The type of the group the corresponding variable belongs to.

- group:

  The group the corresponding variable belongs to. Note: group is nested
  within `indicator.group`.

- par:

  The type of the parameter.

- derived.par:

  Optionally: The derived parameter.

- value:

  The value of the corresponding estimate.

## References

DeMars, C. E. (2011). An analytic comparison of effect sizes for
differential item functioning. Applied Measurement in Education, 24 (3),
189-209. https://doi.org/10.1080/08957347.2011.580255

Monahan, P. O., McHorney, C. A., Stump, T. E. & Perkins, A. J. (2007).
Odds ratio, delta, ETS classification, and standardization measures of
DIF magnitude for binary logistic regression. Journal of Educational and
Behavioral Statistics, 32 (1), 92-109.
https://doi.org/10.3102/1076998606298035

Zieky, M. (1993). Practical questions in the use of DIF statistics in
item development. In P. W. Holland & H. Wainer (Eds.), *Differential
item functioning* (pp. 337-347). Hillsdale, NJ: Lawrence Erlbaum.

## Examples

``` r
# see examples in the help file of defineModel()
```
