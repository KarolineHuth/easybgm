# easybgm 0.5.1

## Bug fixes

* **Changes returned values.** Pairwise estimates from bgms are now placed on
  the edge they belong to. For mixed variable types bgms groups pairwise
  quantities by variable-type pair rather than in triangle order, and easybgm
  placed them by position, so whenever variable types were interleaved the
  estimates landed on the wrong edges. The matrices easybgm builds from them
  are now filled by pair name: `parameters` from `easybgm()`, and
  `parameters`, `parameters_g1`, `parameters_g2` and `overall_estimate` from
  `easybgm_compare()`. Models with a single variable type are unaffected.
* Inclusion Bayes factors (`inc_BF`, and the `MCSE_BF` interval derived from
  it) for bgms fits are now taken from `bgms::extract_inclusion_bf()`. Under a
  stochastic block prior with different within- and between-block
  hyperparameters, easybgm built the prior odds from the estimated posterior
  partition instead of marginalizing over the prior partition, which could
  misstate the Bayes factor by more than an order of magnitude. With bgms
  older than 0.2.0.0, which lacks the extractor, the previous calculation is
  kept.
* **Changes returned values.** A bgms comparison of more than two groups
  (fitted with `group_indicator`) no longer returns `parameters`. It held the average of the bgms contrast
  coefficients of each edge, which depends on the contrast basis, is not a
  group difference, and could have the opposite sign to every pairwise group
  difference. It is replaced by `pairwise_group_differences` (one column per
  pair of groups, e.g. `"group2 - group1"`) and `contrast_coefficients` (the
  coefficients as bgms reports them, labelled `"edge (diffN)"`). `summary()`
  no longer shows the "Average Difference" column, and `plot_network()` stops
  with an error for these fits.
* For a two-group bgms comparison, `parameters` is now documented as group 2
  minus group 1.
* **Changes returned values.** `overall_estimate`, shown as "Across-group
  Estimate" by `summary()` for a bgms comparison fitted with
  `group_indicator`, held the estimate of group 1. It now holds the bgms
  baseline, which is the mean of the group estimates.
* **Changes returned values.** Strength centrality of bgms fits with mixed
  variable types is now computed with each posterior draw placed on the edge
  it belongs to. It was placed by position, which permuted the edges.
* Under a stochastic block prior, the prior inclusion probability that
  `plot_prior_sensitivity()` places on its horizontal axis (`edge.prior`) is
  now taken from `bgms::extract_prior_inclusion_probabilities()`. It was
  computed from the estimated posterior partition. With bgms older than
  0.2.0.0 the previous calculation is kept.
* `summary()` of a bgms fit with mixed variable types now shows each edge its
  own R-hat. The values were placed by position and could belong to another
  edge.
* The Monte Carlo interval of each inclusion Bayes factor (`MCSE_BF`) of a
  bgms fit with mixed variable types now uses that edge's own Monte Carlo
  error. The errors were placed by position, so an interval could combine one
  edge's Bayes factor with another edge's error, and its row label could name
  another edge.
* `summary()` of a bgms comparison now shows each edge its own R-hat. bgms
  reports these in a different edge order than the summary table, so some
  edges (e.g. B-C and A-D) showed each other's R-hat.
* `plot_network()` on a raw bgms comparison of more than two groups now warns
  that the edge weights it shows are contrast coefficients rather than
  pairwise group differences, and may be drawn on the wrong edges.

## Other changes

* The objects returned by `easybgm()` and `easybgm_compare()` for bgms now
  include `packagefit`, the underlying bgms fit, so bgms extractors can be
  called on it without refitting.

# easybgm 0.5.0

## Support for bgms >= 0.2.0.0

* easybgm now works with the S7 fit objects returned by bgms >= 0.2.0.0. This
  switches off the temporary backward-compatibility shim that bgms shipped for
  easybgm 0.4.0, so fits are no longer converted to S3 lists and the per-fit
  compatibility warning is gone.
* `bgms` is now the default fitting package for every data type, including
  `type = "continuous"` and `type = "mixed"`, which previously defaulted to
  BGGM. Pass `package = "BGGM"` or `package = "BDgraph"` to keep the old
  behaviour.
* `type` may now be given as a per-variable character vector, for example
  `type = c("ordinal", "ordinal", "continuous")`.
* Priors may now be given as bgms prior objects (`cauchy_prior()`,
  `normal_prior()`, `beta_prime_prior()`, `bernoulli_prior()`,
  `beta_bernoulli_prior()`, `sbm_prior()`). The older flat arguments still work
  and are translated internally.

## Blume-Capel main effects

* Fits with at least one Blume-Capel variable now return
  `blume_capel_parameters`, a data frame holding the posterior mean, posterior
  standard deviation, 95% credible interval and R-hat of the linear and
  quadratic effect of each Blume-Capel variable, together with the baseline
  category it was fitted with. Unlike the category thresholds of an ordinal
  variable, these two parameters are usually of substantive interest, so they
  are also printed by `summary()` rather than left in the fit object.
* With `save = TRUE`, the posterior draws of those effects are returned in
  `samples_blume_capel`.
* Baseline categories are reported on the scale of the input data. bgms recodes
  discrete scores to start at 0 and shifts the baseline category with them, so
  the value it stores internally can be lower than the one the user supplied.

## Bug fixes

* For Blume-Capel variables, the two columns of `thresholds` were labelled
  `cat (1)` and `cat (2)`, the same headers bgms uses for genuine category
  thresholds. They are in fact the linear and quadratic effect, and are now
  named accordingly. Where Blume-Capel and ordinal variables share one matrix
  the headers cannot describe both, so the per-row meaning is recorded in the
  matrix's `"variable_type"` attribute.
* `print()` on an unsummarised `easybgm` object printed the closing notes twice,
  once from the summary it prints internally and once from its own tail.

* `centrality` for bgms fits was computed from a mis-permuted edge matrix: the
  posterior samples were read back in BGGM's upper-triangle order rather than
  the lower-triangle order bgms uses. Per-node strengths were therefore
  permuted, and `plot_centrality()` reported them under the wrong node labels.
  The ordering is now stated explicitly at each call site.
* `structure` was returned as a complete graph (a matrix of ones) whenever
  `save = FALSE`, which is the default. `plot_structure()` consequently drew a
  fully connected network. It is now the median probability model in both
  branches, matching the documentation.
* The Monte Carlo interval in `MCSE_BF` mixed two estimators: it took the
  binomial variance of the raw indicator average but divided it by the
  effective sample size of the Rao-Blackwellized chain, and attached the result
  to a Rao-Blackwellized Bayes factor. The interval was too wide by up to about
  40%. It is now computed from the Monte Carlo standard error that bgms reports
  for the Rao-Blackwellized inclusion probability.
* `plot_centrality()` and `plot_prior_sensitivity()` failed on lists of raw
  bgms fit objects, because bgms no longer reports `save` among the fit
  arguments. Both now work, and both record the model type correctly.
* `clusterBayesfactor()` failed on a raw bgms fit object with "invalid to use
  names()<- on an S4 object". It now reads the prior and the block posterior
  through the bgms extractor functions, and gives an informative error when the
  fit was not estimated with the Stochastic Block Model prior.
* The legacy `interaction_scale` argument no longer leaks a bgms deprecation
  warning; it is translated to `cauchy_prior(scale)` like `pairwise_scale`.
* Corrected the documented defaults for `interaction_prior` and
  `precision_scale_prior`, and documented `precision_graph_prior` and
  `difference_family`.

## Results that change

Fitting with bgms >= 0.2.0.0 changes several numbers relative to easybgm 0.4.0
with bgms 0.1.6.3. None of these is a bug in either package:

* **Edge weights are about half their former size.** bgms now reports pairwise
  parameters on the association scale, the coefficient entering each
  conditional as `2 * omega * x`, where 0.1.6.3 stored `2 * omega`. This
  affects `parameters`, `samples_posterior`, `centrality`, and every plot drawn
  from them.
* **Edge weights are not on the same scale across fitting packages.** BGGM and
  BDgraph report partial correlations; bgms reports the pairwise association
  parameter. For bgms fits of continuous and mixed data the partial
  correlations and the precision matrix are returned separately, in
  `partial_correlations` and `precision_matrix`, and `summary()` now states
  which scale the reported edge weights are on.
* **Inclusion probabilities and Bayes factors are Rao-Blackwellized.** They no
  longer saturate at 0 or 1 on short chains, so inclusion Bayes factors are
  finite where they used to be `0` or `Inf`, and an edge can cross the
  median-probability threshold differently than before.
* **The default interaction prior changed** from a Cauchy to
  `normal_prior(scale = 1)`, on the new coordinate. A Normal slab has much
  lighter tails than a Cauchy and constrains weakly identified edges more
  tightly.
* **`convergence_parameter` is the classic split-R-hat.** bgms 0.1.6.3 applied
  a degrees-of-freedom adjustment that reported about 1.29 on nearly saturated
  indicators, that is, on the most decisive edges. Those now report near 1. `NA`
  and `Inf` are possible when all chains are identical or stuck.
* **Group comparisons keep every category any group observes.** bgms 0.1.6.3
  merged categories that were not observed in every group, which biased the
  affected variable's pairwise parameters. Results move most where groups have
  unequal category support.

## Other changes

* easybgm continues to support bgms 0.1.6.3, as stated in `DESCRIPTION`. The
  fixes above apply to both bgms versions, and the test suite now exercises
  them on 0.1.6.3 as well as on 0.2.0.0 rather than skipping them.
* The examples and tests now pass `warmup` explicitly. bgms defaults to
  `warmup = 2000` regardless of `iter`, which dominated the runtime of the
  examples. `warmup` is a `bgm()` argument in both supported bgms versions.
* Removed the unused `vdiffr` dependency and the `LazyData` field, and dropped
  some dead version-gating code.
