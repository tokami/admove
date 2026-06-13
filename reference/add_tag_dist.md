# Compute predicted location distributions for a single tag

Propagates the model's predicted spatial location distribution for one
archival tag from its release to recovery time, using either the CTMC
forward-pass (engine 2) or repeated Kalman-filter track simulations
(engine 1). The result is stored in `fit$tag_dist` and consumed by
[`plot_tag_dist()`](https://tokami.github.io/admove/reference/plot_tag_dist.md).
Separating the expensive computation from rendering means plot
aesthetics can be changed without re-running the model.

## Usage

``` r
add_tag_dist(fit, i = 1, dt = 0.5, engine = NULL, xrel0 = NULL, yrel0 = NULL)
```

## Arguments

- fit:

  A fitted object of class `admove`, as returned by
  [`admove()`](https://tokami.github.io/admove/reference/admove.md).

- i:

  Integer index of the tag to use. Default is `1`.

- dt:

  Time step used when simulating tracks (engine 1 only). Default is
  `0.5`.

- engine:

  Optional integer overriding the engine stored in `fit$conf$engine`.
  `1` = Kalman filter, `2` = CTMC.

- xrel0, yrel0:

  Optional coordinates overriding the release location for CTMC-based
  predictions.

## Value

A copy of `fit` with the additional component `$tag_dist`, a list
containing the precomputed distributions and metadata required by
[`plot_tag_dist()`](https://tokami.github.io/admove/reference/plot_tag_dist.md).

## See also

[`plot_tag_dist()`](https://tokami.github.io/admove/reference/plot_tag_dist.md)
