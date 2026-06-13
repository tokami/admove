# Rescale covariate fields

Rescales a covariate field, or one element of a list of covariate
fields, to a specified numeric range.

## Usage

``` r
.rescale_cov(cov, zrange, i = NULL)
```

## Arguments

- cov:

  A covariate field or a list of covariate fields.

- zrange:

  Numeric vector of length 2 giving the target range.

- i:

  Optional index specifying which element of `cov` to rescale if `cov`
  is a list.

## Value

A rescaled covariate field with values mapped to `zrange`.
