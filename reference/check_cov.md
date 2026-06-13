# Check and standardise covariate data

`check_cov()` validates covariate data and ensures that it can be used
by *admove*. If the input is not already of class `admove_cov` or
`admove_cov_list`, the function first tries to convert it using
[`prep_cov()`](https://tokami.github.io/admove/reference/prep_cov.md).

Existing covariate objects are checked for missing values, missing
dimension names, and incompatible dimension-name lengths. Invalid
covariate fields are removed.

## Usage

``` r
check_cov(x, verbose = TRUE)
```

## Arguments

- x:

  Covariate data to check. Can be `NULL`, an object of class
  `admove_cov` or `admove_cov_list`, or another object that can be
  converted with
  [`prep_cov()`](https://tokami.github.io/admove/reference/prep_cov.md).

- verbose:

  Logical; if `TRUE`, print informative messages about removed covariate
  fields. Default: `TRUE`.

## Value

`NULL` if no valid covariate field remains; otherwise an object of class
`admove_cov` or `admove_cov_list`.

## Details

For existing covariate objects, the function removes covariate fields
that:

- contain only missing values,

- have missing or incomplete dimension names, or

- have dimension names whose lengths do not match the array dimensions.
