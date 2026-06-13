# Subset an `admove_cov` object

Subsetting method for objects of class `admove_cov`.

This method preserves the `admove_cov` class and associated attributes
when the result remains a three-dimensional covariate array. If
subsetting returns an object with fewer than three dimensions, the
result is returned as a regular R object without the `admove_cov` class.

In particular, subsetting only along the third dimension (time) returns
a subsetted `admove_cov` object with the corresponding time-related
attributes updated where available.

## Usage

``` r
# S3 method for class 'admove_cov'
x[i, j, k, ..., drop = TRUE]
```

## Arguments

- x:

  An object of class `admove_cov`.

- i:

  Indices for the first dimension.

- j:

  Indices for the second dimension.

- k:

  Indices for the third dimension, typically corresponding to time.

- ...:

  Further indices passed to the underlying array subsetting operation.

- drop:

  Logical; should dimensions of length one be dropped? Default: `TRUE`.

## Value

An object of class `admove_cov` if the subset retains three dimensions;
otherwise a regular subsetted R object.

## Details

The method first performs standard array subsetting on the unclassed
object. If the resulting object still has three dimensions, attributes
other than `dim` and `dimnames` are restored and the original class is
reattached.
