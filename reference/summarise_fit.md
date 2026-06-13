# Summarise a fitted admove model

Summarises the main results of a fitted `admove` model, including
parameter estimates and, where available, associated uncertainty
measures.

## Usage

``` r
summarise_fit(object, CI = 0.95, ...)

# S3 method for class 'admove'
summary(object, ...)
```

## Arguments

- object:

  A fitted model object of class `"admove"`, as returned by
  [`admove()`](https://tokami.github.io/admove/reference/admove.md).

- CI:

  Numeric scalar giving the confidence level used for confidence
  intervals. Default is `0.95`.

- ...:

  Additional arguments passed to internal summary methods.

## Value

A summary object, typically printed for inspection.
