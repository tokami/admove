# Combine `admove_tags` objects

`combine_tags()` combines multiple objects of class `admove_tags` into a
single `admove_tags` object. It is designed to support workflows such as
combining data-storage, mark-resight, and mark-recapture tags into one
unified data set.

All inputs must inherit from `admove_tags` and must have compatible
spatial and temporal reference information.

## Usage

``` r
combine_tags(..., recursive = FALSE)

# S3 method for class 'admove_tags'
c(..., recursive = FALSE)
```

## Arguments

- ...:

  Objects of class `admove_tags`, or a single list containing such
  objects.

- recursive:

  Ignored. Included for compatibility with the generic
  [`base::c()`](https://rdrr.io/r/base/c.html) interface.

## Value

A single object of class `admove_tags` containing all rows from the
input objects.

## Details

Before combining inputs, the function:

- removes `NULL` inputs,

- checks that all remaining inputs inherit from `admove_tags`,

- aligns columns across objects by filling missing columns with `NA`,

- checks that all objects have the same `sref`, and

- checks that all objects have the same `tref`.

If spatial or temporal reference information differs across inputs, the
function stops with an error.

## Examples

``` r
## tags <- combine_tags(dtags, stags, ctags)
## tags <- c(dtags, ctags)
```
