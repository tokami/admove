# Extract data-storage tags

`get_dtags()` extracts data-storage (archival) tags from an object that
directly contains tagging data or from an object with a `tags` element.

## Usage

``` r
get_dtags(x)
```

## Arguments

- x:

  An `admove_tags` object, a data frame with a `tag_type` column, or an
  object containing a `tags` element.

## Value

A subset of the input tagging data containing only data-storage tags
(`tag_type == "d"`).
