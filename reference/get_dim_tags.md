# Get spatial and temporal ranges from tagging data

`get_dim_tags()` extracts the observed temporal range and the spatial
ranges in the `x` and `y` directions from tagging data.

## Usage

``` r
get_dim_tags(tags = NULL)
```

## Arguments

- tags:

  Tagging data as a data frame, an object of class `admove_tags`, or a
  list of tag-specific data frames.

## Value

A list with components:

- trange:

  Range of observed times.

- xrange:

  Range of observed `x` coordinates.

- yrange:

  Range of observed `y` coordinates.
