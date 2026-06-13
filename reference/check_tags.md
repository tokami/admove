# Check and standardise a tag object

`check_tags()` validates tagging data, removes invalid entries,
optionally checks whether observations fall within the spatial and
temporal domain of a grid or data object, and returns a cleaned object
of class `admove_tags`.

The function can be used on raw tagging data, on an `admove_tags`
object, or on higher-level *admove* objects that contain tags.

## Usage

``` r
check_tags(
  x,
  grid = NULL,
  dat = NULL,
  conf = NULL,
  remove_non_recovered_tags = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  Tagging data to check. Can be an object of class `admove_tags`, a data
  frame, a list of tag-specific data frames, or an object that contains
  tags such as `admove_data`, `admove_sim`, or `admove`.

- grid:

  Optional spatial grid used to remove observations outside the spatial
  domain and to assign grid cells.

- dat:

  Optional `admove_data` object used to check whether observations fall
  within the time domain.

- conf:

  Optional configuration list. If provided, only tag types enabled in
  `conf` are retained.

- remove_non_recovered_tags:

  Logical; if `TRUE`, remove tags with fewer than two observations.
  Default: `TRUE`.

- verbose:

  Logical; if `TRUE`, print informative messages about removed entries.
  Default: `TRUE`.

## Value

A cleaned object of class `admove_tags`.

## Details

The function performs several checks, including:

- removal of rows with missing required values in `t`, `x`, or `y`,

- removal of rows with `use = FALSE`,

- removal of observations outside the spatial domain of `grid`,

- removal of observations outside the temporal domain of `dat`, and

- optional removal of tags with fewer than two observations.

Spatial and temporal reference information are attached to the returned
object from `grid` and `dat` when available.
