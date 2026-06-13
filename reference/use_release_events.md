# Group tags by release events

`use_release_events()` groups tags that were released close in space and
time into common release events. This is mainly intended for
conventional mark-recapture tags, where multiple tags may share the same
approximate release location and release time.

The grouping is defined by the supplied spatial grid and time vector.

## Usage

``` r
use_release_events(x, grid, time_cont, tag_types = "c")
```

## Arguments

- x:

  Tagging data as an object of class `admove_tags` or an `admove_data`
  object containing tags.

- grid:

  Spatial grid used to aggregate release locations into common release
  events.

- time_cont:

  Numeric time vector used to aggregate release times into common
  release events.

- tag_types:

  Character vector giving the tag types to group by release event.
  Default: `"c"`.

## Value

An object like the input tagging data, but with selected tags grouped by
common release events.

## Details

For each selected tag type, the first observation of each tag is treated
as the release event. Releases falling into the same spatial grid cell
and the same time interval are grouped together and assigned a common
release event. The remaining observations are then combined under the
new grouped tag id.

## Examples

``` r
## use_release_events(tags, grid = dat$grid, time_cont = dat$time_cont)
```
