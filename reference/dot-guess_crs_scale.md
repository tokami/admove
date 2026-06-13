# Guess a coordinate scaling factor from a spatial reference (sref)

Compute a numeric factor that converts coordinates stored under a given
`admove_sref` to a requested target unit. This extends the previous
CRS-only logic by accounting for `sref$crs_scale`, i.e. the map from CRS
units to stored units.

## Usage

``` r
.guess_crs_scale(sref, units)
```

## Arguments

- sref:

  An `admove_sref` object (or a list with components `crs`, `units`,
  `crs_scale`).

- units:

  Character string describing the desired target units, e.g. `"m"`,
  `"km"`, `"mi"`, `"nmi"`.

## Value

A numeric scalar scaling factor, or `NA_real_` if it cannot be
determined.

## Details

The returned sref is a multiplicative factor `f` such that:
\$\$x\_{\mathrm{target}} = x\_{\mathrm{stored}} \times f\$\$

If CRS units are known (via sf) and both the stored units and target
units can be converted to meters, then: \$\$ f = \frac{\mathrm{meters\\
per\\ stored\\ unit}}{\mathrm{meters\\ per\\ target\\ unit}} \$\$

where meters per stored unit is derived from the CRS units and
`sref$crs_scale`.
