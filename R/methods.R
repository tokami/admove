
## Time reference information -------------------------------------------------------

##' Get spatial reference from an object
##'
##' @description
##' Generic function to extract spatial reference information from an object.
##' \pkg{admove} objects typically store this information in
##' `attr(x, "sref")`.
##'
##' @param x An object.
##' @param ... Further arguments passed to methods (currently unused).
##'
##' @return An object of class `admove_sref`.
##'
##' @name sref
##' @export
sref <- function(x, ...) UseMethod("sref")

##' @rdname sref
##' @export
sref.default <- function(x, ...) {
  if (is.null(x)) return(NULL)
  sp <- attr(x, "sref")
  if (is.null(sp)) stop("Object has no 'sref' attribute.")
  validate_sref(sp)
}

##' Set spatial reference on an object
##'
##' @description
##' Replacement function for [sref()].
##'
##' @param x An object to modify.
##' @param value An object of class `admove_sref`.
##'
##' @return `x` with an updated `"sref"` attribute.
##'
##' @name sref-set
##' @export
`sref<-` <- function(x, value) UseMethod("sref<-")

##' @rdname sref-set
##' @export
`sref<-.default` <- function(x, value) {
  sref <- validate_sref(value)
  attr(x, "sref") <- sref
  x
}




## CRS -------------------------------------------------------------------------------

##' Get coordinate reference system (CRS)
##'
##' @description
##' Generic function to extract the coordinate reference system (CRS) from a
##' supported object.
##'
##' @param x An object from which to extract a CRS.
##' @param ... Further arguments passed to methods.
##'
##' @return
##' A CRS representation. For \pkg{admove} objects, this is typically whatever
##' is stored in the corresponding spatial reference object, for example a WKT
##' string, EPSG code, or another \pkg{sf}-compatible CRS specification.
##'
##' @details
##' `crs()` returns the CRS of the *unscaled* coordinates and deliberately does
##' **not** reflect [units_space()]. Stored coordinates are CRS coordinates times
##' [crs_scale()], so the two are complementary rather than redundant: a
##' metre-based projection whose coordinates are stored in km keeps a metre CRS
##' and `crs_scale = 0.001`. Everything that hands the CRS to \pkg{sf} divides
##' the stored coordinates by `crs_scale` first, so a CRS rewritten to the stored
##' unit would convert twice.
##'
##' ```r
##' crs_aeqd <- paste("+proj=aeqd +lat_0=40 +lon_0=-6 +x_0=0 +y_0=0",
##'                   "+datum=WGS84 +units=m +no_defs")
##' sp <- create_sref(crs = crs_aeqd, units = "km")
##' sf::st_crs(crs(sp))$units_gdal   ## "metre"  -- the projection
##' units_space(sp)                  ## "km"     -- the stored numbers
##' crs_scale(sp)                    ## 0.001    -- links the two
##' ```
##'
##' @seealso [sf::st_crs()], [units_space()], [crs_scale()], [scale_sref()],
##'   [create_sref()]
##'
##' @name crs
##' @export
crs <- function(x, ...) UseMethod("crs")

##' @rdname crs
##' @export
crs.default <- function(x, ...) {
  sref(x)$crs
}

##' @rdname crs
##' @export
crs.admove_sref <- function(x, ...) {
  x$crs
}

##' @rdname crs
##' @export
crs.admove_grid <- function(x, ...) sref(x)$crs

##' @rdname crs
##' @export
crs.admove_cov <- function(x, ...) sref(x)$crs

##' @rdname crs
##' @export
crs.admove_tags <- function(x, ...) sref(x)$crs

##' @rdname crs
##' @export
crs.admove_data <- function(x, ...) sref(x)$crs


##' Set coordinate reference system (CRS)
##'
##' @description
##' Replacement function for [crs()].
##'
##' @param x An object to modify.
##' @param value A CRS specification.
##'
##' @return `x` with updated CRS information.
##'
##' @details
##' Only defined for [create_sref()] objects, where the other fields are
##' recomputed so the spatial reference stays self-consistent. On objects that
##' carry coordinates (`admove_grid`, `admove_cov`, `admove_tags`,
##' `admove_data`) this errors, because changing the metadata alone would leave
##' the coordinates stale: use [add_sref()], [scale_sref()] or
##' [transform_sref()], which move the coordinates too.
##'
##' @seealso [transform_sref()], [add_sref()]
##'
##' @name crs-set
##' @export
`crs<-` <- function(x, value) UseMethod("crs<-")

##' @rdname crs-set
##' @export
`crs<-.default` <- function(x, value) {
  .stop_sref_setter("crs")
}

##' @rdname crs-set
##' @export
`crs<-.admove_sref` <- function(x, value) {

  ## keep the stored-unit label only while it still means the same thing, i.e.
  ## while the new CRS has the same native unit as the old one; otherwise let
  ## create_sref() re-infer it, as add_sref() does on a CRS change
  units <- NA_character_
  if (!.is_na_scalar(x$units) &&
        identical(.infer_sref_units(x$crs, verbose = FALSE),
                  .infer_sref_units(value, verbose = FALSE))) {
    units <- x$units
  }

  create_sref(crs = value, units = units)
}

##' @rdname crs-set
##' @export
`crs<-.admove_grid` <- function(x, value) {
  .stop_sref_setter("crs")
}

##' @rdname crs-set
##' @export
`crs<-.admove_cov` <- function(x, value) {
  .stop_sref_setter("crs")
}

##' @rdname crs-set
##' @export
`crs<-.admove_tags` <- function(x, value) {
  .stop_sref_setter("crs")
}

##' @rdname crs-set
##' @export
`crs<-.admove_data` <- function(x, value) {
  .stop_sref_setter("crs")
}



## units_space -----------------------------------------------------------------------

##' Get spatial units
##'
##' @description
##' Generic function to extract spatial units from a supported object.
##'
##' @param x An object from which to extract spatial units.
##' @param ... Further arguments passed to methods.
##'
##' @return A character string describing the spatial units, such as `"m"`,
##'   `"km"`, or `"degree"`.
##'
##' @name units_space
##' @export
units_space <- function(x, ...) UseMethod("units_space")

##' @rdname units_space
##' @export
units_space.default <- function(x, ...) {
  sref(x)$units
}

##' @rdname units_space
##' @export
units_space.admove_sref <- function(x, ...) {
  x$units
}

##' @rdname units_space
##' @export
units_space.admove_grid <- function(x, ...) sref(x)$units

##' @rdname units_space
##' @export
units_space.admove_cov <- function(x, ...) sref(x)$units

##' @rdname units_space
##' @export
units_space.admove_tags <- function(x, ...) sref(x)$units

##' @rdname units_space
##' @export
units_space.admove_data <- function(x, ...) sref(x)$units



##' Set spatial units
##'
##' @description
##' Replacement function for [units_space()].
##'
##' @param x An object to modify.
##' @param value A character string describing the spatial units.
##'
##' @return `x` with updated spatial units and a matching [crs_scale()].
##'
##' @details
##' Only defined for [create_sref()] objects, where the other fields are
##' recomputed so the spatial reference stays self-consistent. On objects that
##' carry coordinates (`admove_grid`, `admove_cov`, `admove_tags`,
##' `admove_data`) this errors, because changing the metadata alone would leave
##' the coordinates stale: use [add_sref()], [scale_sref()] or
##' [transform_sref()], which move the coordinates too.
##'
##' @seealso [scale_sref()], [add_sref()]
##'
##' @name units_space-set
##' @export
`units_space<-` <- function(x, value) UseMethod("units_space<-")

##' @rdname units_space-set
##' @export
`units_space<-.default` <- function(x, value) {
  .stop_sref_setter("units_space", value)
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_sref` <- function(x, value) {
  create_sref(crs = x$crs, units = value)
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_grid` <- function(x, value) {
  .stop_sref_setter("units_space", value)
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_cov` <- function(x, value) {
  .stop_sref_setter("units_space", value)
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_tags` <- function(x, value) {
  .stop_sref_setter("units_space", value)
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_data` <- function(x, value) {
  .stop_sref_setter("units_space", value)
}



## crs_scale -----------------------------------------------------------------------

##' Get CRS scale
##'
##' @description
##' Generic function to extract the scaling factor relating CRS units to the
##' working spatial units used by an object.
##'
##' @param x An object.
##' @param ... Further arguments passed to methods.
##'
##' @return A numeric scaling factor.
##'
##' @name crs_scale
##' @export
crs_scale <- function(x, ...) UseMethod("crs_scale")

##' @rdname crs_scale
##' @export
crs_scale.default <- function(x, ...) {
  sref(x)$crs_scale
}

##' @rdname crs_scale
##' @export
crs_scale.admove_sref <- function(x, ...) {
  x$crs_scale
}

##' @rdname crs_scale
##' @export
crs_scale.admove_grid <- function(x, ...) sref(x)$crs_scale

##' @rdname crs_scale
##' @export
crs_scale.admove_cov <- function(x, ...) sref(x)$crs_scale

##' @rdname crs_scale
##' @export
crs_scale.admove_tags <- function(x, ...) sref(x)$crs_scale

##' @rdname crs_scale
##' @export
crs_scale.admove_data <- function(x, ...) sref(x)$crs_scale



##' Set CRS scale
##'
##' @description
##' Replacement function for [crs_scale()].
##'
##' @param x An object to modify.
##' @param value A numeric scaling factor.
##'
##' @return `x` with updated CRS scaling information and a matching
##'   [units_space()] label.
##'
##' @details
##' Only defined for [create_sref()] objects, where the other fields are
##' recomputed so the spatial reference stays self-consistent. On objects that
##' carry coordinates (`admove_grid`, `admove_cov`, `admove_tags`,
##' `admove_data`) this errors, because changing the metadata alone would leave
##' the coordinates stale: use [add_sref()], [scale_sref()] or
##' [transform_sref()], which move the coordinates too.
##'
##' @seealso [scale_sref()], [add_sref()]
##'
##' @name crs_scale-set
##' @export
`crs_scale<-` <- function(x, value) UseMethod("crs_scale<-")

##' @rdname crs_scale-set
##' @export
`crs_scale<-.default` <- function(x, value) {
  .stop_sref_setter("crs_scale", value)
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_sref` <- function(x, value) {
  create_sref(crs = x$crs, crs_scale = value)
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_grid` <- function(x, value) {
  .stop_sref_setter("crs_scale", value)
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_cov` <- function(x, value) {
  .stop_sref_setter("crs_scale", value)
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_tags` <- function(x, value) {
  .stop_sref_setter("crs_scale", value)
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_data` <- function(x, value) {
  .stop_sref_setter("crs_scale", value)
}





## Time reference information -------------------------------------------------------


##' Get time reference from an object
##'
##' @description
##' Generic function to extract time reference information from an object.
##' \pkg{admove} objects typically store this information in
##' `attr(x, "tref")`.
##'
##' @param x An object.
##' @param ... Further arguments passed to methods (currently unused).
##'
##' @return An object of class `admove_tref`.
##'
##' @name tref
##' @export
tref <- function(x, ...) UseMethod("tref")

##' @rdname tref
##' @export
tref.default <- function(x, ...) {
  if (is.null(x)) return(NULL)
  tr <- attr(x, "tref")
  if (is.null(tr)) stop("Object has no 'tref' attribute.")
  tr
}


##' Set time reference on an object
##'
##' @description
##' Replacement function for [tref()].
##'
##' @param x An object to modify.
##' @param value An object of class `admove_tref`.
##'
##' @return `x` with an updated `"tref"` attribute.
##'
##' @name tref-set
##' @export
`tref<-` <- function(x, value) UseMethod("tref<-")

##' @rdname tref-set
##' @export
`tref<-.default` <- function(x, value) {
  attr(x, "tref") <- value
  x
}


## origin -----------------------------------------------------------------------------

##' Get temporal origin
##'
##' @description
##' Generic function to extract the temporal origin from a supported object.
##'
##' @param x An object from which to extract the temporal origin.
##' @param ... Further arguments passed to methods.
##'
##' @return A `POSIXct` time origin, or another stored origin representation.
##'
##' @name origin
##' @export
origin <- function(x, ...) UseMethod("origin")

##' @rdname origin
##' @export
origin.default <- function(x, ...) {
  tref(x)$origin
}

##' @rdname origin
##' @export
origin.admove_tref <- function(x, ...) {
  x$origin
}

##' @rdname origin
##' @export
origin.admove_cov <- function(x, ...) tref(x)$origin

##' @rdname origin
##' @export
origin.admove_tags <- function(x, ...) tref(x)$origin

##' @rdname origin
##' @export
origin.admove_data <- function(x, ...) tref(x)$origin




##' Set temporal origin
##'
##' @description
##' Replacement function for [origin()].
##'
##' @param x An object to modify.
##' @param value A temporal origin, typically of class `POSIXct`.
##'
##' @return `x` with updated temporal origin information.
##'
##' @name origin-set
##' @export
`origin<-` <- function(x, value) UseMethod("origin<-")

##' @rdname origin-set
##' @export
`origin<-.default` <- function(x, value) {
  tr <- tref(x)
  tr$origin <- value
  tref(x) <- tr
  x
}

##' @rdname origin-set
##' @export
`origin<-.admove_tref` <- function(x, value) {
  x$origin <- value
  x
}

##' @rdname origin-set
##' @export
`origin<-.admove_cov` <- function(x, value) {
  tr <- sref(x)
  tr$origin <- value
  tref(x) <- tr
  x
}

##' @rdname origin-set
##' @export
`origin<-.admove_tags` <- function(x, value) {
  tr <- sref(x)
  tr$origin <- value
  tref(x) <- tr
  x
}

##' @rdname origin-set
##' @export
`origin<-.admove_data` <- function(x, value) {
  tr <- sref(x)
  tr$origin <- value
  tref(x) <- tr
  x
}


## units_time -------------------------------------------------------------------------

##' Get temporal units
##'
##' @description
##' Generic function to extract temporal units from a supported object.
##'
##' @param x An object from which to extract temporal units.
##' @param ... Further arguments passed to methods.
##'
##' @return A character string describing the temporal units, such as `"day"`,
##'   `"month"`, `"year"`, or `"week"`.
##'
##' @name units_time
##' @export
units_time <- function(x, ...) UseMethod("units_time")

##' @rdname units_time
##' @export
units_time.default <- function(x, ...) {
  tref(x)$units
}

##' @rdname units_time
##' @export
units_time.admove_tref <- function(x, ...) {
  x$units
}

##' @rdname units_time
##' @export
units_time.admove_cov <- function(x, ...) tref(x)$units

##' @rdname units_time
##' @export
units_time.admove_tags <- function(x, ...) tref(x)$units

##' @rdname units_time
##' @export
units_time.admove_data <- function(x, ...) tref(x)$units


##' Set temporal units
##'
##' @description
##' Replacement function for [units_time()].
##'
##' @param x An object to modify.
##' @param value A character string describing the temporal units.
##'
##' @return `x` with updated temporal units.
##'
##' @name units_time-set
##' @export
`units_time<-` <- function(x, value) UseMethod("units_time<-")

##' @rdname units_time-set
##' @export
`units_time<-.default` <- function(x, value) {
  tr <- tref(x)
  tr$units <- value
  tref(x) <- tr
  x
}

##' @rdname units_time-set
##' @export
`units_time<-.admove_tref` <- function(x, value) {
  x$units <- value
  x
}

##' @rdname units_time-set
##' @export
`units_time<-.admove_cov` <- function(x, value) {
  tr <- sref(x)
  tr$units <- value
  tref(x) <- tr
  x
}

##' @rdname units_time-set
##' @export
`units_time<-.admove_tags` <- function(x, value) {
  tr <- sref(x)
  tr$units <- value
  tref(x) <- tr
  x
}

##' @rdname units_time-set
##' @export
`units_time<-.admove_data` <- function(x, value) {
  tr <- sref(x)
  tr$units <- value
  tref(x) <- tr
  x
}



## period -----------------------------------------------------------------------------

##' Get temporal period
##'
##' @description
##' Generic function to extract the temporal period from a supported object.
##'
##' @param x An object from which to extract the period.
##' @param ... Further arguments passed to methods.
##'
##' @return A numeric value giving the number of time steps per year, or another
##'   stored period representation.
##'
##' @name period
##' @export
period <- function(x, ...) UseMethod("period")

##' @rdname period
##' @export
period.default <- function(x, ...) {
  tref(x)$period
}

##' @rdname period
##' @export
period.admove_tref <- function(x, ...) {
  x$period
}

##' @rdname period
##' @export
period.admove_cov <- function(x, ...) tref(x)$period

##' @rdname period
##' @export
period.admove_tags <- function(x, ...) tref(x)$period

##' @rdname period
##' @export
period.admove_data <- function(x, ...) tref(x)$period



##' Set temporal period
##'
##' @description
##' Replacement function for [period()].
##'
##' @param x An object to modify.
##' @param value A period value.
##'
##' @return `x` with updated period information.
##'
##' @name period-set
##' @export
`period<-` <- function(x, value) UseMethod("period<-")


##' @rdname period-set
##' @export
`period<-.default` <- function(x, value) {
  .set_period(x, value)
}


##' @rdname period-set
##' @export
`period<-.admove_tref` <- function(x, value) {
  x$period <- value
  x
}

##' @rdname period-set
##' @export
`period<-.admove_cov` <- function(x, value) {
  .set_period(x, value)
}

##' @rdname period-set
##' @export
`period<-.admove_tags` <- function(x, value) {
  .set_period(x, value)
}

##' @rdname period-set
##' @export
`period<-.admove_data` <- function(x, value) {
  .set_period(x, value)
}


## Set the seasonal period on the time reference of an object. The period is
## part of the time reference, so an object without one has nowhere to store it.
.set_period <- function(x, value) {
  tr <- tref(x)
  if (is.null(tr)) {
    stop("This object has no time reference, so a period cannot be set on it. ",
         "Add one first, e.g. x <- add_tref(x, create_tref(origin = ..., ",
         "units = ..., period = ", deparse(value), ")).", call. = FALSE)
  }
  tr$period <- value
  tref(x) <- tr
  x
}



## Other ----------------------------------------------------------------------------


##' Print admove objects
##'
##' @description
##' Print methods for `admove` classes.
##'
##' @param x An object to print.
##' @param ... Further arguments passed to or from other methods.
##'
##' @return
##' The object is printed to the console. Invisibly returns `x`.
##'
##' @name print-admove
NULL
