
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
##' @seealso [sf::st_crs()]
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
##' @name crs-set
##' @export
`crs<-` <- function(x, value) UseMethod("crs<-")

##' @rdname crs-set
##' @export
`crs<-.default` <- function(x, value) {
  sp <- sref(x)
  sp$crs <- value
  sref(x) <- sp
  x
}

##' @rdname crs-set
##' @export
`crs<-.admove_sref` <- function(x, value) {
  x$crs <- value
  x
}

##' @rdname crs-set
##' @export
`crs<-.admove_grid` <- function(x, value) {
  sp <- sref(x)
  sp$crs <- value
  sref(x) <- sp
  x
}

##' @rdname crs-set
##' @export
`crs<-.admove_cov` <- function(x, value) {
  sp <- sref(x)
  sp$crs <- value
  sref(x) <- sp
  x
}

##' @rdname crs-set
##' @export
`crs<-.admove_tags` <- function(x, value) {
  sp <- sref(x)
  sp$crs <- value
  sref(x) <- sp
  x
}

##' @rdname crs-set
##' @export
`crs<-.admove_data` <- function(x, value) {
  sp <- sref(x)
  sp$crs <- value
  sref(x) <- sp
  x
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
##' @return `x` with updated spatial units.
##'
##' @name units_space-set
##' @export
`units_space<-` <- function(x, value) UseMethod("units_space<-")

##' @rdname units_space-set
##' @export
`units_space<-.default` <- function(x, value) {
  sp <- sref(x)
  sp$units <- value
  sref(x) <- sp
  x
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_sref` <- function(x, value) {
  x$units <- value
  x
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_grid` <- function(x, value) {
  sp <- sref(x)
  sp$units <- value
  sref(x) <- sp
  x
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_cov` <- function(x, value) {
  sp <- sref(x)
  sp$units <- value
  sref(x) <- sp
  x
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_tags` <- function(x, value) {
  sp <- sref(x)
  sp$units <- value
  sref(x) <- sp
  x
}

##' @rdname units_space-set
##' @export
`units_space<-.admove_data` <- function(x, value) {
  sp <- sref(x)
  sp$units <- value
  sref(x) <- sp
  x
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
##' @return `x` with updated CRS scaling information.
##'
##' @name crs_scale-set
##' @export
`crs_scale<-` <- function(x, value) UseMethod("crs_scale<-")

##' @rdname crs_scale-set
##' @export
`crs_scale<-.default` <- function(x, value) {
  sp <- sref(x)
  sp$crs_scale <- value
  sref(x) <- sp
  x
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_sref` <- function(x, value) {
  x$crs_scale <- value
  x
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_grid` <- function(x, value) {
  sp <- sref(x)
  sp$crs_scale <- value
  sref(x) <- sp
  x
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_cov` <- function(x, value) {
  sp <- sref(x)
  sp$crs_scale <- value
  sref(x) <- sp
  x
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_tags` <- function(x, value) {
  sp <- sref(x)
  sp$crs_scale <- value
  sref(x) <- sp
  x
}

##' @rdname crs_scale-set
##' @export
`crs_scale<-.admove_data` <- function(x, value) {
  sp <- sref(x)
  sp$crs_scale <- value
  sref(x) <- sp
  x
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
  tr <- tref(x)
  tr$period <- value
  tref(x) <- tr
  x
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
  tr <- sref(x)
  tr$period <- value
  tref(x) <- tr
  x
}

##' @rdname period-set
##' @export
`period<-.admove_tags` <- function(x, value) {
  tr <- sref(x)
  tr$period <- value
  tref(x) <- tr
  x
}

##' @rdname period-set
##' @export
`period<-.admove_data` <- function(x, value) {
  tr <- sref(x)
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
