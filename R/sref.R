
##' Spatial reference for admove objects
##'
##' A lightweight container that stores the spatial reference information used
##' throughout \pkg{admove}. It bundles:
##' \itemize{
##'   \item \strong{crs}: a CRS specification (typically WKT, PROJ string, or EPSG)
##'   \item \strong{units}: the units of the coordinates stored in the object
##'   \item \strong{crs_scale}: conversion factor from CRS-units to stored units
##' }
##'
##' The intention is that \emph{objects store their coordinates in the chosen
##' display/model units} (e.g. km), while \code{crs} still defines the underlying
##' standard CRS (typically meters). The \code{crs_scale} links the two.
##'
##' @param crs CRS specification. Can be WKT/PROJ string, EPSG integer, or an
##'   \code{sf::crs} input.
##' @param units Character string describing units of stored coordinates, e.g.
##'   \code{"m"}, \code{"km"}, \code{"degree"}.
##' @param crs_scale Numeric scalar conversion factor from CRS units to stored
##'   units. For example, if CRS is meters and stored coordinates are kilometers,
##'   \code{crs_scale = 0.001}.
##'
##' @return An object of class \code{admove_sref}.
##'
##' @examples
##' sp <- create_sref(crs = 32631, units = "km", crs_scale = 0.001)
##' sp
##'
##' @export
create_sref <- function(crs = NA, units = NA_character_, crs_scale = 1) {

  ## crs
  if (is.null(crs)) {
    crs <- NA
  }
  crs <- .normalize_crs_wkt(crs)
  if (requireNamespace("sf", quietly = TRUE)) {
    crs_sf <- try(sf::st_crs(crs), silent = TRUE)
    if (inherits(crs_sf, "try-error") || is.na(crs_sf)) {
      crs_sf <- NA
    }
  } else crs_sf <- NA

  ## units
  if (!is.null(crs_sf) && !is.na(crs_sf) && !is.null(crs_sf$units_gdal)) {
    units_crs <- crs_sf$units_gdal
  } else units_crs <- NA
  if (is.null(units) || is.na(units)) {
    units <- units_crs
  }

  if (!is.na(units) && !is.na(units_crs)) {
    crs_scale <- .in_m(units_crs) / .in_m(units)
  }

  ## crs_scale
  if (is.null(crs_scale) || is.na(crs_scale) || length(crs_scale) == 0) {
    if ((!is.null(crs) && !is.na(crs)) || (!is.na(units) && !is.null(units))) {
      crs_scale <- 1
    } else {
      crs_scale <- NA
    }
  }
  sp <- structure(
    list(crs = crs, units = units, crs_scale = crs_scale),
    class = "admove_sref"
  )
  validate_sref(sp)
}

##' Validate an admove spatial reference
##'
##' Internal helper used to ensure spatial reference objects are consistent.
##'
##' @param x An \code{admove_sref} object.
##'
##' @return The validated \code{admove_sref} (invisibly identical to input).
##'
##' @keywords internal
validate_sref <- function(x) {
  if (!inherits(x, "admove_sref")) {
    stop("'x' must inherit from class 'admove_sref'.")
  }
  ## if (!is.numeric(x$crs_scale) || length(x$crs_scale) != 1 ||
  ##     is.na(x$crs_scale) || !is.finite(x$crs_scale) || x$crs_scale <= 0) {
  ##   stop("'crs_scale' must be a single positive finite number.")
  ## }
  x
}


##' Add or harmonise a spatial reference on an object
##'
##' Attach a spatial reference (\code{sref}) to an object, or harmonize an
##' existing spatial reference with that of another object or sref-like
##' specification.
##'
##' The spatial reference bundles a coordinate reference system (CRS), the units
##' of the stored coordinates, and a unit scaling factor linking CRS units to the
##' stored coordinate units.
##'
##' This function is primarily intended to add or update \emph{metadata}. If the
##' object already has a spatial reference and the requested CRS is the same but
##' \code{crs_scale} differs, the stored coordinates are rescaled via
##' \code{\link{scale_sref}}. If the requested CRS differs, changing the CRS
##' metadata alone would change the meaning of the stored coordinates. Therefore,
##' by default the function throws an error in this case. To preserve the
##' represented locations while adopting the new CRS, set
##' \code{transform_crs = TRUE}; this calls \code{\link{transform_sref}}
##' internally.
##'
##' If a CRS is provided and the \pkg{sf} package is available, the CRS is
##' normalized using \code{\link[sf]{st_crs}} and stored as WKT for stability.
##' If \code{units} are not provided, they are inferred from the target CRS when
##' possible. If \code{crs_scale} is not provided, it is inferred with
##' \code{\link{.guess_crs_scale}} where possible.
##'
##' @param x An object to which spatial reference information should be attached
##'   or whose spatial reference should be harmonized.
##' @param sref Optional spatial reference specification. This can be:
##'   \itemize{
##'   \item an object of class \code{"admove_sref"},
##'   \item a named list containing one or more of \code{crs}, \code{units}, and
##'     \code{crs_scale},
##'   \item another object from which \code{sref(sref)} can be extracted.
##'   }
##'   If \code{NULL}, existing spatial reference information on \code{x} is kept
##'   where possible.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed when
##'   CRS, units, or scaling information are inferred, retained, or transformed.
##' @param transform_crs Logical; if \code{TRUE} and the current CRS of
##'   \code{x} differs from the requested CRS, stored coordinates are transformed
##'   with \code{\link{transform_sref}} so that represented locations remain the
##'   same. If \code{FALSE} (default), differing non-missing CRS values cause an
##'   error.
##'
##' @return \code{x} with updated \code{sref(x)}. Depending on the input and
##'   requested spatial reference, the stored coordinates may also be rescaled
##'   (same CRS, different \code{crs_scale}) or transformed (different CRS and
##'   \code{transform_crs = TRUE}).
##'
##' @details
##' The function distinguishes between three kinds of spatial-reference updates:
##' \enumerate{
##'   \item \strong{Adding missing spatial metadata}: if \code{x} has no spatial
##'     reference, or if parts of it are missing, missing fields are filled from
##'     \code{sref} where possible.
##'   \item \strong{Changing stored coordinate units}: if the CRS is unchanged but
##'     \code{crs_scale} differs, stored coordinates are rescaled using
##'     \code{\link{scale_sref}}.
##'   \item \strong{Changing CRS}: if CRS differ, stored coordinates must be
##'     transformed if the represented locations are to remain unchanged. This is
##'     not done silently; use \code{transform_crs = TRUE} or call
##'     \code{\link{transform_sref}} directly.
##' }
##'
##' @seealso \code{\link{transform_sref}}, \code{\link{scale_sref}},
##'   \code{\link{create_sref}}, \code{\link{sref}},
##'   \code{\link{.guess_crs_scale}}
##'
##' @examples
##' ## Not run:
##' ## Add a new spatial reference
##' grid <- create_grid(cellsize = 5e3, xrange = c(0, 5e4), yrange = c(0, 5e4))
##' grid <- add_sref(grid, sref = list(crs = 32631, units = "m"))
##'
##' ## Same CRS, but store coordinates in km rather than m
##' grid <- add_sref(grid, sref = list(crs = 32631, units = "km",
##'                                      crs_scale = 0.001))
##'
##' ## Adopt the CRS of another object while preserving locations
##' ## tags <- add_sref(tags, sref = other_tags, transform_crs = TRUE)
##' ## End(Not run)
##'
##' @export
add_sref <- function(x, sref = NULL, verbose = TRUE, transform_crs = FALSE) {

  sref0 <- try(sref(x), silent = TRUE)

  sp_target <- .get_sref_target(sref)

  crs_in <- units_in <- crs_scale_in <- NULL
  if (!is.null(sp_target)) {
    crs_in <- sp_target$crs
    units_in <- sp_target$units
    crs_scale_in <- sp_target$crs_scale
  }

  crs0 <- units0 <- NA
  crs_scale0 <- NA_real_

  if (!is.null(sref0) && !inherits(sref0, "try-error")) {
    crs0 <- sref0$crs %||% NA
    units0 <- sref0$units %||% NA_character_
    crs_scale0 <- sref0$crs_scale %||% NA_real_
  }

  ## determine target CRS first
  if (is.null(crs_in)) {
    crs_store <- crs0
  } else {
    crs_store <- .normalize_crs_wkt(crs_in)
  }

  crs_changed <- !.same_crs(crs0, crs_store)

  ## determine target units
  if (is.null(units_in)) {
    if (crs_changed) {
      units <- .infer_sref_units(crs_store, verbose = verbose)
    } else if (!.is_na_scalar(units0)) {
      units <- units0
    } else {
      units <- .infer_sref_units(crs_store, verbose = verbose)
    }
  } else {
    units <- units_in
  }

  ## determine target crs_scale
  if (is.null(crs_scale_in)) {

    units_changed <- !(.is_na_scalar(units0) || .is_na_scalar(units)) &&
      !identical(as.character(units0), as.character(units))

    if (!is.na(crs_scale0) && is.finite(crs_scale0) && !crs_changed && !units_changed) {
      crs_scale <- crs_scale0
    } else {
      crs_scale <- .infer_sref_scale(crs_store, units, verbose = verbose)
    }

  } else {
    crs_scale <- crs_scale_in
  }

  if (is.null(crs_scale) || is.na(crs_scale) || !is.finite(crs_scale)) {
    if (!.is_na_scalar(crs_store) || !.is_na_scalar(units)) {
      crs_scale <- 1
    } else {
      crs_scale <- NA_real_
    }
  }

  ## existing sref present: harmonize / transform / scale
  if (!is.null(sref0) && !inherits(sref0, "try-error")) {

    if (!.is_na_scalar(crs0) && !.is_na_scalar(crs_store) && crs_changed) {

      if (isTRUE(transform_crs)) {
        if (isTRUE(verbose)) {
          message("CRS differs; transforming stored coordinates to the target CRS.")
        }
        x <- transform_sref(
          x,
          crs = crs_store,
          units = units,
          crs_scale = crs_scale,
          verbose = FALSE
        )
      } else {
        stop(
          "CRS differs. Replacing the CRS metadata alone would change the meaning ",
          "of the stored coordinates. Use transform_sref() or ",
          "add_sref(..., transform_crs = TRUE)."
        )
      }

    } else if (!is.na(crs_scale0) && is.finite(crs_scale0) &&
               !is.na(crs_scale) && is.finite(crs_scale) &&
               !isTRUE(all.equal(crs_scale0, crs_scale))) {

      x <- scale_sref(
        x,
        scale = crs_scale / crs_scale0,
        verbose = verbose
      )
    }
  }

  sp_target <- create_sref(crs = crs_store,
                           units = units,
                           crs_scale = crs_scale)
  sref(x) <- sp_target
  x
}


##' Transform stored coordinates to a new CRS while preserving locations
##'
##' Transform the stored coordinates of an object from its current spatial
##' reference to a new CRS, and update \code{sref(x)} accordingly.
##'
##' This function should be used when an object already has a valid spatial
##' reference and you want to adopt a different CRS without changing the
##' underlying locations represented by the coordinates. In other words, it
##' changes the spatial coordinate system, not the represented points.
##'
##' The source coordinates are interpreted using the current \code{sref(x)}. The
##' current \code{crs_scale} is used to convert stored coordinates back to CRS
##' units before transformation. After transformation, coordinates are converted
##' to the requested stored units using the target \code{crs_scale}.
##'
##' The new spatial reference can be supplied directly via \code{crs},
##' \code{units}, and \code{crs_scale}, or indirectly via \code{sref} as a
##' sref-like object or another object from which \code{sref(sref)} can be
##' extracted.
##'
##' At present, this function is implemented for point-based objects with
##' coordinate columns named \code{x}/\code{y}, \code{x0}/\code{y0},
##' \code{x1}/\code{y1}, and so on. It is not implemented for
##' \code{admove_grid} or \code{admove_cov}, because reprojection of regular
##' grids and covariate arrays requires resampling rather than simple coordinate
##' transformation.
##'
##' @param x An object with an existing valid \code{sref} and supported stored
##'   coordinates.
##' @param sref Optional object or sref-like specification providing target
##'   spatial-reference fields. This can be:
##'   \itemize{
##'   \item an object of class \code{"admove_sref"},
##'   \item a named list containing one or more of \code{crs}, \code{units}, and
##'     \code{crs_scale},
##'   \item another object from which \code{sref(sref)} can be extracted.
##'   }
##' @param crs Optional target CRS. Can be an \code{sf::crs} object, an EPSG
##'   code, or a character string. If supplied, this takes precedence over
##'   \code{sref$crs}.
##' @param units Optional target stored coordinate units, e.g. \code{"m"},
##'   \code{"km"}, or \code{"degree"}. If \code{NULL}, they are inferred from the
##'   target CRS where possible.
##' @param crs_scale Optional numeric scalar giving the conversion factor from
##'   target CRS units to target stored units. If \code{NULL}, it is inferred
##'   from \code{crs} and \code{units} where possible.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed.
##'
##' @return \code{x} with transformed coordinates and updated \code{sref(x)}.
##'
##' @details
##' Let \eqn{c_{\mathrm{old}}} denote stored coordinates, \eqn{s_{\mathrm{old}}}
##' the old \code{crs_scale}, and \eqn{\mathcal{T}} the CRS transformation from
##' the old CRS to the new CRS. The transformation is applied as:
##' \deqn{
##' c_{\mathrm{new}} =
##' \mathcal{T}\left(c_{\mathrm{old}} / s_{\mathrm{old}}\right) \times
##' s_{\mathrm{new}},
##' }
##' where coordinates are first converted from stored units to source CRS units,
##' transformed to the target CRS, and then converted from target CRS units to
##' target stored units.
##'
##' If the target CRS is the same as the current CRS, the function falls back to
##' \code{\link{scale_sref}} if only the stored coordinate units differ.
##'
##' @seealso \code{\link{add_sref}}, \code{\link{scale_sref}},
##'   \code{\link{create_sref}}, \code{\link{sref}},
##'   \code{\link{.guess_crs_scale}}
##'
##' @examples
##' ## Not run:
##' ## Reproject point coordinates from lon/lat to a projected CRS
##' ## tags <- transform_sref(tags, crs = 3035, units = "km", crs_scale = 0.001)
##'
##' ## Use the spatial reference of another object
##' ## tags <- transform_sref(tags, sref = other_tags)
##' ## End(Not run)
##'
##' @export
transform_sref <- function(x, sref = NULL, crs = NULL, units = NULL,
                            crs_scale = NULL, verbose = TRUE) {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for CRS transformation.")
  }

  ## browser()

  sp0 <- try(sref(x), silent = TRUE)
  if (inherits(sp0, "try-error") || is.null(sp0)) {
    stop("Input object 'x' must already have a valid sref.")
  }

  crs0 <- sp0$crs %||% NA
  units0 <- sp0$units %||% NA_character_
  crs_scale0 <- sp0$crs_scale %||% NA_real_

  if (.is_na_scalar(crs0)) {
    stop("Cannot transform spatial reference because sref(x)$crs is NA.")
  }

  if (is.na(crs_scale0) || !is.finite(crs_scale0)) {
    stop("Cannot transform spatial reference because sref(x)$crs_scale is missing or invalid.")
  }

  sp_target <- .get_sref_target(sref)

  crs_in <- sp_target$crs %||% NULL
  units_in <- sp_target$units %||% NULL
  crs_scale_in <- sp_target$crs_scale %||% NULL

  if (!is.null(crs)) crs_in <- crs
  crs_store <- .normalize_crs_wkt(crs_in)

  ## TODO: use scale_crs if crs_in and sp_target are NULL/NA?

  if (.is_na_scalar(crs_store)) {
    stop("Provide a valid target 'crs' or an object/sref with sref$crs.")
  }

  crs_changed <- !.same_crs(crs0, crs_store)

  if (is.null(units)) {
    if (!is.null(units_in)) {
      units <- units_in
    } else if (crs_changed) {
      units <- .infer_sref_units(crs_store, verbose = verbose)
    } else if (!.is_na_scalar(units0)) {
      units <- units0
    } else {
      units <- .infer_sref_units(crs_store, verbose = verbose)
    }
  }

  if (is.null(crs_scale)) {
    if (!is.null(crs_scale_in)) {
      crs_scale <- crs_scale_in
    } else if (!crs_changed && !.is_na_scalar(units0) && !.is_na_scalar(units) &&
               identical(as.character(units0), as.character(units)) &&
               is.finite(crs_scale0)) {
      crs_scale <- crs_scale0
    } else {
      crs_scale <- .infer_sref_scale(crs_store, units, verbose = verbose)
    }
  }

  if (is.null(crs_scale) || is.na(crs_scale) || !is.finite(crs_scale)) {
    if (!.is_na_scalar(crs_store) || !.is_na_scalar(units)) {
      crs_scale <- 1
    } else {
      crs_scale <- NA_real_
    }
  }

  sp_new <- create_sref(crs = crs_store, units = units, crs_scale = crs_scale)

  if (!crs_changed) {

    if (!isTRUE(all.equal(crs_scale0, crs_scale))) {
      if (isTRUE(verbose)) {
        message("Target CRS is the same; rescaling stored coordinates only.")
      }
      x <- scale_sref(x, scale = crs_scale / crs_scale0, verbose = verbose)
    } else if (isTRUE(verbose)) {
      message("Target CRS and coordinate scaling are unchanged.")
    }

    sref(x) <- sp_new
    return(x)
  }

  x <- .transform_xy_columns(x, sp_old = sp0, sp_new = sp_new)
  sref(x) <- sp_new

  if (isTRUE(verbose)) {
    message("Transformed stored coordinates to target CRS.")
  }

  x
}


##' Change coordinate units by scaling an admove object
##'
##' Multiply the stored coordinates of an \code{admove_*} object by \code{scale}
##' and update its spatial reference accordingly. The spatial reference is
##' updated by:
##' \itemize{
##'   \item multiplying \code{sref(x)$crs_scale} by \code{scale}
##'   \item updating \code{sref(x)$units} where possible (e.g. m <-> km)
##' }
##'
##' @param x An \code{admove_grid}, \code{admove_*tags}, \code{admove_cov}, or
##'   \code{admove_data} object.
##' @param scale Numeric scalar. Coordinates are multiplied by this factor.
##'   For example, to convert meters to kilometers (numerically), use
##'   \code{scale = 0.001}.
##' @param units units
##' @param verbose Logical; if \code{TRUE}, prints informational messages.
##'
##' @return \code{x} with rescaled coordinates and updated \code{sref(x)}.
##'
##' @examples
##'
##' grid <- create_grid(cellsize = 5e3, xrange = c(0, 5e4), yrange = c(0, 5e4))
##'
##' grid <- add_sref(grid, list(crs = 32631, units = "m", crs_scale = 1))
##'
##' grid <- scale_sref(grid, scale = 0.001)
##'
##' units_space(grid)
##'
##' crs_scale(grid)
##'
##'
##' @export
scale_sref <- function(x, scale = 1, units = NULL, verbose = TRUE) {

  use_units <- ifelse(!is.null(units), TRUE, FALSE)

  if (!use_units && (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0)) {
    stop("'scale' must be a single positive finite number.")
  }

  if (use_units) {
    scale <- .guess_crs_scale(sref(x), units)
  }

  ## ---- rescale coordinates --------------------------------------------------
  if (inherits(x, "admove_grid")) {

    x$xygrid <- x$xygrid * scale
    x$cellsize <- x$cellsize * scale
    x$xgr <- x$xgr * scale
    x$ygr <- x$ygr * scale
    x$xrange <- x$xrange * scale
    x$yrange <- x$yrange * scale

  } else if (inherits(x, "admove_tags")) {

    x$x <- x$x * scale
    x$y <- x$y * scale

  } else if (inherits(x, "admove_cov")) {

    ## assumes dimnames store x/y centers as character
    xcen <- as.numeric(attributes(x)$dimnames[[1]])
    ycen <- as.numeric(attributes(x)$dimnames[[2]])
    xcen <- xcen * scale
    ycen <- ycen * scale
    attributes(x)$dimnames[[1]] <- as.character(xcen)
    attributes(x)$dimnames[[2]] <- as.character(ycen)

  } else if (inherits(x, "admove_data")) {

    browser()
    x$x <- x$x * scale
    x$y <- x$y * scale

  } else {
    stop("Only works for admove_grid, admove_tags, admove_cov, or admove_data.")
  }

  ## ---- update spatial reference --------------------------------------------
  sp <- sref(x)

  ## update numeric map CRS-units -> stored-units
  sp$crs_scale <- sp$crs_scale * scale

  ## update units label where possible
  u <- sp$units
  eps <- 1e-10

  u_split <- strsplit(u, "_x_")[[1]]
  if (length(u_split) == 2) {
    u <- u_split[1]
    scale <- scale * as.numeric(u_split[2])  ## same as updated sp$crs_scale
  } else if (length(u_split) != 1) stop("Not sure what to do.")

  if (isTRUE(abs(scale - 0.001) < eps) &&
        u %in% c("meter", "metre", "m")) {
    sp$units <- "km"
  } else if (isTRUE(abs(scale - 1e3) < eps) &&
               u %in% c("kilometer", "kilometre", "km")) {
    sp$units <- "m"
  } else if (isTRUE(abs(scale - 1) > eps)) {
    sp$units <- paste0(u, "_x_", format(scale, scientific = TRUE))
  }

  sref(x) <- sp
  x
}



sref_equal <- function(a, b, tol = 1e-12) {

  if (!inherits(a, "admove_sref")) {
    a <- sref(a)
  }

  if (!inherits(b, "admove_sref")) {
    b <- sref(b)
  }

  ## if (!inherits(a, "admove_sref") || !inherits(b, "admove_sref")) return(FALSE)

  ## CRS compare: use sf when available, else string compare
  crs_a <- a$crs
  crs_b <- b$crs
  if (requireNamespace("sf", quietly = TRUE)) {
    ca <- try(sf::st_crs(crs_a), silent = TRUE)
    cb <- try(sf::st_crs(crs_b), silent = TRUE)
    if (!(inherits(ca, "try-error") || inherits(cb, "try-error") || is.na(ca) || is.na(cb))) {
      ## compare by WKT (stable)
      crs_a <- ca$wkt
      crs_b <- cb$wkt
    }
  }
  if (!identical(crs_a, crs_b)) return(FALSE)

  if (is.na(a$crs_scale) && is.na(b$crs_scale) &&
        is.na(a$units) && is.na(b$units) &&
        is.na(a$crs) && is.na(b$crs)) return(TRUE)

  ## units + crs_scale
  if (!identical(a$units, b$units)) return(FALSE)
  if (!is.na(a$crs_scale) && !is.na(b$crs_scale) &&
      abs(a$crs_scale - b$crs_scale) > tol) return(FALSE)

  TRUE
}





## Internal functions -----------------------------------------------------------

##' Guess a coordinate scaling factor from a spatial reference (sref)
##'
##' Compute a numeric factor that converts coordinates stored under a given
##' \code{admove_sref} to a requested target unit. This extends the previous
##' CRS-only logic by accounting for \code{sref$crs_scale}, i.e. the map from
##' CRS units to stored units.
##'
##' The returned sref is a multiplicative factor \code{f} such that:
##' \deqn{x_{\mathrm{target}} = x_{\mathrm{stored}} \times f}
##'
##' If CRS units are known (via \pkg{sf}) and both the stored units and target
##' units can be converted to meters, then:
##' \deqn{
##' f = \frac{\mathrm{meters\ per\ stored\ unit}}{\mathrm{meters\ per\ target\ unit}}
##' }
##'
##' where meters per stored unit is derived from the CRS units and
##' \code{sref$crs_scale}.
##'
##' @param sref An \code{admove_sref} object (or a list with components
##'   \code{crs}, \code{units}, \code{crs_scale}).
##' @param units Character string describing the desired target units, e.g.
##'   \code{"m"}, \code{"km"}, \code{"mi"}, \code{"nmi"}.
##'
##' @return A numeric scalar scaling factor, or \code{NA_real_} if it cannot be
##'   determined.
##'
.guess_crs_scale <- function(sref, units) {

  if (is.null(sref) || all(is.na(sref))) stop("'sref' must be provided.")
  if (is.null(units) || !length(units)) stop("'units' must be provided.")

  crs <- sref$crs
  units_stored <- sref$units
  crs_scale <- sref$crs_scale

  if (is.null(crs_scale) || length(crs_scale) != 1 ||
        !is.finite(crs_scale) || crs_scale <= 0) {
    return(NA_real_)
  }

  crs_sf <- NULL
  if (!is.null(crs) && requireNamespace("sf", quietly = TRUE)) {
    crs_sf <- try(sf::st_crs(crs), silent = TRUE)
    if (inherits(crs_sf, "try-error") || is.na(crs_sf)) crs_sf <- NULL
  }

  if (units == "degree") {
    if (!is.null(crs_sf) && !is.null(crs_sf$units_gdal) && crs_sf$units_gdal == "degree") {
      return(1)
    } else {
      return(NA_real_)
    }
  }

  ## target units (meters per target unit)
  target_m_per_unit <- .in_m(units)
  if (is.na(target_m_per_unit)) return(NA_real_)

  ## stored units: we derive meters per stored unit from CRS units and crs_scale
  ## crs_scale maps CRS-units -> stored-units, so:
  ## stored_units = crs_units * crs_scale
  ## => 1 stored unit = (1 / crs_scale) CRS units
  units_crs <- NA_character_
  if (!is.null(crs_sf) && !is.null(crs_sf$units_gdal)) units_crs <- crs_sf$units_gdal
  crs_m_per_unit <- .in_m(units_crs)

  if (!is.na(crs_m_per_unit)) {
    stored_m_per_unit <- crs_m_per_unit / crs_scale
  } else {
    ## fallback: if CRS units are unknown, try to use the stored unit label directly
    stored_m_per_unit <- .in_m(units_stored)
    if (is.na(stored_m_per_unit)) return(NA_real_)
  }

  ## return factor to multiply stored coordinates to get target coordinates
  return(stored_m_per_unit / target_m_per_unit)
}

.normalise_unit <- function(u) {
  if (is.null(u) || is.na(u) || !nzchar(u)) return(NA_character_)
  u <- tolower(trimws(as.character(u)))
  u <- gsub("\\.", "", u)
  u <- gsub("s$", "", u)  ## crude plural removal (meters -> meter)
  return(u)
}

## conversion: unit -> meters per unit
.in_m <- function(u) {
  u <- .normalise_unit(u)

  res <- NA

  idx <- grep("_x_", u)
  if (length(idx) > 0) {
    tmp <- strsplit(u, "_x_")[[1]]
    u <- tmp[1]
    res2 <- 1/as.numeric(tmp[2])
  } else {
    res2 <- 1
  }

  ## metric
  if (u %in% c("m", "meter", "metre")) res <- 1
  if (u %in% c("km", "kilometer", "kilometre")) res <- 1e3
  if (u %in% c("cm", "centimeter", "centimetre")) res <- 1e-2
  if (u %in% c("mm", "millimeter", "millimetre")) res <- 1e-3
  if (u %in% c("um", "micrometer", "micrometre")) res <- 1e-6

  ## imperial / US customary
  if (u %in% c("in", "inch")) res <- 0.0254
  if (u %in% c("ft", "foot")) res <- 0.3048
  if (u %in% c("yd", "yard")) res <- 0.9144
  if (u %in% c("mi", "mile", "statutemile")) res <- 1609.344

  ## nautical
  if (u %in% c("nmi", "nauticalmile", "nautical_mile", "nmile")) res <- 1852

  res * res2
}


.get_sref_target <- function(sref) {

  if (is.null(sref)) return(NULL)

  if (inherits(sref, "admove_sref") ||
      any(names(sref) %in% c("crs", "units", "crs_scale"))) {
    return(sref)
  }

  sp <- try(sref(sref), silent = TRUE)
  if (inherits(sp, "try-error")) return(NULL)

  sp
}

.normalize_crs_wkt <- function(crs) {

  if (is.null(crs) || .is_na_scalar(crs)) return(NA)

  if (requireNamespace("sf", quietly = TRUE)) {
    crs_sf <- try(sf::st_crs(crs), silent = TRUE)
    if (!inherits(crs_sf, "try-error") && !is.na(crs_sf)) {
      return(crs_sf$wkt)
    }
  }

  crs
}

.same_crs <- function(x, y) {

  if (.is_na_scalar(x) && .is_na_scalar(y)) return(TRUE)
  if (.is_na_scalar(x) || .is_na_scalar(y)) return(FALSE)

  if (requireNamespace("sf", quietly = TRUE)) {
    x_sf <- try(sf::st_crs(x), silent = TRUE)
    y_sf <- try(sf::st_crs(y), silent = TRUE)

    if (!inherits(x_sf, "try-error") && !inherits(y_sf, "try-error") &&
        !is.na(x_sf) && !is.na(y_sf)) {
      return(identical(x_sf$wkt, y_sf$wkt))
    }
  }

  identical(as.character(x), as.character(y))
}

.infer_sref_units <- function(crs, verbose = TRUE) {

  if (is.null(crs) || .is_na_scalar(crs)) return(NA_character_)

  if (requireNamespace("sf", quietly = TRUE)) {
    crs_sf <- try(sf::st_crs(crs), silent = TRUE)
    if (!inherits(crs_sf, "try-error") && !is.na(crs_sf) &&
        !is.null(crs_sf$units_gdal) && !is.na(crs_sf$units_gdal)) {
      return(crs_sf$units_gdal)
    }
  }

  if (isTRUE(verbose)) {
    message("Could not infer units from CRS. Consider providing spatial 'units' explicitly.")
  }

  NA_character_
}

.infer_sref_scale <- function(crs, units, verbose = TRUE) {

  if (is.null(crs) || .is_na_scalar(crs) || is.null(units) || .is_na_scalar(units)) {
    return(NA_real_)
  }

  out <- try(
    .guess_crs_scale(
      sref = list(crs = crs, units = units, crs_scale = 1),
      units = units
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error") || !is.numeric(out) || length(out) != 1L ||
      !is.finite(out) || is.na(out)) {
    if (isTRUE(verbose)) message("Could not infer 'crs_scale'.")
    return(NA_real_)
  }

  out
}

.xy_pair_names <- function(x) {

  nms <- names(x)
  if (is.null(nms)) {
    return(data.frame(x = character(0), y = character(0), stringsAsFactors = FALSE))
  }

  xcols <- grep("^x[0-9]*$", nms, value = TRUE)
  if (!length(xcols)) {
    return(data.frame(x = character(0), y = character(0), stringsAsFactors = FALSE))
  }

  suffix <- sub("^x", "", xcols)
  ycols <- paste0("y", suffix)
  keep <- ycols %in% nms

  data.frame(
    x = xcols[keep],
    y = ycols[keep],
    stringsAsFactors = FALSE
  )
}

.transform_xy_columns <- function(x, sp_old, sp_new) {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for CRS transformation.")
  }

  if (inherits(x, "admove_grid")) {
    stop("transform_sref() is not implemented for 'admove_grid'. Grid reprojection requires resampling.")
  }

  if (inherits(x, "admove_cov")) {
    stop("transform_sref() is not implemented for 'admove_cov'. Covariate reprojection requires resampling.")
  }

  pairs <- .xy_pair_names(x)

  if (!nrow(pairs)) {
    stop(
      "No supported x/y coordinate columns found. ",
      "Currently supported names are x/y, x0/y0, x1/y1, ..."
    )
  }

  for (i in seq_len(nrow(pairs))) {

    xc <- pairs$x[i]
    yc <- pairs$y[i]

    xv <- x[[xc]]
    yv <- x[[yc]]

    ok <- is.finite(xv) & is.finite(yv)
    if (!any(ok)) next

    pts <- data.frame(
      x = xv[ok] / sp_old$crs_scale,
      y = yv[ok] / sp_old$crs_scale
    )

    sf_pts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = sp_old$crs)
    sf_pts <- sf::st_transform(sf_pts, crs = sp_new$crs)
    cc <- sf::st_coordinates(sf_pts)

    xv_new <- xv
    yv_new <- yv

    xv_new[ok] <- cc[, 1] * sp_new$crs_scale
    yv_new[ok] <- cc[, 2] * sp_new$crs_scale

    x[[xc]] <- xv_new
    x[[yc]] <- yv_new
  }

  x
}


.extract_wkt_name <- function(wkt, key) {
  m <- regexec(paste0(key, '\\["([^"]+)"'), wkt)
  mm <- regmatches(wkt, m)[[1]]
  if (length(mm) >= 2) mm[2] else NA_character_
}

.extract_wkt_name <- function(wkt, key) {
  if (is.null(wkt) || length(wkt) != 1L || is.na(wkt) || !nzchar(wkt)) {
    return(NA_character_)
  }

  m <- regexec(paste0(key, '\\["([^"]+)"'), wkt)
  mm <- regmatches(wkt, m)[[1]]

  if (length(mm) >= 2L) mm[2] else NA_character_
}


.format_sref_short <- function(x, labw = 15) {
  stopifnot(inherits(x, "admove_sref"))

  ## Helper: test whether a scalar field is effectively missing
  is_missing1 <- function(z) {
    is.null(z) || length(z) != 1L || is.na(z) || !nzchar(as.character(z))
  }

  ## Extract raw fields safely
  crs_raw <- if (!is_missing1(x$crs)) as.character(x$crs) else NA_character_
  units_raw <- if (!is_missing1(x$units)) as.character(x$units) else NA_character_
  scale_raw <- if (!is_missing1(x$crs_scale)) as.numeric(x$crs_scale) else NA_real_

  ## Case 1: CRS not specified at all
  if (is.na(crs_raw)) {
    out <- c(
      sprintf(paste0("  %-", labw, "s %s"), "crs:", "not specified"),
      if (!is.na(units_raw))
        sprintf(paste0("  %-", labw, "s %s"), "stored units:", units_raw),
      if (!is.na(scale_raw))
        sprintf(paste0("  %-", labw, "s %s"), "crs scale:",
                format(signif(scale_raw, 4), trim = TRUE))
    )
    return(out)
  }

  ## Defaults from WKT only (works without sf)
  wkt <- crs_raw

  crs_name <- .extract_wkt_name(wkt, "PROJCRS")
  method <- .extract_wkt_name(wkt, "METHOD")
  datum <- .extract_wkt_name(wkt, "DATUM")
  crs_units <- .extract_wkt_name(wkt, "LENGTHUNIT")
  epsg <- NA_integer_

  if (is.na(datum)) {
    datum <- .extract_wkt_name(wkt, "BASEGEOGCRS")
  }

  ## If PROJCRS is "unknown", prefer the projection method
  if (!is.na(crs_name) && tolower(crs_name) == "unknown") {
    crs_name <- NA_character_
  }
  if (is.na(crs_name)) {
    crs_name <- method
  }
  if (is.na(crs_name) || !nzchar(crs_name)) {
    crs_name <- "specified CRS"
  }

  ## Optional enrichment via sf, but never fail if sf is absent or errors
  if (requireNamespace("sf", quietly = TRUE)) {
    crs_obj <- tryCatch(sf::st_crs(crs_raw), error = function(e) NULL)

    if (!is.null(crs_obj)) {
      nm <- tryCatch(crs_obj$Name, error = function(e) NA_character_)
      if (!is.na(nm) && nzchar(nm) && tolower(nm) != "unknown") {
        crs_name <- nm
      }

      ep <- tryCatch(crs_obj$epsg, error = function(e) NA_integer_)
      if (!is.na(ep)) {
        epsg <- ep
      }

      un <- tryCatch(crs_obj$units_gdal, error = function(e) NA_character_)
      if (!is.na(un) && nzchar(un)) {
        crs_units <- un
      }
    }
  }

  ## Identify CRS type
  crs_id <- if (!is.na(epsg)) paste0("EPSG:", epsg) else "custom"

  ## Scale text
  scale_txt <- NA_character_
  if (!is.na(scale_raw)) {
    if (!is.na(crs_units) && !is.na(units_raw)) {
      scale_txt <- sprintf("1 %s = %s %s",
                           crs_units,
                           format(signif(scale_raw, 4), trim = TRUE),
                           units_raw)
    } else {
      scale_txt <- format(signif(scale_raw, 4), trim = TRUE)
    }
  }

  ## Build output
  crs_line <- if (!is.na(crs_id)) {
    sprintf(paste0("  %-", labw, "s %s [%s]"), "crs:", crs_name, crs_id)
  } else {
    sprintf(paste0("  %-", labw, "s %s"), "crs:", crs_name)
  }

  out <- c(
    crs_line,
    if (!is.na(datum))
      sprintf(paste0("  %-", labw, "s %s"), "datum:", datum),
    if (!is.na(crs_units))
      sprintf(paste0("  %-", labw, "s %s"), "crs units:", crs_units),
    if (!is.na(units_raw))
      sprintf(paste0("  %-", labw, "s %s"), "stored units:", units_raw),
    if (!is.na(scale_txt))
      sprintf(paste0("  %-", labw, "s %s"), "crs scale:", scale_txt)
  )

  out
}



## s3 methods ----------------------------------------------------------------------

##' sf-compatible CRS accessor for admove objects
##'
##' S3 method for \code{\link[sf]{st_crs}} that returns the CRS stored in
##' \code{sref(x)$crs} as an \code{sf::crs} object.
##'
##' @param x An object with spatial reference (\code{sref}).
##' @param ... Further arguments (unused).
##'
##' @return An \code{sf::crs} object.
##'
##' @details
##' Requires the \pkg{sf} package to be installed.
##'
##' @name admove-st_crs
NULL

##' @rdname admove-st_crs
##' @export
st_crs.admove_grid <- function(x, ...) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for st_crs(). Please install it.")
  }
  sf::st_crs(crs(x))
}

##' @rdname admove-st_crs
##' @export
st_crs.admove_cov <- function(x, ...) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for st_crs(). Please install it.")
  }
  sf::st_crs(crs(x))
}

##' @rdname admove-st_crs
##' @export
st_crs.admove_tags <- function(x, ...) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for st_crs(). Please install it.")
  }
  sf::st_crs(crs(x))
}
