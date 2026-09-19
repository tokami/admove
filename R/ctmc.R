## Main functions ---------------------------------------------------------------------

##' Group tags into release events for the CTMC engine
##'
##' @description
##' `set_release_events()` records which tags the CTMC engine may treat as one
##' release event. Tags released in the same place at the same time share a
##' single forward pass through the likelihood instead of repeating an identical
##' one each, which is where nearly all of the cost of a conventional-tag CTMC
##' fit sits.
##'
##' @param dat A data list of class `admove_data`, as produced by [setup_data()].
##' @param tag_types Tag types eligible for grouping, as the one-letter codes.
##'   Default `"c"` (mark-recapture).
##' @param t_tol Time tolerance for merging releases, in the model's time units
##'   (see [units_time()]). `0` (the default) merges only exactly simultaneous
##'   releases.
##' @param dist_tol Distance tolerance for merging releases, in the spatial
##'   units of the reference system (see [units_space()]). `0` (the default)
##'   merges only releases already in the same grid cell.
##' @param verbose Logical; if `TRUE`, informative messages may be printed.
##'
##' @details
##' Grouping is only valid for tags whose distribution is never conditioned on
##' an observation, i.e. tag types with `conf$do_update` set to `FALSE`. With
##' the default configuration that is mark-recapture tags: their position is a
##' pure forecast from release, so it depends on the release alone and every tag
##' released together follows the same trajectory. Archival and mark-resight
##' tags are updated at each observation and diverge immediately, so they are
##' never grouped whatever is set here. A model may mix both; the split is per
##' tag, not per model.
##'
##' **Calling this function is optional.** With the defaults
##' (`conf$ctmc_groups = "auto"`) the engine already groups exactly, on the
##' release cell and release time. It is needed only to widen the grouping with
##' `t_tol`/`dist_tol`, or to restrict `tag_types`.
##'
##' Exact grouping is free: it changes nothing but the arithmetic order, because
##' the tags really do share a trajectory. **Non-zero tolerances are an
##' approximation** — releases that differ are treated as if they did not, so
##' the likelihood changes. The function therefore warns, naming how many tags
##' were merged and the largest displacement it applied, so the size of the
##' approximation is visible rather than implied.
##'
##' What is stored is a policy, not a materialised grouping: [admove()] runs
##' [check_tags()] again, which can drop rows, recompute `ic` and rename ids, so
##' anything resolved earlier would go stale silently. Ids that cannot be found
##' at fit time fall back to exact grouping, with one message.
##'
##' @return
##' `dat` with a `release_events` element of class `"admove_release_events"`.
##'
##' @examples
##' dat <- set_release_events(skjepo$sim$dat)
##'
##' @export
set_release_events <- function(dat, tag_types = "c", t_tol = 0, dist_tol = 0,
                               verbose = TRUE) {

  .check_class(dat, "admove_data")

  if (!is.character(tag_types) || length(tag_types) == 0 ||
        !all(tag_types %in% c("d", "s", "c", "a"))) {
    stop("'tag_types' must be one or more of \"d\", \"s\", \"c\" and \"a\".",
         call. = FALSE)
  }
  .check_pos_or_zero <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0) {
      stop("'", nm, "' must be a single finite number >= 0.", call. = FALSE)
    }
  }
  .check_pos_or_zero(t_tol, "t_tol")
  .check_pos_or_zero(dist_tol, "dist_tol")

  exact <- t_tol == 0 && dist_tol == 0
  keys <- NULL

  if (!exact) {

    tags <- dat$tags
    if (is.null(tags) || nrow(tags) == 0) {
      stop("'dat' holds no tags to group.", call. = FALSE)
    }
    rel <- do.call(rbind, lapply(split(tags, tags$id), function(z) z[1, ]))
    rel <- rel[as.character(rel$tag_type) %in% tag_types, , drop = FALSE]

    if (nrow(rel) > 0) {
      ## Greedy single pass in time order: each release joins the first open
      ## cluster within both tolerances, otherwise opens one. Order-dependent
      ## by construction, which is why the displacements are reported rather
      ## than assumed small.
      o <- order(rel$t)
      rel <- rel[o, , drop = FALSE]
      cl <- integer(nrow(rel))
      ct <- numeric(0); cx <- numeric(0); cy <- numeric(0)
      dmax <- 0; tmax <- 0
      for (k in seq_len(nrow(rel))) {
        hit <- 0L
        if (length(ct) > 0) {
          dt <- abs(rel$t[k] - ct)
          dd <- sqrt((rel$x[k] - cx)^2 + (rel$y[k] - cy)^2)
          ok <- which(dt <= t_tol & dd <= dist_tol)
          if (length(ok) > 0) {
            hit <- ok[which.min(dd[ok])]
            tmax <- max(tmax, dt[hit])
            dmax <- max(dmax, dd[hit])
          }
        }
        if (hit == 0L) {
          ct <- c(ct, rel$t[k]); cx <- c(cx, rel$x[k]); cy <- c(cy, rel$y[k])
          hit <- length(ct)
        }
        cl[k] <- hit
      }
      keys <- stats::setNames(paste0("re", cl), as.character(rel$id))

      nmerged <- sum(table(cl) > 1L)
      ntags_m <- sum(cl %in% as.integer(names(table(cl))[table(cl) > 1L]))
      if (nmerged > 0) {
        warning("set_release_events(): merged ", ntags_m, " tags into ",
                nmerged, " release event(s) that are not identical, moving ",
                "releases by up to ", signif(tmax, 3), " in time and ",
                signif(dmax, 3), " in distance. This changes the likelihood; ",
                "t_tol = 0 and dist_tol = 0 group only tags that really do ",
                "share a release.", call. = FALSE)
      }
    }
  }

  if (verbose && !any(tag_types %in% "c")) {
    message("tag_types does not include \"c\". With the default ",
            "conf$do_update, mark-recapture tags are the only ones that can ",
            "be grouped, so nothing will be.")
  }

  dat$release_events <- structure(list(tag_types = tag_types,
                                       t_tol = t_tol,
                                       dist_tol = dist_tol,
                                       exact = exact,
                                       keys = keys),
                                  class = "admove_release_events")
  return(dat)
}


##' @rdname print-admove
##' @method print admove_release_events
##' @export
print.admove_release_events <- function(x, ...) {
  cat("Release-event grouping for the CTMC engine\n")
  cat("  tag types :", paste(x$tag_types, collapse = ", "), "\n")
  if (isTRUE(x$exact)) {
    cat("  grouping  : exact (same release cell and time)\n")
  } else {
    cat("  grouping  : approximate, t_tol =", x$t_tol,
        " dist_tol =", x$dist_tol, "\n")
    cat("  tags mapped:", length(x$keys), "\n")
  }
  invisible(x)
}


## Internal functions ---------------------------------------------------------------------

## The CTMC engine of nll(). Extracted from the engine branch so that the
## per-tag pass and the release-event pass can share the same leaves
## (.ctmc_generator, .ctmc_step, .ctmc_obs_dist) -- two loops evaluating the
## same model is a drift risk, and sharing the leaves is what contains it.
##
## Every function here re-declares the three RTMB overloads: they are local to
## the function that declares them, and R's JIT drops them otherwise. Leaving
## one out fails loudly, but only on a path a small fixture may not reach.


## The generator at dt = 1: the instantaneous rate matrix, mass-balanced.
##
## Every off-diagonal entry is exactly linear in dt -- the drift terms carry it
## as a factor, the diffusion as hD = D * dt, fill_inst_mat() applies either
## pos(v) = 0.5*(v + abs(v)) (positively homogeneous) or 0.5*v, and the diagonal
## is a row sum of those. So Mstar(slice, dt) == dt * Q(slice) and the step is
## applied in .ctmc_step() instead. Do not fold dt back in here: it is what
## makes the cache key the slice alone, and what makes a subdivided step exact.
## See dev/code_notes.org, "Caching the CTMC generator".
##
## Depends on the parameters and on t only through the covariate and spline
## slices t falls in -- never on the tag.
.ctmc_generator <- function(ctx, t, template) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  "diag<-" <- RTMB::ADoverload("diag<-")

  dat <- ctx$dat
  nc <- ctx$nc
  xygrid <- ctx$xygrid
  nextTo <- ctx$nextTo
  next_dist <- ctx$next_dist

  ## Set to zero
  if (!is.null(template)) {
    Zstar <- Astar <- Dstar <- template
    Zstar@x[] <- Astar@x[] <- Dstar@x[] <- 0
  } else {
    Zstar <- Astar <- Dstar <- RTMB::matrix(0, nc, nc)
  }

  ## taxis
  if (dat$use_taxis) {
    move <- ctx$kappa * ctx$habi$tax$grad(xygrid, t)  ## distance / time
    Zstar <- fill_inst_mat(Zstar, move, nextTo, next_dist, dat$drift_scheme)
  }

  ## advection
  if (dat$use_advection) {
    move <- cbind(ctx$habi$adv_x$val(xygrid, t),
                  ctx$habi$adv_y$val(xygrid, t))  ## distance / time
    Astar <- fill_inst_mat(Astar, move, nextTo, next_dist, dat$drift_scheme)
  }

  ## diffusion
  D <- exp(ctx$habi$dif$val(xygrid, t)) ## distance^2 / time
  for (k in 1:4) {
    j <- k + 1
    ind <- which(!is.na(nextTo[, j]))
    Dstar[cbind(ind, nextTo[ind, j])] <- D[ind] / next_dist[k]^2
    ## (distance^2 / time) / (distance^2) = 1 / time
  }

  ## Movement rates
  Mstar <- Zstar + Astar + Dstar

  ## Mass balance
  Mstar[cbind(1:nc, 1:nc)] <- 0
  Mstar[cbind(1:nc, 1:nc)] <- -RTMB::rowSums(Mstar)

  return(Mstar)
}


## Resolve the lattice / grouping mode. NULL-safe: a conf built before this
## setting existed reaches nll() without it, and admove() never calls
## check_conf(). Unlike drift_scheme, the fallback here is the NEW behaviour --
## an old object then refits faster without being edited, and reproducing the
## old lattice is an explicit "off".
.ctmc_mode <- function(dat) {
  if (is.null(dat$ctmc_groups)) return("auto")
  return(dat$ctmc_groups)
}


## Every covariate and spline slice boundary on the absolute time axis, so that
## no lattice step spans two slices and the generator is exactly constant
## within a step. That is what makes subdividing a step exact -- with Q constant,
## exp(Q*(a+b)) == exp(Q*a) %*% exp(Q*b) -- and hence what lets tags released
## together share one pass. See dev/code_notes.org, "The aligned CTMC lattice".
.ctmc_breaks <- function(dat) {

  tv <- c(dat$time_cov, dat$time_spline)
  if (length(tv) == 0) return(numeric(0))

  ## nll() is handed the split() list, but this is also callable on the
  ## data.frame form setup_data() stores
  tt <- if (is.data.frame(dat$tags)) {
    dat$tags$t
  } else {
    unlist(lapply(dat$tags, function(z) z$t))
  }
  tt <- tt[is.finite(tt)]
  if (length(tt) == 0) return(numeric(0))
  tmin <- min(tt)
  tmax <- max(tt)

  seas <- c(dat$seasonal_cov, dat$seasonal_spline)
  per <- dat$period
  wrap <- !is.null(per) && length(per) == 1L && is.finite(per) && per > 0
  kseq <- if (wrap) floor(tmin / per):ceiling(tmax / per) else 0

  br <- numeric(0)
  any_seasonal <- FALSE
  for (i in seq_along(tv)) {
    v <- as.numeric(tv[[i]])
    if (length(v) == 0) next
    if (isTRUE(seas[i]) && wrap) {
      ## seasonal values are phases in [0, period): replicate them over the
      ## periods the tags span
      v <- as.numeric(outer(v, kseq * per, "+"))
      any_seasonal <- TRUE
    }
    br <- c(br, v)
  }

  ## t %% period is discontinuous at the period boundaries themselves, and
  ## seasonal_cov carries no guarantee that 0 is among its values
  if (any_seasonal) br <- c(br, kseq * per)

  br <- br[is.finite(br)]
  return(sort(unique(br)))
}


## The covariate and spline slice indices t falls in, as one key. Taken from
## the habi object itself rather than from the breakpoints, so it cannot drift
## from what val()/grad() actually read. All four habi objects are built with
## the same time arguments, so one of them serves.
.ctmc_slice_key <- function(habi, t) {
  paste(habi$dif$slice(t), collapse = ",")
}


## Memoise the dt = 1 generator on the slice t falls in.
##
## One cache per taping, and it must stay that way. A cached AD object holds
## handles into the tape that was active when it was built; replaying it under
## a different tape returns a silently wrong number -- no error, no warning.
## This is the same invariant the spline derivative tape depends on; see
## dev/code_notes.org, "Two invariants the cache depends on".
##
## The live trap here is that tmb_all <- c(dat, conf) is captured by the closure
## MakeADFun() tapes and persists across every obj$fn() call and every retape.
## Writing the cache into dat, into a package-level memo, or into anything else
## that outlives this call would therefore look like it works and be wrong.
.ctmc_gen_cache <- function(ctx, template) {

  cache <- new.env(parent = emptyenv())

  Qof <- function(t) {
    key <- .ctmc_slice_key(ctx$habi, t)
    Q <- cache[[key]]
    if (!is.null(Q)) return(Q)
    Q <- .ctmc_generator(ctx, t, template)
    assign(key, Q, envir = cache)
    return(Q)
  }

  list(Q = Qof, env = cache)
}


## Propagate the cell distribution over one step of length dt.
.ctmc_step <- function(ctx, Q, dt, last_dist) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  ## counters is an environment, so this counts across the whole call
  ctx$counters$nstep <- ctx$counters$nstep + 1L

  Mstar <- dt * Q

  if (identical(ctx$dat$ctmc_method, "expav")) {

    pred_dist <- as.vector(RTMB::expAv(Mstar,
                               last_dist,
                               transpose = TRUE,
                               uniformization = TRUE,
                               rescale_freq = 1,
                               trace = FALSE))

  } else {

    M <- Matrix::expm(Mstar)
    pred_dist <- as.vector(RTMB::matrix(last_dist, 1, ctx$nc) %*% M)

  }

  return(pred_dist)
}


## The distribution one observation row is scored against: an exact-cell
## indicator, or the cell-integrated normal when that row carries observation
## error. Depends on the row only, never on the tag's history, which is what
## lets members of a release event be scored off one shared forward pass.
.ctmc_obs_dist <- function(ctx, tag, ind_obs_j, ind_tt, ev, last_ev) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  dat <- ctx$dat
  xygrid <- ctx$xygrid
  cs <- ctx$cs
  sdO <- ctx$sdO

  ## default
  this_dist <- rep(0, nrow(xygrid))
  this_dist[tag$ic[ind_obs_j]] <- 1

  ## obs uncertainty
  ## see the note in the KF branch: "all but the last observation" is
  ## a comparison against the last observation event, not against the
  ## number of integration time points
  if ((dat$obs_var_type[ind_tt] == 1 && ev[ind_obs_j] != last_ev) ||
        dat$obs_var_type[ind_tt] == 2 ||
        dat$obs_var_type[ind_tt] == 3) {

    xLo <- xygrid[,1] - cs[1] / 2
    xUp <- xygrid[,1] + cs[1] / 2
    yLo <- xygrid[,2] - cs[2] / 2
    yUp <- xygrid[,2] + cs[2] / 2

    xObs <- tag$x[ind_obs_j]
    yObs <- tag$y[ind_obs_j]
    if (dat$obs_var_type[ind_tt] == 3) {
      sdx <- tag$sdx[ind_obs_j]
      sdy <- tag$sdy[ind_obs_j]
    } else {
      sdx <- sdO[1,ind_tt]
      sdy <- sdO[2,ind_tt]
    }
    px <- RTMB::pnorm(xUp, mean = xObs, sd = sdx) -
      RTMB::pnorm(xLo, mean = xObs, sd = sdx)
    py <- RTMB::pnorm(yUp, mean = yObs, sd = sdy) -
      RTMB::pnorm(yLo, mean = yObs, sd = sdy)
    pxy <- px * py
    pxy <- pxy / sum(pxy)

    if (all(!is.na(pxy))) {
      this_dist <- pxy
    }
  }

  return(this_dist)
}


## Tag types eligible for release-event grouping. Only a policy, so it stays
## correct across check_tags(), which runs inside admove() and can drop rows,
## recompute ic and rename ids -- anything materialised earlier goes stale
## silently.
.ctmc_group_types <- function(dat) {
  re <- dat$release_events
  if (!is.null(re) && !is.null(re$tag_types)) return(re$tag_types)
  return("c")
}


## The release-event key of tag i, or NULL when the tag cannot share a pass.
##
## A tag is groupable only if its distribution is never conditioned on an
## observation: one update and its state leaves its group-mates', so the shared
## pass is invalid from that point on. That is a property of the tag TYPE via
## conf$do_update, tested over every row rather than row 1.
##
## Read dat$do_update, never tag$update: check_tags() writes that column
## (R/tags.R) and nll() has never read it, so a rule written against it would
## ignore the configuration and silently group updating tags.
.ctmc_groupable <- function(ctx, i) {

  dat <- ctx$dat
  tag <- dat$tags[[i]]

  if (nrow(tag) < 2L) return(NULL)

  ## a deactivated or unlocated release makes the whole tag a no-op or an
  ## error in the per-tag pass; leave both there, unchanged
  if (isTRUE(tag$use[1] == 0)) return(NULL)
  if (is.na(tag$ic[1])) return(NULL)

  tt <- as.integer(tag$tag_type)
  if (any(is.na(tt))) return(NULL)

  ## %in% FALSE is NA-safe: an out-of-range tag type gives NA and is not
  ## grouped, which is the safe direction
  if (!all(dat$do_update[tt] %in% FALSE)) return(NULL)

  if (!all(as.character(tag$tag_type) %in% .ctmc_group_types(dat))) return(NULL)

  ## dt_min is per tag (the median step enters it), so members of one event
  ## must agree on it or they do not share a lattice -- hence it is in the key
  dt_min <- min(dat$min_dt, median(diff(sort(unique(tag$t)))))
  if (!is.finite(dt_min) || dt_min <= 0) return(NULL)

  ## a widened grouping supplies its own key; dt_min stays in it either way,
  ## since members that disagree on it do not share a lattice
  re <- dat$release_events
  if (!is.null(re) && !is.null(re$keys)) {
    k <- re$keys[[as.character(tag$id[1])]]
    if (!is.null(k)) return(paste(k, format(dt_min, digits = 15)))
  }

  return(paste(tag$ic[1],
               format(tag$t[1], digits = 15),
               format(dt_min, digits = 15)))
}


## Split the tags into release events and everything else. Exact by default:
## same release cell, same release time, same dt_min.
.ctmc_partition <- function(ctx) {

  ntags <- length(ctx$dat$tags)

  keys <- vapply(seq_len(ntags), function(i) {
    k <- .ctmc_groupable(ctx, i)
    if (is.null(k)) NA_character_ else k
  }, character(1))

  singles <- which(is.na(keys))
  groups <- list()

  ok <- which(!is.na(keys))
  if (length(ok) > 0) {
    for (m in split(ok, keys[ok])) {
      ## a group of one costs exactly what the per-tag pass costs
      if (length(m) < 2L) {
        singles <- c(singles, m)
      } else {
        groups[[length(groups) + 1L]] <- list(members = m)
      }
    }
  }

  return(list(groups = groups, singles = sort(singles)))
}


## One tag, propagated on its own lattice. The reference implementation: the
## release-event pass must reproduce it exactly for the tags it takes over.
.ctmc_tag_pass <- function(ctx, i, loglik_tags, gen, time_mode, breaks) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  dat <- ctx$dat
  tag <- dat$tags[[i]]

  ## ambiguous final observation -- see the KF branch above
  ev <- .tag_events(tag)
  last_ev <- max(ev)
  amb <- sum(ev == last_ev) > 1L
  amb_lw <- NULL

  ## unique(): candidate positions may share a time
  dt_min <- min(dat$min_dt, median(diff(sort(unique(tag$t)))))
  if (time_mode == "fill_gaps" && (!is.finite(dt_min) || dt_min <= 0)) {
    return(loglik_tags)
  }

  out <- .build_time_breaks(tag$t,
                            mode = time_mode,
                            dt_min = dt_min,
                            dt = dat$dt,
                            eps = 0.1,
                            breaks = breaks)

  ts <- out$ts
  dts <- out$dts
  nts <- out$nts
  observed <- out$observed
  if (length(observed) == 0) return(loglik_tags)

  ind_tag_type <- as.integer(tag$tag_type)

  if (tag$use[1] == 0) return(loglik_tags)

  ## Distribution probability
  last_dist <- rep(0, ctx$nc)
  last_dist[tag$ic[1]] <- 1

  if (nts < 2) stop("Something went wrong (nts < 2).")

  ## Loop over time
  for (t in 2:nts) {

    dt <- dts[t-1]

    ## dist prob after move
    pred_dist <- .ctmc_step(ctx, gen$Q(ts[t-1]), dt, last_dist)

    if (t %in% observed) {

      ind_obs <- which(observed == t) + 1

      ## multiple observation in same time
      for (j in seq_along(ind_obs)) {
        ind_obs_j <- ind_obs[j]
        ind_tt <- ind_tag_type[ind_obs_j]

        if (is.na(tag$ic[ind_obs_j]) || tag$use[ind_obs_j] == 0) {
          last_dist <- pred_dist
          next()
        }

        this_dist <- .ctmc_obs_dist(ctx, tag, ind_obs_j, ind_tt, ev, last_ev)

        update_dist <- pred_dist * this_dist

        if (amb && ev[ind_obs_j] == last_ev) {

          ## mixture over the candidate positions; no update, exactly as in
          ## the KF branch (the ambiguous event is the last one)
          term <- log(tag$prob[ind_obs_j]) + log(sum(update_dist))
          amb_lw <- if (is.null(amb_lw)) term else c(amb_lw, term)
          last_dist <- pred_dist

        } else {

          ## likelihood
          loglik_tags[i] <- loglik_tags[i] + log(sum(update_dist))

          ## update
          if (isTRUE(dat$do_update[ind_tt])) {
            last_dist <- update_dist / sum(update_dist)
          } else {
            last_dist <- pred_dist
          }
        }
      }
    } else {
      last_dist <- pred_dist
    }
  }

  ## log sum_k prob_k * density_k, accumulated stably (see .logsumexp_ad)
  if (amb && !is.null(amb_lw)) {
    loglik_tags[i] <- loglik_tags[i] + .logsumexp_ad(amb_lw)
  }

  return(loglik_tags)
}


## One release event: a single forward pass serving every member.
##
## Valid only because no member is ever updated (.ctmc_groupable()), so the
## distribution depends on the release alone and each member just reads its own
## likelihood off the shared trajectory at its own observation times.
##
## The lattice is built from the UNION of the members' observation times, so it
## refines the lattice each member would get alone. That refinement is exact
## rather than approximate only because the lattice is aligned to the slice
## boundaries: with Q constant across a step, exp(Q*(a+b)) == exp(Q*a) exp(Q*b).
## See dev/code_notes.org, "The aligned CTMC lattice".
.ctmc_event_pass <- function(ctx, grp, loglik_tags, gen, time_mode, breaks) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  dat <- ctx$dat
  members <- grp$members
  nm <- length(members)

  ## per-member bookkeeping, all off the tape
  info <- lapply(members, function(i) {
    tg <- dat$tags[[i]]
    ev <- .tag_events(tg)
    list(tag = tg, ev = ev, last_ev = max(ev),
         amb = sum(ev == max(ev)) > 1L,
         tt = as.integer(tg$tag_type))
  })

  first <- info[[1]]$tag
  dt_min <- min(dat$min_dt, median(diff(sort(unique(first$t)))))

  ## the union of every member's times: the shared lattice
  t_all <- sort(unique(unlist(lapply(info, function(z) z$tag$t))))

  out <- .build_time_breaks(t_all,
                            mode = time_mode,
                            dt_min = dt_min,
                            dt = dat$dt,
                            eps = 0.1,
                            breaks = breaks)

  ts <- out$ts
  dts <- out$dts
  nts <- out$nts
  if (nts < 2) stop("Something went wrong (nts < 2).")

  ## Which rows are scored at which lattice index. Sorted by (k, member, row)
  ## so the floating-point summation order is reproducible.
  obs <- list()
  for (m in seq_len(nm)) {
    tg <- info[[m]]$tag
    for (r in seq_len(nrow(tg))[-1]) {
      obs[[length(obs) + 1L]] <- c(k = which.min(abs(ts - tg$t[r])),
                                   m = m, r = r)
    }
  }
  obs <- do.call(rbind, obs)
  obs <- obs[order(obs[, "k"], obs[, "m"], obs[, "r"]), , drop = FALSE]

  by_k <- vector("list", nts)
  for (q in seq_len(nrow(obs))) {
    kk <- obs[q, "k"]
    by_k[[kk]] <- c(by_k[[kk]], q)
  }

  ## Distribution probability
  last_dist <- rep(0, ctx$nc)
  last_dist[first$ic[1]] <- 1

  ## one slot per member: a shared vector would pool the mixtures across
  ## members and give a finite, plausible, wrong likelihood
  amb_lw <- vector("list", nm)

  for (t in 2:nts) {

    ## dist prob after move
    pred_dist <- .ctmc_step(ctx, gen$Q(ts[t-1]), dts[t-1], last_dist)

    for (q in by_k[[t]]) {

      m <- obs[q, "m"]
      r <- obs[q, "r"]
      z <- info[[m]]
      tg <- z$tag
      ind_tt <- z$tt[r]

      if (is.na(tg$ic[r]) || tg$use[r] == 0) next()

      this_dist <- .ctmc_obs_dist(ctx, tg, r, ind_tt, z$ev, z$last_ev)

      update_dist <- pred_dist * this_dist

      if (z$amb && z$ev[r] == z$last_ev) {

        ## mixture over the candidate positions; no update, exactly as in
        ## the per-tag pass
        term <- log(tg$prob[r]) + log(sum(update_dist))
        amb_lw[[m]] <- if (is.null(amb_lw[[m]])) term else c(amb_lw[[m]], term)

      } else {

        ## straight into the member's own entry: identical to the statement in
        ## .ctmc_tag_pass(), so there is no index mapping to get wrong
        loglik_tags[members[m]] <- loglik_tags[members[m]] +
          log(sum(update_dist))
      }
    }

    ## never updated: groupability guarantees no member conditions the state
    last_dist <- pred_dist
  }

  ## log sum_k prob_k * density_k, accumulated stably (see .logsumexp_ad)
  for (m in seq_len(nm)) {
    if (info[[m]]$amb && !is.null(amb_lw[[m]])) {
      loglik_tags[members[m]] <- loglik_tags[members[m]] +
        .logsumexp_ad(amb_lw[[m]])
    }
  }

  return(loglik_tags)
}


## Entry point from nll(). Returns loglik_tags with the CTMC contributions added.
.ctmc_loglik <- function(ctx, loglik_tags) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  dat <- ctx$dat

  time_mode <- ifelse(is.null(dat$dt) || is.na(dat$dt),
                      "fill_gaps", "fixed_dt")

  template <- if (identical(dat$ctmc_method, "expav")) {
    make_mstar_template(ctx$nextTo, ad = TRUE)
  } else {
    NULL
  }

  ## Plain integers, not AD: only non-AD code writes them, so they are free.
  ## They make the cache and the partition testable on their mechanism rather
  ## than on wall time, which cannot be asserted without flakiness.
  ctx$counters <- new.env(parent = emptyenv())
  ctx$counters$nstep <- 0L

  gen <- .ctmc_gen_cache(ctx, template)

  ## "off" reproduces the per-tag lattice: .build_time_breaks() falls back to
  ## build_time() when there are no breaks to insert
  mode <- .ctmc_mode(dat)
  breaks <- if (identical(mode, "off")) numeric(0) else .ctmc_breaks(dat)

  part <- if (identical(mode, "auto")) {
    .ctmc_partition(ctx)
  } else {
    list(groups = list(), singles = seq_len(length(dat$tags)))
  }

  ctx$counters$ngroup <- length(part$groups)
  ctx$counters$ngrouped <- sum(vapply(part$groups,
                                      function(g) length(g$members),
                                      integer(1)))

  if (isTRUE(dat$dbg)) {
    ## every tag in exactly one place, and every tag somewhere: an index that
    ## drifts here moves a likelihood term onto the wrong tag without changing
    ## the total, so only per-tag output would show it
    all_ix <- sort(c(unlist(lapply(part$groups, function(g) g$members)),
                     part$singles))
    stopifnot(identical(as.integer(all_ix), seq_len(length(dat$tags))))
  }

  for (g in part$groups) {
    loglik_tags <- .ctmc_event_pass(ctx, g, loglik_tags, gen, time_mode, breaks)
  }

  for (i in part$singles) {
    loglik_tags <- .ctmc_tag_pass(ctx, i, loglik_tags, gen, time_mode, breaks)
  }

  return(list(loglik_tags = loglik_tags,
              ngen = length(ls(gen$env)),
              nstep = ctx$counters$nstep,
              ngroup = ctx$counters$ngroup,
              ngrouped = ctx$counters$ngrouped))
}
