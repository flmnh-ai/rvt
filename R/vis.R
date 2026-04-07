# ============================================================================
# Visualization algorithms: slope, hillshade, SVF, openness, sLRM
#
# Faithful port of rvt-py (github.com/EarthObservation/RVT_py)
# Original: Kokalj, Maroh, Ostir, Zaksek, Coz (ZRC SAZU / Univ. Ljubljana)
#
# References:
#   Kokalj & Somrak 2019. Remote Sensing 11(7): 747.
#   Zaksek, Ostir & Kokalj 2011. Remote Sensing 3: 398-415.
# ============================================================================


# ---- Terrain presets -------------------------------------------------------

#' Terrain visualization presets
#'
#' Returns a named list of default parameters for a given terrain type.
#' These match the presets in rvt-py's `default_terrains_settings.json`,
#' with the addition of `ve_default` (suggested vertical exaggeration).
#'
#' @param name One of `"general"`, `"flat"`, `"steep"`.
#' @return Named list with elements: `hs_sun_el`, `svf_r_max`, `svf_noise`,
#'   `slope` (2-element range), `svf` (2-element range), `opns` (2-element
#'   range), `ve_default`.
#' @export
#' @examples
#' rvt_terrain("flat")
#' rvt_terrain("general")$svf_r_max
rvt_terrain <- function(name = c("general", "flat", "steep")) {
  name <- match.arg(name)
  switch(name,
    general = list(hs_sun_el = 35,  svf_r_max = 10L, svf_noise = 0L,
                   slope = c(0, 50),  svf = c(0.7, 1),   opns = c(68, 93),
                   ve_default = 1),
    flat    = list(hs_sun_el = 15,  svf_r_max = 20L, svf_noise = 3L,
                   slope = c(0, 15),  svf = c(0.9, 1),   opns = c(85, 93),
                   ve_default = 3),
    steep   = list(hs_sun_el = 55,  svf_r_max = 10L, svf_noise = 0L,
                   slope = c(0, 60),  svf = c(0.55, 1),  opns = c(55, 95),
                   ve_default = 1)
  )
}


# ---- Slope & Aspect -------------------------------------------------------

#' Slope and aspect from a DEM matrix
#'
#' Computes slope and aspect using central differences on a 3x3 neighbourhood.
#'
#' @param dem Numeric matrix of elevations.
#' @param res_x Pixel size in X direction (metres).
#' @param res_y Pixel size in Y direction (metres).
#' @param units Output units: `"radian"` (default), `"degree"`, or `"percent"`.
#' @param ve_factor Vertical exaggeration applied before computation.
#' @return List with `$slope` and `$aspect` matrices.
#' @export
rvt_slope_aspect <- function(dem, res_x = 1, res_y = 1,
                             units = "radian", ve_factor = 1) {
  dem <- dem * ve_factor
  dem_pad <- pad_edge(dem, 1L)
  dzdx <- ((roll_fill(dem_pad, 1L, 2L) - roll_fill(dem_pad, -1L, 2L)) / 2) / res_x
  dzdy <- ((roll_fill(dem_pad, -1L, 1L) - roll_fill(dem_pad, 1L, 1L)) / 2) / res_y
  tan_slope <- sqrt(dzdx^2 + dzdy^2)
  slope_out <- switch(units,
    percent = tan_slope * 100,
    degree  = atan(tan_slope) * (180 / pi),
    radian  = atan(tan_slope)
  )
  dzdy_safe <- dzdy; dzdy_safe[dzdy_safe == 0] <- 1e-8
  aspect_out <- atan2(dzdx, dzdy_safe)
  if (units == "degree") aspect_out <- aspect_out * (180 / pi)
  rows <- 2:(nrow(dem_pad) - 1); cols <- 2:(ncol(dem_pad) - 1)
  slope_out <- slope_out[rows, cols]; aspect_out <- aspect_out[rows, cols]
  na_mask <- is.na(dem)
  slope_out[na_mask] <- NA_real_; aspect_out[na_mask] <- NA_real_
  list(slope = slope_out, aspect = aspect_out)
}


# ---- Hillshade ------------------------------------------------------------

#' Analytical hillshade from a DEM matrix
#'
#' @param dem Numeric matrix of elevations.
#' @param res_x,res_y Pixel size in metres.
#' @param sun_azimuth Sun azimuth in degrees (default 315).
#' @param sun_elevation Sun elevation in degrees (default 35).
#' @param ve_factor Vertical exaggeration.
#' @param slope,aspect Pre-computed slope/aspect in radians (optional).
#' @return Numeric matrix of hillshade values \[0, 1\].
#' @export
rvt_hillshade <- function(dem, res_x, res_y,
                          sun_azimuth = 315, sun_elevation = 35,
                          ve_factor = 1, slope = NULL, aspect = NULL) {
  if (is.null(slope) || is.null(aspect)) {
    sa <- rvt_slope_aspect(dem, res_x, res_y, "radian", ve_factor)
    slope <- sa$slope; aspect <- sa$aspect
  }
  sz <- pi / 2 - sun_elevation * pi / 180
  sa <- sun_azimuth * pi / 180
  hs <- cos(sz) * cos(slope) + sin(sz) * sin(slope) * cos(aspect - sa)
  hs[!is.na(hs) & hs < 0] <- 0
  hs
}


# ---- SVF / Openness -------------------------------------------------------

#' Core SVF/Openness computation (internal, optimised)
#'
#' @param height Padded height matrix (elevation * ve / resolution).
#' @param r_max Maximum search radius in pixels.
#' @param r_min Minimum search radius in pixels.
#' @param n_dir Number of search directions.
#' @param compute_svf Compute sky-view factor?
#' @param compute_opns Compute positive openness?
#' @param verbose Print progress?
#' @return List with `$svf` and/or `$opns` matrices.
#' @keywords internal
rvt_svf_compute <- function(height, r_max, r_min, n_dir,
                            compute_svf = TRUE, compute_opns = FALSE,
                            verbose = FALSE) {
  nr_pad <- nrow(height); nc_pad <- ncol(height)
  nr_out <- nr_pad - 2L * r_max
  nc_out <- nc_pad - 2L * r_max
  moves <- horizon_shift_vector(n_dir, r_max, r_min)

  out_rows <- (r_max + 1L):(r_max + nr_out)
  out_cols <- (r_max + 1L):(r_max + nc_out)

  height_center <- height[out_rows, out_cols]

  has_na <- anyNA(height)

  if (compute_svf)  svf_acc  <- matrix(0, nr_out, nc_out)
  if (compute_opns) opns_acc <- matrix(0, nr_out, nc_out)

  for (i_dir in seq_len(n_dir)) {
    if (verbose && (i_dir %% max(1L, n_dir %/% 4L) == 1L))
      message(sprintf("  Horizon: direction %d / %d", i_dir, n_dir))

    mv <- moves[[i_dir]]
    n_steps <- length(mv$distance)

    max_slope <- matrix(-1000, nr_out, nc_out)

    for (i_rad in seq_len(n_steps)) {
      dy <- mv$shift_row[i_rad]
      dx <- mv$shift_col[i_rad]
      d  <- mv$distance[i_rad]

      shifted <- height[out_rows - dy, out_cols - dx]

      if (has_na) {
        max_slope <- pmax(max_slope, (shifted - height_center) / d, na.rm = TRUE)
      } else {
        max_slope <- pmax(max_slope, (shifted - height_center) / d)
      }
    }

    max_slope <- atan(max_slope)

    if (compute_svf)
      svf_acc <- svf_acc + (1 - sin(pmax(max_slope, 0)))
    if (compute_opns)
      opns_acc <- opns_acc + max_slope
  }

  result <- list()
  if (compute_svf)  result$svf  <- svf_acc / n_dir
  if (compute_opns) result$opns <- (pi / 2 - opns_acc / n_dir) * 180 / pi
  result
}


#' Sky-view factor and positive openness from a DEM matrix
#'
#' Computes the sky-view factor (proportion of visible hemisphere) and/or
#' positive openness for each pixel by tracing horizon angles in multiple
#' directions.
#'
#' @param dem Numeric matrix of elevations.
#' @param resolution Pixel size in metres (single value; assumed square).
#' @param n_dir Number of search directions (default 16).
#' @param r_max Search radius in pixels (default 10).
#' @param noise Noise removal level 0--3 (default 0). Higher values skip
#'   near-pixel samples to reduce noise from micro-topography.
#' @param ve_factor Vertical exaggeration.
#' @param compute_svf Compute sky-view factor? (default TRUE)
#' @param compute_opns Compute positive openness? (default FALSE)
#' @param verbose Print progress?
#' @return List with `$svf` and/or `$opns` matrices.
#' @export
rvt_sky_view_factor <- function(dem, resolution,
                                n_dir = 16L, r_max = 10L, noise = 0L,
                                ve_factor = 1,
                                compute_svf = TRUE, compute_opns = FALSE,
                                verbose = FALSE) {
  stopifnot(noise %in% 0:3)
  sc_r_min <- c(0, 10, 20, 40)
  r_min <- max(round(r_max * sc_r_min[noise + 1L] * 0.01), 1)
  na_mask <- is.na(dem)
  height <- dem * ve_factor / resolution
  height <- pad_reflect(height, r_max)
  result <- rvt_svf_compute(height, r_max, r_min, n_dir,
                            compute_svf, compute_opns, verbose)
  for (nm in names(result)) result[[nm]][na_mask] <- NA_real_
  result
}


# ---- Simple Local Relief Model ---------------------------------------------

#' Simple Local Relief Model (sLRM)
#'
#' Subtracts a mean-filtered DEM from the original to isolate local elevation
#' anomalies. Positive values indicate locally high areas (mounds, ridges);
#' negative values indicate locally low areas (pits, depressions).
#'
#' Uses a fast box filter via cumulative sums (no external dependencies).
#'
#' @param dem Numeric matrix of elevations.
#' @param radius_cell Radius of the mean filter in pixels (default 20).
#' @param ve_factor Vertical exaggeration (default 1).
#' @return Numeric matrix of local relief values (centred near 0).
#' @export
rvt_slrm <- function(dem, radius_cell = 20L, ve_factor = 1) {
  dem <- dem * ve_factor
  nr <- nrow(dem); nc <- ncol(dem)

  # Mean filter via cumulative sums (fast, no dependencies)
  pad <- pad_edge(dem, radius_cell)
  pad[is.na(pad)] <- 0
  nr_p <- nrow(pad); nc_p <- ncol(pad)

  # Count non-NA cells for proper averaging
  valid <- pad_edge(!is.na(dem) * 1, radius_cell)

  # Row-wise cumsum then col-wise cumsum for box filter
  cs <- apply(pad, 2, cumsum)
  cs <- rbind(rep(0, nc_p), cs)
  box <- cs[seq(2 * radius_cell + 2L, nr_p + 1L), ] -
         cs[seq(1L, nr_p - 2L * radius_cell), ]

  cs_v <- apply(valid, 2, cumsum)
  cs_v <- rbind(rep(0, ncol(valid)), cs_v)
  box_v <- cs_v[seq(2 * radius_cell + 2L, nrow(valid) + 1L), ] -
           cs_v[seq(1L, nrow(valid) - 2L * radius_cell), ]

  cs2 <- t(apply(box, 1, cumsum))
  cs2 <- cbind(rep(0, nrow(cs2)), cs2)
  box_sum <- cs2[, seq(2 * radius_cell + 2L, ncol(cs2))] -
             cs2[, seq(1L, ncol(cs2) - 2L * radius_cell - 1L)]

  cs2_v <- t(apply(box_v, 1, cumsum))
  cs2_v <- cbind(rep(0, nrow(cs2_v)), cs2_v)
  box_cnt <- cs2_v[, seq(2 * radius_cell + 2L, ncol(cs2_v))] -
             cs2_v[, seq(1L, ncol(cs2_v) - 2L * radius_cell - 1L)]

  box_cnt[box_cnt == 0] <- 1
  mean_dem <- box_sum / box_cnt

  out <- dem - mean_dem[seq_len(nr), seq_len(nc)]
  out[is.na(dem)] <- NA_real_
  out
}
