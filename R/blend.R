# ============================================================================
# Blend modes and VAT composite
#
# Blend mode formulas match rvt-py's blend_func.py exactly.
# VAT blend order matches rvt-py's default_blender_combinations.json
# "Archaeological" preset.
# ============================================================================


# ---- Blend modes (internal) ------------------------------------------------

blend_normal    <- function(a, bg) a
blend_multiply  <- function(a, bg) a * bg
blend_screen    <- function(a, bg) 1 - (1 - a) * (1 - bg)

blend_soft_light <- function(a, bg) {
  result <- bg
  lo <- which(a < 0.5 & !is.na(a))
  hi <- which(a >= 0.5 & !is.na(a))
  result[lo] <- 2 * bg[lo] * a[lo] + bg[lo]^2 * (1 - 2 * a[lo])
  result[hi] <- 2 * bg[hi] * (1 - a[hi]) + sqrt(bg[hi]) * (2 * a[hi] - 1)
  result
}

blend_overlay <- function(a, bg) {
  result <- bg
  hi <- which(bg > 0.5 & !is.na(bg))
  lo <- which(bg <= 0.5 & !is.na(bg))
  result[hi] <- 1 - (1 - 2 * (bg[hi] - 0.5)) * (1 - a[hi])
  result[lo] <- 2 * bg[lo] * a[lo]
  result
}

blend_luminosity <- function(a, bg) {
  if (is.list(a)) a <- 0.3 * a[[1]] + 0.59 * a[[2]] + 0.11 * a[[3]]
  if (!is.list(bg)) return(a)
  d <- a - (0.3 * bg[[1]] + 0.59 * bg[[2]] + 0.11 * bg[[3]])
  r <- bg[[1]] + d; g <- bg[[2]] + d; b <- bg[[3]] + d
  min_c <- pmin(r, g, b); max_c <- pmax(r, g, b)
  lo <- which(min_c < 0); hi <- which(max_c > 1)
  if (length(lo) > 0) {
    den <- a[lo] - min_c[lo]; den[den == 0] <- 1e-10
    r[lo] <- a[lo] + (r[lo] - a[lo]) * a[lo] / den
    g[lo] <- a[lo] + (g[lo] - a[lo]) * a[lo] / den
    b[lo] <- a[lo] + (b[lo] - a[lo]) * a[lo] / den
  }
  if (length(hi) > 0) {
    den <- max_c[hi] - a[hi]; den[den == 0] <- 1e-10
    r[hi] <- a[hi] + (r[hi] - a[hi]) * (1 - a[hi]) / den
    g[hi] <- a[hi] + (g[hi] - a[hi]) * (1 - a[hi]) / den
    b[hi] <- a[hi] + (b[hi] - a[hi]) * (1 - a[hi]) / den
  }
  list(r, g, b)
}


# ---- Exported blend utilities ----------------------------------------------

#' Dispatch a blend mode by name
#'
#' Applies a named blend mode to an active layer and background layer.
#' Supported modes: `"normal"`, `"multiply"`, `"screen"`, `"overlay"`,
#' `"soft_light"`, `"luminosity"`.
#'
#' @param mode Character string naming the blend mode.
#' @param active Numeric matrix (the layer being blended in).
#' @param bg Numeric matrix or RGB list (the background).
#' @return Blended result (matrix or RGB list for luminosity).
#' @export
blend_dispatch <- function(mode, active, bg) {
  switch(tolower(mode),
    "normal"     = blend_normal(active, bg),
    "multiply"   = blend_multiply(active, bg),
    "screen"     = blend_screen(active, bg),
    "overlay"    = blend_overlay(active, bg),
    "soft_light" = blend_soft_light(active, bg),
    "luminosity" = blend_luminosity(active, bg),
    stop("Unknown blend mode: ", mode)
  )
}

#' Apply opacity to a blended result
#'
#' Linear interpolation between blended and background layers.
#'
#' @param blended The blended layer.
#' @param bg The background layer.
#' @param opacity Opacity 0--100 (or 0--1).
#' @return Composited result.
#' @keywords internal
apply_opacity <- function(blended, bg, opacity) {
  if (opacity > 1) opacity <- opacity / 100
  blended * opacity + bg * (1 - opacity)
}

#' Normalize a visualization layer for blending
#'
#' Clips values to \[min, max\], scales to \[0, 1\], and inverts slope layers
#' (so dark = steep, matching the rvt-py convention).
#'
#' @param image Numeric matrix.
#' @param vis_type Character label (e.g. `"slope"`, `"hillshade"`, `"svf"`).
#' @param min_val,max_val Normalization range.
#' @return Numeric matrix in \[0, 1\].
#' @export
normalize_vis <- function(image, vis_type, min_val, max_val) {
  norm <- normalize_lin(image, min_val, max_val)
  if (tolower(vis_type) %in% c("slope gradient", "slope"))
    norm <- 1 - norm
  norm
}


# ---- VAT composite --------------------------------------------------------

#' Visualization for Archaeological Topography (VAT)
#'
#' Blends hillshade, slope, positive openness, and sky-view factor into a
#' single composite following Kokalj & Somrak (2019). Layer order and blend
#' modes match the rvt-py "Archaeological" preset.
#'
#' @param dem Numeric matrix of elevations.
#' @param res_x Pixel size in X (metres).
#' @param res_y Pixel size in Y (metres).
#' @param terrain Terrain preset: `"general"`, `"flat"`, or `"steep"`.
#' @param svf_n_dir SVF search directions (default 16).
#' @param svf_r_max SVF search radius in pixels (`NULL` = terrain preset).
#' @param svf_noise SVF noise level 0--3 (`NULL` = terrain preset).
#' @param hs_sun_el Hillshade sun elevation in degrees (`NULL` = terrain preset).
#' @param ve_factor Vertical exaggeration. If `NULL` (default), uses the
#'   terrain preset's `ve_default` (e.g. 3 for `"flat"`). This differs from
#'   rvt-py, which always defaults to 1 regardless of preset.
#' @param verbose Print progress messages? (default TRUE)
#' @return Numeric matrix of blended values \[0, 1\].
#' @export
#' @examples
#' # Synthetic 100x100 DEM with a mound
#' dem <- outer(1:100, 1:100, function(x, y)
#'   5 * exp(-((x - 50)^2 + (y - 50)^2) / 200))
#' result <- rvt_vat(dem, 1, 1, terrain = "flat", verbose = FALSE)
rvt_vat <- function(dem, res_x, res_y,
                    terrain = "general",
                    svf_n_dir = 16L,
                    svf_r_max = NULL, svf_noise = NULL, hs_sun_el = NULL,
                    ve_factor = NULL, verbose = TRUE) {
  t <- rvt_terrain(terrain)
  svf_r_max <- svf_r_max %||% t$svf_r_max
  svf_noise <- svf_noise %||% t$svf_noise
  hs_sun_el <- hs_sun_el %||% t$hs_sun_el
  ve_factor <- ve_factor %||% t$ve_default
  t0 <- proc.time()

  if (verbose) message("[1/4] Slope & aspect...")
  sa_rad <- rvt_slope_aspect(dem, res_x, res_y, "radian", ve_factor)
  slope_deg <- sa_rad$slope * (180 / pi)

  if (verbose) message("[2/4] Hillshade...")
  hs <- rvt_hillshade(dem, res_x, res_y, 315, hs_sun_el, ve_factor,
                      sa_rad$slope, sa_rad$aspect)

  if (verbose) message("[3/4] SVF & openness...")
  svf_opns <- rvt_sky_view_factor(dem, res_x, svf_n_dir, svf_r_max, svf_noise,
                                  ve_factor, TRUE, TRUE, verbose)

  if (verbose) message("[4/4] Blending...")
  n_hs   <- normalize_vis(hs,            "hillshade", 0,          1)
  n_slp  <- normalize_vis(slope_deg,     "slope",     t$slope[1], t$slope[2])
  n_opns <- normalize_vis(svf_opns$opns, "openness",  t$opns[1],  t$opns[2])
  n_svf  <- normalize_vis(svf_opns$svf,  "svf",       t$svf[1],   t$svf[2])

  # Blend bottom-to-top (same order as rvt-py VAT)
  rendered <- n_hs
  rendered <- apply_opacity(blend_dispatch("luminosity", n_slp, rendered),
                            rendered, 50)
  rendered <- apply_opacity(blend_dispatch("overlay", n_opns, rendered),
                            rendered, 50)
  rendered <- apply_opacity(blend_dispatch("multiply", n_svf, rendered),
                            rendered, 25)
  rendered[!is.na(rendered) & rendered < 0] <- 0
  rendered[!is.na(rendered) & rendered > 1] <- 1

  if (verbose) message(sprintf("Done. (%.1fs)", (proc.time() - t0)[3]))
  rendered
}


#' VAT Combined (general + flat terrain blend)
#'
#' Computes VAT for both general and flat terrain presets and alpha-blends
#' them together. Useful for terrain that has both broad relief and subtle
#' flat-area features.
#'
#' @param dem Numeric matrix of elevations.
#' @param res_x,res_y Pixel size in metres.
#' @param general_opacity Opacity of the general-terrain layer (0--100,
#'   default 50).
#' @param svf_n_dir SVF search directions.
#' @param ve_factor Vertical exaggeration. If `NULL` (default), each sub-call
#'   uses its terrain preset's `ve_default`.
#' @param verbose Print progress?
#' @return Numeric matrix of blended values \[0, 1\].
#' @export
rvt_vat_combined <- function(dem, res_x, res_y,
                             general_opacity = 50,
                             svf_n_dir = 16L, ve_factor = NULL,
                             verbose = TRUE) {
  if (verbose) message("=== VAT General ===")
  vat_gen <- rvt_vat(dem, res_x, res_y, "general", svf_n_dir,
                     ve_factor = ve_factor, verbose = verbose)
  if (verbose) message("=== VAT Flat ===")
  vat_flat <- rvt_vat(dem, res_x, res_y, "flat", svf_n_dir,
                      ve_factor = ve_factor, verbose = verbose)
  if (verbose) message("=== Combining ===")
  out <- apply_opacity(vat_gen, vat_flat, general_opacity)
  out[!is.na(out) & out < 0] <- 0
  out[!is.na(out) & out > 1] <- 1
  out
}
