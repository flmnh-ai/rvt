# ============================================================================
# terra SpatRaster convenience wrappers
#
# These require the terra package (in Suggests, not Imports).
# All computation is done by the matrix-based functions in vis.R / blend.R;
# these just handle raster <-> matrix conversion and CRS/extent bookkeeping.
# ============================================================================

.check_terra <- function() {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for SpatRaster methods. ",
         "Install it with install.packages('terra').", call. = FALSE)
  }
}

.rast_to_mat <- function(r) {
  terra::as.matrix(r[[1]], wide = TRUE)
}

.mat_to_rast <- function(m, template) {
  r <- terra::rast(m)
  terra::ext(r) <- terra::ext(template)
  terra::crs(r) <- terra::crs(template)
  r
}


#' VAT from a SpatRaster DEM
#'
#' Convenience wrapper around [rvt_vat()] that accepts and returns
#' [terra::SpatRaster] objects.
#'
#' @param dem A [terra::SpatRaster] (single layer).
#' @param terrain Terrain preset: `"general"`, `"flat"`, or `"steep"`.
#' @param ... Additional arguments passed to [rvt_vat()].
#' @return A [terra::SpatRaster] with values in \[0, 1\].
#' @export
vat <- function(dem, terrain = "general", ...) {
  .check_terra()
  m <- .rast_to_mat(dem); res <- terra::res(dem)
  .mat_to_rast(rvt_vat(m, res[1], res[2], terrain, ...), dem)
}

#' VAT Combined from a SpatRaster DEM
#'
#' @param dem A [terra::SpatRaster] (single layer).
#' @param general_opacity Opacity of the general-terrain layer (0--100).
#' @param ... Additional arguments passed to [rvt_vat_combined()].
#' @return A [terra::SpatRaster] with values in \[0, 1\].
#' @export
vat_combined <- function(dem, general_opacity = 50, ...) {
  .check_terra()
  m <- .rast_to_mat(dem); res <- terra::res(dem)
  .mat_to_rast(rvt_vat_combined(m, res[1], res[2], general_opacity, ...), dem)
}

#' Sky-view factor from a SpatRaster DEM
#'
#' @param dem A [terra::SpatRaster] (single layer).
#' @param n_dir Number of search directions (default 16).
#' @param r_max Search radius in pixels (default 10).
#' @param noise Noise removal level 0--3 (default 0).
#' @param ... Additional arguments passed to [rvt_sky_view_factor()].
#' @return A [terra::SpatRaster].
#' @export
svf <- function(dem, n_dir = 16L, r_max = 10L, noise = 0L, ...) {
  .check_terra()
  m <- .rast_to_mat(dem); res <- terra::res(dem)[1]
  r <- rvt_sky_view_factor(m, res, n_dir, r_max, noise, compute_svf = TRUE, ...)
  .mat_to_rast(r$svf, dem)
}

#' Positive openness from a SpatRaster DEM
#'
#' @param dem A [terra::SpatRaster] (single layer).
#' @param n_dir Number of search directions (default 16).
#' @param r_max Search radius in pixels (default 10).
#' @param noise Noise removal level 0--3 (default 0).
#' @param ... Additional arguments passed to [rvt_sky_view_factor()].
#' @return A [terra::SpatRaster] with openness in degrees.
#' @export
positive_openness <- function(dem, n_dir = 16L, r_max = 10L, noise = 0L, ...) {
  .check_terra()
  m <- .rast_to_mat(dem); res <- terra::res(dem)[1]
  r <- rvt_sky_view_factor(m, res, n_dir, r_max, noise,
                           compute_svf = FALSE, compute_opns = TRUE, ...)
  .mat_to_rast(r$opns, dem)
}

#' Simple Local Relief Model from a SpatRaster DEM
#'
#' @param dem A [terra::SpatRaster] (single layer).
#' @param radius_cell Mean filter radius in pixels (default 20).
#' @param ve_factor Vertical exaggeration (default 1).
#' @return A [terra::SpatRaster].
#' @export
slrm <- function(dem, radius_cell = 20L, ve_factor = 1) {
  .check_terra()
  m <- .rast_to_mat(dem)
  .mat_to_rast(rvt_slrm(m, radius_cell, ve_factor), dem)
}
