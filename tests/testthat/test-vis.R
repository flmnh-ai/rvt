# Basic tests for visualization functions

# Synthetic DEM: Gaussian mound on a flat plane
make_test_dem <- function(n = 50) {
  outer(1:n, 1:n, function(x, y)
    10 * exp(-((x - n/2)^2 + (y - n/2)^2) / (n * 2)))
}

test_that("rvt_slope_aspect returns correct dimensions", {
  dem <- make_test_dem()
  sa <- rvt_slope_aspect(dem, 1, 1)
  expect_equal(dim(sa$slope), dim(dem))
  expect_equal(dim(sa$aspect), dim(dem))
})

test_that("slope is non-negative", {
  sa <- rvt_slope_aspect(make_test_dem(), 1, 1, "radian")
  expect_true(all(sa$slope >= 0, na.rm = TRUE))
})

test_that("slope units convert correctly", {
  dem <- make_test_dem()
  sa_rad <- rvt_slope_aspect(dem, 1, 1, "radian")
  sa_deg <- rvt_slope_aspect(dem, 1, 1, "degree")
  expect_equal(sa_rad$slope * 180 / pi, sa_deg$slope, tolerance = 1e-10)
})

test_that("rvt_hillshade returns values in [0, 1]", {
  hs <- rvt_hillshade(make_test_dem(), 1, 1)
  expect_true(all(hs >= 0 & hs <= 1, na.rm = TRUE))
})

test_that("rvt_sky_view_factor returns correct dimensions", {
  dem <- make_test_dem(30)
  result <- rvt_sky_view_factor(dem, 1, n_dir = 8L, r_max = 5L,
                                compute_svf = TRUE, compute_opns = TRUE)
  expect_equal(dim(result$svf), dim(dem))
  expect_equal(dim(result$opns), dim(dem))
})

test_that("SVF values are in [0, 1]", {
  dem <- make_test_dem(30)
  result <- rvt_sky_view_factor(dem, 1, n_dir = 8L, r_max = 5L)
  expect_true(all(result$svf >= 0 & result$svf <= 1, na.rm = TRUE))
})

test_that("rvt_slrm returns correct dimensions", {
  dem <- make_test_dem(50)
  slrm_out <- rvt_slrm(dem, radius_cell = 10L)
  expect_equal(dim(slrm_out), dim(dem))
})

test_that("sLRM of flat surface is ~zero", {
  flat <- matrix(5, 50, 50)
  slrm_out <- rvt_slrm(flat, radius_cell = 10L)
  expect_true(all(abs(slrm_out) < 1e-10))
})

test_that("sLRM detects a mound as positive", {
  dem <- make_test_dem(50)
  slrm_out <- rvt_slrm(dem, radius_cell = 15L)
  # Centre of the mound should be positive (locally high)
  expect_true(slrm_out[25, 25] > 0)
})

test_that("rvt_terrain returns valid presets", {
  for (name in c("general", "flat", "steep")) {
    t <- rvt_terrain(name)
    expect_true(is.list(t))
    expect_true(all(c("hs_sun_el", "svf_r_max", "svf_noise",
                       "slope", "svf", "opns", "ve_default") %in% names(t)))
  }
})

test_that("rvt_vat returns values in [0, 1]", {
  dem <- make_test_dem(40)
  result <- rvt_vat(dem, 1, 1, terrain = "flat",
                    svf_n_dir = 8L, verbose = FALSE)
  expect_equal(dim(result), dim(dem))
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("NAs are preserved through pipeline", {
  dem <- make_test_dem(40)
  dem[20, 20] <- NA
  result <- rvt_vat(dem, 1, 1, terrain = "general",
                    svf_n_dir = 8L, verbose = FALSE)
  expect_true(is.na(result[20, 20]))
})
