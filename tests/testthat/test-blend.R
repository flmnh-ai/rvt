# Tests for blend modes

test_that("blend_dispatch supports all modes", {
  a <- matrix(0.5, 3, 3)
  bg <- matrix(0.7, 3, 3)
  for (mode in c("normal", "multiply", "screen", "overlay", "soft_light")) {
    result <- blend_dispatch(mode, a, bg)
    expect_true(is.matrix(result))
    expect_equal(dim(result), c(3, 3))
  }
})

test_that("blend_dispatch errors on unknown mode", {
  expect_error(blend_dispatch("bogus", matrix(1), matrix(1)),
               "Unknown blend mode")
})

test_that("normal blend returns active layer", {
  a <- matrix(0.3, 2, 2)
  bg <- matrix(0.9, 2, 2)
  expect_equal(blend_dispatch("normal", a, bg), a)
})

test_that("multiply blend is commutative", {
  a <- matrix(c(0.2, 0.5, 0.8, 1.0), 2, 2)
  bg <- matrix(c(0.3, 0.6, 0.9, 0.1), 2, 2)
  expect_equal(blend_dispatch("multiply", a, bg),
               blend_dispatch("multiply", bg, a))
})

test_that("screen blend of black returns background", {
  bg <- matrix(0.5, 2, 2)
  result <- blend_dispatch("screen", matrix(0, 2, 2), bg)
  expect_equal(result, bg)
})

test_that("normalize_vis inverts slope", {
  slope <- matrix(c(0, 25, 50), 1, 3)
  norm <- normalize_vis(slope, "slope", 0, 50)
  expect_equal(norm[1, 1], 1)  # 0 degrees = max brightness
  expect_equal(norm[1, 3], 0)  # 50 degrees = min brightness
})

test_that("normalize_vis does not invert non-slope", {
  vals <- matrix(c(0, 0.5, 1), 1, 3)
  norm <- normalize_vis(vals, "hillshade", 0, 1)
  expect_equal(norm[1, 1], 0)
  expect_equal(norm[1, 3], 1)
})
