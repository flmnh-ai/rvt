# Internal utility functions (not exported)

roll_2d <- function(m, dy, dx) {
  nr <- nrow(m); nc <- ncol(m)
  m[((seq_len(nr) - 1L - dy) %% nr) + 1L,
    ((seq_len(nc) - 1L - dx) %% nc) + 1L]
}

pad_reflect <- function(m, pad) {
  nr <- nrow(m); nc <- ncol(m)
  m[c(seq.int(pad + 1L, 2L), seq_len(nr), seq.int(nr - 1L, nr - pad)),
    c(seq.int(pad + 1L, 2L), seq_len(nc), seq.int(nc - 1L, nc - pad))]
}

pad_edge <- function(m, pad) {
  nr <- nrow(m); nc <- ncol(m)
  m[c(rep(1L, pad), seq_len(nr), rep(nr, pad)),
    c(rep(1L, pad), seq_len(nc), rep(nc, pad))]
}

roll_fill <- function(m, shift, axis) {
  rolled <- if (axis == 1L) roll_2d(m, shift, 0L) else roll_2d(m, 0L, shift)
  na_mask <- is.na(rolled)
  rolled[na_mask] <- m[na_mask]
  rolled
}

normalize_lin <- function(x, minimum, maximum) {
  ok <- !is.na(x)
  x[ok & x > maximum] <- maximum
  x[ok & x < minimum] <- minimum
  x <- (x - minimum) / (maximum - minimum)
  x[ok & x > 1] <- 1
  x[ok & x < 0] <- 0
  x
}

horizon_shift_vector <- function(n_dir = 16L, radius_pixels = 10L,
                                 min_radius = 1) {
  angles <- (2 * pi / n_dir) * (0:(n_dir - 1L))
  dir_x <- cos(angles); dir_y <- sin(angles)
  scale <- 3
  n_radii <- as.integer((radius_pixels - min_radius) * scale) + 1L
  radii <- seq(0, n_radii - 1L) / scale + min_radius
  lapply(seq_len(n_dir), function(i) {
    x_int <- round(dir_x[i] * radii); y_int <- round(dir_y[i] * radii)
    coords <- unique(complex(real = x_int, imaginary = y_int))
    sx <- as.integer(Re(coords)); sy <- as.integer(Im(coords))
    dist <- sqrt(as.numeric(sx)^2 + as.numeric(sy)^2)
    ord <- order(dist); valid <- dist[ord] > 0; ord <- ord[valid]
    list(shift_row = sx[ord], shift_col = sy[ord], distance = dist[ord])
  })
}
