# rvt

Pure-R implementation of the [Relief Visualization Toolbox](https://github.com/EarthObservation/RVT_py) (RVT) algorithms for terrain visualization from LiDAR data. Faithful port of rvt-py (Kokalj et al.) with no Python dependency.

## Installation

```r
# install.packages("remotes")
remotes::install_github("flmnh-ai/rvt")
```

## What it does

Computes terrain visualizations from a DEM matrix or `terra::SpatRaster`:

- **Hillshade** — analytical hillshading with configurable sun position
- **Slope & aspect** — central-difference gradient in radians, degrees, or percent
- **Sky-view factor (SVF)** — proportion of visible hemisphere per pixel
- **Positive openness** — angular measure of terrain exposure
- **Simple local relief model (sLRM)** — mean-filter subtraction to isolate local anomalies
- **VAT composite** — the full Visualization for Archaeological Topography blend (Kokalj & Somrak 2019)

All core functions operate on plain numeric matrices with zero dependencies. `terra` is in Suggests for optional `SpatRaster` convenience wrappers.

## Usage

```r
library(rvt)

# Synthetic DEM with a mound
dem <- outer(1:100, 1:100, function(x, y)
  5 * exp(-((x - 50)^2 + (y - 50)^2) / 200))

# Full VAT composite
result <- rvt_vat(dem, res_x = 1, res_y = 1, terrain = "flat")

# Individual layers
sa <- rvt_slope_aspect(dem, 1, 1, units = "degree")
hs <- rvt_hillshade(dem, 1, 1, sun_elevation = 35)
svf_result <- rvt_sky_view_factor(dem, resolution = 1, r_max = 10L)
slrm_result <- rvt_slrm(dem, radius_cell = 20L)
```

### With terra

```r
library(terra)
library(rvt)

r <- rast("dem.tif")
vat_raster <- vat(r, terrain = "flat")
svf_raster <- svf(r, r_max = 15L)
plot(vat_raster)
```

### Terrain presets

Three presets match rvt-py's defaults, controlling sun elevation, SVF radius, noise filtering, and normalization ranges:

| Preset | Best for | VE default |
|--------|----------|-----------|
| `"general"` | Mixed terrain | 1 |
| `"flat"` | Subtle features, floodplains | 3 |
| `"steep"` | Mountainous terrain | 1 |

```r
rvt_terrain("flat")
#> $hs_sun_el [1] 15
#> $svf_r_max [1] 20
#> ...
```

## Blend modes

The VAT composite uses a specific layer stack with blend modes matching the rvt-py "Archaeological" preset. You can also use `blend_dispatch()` directly:

```r
blend_dispatch("multiply", layer_a, layer_b)
# Supported: normal, multiply, screen, overlay, soft_light, luminosity
```

## Differences from rvt-py

This package is a faithful port with two intentional improvements:

- **Automatic vertical exaggeration.** `rvt_vat()` and `rvt_vat_combined()` default `ve_factor = NULL`, which applies the terrain preset's recommended VE (e.g. 3 for `"flat"`). In rvt-py, VE always defaults to 1 regardless of preset. Pass `ve_factor = 1` explicitly to match rvt-py behaviour.
- **Square pixels assumed for SVF/openness.** `rvt_sky_view_factor()` takes a single `resolution` parameter, same as rvt-py. Non-square pixels are not supported for horizon tracing.

## References

- Kokalj, Z. & Somrak, M. (2019). Why not a single image? Combining visualizations to facilitate fieldwork and on-screen mapping. *Remote Sensing*, 11(7), 747.
- Zakšek, K., Oštir, K. & Kokalj, Ž. (2011). Sky-view factor as a relief visualization technique. *Remote Sensing*, 3, 398–415.
- [rvt-py](https://github.com/EarthObservation/RVT_py) — the original Python implementation by ZRC SAZU / University of Ljubljana.

## License

Apache License 2.0
