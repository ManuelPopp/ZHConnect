#!/usr/bin/env/Rscript --vanilla
#>------------------------------------<
#
# Script name: rasterize_core_areas
#
# Author: Manuel R. Popp
# Email: manuel.popp@wsl.ch
#
# Date Created: 2023-10-04
#
# ---------------------------
#
# Description: Rasterise the "Kerngebiete" shapefiles.
#
# Notes: The current version includes areas from neighbouring cantons within a
# given buffer area.
#
#>----------------------------------------------------------------------------<|
#> Install/load packages
packages <- c("terra", "smoothr", "lwgeom", "dplyr")

for(i in 1:NROW(packages)){
  if(!require(
    packages[i], character.only = TRUE
  )){
    install.packages(
      packages[i],
      dependencies = TRUE
    )
    require(
      packages[i],
      character.only = TRUE
    )
  }
}

#>----------------------------------------------------------------------------<|
#> Load data
pa_file <- "high_suitability.shp"
main_dir <- "L:/poppman/shared/dami"
min_polygon_area <- 625

vector_file <- file.path(
  main_dir, "dat", "pat", pa_file
  )

template_file <- file.path(
  main_dir, "dat", "res", "resistance_map.tif"
  )

# Remove road verges is now default
# if (remove_road_verges) {
#   habitat_type <- paste0(habitat_type, "_noTBA")
# }

output_file <- sub(".shp", ".tif", vector_file)

polygons_source <- terra::vect(vector_file)
polygons_source <- terra::disagg(polygons_source)

cat(
  "Number of patches within ZH:",
  length(polygons_source),
  "\n"
  )

#>----------------------------------------------------------------------------<|
#> Clean up source polygon layer
min_area <- units::set_units(min_polygon_area, m^2)
polygons <- smoothr::drop_crumbs(
  polygons_source,
  threshold = min_area,
  drop_empty = TRUE
  )

patchsizes <- terra::expanse(polygons, unit = "ha")

#>----------------------------------------------------------------------------<|
#> Load template raster
template <- terra::rast(template_file)

#>----------------------------------------------------------------------------<|
#> Rasterise polygons
raster <- terra::rasterize(
  polygons,
  template,
  touches = TRUE,
  filename = output_file,
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2S",
    gdal = c("COMPRESS=DEFLATE", "TFW=YES")
    )
  )

vector <- terra::as.polygons(raster)

sizes_25 <- terra::expanse(terra::disagg(vector))

sink(sub(".shp", ".info", vector_file))
cat("Info on ", vector_file, "\n")
cat("Protected areas within the original file:\n")
cat("Mean patchsize:", mean(patchsizes), " ha\n")
cat("Median patchsize:", median(patchsizes), "ha\n")
cat("Number of patches:", length(patchsizes), "\n\n")
cat("Protected areas within the exported raster file:\n")
cat(
  "Mean patch size (25*25 m):", mean(
    sizes_25[sizes_25 < 1000000000]
    ), "ha\n"
  )
cat(
  "Median patch size (25*25 m):", median(
    sizes_25[sizes_25 < 1000000000]
  ), "ha\n"
)
cat("Number of patches (25*25 m):", length(sizes_25), "\n")
sink()