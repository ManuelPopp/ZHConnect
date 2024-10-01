#!/usr/bin/env/Rscript --vanilla
#>------------------------------------<
#
# Script name: create_patchsets
#
# Author: Manuel R. Popp
# Email: manuel.popp@wsl.ch
#
# Date Created: 2023-09-04
#
# ---------------------------
#
# Description: Add patches to the existing infrastructure until target area is
# reached. This procedure is used to create samples (hypothetical future
# ecological infrastructure to generate graphs and calculate metrics on.)
#
# Notes: -
#
#>----------------------------------------------------------------------------<|
#> Install/load packages
packages <- c("tools", "optparse", "stringr", "dplyr", "terra", "tidyterra")

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
#> Settings
CLASSVALUE <- 1

#>----------------------------------------------------------------------------<|
#> Options
option_list = list(
  make_option(c("-t", "--habitat"), type = "character",
              help = "Habitat type."
  ),
  make_option(c("-n", "--n_samples"), type = "numeric", default = 5000,
              help = "Number subsetsets to generate."
  ),
  make_option(c("-p", "--patchfile"), type = "character",
              default = "patches_1.3_ha_dj.tif",
              help = "GeoTIFF encoding the patches."
  )
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
#opt <- list(area = 530, habitat = "test", n_samples = 3300, patchfile = "Subset_1ha.tif")

#>----------------------------------------------------------------------------<|
#> Functions
get_file_location <-  function(){
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value,
                    into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)

  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }else if (length(this_file < 3)) {# quite arbitrary; might produce unexpected behaviour
    this_file <- getwd()
  }
  
  return(this_file)
}

#>----------------------------------------------------------------------------<|
#> Settings
## Parameters
### Habitat type
HABITATTYPE <- opt$habitat
NSAMPLES <- opt$n_samples
PATCHFILE <- opt$patchfile

## Directories
dir_this_file <- get_file_location()
if (Sys.info()["sysname"] == "Windows") {
  dir_this_file <- "L:/poppman/shared/dami/rsc"
} else {
  dir_this_file <- "/lud11/poppman/shared/dami/rsc"
}

dir_main <- dirname(dir_this_file)

dir_dat <- file.path(dir_main, "dat")
dir_pat <- file.path(dir_dat, "pat")
dir_spl <- file.path(dir_dat, "spl", HABITATTYPE)
dir.create(dir_spl, showWarnings = TRUE, recursive = TRUE, mode = "0777")

#>----------------------------------------------------------------------------<|
#> Select additional patches
#> Patches are drawn in sets NSETS times.

## Load patches raster
patches_file <- file.path(dir_dat, "pat", PATCHFILE)
patches_raster <- terra::rast(patches_file)
print(paste("Patches:", patches_file))

## Load existing infrastructure
if (tools::file_ext(HABITATTYPE) == "") {
  print(
  paste(
    "Current infrastructure",
    file.path(
      dir_pat, paste0(HABITATTYPE, ".tif")
      )
    )
  )

  infrastructure <- terra::rast(
    file.path(
      dir_pat, paste0(HABITATTYPE, ".tif")
    )
  )
} else if (tools::file_ext(HABITATTYPE) == "tif") {
  print(
  paste(
    "Current infrastructure",
    file.path(
      dir_pat, HABITATTYPE
      )
    )
  )

  infrastructure <- terra::rast(
    file.path(
      dir_pat, HABITATTYPE
    )
  )
} else {
  print(
  paste(
    "Current infrastructure",
    file.path(
      dir_pat, HABITATTYPE
      )
    )
  )

  infrastructure <- terra::vect(
    file.path(
      dir_pat, HABITATTYPE
    )
  ) %>%
  terra::rasterize(patches_raster)
}

# Load infrastructure as shapefile
if (tools::file_ext(HABITATTYPE) == "shp") {
  infrastructure_shp <- terra::vect(
    file.path(
      dir_pat, HABITATTYPE
    )
  )
} else {
  infrastructure_shp <- terra::vect(
    file.path(
      dir_pat, paste0(HABITATTYPE, ".shp")
    )
  )
}

CUMAREA <- sum(terra::expanse(infrastructure_shp, unit = "ha"))
print(
  paste0(
    "Settings: Habitat=",
    HABITATTYPE, " Area=",
    CUMAREA,
    " N samples=",
    NSAMPLES,
    " Patch file=", PATCHFILE
  )
)

#HABITATTYPE <- tools::file_path_sans_ext(HABITATTYPE)

# Patches to polygons
patches_polygon <- terra::as.polygons(patches_raster)
names(patches_polygon) <- "POLYID"
print(paste("Found", length(patches_polygon), "polygons."))

polygon_ids <- patches_polygon$POLYID
polygon_areas <- terra::expanse(patches_polygon, unit = "ha")
polygons <- data.frame(id = polygon_ids, area = polygon_areas)

## Draw samples
sample_meta <- list()
remaining <- polygons
n <- 0

pb = txtProgressBar(min = 0, max = NSAMPLES, initial = 0)

for (i in 1:NSAMPLES) {
  pol_ids <- c()
  pol_areas <- c()
  cum_area <- 0

  ## Collect sample
  #print("Collecting patches...")
  #flush.console()

  while (cum_area < CUMAREA) {
    random_rownum <- sample(
      seq(1, nrow(remaining))[which(!(remaining$id %in% pol_ids))],
      1, replace = FALSE
      )
    row <- remaining[random_rownum, ]

    if (length(pol_ids) > 0) {
      if (!(row$id %in% pol_ids)) {
        pol_ids <- c(pol_ids, row$id)
        pol_areas <- c(pol_areas, row$area)
      }
    }else{
      pol_ids <- row$id
      pol_areas <- row$area
    }

    cum_area <- sum(pol_areas)
    #print(paste("Cumulative area =", cum_area))
    
    remaining <- remaining[-random_rownum,]
    
    if (nrow(remaining) == 0) {
      n <- n + 1
      print(paste("Entire patch set was sampled", n, "times."))
      remaining <- polygons
    }
  }

  ## Write polygon subset
  print(paste("Subset", i, "contains", length(pol_ids), "patches."))
  
  polygon_subset <- patches_polygon[patches_polygon$POLYID %in% pol_ids]

  sample_path <- file.path(dir_spl, paste0("sample_", as.character(i), ".tif"))

  tmp <- terra::mask(
    infrastructure,
    mask = polygon_subset,
    inverse = TRUE,
    updatevalue = CLASSVALUE,
    touches = FALSE
  )
  
  terra::writeRaster(
    tmp,
    filename = sample_path,
    datatype = "INT2S",
    gdal = "TFW=YES",
    overwrite = TRUE
    )

  sample_meta[[i]] <- c(
    Path = sample_path,
    n_polygons = length(pol_ids),
    area = cum_area,
    polygon_ids = paste(pol_ids, collapse = ":")
    )
  
  setTxtProgressBar(pb, i)
}

meta <- as.data.frame(
  do.call(rbind, sample_meta)
)

# Write metadata to file (Unix eol is important for reading in bash!)
con <- file(file.path(dir_spl, "patch_meta.csv"))

if( !isOpen(con = con, rw = "wb") ) {
  open(con, open = "wb" )
}

write(x = paste(colnames(meta), collapse = ","), file = con)
for (i in 1:nrow( meta )) {
  write(paste(meta[i, ], collapse = ","), file = con)
}

close(con)

terra::writeRaster(
  infrastructure,
  file.path(dir_spl, paste0("Current_Infrastructure.tif")),
  datatype = "INT2S",
  gdal = "TFW=YES",
  overwrite = TRUE
)

fileConn <-file(file.path(dir_spl, "patch_meta.info"))
writeLines(c(
  paste("Habitat:", HABITATTYPE),
  paste("Cumulative_area:", CUMAREA),
  paste("N_samples:", NSAMPLES),
  paste("Patchfile:", PATCHFILE)
), fileConn)
close(fileConn)